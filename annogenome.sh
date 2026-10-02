#!/usr/bin/env bash
###############################################################################
# genome_annotation_pipeline.sh  Gaurav Sablok gsablok@proton.me
#
# End-to-end eukaryotic genome annotation pipeline:
#   1. Genome QC / indexing
#   2. Repeat identification & masking      (RepeatModeler2 + RepeatMasker)
#   3. RNA-seq alignment                    (HISAT2 -> sorted BAM)
#   4. Evidence-based gene prediction       (BRAKER1 / BRAKER2 / BRAKER3)
#   5. Merge independent BRAKER runs        (TSEBRA)   [only if both RNA+protein
#                                                        were run separately]
#   6. Annotation completeness check        (BUSCO, protein mode)
#   7. Extract final outputs                (GFF3, CDS, proteins, stats)
#
# Requires (on PATH or as conda envs, see ENV VARS below):
#   hisat2, hisat2-build, samtools, RepeatModeler, RepeatMasker, BuildDatabase,
#   braker.pl (BRAKER2/3), augustus, gffread, TSEBRA, busco, AGAT (optional),
#   GeneMark-ETP / GeneMark-ES (BRAKER dependency), ProtHint (bundled w/ BRAKER)
#
# Usage:
#   ./genome_annotation_pipeline.sh -g genome.fa -s "Species_name" \
#       -1 rnaseq_R1.fastq.gz -2 rnaseq_R2.fastq.gz \
#       -p proteins.fa -o /path/to/outdir -t 16
#
# Any single evidence type (RNA-seq only, protein only, or both) is supported;
# the script auto-selects BRAKER1 (RNA only), BRAKER2 (protein only), or
# BRAKER3 (RNA+protein together) and runs the recommended single BRAKER3
# call when both are present -- toggle RUN_SEPARATE_AND_MERGE=1 below if you
# instead want to run BRAKER1 + BRAKER2 separately and merge with TSEBRA.
###############################################################################

set -euo pipefail
IFS=$'\n\t'

###############################################################################
# 0. DEFAULTS / ARGUMENT PARSING
###############################################################################
GENOME=""
SPECIES=""
READS_1=""
READS_2=""
READS_SE=""            # optional single-end RNA-seq
PROTEINS=""
OUTDIR="$(pwd)/annotation_run"
THREADS=8
BUSCO_LINEAGE="eukaryota_odb10"
RUN_SEPARATE_AND_MERGE=0   # 1 = run BRAKER1+BRAKER2 separately, merge with TSEBRA
SKIP_REPEATMASK=0
FUNCTIONAL_ANNOT=0         # 1 = attempt InterProScan/eggNOG-mapper if installed

usage() {
  cat <<EOF
Usage: $0 -g genome.fa -s "Genus species" [options]

Required:
  -g FILE     Genome FASTA
  -s STRING   Species identifier (no spaces; used for AUGUSTUS species name)

Evidence (at least one of RNA-seq or protein required):
  -1 FILE     RNA-seq forward reads (fastq[.gz]) - paired
  -2 FILE     RNA-seq reverse reads (fastq[.gz]) - paired
  -U FILE     RNA-seq single-end reads (fastq[.gz])
  -p FILE     Protein evidence FASTA (e.g. related species proteomes / OrthoDB)

Options:
  -o DIR      Output directory (default: ./annotation_run)
  -t INT      Threads (default: 8)
  -b STRING   BUSCO lineage dataset (default: eukaryota_odb10)
  -m          Skip repeat masking (genome already masked)
  -M          Run BRAKER1 + BRAKER2 separately and merge via TSEBRA
              (default: single combined BRAKER3 run when both evidence given)
  -F          Attempt functional annotation (InterProScan / eggNOG-mapper)
  -h          Show this help

Example:
  $0 -g genome.fa -s Drosophila_melanogaster \\
     -1 rnaseq_1.fq.gz -2 rnaseq_2.fq.gz -p refseq_proteins.fa \\
     -o ./dmel_annotation -t 24
EOF
  exit 1
}

while getopts "g:s:1:2:U:p:o:t:b:mMFh" opt; do
  case $opt in
    g) GENOME="$OPTARG" ;;
    s) SPECIES="$OPTARG" ;;
    1) READS_1="$OPTARG" ;;
    2) READS_2="$OPTARG" ;;
    U) READS_SE="$OPTARG" ;;
    p) PROTEINS="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    b) BUSCO_LINEAGE="$OPTARG" ;;
    m) SKIP_REPEATMASK=1 ;;
    M) RUN_SEPARATE_AND_MERGE=1 ;;
    F) FUNCTIONAL_ANNOT=1 ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "$GENOME" || -z "$SPECIES" ]] && { echo "ERROR: -g and -s are required."; usage; }
[[ -z "$READS_1" && -z "$READS_SE" && -z "$PROTEINS" ]] && {
  echo "ERROR: provide RNA-seq reads (-1/-2 or -U) and/or protein evidence (-p)."; usage; }
[[ ! -f "$GENOME" ]] && { echo "ERROR: genome file not found: $GENOME"; exit 1; }

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"
LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"; }
run() { log "CMD: $*"; eval "$@" >>"$LOGFILE" 2>&1; }
need() { command -v "$1" >/dev/null 2>&1 || { log "ERROR: required tool '$1' not found on PATH."; exit 1; }; }

log "=============================================="
log " Genome Annotation Pipeline starting"
log " Genome:    $GENOME"
log " Species:   $SPECIES"
log " Outdir:    $OUTDIR"
log " Threads:   $THREADS"
log "=============================================="

###############################################################################
# 1. GENOME QC / PREP
###############################################################################
GENOME_DIR="$OUTDIR/00_genome"
mkdir -p "$GENOME_DIR"
CLEAN_GENOME="$GENOME_DIR/genome.fa"

log "STEP 1: Preparing and QC'ing genome FASTA"
need samtools

# Simplify headers (BRAKER/AUGUSTUS are picky about long/odd FASTA headers)
awk '{ if ($0 ~ /^>/) { print $1 } else { print $0 } }' "$GENOME" > "$CLEAN_GENOME"
run "samtools faidx '$CLEAN_GENOME'"

N_SEQS=$(grep -c '^>' "$CLEAN_GENOME")
GENOME_SIZE=$(awk '{sum+=$2} END{print sum}' "$CLEAN_GENOME.fai")
log "Genome contains $N_SEQS sequences, total length $GENOME_SIZE bp"

###############################################################################
# 2. REPEAT IDENTIFICATION & MASKING
###############################################################################
MASKED_GENOME="$CLEAN_GENOME"

if [[ "$SKIP_REPEATMASK" -eq 0 ]]; then
  log "STEP 2: Repeat modeling and masking (RepeatModeler2 + RepeatMasker)"
  need BuildDatabase
  need RepeatModeler
  need RepeatMasker

  REPEAT_DIR="$OUTDIR/01_repeats"
  mkdir -p "$REPEAT_DIR"
  cd "$REPEAT_DIR"

  run "BuildDatabase -name '${SPECIES}_db' '$CLEAN_GENOME'"
  run "RepeatModeler -database '${SPECIES}_db' -threads $THREADS -LTRStruct"

  DENOVO_LIB="$REPEAT_DIR/${SPECIES}_db-families.fa"
  if [[ ! -s "$DENOVO_LIB" ]]; then
    log "WARNING: RepeatModeler library not found, falling back to RepeatMasker's built-in repeat db"
    run "RepeatMasker -pa $THREADS -xsmall -gff -dir '$REPEAT_DIR' '$CLEAN_GENOME'"
  else
    run "RepeatMasker -pa $THREADS -xsmall -gff -lib '$DENOVO_LIB' -dir '$REPEAT_DIR' '$CLEAN_GENOME'"
  fi

  MASKED_GENOME="$REPEAT_DIR/$(basename "$CLEAN_GENOME").masked"
  if [[ ! -s "$MASKED_GENOME" ]]; then
    log "ERROR: RepeatMasker did not produce a masked genome; check $LOGFILE"
    exit 1
  fi
  cd "$OUTDIR"
  log "Soft-masked genome: $MASKED_GENOME"
else
  log "STEP 2: Skipped (assuming genome is already masked, -m given)"
fi

###############################################################################
# 3. RNA-SEQ ALIGNMENT (HISAT2)
###############################################################################
RNASEQ_BAM=""
if [[ -n "$READS_1" || -n "$READS_SE" ]]; then
  log "STEP 3: Aligning RNA-seq reads with HISAT2"
  need hisat2
  need hisat2-build
  need samtools

  ALIGN_DIR="$OUTDIR/02_rnaseq_alignment"
  mkdir -p "$ALIGN_DIR"
  HISAT_INDEX="$ALIGN_DIR/${SPECIES}_hisat2_index"

  run "hisat2-build -p $THREADS '$MASKED_GENOME' '$HISAT_INDEX'"

  SAM_OUT="$ALIGN_DIR/rnaseq.sam"
  if [[ -n "$READS_1" && -n "$READS_2" ]]; then
    run "hisat2 -p $THREADS --dta -x '$HISAT_INDEX' -1 '$READS_1' -2 '$READS_2' -S '$SAM_OUT'"
  elif [[ -n "$READS_SE" ]]; then
    run "hisat2 -p $THREADS --dta -x '$HISAT_INDEX' -U '$READS_SE' -S '$SAM_OUT'"
  fi

  RNASEQ_BAM="$ALIGN_DIR/rnaseq.sorted.bam"
  run "samtools sort -@ $THREADS -o '$RNASEQ_BAM' '$SAM_OUT'"
  run "samtools index '$RNASEQ_BAM'"
  rm -f "$SAM_OUT"
  log "Sorted, indexed RNA-seq BAM: $RNASEQ_BAM"
else
  log "STEP 3: Skipped (no RNA-seq reads supplied)"
fi

###############################################################################
# 4. BRAKER GENE PREDICTION
###############################################################################
need braker.pl
BRAKER_ROOT="$OUTDIR/03_braker"
mkdir -p "$BRAKER_ROOT"

run_braker() {
  # $1 = run name, $2 = extra braker.pl args
  local name="$1"
  local extra_args="$2"
  local dir="$BRAKER_ROOT/$name"
  mkdir -p "$dir"
  log "Running BRAKER ($name) -> $dir"
  run "braker.pl \
      --genome='$MASKED_GENOME' \
      --species='${SPECIES}_${name}' \
      --workingdir='$dir' \
      --threads=$THREADS \
      --gff3 \
      --useexisting=false \
      $extra_args"
}

FINAL_GFF=""

if [[ "$RUN_SEPARATE_AND_MERGE" -eq 1 && -n "$RNASEQ_BAM" && -n "$PROTEINS" ]]; then
  log "STEP 4: Running BRAKER1 (RNA-seq) and BRAKER2 (protein) separately, then merging with TSEBRA"

  run_braker "rnaseq"  "--bam='$RNASEQ_BAM'"
  run_braker "protein" "--prot_seq='$PROTEINS' --epmode"

  need tsebra.py
  MERGE_DIR="$OUTDIR/04_tsebra_merge"
  mkdir -p "$MERGE_DIR"

  RNASEQ_GTF="$BRAKER_ROOT/rnaseq/braker.gtf"
  PROTEIN_GTF="$BRAKER_ROOT/protein/braker.gtf"
  RNASEQ_HINTS="$BRAKER_ROOT/rnaseq/hintsfile.gff"
  PROTEIN_HINTS="$BRAKER_ROOT/protein/hintsfile.gff"

  run "tsebra.py \
        -g '$RNASEQ_GTF','$PROTEIN_GTF' \
        -c /usr/share/tsebra/config/default.cfg \
        -e '$RNASEQ_HINTS','$PROTEIN_HINTS' \
        -o '$MERGE_DIR/braker_combined.gtf'"

  need gffread
  FINAL_GFF="$MERGE_DIR/braker_combined.gff3"
  run "gffread '$MERGE_DIR/braker_combined.gtf' -o '$FINAL_GFF'"

elif [[ -n "$RNASEQ_BAM" && -n "$PROTEINS" ]]; then
  log "STEP 4: Running combined BRAKER3 (RNA-seq + protein evidence in one run)"
  run_braker "braker3" "--bam='$RNASEQ_BAM' --prot_seq='$PROTEINS'"
  FINAL_GFF="$BRAKER_ROOT/braker3/braker.gff3"

elif [[ -n "$RNASEQ_BAM" ]]; then
  log "STEP 4: Running BRAKER1 (RNA-seq evidence only)"
  run_braker "braker1" "--bam='$RNASEQ_BAM'"
  FINAL_GFF="$BRAKER_ROOT/braker1/braker.gff3"

elif [[ -n "$PROTEINS" ]]; then
  log "STEP 4: Running BRAKER2 (protein evidence only)"
  run_braker "braker2" "--prot_seq='$PROTEINS' --epmode"
  FINAL_GFF="$BRAKER_ROOT/braker2/braker.gff3"
fi

if [[ -z "$FINAL_GFF" || ! -s "$FINAL_GFF" ]]; then
  log "ERROR: BRAKER did not produce a final GFF3. Check logs under $BRAKER_ROOT"
  exit 1
fi
log "BRAKER gene models: $FINAL_GFF"

###############################################################################
# 5. TRAINED AUGUSTUS SPECIES PARAMETERS
###############################################################################
# BRAKER trains a new AUGUSTUS species model automatically (named
# ${SPECIES}_<run>), stored in AUGUSTUS_CONFIG_PATH/species/. This lets you
# re-run AUGUSTUS standalone later, e.g. on a new assembly of the same
# species, without re-running the full BRAKER pipeline:
#
#   augustus --species=${SPECIES}_braker3 --gff3=on new_genome.fa > new_predictions.gff3
#
log "STEP 5: AUGUSTUS species parameters trained by BRAKER are available under \$AUGUSTUS_CONFIG_PATH/species/${SPECIES}_*"
log "         Re-use them with: augustus --species=<trained_name> --gff3=on <genome.fa>"

###############################################################################
# 6. EXTRACT PROTEIN / CDS / TRANSCRIPT SEQUENCES
###############################################################################
log "STEP 6: Extracting protein, CDS, and transcript FASTA from final gene set"
need gffread
FINAL_DIR="$OUTDIR/05_final_annotation"
mkdir -p "$FINAL_DIR"

cp "$FINAL_GFF" "$FINAL_DIR/${SPECIES}.annotation.gff3"
run "gffread '$FINAL_DIR/${SPECIES}.annotation.gff3' -g '$MASKED_GENOME' \
      -y '$FINAL_DIR/${SPECIES}.proteins.fa' \
      -x '$FINAL_DIR/${SPECIES}.cds.fa' \
      -w '$FINAL_DIR/${SPECIES}.transcripts.fa'"

N_GENES=$(awk '$3=="gene"' "$FINAL_DIR/${SPECIES}.annotation.gff3" | wc -l)
N_MRNA=$(awk '$3=="mRNA"' "$FINAL_DIR/${SPECIES}.annotation.gff3" | wc -l)
log "Final gene set: $N_GENES genes, $N_MRNA mRNAs"

###############################################################################
# 7. BUSCO COMPLETENESS ASSESSMENT
###############################################################################
if command -v busco >/dev/null 2>&1; then
  log "STEP 7: Running BUSCO (protein mode, lineage=$BUSCO_LINEAGE) on predicted proteins"
  BUSCO_DIR="$OUTDIR/06_busco"
  mkdir -p "$BUSCO_DIR"
  cd "$BUSCO_DIR"
  run "busco -i '$FINAL_DIR/${SPECIES}.proteins.fa' \
        -l '$BUSCO_LINEAGE' \
        -o '${SPECIES}_busco' \
        -m protein \
        -c $THREADS -f"
  cd "$OUTDIR"
  log "BUSCO results: $BUSCO_DIR/${SPECIES}_busco"
else
  log "STEP 7: Skipped (busco not found on PATH) -- install with 'conda install -c bioconda busco' to validate completeness"
fi

###############################################################################
# 8. OPTIONAL FUNCTIONAL ANNOTATION
###############################################################################
if [[ "$FUNCTIONAL_ANNOT" -eq 1 ]]; then
  log "STEP 8: Functional annotation requested (-F)"
  if command -v interproscan.sh >/dev/null 2>&1; then
    FUNC_DIR="$OUTDIR/07_functional"
    mkdir -p "$FUNC_DIR"
    run "interproscan.sh -i '$FINAL_DIR/${SPECIES}.proteins.fa' \
          -f TSV,GFF3 -goterms -pa -cpu $THREADS -d '$FUNC_DIR'"
    log "InterProScan output written to $FUNC_DIR"
  else
    log "WARNING: -F given but interproscan.sh not found on PATH; skipping functional annotation"
  fi
  if command -v emapper.py >/dev/null 2>&1; then
    FUNC_DIR="$OUTDIR/07_functional"
    mkdir -p "$FUNC_DIR"
    run "emapper.py -i '$FINAL_DIR/${SPECIES}.proteins.fa' \
          --cpu $THREADS -o '${SPECIES}_eggnog' --output_dir '$FUNC_DIR'"
    log "eggNOG-mapper output written to $FUNC_DIR"
  fi
else
  log "STEP 8: Skipped (pass -F to attempt InterProScan/eggNOG-mapper functional annotation)"
fi

###############################################################################
# 9. SUMMARY
###############################################################################
SUMMARY="$OUTDIR/PIPELINE_SUMMARY.txt"
{
  echo "Genome Annotation Pipeline Summary"
  echo "==================================="
  echo "Date:              $(date)"
  echo "Species:           $SPECIES"
  echo "Input genome:      $GENOME"
  echo "Masked genome:     $MASKED_GENOME"
  echo "RNA-seq BAM:       ${RNASEQ_BAM:-none}"
  echo "Protein evidence:  ${PROTEINS:-none}"
  echo "Final GFF3:        $FINAL_DIR/${SPECIES}.annotation.gff3"
  echo "Proteins FASTA:    $FINAL_DIR/${SPECIES}.proteins.fa"
  echo "CDS FASTA:         $FINAL_DIR/${SPECIES}.cds.fa"
  echo "Transcripts FASTA: $FINAL_DIR/${SPECIES}.transcripts.fa"
  echo "Gene count:        $N_GENES"
  echo "mRNA count:        $N_MRNA"
  echo "BUSCO dir:         $OUTDIR/06_busco (if run)"
  echo "Full log:          $LOGFILE"
} | tee "$SUMMARY"

log "=============================================="
log " Pipeline complete. See $SUMMARY"
log "=============================================="
