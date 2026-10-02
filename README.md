# annogenome 

- my genome annotation pipeline releasing today.
- Entire pipeline using Austus, Braker, TSEBRA, HISAT2 and annotates complete genome with evidence. 
- Evaluates completeness. 


```
./genome_annotation_pipeline.sh -g genome.fa -s Species_name \
  -1 rnaseq_R1.fq.gz -2 rnaseq_R2.fq.gz -p proteins.fa \
  -o ./annotation_out -t 24

```


```
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


```

Gaurav Sablok \
gsablok@proton.me
