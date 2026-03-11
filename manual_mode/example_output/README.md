# Example Output Files

These files show what to expect at each QC checkpoint. Compare your results against these examples.

| File | Description | What to look for |
|------|-------------|------------------|
| `fastp_qc_good.txt` | Key fastp metrics from a passing sample | >80% survival, >85% Q30 |
| `fastp_qc_bad.txt` | Key fastp metrics from a failing sample | <70% survival or Q30 |
| `flagstat_example.txt` | Annotated samtools flagstat output | Mapping rate interpretation |
| `annotation_example.tsv` | 10 rows of normal variant annotation | Synonymous + missense, 0 frameshifts |
| `annotation_with_frameshift.tsv` | THE KEY LESSON — what bad data looks like | Multiple frameshifts = artifact |
| `variant_classification_example.tsv` | Pass 1 vs Pass 2 comparison | FIXED vs INTRAHOST categories |
| `expected_qc_summary.txt` | Full QC output for a good sample | All checkpoints passing |
