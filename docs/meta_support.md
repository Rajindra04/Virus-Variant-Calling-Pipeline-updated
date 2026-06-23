Metagenomic mode (host removal, assembly, and classification)

This document explains how to use the pipeline's new --meta mode and the external tools and references you need.

1) Overview
- When you run the pipeline with --meta, the following steps are performed after fastp trimming and QC:
  a) Host removal: trimmed reads are mapped to a combined host reference (human + mosquito) using bwa-mem2. Reads where both mates are unmapped are retained (host-depleted paired FASTQs).
  b) De novo assembly: host-depleted reads are assembled (MEGAHIT or metaSPAdes) to produce contigs enriched for non-host/viral sequences.
  c) Classification: contigs are classified using Kraken2 (if kraken2_db is available). Reports are written under the sample assembly directory.
  d) Downstream mapping: host-depleted FASTQs are used for the standard mapping and variant-calling workflow.

2) Preparing host reference
- Concatenate human and mosquito FASTA files to build a combined host reference:

  cat human.fasta mosquito.fasta > references/host_combined.fasta

- Index the combined reference for bwa-mem2:

  bwa-mem2 index references/host_combined.fasta

- Note: keep the original separate FASTAs if you want to compute per-host mapping statistics.

3) Installing required tools
- The pipeline already requires bwa-mem2 and samtools. For metagenomic mode, install:
  - megahit (or metaSPAdes)
  - kraken2 (optional, for classification)

Using conda (recommended):

  conda install -c bioconda megahit kraken2 spades

Kraken2 databases are large and must be downloaded separately. See Kraken2 documentation: https://github.com/DerrickWood/kraken2

4) Config
- Add a `meta` section to your existing YAML config (see configs/meta_example.yaml in this branch).
- Required field when using --meta: meta.host_reference
- Optional fields: assembler, min_contig_len, kraken2_db, run_read_level_classification

5) Resource considerations
- Assemblers and Kraken2 DBs can be memory- and disk-intensive. Run on a machine with sufficient RAM and disk space.
- For large datasets, consider skipping assembly (set assembler to null) and using read-level classification instead (kraken2 on reads) to get taxonomic composition faster.

6) Output
- host-depleted FASTQs: {sample}_output/{sample}_host_depleted_1.fastq.gz and _2.fastq.gz
- assembly dir: {sample}_output/assembly_{sample}/ (contains contigs and kraken2 reports when run)
- Logs include the assembler and classifier commands for provenance.

7) Example run

  run_pipeline \
    --samplesheet samples.tsv \
    --reference references/NC_001477.1.fasta \
    --config configs/config_with_meta.yaml \
    --threads 8 \
    --meta

