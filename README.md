---

## License

This project is licensed under the terms of the MIT License. See the [LICENSE](LICENSE) file for comprehensive details.

## Contact

For production pipeline questions, environment bugs,The updated file text is already written completely inside a copyable **Markdown block** above! 

However, if you want the absolute raw, unformatted Markdown text (without any outer code block container) so you can copy and paste it directly into your `README.md` file, here it is:

---

# Virus Variant Calling Pipeline

This pipeline processes paired-end FASTQ files to perform variant calling and generate consensus sequences for viral genomes, specifically optimized for Dengue virus (DENV1-3). It integrates a robust suite of bioinformatics tools to map reads, track quality control profiles, process structural alignments, and annotate discovered variants.

![Virus Variant Calling Pipeline](docs/pipeline_figure.png)

## Manual Mode Tutorial

**New users start here:** The [Manual Mode Tutorial](manual_mode/TUTORIAL.md) walks you through every individual execution stage with QC checkpoints, strict parameter explanations, and core biological context. It teaches you how to evaluate whether output trends are reasonable before publishing downstream results.

## Table of Contents
- [Manual Mode Tutorial](#manual-mode-tutorial)
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Memory Requirements](#memory-requirements)
- [Parallelization Performance](#parallelization-performance)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Overview

The workflow processes raw read sequences through the following sequential modules:
1. **Create Sample Sheet**: Automatically scans your data and constructs a clean index map (`samplesheet.tsv`) containing absolute file locations.
2. **Map Reads**: Performs quality trimming via `fastp`, generates absolute quality distribution plots with `FastQC`, and aligns reads using `bwa-mem2`.
3. **SAM to BAM Conversion**: Compresses raw sequence layouts into sorted, indexed BAM files using `samtools`.
4. **SnpEff Database Creation**: Dynamically assembles a localized annotation infrastructure directly from structural GenBank records.
5. **Variant Calling and Consensus**: Executes modern variant discovery using `bcftools` (and optionally `GATK`), applies strict quality filters, and constructs complete `.fasta` consensus sequences.
6. **Summarization**: Collates coverage depths, consensus assets, and SnpEff annotations into polished cross-sample spreadsheets.

---

## Prerequisites

- **Operating System**: Linux (Tested extensively on Ubuntu) or macOS.
- **Environment Engine**: Miniconda or Anaconda.
- **Input Structures**:
  - Paired-end FASTQ files (e.g., `sample_R1_001.fastq.gz` / `sample_R2_001.fastq.gz`).
  - Reference genome sequence in standard FASTA layout.
  - Corresponding NCBI GenBank feature flatfile (e.g., `NC_001477.1.gb`).

### ⚠️ Critical Runtime Dependency: Java Environments
Variant calling and genomic structural annotation rely on distinct execution engines that depend on different Java runtimes. `GATK` pipelines operate optimally under **Java 8 or 17**, whereas modern instances of `SnpEff` require **Java 11 to 21**. 

To prevent systemic environment conflicts, you can explicitly route execution to explicit paths using the following runtime flags:
* `--gatk_java`: Absolute path pointing directly to your GATK-compatible Java binary.
* `--snpeff_java`: Absolute path pointing directly to your SnpEff-compatible Java binary.

---

## Installation

### 1. Clone the Source Repository
```bash
git clone [https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline-updated.git](https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline-updated.git)
cd Virus-Variant-Calling-Pipeline-updated
```
2. Configure Environment Contexts
Build the foundational conda infrastructure and map editable local installation bindings:

```Bash
conda env create -f environment.yml
conda activate dengue_pipeline
pip install -r requirements.txt
pip install -e .
```
3. Provision Native System Java Runtimes (If Required)
If your default conda environment path cannot handle multi-version execution natively, install the required JDK packages globally on your machine:

On Ubuntu / Debian systems:

```Bash
sudo apt-get update
sudo apt-get install openjdk-17-jdk openjdk-21-jdk

# Check and verify valid path mappings on your local system
update-alternatives --list java
```
On macOS environments (via Homebrew):

```Bash
brew install openjdk@17 openjdk@21
```

# Locate executable binary points
```
ls /usr/local/opt/openjdk@17/bin/java
ls /opt/homebrew/opt/openjdk@21/bin/java
```
Usage
1. Structure Working Materials
Place raw sequencing datasets inside the fastq/ directory space.

Verify target references and feature annotations exist within the references/ workspace.

2. Primary Execution Commands
Standard Execution (Explicit Dual-Java Paths & Trimming Active)
```Bash
run_pipeline \
  --input_dir fastq/ \
  --reference_fasta references/NC_001477.1.fasta \
  --primer_bed primers/denv1_primers.bed \
  --genbank_file NC_001477.1.gb \
  --output_dir output/ \
  --config configs/denv1.yaml \
  --gatk_java /usr/lib/jvm/java-17-openjdk-amd64/bin/java \
  --snpeff_java /usr/lib/jvm/java-21-openjdk-amd64/bin/java \
  --parallel
```
Optional Flag Adaptations
Skipping Primer Trimming: Amplicon primer sequence removal is optional. Pass a valid target coordinates file via --primer_bed to perform ivar trim. If you omit this flag entirely, the system skips trimming and safely parses raw BAM outputs. This is ideal for untargeted metagenomic datasets:

```Bash
# Omitting --primer_bed completely avoids amplicon processing modifications
run_pipeline --input_dir fastq/ --reference_fasta ref.fa --output_dir output/ --config configs/denv1.yaml
```
Config-Driven Annotations: If you lack a structural .gb file asset, pass --annotation_mode config to run lightweight variant filtering using configuration rules rather than deep SnpEff database lookups.

Multithreading Throughput: Pass the --parallel flag to automatically check system architectures and forward optimal computing blocks to intensive core dependencies (fastp, fastqc, bwa-mem2, and samtools).

```Directory Structure
Virus-Variant-Calling-Pipeline-updated/
├── fastq/                    # Raw paired-end sequencing inputs (User-provided)
├── references/               # Fixed reference genome configurations
│   ├── NC_001477.1.fasta     # Dengue Virus Type 1
│   ├── NC_001474.2.fasta     # Dengue Virus Type 2
│   └── NC_001475.2.fasta     # Dengue Virus Type 3
├── configs/                  # Host configuration definitions
│   ├── denv1.yaml
│   ├── denv2.yaml
│   └── denv3.yaml
├── NC_001477.1.gb            # SnpEff reference feature assets
├── output/                   # Auto-generated pipeline deliverables
│   ├── sam_files/            # Clean sequence alignment records (.sam)
│   ├── trimmed_fastq/        # Intermediate processed reads and fastp summaries
│   ├── fastqc_reports/       # Visual read distribution diagnostics
│   ├── bam_files/            # Compressed, sorted alignment maps (.bam)
│   ├── vcf_files/            # Discovered variants (.vcf.gz)
│   └── consensus_sequences/  # Extracted high-confidence viral genomes (.fasta)
├── virus_pipeline/           # Pipeline architectural source files
│   ├── create_samplesheet.py
│   ├── map_reads.py
│   ├── samtobamdenv.py
│   └── variant_calling_consensus.py
├── environment.yml           # Conda environment definition file
├── requirements.txt          # Python runtime requirements
├── setup.py                  # Local package definition file
└── run_pipeline.py           # Primary master pipeline entry point
```
Memory Requirements
The variant processing workflow can be memory-intensive, especially during GATK HaplotypeCaller routines. By default, the system caps Java heap execution allocation flags at 4 GB.

If execution hangs or terminates unexpectedly, check your total system capacity (free -h).

Systems with 8 GB+ RAM: Recommended configuration. Run standard commands.

Systems with limited resources (4 GB - 6 GB RAM): Pass a modified heap limitation directive via the CLI command or alter the YAML configuration:

```Bash
run_pipeline --input_dir fastq/ [options...] --gatk_memory 2g
```
Parallelization Performance
When running with the --parallel flag, the pipeline optimizes performance using a two-tier execution structure:

Tool-Level Parallelism: The pipeline forwards your thread pool to underlying multi-threaded tasks. Tools like fastp, FastQC, and bwa-mem2 utilize parallel threads simultaneously to process your sequence data as fast as your hardware allows.

Sample-Level Linear Progression: Individual samples are processed sequentially. This prevents disk I/O write bottlenecks, keeping your drive clear of read/write collisions while still processing individual data files at maximum computing speed.

Troubleshooting
Mismatched File Errors (KeyError: 'fastq_1')
If the execution halts at Step 2 with an error indicating column indices are missing, your sample sheet formatting is misaligned. Ensure that create_samplesheet.py writes headers matching exactly what map_reads.py checks for: sample_name, fastq_1, and fastq_2. Run head -n 2 output/samplesheet.tsv to verify.

Local Updates Not Registering
If you alter configuration scripts or adjust processing step files locally, but the global command execution runs legacy cached scripts, force an ecosystem update clear out:

```Bash
rm -rf *.egg-info build/ dist/
pip install --force-reinstall -e .
```
License
This project is licensed under the terms of the MIT License. See the LICENSE file for comprehensive details.

Contact
For production pipeline questions, environment bugs, or performance optimization requests, please open an issue in the tracker:

GitHub Repository Main Maintainer: Rajindra04
