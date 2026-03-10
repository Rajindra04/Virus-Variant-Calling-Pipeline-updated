# Virus Variant Calling Pipeline

This pipeline processes paired-end FASTQ files to perform variant calling and generate consensus sequences for viral genomes, specifically designed for Dengue virus (DENV1). It uses a series of bioinformatics tools to map reads, convert SAM to BAM, call variants, annotate them with SnpEff, and summarize results.

### Modules 1–6: Core Pipeline
![Pipeline Modules 1-6](docs/pipeline_modules_1-6.png)

### Module 7: Self-Reference Variant Calling
![Pipeline Module 7](docs/pipeline_module_7.png)

## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Contact](#contact)

## Overview

The pipeline performs the following steps:
1. **Create Sample Sheet**: Generates a `samplesheet.tsv` from FASTQ files.
2. **Map Reads**: Trims reads with `fastp`, runs `FastQC`, and maps reads to a reference using `bwa-mem2`.
3. **SAM to BAM Conversion**: Converts SAM files to sorted and indexed BAM files using `samtools`.
4. **SnpEff Database Creation**: Builds a SnpEff database from a GenBank file.
5. **Variant Calling and Consensus**: Performs variant calling with `ivar` and `GATK`, generating consensus sequences.
6. **Summarization**: Summarizes coverage, consensus FASTA, and SnpEff annotations.

## Prerequisites

- **Operating System**: Linux (tested on Ubuntu) or macOS.
- **Conda**: Miniconda or Anaconda installed.
- **Input Files**:
  - Paired-end FASTQ files (e.g., `fastq/D1-1_S1_L001_R1_001.fastq.gz`, `fastq/D1-1_S1_L001_R2_001.fastq.gz`).
  - Reference FASTA file (e.g., `references/denv1.fasta`).
  - GenBank file for SnpEff (e.g., `denv1.gb`, downloadable from NCBI: `NC_001477.1`).

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline.git
   cd Virus-Variant-Calling-Pipeline
   ```

2. **Set Up Conda Environment**:
   Create the `dengue_pipeline` environment using the provided `environment.yml`:
   ```bash
   conda env create -f environment.yml
   conda activate dengue_pipeline
   conda install -c conda-forge openblas
   pip install -r requirements.txt
   pip install .
   ```
   Install Pip manually if above commands does not work and run the commands again
   ```bash
   sudo apt-get apt-get install pip3
   ```

4. **Verify Tools**:
   Ensure all required tools are installed:
   ```bash
   which bwa-mem2 samtools fastp fastqc gatk snpeff snpsift ivar bcftools
   python --version  # Should output Python 3.11.x
   ```

5. **Download GenBank File (if not provided)**:
   ```bash
   wget -O denv1.gb "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_001477.1&rettype=gb&retmode=text"
   ```

## Usage

1. **Prepare Input Files**:
   - Place paired-end FASTQ files in a `fastq/` directory.
   - Ensure the reference FASTA (`references/denv1.fasta`) and GenBank file (`denv1.gb`) are in the repository root or a specified directory.

2. **Run the Pipeline with test data**:
   ```bash
   run_pipeline \
     --input_dir fastq/ \
     --reference_fasta references/denv1.fasta \
     --genbank_file denv1.gb \
     --output_dir output/ \
     --database_name denv1
   ```

3. **Output Files**:
   - `output/samplesheet.tsv`: Sample sheet with FASTQ file paths.
   - `output/sam_files/*.sam`: SAM files from read mapping.
   - `output/*.sorted.bam`: Sorted and indexed BAM files.
   - `output/*.vcf`: Variant call files.
   - `output/*_snpEff_genes.txt`: SnpEff annotation summaries.
   - `output/summary_table.csv`: Summary of variants and annotations.
   - `output/chart_data.json`: Data for visualization.

## Directory Structure

```
Virus-Variant-Calling-Pipeline/
├── fastq/                    # Input FASTQ files
├── references/               # Reference FASTA file (e.g., denv1.fasta)
├── denv1.gb                  # GenBank file for SnpEff
├── output/                   # Output directory
│   ├── sam_files/            # SAM files from map_reads.py
│   ├── *.sorted.bam          # Sorted BAM files
│   ├── *.vcf                 # Variant call files
│   ├── *_snpEff_genes.txt    # SnpEff annotation files
│   ├── summary_table.csv     # Summary table
│   └── chart_data.json       # Visualization data
├── virus_pipeline/           # Pipeline scripts
│   ├── create_samplesheet.py
│   ├── map_reads.py
│   ├── samtobamdenv.py
│   ├── create_snpeff_database.py
│   ├── sam2consensus_test2_ivar.py
│   ├── variant_calling_consensus.py
│   ├── summarize_result.py
│   └── summarize_snpEff.py
├── environment.yml           # Conda environment file
├── requirements.txt          # Additional Python dependencies (if needed)
├── conda-recipe/             # Conda package recipe
│   └── meta.yaml
└── run_pipeline.py           # Main pipeline script
```

## Troubleshooting

- **No BAM Files in `output/`**:
  - Check if SAM files exist in `output/sam_files/`:
    ```bash
    ls -l output/sam_files/
    ```
  - Run `map_reads.py` manually to verify SAM file generation:
    ```bash
    python virus_pipeline/map_reads.py --samplesheet output/samplesheet.tsv --reference references/denv1.fasta
    ```
  - Ensure `samtobamdenv.py` finds SAM files:
    ```bash
    python virus_pipeline/samtobamdenv.py --input_dir output/sam_files --reference_fasta references/denv1.fasta --output_dir output
    ```

- **SnpEff Database Errors**:
  - Verify the GenBank file (`denv1.gb`) is valid and matches the reference FASTA.
  - Check SnpEff logs in `output/snpEff.config`.

- **Dependency Issues**:
  - Ensure Conda channels are configured:
    ```bash
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --set channel_priority strict
    ```
  - Recreate the environment if needed:
    ```bash
    conda env remove -n dengue_pipeline
    conda env create -f environment.yml
    ```

- **File Permission Issues**:
  ```bash
  chmod -R u+rwX output/
  ```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For issues or questions, please contact the maintainer at:
- GitHub: [Rajindra04](https://github.com/Rajindra04)
- Email: [Add your email or preferred contact method]
