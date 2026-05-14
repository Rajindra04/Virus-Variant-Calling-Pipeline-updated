# Detailed Installation Guide

## Table of Contents
1. [Quick Start](#quick-start)
2. [System Requirements](#system-requirements)
3. [Step-by-Step Installation](#step-by-step-installation)
4. [Java Configuration](#java-configuration)
5. [Verification](#verification)
6. [Common Issues and Solutions](#common-issues-and-solutions)

---

## Quick Start

For users with working conda and basic familiarity:

```bash
git clone https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline-updated.git
cd Virus-Variant-Calling-Pipeline-updated
conda env create -f environment.yml
conda activate dengue_pipeline
pip install -r requirements.txt
pip install .
java -version  # Verify Java is available
````
System Requirements
Operating System
Linux (Ubuntu 18.04+, CentOS 7+) or macOS 10.14+
Windows: Use WSL2 (Windows Subsystem for Linux)
Hardware
Minimum: 8 GB RAM, 20 GB disk space
Recommended: 16+ GB RAM, 50+ GB disk space (for multiple samples)
Software
Conda (Miniconda or Anaconda)
Git
Java 11+ (will be installed via conda)
Step-by-Step Installation
Step 1: Install Miniconda (if not already installed)
Linux:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
macOS:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```
Step 2: Clone Repository
```bash
git clone https://github.com/Rajindra04/Virus-Variant-Calling-Pipeline-updated.git
cd Virus-Variant-Calling-Pipeline-updated
```
Step 3: Create Conda Environment
```bash
conda env create -f environment.yml
conda activate dengue_pipeline
Step 4: Install Python Package
```
```bash
pip install -r requirements.txt
pip install .
```
Step 5: Verify Installation
```bash
which bwa-mem2 samtools fastp fastqc gatk snpeff snpsift ivar bcftools
python --version  # Should output Python 3.11.x
java -version    # Should output Java 11 or higher
```
Java Configuration
Default Setup (Recommended)
The conda environment provides Java 11+:

```bash
java -version
```
Custom Java Paths (If Needed)
```bash
run_pipeline \
  --input_dir fastq/ \
  --reference_fasta references/NC_001477.1.fasta \
  --genbank_file NC_001477.1.gb \
  --output_dir output/ \
  --config configs/denv1.yaml \
  --gatk_java /usr/lib/jvm/java-17-openjdk-amd64/bin/java \
  --snpeff_java /usr/lib/jvm/java-11-openjdk-amd64/bin/java
```
Common Issues
Issue: "Unsupported major.minor version"
bash
# Use correct Java version
```
run_pipeline ... --gatk_java /path/to/java-17/bin/java
```
Issue: "conda: command not found"
```bash
exec bash
source ~/miniconda3/bin/activate
Issue: Environment creation fails
```
```bash
conda clean --all
conda env remove -n dengue_pipeline
conda env create -f environment.yml -v
```


