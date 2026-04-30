# Virus Variant Calling Pipeline — Manual Mode Tutorial

## Why This Tutorial Exists

This pipeline was previously run as a single automated command. That approach produced results containing **multiple frameshift variants per genome** — a biological impossibility for Dengue virus, which encodes all its proteins in a single polyprotein reading frame. A real frameshift would destroy every downstream protein, producing a non-viable virus. Multiple frameshifts in annotated output almost certainly indicate sequencing artifacts, mapping errors, or incorrect variant calling parameters.

**This tutorial teaches you to run the pipeline step by step**, with QC checkpoints at every stage, so you can:

1. **Understand** what each tool does and why
2. **Evaluate** whether output is biologically reasonable
3. **Catch problems** before they propagate to final results
4. **Never publish** biologically impossible results

---

## Prerequisites

### Software
- **Conda** (Miniconda or Anaconda)
- All bioinformatics tools are installed via `environment.yml` — no manual tool installation needed

### Disk Space
- ~500 MB for raw FASTQ downloads (all 40 samples)
- ~5 GB peak working space during processing
- ~2 GB final output

### Network
- Internet access for SRA downloads (Step 1 only)

### Time
- ~1-2 hours for a single sample end-to-end
- ~1-2 days for all 40 samples (can run in parallel)

---

## Quick Start: One Sample End-to-End

```bash
# Clone the repo and enter it
git clone https://github.com/mihinduk/Virus-Variant-Calling-Pipeline-updated.git
cd Virus-Variant-Calling-Pipeline-updated

# Set the serotype (denv1, denv2, or denv3)
export SEROTYPE=denv2

# Run setup (creates conda env, checks tools, indexes reference)
cd manual_mode
bash 00_setup.sh

# Download one test sample
bash 01_fetch_data.sh SRR35818859

# Create samplesheet
bash 02_create_samplesheet.sh

# Process the sample (each step checks the previous step's output)
bash 03_qc_trimming.sh SRR35818859
bash 04_mapping.sh SRR35818859
bash 05_sam_to_bam.sh SRR35818859
bash 06_dedup_primer_trim.sh SRR35818859

# Build annotation database (once per serotype)
bash 07_snpeff_database.sh

# Continue processing
bash 08_coverage_consensus.sh SRR35818859
bash 09_variant_calling.sh SRR35818859
bash 10_annotation.sh SRR35818859        # ← THE CRITICAL CHECKPOINT

# Pass 2: self-reference analysis
bash 11_self_reference_prep.sh SRR35818859
bash 12_self_reference_calling.sh SRR35818859
bash 13_variant_comparison.sh SRR35818859

# Protein extraction (processes all samples in output directory)
bash 14_protein_extraction.sh
```

---

## Module Guide

### Module 0: Setup (`00_setup.sh`)

**Run first.** Sets environment variables, creates directories, checks tools.

Edit the `SEROTYPE` variable at the top (or `export SEROTYPE=denv2` before running). All downstream scripts source this file automatically.

```bash
# Process DENV1 instead of the default DENV2:
export SEROTYPE=denv1
bash 00_setup.sh
```
Run first. Sets environment variables, creates directories, checks tools, and defines Java paths.

Edit the top of 00_setup.sh to permanently set your Java paths and memory limits:
How to set custom paths:

```
# Example: Pointing to specific installations
export GATK_JAVA=/usr/lib/jvm/java-8-openjdk/bin/java
export SNPEFF_JAVA=/usr/lib/jvm/java-11-openjdk/bin/java
```
# Then run your scripts as usual
bash 09_variant_calling.sh SRR35818859


### Module 1: Fetch Data (`01_fetch_data.sh`)

Downloads raw FASTQ files from NCBI SRA (BioProject PRJNA1346681).

```bash
bash 01_fetch_data.sh SRR35818859          # One sample
bash 01_fetch_data.sh --serotype denv1     # All DENV1 (6 samples)
bash 01_fetch_data.sh --all                # All 40 samples
```

**QC**: Verifies both R1 and R2 files exist and are non-empty.

### Module 2: Sample Sheet (`02_create_samplesheet.sh`)

Pairs R1/R2 files into a tab-delimited samplesheet. Run once after downloading.

**QC**: Verifies all file pairs exist and samplesheet is non-empty.

### Module 3: QC & Trimming (`03_qc_trimming.sh`)

**Tool**: `fastp` removes low-quality bases, adapters, and short reads. `FastQC` generates visual reports.

**Key parameters** (from config YAML):
- Q20 quality threshold (1% error rate)
- 50 bp minimum read length
- Sliding window quality trimming
- Auto-detect adapters for paired-end

**QC checkpoint**:
| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Reads surviving | >80% | 70-80% | <70% |
| Q30 rate | >85% | 80-85% | <70% |

### Module 4: Mapping (`04_mapping.sh`)

**Tool**: `bwa-mem2` aligns reads to the reference genome.

**Why bwa-mem2**: Faster successor to bwa-mem with identical accuracy. Uses the FM-index algorithm for efficient alignment.

**Key parameter**: Read group tags (`@RG`) — required by GATK for sample identification.

### Module 5: SAM to BAM (`05_sam_to_bam.sh`)

**Tools**: `samtools view`, `sort`, `index`, `flagstat`

Converts text SAM → compressed BAM, filters out:
- Unmapped reads (flag 0x4)
- Secondary alignments (flag 0x100)
- Supplementary alignments (flag 0x800)
- Reads with MAPQ < 20

**QC checkpoint**:
| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Mapping rate | >80% | 50-80% | <50% |

A mapping rate below 50% usually means wrong reference or heavy host contamination.

### Module 6: Deduplication & Primer Trimming (`06_dedup_primer_trim.sh`)

**Tools**: `samtools markdup`, `ivar trim`

1. **Deduplication**: Removes PCR duplicate reads that inflate coverage artificially
2. **Primer trimming**: Removes synthetic primer sequences from read ends

**Why primer trimming matters**: Primers are not from the sample. Without trimming, mutations in primer binding sites become invisible, and primer-derived bases create false reference calls.

**QC checkpoint**:
| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Duplication rate | <50% | 50-80% | >80% |

### Module 7: snpEff Database (`07_snpeff_database.sh`)

**Run ONCE per serotype**, not per sample. Builds a local annotation database from the GenBank file.

If you've already built the database for DENV2, skip this when processing additional DENV2 samples.

### Module 8: Coverage & Consensus (`08_coverage_consensus.sh`)

**Tools**: `samtools depth`, `samtools mpileup`, `ivar consensus`

- Calculates per-position read depth
- Generates consensus sequence (majority base at each position)
- Positions below 20x depth → "N" (unknown)

**QC checkpoint**:
| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Genome at ≥20x | >90% | 70-90% | <70% |
| N content | <5% | 5-20% | >20% |

### Module 9: Variant Calling (`09_variant_calling.sh`)

**Tool**: GATK HaplotypeCaller with hard filtering

**Key parameters**:
- Ploidy = 1 (haploid virus)
- Four hard filters: QD, strand bias, mapping quality, depth
- Java heap capped at `GATK_MEMORY` (default 4g)

**If GATK appears to hang**: Your VM may not have enough memory. GATK's JVM will thrash into swap on low-memory systems. Check with `free -h` (Linux) or `sysctl hw.memsize` (macOS). To reduce memory usage:
```bash
GATK_MEMORY=2g bash 09_variant_calling.sh SRR35818859
```
Or edit `GATK_MEMORY` in `00_setup.sh`. Minimum recommended system RAM: 8 GB.

**QC checkpoint**:
| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| Variant count | 5-50 | 50-200 | >200 or 0 |

### Module 10: Annotation (`10_annotation.sh`) — THE CRITICAL MODULE

**Tool**: `annotate_from_config.py`

Maps each variant to its gene, determines effect (synonymous, missense, frameshift).

#### The Frameshift Rule

> **Dengue virus encodes a single polyprotein. A real frameshift would destroy all downstream proteins and produce a non-viable virus. If you see >1 frameshift in your annotation, your data has artifacts. Do not publish.**

**QC checkpoint**:
| Metric | PASS | WARN | FAIL |
|--------|------|------|------|
| **Frameshifts** | **0** | **1** | **>1 = ARTIFACTS** |

Compare your output against:
- `example_output/annotation_example.tsv` — what GOOD data looks like
- `example_output/annotation_with_frameshift.tsv` — what BAD data looks like

### Modules 11-13: Self-Reference (Pass 2)

**Purpose**: Separate fixed lineage mutations from intrahost variants.

- **11**: Fills Ns in consensus with reference bases, indexes as self-reference
- **12**: Re-maps reads to self-reference, calls variants at ≥3% frequency with ivar
- **13**: Compares Pass 1 vs Pass 2 to classify variants as FIXED_SHARED, FIXED_UNIQUE, MIXED, or INTRAHOST

**Screen recommended** for Module 12 (10-30 min per sample). See `lib/screen_tutorial.md`.

### Module 14: Protein Extraction (`14_protein_extraction.sh`)

Extracts per-protein FASTA files (E, NS1, NS3, NS5, etc.) from consensus sequences. **Indel-aware**: adjusts protein boundaries based on VCF indels.

---

## Intermediate Batching: Multiple Samples

After validating one sample, process a full serotype:

```bash
export SEROTYPE=denv1
bash 00_setup.sh

# Download all DENV1 samples
bash 01_fetch_data.sh --serotype denv1
bash 02_create_samplesheet.sh

# Process each sample through the pipeline
for SAMPLE in SRR35818875 SRR35818879 SRR35818882 SRR35818886 SRR35818888 SRR35818898; do
    echo "=== Processing ${SAMPLE} ==="
    bash 03_qc_trimming.sh "$SAMPLE"
    bash 04_mapping.sh "$SAMPLE"
    bash 05_sam_to_bam.sh "$SAMPLE"
    bash 06_dedup_primer_trim.sh "$SAMPLE"
    bash 08_coverage_consensus.sh "$SAMPLE"
    bash 09_variant_calling.sh "$SAMPLE"
    bash 10_annotation.sh "$SAMPLE"
    bash 11_self_reference_prep.sh "$SAMPLE"
    bash 12_self_reference_calling.sh "$SAMPLE"
    bash 13_variant_comparison.sh "$SAMPLE"
done

# Build snpEff database once (before or during the loop)
bash 07_snpeff_database.sh

# Extract proteins for all samples at once
bash 14_protein_extraction.sh
```

For parallel processing on a multi-core machine, use `screen` sessions (see `lib/screen_tutorial.md`).

---

## QC Reference Card

| Checkpoint | Script | Metric | PASS | WARN | FAIL |
|-----------|--------|--------|------|------|------|
| Read QC | 03 | Survival rate | >80% | 70-80% | <70% |
| Read QC | 03 | Q30 rate | >85% | 80-85% | <70% |
| Mapping | 05 | Mapping rate | >80% | 50-80% | <50% |
| Dedup | 06 | Duplication | <50% | 50-80% | >80% |
| Coverage | 08 | Genome ≥20x | >90% | 70-90% | <70% |
| Coverage | 08 | N content | <5% | 5-20% | >20% |
| Variants | 09 | Count | 5-50 | 50-200 | >200 or 0 |
| **Annotation** | **10** | **Frameshifts** | **0** | **1** | **>1** |

---

## Troubleshooting

### "command not found" for bioinformatics tools
```bash
# Activate the conda environment
conda activate dengue_pipeline
# Or on HTCF:
source /ref/sahlab/software/anaconda3/bin/activate
conda activate dengue_pipeline
```

### "No FASTQ files found"
Check that files are in the `fastq/` directory with the naming pattern `{SAMPLE}_1.fastq.gz` and `{SAMPLE}_2.fastq.gz`.

### Very low mapping rate (<10%)
You likely have the wrong reference genome. BLAST a few reads against NCBI to confirm the correct serotype.

### 0 variants called
- Check that the BAM has alignments: `samtools view -c sample_analysis.bam`
- Check coverage: is the genome actually covered?
- Verify the reference matches the sample serotype

### GATK errors about read groups
Ensure you used the `-R` flag in `bwa-mem2 mem` (Module 4 does this automatically).

### "Multiple frameshifts detected"
**This is the key lesson of this tutorial.** Do NOT proceed. Investigate:
1. Look at the indels in the VCF — are they near primer binding sites?
2. Check coverage at indel positions — low coverage = unreliable
3. View the region in IGV — do the reads actually support the indel?
4. Most likely cause: sequencing errors in homopolymer regions amplified by low stringency filtering

### Out of disk space
- SAM files are the largest intermediates — scripts 05 deletes them after conversion
- Remove intermediate BAMs after each sample completes
- Use `du -sh output/` to monitor space usage

---

## Adapting for Other Viruses

See [15_adapt_to_new_virus.md](15_adapt_to_new_virus.md) for a step-by-step guide to using this pipeline with non-DENV viruses.
