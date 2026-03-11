# Adapting the Pipeline for a New Virus

This guide walks you through modifying the manual mode tutorial to work with a non-DENV virus.

## Step 1: Obtain Reference Files

You need three files for any new virus:

| File | Source | Format |
|------|--------|--------|
| Reference FASTA | NCBI RefSeq | `.fasta` with one sequence |
| GenBank annotation | NCBI Nucleotide | `.gb` with gene/CDS features |
| Primer BED | Your primer scheme | 6-column BED |

### Download from NCBI

```bash
# Example: Zika virus (NC_012532.1)
# Reference FASTA
efetch -db nucleotide -id NC_012532.1 -format fasta > references/NC_012532.1.fasta

# GenBank annotation
efetch -db nucleotide -id NC_012532.1 -format gb > NC_012532.1.gb
```

### Primer BED format

```
# chrom          start  end    name           pool  strand
NC_012532.1      10     30     ZIKV_1_LEFT    1     +
NC_012532.1      380    400    ZIKV_1_RIGHT   1     -
```

If you don't have a primer scheme, you can still run the pipeline — just skip primer trimming steps (remove the `ivar trim` calls in scripts 06 and 12).

## Step 2: Create a Config YAML

Copy an existing config and modify:

```bash
cp configs/denv2.yaml configs/zikv.yaml
```

### What to change

| Field | DENV2 value | What to set |
|-------|-------------|-------------|
| `virus_name` | "dengue virus type 2" | Your virus name |
| `expected_genome_size` | 10723 | Your genome size (bp) |
| `database_name` | "NC_001474.2" | Your RefSeq accession |
| `transcript_annotations` | DENV2 proteins | Your virus's proteins |
| `gene_coordinates` | DENV2 coordinates | Your gene/protein coordinates |

### Getting gene coordinates

Parse from the GenBank file:

```bash
# Quick extraction of CDS and mat_peptide coordinates
python3 -c "
from Bio import SeqIO
for feat in SeqIO.read('YOUR_VIRUS.gb', 'genbank').features:
    if feat.type in ('CDS', 'mat_peptide'):
        name = feat.qualifiers.get('product', ['unknown'])[0]
        print(f'{name}: {feat.location.start+1}-{feat.location.end}')
"
```

Or manually inspect the GenBank file — look for `CDS` and `mat_peptide` features.

## Step 3: Update 00_setup.sh

Add a new case to the serotype selector:

```bash
zikv)
    REFERENCE_FASTA="${PIPELINE_ROOT}/references/NC_012532.1.fasta"
    GENBANK_FILE="${PIPELINE_ROOT}/NC_012532.1.gb"
    PRIMER_BED="${PIPELINE_ROOT}/primers/zikv_primers.bed"
    CONFIG_YAML="${PIPELINE_ROOT}/configs/zikv.yaml"
    DATABASE_NAME="NC_012532.1"
    GENOME_SIZE=10794
    ;;
```

Then set: `SEROTYPE=zikv`

## Step 4: Adjust QC Thresholds

QC thresholds depend on the virus biology:

### Frameshift tolerance by virus type

| Virus type | Example | Expected frameshifts | Reason |
|------------|---------|---------------------|--------|
| **Single polyprotein** | DENV, ZIKV, WNV, YFV | **0-1** | One reading frame for all proteins |
| **Segmented** | Influenza, Bunyaviruses | **0-1 per segment** | Each segment has its own ORF(s) |
| **Multiple ORFs** | HIV, Coronaviruses | **0-1 per ORF** | Separate reading frames |
| **Large DNA viruses** | HSV, CMV | **More tolerant** | Many independent genes |

### Variant count expectations

| Genome size | Expected variants vs reference |
|-------------|-------------------------------|
| ~10 kb (DENV, ZIKV) | 5-50 |
| ~30 kb (SARS-CoV-2) | 15-100 |
| ~150 kb (HSV) | 50-500 |

### Coverage thresholds

Coverage thresholds are generally the same regardless of virus:
- **≥90% at ≥20x** is the target for confident consensus calling
- Amplicon sequencing should achieve this for most samples
- Metagenomic sequencing may need lower thresholds (e.g., 70% at ≥10x)

## Step 5: Validate with a Single Sample

Before processing a batch:

1. Run one sample through the entire pipeline (scripts 03-14)
2. Verify all QC checkpoints pass
3. Check the annotation TSV manually:
   - Are gene names correct?
   - Are coordinate ranges right?
   - Are synonymous/missense calls reasonable?
4. Spot-check a few variants in IGV
5. Only then process the full dataset

## Common Pitfalls

### Wrong reference genome
- **Symptom**: Very low mapping rate (<50%), many variants (>200)
- **Fix**: BLAST a few reads against NCBI to confirm the correct reference

### Missing primer BED
- **Symptom**: Primer-derived bases appear as reference calls
- **Fix**: Either obtain the correct primer scheme or skip primer trimming and note the limitation

### Incorrect gene coordinates
- **Symptom**: Annotations don't match known proteins, off-by-one errors
- **Fix**: Verify coordinates against the GenBank file; remember GenBank uses 1-based coordinates

### Genome with multiple segments
- **Symptom**: Only one segment gets processed
- **Fix**: Process each segment as a separate "reference" — run the pipeline once per segment
