#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 0: Environment Setup (Updated for Multi-Java Support)
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

# ─── DETECT HOW WE WERE CALLED ────────────────────────────────
SETUP_SOURCE_ONLY=false
if [[ "${1:-}" == "--source-only" ]]; then
    SETUP_SOURCE_ONLY=true
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ═══════════════════════════════════════════════════════════════
# EDIT THIS BLOCK
# ═══════════════════════════════════════════════════════════════

SEROTYPE="${SEROTYPE:-denv2}"

case "${SEROTYPE}" in
    denv1)
        REFERENCE_FASTA="${PIPELINE_ROOT}/references/NC_001477.1.fasta"
        GENBANK_FILE="${PIPELINE_ROOT}/NC_001477.1.gb"
        PRIMER_BED="${PIPELINE_ROOT}/primers/denv1_primers.bed"
        CONFIG_YAML="${PIPELINE_ROOT}/configs/denv1.yaml"
        DATABASE_NAME="NC_001477.1"
        GENOME_SIZE=10735
        ;;
    denv2)
        REFERENCE_FASTA="${PIPELINE_ROOT}/references/NC_001474.2.fasta"
        GENBANK_FILE="${PIPELINE_ROOT}/NC_001474.2.gb"
        PRIMER_BED="${PIPELINE_ROOT}/primers/denv2_primers.bed"
        CONFIG_YAML="${PIPELINE_ROOT}/configs/denv2.yaml"
        DATABASE_NAME="NC_001474.2"
        GENOME_SIZE=10723
        ;;
    denv3)
        REFERENCE_FASTA="${PIPELINE_ROOT}/references/NC_001475.2.fasta"
        GENBANK_FILE="${PIPELINE_ROOT}/NC_001475.2.gb"
        PRIMER_BED="${PIPELINE_ROOT}/primers/denv3_primers.bed"
        CONFIG_YAML="${PIPELINE_ROOT}/configs/denv3.yaml"
        DATABASE_NAME="NC_001475.2"
        GENOME_SIZE=10707
        ;;
    *)
        echo "ERROR: Unknown serotype '${SEROTYPE}'."
        exit 1
        ;;
esac

THREADS="${THREADS:-4}"
GATK_MEMORY="${GATK_MEMORY:-4g}"

# Directory layout
FASTQ_DIR="${PIPELINE_ROOT}/fastq"
OUTPUT_DIR="${PIPELINE_ROOT}/output/${SEROTYPE}"
SAM_DIR="${OUTPUT_DIR}/sam_files"
BAM_DIR="${OUTPUT_DIR}"
VCF_DIR="${OUTPUT_DIR}"
CONSENSUS_DIR="${OUTPUT_DIR}"
ANNOTATION_DIR="${OUTPUT_DIR}"
PASS2_DIR="${OUTPUT_DIR}/pass2"

# Conda environment names
CONDA_ENV_NAME="dengue_pipeline"
JAVA8_ENV="java8_env"   # Specifically for GATK
JAVA11_ENV="java11_env" # Specifically for SnpEff

# ═══════════════════════════════════════════════════════════════
# Export all variables
# ═══════════════════════════════════════════════════════════════
export SEROTYPE PIPELINE_ROOT SCRIPT_DIR
export REFERENCE_FASTA GENBANK_FILE PRIMER_BED CONFIG_YAML
export DATABASE_NAME GENOME_SIZE THREADS GATK_MEMORY
export FASTQ_DIR OUTPUT_DIR SAM_DIR BAM_DIR VCF_DIR
export CONSENSUS_DIR ANNOTATION_DIR PASS2_DIR

# ─── If sourced with --source-only, stop here ─────────────────
if [ "$SETUP_SOURCE_ONLY" = true ]; then
    # We must still ensure the Java paths are exported when sourced
    export GATK_JAVA="$(conda info --base)/envs/${JAVA8_ENV}/bin/java"
    export SNPEFF_JAVA="$(conda info --base)/envs/${JAVA11_ENV}/bin/java"
    return 0 2>/dev/null || exit 0
fi

source "${SCRIPT_DIR}/lib/colors.sh"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Manual Mode Setup — ${SEROTYPE^^}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"

# ─── STEP 1: Conda and Java Environments ───────────────────────
echo -e "${BLUE}[1/4]${RESET} Setting up Conda and Java environments..."

# Create Main Env
if ! conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo -e "  Creating main environment..."
    conda env create -f "${PIPELINE_ROOT}/environment.yml"
fi

# Create Specific Java 8 Env (for GATK)
if ! conda env list | grep -q "^${JAVA8_ENV} "; then
    echo -e "  Installing Java 8 for GATK..."
    conda create -n "$JAVA8_ENV" openjdk=8 -y
fi

# Create Specific Java 11 Env (for SnpEff)
if ! conda env list | grep -q "^${JAVA11_ENV} "; then
    echo -e "  Installing Java 11 for SnpEff..."
    conda create -n "$JAVA11_ENV" openjdk=11 -y
fi

# Set Java Paths
CONDA_BASE=$(conda info --base)
export GATK_JAVA="${CONDA_BASE}/envs/${JAVA8_ENV}/bin/java"
export SNPEFF_JAVA="${CONDA_BASE}/envs/${JAVA11_ENV}/bin/java"

# Activate Main Env
source activate "${CONDA_ENV_NAME}"

print_pass "Main environment active"
print_pass "GATK Java (v8) path: $GATK_JAVA"
print_pass "SnpEff Java (v11) path: $SNPEFF_JAVA"

# ─── STEP 2: Verify tools ─────────────────────────────────────
echo ""
echo -e "${BLUE}[2/4]${RESET} Verifying tools..."

REQUIRED_TOOLS=(bwa-mem2 samtools fastp fastqc gatk ivar bcftools python3)
all_ok=true

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "$tool" &>/dev/null; then
        print_pass "$tool found"
    else
        print_fail "$tool NOT FOUND"
        all_ok=false
    fi
done

# Verify Java executables actually exist at the exported paths
if [ ! -x "$GATK_JAVA" ]; then print_fail "GATK Java not found"; all_ok=false; fi
if [ ! -x "$SNPEFF_JAVA" ]; then print_fail "SnpEff Java not found"; all_ok=false; fi

if [ "$all_ok" = false ]; then exit 1; fi

# ─── STEP 3: Verify reference files ───────────────────────────
echo ""
echo -e "${BLUE}[3/4]${RESET} Verifying reference files for ${SEROTYPE}..."

for f in "$REFERENCE_FASTA" "$GENBANK_FILE" "$PRIMER_BED" "$CONFIG_YAML"; do
    if [ -f "$f" ]; then print_pass "$(basename "$f")"; else print_fail "Missing: $f"; all_ok=false; fi
done

# Sequence Dictionary Check
if [ ! -f "${REFERENCE_FASTA%.fasta}.dict" ] && [ ! -f "${REFERENCE_FASTA}.dict" ]; then
    print_warn "Creating sequence dictionary using GATK Java..."
    # Use the specific GATK Java path here
    "$GATK_JAVA" -Xmx2g -jar $(which gatk) CreateSequenceDictionary -R "$REFERENCE_FASTA" 2>/dev/null || \
    gatk CreateSequenceDictionary -R "$REFERENCE_FASTA"
fi

# ─── STEP 4: Create directories ───────────────────────────────
echo ""
echo -e "${BLUE}[4/4]${RESET} Creating output directories..."

for dir in "$FASTQ_DIR" "$OUTPUT_DIR" "$SAM_DIR" "$PASS2_DIR"; do
    mkdir -p "$dir"
    print_pass "$(basename "$dir")/ ready"
done

echo ""
echo -e "${GREEN}${BOLD}  Setup complete!${RESET}"
echo -e "  Next step: ${CYAN}bash 01_fetch_data.sh${RESET}"
