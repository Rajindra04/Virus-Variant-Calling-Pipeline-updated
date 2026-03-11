#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 0: Environment Setup
# ═══════════════════════════════════════════════════════════════
# TOOL:    conda, bash
# INPUT:   None (first script to run)
# OUTPUT:  Environment variables, directory structure
# TIME:    ~1 minute (or ~10 min if creating conda env)
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Sets up all paths, creates output directories, activates the
#   conda environment, and verifies every tool is available.
#
# WHY IT MATTERS:
#   Every subsequent script sources this file. If paths are wrong
#   or tools are missing, you'll get cryptic errors mid-pipeline.
#   Catching problems here saves hours of debugging.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

# ─── DETECT HOW WE WERE CALLED ────────────────────────────────
# Other scripts call: source 00_setup.sh --source-only
# Running directly:   bash 00_setup.sh  (full setup + checks)
SETUP_SOURCE_ONLY=false
if [[ "${1:-}" == "--source-only" ]]; then
    SETUP_SOURCE_ONLY=true
fi

# ─── PIPELINE ROOT (auto-detected) ────────────────────────────
# This is the directory containing the git repo
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ═══════════════════════════════════════════════════════════════
# EDIT THIS BLOCK — adjust for your project
# ═══════════════════════════════════════════════════════════════

# Which serotype are you processing? (denv1, denv2, or denv3)
SEROTYPE="${SEROTYPE:-denv2}"

# Reference files (auto-set from serotype)
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
        echo "ERROR: Unknown serotype '${SEROTYPE}'. Use denv1, denv2, or denv3."
        exit 1
        ;;
esac

# Processing parameters
THREADS="${THREADS:-4}"

# Directory layout
FASTQ_DIR="${PIPELINE_ROOT}/fastq"
OUTPUT_DIR="${PIPELINE_ROOT}/output/${SEROTYPE}"
SAM_DIR="${OUTPUT_DIR}/sam_files"
BAM_DIR="${OUTPUT_DIR}"
VCF_DIR="${OUTPUT_DIR}"
CONSENSUS_DIR="${OUTPUT_DIR}"
ANNOTATION_DIR="${OUTPUT_DIR}"
PASS2_DIR="${OUTPUT_DIR}/pass2"

# Conda environment name (from environment.yml)
CONDA_ENV_NAME="dengue_pipeline"

# ═══════════════════════════════════════════════════════════════
# Export all variables for child scripts
# ═══════════════════════════════════════════════════════════════
export SEROTYPE PIPELINE_ROOT SCRIPT_DIR
export REFERENCE_FASTA GENBANK_FILE PRIMER_BED CONFIG_YAML
export DATABASE_NAME GENOME_SIZE THREADS
export FASTQ_DIR OUTPUT_DIR SAM_DIR BAM_DIR VCF_DIR
export CONSENSUS_DIR ANNOTATION_DIR PASS2_DIR
export CONDA_ENV_NAME

# ─── If sourced with --source-only, stop here ─────────────────
if [ "$SETUP_SOURCE_ONLY" = true ]; then
    return 0 2>/dev/null || exit 0
fi

# ═══════════════════════════════════════════════════════════════
# FULL SETUP MODE (run directly)
# ═══════════════════════════════════════════════════════════════

source "${SCRIPT_DIR}/lib/colors.sh"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Manual Mode Setup — ${SEROTYPE^^}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── STEP 1: Check conda environment ──────────────────────────
echo -e "${BLUE}[1/4]${RESET} Checking conda environment..."

if ! command -v conda &>/dev/null; then
    echo -e "${YELLOW}  conda not found in PATH. Attempting to activate...${RESET}"
    # HTCF-specific activation
    if [ -f /ref/sahlab/software/anaconda3/bin/activate ]; then
        source /ref/sahlab/software/anaconda3/bin/activate
    else
        print_fail "conda not found. Install Miniconda or Anaconda first."
        exit 1
    fi
fi

# Check if environment exists
if ! conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo -e "${YELLOW}  Environment '${CONDA_ENV_NAME}' not found. Creating it...${RESET}"
    echo -e "  This will take ~10 minutes."
    conda env create -f "${PIPELINE_ROOT}/environment.yml"
fi

# Activate
if [[ "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV_NAME}" ]]; then
    echo -e "  Activating ${CONDA_ENV_NAME}..."
    conda activate "${CONDA_ENV_NAME}" 2>/dev/null || \
        source activate "${CONDA_ENV_NAME}" 2>/dev/null || {
        print_fail "Could not activate conda env '${CONDA_ENV_NAME}'"
        exit 1
    }
fi
print_pass "Conda environment '${CONDA_ENV_NAME}' active"

# ─── STEP 2: Verify tools ─────────────────────────────────────
echo ""
echo -e "${BLUE}[2/4]${RESET} Verifying tools..."

REQUIRED_TOOLS=(bwa-mem2 samtools fastp fastqc gatk ivar bcftools python3)
OPTIONAL_TOOLS=(snpEff snpSift prefetch fasterq-dump)
all_ok=true

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "$tool" &>/dev/null; then
        print_pass "$tool found: $(command -v "$tool")"
    else
        print_fail "$tool NOT FOUND — install via conda"
        all_ok=false
    fi
done

for tool in "${OPTIONAL_TOOLS[@]}"; do
    if command -v "$tool" &>/dev/null; then
        print_pass "$tool found (optional)"
    else
        print_warn "$tool not found (optional — needed for some steps)"
    fi
done

if [ "$all_ok" = false ]; then
    echo ""
    print_fail "Missing required tools. Run: conda env create -f ${PIPELINE_ROOT}/environment.yml"
    exit 1
fi

# ─── STEP 3: Verify reference files ───────────────────────────
echo ""
echo -e "${BLUE}[3/4]${RESET} Verifying reference files for ${SEROTYPE}..."

for f in "$REFERENCE_FASTA" "$GENBANK_FILE" "$PRIMER_BED" "$CONFIG_YAML"; do
    if [ -f "$f" ]; then
        print_pass "$(basename "$f")"
    else
        print_fail "Missing: $f"
        all_ok=false
    fi
done

# Check if reference is indexed
if [ ! -f "${REFERENCE_FASTA}.fai" ]; then
    print_warn "Reference FASTA not indexed — will index now"
    samtools faidx "$REFERENCE_FASTA"
    print_pass "Index created: ${REFERENCE_FASTA}.fai"
fi

if [ ! -f "${REFERENCE_FASTA%.fasta}.dict" ] && [ ! -f "${REFERENCE_FASTA}.dict" ]; then
    print_warn "No sequence dictionary — will create now"
    gatk CreateSequenceDictionary -R "$REFERENCE_FASTA" 2>/dev/null
    print_pass "Sequence dictionary created"
fi

# Check bwa-mem2 index
if ! ls "${REFERENCE_FASTA}".* 2>/dev/null | grep -q "bwt"; then
    print_warn "bwa-mem2 index not found — will index now (takes ~1 min)"
    bwa-mem2 index "$REFERENCE_FASTA" 2>/dev/null
    print_pass "bwa-mem2 index created"
fi

# ─── STEP 4: Create directories ───────────────────────────────
echo ""
echo -e "${BLUE}[4/4]${RESET} Creating output directories..."

for dir in "$FASTQ_DIR" "$OUTPUT_DIR" "$SAM_DIR" "$PASS2_DIR"; do
    mkdir -p "$dir"
    print_pass "$(basename "$dir")/ ready"
done

# ─── Summary ──────────────────────────────────────────────────
echo ""
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${GREEN}${BOLD}  Setup complete!${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""
echo -e "  Serotype:     ${BOLD}${SEROTYPE^^}${RESET}"
echo -e "  Reference:    $(basename "$REFERENCE_FASTA")"
echo -e "  Primers:      $(basename "$PRIMER_BED")"
echo -e "  Config:       $(basename "$CONFIG_YAML")"
echo -e "  Output:       ${OUTPUT_DIR}"
echo -e "  Threads:      ${THREADS}"
echo ""
echo -e "  Next step:  ${CYAN}bash 01_fetch_data.sh${RESET}"
echo ""
