#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 2: Create Sample Sheet
# ═══════════════════════════════════════════════════════════════
# TOOL:    python (create_samplesheet.py)
# INPUT:   FASTQ files in fastq/ directory
# OUTPUT:  samplesheet.tsv with columns: sample_name, read1, read2
# TIME:    <1 minute
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Scans the fastq/ directory and matches R1/R2 pairs by name.
#   Produces a tab-delimited file that tells every downstream
#   step which files belong to which sample.
#
# WHY IT MATTERS:
#   Mismatched pairs (R1 from sample A with R2 from sample B)
#   would produce nonsense alignments. The samplesheet is your
#   single source of truth for file pairing.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 2: Create Sample Sheet${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

SAMPLESHEET="${OUTPUT_DIR}/samplesheet.tsv"

# ─── Check FASTQ directory ────────────────────────────────────
if [ ! -d "$FASTQ_DIR" ] || [ -z "$(ls -A "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null)" ]; then
    print_fail "No .fastq.gz files found in ${FASTQ_DIR}"
    echo "  Run 01_fetch_data.sh first."
    exit 1
fi

# Count FASTQ files
N_FILES=$(ls "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
print_info "Found ${N_FILES} FASTQ files in ${FASTQ_DIR}"

# ─── STEP 1: Run create_samplesheet.py ────────────────────────
# The pipeline's built-in samplesheet creator handles SRA naming
# conventions (_1.fastq.gz / _2.fastq.gz).
#
# PARAMETERS:
#   directory
#       WHAT: Path to folder with FASTQ files.
#   sample_sheet_file
#       WHAT: Output TSV path.
#   --sample_names (optional)
#       WHAT: Override auto-detected sample names.
#       WHY: SRA accessions (SRR...) are the default names. You
#            can provide comma-separated friendlier names if needed.
#       EXAMPLE: --sample_names "Patient1,Patient2,Patient3"

mkdir -p "$(dirname "$SAMPLESHEET")"

echo "  Running create_samplesheet.py..."
python3 "${PIPELINE_ROOT}/virus_pipeline/create_samplesheet.py" \
    "${FASTQ_DIR}" \
    "${SAMPLESHEET}"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: samplesheet has >0 rows, all files exist       ║
# ║  ✗ FAIL if: empty samplesheet or missing files             ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
print_header "Sample Sheet"

if [ ! -f "$SAMPLESHEET" ]; then
    print_fail "Samplesheet was not created"
    exit 1
fi

N_SAMPLES=$(tail -n+2 "$SAMPLESHEET" | wc -l | tr -d ' ')
print_info "Samples found: ${N_SAMPLES}"

if [ "$N_SAMPLES" -eq 0 ]; then
    print_fail "Samplesheet is empty — check FASTQ naming convention"
    exit 1
fi

# Verify all files referenced in the samplesheet exist
ALL_EXIST=true
while IFS=$'\t' read -r name r1 r2; do
    if [ ! -f "$r1" ]; then
        print_fail "R1 missing: $r1"
        ALL_EXIST=false
    fi
    if [ ! -f "$r2" ]; then
        print_fail "R2 missing: $r2"
        ALL_EXIST=false
    fi
done < <(tail -n+2 "$SAMPLESHEET")

if [ "$ALL_EXIST" = true ]; then
    print_pass "All ${N_SAMPLES} sample pairs verified"
else
    print_fail "Some FASTQ files referenced in samplesheet are missing"
    exit 1
fi

# Show the samplesheet
echo ""
echo -e "${BOLD}  Sample sheet contents:${RESET}"
column -t -s$'\t' "$SAMPLESHEET" | head -20 | sed 's/^/    /'

echo ""
echo -e "  Saved to: ${SAMPLESHEET}"
echo ""
echo -e "  Next step:  ${CYAN}bash 03_qc_trimming.sh <sample_name>${RESET}"
echo -e "  Example:    ${CYAN}bash 03_qc_trimming.sh SRR35818859${RESET}"
echo ""
