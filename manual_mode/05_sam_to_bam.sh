#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 5: SAM to BAM Conversion & Quality Filtering
# ═══════════════════════════════════════════════════════════════
# TOOL:    samtools view, samtools sort, samtools index, samtools flagstat
# INPUT:   SAM file from mapping
# OUTPUT:  Sorted, filtered, indexed BAM + flagstat report
# TIME:    ~2-3 minutes per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Converts the text SAM to compressed binary BAM, removes
#   unmapped/low-quality/supplementary reads, sorts by genomic
#   position, and creates an index for random access.
#
# WHY IT MATTERS:
#   SAM files are huge text files — BAM is ~5x smaller. Sorting
#   and indexing is required by every downstream tool (ivar,
#   GATK, IGV). Quality filtering removes reads that mapped
#   poorly, which would otherwise create false variants.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 5: SAM → BAM — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
SAM_FILE="${SAM_DIR}/${SAMPLE}.sam"
if [ ! -f "$SAM_FILE" ]; then
    print_fail "SAM file not found: $SAM_FILE"
    echo "  Run 04_mapping.sh ${SAMPLE} first."
    exit 1
fi

# ─── Output files ─────────────────────────────────────────────
RAW_BAM="${OUTPUT_DIR}/${SAMPLE}_raw.bam"
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE}.sorted.bam"
FLAGSTAT="${OUTPUT_DIR}/${SAMPLE}.flagstat"

# ─── STEP 1: Convert SAM to BAM with quality filter ───────────
# PARAMETERS:
#   -q 20
#       WHAT: Minimum mapping quality (MAPQ). 20 means ≤1% chance
#             the read is mapped to the wrong location.
#       DEFAULT: 20 is standard for viral amplicon data.
#       HOW TO CHANGE: Lower to 10 for divergent references; never 0.
#
#   -F 4
#       WHAT: Exclude unmapped reads (SAM flag 0x4).
#       WHY: Unmapped reads have no position — useless for variants.
#
#   -F 256
#       WHAT: Exclude secondary alignments (0x100).
#       WHY: Keep only the best alignment per read.
#
#   -F 2048
#       WHAT: Exclude supplementary alignments (0x800).
#       WHY: These are chimeric read parts — artifacts in amplicon data.

echo -e "${BLUE}[1/4]${RESET} Filtering and converting SAM → BAM..."
samtools view \
    -bS \
    -q 20 \
    -F 4 -F 256 -F 2048 \
    "$SAM_FILE" \
    > "$RAW_BAM"

# ─── STEP 2: Sort by position ─────────────────────────────────
# PARAMETERS:
#   -@ THREADS
#       WHAT: Compression threads for writing BAM.
#   Output is position-sorted — required for indexing.

echo -e "${BLUE}[2/4]${RESET} Sorting BAM by position..."
samtools sort \
    -@ "${THREADS}" \
    -o "$SORTED_BAM" \
    "$RAW_BAM"

# ─── STEP 3: Index BAM ────────────────────────────────────────
# Creates .bai index for random access. Required by GATK, ivar, IGV.

echo -e "${BLUE}[3/4]${RESET} Indexing BAM..."
samtools index "$SORTED_BAM"

# ─── STEP 4: Generate flagstat ─────────────────────────────────
# Flagstat gives a one-page summary of mapping statistics.

echo -e "${BLUE}[4/4]${RESET} Generating mapping statistics..."
samtools flagstat "$SORTED_BAM" > "$FLAGSTAT"

# Clean up intermediate files
rm -f "$RAW_BAM"
rm -f "$SAM_FILE"
print_info "Removed intermediate SAM and unsorted BAM"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: mapping rate >80%                              ║
# ║  ⚠ WARN if: 50-80%                                         ║
# ║  ✗ FAIL if: <50% → wrong reference? contamination?         ║
# ╚══════════════════════════════════════════════════════════════╝
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_mapping "$FLAGSTAT"

echo ""
echo -e "  BAM file:  ${SORTED_BAM}"
echo -e "  Flagstat:  ${FLAGSTAT}"
echo ""
echo -e "  Next step:  ${CYAN}bash 06_dedup_primer_trim.sh ${SAMPLE}${RESET}"
echo ""
