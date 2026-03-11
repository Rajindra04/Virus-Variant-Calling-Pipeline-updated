#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 6: Deduplication & Primer Trimming
# ═══════════════════════════════════════════════════════════════
# TOOL:    samtools markdup, ivar trim
# INPUT:   Sorted BAM from previous step
# OUTPUT:  Deduplicated, primer-trimmed BAM
# TIME:    ~2-5 minutes per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   1. Marks and removes PCR duplicate reads (same fragment
#      sequenced multiple times due to PCR amplification).
#   2. Trims primer sequences from aligned read ends using
#      the primer BED file.
#
# WHY IT MATTERS:
#   PCR duplicates inflate coverage artificially and bias allele
#   frequencies. Primer sequences are NOT from the sample — they
#   are synthetic oligonucleotides. If not trimmed, mutations in
#   primer binding sites become invisible, and primer-derived
#   bases create false reference calls.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 6: Dedup & Primer Trim — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE}.sorted.bam"
if [ ! -f "$SORTED_BAM" ]; then
    print_fail "Sorted BAM not found: $SORTED_BAM"
    echo "  Run 05_sam_to_bam.sh ${SAMPLE} first."
    exit 1
fi

# ─── Output files ─────────────────────────────────────────────
NAMESORT_BAM="${OUTPUT_DIR}/${SAMPLE}_namesort.bam"
FIXMATE_BAM="${OUTPUT_DIR}/${SAMPLE}_fixmate.bam"
POSORT_BAM="${OUTPUT_DIR}/${SAMPLE}_posort.bam"
DEDUP_BAM="${OUTPUT_DIR}/${SAMPLE}_dedup.bam"
DEDUP_STATS="${OUTPUT_DIR}/${SAMPLE}_dedup_stats.txt"
ANALYSIS_BAM="${OUTPUT_DIR}/${SAMPLE}_analysis.bam"

# ─── STEP 1: Prepare for deduplication ────────────────────────
# samtools markdup needs mate information from fixmate.
# Workflow: name-sort → fixmate → position-sort → markdup

echo -e "${BLUE}[1/4]${RESET} Name-sorting for fixmate..."
samtools sort -n -@ "${THREADS}" -o "$NAMESORT_BAM" "$SORTED_BAM"

echo -e "${BLUE}[2/4]${RESET} Adding mate score tags..."
samtools fixmate -m "$NAMESORT_BAM" "$FIXMATE_BAM"

echo -e "${BLUE}[3/4]${RESET} Position-sorting..."
samtools sort -@ "${THREADS}" -o "$POSORT_BAM" "$FIXMATE_BAM"

# ─── STEP 2: Mark and remove duplicates ───────────────────────
# PARAMETERS:
#   -r
#       WHAT: Remove duplicates (don't just mark them).
#       DEFAULT: Remove — we want clean data for variant calling.
#       HOW TO CHANGE: Drop -r to keep duplicates marked but present.
#
#   -f STATS_FILE
#       WHAT: Write deduplication statistics to a file.

echo -e "${BLUE}[4/4]${RESET} Removing PCR duplicates..."
samtools markdup \
    -r \
    -f "$DEDUP_STATS" \
    "$POSORT_BAM" \
    "$DEDUP_BAM"

samtools index "$DEDUP_BAM"

# ─── STEP 3: Primer trimming with ivar ────────────────────────
# PARAMETERS:
#   -b PRIMER_BED
#       WHAT: BED file with primer coordinates.
#       FORMAT: chrom start end name pool strand
#       WHY: Tells ivar exactly which bases are primer-derived.
#
#   -q 20
#       WHAT: Minimum quality for trimmed reads.
#       DEFAULT: Q20 matches our fastp threshold.
#
#   -m 50
#       WHAT: Minimum length after trimming.
#       DEFAULT: 50 bp — same as fastp. Reads shorter than this
#              after primer removal are discarded.
#
#   -s 4
#       WHAT: Sliding window size for quality trimming.
#       DEFAULT: 4 — standard setting.
#
#   -e
#       WHAT: Include reads that don't overlap any primer.
#       WHY: Some valid reads may fall between primer pairs.
#            Dropping them loses real data.

echo ""
echo -e "${BLUE}  Trimming primers...${RESET}"
echo "  Primer BED: $(basename "$PRIMER_BED")"

ivar trim \
    -i "$DEDUP_BAM" \
    -b "$PRIMER_BED" \
    -p "${OUTPUT_DIR}/${SAMPLE}_trimmed" \
    -q 20 \
    -m 50 \
    -s 4 \
    -e \
    2>&1 | tail -3

# Sort and index the primer-trimmed BAM
samtools sort -@ "${THREADS}" -o "$ANALYSIS_BAM" "${OUTPUT_DIR}/${SAMPLE}_trimmed.bam"
samtools index "$ANALYSIS_BAM"

# Clean up intermediate files
rm -f "$NAMESORT_BAM" "$FIXMATE_BAM" "$POSORT_BAM" "$DEDUP_BAM" "${DEDUP_BAM}.bai"
rm -f "${OUTPUT_DIR}/${SAMPLE}_trimmed.bam"
print_info "Removed intermediate BAM files"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: duplication rate <50%                           ║
# ║  ⚠ WARN if: 50-80%                                         ║
# ║  ✗ FAIL if: >80% → library complexity very low             ║
# ╚══════════════════════════════════════════════════════════════╝
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_dedup "$DEDUP_STATS"

echo ""
echo -e "  Analysis BAM: ${ANALYSIS_BAM}"
echo -e "  Dedup stats:  ${DEDUP_STATS}"
echo ""
echo -e "  Next step:  ${CYAN}bash 07_snpeff_database.sh${RESET}  (run ONCE per serotype)"
echo ""
