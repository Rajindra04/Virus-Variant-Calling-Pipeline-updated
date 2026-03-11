#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 4: Read Mapping
# ═══════════════════════════════════════════════════════════════
# TOOL:    bwa-mem2
# INPUT:   Trimmed paired FASTQs, indexed reference FASTA
# OUTPUT:  SAM file (aligned reads)
# TIME:    ~5-10 minutes per sample
# SCREEN:  Recommended
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Aligns each sequencing read to the reference genome,
#   finding the best position where it matches. The output
#   is a SAM file — a text record of where each read maps.
#
# WHY IT MATTERS:
#   Mapping is the foundation of all variant calling. Reads
#   mapped to the wrong position produce false variants.
#   bwa-mem2 is the successor to bwa-mem — faster, same
#   accuracy. The reference must be indexed first (done in
#   00_setup.sh).
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 4: Read Mapping — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input files ──────────────────────────────────────────────
TRIM_R1="${OUTPUT_DIR}/${SAMPLE}_trimmed_1.fastq.gz"
TRIM_R2="${OUTPUT_DIR}/${SAMPLE}_trimmed_2.fastq.gz"

for f in "$TRIM_R1" "$TRIM_R2"; do
    if [ ! -f "$f" ]; then
        print_fail "Trimmed FASTQ not found: $f"
        echo "  Run 03_qc_trimming.sh ${SAMPLE} first."
        exit 1
    fi
done

# ─── Output ───────────────────────────────────────────────────
mkdir -p "${SAM_DIR}"
SAM_FILE="${SAM_DIR}/${SAMPLE}.sam"

# ─── STEP 1: bwa-mem2 mem ─────────────────────────────────────
# PARAMETERS:
#   -t THREADS
#       WHAT: Number of parallel threads for alignment.
#       DEFAULT: 4. Increase to 8-16 on nodes with more cores.
#
#   -R "@RG\tID:...\tSM:...\tPL:ILLUMINA"
#       WHAT: Read group tag embedded in every alignment.
#       WHY: GATK requires read groups. SM (sample) tag groups
#            reads from the same biological sample. PL (platform)
#            tells GATK what error model to expect.
#       HOW TO CHANGE: Always set SM to the sample name.

echo -e "${BLUE}[1/1]${RESET} Mapping reads with bwa-mem2..."
echo "  Reference: $(basename "$REFERENCE_FASTA")"
echo "  R1: $(basename "$TRIM_R1")"
echo "  R2: $(basename "$TRIM_R2")"
echo ""

bwa-mem2 mem \
    -t "${THREADS}" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "${REFERENCE_FASTA}" \
    "${TRIM_R1}" \
    "${TRIM_R2}" \
    > "${SAM_FILE}" \
    2> "${OUTPUT_DIR}/${SAMPLE}_bwa.log"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: SAM file exists and is non-empty               ║
# ║  ✗ FAIL if: empty or missing SAM                           ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
if [ -s "$SAM_FILE" ]; then
    SAM_SIZE=$(du -h "$SAM_FILE" | cut -f1)
    ALIGN_COUNT=$(grep -cv "^@" "$SAM_FILE" || echo 0)
    print_pass "SAM file created: ${SAM_SIZE} (${ALIGN_COUNT} alignment records)"
else
    print_fail "SAM file is empty or missing"
    echo "  Check log: ${OUTPUT_DIR}/${SAMPLE}_bwa.log"
    exit 1
fi

echo ""
echo -e "  Next step:  ${CYAN}bash 05_sam_to_bam.sh ${SAMPLE}${RESET}"
echo ""
