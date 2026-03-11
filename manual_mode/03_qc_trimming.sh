#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 3: Quality Control & Read Trimming
# ═══════════════════════════════════════════════════════════════
# TOOL:    fastp, FastQC
# INPUT:   Raw paired FASTQ files
# OUTPUT:  Trimmed FASTQs, fastp JSON/HTML report, FastQC reports
# TIME:    ~2-5 minutes per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Removes low-quality bases from the ends of reads, trims
#   adapter sequences, and filters out reads that are too short.
#   FastQC then generates a visual quality report.
#
# WHY IT MATTERS:
#   Low-quality bases cause false variant calls. Adapter
#   contamination causes false insertions. Short reads map
#   ambiguously. Trimming is the most impactful QC step.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>  (e.g., SRR35818859)}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 3: QC & Trimming — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Locate input files ───────────────────────────────────────
R1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"

for f in "$R1" "$R2"; do
    if [ ! -f "$f" ]; then
        print_fail "Input file not found: $f"
        exit 1
    fi
done
print_pass "Input files found"

# ─── Output files ─────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR}"
TRIM_R1="${OUTPUT_DIR}/${SAMPLE}_trimmed_1.fastq.gz"
TRIM_R2="${OUTPUT_DIR}/${SAMPLE}_trimmed_2.fastq.gz"
FASTP_JSON="${OUTPUT_DIR}/${SAMPLE}_fastp.json"
FASTP_HTML="${OUTPUT_DIR}/${SAMPLE}_fastp.html"

# ─── STEP 1: fastp ────────────────────────────────────────────
# PARAMETERS (from config YAML):
#   --qualified_quality_phred 20
#       WHAT: Bases below Q20 (1% error rate) are "unqualified."
#       DEFAULT: Q20 is standard for amplicon sequencing.
#       HOW TO CHANGE: Lower to Q15 for low-quality runs; never below Q10.
#
#   --length_required 50
#       WHAT: Discard reads shorter than 50 bp after trimming.
#       DEFAULT: 50 bp is minimum for unique mapping to a ~10.7 kb genome.
#       HOW TO CHANGE: Can lower to 30 for very short amplicons.
#
#   --cut_front / --cut_tail
#       WHAT: Sliding window quality trimming from both ends.
#       DEFAULT: Enabled — read ends are always lowest quality.
#
#   --cut_window_size 4 / --cut_mean_quality 20
#       WHAT: 4-base sliding window; trim if mean quality drops below Q20.
#       DEFAULT: Standard Trimmomatic-equivalent settings.
#
#   --detect_adapter_for_pe
#       WHAT: Auto-detect adapter sequences for paired-end data.
#       DEFAULT: Enabled — no need to specify adapter sequences manually.
#
#   --correction
#       WHAT: Correct mismatched bases in overlapping paired-end regions.
#       DEFAULT: Enabled — improves accuracy at no cost.
#
#   --overlap_len_require 30
#       WHAT: Minimum overlap length to attempt correction.
#       DEFAULT: 30 bp prevents false corrections.

echo -e "${BLUE}[1/2]${RESET} Running fastp..."
fastp \
    --in1 "$R1" \
    --in2 "$R2" \
    --out1 "$TRIM_R1" \
    --out2 "$TRIM_R2" \
    --json "$FASTP_JSON" \
    --html "$FASTP_HTML" \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --detect_adapter_for_pe \
    --correction \
    --overlap_len_require 30 \
    --thread "${THREADS}" \
    2>&1 | tail -5

# ─── STEP 2: FastQC on trimmed reads ──────────────────────────
# PARAMETERS:
#   --outdir: Where to write HTML reports.
#   --threads: Parallel processing.
#
# FastQC reports are for VISUAL inspection — open the HTML in a browser.
# Key things to check:
#   • "Per base sequence quality" — should be green (Q>28) across all positions
#   • "Adapter Content" — should be flat (near 0%) after trimming
#   • "Sequence Length Distribution" — should peak near your amplicon size

echo -e "${BLUE}[2/2]${RESET} Running FastQC on trimmed reads..."
fastqc \
    --outdir "${OUTPUT_DIR}" \
    --threads "${THREADS}" \
    "$TRIM_R1" "$TRIM_R2" \
    2>&1 | tail -2

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: >80% reads surviving AND Q30 >85%              ║
# ║  ⚠ WARN if: 70-80% surviving OR Q30 80-85%                 ║
# ║  ✗ FAIL if: <70% surviving OR Q30 <70% → do not continue   ║
# ╚══════════════════════════════════════════════════════════════╝
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_fastp "$FASTP_JSON"

echo ""
echo -e "  Reports: ${FASTP_HTML}"
echo -e "           ${OUTPUT_DIR}/${SAMPLE}_trimmed_1_fastqc.html"
echo ""
echo -e "  Next step:  ${CYAN}bash 04_mapping.sh ${SAMPLE}${RESET}"
echo ""
