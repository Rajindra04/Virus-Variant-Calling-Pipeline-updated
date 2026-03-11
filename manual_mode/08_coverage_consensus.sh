#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 8: Coverage Analysis & Consensus Calling
# ═══════════════════════════════════════════════════════════════
# TOOL:    samtools depth, samtools mpileup, ivar consensus
# INPUT:   Analysis BAM (deduplicated, primer-trimmed)
# OUTPUT:  Coverage file, consensus FASTA, coverage plot
# TIME:    ~2-5 minutes per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   1. Calculates per-position read depth across the genome.
#   2. Generates a consensus sequence — the "best guess" of
#      this sample's actual genome, base by base.
#
# WHY IT MATTERS:
#   Coverage tells you which regions of the genome have enough
#   data to make confident calls. Positions with <20x coverage
#   become "N" (unknown) in the consensus. The consensus is
#   used for phylogenetics and as the self-reference in Pass 2.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 8: Coverage & Consensus — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
ANALYSIS_BAM="${OUTPUT_DIR}/${SAMPLE}_analysis.bam"
if [ ! -f "$ANALYSIS_BAM" ]; then
    print_fail "Analysis BAM not found: $ANALYSIS_BAM"
    echo "  Run 06_dedup_primer_trim.sh ${SAMPLE} first."
    exit 1
fi

# ─── Output files ─────────────────────────────────────────────
COVERAGE="${OUTPUT_DIR}/${SAMPLE}_coverage.txt"
CONSENSUS="${OUTPUT_DIR}/${SAMPLE}.fa"

# ─── STEP 1: Calculate coverage ───────────────────────────────
# PARAMETERS:
#   -a
#       WHAT: Output ALL positions, including zero-depth.
#       WHY: We need to know where there is NO coverage.
#   -q 20
#       WHAT: Minimum base quality to count a base.
#       DEFAULT: Q20 — consistent with other quality thresholds.
#   -Q 20
#       WHAT: Minimum mapping quality.

echo -e "${BLUE}[1/3]${RESET} Calculating per-position coverage..."
samtools depth \
    -a \
    -q 20 \
    -Q 20 \
    "$ANALYSIS_BAM" \
    > "$COVERAGE"

# Quick coverage summary
MEAN_DEPTH=$(awk '{sum+=$3} END {printf "%.0f", sum/NR}' "$COVERAGE")
ZERO_DEPTH=$(awk '$3==0 {count++} END {print count+0}' "$COVERAGE")
print_info "Mean depth: ${MEAN_DEPTH}x"
print_info "Zero-depth positions: ${ZERO_DEPTH}"

# ─── STEP 2: Generate consensus with ivar ─────────────────────
# PARAMETERS:
#   samtools mpileup:
#     -d 10000
#         WHAT: Maximum read depth to consider per position.
#         DEFAULT: 10000 — high enough to not truncate amplicon data.
#     -Q 0
#         WHAT: Minimum base quality for mpileup. Set to 0 here
#               because ivar does its own quality filtering.
#     -q 0
#         WHAT: Minimum mapping quality for mpileup. Same reason.
#     -A
#         WHAT: Include anomalous read pairs.
#
#   ivar consensus:
#     -q 20
#         WHAT: Minimum quality score to count a base.
#         DEFAULT: Q20 — 1% error rate threshold.
#     -t 0.5
#         WHAT: Minimum frequency to call a consensus base.
#         DEFAULT: 0.5 (majority rule). A base needs >50% of reads
#                  to be called. Below 50%, it becomes ambiguous.
#         HOW TO CHANGE: Lower to 0.3 for mixed infections (NOT
#                        recommended for standard analysis).
#     -m 20
#         WHAT: Minimum depth to make a call. Below this → "N".
#         DEFAULT: 20x is standard for confident base calling.
#         HOW TO CHANGE: Lower to 10x for low-coverage samples,
#                        but results are less reliable.
#     -n N
#         WHAT: Character for positions below minimum depth.
#         DEFAULT: "N" — standard ambiguity code.

echo -e "${BLUE}[2/3]${RESET} Generating consensus sequence..."
samtools mpileup \
    -d 10000 \
    -Q 0 \
    -q 0 \
    -A \
    --reference "$REFERENCE_FASTA" \
    "$ANALYSIS_BAM" \
    | ivar consensus \
        -q 20 \
        -t 0.5 \
        -m 20 \
        -n N \
        -p "${OUTPUT_DIR}/${SAMPLE}" \
    2>&1 | tail -3

# Rename ivar output to .fa
if [ -f "${OUTPUT_DIR}/${SAMPLE}.fa" ]; then
    print_pass "Consensus FASTA created"
elif [ -f "${OUTPUT_DIR}/${SAMPLE}.consensus.fa" ]; then
    mv "${OUTPUT_DIR}/${SAMPLE}.consensus.fa" "${OUTPUT_DIR}/${SAMPLE}.fa"
    print_pass "Consensus FASTA created"
fi

# ─── STEP 3: Generate coverage plot ───────────────────────────
echo -e "${BLUE}[3/3]${RESET} Generating coverage plot..."
python3 -c "
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

positions, depths = [], []
with open('${COVERAGE}') as f:
    for line in f:
        cols = line.strip().split('\t')
        positions.append(int(cols[1]))
        depths.append(int(cols[2]))

fig, ax = plt.subplots(figsize=(12, 5))
ax.fill_between(positions, depths, alpha=0.4, color='steelblue')
ax.plot(positions, depths, linewidth=0.5, color='steelblue')
ax.axhline(y=20, color='red', linestyle='--', linewidth=1, label='20X threshold')
ax.set_yscale('log')
ax.set_ylim(bottom=0.5)
ax.set_xlabel('Genome Position (bp)')
ax.set_ylabel('Read Depth (log scale)')
ax.set_title('Coverage: ${SAMPLE}')
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)
fig.savefig('${OUTPUT_DIR}/${SAMPLE}_coverage.png', dpi=150, bbox_inches='tight')
plt.close()
print('  Plot saved.')
"
print_pass "Coverage plot: ${OUTPUT_DIR}/${SAMPLE}_coverage.png"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: ≥90% genome at ≥20x AND <5% N content         ║
# ║  ⚠ WARN if: 70-90% at ≥20x OR 5-20% N                     ║
# ║  ✗ FAIL if: <70% at ≥20x OR >20% N → do not continue      ║
# ╚══════════════════════════════════════════════════════════════╝
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_coverage "$COVERAGE" "${OUTPUT_DIR}/${SAMPLE}.fa" "$GENOME_SIZE"

echo ""
echo -e "  Consensus:     ${OUTPUT_DIR}/${SAMPLE}.fa"
echo -e "  Coverage:      ${COVERAGE}"
echo -e "  Coverage plot: ${OUTPUT_DIR}/${SAMPLE}_coverage.png"
echo ""
echo -e "  Next step:  ${CYAN}bash 09_variant_calling.sh ${SAMPLE}${RESET}"
echo ""
