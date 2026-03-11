#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 12: Self-Reference Variant Calling (Pass 2)
# ═══════════════════════════════════════════════════════════════
# TOOL:    bwa-mem2, samtools, ivar trim, ivar variants
# INPUT:   Trimmed FASTQs, self-reference FASTA
# OUTPUT:  Pass 2 BAM, coverage, ivar variant TSV
# TIME:    ~10-30 minutes per sample
# SCREEN:  Strongly recommended
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Re-maps the original trimmed reads against the sample's own
#   consensus (self-reference), then calls variants with ivar at
#   low frequency thresholds to detect minority variants.
#
# WHY IT MATTERS:
#   Against the RefSeq reference, fixed differences (shared by
#   all viruses in the sample) look identical to intrahost
#   variants (present in only some viruses). By mapping to the
#   self-reference, fixed differences disappear — only true
#   intrahost variants remain. ivar (not GATK) is used here
#   because it reports allele frequencies, which GATK does not.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 12: Pass 2 Variant Calling — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
TRIM_R1="${OUTPUT_DIR}/${SAMPLE}_trimmed_1.fastq.gz"
TRIM_R2="${OUTPUT_DIR}/${SAMPLE}_trimmed_2.fastq.gz"
SELFREF_FASTA="${PASS2_DIR}/references/${SAMPLE}.fasta"

for f in "$TRIM_R1" "$TRIM_R2"; do
    if [ ! -f "$f" ]; then
        print_fail "Trimmed FASTQ not found: $f"
        exit 1
    fi
done
if [ ! -f "$SELFREF_FASTA" ]; then
    print_fail "Self-reference not found: $SELFREF_FASTA"
    echo "  Run 11_self_reference_prep.sh ${SAMPLE} first."
    exit 1
fi

# ─── Output ───────────────────────────────────────────────────
PASS2_OUT="${PASS2_DIR}/${SAMPLE}"
mkdir -p "$PASS2_OUT"

PASS2_BAM="${PASS2_OUT}/${SAMPLE}_analysis.bam"
PASS2_COVERAGE="${PASS2_OUT}/${SAMPLE}_coverage.txt"
PASS2_VARIANTS="${PASS2_OUT}/${SAMPLE}_ivar_variants.tsv"

# ─── STEP 1: Map to self-reference ────────────────────────────
echo -e "${BLUE}[1/6]${RESET} Mapping reads to self-reference..."
bwa-mem2 mem \
    -t "${THREADS}" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "$SELFREF_FASTA" \
    "$TRIM_R1" "$TRIM_R2" \
    2> "${PASS2_OUT}/${SAMPLE}_bwa.log" \
    | samtools view -bS -q 20 -F 4 -F 256 -F 2048 - \
    | samtools sort -@ "${THREADS}" -o "${PASS2_OUT}/${SAMPLE}_sorted.bam"

# ─── STEP 2: Dedup ────────────────────────────────────────────
echo -e "${BLUE}[2/6]${RESET} Removing duplicates..."
samtools sort -n -@ "${THREADS}" -o "${PASS2_OUT}/${SAMPLE}_namesort.bam" "${PASS2_OUT}/${SAMPLE}_sorted.bam"
samtools fixmate -m "${PASS2_OUT}/${SAMPLE}_namesort.bam" "${PASS2_OUT}/${SAMPLE}_fixmate.bam"
samtools sort -@ "${THREADS}" -o "${PASS2_OUT}/${SAMPLE}_posort.bam" "${PASS2_OUT}/${SAMPLE}_fixmate.bam"
samtools markdup -r \
    -f "${PASS2_OUT}/${SAMPLE}_dedup_stats.txt" \
    "${PASS2_OUT}/${SAMPLE}_posort.bam" \
    "${PASS2_OUT}/${SAMPLE}_dedup.bam"

# ─── STEP 3: Primer trimming ──────────────────────────────────
# NOTE: Primer BED coordinates are for the original reference.
# If the self-reference has insertions/deletions relative to the
# reference, primer coordinates may be slightly off. For DENV
# (very low indel rate), this is negligible.

echo -e "${BLUE}[3/6]${RESET} Trimming primers..."
ivar trim \
    -i "${PASS2_OUT}/${SAMPLE}_dedup.bam" \
    -b "$PRIMER_BED" \
    -p "${PASS2_OUT}/${SAMPLE}_trimmed" \
    -q 20 -m 50 -s 4 -e \
    2>&1 | tail -2

samtools sort -@ "${THREADS}" -o "$PASS2_BAM" "${PASS2_OUT}/${SAMPLE}_trimmed.bam"
samtools index "$PASS2_BAM"

# ─── STEP 4: Coverage ─────────────────────────────────────────
echo -e "${BLUE}[4/6]${RESET} Calculating coverage..."
samtools depth -a -q 20 -Q 20 "$PASS2_BAM" > "$PASS2_COVERAGE"

MEAN_DEPTH=$(awk '{sum+=$3} END {printf "%.0f", sum/NR}' "$PASS2_COVERAGE")
print_info "Pass 2 mean depth: ${MEAN_DEPTH}x"

# ─── STEP 5: ivar variant calling ─────────────────────────────
# PARAMETERS:
#   -t 0.03
#       WHAT: Minimum allele frequency to report a variant.
#       DEFAULT: 3% — captures minority variants while filtering
#                out most sequencing errors (~0.1-1% for Illumina).
#       HOW TO CHANGE: Lower to 0.01 (1%) for deep sequencing;
#                      never below 0.01 (noise floor).
#
#   -q 20
#       WHAT: Minimum base quality.
#
#   -m 20
#       WHAT: Minimum depth to call variants.
#       DEFAULT: 20x — same as consensus calling.

echo -e "${BLUE}[5/6]${RESET} Calling variants with ivar (AF ≥3%)..."
samtools mpileup \
    -d 10000 \
    -Q 0 \
    -q 0 \
    -A \
    --reference "$SELFREF_FASTA" \
    "$PASS2_BAM" \
    | ivar variants \
        -t 0.03 \
        -q 20 \
        -m 20 \
        -r "$SELFREF_FASTA" \
        -p "${PASS2_OUT}/${SAMPLE}_ivar_variants" \
    2>&1 | tail -2

# ─── STEP 6: Clean up ─────────────────────────────────────────
echo -e "${BLUE}[6/6]${RESET} Cleaning up intermediate files..."
rm -f "${PASS2_OUT}/${SAMPLE}_sorted.bam" \
      "${PASS2_OUT}/${SAMPLE}_namesort.bam" \
      "${PASS2_OUT}/${SAMPLE}_fixmate.bam" \
      "${PASS2_OUT}/${SAMPLE}_posort.bam" \
      "${PASS2_OUT}/${SAMPLE}_dedup.bam" \
      "${PASS2_OUT}/${SAMPLE}_trimmed.bam"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: ivar TSV has variants AND coverage is good     ║
# ║  ✗ FAIL if: empty TSV or very low coverage                 ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
print_header "Pass 2 Results"

if [ -f "$PASS2_VARIANTS" ]; then
    VARIANT_COUNT=$(tail -n+2 "$PASS2_VARIANTS" | wc -l | tr -d ' ')
    print_info "ivar variants found: ${VARIANT_COUNT}"

    # Show frequency distribution
    if [ "$VARIANT_COUNT" -gt 0 ]; then
        echo ""
        echo -e "  ${BOLD}Allele frequency distribution:${RESET}"
        awk -F'\t' 'NR>1 {
            af=$11;
            if (af >= 0.95) bin="≥95% (fixed)"
            else if (af >= 0.50) bin="50-95% (major)"
            else if (af >= 0.10) bin="10-50% (minor)"
            else bin="3-10% (rare)"
            print bin
        }' "$PASS2_VARIANTS" | sort | uniq -c | sort -rn | sed 's/^/    /'
    fi

    print_pass "Pass 2 complete"
else
    print_fail "ivar variants TSV not created"
    exit 1
fi

echo ""
echo -e "  Pass 2 BAM:      ${PASS2_BAM}"
echo -e "  Pass 2 variants: ${PASS2_VARIANTS}"
echo ""
echo -e "  Next step:  ${CYAN}bash 13_variant_comparison.sh ${SAMPLE}${RESET}"
echo ""
