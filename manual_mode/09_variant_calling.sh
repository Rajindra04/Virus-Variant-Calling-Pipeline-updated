#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 9: Variant Calling
# ═══════════════════════════════════════════════════════════════
# TOOL:    GATK HaplotypeCaller, GATK VariantFiltration, GATK SelectVariants
# INPUT:   Analysis BAM, reference FASTA
# OUTPUT:  Raw VCF, filtered VCF (PASS variants only)
# TIME:    ~5-15 minutes per sample
# SCREEN:  Recommended
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Identifies positions where this sample differs from the
#   reference genome. GATK examines the pileup of reads at
#   each position and uses a probabilistic model to distinguish
#   real variants from sequencing errors.
#
# WHY IT MATTERS:
#   This is the core scientific output — the list of mutations.
#   False positives (calling errors as mutations) and false
#   negatives (missing real mutations) both compromise results.
#   The filtering step removes low-confidence calls.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 9: Variant Calling — ${SAMPLE}${RESET}"
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
RAW_VCF="${OUTPUT_DIR}/${SAMPLE}_raw.vcf"
FILTERED_VCF="${OUTPUT_DIR}/${SAMPLE}_filtered.vcf"
PASS_VCF="${OUTPUT_DIR}/${SAMPLE}_pass.vcf"

# ─── STEP 1: GATK HaplotypeCaller ─────────────────────────────
# PARAMETERS:
#   -ploidy 1
#       WHAT: Tells GATK this is a haploid genome.
#       WHY: Viruses have one copy of each gene (no diploid
#            heterozygosity). Without this, GATK would model
#            heterozygous calls that don't exist in viruses.
#       HOW TO CHANGE: Never change for single-stranded RNA viruses.
#
#   --standard-min-confidence-threshold-for-calling 30
#       WHAT: Minimum phred-scaled confidence to emit a variant.
#       DEFAULT: 30 means 99.9% confidence the variant is real.
#       HOW TO CHANGE: Lower to 20 for exploratory analysis;
#                      raise to 50 for very conservative calling.
#
#   --min-base-quality-score 20
#       WHAT: Ignore bases below Q20 when evaluating variants.
#       DEFAULT: Q20 — consistent with all other steps.

echo -e "${BLUE}[1/3]${RESET} Running GATK HaplotypeCaller..."
echo "  This may take 5-15 minutes..."
echo ""

GATK_MEMORY="${GATK_MEMORY:-4g}"
GATK_LOG="${OUTPUT_DIR}/${SAMPLE}_haplotypecaller.log"
gatk --java-options "-Xmx${GATK_MEMORY}" HaplotypeCaller \
    -R "$REFERENCE_FASTA" \
    -I "$ANALYSIS_BAM" \
    -O "$RAW_VCF" \
    -ploidy 1 \
    --standard-min-confidence-threshold-for-calling 30 \
    --min-base-quality-score 20 \
    2>"$GATK_LOG" || { cat "$GATK_LOG" >&2; exit 1; }
grep -E "^(INFO|WARN)" "$GATK_LOG" | tail -5 || true

# ─── STEP 2: Apply variant filters ────────────────────────────
# PARAMETERS (hard filters from config YAML):
#   QD < 2.0
#       WHAT: Quality by Depth. Low QD = low-confidence variant
#             inflated by high depth. Removes systematic errors.
#
#   FS > 60.0
#       WHAT: Fisher Strand bias. High FS = variant only seen on
#             one strand. Real variants appear on both strands.
#
#   MQ < 40.0
#       WHAT: Root Mean Square Mapping Quality. Low MQ = reads
#             mapping ambiguously. Unreliable variant position.
#
#   DP < 20
#       WHAT: Total depth below 20x. Not enough evidence to
#             confidently call a variant.

echo -e "${BLUE}[2/3]${RESET} Applying quality filters..."
gatk --java-options "-Xmx${GATK_MEMORY}" VariantFiltration \
    -R "$REFERENCE_FASTA" \
    -V "$RAW_VCF" \
    -O "$FILTERED_VCF" \
    --filter-expression "QD < 2.0" --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "StrandBias" \
    --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
    --filter-expression "DP < 20" --filter-name "LowDepth" \
    2>&1 | tail -2

# ─── STEP 3: Select only PASS variants ────────────────────────
echo -e "${BLUE}[3/3]${RESET} Selecting PASS variants..."
gatk --java-options "-Xmx${GATK_MEMORY}" SelectVariants \
    -R "$REFERENCE_FASTA" \
    -V "$FILTERED_VCF" \
    -O "$PASS_VCF" \
    --exclude-filtered \
    2>&1 | tail -2

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: 5-50 variants                                  ║
# ║  ⚠ WARN if: 50-200 variants                                ║
# ║  ✗ FAIL if: >200 or 0 variants → contamination?            ║
# ╚══════════════════════════════════════════════════════════════╝
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_variants "$PASS_VCF"

# Show variant summary
echo ""
echo -e "  ${BOLD}Variant summary:${RESET}"
RAW_COUNT=$(grep -cv "^#" "$RAW_VCF" || echo 0)
FILT_COUNT=$(grep -cv "^#" "$FILTERED_VCF" || echo 0)
PASS_COUNT=$(grep -cv "^#" "$PASS_VCF" || echo 0)
echo "    Raw variants:      ${RAW_COUNT}"
echo "    After filtering:   ${PASS_COUNT} PASS, $((FILT_COUNT - PASS_COUNT)) filtered out"

# Show SNP vs indel breakdown
SNP_COUNT=$(grep -v "^#" "$PASS_VCF" | awk 'length($4)==1 && length($5)==1' | wc -l | tr -d ' ')
INDEL_COUNT=$((PASS_COUNT - SNP_COUNT))
echo "    SNPs:              ${SNP_COUNT}"
echo "    Indels:            ${INDEL_COUNT}"

if [ "$INDEL_COUNT" -gt 2 ]; then
    echo ""
    print_warn "${INDEL_COUNT} indels detected — review carefully in Module 10"
    echo "  Indels in a single-polyprotein virus often indicate artifacts."
fi

echo ""
echo -e "  Raw VCF:      ${RAW_VCF}"
echo -e "  Filtered VCF: ${PASS_VCF}"
echo ""
echo -e "  Next step:  ${CYAN}bash 10_annotation.sh ${SAMPLE}${RESET}"
echo ""
