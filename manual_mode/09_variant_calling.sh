#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 9: Variant Calling
# ═══════════════════════════════════════════════════════════════
# TOOL:    GATK HaplotypeCaller, GATK VariantFiltration, GATK SelectVariants
# INPUT:   Analysis BAM, reference FASTA
# OUTPUT:  Raw VCF, filtered VCF (PASS variants only)
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

# ─── Verification ─────────────────────────────────────────────
# Ensure we have our specific Java path from setup
if [ -z "${GATK_JAVA:-}" ]; then
    print_fail "GATK_JAVA path not set. Please run 00_setup.sh first."
    exit 1
fi

ANALYSIS_BAM="${OUTPUT_DIR}/${SAMPLE}_analysis.bam"
if [ ! -f "$ANALYSIS_BAM" ]; then
    print_fail "Analysis BAM not found: $ANALYSIS_BAM"
    exit 1
fi

# Locate the GATK jar file from the conda environment
GATK_JAR=$(which gatk)

# ─── Output files ─────────────────────────────────────────────
RAW_VCF="${OUTPUT_DIR}/${SAMPLE}_raw.vcf"
FILTERED_VCF="${OUTPUT_DIR}/${SAMPLE}_filtered.vcf"
PASS_VCF="${OUTPUT_DIR}/${SAMPLE}_pass.vcf"

# ─── STEP 1: GATK HaplotypeCaller ─────────────────────────────
echo -e "${BLUE}[1/3]${RESET} Running GATK HaplotypeCaller..."
echo "  Using Java: $GATK_JAVA"

GATK_MEMORY="${GATK_MEMORY:-4g}"
GATK_LOG="${OUTPUT_DIR}/${SAMPLE}_haplotypecaller.log"

# Updated execution: Using specific Java path + the gatk wrapper/jar
"$GATK_JAVA" "-Xmx${GATK_MEMORY}" -jar "$GATK_JAR" HaplotypeCaller \
    -R "$REFERENCE_FASTA" \
    -I "$ANALYSIS_BAM" \
    -O "$RAW_VCF" \
    -ploidy 1 \
    --standard-min-confidence-threshold-for-calling 30 \
    --min-base-quality-score 20 \
    2>"$GATK_LOG" || { cat "$GATK_LOG" >&2; exit 1; }

grep -E "^(INFO|WARN)" "$GATK_LOG" | tail -5 || true

# ─── STEP 2: Apply variant filters ────────────────────────────
echo -e "${BLUE}[2/3]${RESET} Applying quality filters..."
"$GATK_JAVA" "-Xmx${GATK_MEMORY}" -jar "$GATK_JAR" VariantFiltration \
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
"$GATK_JAVA" "-Xmx${GATK_MEMORY}" -jar "$GATK_JAR" SelectVariants \
    -R "$REFERENCE_FASTA" \
    -V "$FILTERED_VCF" \
    -O "$PASS_VCF" \
    --exclude-filtered \
    2>&1 | tail -2

# ─── Summary ──────────────────────────────────────────────────
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_variants "$PASS_VCF"

echo ""
echo -e "  ${BOLD}Variant summary for ${SAMPLE}:${RESET}"
RAW_COUNT=$(grep -cv "^#" "$RAW_VCF" || echo 0)
PASS_COUNT=$(grep -cv "^#" "$PASS_VCF" || echo 0)
SNP_COUNT=$(grep -v "^#" "$PASS_VCF" | awk 'length($4)==1 && length($5)==1' | wc -l | tr -d ' ')
INDEL_COUNT=$((PASS_COUNT - SNP_COUNT))

echo "    Raw calls:       ${RAW_COUNT}"
echo "    PASS variants:   ${PASS_COUNT}"
echo "    SNPs:            ${SNP_COUNT}"
echo "    Indels:          ${INDEL_COUNT}"

if [ "$INDEL_COUNT" -gt 2 ]; then
    print_warn "High indel count (${INDEL_COUNT}). Possible sequencing artifacts."
fi

echo ""
echo -e "  Next step:  ${CYAN}bash 10_annotation.sh ${SAMPLE}${RESET}"
