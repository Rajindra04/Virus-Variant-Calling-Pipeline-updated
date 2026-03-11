#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 13: Variant Comparison (Pass 1 vs Pass 2)
# ═══════════════════════════════════════════════════════════════
# TOOL:    compare_variants.py
# INPUT:   Pass 1 annotations TSV, Pass 2 ivar variants TSV
# OUTPUT:  Variant classification TSV (FIXED/INTRAHOST)
# TIME:    <1 minute per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Compares variants found in Pass 1 (vs reference) with
#   variants found in Pass 2 (vs self-reference) to classify
#   each variant as:
#
#   FIXED_SHARED:  Present at ≥95% frequency AND found in other
#                  samples of the same serotype. These are
#                  lineage-defining mutations.
#
#   FIXED_UNIQUE:  Present at ≥95% frequency but unique to this
#                  sample. Could be real or could be artifacts.
#
#   MIXED:         Present at high frequency in Pass 1 but
#                  confirmed as mixed in Pass 2.
#
#   INTRAHOST:     Low-frequency variants — true minority
#                  populations within the host.
#
# WHY IT MATTERS:
#   Without this classification, you cannot distinguish between
#   mutations that define the outbreak lineage (shared by all
#   samples) and mutations unique to individual patients
#   (potential intrahost evolution or transmission bottleneck
#   effects).
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 13: Variant Comparison — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
ANNOTATION_TSV="${OUTPUT_DIR}/${SAMPLE}_annotations.tsv"
PASS2_VARIANTS="${PASS2_DIR}/${SAMPLE}/${SAMPLE}_ivar_variants.tsv"

if [ ! -f "$ANNOTATION_TSV" ]; then
    print_fail "Pass 1 annotations not found: $ANNOTATION_TSV"
    echo "  Run 10_annotation.sh ${SAMPLE} first."
    exit 1
fi
if [ ! -f "$PASS2_VARIANTS" ]; then
    print_fail "Pass 2 variants not found: $PASS2_VARIANTS"
    echo "  Run 12_self_reference_calling.sh ${SAMPLE} first."
    exit 1
fi

# ─── Create sample list for compare_variants.py ───────────────
# The comparison script expects a TSV with SAMPLE and SEROTYPE columns.
# For single-sample mode, we create a minimal file.
SAMPLE_LIST="${OUTPUT_DIR}/${SAMPLE}_sample_list.tsv"
echo -e "SAMPLE\tSEROTYPE" > "$SAMPLE_LIST"
echo -e "${SAMPLE}\t${SEROTYPE}" >> "$SAMPLE_LIST"

# ─── Output ───────────────────────────────────────────────────
COMPARISON_DIR="${OUTPUT_DIR}/variant_comparison"
mkdir -p "$COMPARISON_DIR"

# ─── STEP 1: Run variant comparison ───────────────────────────
# PARAMETERS:
#   --pass1_dir: Directory containing *_annotations.tsv files
#   --pass2_dir: Directory containing sample subdirectories
#                with *_ivar_variants.tsv
#   --sample_list: TSV mapping samples to serotypes
#   --output_dir: Where to write classification results

echo -e "${BLUE}[1/1]${RESET} Comparing Pass 1 vs Pass 2 variants..."
python3 "${PIPELINE_ROOT}/virus_pipeline/compare_variants.py" \
    --pass1_dir "$OUTPUT_DIR" \
    --pass2_dir "$PASS2_DIR" \
    --sample_list "$SAMPLE_LIST" \
    --output_dir "$COMPARISON_DIR"

CLASSIFICATION="${COMPARISON_DIR}/${SAMPLE}_variant_classification.tsv"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: classification file exists with variants        ║
# ║  ✗ FAIL if: no output produced                             ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
print_header "Variant Classification"

if [ -f "$CLASSIFICATION" ]; then
    TOTAL=$(tail -n+2 "$CLASSIFICATION" | wc -l | tr -d ' ')
    print_pass "Classification complete: ${TOTAL} variants classified"

    # Show classification breakdown
    echo ""
    echo -e "  ${BOLD}Classification summary:${RESET}"
    awk -F'\t' 'NR>1 {print $NF}' "$CLASSIFICATION" | sort | uniq -c | sort -rn | sed 's/^/    /'

    echo ""
    echo -e "  ${BOLD}What these mean:${RESET}"
    echo "    FIXED_SHARED  — Lineage-defining: all samples in the outbreak share this"
    echo "    FIXED_UNIQUE  — Sample-specific: only this patient has this mutation"
    echo "    MIXED         — High-frequency but not fixed (mixed infection?)"
    echo "    INTRAHOST     — Minority variant: present in <95% of viruses in this sample"
else
    print_fail "Classification file not created"
    exit 1
fi

# Clean up temp sample list
rm -f "$SAMPLE_LIST"

echo ""
echo -e "  Classification: ${CLASSIFICATION}"
echo ""
echo -e "  Next step:  ${CYAN}bash 14_protein_extraction.sh ${SAMPLE}${RESET}"
echo ""
