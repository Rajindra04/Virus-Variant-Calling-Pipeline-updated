#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 10: Variant Annotation — THE CRITICAL MODULE
# ═══════════════════════════════════════════════════════════════
# TOOL:    annotate_from_config.py (config-based annotation)
# INPUT:   Filtered VCF (_pass.vcf), reference FASTA, config YAML
# OUTPUT:  Annotation TSV with gene names, effects, protein changes
# TIME:    <1 minute per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Takes each variant from the VCF and determines:
#   - Which gene/protein it falls in
#   - Whether it changes the amino acid (missense vs synonymous)
#   - Whether it shifts the reading frame (frameshift)
#   - The HGVS notation (e.g., p.Ala123Val)
#
# WHY IT MATTERS:
#   This is where you discover biologically meaningful mutations.
#   More importantly, this is where you catch ARTIFACTS.
#
#   ╔══════════════════════════════════════════════════════════╗
#   ║  THE KEY LESSON:                                        ║
#   ║                                                          ║
#   ║  Dengue virus has ONE polyprotein reading frame.         ║
#   ║  A real frameshift would destroy every protein           ║
#   ║  downstream of the indel → non-viable virus.            ║
#   ║                                                          ║
#   ║  If you see >1 frameshift in your annotation:           ║
#   ║  YOUR DATA HAS ARTIFACTS. DO NOT PUBLISH.               ║
#   ╚══════════════════════════════════════════════════════════╝
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 10: Variant Annotation — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
PASS_VCF="${OUTPUT_DIR}/${SAMPLE}_pass.vcf"
if [ ! -f "$PASS_VCF" ]; then
    print_fail "Filtered VCF not found: $PASS_VCF"
    echo "  Run 09_variant_calling.sh ${SAMPLE} first."
    exit 1
fi

# ─── STEP 1: Run config-based annotation ──────────────────────
# This uses the pipeline's annotate_from_config.py, which reads
# gene coordinates from the config YAML (no snpEff database needed).
#
# WHAT IT PRODUCES (24-column TSV):
#   CHROM POS ID REF ALT QUAL FILTER Total_Depth Allele_Frequency
#   Strand_Bias Allelic_Depths EFFECT PUTATIVE_IMPACT GENE_NAME
#   GENE_ID FEATURE_TYPE FEATURE_ID TRANSCRIPT_TYPE HGVSc HGVSp
#   cDNA_POSITION_AND_LENGTH CDS_POSITION_AND_LENGTH
#   PROTEIN_POSITION_AND_LENGTH ERROR
#
# EFFECT column values:
#   synonymous_variant     — same amino acid (silent)
#   missense_variant       — different amino acid
#   frameshift_variant     — indel NOT divisible by 3 → SHIFTS READING FRAME
#   inframe_insertion      — insertion divisible by 3 → adds amino acids
#   inframe_deletion       — deletion divisible by 3 → removes amino acids
#   stop_gained            — creates premature stop codon
#   upstream_gene_variant  — before the coding region
#   downstream_gene_variant — after the coding region

echo -e "${BLUE}[1/1]${RESET} Annotating variants..."
python3 "${PIPELINE_ROOT}/virus_pipeline/annotate_from_config.py" \
    --vcf "$PASS_VCF" \
    --reference "$REFERENCE_FASTA" \
    --config "$CONFIG_YAML" \
    --sample_name "$SAMPLE" \
    --output_dir "$OUTPUT_DIR"

ANNOTATION_TSV="${OUTPUT_DIR}/${SAMPLE}_annotations.tsv"

if [ ! -f "$ANNOTATION_TSV" ]; then
    print_fail "Annotation TSV was not created"
    exit 1
fi

# ─── Show annotation summary ──────────────────────────────────
echo ""
echo -e "  ${BOLD}Annotation Summary:${RESET}"

# Count variants by effect
while IFS= read -r line; do
    echo "    $line"
done < <(awk -F'\t' 'NR>1 {print $12}' "$ANNOTATION_TSV" | sort | uniq -c | sort -rn)

# Count variants by gene
echo ""
echo -e "  ${BOLD}Variants by Gene:${RESET}"
while IFS= read -r line; do
    echo "    $line"
done < <(awk -F'\t' 'NR>1 {print $14}' "$ANNOTATION_TSV" | sort | uniq -c | sort -rn)

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT — MOST IMPORTANT IN THE ENTIRE PIPELINE      ║
# ║  ✓ PASS if: 0 frameshifts                                  ║
# ║  ⚠ WARN if: 1 frameshift (review manually)                 ║
# ║  ✗ FAIL if: >1 frameshift → ARTIFACTS — DO NOT CONTINUE    ║
# ╚══════════════════════════════════════════════════════════════╝
source "${SCRIPT_DIR}/lib/qc_check.sh"
qc_check_frameshifts "$ANNOTATION_TSV"

echo ""
echo -e "  Annotation: ${ANNOTATION_TSV}"
echo ""
echo -e "  ${BOLD}How to read the annotation file:${RESET}"
echo "    column -t -s\$'\\t' ${ANNOTATION_TSV} | less -S"
echo ""
echo -e "  ${BOLD}Compare with example outputs:${RESET}"
echo "    Good: manual_mode/example_output/annotation_example.tsv"
echo "    Bad:  manual_mode/example_output/annotation_with_frameshift.tsv"
echo ""
echo -e "  Next step:  ${CYAN}bash 11_self_reference_prep.sh ${SAMPLE}${RESET}"
echo ""
