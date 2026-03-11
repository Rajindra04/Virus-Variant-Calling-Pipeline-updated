#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 14: Protein Sequence Extraction
# ═══════════════════════════════════════════════════════════════
# TOOL:    extract_proteins.py
# INPUT:   Consensus FASTA(s), filtered VCF(s), config YAML
# OUTPUT:  Per-protein FASTA files, protein summary TSV
# TIME:    <1 minute
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Extracts individual protein sequences from the consensus
#   genome using gene coordinate definitions from the config.
#   For DENV, the polyprotein is cleaved into ~14 mature
#   proteins (Capsid, prM, E, NS1, NS2A, NS2B, NS3, NS4A, 2K,
#   NS4B, NS5). This step produces one multi-FASTA file per
#   protein, with sequences from all processed samples.
#
# WHY IT MATTERS:
#   Per-protein FASTAs are needed for protein-level phylogenetics,
#   structural analysis, and comparing specific proteins across
#   samples. The script is INDEL-AWARE: if upstream insertions
#   or deletions exist in the VCF, it adjusts protein boundaries
#   before slicing the consensus. Without this correction, an
#   upstream insertion would shift all downstream protein
#   boundaries by the insertion length.
#
# IMPORTANT: This step requires BOTH consensus .fa files AND
#   _filtered.vcf files in the same directory (consensus_dir).
#   The VCF is used to detect indels for coordinate adjustment.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 14: Protein Extraction — ${SEROTYPE^^}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
# This step processes ALL samples in the output directory at once.
# Both .fa (consensus) and _filtered.vcf files must be present.

CONSENSUS_COUNT=$(ls "${OUTPUT_DIR}"/*.fa 2>/dev/null | wc -l | tr -d ' ')
VCF_COUNT=$(ls "${OUTPUT_DIR}"/*_filtered.vcf 2>/dev/null | wc -l | tr -d ' ')

if [ "$CONSENSUS_COUNT" -eq 0 ]; then
    print_fail "No consensus FASTA files found in ${OUTPUT_DIR}"
    echo "  Run 08_coverage_consensus.sh for each sample first."
    exit 1
fi

print_info "Consensus FASTAs found: ${CONSENSUS_COUNT}"
print_info "Filtered VCFs found:    ${VCF_COUNT}"

if [ "$VCF_COUNT" -eq 0 ]; then
    print_warn "No filtered VCFs found — indel-aware extraction disabled"
    echo "  Protein boundaries may be inaccurate if indels exist."
fi

# ─── Output ───────────────────────────────────────────────────
PROTEIN_DIR="${OUTPUT_DIR}/proteins"
mkdir -p "$PROTEIN_DIR"

# ─── STEP 1: Run protein extraction ───────────────────────────
# PARAMETERS:
#   --consensus_dir
#       WHAT: Directory containing both .fa consensus files AND
#             _filtered.vcf files. The script pairs them by
#             sample name.
#       IMPORTANT: Must contain BOTH file types for indel-aware
#                  extraction to work.
#
#   --config_tsv
#       WHAT: Config YAML with gene_coordinates section (or a
#             config review TSV with [TRANSCRIPT_ANNOTATIONS]).
#       WHY: Defines where each protein starts and ends in the
#            reference genome coordinates.
#
#   --reference
#       WHAT: Reference FASTA for codon translation comparison.
#
#   --output_dir
#       WHAT: Where to write per-protein FASTAs and summary.

echo -e "${BLUE}[1/1]${RESET} Extracting protein sequences..."
echo "  Config:       $(basename "$CONFIG_YAML")"
echo "  Reference:    $(basename "$REFERENCE_FASTA")"
echo "  Consensus:    ${OUTPUT_DIR}/ (${CONSENSUS_COUNT} samples)"
echo ""

python3 "${PIPELINE_ROOT}/virus_pipeline/extract_proteins.py" \
    --consensus_dir "$OUTPUT_DIR" \
    --config_tsv "$CONFIG_YAML" \
    --reference "$REFERENCE_FASTA" \
    --output_dir "$PROTEIN_DIR"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: protein FASTAs and summary created             ║
# ║  ✗ FAIL if: no output files                                ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
print_header "Protein Extraction"

PROTEIN_FASTAS=$(ls "${PROTEIN_DIR}"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
SUMMARY_FILE="${PROTEIN_DIR}/protein_summary.tsv"

if [ "$PROTEIN_FASTAS" -gt 0 ]; then
    print_pass "${PROTEIN_FASTAS} protein FASTA files created"
    echo ""
    echo -e "  ${BOLD}Proteins extracted:${RESET}"
    for f in "${PROTEIN_DIR}"/*.fasta; do
        SEQS=$(grep -c "^>" "$f" || echo 0)
        echo "    $(basename "$f" .fasta): ${SEQS} sequences"
    done
else
    print_fail "No protein FASTA files created"
    exit 1
fi

if [ -f "$SUMMARY_FILE" ]; then
    echo ""
    print_pass "Protein summary: ${SUMMARY_FILE}"
    echo ""
    echo -e "  ${BOLD}Preview (first 5 rows):${RESET}"
    head -6 "$SUMMARY_FILE" | column -t -s$'\t' | sed 's/^/    /'
else
    print_warn "protein_summary.tsv not created"
fi

echo ""
echo -e "  Protein FASTAs: ${PROTEIN_DIR}/"
echo ""
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${GREEN}${BOLD}  Pipeline complete for ${SEROTYPE^^}!${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""
echo -e "  All outputs are in: ${OUTPUT_DIR}"
echo ""
echo -e "  ${BOLD}What to do next:${RESET}"
echo "    1. Review annotation files for frameshifts"
echo "    2. Check variant classification (Pass 1 vs Pass 2)"
echo "    3. Use protein FASTAs for phylogenetic analysis"
echo "    4. Process remaining serotypes (change SEROTYPE in 00_setup.sh)"
echo ""
