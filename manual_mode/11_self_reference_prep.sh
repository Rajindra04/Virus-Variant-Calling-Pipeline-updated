#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 11: Self-Reference Preparation (Pass 2 Setup)
# ═══════════════════════════════════════════════════════════════
# TOOL:    fill_consensus_ns.py, bwa-mem2 index, samtools faidx,
#          gatk CreateSequenceDictionary
# INPUT:   Consensus FASTA from Pass 1, reference FASTA
# OUTPUT:  N-filled, indexed self-reference genome
# TIME:    ~2-3 minutes per sample
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Takes the consensus sequence from Pass 1 and prepares it as
#   a custom reference genome for Pass 2. Positions called as "N"
#   (insufficient coverage) are filled with the original reference
#   base, then the genome is indexed for alignment.
#
# WHY IT MATTERS:
#   Pass 1 maps reads against a reference (e.g., NC_001474.2)
#   that may differ from your sample at dozens of positions.
#   These fixed differences create noise when looking for
#   intrahost variants. By re-mapping to the sample's OWN
#   consensus, you separate:
#     • FIXED variants (sample differs from reference)
#     • INTRAHOST variants (minority populations within sample)
#   This two-pass approach is essential for accurate intrahost
#   diversity analysis.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

SAMPLE="${1:?Usage: $0 <sample_name>}"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 11: Self-Reference Prep — ${SAMPLE}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Input ────────────────────────────────────────────────────
CONSENSUS="${OUTPUT_DIR}/${SAMPLE}.fa"
if [ ! -f "$CONSENSUS" ]; then
    print_fail "Consensus FASTA not found: $CONSENSUS"
    echo "  Run 08_coverage_consensus.sh ${SAMPLE} first."
    exit 1
fi

# ─── Output ───────────────────────────────────────────────────
SELFREF_DIR="${PASS2_DIR}/references"
mkdir -p "$SELFREF_DIR"
SELFREF_FASTA="${SELFREF_DIR}/${SAMPLE}.fasta"

# ─── STEP 1: Fill Ns with reference bases ─────────────────────
# PARAMETERS:
#   --consensus: The Pass 1 consensus (may have N positions).
#   --reference: The original reference (NC_001474.2 etc).
#   --output: The "filled" FASTA — no Ns, ready for indexing.
#   --sample_name: Used as the FASTA header.
#
# WHY FILL Ns:
#   bwa-mem2 cannot align reads to "N" bases — they would create
#   unmappable gaps in the reference. Filling with the original
#   reference base is the best guess (the sample probably matches
#   the reference at low-coverage positions).

echo -e "${BLUE}[1/4]${RESET} Filling N positions with reference bases..."
python3 "${PIPELINE_ROOT}/virus_pipeline/fill_consensus_ns.py" \
    --consensus "$CONSENSUS" \
    --reference "$REFERENCE_FASTA" \
    --output "$SELFREF_FASTA" \
    --sample_name "$SAMPLE"

# Report N filling
TOTAL_BASES=$(grep -v "^>" "$CONSENSUS" | tr -d '\n' | wc -c | tr -d ' ')
N_COUNT=$(grep -v "^>" "$CONSENSUS" | tr -d '\n' | tr -cd 'Nn' | wc -c | tr -d ' ')
N_PCT=$(echo "scale=1; $N_COUNT * 100 / $TOTAL_BASES" | bc)
print_info "Filled ${N_COUNT} N positions (${N_PCT}% of genome)"

# ─── STEP 2: Index with samtools ──────────────────────────────
echo -e "${BLUE}[2/4]${RESET} Creating samtools index..."
samtools faidx "$SELFREF_FASTA"

# ─── STEP 3: Index with bwa-mem2 ──────────────────────────────
echo -e "${BLUE}[3/4]${RESET} Creating bwa-mem2 index..."
bwa-mem2 index "$SELFREF_FASTA" 2>&1 | tail -1

# ─── STEP 4: Create GATK sequence dictionary ──────────────────
echo -e "${BLUE}[4/4]${RESET} Creating sequence dictionary..."
gatk CreateSequenceDictionary \
    -R "$SELFREF_FASTA" \
    2>&1 | tail -1

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: all index files created                        ║
# ║  ✗ FAIL if: any index missing                              ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
print_header "Self-Reference"

ALL_OK=true
for ext in ".fasta" ".fasta.fai" ".dict"; do
    EXPECTED="${SELFREF_DIR}/${SAMPLE}${ext}"
    if [ "$ext" = ".dict" ]; then
        EXPECTED="${SELFREF_FASTA%.fasta}.dict"
    fi
    if [ -f "$EXPECTED" ]; then
        print_pass "$(basename "$EXPECTED")"
    else
        print_fail "Missing: $(basename "$EXPECTED")"
        ALL_OK=false
    fi
done

# Check bwa-mem2 index files
if ls "${SELFREF_FASTA}"*.bwt* &>/dev/null; then
    print_pass "bwa-mem2 index files"
else
    print_fail "bwa-mem2 index files missing"
    ALL_OK=false
fi

if [ "$ALL_OK" = false ]; then
    exit 1
fi

echo ""
echo -e "  Self-reference: ${SELFREF_FASTA}"
echo ""
echo -e "  Next step:  ${CYAN}bash 12_self_reference_calling.sh ${SAMPLE}${RESET}"
echo ""
