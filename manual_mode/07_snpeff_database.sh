#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 7: Build snpEff Database
# ═══════════════════════════════════════════════════════════════
# TOOL:    snpEff
# INPUT:   GenBank file (.gb)
# OUTPUT:  snpEff database for the reference genome
# TIME:    ~1-2 minutes
# SCREEN:  Not required
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Builds a local snpEff annotation database from a GenBank
#   file. This database maps genome coordinates to gene names,
#   protein products, and reading frames — so variant calls can
#   be annotated as "missense in NS3" instead of just "A→G at
#   position 4750."
#
# WHY IT MATTERS:
#   Without annotation, you have a list of positions and base
#   changes. With annotation, you know which protein is affected
#   and whether the amino acid changed. This is essential for
#   biological interpretation.
#
# NOTE: Run this ONCE per serotype, not per sample. If you
#   already built the database for DENV2, skip this when
#   processing additional DENV2 samples.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Module 7: Build snpEff Database — ${DATABASE_NAME}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo ""

# ─── Check if database already exists ─────────────────────────
SNPEFF_DATA_DIR=$(python3 -c "
import subprocess, os
result = subprocess.run(['snpEff', '-version'], capture_output=True, text=True)
# Find snpEff data directory
import glob
candidates = glob.glob(os.path.expanduser('~') + '/**/snpEff/data', recursive=True)
if not candidates:
    # Check conda env
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    candidates = glob.glob(conda_prefix + '/share/snpeff-*/data')
if candidates:
    print(candidates[0])
else:
    print('')
" 2>/dev/null || echo "")

# Alternative: find via snpEff config
if [ -z "$SNPEFF_DATA_DIR" ]; then
    SNPEFF_CONFIG=$(find "${CONDA_PREFIX:-/usr}" -name "snpEff.config" 2>/dev/null | head -1)
    if [ -n "$SNPEFF_CONFIG" ]; then
        SNPEFF_DATA_DIR="$(dirname "$SNPEFF_CONFIG")/data"
    fi
fi

if [ -d "${SNPEFF_DATA_DIR}/${DATABASE_NAME}" ]; then
    print_pass "Database '${DATABASE_NAME}' already exists"
    echo "  Location: ${SNPEFF_DATA_DIR}/${DATABASE_NAME}"
    echo "  Skipping build. Delete the directory to force rebuild."
    exit 0
fi

# ─── Verify GenBank file ──────────────────────────────────────
if [ ! -f "$GENBANK_FILE" ]; then
    print_fail "GenBank file not found: $GENBANK_FILE"
    exit 1
fi
print_pass "GenBank file: $(basename "$GENBANK_FILE")"

# ─── STEP 1: Build database using pipeline's Python script ────
# The pipeline includes create_snpeff_database.py which handles
# all the snpEff configuration and database building.
#
# WHAT IT DOES:
#   1. Creates a directory in snpEff's data folder
#   2. Copies the GenBank file as genes.gbk
#   3. Adds the genome entry to snpEff.config
#   4. Runs snpEff build -genbank

echo -e "${BLUE}[1/1]${RESET} Building snpEff database..."
echo "  Database name: ${DATABASE_NAME}"
echo "  GenBank file:  $(basename "$GENBANK_FILE")"
echo ""

python3 "${PIPELINE_ROOT}/virus_pipeline/create_snpeff_database.py" \
    --genbank_file "$GENBANK_FILE" \
    --database_name "$DATABASE_NAME"

# ╔══════════════════════════════════════════════════════════════╗
# ║  QC CHECKPOINT                                              ║
# ║  ✓ PASS if: snpEff dump shows genes and proteins            ║
# ║  ✗ FAIL if: database not found or empty                    ║
# ╚══════════════════════════════════════════════════════════════╝

echo ""
print_header "snpEff Database"

# Verify by dumping database info
if snpEff dump "$DATABASE_NAME" 2>/dev/null | head -5 | grep -q "Genome"; then
    GENE_COUNT=$(snpEff dump "$DATABASE_NAME" 2>/dev/null | grep -c "Gene:" || echo 0)
    print_pass "Database '${DATABASE_NAME}' built successfully (${GENE_COUNT} genes)"
else
    # Try a simpler check
    if snpEff databases 2>/dev/null | grep -q "$DATABASE_NAME"; then
        print_pass "Database '${DATABASE_NAME}' registered in snpEff"
    else
        print_warn "Could not verify database — try: snpEff dump ${DATABASE_NAME}"
    fi
fi

echo ""
echo -e "  This database is reused for ALL ${SEROTYPE^^} samples."
echo ""
echo -e "  Next step:  ${CYAN}bash 08_coverage_consensus.sh <sample_name>${RESET}"
echo ""
