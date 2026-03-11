#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Module 1: Fetch Data from SRA
# ═══════════════════════════════════════════════════════════════
# TOOL:    prefetch, fasterq-dump, gzip
# INPUT:   SRR accession number(s)
# OUTPUT:  Paired FASTQ files in fastq/ directory
# TIME:    ~5-10 minutes per sample (depends on network)
# SCREEN:  Recommended for batch downloads
# ═══════════════════════════════════════════════════════════════
#
# WHAT THIS DOES:
#   Downloads raw sequencing reads from NCBI's Sequence Read
#   Archive (SRA). Each sample has two files (R1 and R2) because
#   these are paired-end reads — the same DNA fragment was
#   sequenced from both ends.
#
# WHY IT MATTERS:
#   These are the raw data. Everything downstream depends on
#   having complete, uncorrupted FASTQ files. Always verify
#   both files exist and are non-empty after download.
#
# ═══════════════════════════════════════════════════════════════

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_setup.sh" --source-only
source "${SCRIPT_DIR}/lib/colors.sh"

# ═══════════════════════════════════════════════════════════════
# All 40 samples from BioProject PRJNA1346681
# BioSample accessions: SAMN52839023–SAMN52839062
#
# To regenerate this list:
#   esearch -db sra -query PRJNA1346681 | \
#     efetch -format runinfo | cut -d',' -f1 | tail -n+2
# ═══════════════════════════════════════════════════════════════

# DENV1 samples (6)
DENV1_SAMPLES=(
    SRR35818875
    SRR35818879
    SRR35818882
    SRR35818886
    SRR35818888
    SRR35818898
)

# DENV2 samples (30)
DENV2_SAMPLES=(
    SRR35818859
    SRR35818860
    SRR35818861
    SRR35818863
    SRR35818864
    SRR35818865
    SRR35818866
    SRR35818867
    SRR35818868
    SRR35818869
    SRR35818870
    SRR35818871
    SRR35818872
    SRR35818874
    SRR35818876
    SRR35818877
    SRR35818878
    SRR35818880
    SRR35818881
    SRR35818883
    SRR35818884
    SRR35818885
    SRR35818887
    SRR35818889
    SRR35818890
    SRR35818891
    SRR35818892
    SRR35818893
    SRR35818895
    SRR35818896
)

# DENV3 samples (4)
DENV3_SAMPLES=(
    SRR35818862
    SRR35818873
    SRR35818894
    SRR35818897
)

# All samples
ALL_SAMPLES=("${DENV1_SAMPLES[@]}" "${DENV2_SAMPLES[@]}" "${DENV3_SAMPLES[@]}")

# ─── Parse arguments ──────────────────────────────────────────
usage() {
    echo "Usage: $0 [SAMPLE_ID | --all | --serotype denv1|denv2|denv3]"
    echo ""
    echo "Examples:"
    echo "  $0 SRR35818859          # Download one sample"
    echo "  $0 --serotype denv1     # Download all DENV1 samples"
    echo "  $0 --all                # Download all 40 samples (~485 MB)"
    echo ""
    echo "If no argument given, downloads the first DENV2 sample as a test."
    exit 1
}

SAMPLES_TO_FETCH=()

case "${1:-}" in
    --all)
        SAMPLES_TO_FETCH=("${ALL_SAMPLES[@]}")
        ;;
    --serotype)
        case "${2:-}" in
            denv1) SAMPLES_TO_FETCH=("${DENV1_SAMPLES[@]}") ;;
            denv2) SAMPLES_TO_FETCH=("${DENV2_SAMPLES[@]}") ;;
            denv3) SAMPLES_TO_FETCH=("${DENV3_SAMPLES[@]}") ;;
            *) echo "ERROR: Unknown serotype '${2:-}'. Use denv1, denv2, or denv3."; exit 1 ;;
        esac
        ;;
    --help|-h)
        usage
        ;;
    "")
        # Default: one test sample
        SAMPLES_TO_FETCH=("SRR35818859")
        echo -e "${YELLOW}No sample specified — downloading SRR35818859 as a test.${RESET}"
        echo ""
        ;;
    SRR*)
        SAMPLES_TO_FETCH=("$1")
        ;;
    *)
        echo "ERROR: Unrecognized argument '$1'"
        usage
        ;;
esac

# ─── Check tools ──────────────────────────────────────────────
for tool in prefetch fasterq-dump; do
    if ! command -v "$tool" &>/dev/null; then
        echo -e "${RED}ERROR: $tool not found.${RESET}"
        echo "Install SRA Toolkit: conda install -c bioconda sra-tools"
        exit 1
    fi
done

# ─── Create output directory ──────────────────────────────────
mkdir -p "${FASTQ_DIR}"

echo -e "${BOLD}Downloading ${#SAMPLES_TO_FETCH[@]} sample(s) to ${FASTQ_DIR}${RESET}"
echo ""

# ─── Download each sample ─────────────────────────────────────
SUCCESS=0
FAILED=0

for SRR in "${SAMPLES_TO_FETCH[@]}"; do
    echo -e "${CYAN}── ${SRR} ──────────────────────────────────────${RESET}"

    # Skip if already downloaded
    if [ -f "${FASTQ_DIR}/${SRR}_1.fastq.gz" ] && [ -f "${FASTQ_DIR}/${SRR}_2.fastq.gz" ]; then
        print_pass "Already exists — skipping"
        SUCCESS=$((SUCCESS + 1))
        continue
    fi

    # STEP 1: prefetch — download SRA file to local cache
    # PARAMETERS:
    #   --max-size 5G
    #       WHAT: Maximum file size to download.
    #       DEFAULT: Some SRA files can be large; 5G is generous for amplicon data.
    echo "  Prefetching..."
    if ! prefetch --max-size 5G "$SRR" 2>&1 | tail -1; then
        print_fail "prefetch failed for ${SRR}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # STEP 2: fasterq-dump — convert SRA to FASTQ
    # PARAMETERS:
    #   --split-files
    #       WHAT: Separate paired reads into _1.fastq and _2.fastq.
    #       WHY: Paired-end data must stay in two synchronized files.
    #   --outdir
    #       WHAT: Where to write FASTQ files.
    #   --threads
    #       WHAT: Parallel decompression threads.
    echo "  Converting to FASTQ..."
    if ! fasterq-dump --split-files \
            --outdir "${FASTQ_DIR}" \
            --threads "${THREADS}" \
            "$SRR" 2>&1 | tail -1; then
        print_fail "fasterq-dump failed for ${SRR}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # STEP 3: Compress FASTQ files
    # Raw FASTQ files are ~5x larger than gzipped.
    echo "  Compressing..."
    gzip -f "${FASTQ_DIR}/${SRR}_1.fastq" 2>/dev/null || true
    gzip -f "${FASTQ_DIR}/${SRR}_2.fastq" 2>/dev/null || true

    # ╔══════════════════════════════════════════════╗
    # ║  QC CHECKPOINT                               ║
    # ║  ✓ PASS if: both files exist and non-empty   ║
    # ║  ✗ FAIL if: missing or zero-size files       ║
    # ╚══════════════════════════════════════════════╝
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"

    if [ -s "$R1" ] && [ -s "$R2" ]; then
        R1_SIZE=$(du -h "$R1" | cut -f1)
        R2_SIZE=$(du -h "$R2" | cut -f1)
        print_pass "${SRR}: R1=${R1_SIZE}, R2=${R2_SIZE}"
        SUCCESS=$((SUCCESS + 1))
    else
        print_fail "${SRR}: Missing or empty FASTQ files!"
        FAILED=$((FAILED + 1))
    fi
    echo ""
done

# ─── Summary ──────────────────────────────────────────────────
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"
echo -e "  Downloaded: ${GREEN}${SUCCESS}${RESET}  Failed: ${RED}${FAILED}${RESET}"
echo -e "${BOLD}═══════════════════════════════════════════════════════════════${RESET}"

if [ "$FAILED" -gt 0 ]; then
    echo -e "${YELLOW}  Re-run with failed sample IDs to retry.${RESET}"
fi

echo ""
echo -e "  Next step:  ${CYAN}bash 02_create_samplesheet.sh${RESET}"
echo ""
