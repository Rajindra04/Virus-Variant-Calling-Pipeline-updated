#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# QC Checkpoint Functions
# ═══════════════════════════════════════════════════════════════
# Source this file to get QC validation functions.
# Each function parses tool output, prints PASS/WARN/FAIL,
# and exits non-zero on FAIL.
#
# Usage:
#   source lib/qc_check.sh
#   qc_check_fastp  output_dir/sample_fastp.json
#   qc_check_mapping  output_dir/sample.flagstat
#   ...
# ═══════════════════════════════════════════════════════════════

SCRIPT_DIR_QC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR_QC}/colors.sh"

# ───────────────────────────────────────────────────────────────
# qc_check_fastp <fastp_json>
#   PASS: >80% reads surviving AND Q30 >85%
#   WARN: 70-80% surviving OR Q30 80-85%
#   FAIL: <70% surviving OR Q30 <70%
# ───────────────────────────────────────────────────────────────
qc_check_fastp() {
    local json="${1:?Usage: qc_check_fastp <fastp.json>}"
    print_header "fastp QC"

    if [ ! -f "$json" ]; then
        print_fail "fastp JSON not found: $json"
        return 1
    fi

    # Extract metrics using python (available in conda env)
    local metrics
    metrics=$(python3 -c "
import json, sys
with open('${json}') as f:
    d = json.load(f)
before = d['summary']['before_filtering']['total_reads']
after  = d['summary']['after_filtering']['total_reads']
q30    = d['summary']['after_filtering']['q30_rate'] * 100
surv   = (after / before * 100) if before > 0 else 0
print(f'{surv:.1f} {q30:.1f} {before} {after}')
")

    local surv q30 before after
    read -r surv q30 before after <<< "$metrics"

    print_info "Reads before filtering: ${before}"
    print_info "Reads after filtering:  ${after}"
    print_info "Surviving: ${surv}%"
    print_info "Q30 rate:  ${q30}%"

    local status=0

    # Check survival rate
    if (( $(echo "$surv < 70" | bc -l) )); then
        print_fail "Survival rate ${surv}% < 70% — too many reads lost"
        status=1
    elif (( $(echo "$surv < 80" | bc -l) )); then
        print_warn "Survival rate ${surv}% is between 70-80%"
    else
        print_pass "Survival rate ${surv}% ≥ 80%"
    fi

    # Check Q30
    if (( $(echo "$q30 < 70" | bc -l) )); then
        print_fail "Q30 rate ${q30}% < 70% — very low base quality"
        status=1
    elif (( $(echo "$q30 < 85" | bc -l) )); then
        print_warn "Q30 rate ${q30}% is between 70-85%"
    else
        print_pass "Q30 rate ${q30}% ≥ 85%"
    fi

    return $status
}

# ───────────────────────────────────────────────────────────────
# qc_check_mapping <flagstat_file>
#   PASS: mapping rate >80%
#   WARN: 50-80%
#   FAIL: <50%
# ───────────────────────────────────────────────────────────────
qc_check_mapping() {
    local flagstat="${1:?Usage: qc_check_mapping <flagstat_file>}"
    print_header "Mapping QC"

    if [ ! -f "$flagstat" ]; then
        print_fail "Flagstat file not found: $flagstat"
        return 1
    fi

    # Parse mapping rate from samtools flagstat output
    local total mapped pct
    total=$(grep "in total" "$flagstat" | head -1 | awk '{print $1}')
    mapped=$(grep "mapped (" "$flagstat" | head -1 | awk '{print $1}')

    if [ "$total" -eq 0 ] 2>/dev/null; then
        print_fail "No reads in flagstat file"
        return 1
    fi

    pct=$(echo "scale=1; $mapped * 100 / $total" | bc)

    print_info "Total reads:  ${total}"
    print_info "Mapped reads: ${mapped}"
    print_info "Mapping rate: ${pct}%"

    if (( $(echo "$pct < 50" | bc -l) )); then
        print_fail "Mapping rate ${pct}% < 50% — wrong reference or contamination?"
        return 1
    elif (( $(echo "$pct < 80" | bc -l) )); then
        print_warn "Mapping rate ${pct}% is between 50-80%"
        return 0
    else
        print_pass "Mapping rate ${pct}% ≥ 80%"
        return 0
    fi
}

# ───────────────────────────────────────────────────────────────
# qc_check_dedup <dedup_stats_file>
#   PASS: duplication rate <50%
#   WARN: 50-80%
#   FAIL: >80%
# ───────────────────────────────────────────────────────────────
qc_check_dedup() {
    local stats="${1:?Usage: qc_check_dedup <markdup_stats_file>}"
    print_header "Deduplication QC"

    if [ ! -f "$stats" ]; then
        print_fail "Dedup stats file not found: $stats"
        return 1
    fi

    # samtools markdup -s writes stats to stderr, captured to file
    # Look for "DUPLICATE" line or parse the summary
    local total dups pct
    # Try parsing samtools markdup stats format
    if grep -q "DUPLICATE TOTAL" "$stats"; then
        # samtools markdup -f format
        total=$(grep "^EXAMINED" "$stats" | awk '{print $2}')
        dups=$(grep "^DUPLICATE TOTAL" "$stats" | awk '{print $2}')
    else
        # Fallback: parse from "duplicate" percentage line
        total=$(grep -i "read" "$stats" | head -1 | awk '{print $1}')
        dups=$(grep -i "duplicate" "$stats" | head -1 | awk '{print $1}')
    fi

    if [ -z "$total" ] || [ "$total" -eq 0 ] 2>/dev/null; then
        print_warn "Could not parse dedup stats — check file format"
        return 0
    fi

    pct=$(echo "scale=1; $dups * 100 / $total" | bc)

    print_info "Total reads examined: ${total}"
    print_info "Duplicates found:     ${dups}"
    print_info "Duplication rate:     ${pct}%"

    if (( $(echo "$pct > 80" | bc -l) )); then
        print_fail "Duplication rate ${pct}% > 80% — library complexity very low"
        return 1
    elif (( $(echo "$pct > 50" | bc -l) )); then
        print_warn "Duplication rate ${pct}% is between 50-80%"
        return 0
    else
        print_pass "Duplication rate ${pct}% < 50%"
        return 0
    fi
}

# ───────────────────────────────────────────────────────────────
# qc_check_coverage <coverage_file> <consensus_fasta> <genome_size>
#   PASS: ≥90% genome at ≥20x AND <5% N content
#   WARN: 70-90% at ≥20x OR 5-20% N
#   FAIL: <70% at ≥20x OR >20% N
# ───────────────────────────────────────────────────────────────
qc_check_coverage() {
    local coverage="${1:?Usage: qc_check_coverage <coverage.txt> <consensus.fa> <genome_size>}"
    local consensus="${2:?}"
    local genome_size="${3:?}"
    print_header "Coverage QC"

    if [ ! -f "$coverage" ]; then
        print_fail "Coverage file not found: $coverage"
        return 1
    fi
    if [ ! -f "$consensus" ]; then
        print_fail "Consensus FASTA not found: $consensus"
        return 1
    fi

    local status=0

    # Calculate % genome at ≥20x
    local bases_20x pct_20x
    bases_20x=$(awk '$3 >= 20 {count++} END {print count+0}' "$coverage")
    pct_20x=$(echo "scale=1; $bases_20x * 100 / $genome_size" | bc)

    print_info "Genome size:    ${genome_size} bp"
    print_info "Bases at ≥20x:  ${bases_20x}"
    print_info "Genome at ≥20x: ${pct_20x}%"

    if (( $(echo "$pct_20x < 70" | bc -l) )); then
        print_fail "Genome coverage ${pct_20x}% < 70% at ≥20x"
        status=1
    elif (( $(echo "$pct_20x < 90" | bc -l) )); then
        print_warn "Genome coverage ${pct_20x}% is between 70-90% at ≥20x"
    else
        print_pass "Genome coverage ${pct_20x}% ≥ 90% at ≥20x"
    fi

    # Calculate N content in consensus
    local total_bases n_count n_pct
    total_bases=$(grep -v "^>" "$consensus" | tr -d '\n' | wc -c | tr -d ' ')
    n_count=$(grep -v "^>" "$consensus" | tr -d '\n' | tr -cd 'Nn' | wc -c | tr -d ' ')
    n_pct=$(echo "scale=1; $n_count * 100 / $total_bases" | bc)

    print_info "Consensus length: ${total_bases} bp"
    print_info "N positions:      ${n_count}"
    print_info "N content:        ${n_pct}%"

    if (( $(echo "$n_pct > 20" | bc -l) )); then
        print_fail "N content ${n_pct}% > 20% — too many uncalled bases"
        status=1
    elif (( $(echo "$n_pct > 5" | bc -l) )); then
        print_warn "N content ${n_pct}% is between 5-20%"
    else
        print_pass "N content ${n_pct}% < 5%"
    fi

    return $status
}

# ───────────────────────────────────────────────────────────────
# qc_check_variants <filtered_vcf>
#   PASS: 5-50 variants
#   WARN: 50-200 variants
#   FAIL: >200 or 0 variants
# ───────────────────────────────────────────────────────────────
qc_check_variants() {
    local vcf="${1:?Usage: qc_check_variants <filtered.vcf>}"
    print_header "Variant Count QC"

    if [ ! -f "$vcf" ]; then
        print_fail "VCF file not found: $vcf"
        return 1
    fi

    local count
    count=$(grep -v "^#" "$vcf" | wc -l | tr -d ' ')

    print_info "Variants passing filters: ${count}"

    if [ "$count" -eq 0 ]; then
        print_fail "0 variants found — no variation from reference?"
        return 1
    elif [ "$count" -gt 200 ]; then
        print_fail "${count} variants > 200 — possible contamination or wrong reference"
        return 1
    elif [ "$count" -gt 50 ]; then
        print_warn "${count} variants is between 50-200"
        return 0
    else
        print_pass "${count} variants in expected range (5-50)"
        return 0
    fi
}

# ───────────────────────────────────────────────────────────────
# qc_check_frameshifts <annotation_tsv>
#   PASS: 0 frameshifts
#   WARN: 1 frameshift
#   FAIL: >1 frameshift — THE KEY LESSON
# ───────────────────────────────────────────────────────────────
qc_check_frameshifts() {
    local tsv="${1:?Usage: qc_check_frameshifts <annotations.tsv>}"
    print_header "Frameshift QC (CRITICAL)"

    if [ ! -f "$tsv" ]; then
        print_fail "Annotation TSV not found: $tsv"
        return 1
    fi

    local count
    count=$(awk -F'\t' '$12 ~ /frameshift/' "$tsv" | wc -l | tr -d ' ')

    print_info "Frameshift variants found: ${count}"

    if [ "$count" -gt 1 ]; then
        echo ""
        echo -e "${RED}${BOLD}  ╔══════════════════════════════════════════════════════════════╗${RESET}"
        echo -e "${RED}${BOLD}  ║                    *** CRITICAL WARNING ***                  ║${RESET}"
        echo -e "${RED}${BOLD}  ║                                                              ║${RESET}"
        echo -e "${RED}${BOLD}  ║  ${count} frameshift variants detected!                           ║${RESET}"
        echo -e "${RED}${BOLD}  ║                                                              ║${RESET}"
        echo -e "${RED}${BOLD}  ║  Dengue virus encodes a SINGLE polyprotein.                  ║${RESET}"
        echo -e "${RED}${BOLD}  ║  A real frameshift would destroy ALL downstream proteins     ║${RESET}"
        echo -e "${RED}${BOLD}  ║  and produce a non-viable virus.                             ║${RESET}"
        echo -e "${RED}${BOLD}  ║                                                              ║${RESET}"
        echo -e "${RED}${BOLD}  ║  Multiple frameshifts almost certainly indicate:              ║${RESET}"
        echo -e "${RED}${BOLD}  ║    • Sequencing artifacts (homopolymer errors)                ║${RESET}"
        echo -e "${RED}${BOLD}  ║    • Mapping errors near primer binding sites                 ║${RESET}"
        echo -e "${RED}${BOLD}  ║    • Incorrect variant calling parameters                    ║${RESET}"
        echo -e "${RED}${BOLD}  ║                                                              ║${RESET}"
        echo -e "${RED}${BOLD}  ║  DO NOT publish these results without manual review!         ║${RESET}"
        echo -e "${RED}${BOLD}  ╚══════════════════════════════════════════════════════════════╝${RESET}"
        echo ""
        echo -e "  Frameshifts found:"
        awk -F'\t' '$12 ~ /frameshift/ {printf "    POS=%s  REF=%s  ALT=%s  GENE=%s  EFFECT=%s\n", $2, $4, $5, $14, $12}' "$tsv"
        echo ""
        return 1
    elif [ "$count" -eq 1 ]; then
        print_warn "1 frameshift found — review manually (could be real but rare)"
        awk -F'\t' '$12 ~ /frameshift/ {printf "    POS=%s  REF=%s  ALT=%s  GENE=%s\n", $2, $4, $5, $14}' "$tsv"
        return 0
    else
        print_pass "0 frameshifts — consistent with viable single-polyprotein virus"
        return 0
    fi
}
