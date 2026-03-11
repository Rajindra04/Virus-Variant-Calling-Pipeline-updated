#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Terminal color definitions for the manual mode tutorial
# ═══════════════════════════════════════════════════════════════
# Source this file to get colored output in QC checkpoints.
# All colors reset with ${RESET}.
# ═══════════════════════════════════════════════════════════════

# Only enable colors if stdout is a terminal
if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    CYAN='\033[0;36m'
    BOLD='\033[1m'
    RESET='\033[0m'
else
    RED=''
    GREEN=''
    YELLOW=''
    BLUE=''
    CYAN=''
    BOLD=''
    RESET=''
fi

# Convenience functions
print_pass() { echo -e "${GREEN}${BOLD}  ✓ PASS${RESET}: $1"; }
print_warn() { echo -e "${YELLOW}${BOLD}  ⚠ WARN${RESET}: $1"; }
print_fail() { echo -e "${RED}${BOLD}  ✗ FAIL${RESET}: $1"; }
print_info() { echo -e "${BLUE}  ℹ INFO${RESET}: $1"; }
print_header() {
    echo ""
    echo -e "${CYAN}${BOLD}╔══════════════════════════════════════════════╗${RESET}"
    echo -e "${CYAN}${BOLD}║  QC CHECKPOINT: $1${RESET}"
    echo -e "${CYAN}${BOLD}╚══════════════════════════════════════════════╝${RESET}"
}
