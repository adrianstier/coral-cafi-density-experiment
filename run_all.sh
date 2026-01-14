#!/bin/bash
# ==============================================================================
# Master Script: Run Complete CAFI Analysis Pipeline
# ==============================================================================
#
# Purpose: Execute all core analysis scripts in proper order
# Runtime: ~5-10 minutes on modern laptop
# Requirements: R >= 4.3.0 with required packages (see scripts/MRB/1.libraries.R)
#
# Usage:
#   chmod +x run_all.sh
#   ./run_all.sh
#
# For quick analysis (core scripts only):
#   ./run_all.sh --quick
#
# ==============================================================================

set -e  # Exit immediately if any command fails
set -u  # Exit if undefined variable is used

# Color output for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Navigate to repository root
cd "$(dirname "$0")"

# Function to print colored status messages
print_status() {
    echo -e "${BLUE}==>${NC} $1"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Function to run R script with timing and error handling
run_r_script() {
    local script=$1
    local description=$2

    print_status "Running: $description"
    echo "  Script: $script"

    START=$(date +%s)

    if Rscript "$script"; then
        END=$(date +%s)
        DIFF=$((END - START))
        print_success "Completed in ${DIFF}s: $description"
    else
        print_error "FAILED: $description"
        print_error "Check output above for error details"
        exit 1
    fi

    echo ""
}

# Parse command line arguments
QUICK_MODE=false
if [ $# -gt 0 ]; then
    if [ "$1" = "--quick" ]; then
        QUICK_MODE=true
        print_warning "Running in QUICK mode (core scripts only)"
        echo ""
    fi
fi

# ==============================================================================
# Main Analysis Pipeline
# ==============================================================================

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  CAFI Analysis Pipeline - Stier et al. 2026                    ║"
echo "║  DOI: 10.5281/zenodo.18239647                                  ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

TOTAL_START=$(date +%s)

# Step 1: Load Libraries and Dependencies
print_status "Step 1/5: Loading R packages and dependencies"
run_r_script "scripts/MRB/1.libraries.R" "Package loading and verification"

if [ "$QUICK_MODE" = false ]; then
    # Step 2: CAFI Community Analyses (Optional for quick mode)
    print_status "Step 2/5: CAFI community analyses"
    run_r_script "scripts/MRB/3.abundance.R" "CAFI abundance analysis"
    run_r_script "scripts/MRB/4d.diversity.R" "CAFI diversity metrics"
    run_r_script "scripts/MRB/5.fishes.R" "Fish community analysis"
    run_r_script "scripts/MRB/12.nmds_permanova_cafi.R" "Community composition (NMDS & PERMANOVA)"
else
    print_warning "Skipping CAFI community analyses (quick mode)"
    echo ""
fi

# Step 3: Coral Growth Analysis (REQUIRED)
print_status "Step 3/5: Coral growth analysis"
run_r_script "scripts/MRB/6.coral-growth.R" "Allometric growth models"

# Step 4: Coral Physiology and Performance (REQUIRED)
print_status "Step 4/5: Coral physiology and performance"
run_r_script "scripts/MRB/7.coral-physiology.R" "Physiological metrics and integrated performance"

# Step 5: CAFI-Coral Relationships (REQUIRED)
print_status "Step 5/5: CAFI-coral feedbacks"
run_r_script "scripts/MRB/8.coral-caffi.R" "Community-performance relationships"

# Compile Final Statistics
print_status "Compiling manuscript statistics"
run_r_script "scripts/MRB/14.compile-manuscript-statistics.R" "Statistical summaries"

# ==============================================================================
# Summary
# ==============================================================================

TOTAL_END=$(date +%s)
TOTAL_DIFF=$((TOTAL_END - TOTAL_START))
MINUTES=$((TOTAL_DIFF / 60))
SECONDS=$((TOTAL_DIFF % 60))

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║  Analysis Complete!                                            ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
print_success "Total runtime: ${MINUTES}m ${SECONDS}s"
echo ""

# Verify key outputs exist
print_status "Verifying outputs..."

OUTPUTS=(
    "output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.csv"
    "output/MRB/figures/coral/growth_treatment_comparison.png"
    "output/MRB/figures/coral/physiology_by_treatment.png"
    "output/MRB/figures/coral-cafi/cafi_community_vs_coral_performance.png"
)

ALL_PRESENT=true
for file in "${OUTPUTS[@]}"; do
    if [ -f "$file" ]; then
        print_success "Found: $file"
    else
        print_warning "Missing: $file (may be normal if script was modified)"
        ALL_PRESENT=false
    fi
done

echo ""

if [ "$ALL_PRESENT" = true ]; then
    print_success "All expected outputs present"
else
    print_warning "Some outputs missing - check individual script logs"
fi

echo ""
echo "Next steps:"
echo "  1. Review figures in: output/MRB/figures/"
echo "  2. Check statistics in: output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.csv"
echo "  3. Verify session info: output/MRB/objects/sessionInfo_*.txt"
echo ""
print_status "For full documentation, see: docs/REPRODUCIBILITY_GUIDE.md"
echo ""
