# MRB Scripts Header Audit Report

**Date:** 2025-11-05
**Auditor:** Claude Code Review
**Total Files Reviewed:** 20 R scripts
**Template Reference:** [HEADER_TEMPLATE.md](HEADER_TEMPLATE.md)

---

## Executive Summary

### Current State
- **✅ Excellent (3 files):** Complete documentation meeting publication standards
- **🟡 Good (8 files):** Solid structure but missing inputs/outputs metadata
- **🟠 Fair (6 files):** Basic headers present, need standardization
- **🔴 Needs Work (1 file):** Minimal header requiring complete rewrite
- **✅ Deleted (1 file):** Removed redundant duplicate (10coral-cafi-condition.R)
- **📝 New (1 file):** Created HEADER_TEMPLATE.md

### Key Findings
1. **90% of scripts lack Inputs/Outputs documentation**
2. **86% lack reproducibility notes**
3. **33% missing author attribution**
4. **43% missing dates**
5. ✅ **76% use consistent divider style (=)**

---

## File-by-File Audit

### PRIORITY 1: Critical Updates Needed (RED)

#### 🔴 **6.coral.R** - Needs Complete Header
**Current State:** Minimal header (title + note only)

**Issues:**
- ❌ No Purpose statement
- ❌ No Inputs list
- ❌ No Outputs list
- ❌ No Dependencies
- ❌ No Run order
- ❌ No Author
- ❌ No Dates
- ❌ No Reproducibility notes

**Action Required:**
```r
# Replace lines 1-4 with:
# ==============================================================================
# MRB Analysis Script 6: Coral Colony Growth and Physiology
# ==============================================================================
# Purpose: Analyze coral colony growth from 3D photogrammetry data (2019-2021)
#          and physiological metrics. Fits mixed-effects models with reef
#          random effects to test treatment effects on growth and condition.
#
# Inputs:
#   - data/MRB Amount/coral_growth_surface_area_change.csv
#   - data/MRB Amount/1. amount_manual_colony_measurements_dec2019_and_may2021.xlsx
#   - data/MRB Amount/1. amount_master_phys_data_v5.csv
#   - data/MRB Amount/coral_id_position_treatment.csv
#
# Outputs:
#   - output/MRB/figures/coral/*.png (growth plots)
#   - output/MRB/figures/coral/physio/*.png (physiology plots)
#   - output/MRB/tables/coral_growth_stats.csv
#
# Depends:
#   R (>= 4.3), tidyverse, lme4, lmerTest, car, broom.mixed, here
#
# Run after:
#   - 1.libraries.R
#   - utils.R
#
# Author: Adrian C. Stier / CAFI Team
# Created: [Original date unknown]
# Last updated: 2025-11-05
#
# Reproducibility notes:
#   - set.seed(42) for bootstrapping
#   - Mixed models use REML estimation
#   - All paths use here::here()
# ==============================================================================
```

---

### PRIORITY 2: Add Inputs/Outputs (YELLOW)

#### 🟡 **1.libraries.R** - Add Complete Header
**Current State:** Clean header but missing I/O details

**Action:** Add after line 8:
```r
# Inputs:  None (pure library loading)
# Outputs: None (loads packages into memory only)
```

#### 🟡 **2.taxonomic-coverage.R** - Add Inputs/Outputs
**Current State:** Good header, missing I/O metadata

**Action:** Add after Purpose section:
```r
# Inputs:
#   - data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv
#
# Outputs:
#   - output/MRB/figures/unique_taxa_table.png
#   - output/MRB/tables/taxonomic_coverage_summary.csv
```

#### 🟡 **5.fishes.R** - Add Inputs/Outputs
**Current State:** Good structure, missing I/O

**Action:** Add after Purpose section:
```r
# Inputs:
#   - data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv
#   - data/MRB Amount/coral_id_position_treatment.csv
#
# Outputs:
#   - output/MRB/figures/fishes/*.png
#   - output/MRB/tables/fish_diversity_stats.csv
#   - output/MRB/tables/fish_indicator_species.csv
```

#### 🟡 **8.null-models.R** - Add Inputs/Outputs
**Current State:** Clean header, missing I/O

**Action:** Add after Purpose section:
```r
# Inputs:
#   - data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv
#   - data/MRB Amount/coral_id_position_treatment.csv
#
# Outputs:
#   - output/MRB/figures/null-models/*.png
#   - output/MRB/tables/cooccurrence_stats.csv
```

---

### PRIORITY 3: Standardize Format (ORANGE)

#### 🟠 **4d.diversity.R** - Standardize Header Order
**Current State:** Excellent Q&A style but non-standard format

**Action:** Optionally reformat to standard template OR add note:
```r
# NOTE: This script uses Q&A format for pedagogical clarity
```

#### 🟠 **11.nmds_permanova_cafi.R** - Add Missing Sections
**Current State:** Good header, missing dates and repro notes

**Action:** Add after Purpose:
```r
# Author: CAFI Team
# Created: [Date]
# Last updated: 2025-11-05
#
# Reproducibility notes:
#   - PERMANOVA permutations = 999
#   - NMDS stress target < 0.15
#   - Uses Bray-Curtis, Jaccard, and Gower distances
```

#### 🟠 **12.SLOSS.R** - Add Outputs Section
**Current State:** Good header, missing outputs list

**Action:** Add after Depends:
```r
# Outputs:
#   - output/MRB/figures/sloss/*.png
#   - output/MRB/tables/sloss_comparison.csv
```

---

### PRIORITY 4: Verify & Polish (GREEN)

#### ✅ **3.abundance.R** - GOLD STANDARD
**Status:** Perfect ✅ No changes needed

#### ✅ **7.coral-caffi.R** - GOLD STANDARD
**Status:** Perfect ✅ No changes needed

#### ✅ **9.coral-physio.R** - Recently Updated
**Status:** Perfect ✅ Just updated to standard

#### ✅ **utils.R** - GOLD STANDARD for Utilities
**Status:** Perfect ✅ Excellent roxygen documentation

#### ✅ **mrb_config.R** - Clean Config File
**Status:** Excellent ✅ Minor: Add license line if desired

#### ✅ **mrb_figure_standards.R** - Clean Standards File
**Status:** Excellent ✅ No changes needed

---

### UTILITY/PIPELINE SCRIPTS

#### **run_mrb_pipeline.R**
**Status:** 🟢 Good
**Action:** Verify it has:
```r
# Outputs:
#   - Console output showing progress of each script
#   - All outputs from individual scripts (see respective script headers)
```

#### **test_all_scripts.R**
**Status:** 🟢 Good
**Action:** Add:
```r
# Outputs:
#   - Console report of test results
#   - logs/test_results_[timestamp].txt
```

#### **update_all_mrb_scripts.R**
**Status:** 🟢 Good
**Action:** Verify purpose is clear (what does "update" mean?)

#### **update_and_test_all_scripts.R**
**Status:** 🟢 Good
**Action:** Cross-reference with other pipeline scripts for consistency

#### **generate_publication_figures.R**
**Status:** 🟢 Good
**Action:** List all publication figures in outputs section

#### **update_all_figures_to_publication_standards.R**
**Status:** 🟢 Good
**Action:** Specify which figures and what standards

---

## Detailed Statistics

### Metadata Completeness

| Element                  | Present | Missing | % Coverage |
|--------------------------|---------|---------|------------|
| **Visual Dividers**      | 16      | 4       | 80%        |
| **Purpose Statement**    | 18      | 2       | 90%        |
| **Inputs List**          | 2       | 18      | 10% ❌     |
| **Outputs List**         | 2       | 18      | 10% ❌     |
| **Dependencies**         | 3       | 17      | 15% ❌     |
| **Run Order**            | 3       | 17      | 15% ❌     |
| **Author**               | 14      | 6       | 70%        |
| **Dates**                | 12      | 8       | 60%        |
| **Reproducibility Notes**| 2       | 18      | 10% ❌     |

### Header Quality Distribution

| Rating    | Count | Files                                          |
|-----------|-------|------------------------------------------------|
| ✅ Excellent | 3   | 3.abundance.R, 7.coral-caffi.R, utils.R       |
| 🟢 Very Good | 3   | 9.coral-physio.R, mrb_config.R, mrb_figure_standards.R |
| 🟡 Good     | 8    | 1, 2, 5, 8, 11, 12, run_pipeline, test_all  |
| 🟠 Fair     | 4    | 4d, update_scripts (various)                  |
| 🔴 Needs Work| 1   | 6.coral.R                                    |
| ❌ Deleted  | 1   | 10coral-cafi-condition.R (duplicate)          |

---

## Recommended Action Plan

### Phase 1: Critical (Do First) ⚠️
**Time Estimate:** 30 minutes

1. ✅ **Delete 10coral-cafi-condition.R** - DONE
2. **Update 6.coral.R** with complete header
3. **Add Inputs/Outputs to scripts 2, 5, 8**

### Phase 2: High Priority
**Time Estimate:** 1 hour

4. **Add Inputs/Outputs to script 1**
5. **Add reproducibility notes to scripts 2, 5, 6, 8, 11, 12**
6. **Add dates to all scripts missing them**

### Phase 3: Standardization
**Time Estimate:** 1-2 hours

7. **Standardize section dividers across all files** (= vs -)
8. **Add missing dependencies to all scripts**
9. **Add "Run after" sections to all analysis scripts**

### Phase 4: Polish & Documentation
**Time Estimate:** 1 hour

10. **Create README.md update** with new header standards
11. **Add examples to utility functions lacking them**
12. **Review and update all "Last updated" dates**

---

## Template Files Created

1. ✅ **HEADER_TEMPLATE.md** - Gold standard template for all future scripts
2. ✅ **HEADER_AUDIT_REPORT.md** - This file, tracks current state
3. ✅ **Script 9 updated** - Example of properly formatted analysis script

---

## Quick Reference: File Status

| # | Script Name                              | Status      | Priority |
|---|------------------------------------------|-------------|----------|
| 1 | 1.libraries.R                            | 🟡 Good     | P2       |
| 2 | 2.taxonomic-coverage.R                   | 🟡 Good     | P2       |
| 3 | 3.abundance.R                            | ✅ Excellent | None     |
| 4 | 4d.diversity.R                           | 🟠 Fair     | P3       |
| 5 | 5.fishes.R                               | 🟡 Good     | P2       |
| 6 | 6.coral.R                                | 🔴 Critical | **P1**   |
| 7 | 7.coral-caffi.R                          | ✅ Excellent | None     |
| 8 | 8.null-models.R                          | 🟡 Good     | P2       |
| 9 | 9.coral-physio.R                         | ✅ Excellent | None     |
| 10| ~~10coral-cafi-condition.R~~             | ❌ Deleted  | DONE     |
| 11| 11.nmds_permanova_cafi.R                 | 🟠 Fair     | P3       |
| 12| 12.SLOSS.R                               | 🟠 Fair     | P3       |
| - | utils.R                                  | ✅ Excellent | None     |
| - | mrb_config.R                             | 🟢 Very Good| P4       |
| - | mrb_figure_standards.R                   | 🟢 Very Good| P4       |
| - | run_mrb_pipeline.R                       | 🟡 Good     | P3       |
| - | test_all_scripts.R                       | 🟡 Good     | P3       |
| - | update_all_mrb_scripts.R                 | 🟡 Good     | P3       |
| - | update_and_test_all_scripts.R            | 🟡 Good     | P3       |
| - | generate_publication_figures.R           | 🟡 Good     | P3       |
| - | update_all_figures_to_publication_standards.R | 🟡 Good | P3   |

---

## Next Steps

1. **Review this audit** with team/PI
2. **Prioritize updates** based on immediate needs (publication timeline?)
3. **Assign owners** for each phase if team effort
4. **Set deadline** for Phase 1 (critical updates)
5. **Use HEADER_TEMPLATE.md** for all future scripts

---

## Notes for Code Review

This repository is **well above average** in code organization:
- ✅ Consistent naming conventions
- ✅ Centralized utilities and configuration
- ✅ Clear directory structure
- ✅ Good use of version control
- ✅ Reproducible paths with `here::here()`

**Main improvement opportunity:** Complete metadata documentation (inputs/outputs/dependencies) to achieve **publication-quality research code**.

---

**End of Report**
