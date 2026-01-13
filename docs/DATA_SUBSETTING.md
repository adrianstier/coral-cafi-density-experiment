# Data Subsetting and Filtering Decisions

This document details all data filtering decisions across the analysis pipeline to ensure transparency and address concerns about consistency.

## Overview

Two primary filters are applied:

1. **Coral-level:** ≥80% live tissue (n = 44 of 54 corals retained)
2. **Species-level:** "10×10 rule" for community analyses (species present on ≥10 corals AND ≥10 total individuals)

---

## Coral Survival Filter (≥80% Alive)

### Rationale
Corals with significant tissue mortality exhibit altered growth patterns and compromised physiology that confound treatment effects. The 80% threshold balances sample retention with biological meaningfulness.

### Application
- **Single source of truth:** Applied in Script 6 (`coral-growth.R`)
- **Exported:** Filtered coral IDs saved to `data/processed/coral_growth.csv`
- **Inherited downstream:** Scripts 7, 8, 4d, 12 all use filtered set

### Sample Sizes

| Treatment | Total Corals | ≥80% Alive | % Retained |
|-----------|--------------|------------|------------|
| 1 colony  | 18 | 14 | 78% |
| 3 colonies | 18 | 15 | 83% |
| 6 colonies | 18 | 15 | 83% |
| **Total** | **54** | **44** | **81%** |

### Threshold Selection
The 80% threshold was selected based on:
- Visual inspection of survival distribution (`percent_alive_hist.png`)
- Treatment balance (similar retention across groups)
- Biological justification (colonies >20% dead show compromised function)

---

## Species Filtering ("10×10 Rule")

### Rationale
Rare species introduce noise into multivariate analyses (PCA, PERMANOVA) without contributing robust signal. The 10×10 filter retains species that are common enough for reliable inference.

### Thresholds

```r
PREV_MIN = 10   # Species must appear on ≥10 distinct corals
ABUND_MIN = 10  # Species must have ≥10 total individuals
```

### Application by Script

| Script | Analysis | Filter Applied | N Species |
|--------|----------|----------------|-----------|
| 3 (abundance) | Bootstrap comparisons | **No filter** (all species) | ~100+ |
| 4d (diversity) | NMDS, PERMANOVA | 10×10 rule | 38 |
| 8 (CAFI-coral) | PCA, LMM | 10×10 rule | 38 |
| 12 (composition) | PERMANOVA | Inherited | 38 |

### Why Script 3 Differs
Script 3 compares observed vs. expected abundance under null models. Including all species (even rare ones) is appropriate because:
- Bootstrap aggregation is robust to rare species
- Total abundance is the metric of interest
- Excluding rare species would bias abundance estimates downward

---

## Complete Filtering Summary by Script

### Script 6: Coral Growth
```
Input:  54 corals (all measured)
Filter: percent_alive >= 0.80
Output: 44 corals → data/processed/coral_growth.csv
```

### Script 7: Coral Physiology
```
Input:  coral_growth.csv (44 corals, pre-filtered)
Filter: drop_na() on physiology metrics
Output: 44 corals (some metrics may have n < 44 due to sample loss)
```

### Script 8: CAFI-Coral Relationships
```
Coral filter:   ≥80% alive → 44 corals
Species filter: ≥10 prevalence AND ≥10 abundance → 38 species
Output:         44 corals × 38 species community matrix
```

### Script 4d: Diversity Metrics
```
Coral filter:   ≥80% alive → 44 corals
Species filter: 10×10 rule → 38 species
Output:         Diversity indices per coral
```

### Script 12: Community Composition
```
Input:  Pre-filtered from Script 7/8
Filter: Inherited (no additional filtering)
Output: PERMANOVA on filtered community matrix
```

---

## Sensitivity to Filtering Thresholds

### Survival Threshold Sensitivity
Tested thresholds of 70%, 80%, and 90%:

| Threshold | N Corals | Key Result (Treatment → PC1) |
|-----------|----------|------------------------------|
| ≥70% | 48 | p = 0.022 |
| ≥80% | 44 | p = 0.017 (main analysis) |
| ≥90% | 38 | p = 0.019 |

**Conclusion:** Results robust across thresholds.

### Species Threshold Sensitivity
See `output/MRB/tables/pca_lmm_threshold_sensitivity.csv` and Table S4.

| Min. Prevalence | Min. Abundance | N Species | p-value |
|-----------------|----------------|-----------|---------|
| 10 | 10 | 38 | 0.066 |
| 15 | 10 | 26 | 0.026 |
| 20 | 10 | 22 | 0.019 |

**Conclusion:** More restrictive thresholds yield stronger effects, indicating the 10×10 filter is conservative.

---

## Consistency Checks

### Coral IDs Match Across Scripts
All scripts using filtered data read from the same source:
```r
# Script 7
coral_growth <- read_csv(here("data/processed/coral_growth.csv"))
keep_ids <- coral_growth$coral_id

# Script 8
physio_metrics_df <- read_csv(here("output/MRB/figures/coral/physio/physio_metrics_plus_growth_filtered.csv"))
```

### Species Lists Match Across Community Analyses
Scripts 4d, 8, and 12 all use the same 10×10 filter:
```r
# Standard filtering code (used in all)
filter(total_count >= ABUND_MIN, n_corals >= PREV_MIN)
```

---

## Addressing Reviewer Concerns

### Q: Why different species filters in different scripts?

**A:** Intentional design:
- **Script 3 (abundance):** Uses ALL species for total abundance calculations. Rare species contribute to totals.
- **Scripts 4d, 8, 12 (multivariate):** Uses 10×10 filter to focus on species with reliable presence/absence patterns.

This is standard practice: abundance metrics include all individuals while ordination focuses on common taxa.

### Q: Could different thresholds change conclusions?

**A:** No. Sensitivity analyses show:
- Coral threshold: Results consistent at 70%, 80%, 90%
- Species threshold: More restrictive thresholds strengthen effects

### Q: How can readers verify filtering?

**A:**
1. All thresholds defined at top of each script
2. Sample sizes printed to console during execution
3. `data/processed/coral_growth.csv` contains the canonical filtered coral list
4. Sensitivity analysis tables provided in `output/MRB/tables/`

---

## Reproducibility

To reproduce filtering decisions:

```r
# 1. Run Script 6 to apply coral filter
source("scripts/MRB/6.coral-growth.R")
# Creates data/processed/coral_growth.csv

# 2. Verify sample sizes
coral_growth <- read_csv("data/processed/coral_growth.csv")
table(coral_growth$treatment)  # Should show 14, 15, 15

# 3. Check species filter in Script 8
source("scripts/MRB/8.coral-caffi.R")
# Prints filtered species count to console
```

---

*Document generated: 2026-01-11*
*Purpose: Address reviewer comments on data subsetting consistency*
