# MRB Figure Reorganization & Color Standardization Plan

**Date:** 2025-11-05
**Status:** 🚧 IN PROGRESS

---

## Executive Summary

This document outlines the complete reorganization of the MRB figure output directory structure and color standardization across all R scripts. The current system has:

1. **Inconsistent directory organization** (figures in both `output/diversity` and `output/MRB/figures/`)
2. **Color scheme inconsistencies** in 3 high-priority scripts
3. **Unintuitive naming** that makes it difficult to find figures

---

## Current Issues

### 1. Directory Structure Problems

**Current fragmented structure:**
```
output/
├── diversity/              ← 14 figures (should be in output/MRB/figures/)
├── figures/                ← 2 figures (unclear purpose)
└── MRB/
    └── figures/
        ├── abundance/      ← 6 figures
        ├── combined/       ← 35 figures (too many, unclear grouping)
        ├── coral/          ← 6 figures
        │   └── physio/     ← 9 figures
        ├── fishes/         ← 4 figures
        ├── null_models/    ← 2 figures
        ├── ordination/     ← 5 figures
        └── taxonomic_coverage/ ← 4 figures
```

**Problems:**
- Figures split between `output/diversity` and `output/MRB/figures/` with no clear logic
- `output/figures` at root level creates confusion
- `output/MRB/figures/combined` is a catch-all with 35 figures
- No clear mapping between scripts and figure locations

### 2. Color Scheme Issues

**Scripts with inconsistent colors (from Task agent analysis):**

| Priority | Script | Issue | Figures Affected |
|----------|--------|-------|------------------|
| **HIGH** | `12.nmds_permanova_cafi.R` | Custom RColorBrewer "Dark2" palette instead of TREATMENT_COLORS | 5 ordination figures |
| **HIGH** | `10.coral-physio.R` | Doesn't source mrb_figure_standards.R, uses theme_minimal | 35 combined figures |
| **MEDIUM** | `13.SLOSS.R` | Doesn't source mrb_figure_standards.R | Not yet generated |

**Color palette comparison:**
- ✅ **Standard (correct):** Treatment 1=#E69F00 (Orange), Treatment 3=#56B4E9 (Sky Blue), Treatment 6=#009E73 (Green)
- ❌ **RColorBrewer Dark2 (wrong):** #1B9E77 (Teal), #D95F02 (Orange-Red), #7570B3 (Purple)

---

## Proposed New Directory Structure

### Intuitive, Script-Aligned Organization

```
output/MRB/
├── 01_taxonomic_coverage/      # From script 2
│   ├── abundance_by_resolution.png
│   ├── species_richness_by_order.png
│   ├── taxonomic_completeness.png
│   └── unique_taxa_table.png
│
├── 02_abundance/                # From script 3
│   ├── community_abundance_observed_bootstrap.png
│   ├── community_total_abundance.png
│   ├── focal_order_species_np.png
│   ├── rarefied_richness.png
│   ├── species_richness.png
│   └── top30_observed_vs_expected_all.png
│
├── 03_diversity/                # From script 4d
│   ├── alpha/
│   │   ├── 04_alpha_diversity.png
│   │   └── 05_rank_abundance.png
│   ├── beta/
│   │   ├── 06_nmds.png
│   │   ├── 14_coral_nmds_centroids_spiders.png
│   │   ├── 15_between_reef_distance_boxplots.png
│   │   ├── 15_pairwise_permanova_heatmaps.png
│   │   ├── 16_gower_2x2_coral_reef_abund_prop.png
│   │   ├── 16_gower_2x3_with_density.png
│   │   ├── 16_reef_abundance_vs_density.png
│   │   └── 16_vertical_2panel_nmds_prop_plus_top15_density.png
│   ├── methods/
│   │   ├── section17_bars_technique_metric.png
│   │   ├── section17_hist_pvalues.png
│   │   └── section17_sensitivity_curves.png
│   └── supplemental/
│       └── nmds_gower_density_square.png
│
├── 04_fishes/                   # From script 5
│   ├── fish_diversity_metrics.png
│   ├── fish_nmds_ordination.png
│   ├── fish_species_accumulation.png
│   └── top_fish_species_abundance.png
│
├── 05_coral_growth/             # From script 6
│   ├── ANCOVA_Init_vs_Final_Volume_by_Treatment.png
│   ├── DeltaVolume_vs_SA_ANCOVA.png
│   ├── SA_Scaled_Volume_Growth_by_Treatment.png
│   ├── SizeCorrected_Volume_Growth_by_Treatment.png
│   ├── percent_alive_by_treatment.png
│   └── percent_alive_hist.png
│
├── 06_coral_physiology/         # From script 7
│   ├── pc1_loadings_and_scores_paired.png
│   ├── pca_scores_by_treatment.png
│   ├── physio_by_treatment.png
│   ├── physio_growth_pairs.png
│   ├── physio_growth_pca_biplot.png
│   ├── physio_growth_pca_loadings.png
│   ├── physio_growth_pca_scree.png
│   ├── univariate_anova_treatment_effects.png
│   └── univariate_metric_by_treatment.png
│
├── 07_community_coral_integration/ # From script 8 (8.coral-caffi.R)
│   ├── pca_analysis/
│   │   ├── PC1_loadings_COMM_HELLINGER.png
│   │   ├── PC1_loadings_COMM_RAW.png
│   │   ├── PC1_loadings_COMM_SQRT.png
│   │   ├── PC1_loadings_COMM_SQRT_CS.png
│   │   ├── PC1_loadings_COMM_SQRT_CS_aligned.png
│   │   ├── PC1_loadings_COMM_SQRT_aligned.png
│   │   ├── PCA_9panel_overview.png
│   │   ├── PCA_comm_vs_cond_PC1.png
│   │   ├── SQRT_CS_PC1loadings_and_scatter.png
│   │   ├── scree_COMM_HELLINGER.png
│   │   ├── scree_COMM_RAW.png
│   │   ├── scree_COMM_SQRT.png
│   │   ├── scree_COMM_SQRT_CS.png
│   │   └── scree_COND_z.png
│   ├── publication_panels/
│   │   ├── Final_3panel.png
│   │   ├── Loadings_panel.png
│   │   ├── PCA_LOADINGS_HELLINGER_3panel_clean.png
│   │   ├── PCA_LOADINGS_RAW_2panel_clean.png
│   │   ├── PCA_LOADINGS_RAW_3panel_clean.png
│   │   └── panel_6b_oriented.png
│   ├── species_analysis/
│   │   ├── LMM_species_top20_coefplot.png
│   │   ├── species_corr_PC1_vs_RDA1_scatter.png
│   │   ├── species_faceted_LMM_lines_rawX_sqrtAxis.png
│   │   ├── species_faceted_LMM_lines_rawX_sqrtAxis_wide_3x2.png
│   │   ├── species_performance_corr_heatmap_sqrt_byGroup_orderedByPC1.png
│   │   ├── species_stability.png
│   │   ├── species_top20_LMM_panels.png
│   │   ├── species_vs_commPC1_corr.png
│   │   └── species_vs_RDA1_corr.png
│   ├── model_comparison/
│   │   ├── LMM_COMM_vs_COND_comparison.png
│   │   ├── combined_sensitivity_overview.png
│   │   ├── combined_species_corr_panel.png
│   │   ├── corr_vs_adjR2_scatter.png
│   │   ├── multi_metric_heatmaps.png
│   │   └── sample_taxa_richness.png
│   └── networks/
│       ├── sem_network_all_species_r030_reefAdjusted.png
│       └── sem_network_physio_panel.png
│
├── 08_null_models/              # From script 9
│   ├── bootstrap_distributions.png
│   └── global_null_distributions.png
│
├── 09_ordination/               # From script 12
│   ├── NMDS_bray.png
│   ├── NMDS_gower.png
│   ├── NMDS_jaccard.png
│   ├── subsample_prop_sig.png
│   └── subsample_pvalue_distributions.png
│
└── 10_sloss/                    # From script 13 (when generated)
    └── (figures TBD)
```

### Key Improvements

1. **Sequential numbering** matches script order (01-10)
2. **Descriptive names** clearly indicate content
3. **Logical sub-grouping** within complex analyses (e.g., diversity split into alpha/beta/methods)
4. **Script-to-directory mapping** is obvious
5. **Publication vs analysis** figures separated where appropriate

---

## Color Standardization Fixes

### Script 12: `12.nmds_permanova_cafi.R`

**Location:** Lines 132-143, 219-221

**Current code (WRONG):**
```r
# Lines 132-143: Custom color palette
if (exists("treatment_colors", inherits = FALSE)) {
  pal <- treatment_colors
  miss <- setdiff(levels(meta_df$treatment), names(pal))
  if (length(miss)) {
    pal <- c(pal, setNames(RColorBrewer::brewer.pal(max(3, length(miss)), "Dark2")[seq_along(miss)], miss))
  }
} else {
  levs <- levels(meta_df$treatment)
  pal  <- setNames(RColorBrewer::brewer.pal(max(3, length(levs)), "Dark2"), levs)
}

# Line 219: Manual color scale
scale_color_manual(values = pal) +

# Lines 220-221: Non-standard theme
theme_minimal(base_size = 14) +
```

**Fix:**
```r
# Lines 132-143: REMOVE custom palette code entirely

# Line 219: Replace with
scale_color_treatment() +

# Lines 220-221: Replace with
theme_publication() +
```

---

### Script 10: `10.coral-physio.R`

**Location:** Missing source statement, 18 instances of theme_minimal()

**Current code (WRONG):**
```r
# Missing at top of script
# (no source statement for mrb_figure_standards.R)

# Throughout script (18 instances):
theme_minimal(base_size = 14)
```

**Fix:**
```r
# Add after line 51 (after sourcing utils.R):
source("scripts/MRB/mrb_figure_standards.R")

# Replace all 18 instances of theme_minimal() with:
theme_publication()
```

**Locations to fix:**
Lines: 180, 193, 265, 278, 291, 370, 392, 426, 646, 677, 698, 716, 788, 807, 827, 914, 930, 1024

---

### Script 13: `13.SLOSS.R`

**Location:** Missing source statement

**Current code (WRONG):**
```r
# Lines 12-14
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")
# Missing: source("scripts/MRB/mrb_figure_standards.R")
```

**Fix:**
```r
# Add after line 14:
source("scripts/MRB/mrb_figure_standards.R")
```

**Note:** This script's custom SL/SS colors are acceptable as they're for categorical comparisons, not treatments.

---

## Implementation Plan

### Phase 1: Fix Scripts (Estimated: 15 minutes)

1. ✅ **Script 12:** Remove custom colors, use scale_color_treatment(), theme_publication()
2. ✅ **Script 10:** Add mrb_figure_standards.R source, replace all theme_minimal()
3. ✅ **Script 13:** Add mrb_figure_standards.R source

### Phase 2: Backup & Clean (Estimated: 5 minutes)

1. Create backup of current output directory
2. Document all current figure locations in inventory

### Phase 3: Reorganize Directory Structure (Estimated: 10 minutes)

1. Create new directory structure
2. Update all scripts to write to new locations
3. Create mapping document for old→new paths

### Phase 4: Regenerate All Figures (Estimated: 30 minutes)

1. Delete old figures
2. Run full pipeline: `source("scripts/MRB/run_mrb_pipeline.R")`
3. Verify all figures generated in correct locations

### Phase 5: Verification (Estimated: 10 minutes)

1. Count figures in each directory
2. Spot-check color consistency in key figures
3. Update documentation

**Total Estimated Time:** ~70 minutes

---

## Success Criteria

✅ All scripts source `mrb_figure_standards.R`
✅ All scripts use `scale_color_treatment()` / `scale_fill_treatment()` for treatment colors
✅ All scripts use `theme_publication()` instead of `theme_minimal()`
✅ Figures organized in intuitive, script-aligned directory structure
✅ Easy to find figures by analysis type
✅ All figures use consistent TREATMENT_COLORS:
   - Treatment 1 = #E69F00 (Orange)
   - Treatment 3 = #56B4E9 (Sky Blue)
   - Treatment 6 = #009E73 (Green)

---

## Migration Guide for Users

**Finding your figures after reorganization:**

| Old Location | New Location | Script |
|-------------|--------------|--------|
| `output/diversity/*.png` | `output/MRB/03_diversity/` | 4d |
| `output/figures/sem_*.png` | `output/MRB/07_community_coral_integration/networks/` | 8 |
| `output/MRB/figures/combined/*.png` | `output/MRB/07_community_coral_integration/` (subdivided) | 8, 10 |
| `output/MRB/figures/coral/*.png` | `output/MRB/05_coral_growth/` | 6 |
| `output/MRB/figures/coral/physio/*.png` | `output/MRB/06_coral_physiology/` | 7 |
| `output/MRB/figures/ordination/*.png` | `output/MRB/09_ordination/` | 12 |

---

## Next Steps

1. Review and approve this plan
2. Execute Phase 1: Fix scripts
3. Execute Phase 2-5: Reorganize and regenerate
4. Update all documentation and README files

---

**Document Status:** Draft for review
**Last Updated:** 2025-11-05
