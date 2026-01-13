# Codebase Architecture Audit
## Stier et al. 2025 - CAFI 136 MRB Analysis

**Audit Date:** 2026-01-11
**Status:** Cleaned and standardized for publication

---

## 1. Repository Structure

```
Stier-2025-CAFI136-MRB-AMOUNT/
├── data/
│   ├── MRB Amount/              # Raw experimental data
│   ├── processed/               # Generated data files
│   └── README.md                # Data dictionary
├── scripts/
│   └── MRB/
│       ├── 1-14.*.R             # Core analysis pipeline (numbered)
│       ├── utils.R              # Shared utility functions
│       ├── mrb_figure_standards.R  # Publication figure standards
│       ├── mrb_config.R         # Configuration (sources standards)
│       ├── agent_generated/     # Agent-generated supporting scripts
│       └── archive/             # Backup files & old documentation
├── output/
│   └── MRB/
│       ├── figures/             # Generated figures
│       ├── tables/              # Generated tables
│       └── archive/             # Archived outputs
├── agents/                      # Python agent system
├── README.md                    # Main project overview
├── CITATION.cff                 # Citation metadata
├── LICENSE                      # CC-BY 4.0
└── DATA_AVAILABILITY.md         # Data access statement
```

---

## 2. Analysis Pipeline

### Core Scripts (Execute in Order)

| # | Script | Purpose | Key Outputs |
|---|--------|---------|-------------|
| 1 | `1.libraries.R` | Load all R packages | - |
| 2 | `2.taxonomic-coverage.R` | Taxonomic resolution analysis | Taxonomic tables |
| 3 | `3.abundance.R` | CAFI abundance scaling | Abundance figures |
| 4d | `4d.diversity.R` | Alpha/beta diversity, PERMANOVA | Diversity metrics |
| 5 | `5.fishes.R` | Fish community patterns | Fish-specific analyses |
| 6 | `6.coral-growth.R` | Coral growth (3D photogrammetry) | `coral_growth.csv` |
| 7 | `7.coral-physiology.R` | Physiology + integrated performance | Performance PCA |
| 8 | `8.coral-caffi.R` | CAFI-coral community relationships | Community PCAs |
| 9 | `9.null-models.R` | Species co-occurrence null models | Null model results |
| 10 | `10.coral-physio.R` | Advanced ordination (PCA/PCoA/RDA) | Ordination figures |
| 12 | `12.nmds_permanova_cafi.R` | NMDS + PERMANOVA | NMDS plots |
| 13 | `13.SLOSS.R` | SLOSS richness resampling | SLOSS comparisons |
| 14 | `14.compile-manuscript-statistics.R` | Compile all statistics | Stats tables |

### Utility Scripts

| Script | Purpose |
|--------|---------|
| `utils.R` | Data loading, helper functions, themes, shared constants (`ALIVE_THRESH`, `strip_fe()`) |
| `mrb_figure_standards.R` | Publication colors, themes, figure specs, accent colors |
| `mrb_config.R` | Configuration wrapper (sources standards) |

---

## 3. Data Conventions

### Coral ID Format

**Standard:** POC format (e.g., `POC61`, `POC62`)
- Raw data uses `FE-POC61` format
- All scripts strip the `FE-` prefix using `strip_fe()` helper
- Processed data (`coral_growth.csv`) uses POC format

```r
# Standard pattern in all scripts:
strip_fe <- function(x) stringr::str_remove(x, "^FE-")
```

### Species Filtering (10×10 Rule)

**Standard:** Include species with:
- ≥10 corals with presence (prevalence)
- ≥10 total individuals (abundance)

```r
# Standard constants:
PREV_MIN  <- 10
ABUND_MIN <- 10
```

### Coral Survival Threshold

**Standard:** ≥80% tissue alive

```r
ALIVE_THRESH <- 0.80
```

---

## 4. Color and Theme Standards

### Treatment Colors (Colorblind-Friendly)

**Single source of truth:** `mrb_figure_standards.R`

```r
TREATMENT_COLORS <- c(
  "1" = "#E69F00",   # Orange - Single coral
  "3" = "#56B4E9",   # Sky Blue - Three corals
  "6" = "#009E73"    # Green - Six corals
)
```

### Theme Functions

| Function | Use Case |
|----------|----------|
| `theme_publication()` | Main publication figures |
| `theme_multipanel()` | Multi-panel figures |
| `theme_ordination()` | NMDS/PCA plots |
| `theme_heatmap()` | Heatmaps |

### Color Scale Functions

```r
scale_color_treatment()  # Apply treatment colors
scale_fill_treatment()   # Apply treatment fills
```

---

## 5. Issues Fixed (2026-01-11)

### Critical Fixes (First Sweep)

1. **Treatment Colors in `generate_publication_figures.R`**
   - Was: `#FFD92F`, `#8DA0CB`, `#66C2A5` (wrong colors)
   - Fixed to: `#E69F00`, `#56B4E9`, `#009E73` (standard colors)

2. **Coral ID Consistency**
   - `comprehensive_sensitivity_analysis.R` now strips `FE-` prefix to match pipeline
   - All scripts use consistent POC format

3. **Duplicate Color Definitions**
   - Removed TREATMENT_COLORS redefinition from `mrb_config.R`
   - Single source of truth in `mrb_figure_standards.R`

### Second Sweep Fixes (2026-01-11)

4. **Shared Constants in `utils.R`**
   - Added `ALIVE_THRESH <- 0.80` as shared constant
   - Added `strip_fe()` function for consistent coral ID handling
   - All scripts can now source utils.R for these shared definitions

5. **Removed Duplicate ALIVE_THRESH in `8.coral-caffi.R`**
   - Was defined twice (lines 151 and 2161)
   - Now references shared constant from utils.R

6. **Fixed Hardcoded Colors in `10.coral-physio.R`**
   - Replaced hardcoded hex colors with `ACCENT_COLOR` and `ACCENT_COLOR_ALT`
   - Added `ACCENT_COLOR` and `ACCENT_COLOR_ALT` to `mrb_figure_standards.R`

7. **Standardized Output Directory in `generate_publication_figures.R`**
   - Was: `output/figures_mrb/publication/`
   - Fixed to: `output/MRB/figures/publication/` (matches standard structure)
   - Now sources `mrb_figure_standards.R` for TREATMENT_COLORS

8. **Archived Broken Agent Scripts**
   - `analysis_test_abundance.R` → archive (missing data source)
   - `figure_total_abundance_by_treatment.R` → archive (undefined function)
   - `treatment_effect_abundance.R` → archive (circular sourcing)

### Third Sweep Fixes (2026-01-11)

9. **Removed Duplicate `strip_fe()` Definitions**
   - Removed from scripts 4d, 6, 7, 8, 12
   - All now use shared function from `utils.R`

10. **Standardized ALIVE_THRESH Naming**
    - Changed lowercase `alive_thresh` to uppercase `ALIVE_THRESH` in scripts 6, 7
    - Removed duplicate definitions from scripts 4d, 8
    - All scripts now use shared constant from `utils.R`

11. **Added Missing source() Calls**
    - Script 14: Added sourcing of `mrb_figure_standards.R` and `utils.R`

12. **Fixed Hardcoded Colors in Script 9**
    - Replaced `#D55E00`, `#0072B2` with `ACCENT_COLOR_ALT`, `ACCENT_COLOR`
    - Uses `TREATMENT_COLORS` for fill values

13. **Removed Duplicate Theme/Save Wrapper Functions**
    - Added `show_and_save()` to `utils.R` as centralized function
    - Removed duplicate `save_both`, `show_and_save` definitions from scripts 4d, 6, 7, 8, 12
    - All scripts now use shared functions from `utils.R`

### Fourth Sweep Fixes (2026-01-12)

14. **Added Missing Package to `1.libraries.R`**
    - Added `RColorBrewer` library (required by script 12)

15. **Removed Redundant Library Loading**
    - Script 8: Removed 25-line package loading block (lines 78-104)
    - Script 12: Removed 15-line suppressPackageStartupMessages block (lines 20-34)
    - All now rely on `source("scripts/MRB/1.libraries.R")`

16. **Fixed File Header in Script 3**
    - Changed from `05_MRB_community_abundance_analysis.R` to `3.abundance.R`
    - Updated "Run after" to reference current script names

17. **Added Standard Source Calls to Script 3**
    - Added `source("scripts/MRB/1.libraries.R")` and `source("scripts/MRB/utils.R")`
    - Removed custom package loading block

18. **Fixed Hardcoded Colors in `snail_corallivore_analysis.R`**
    - Replaced hardcoded `cols_trt` with `TREATMENT_COLORS` constant
    - Added standard source calls for 1.libraries.R, utils.R, mrb_figure_standards.R

19. **Standardized Agent-Generated Scripts**
    - `figureS1_physiology_reformatted.R`: Replaced library calls with standard source pattern
    - `tableS4_sensitivity_analysis.R`: Added standard source calls
    - `verify_hochberg_pvalues.R`: Added standard source calls

### Files Archived

Moved to `scripts/MRB/archive/`:
- `10.coral-physio.R.bak`, `10.coral-physio.R.bak2` (backup files)

Moved to `scripts/MRB/archive/documentation/`:
- `COLOR_STANDARDIZATION_FIXES.md`
- `FIGURE_16_IMPROVEMENTS.md`
- `FIGURE_COLOR_FIX_SUMMARY.md`
- `FIGURE_REORGANIZATION_PLAN.md`
- `HEADER_AUDIT_REPORT.md`
- `HEADER_UPDATE_SUMMARY.md`
- `PLOT_STANDARDIZATION_PLAN.md`
- `SCRIPT_REORGANIZATION_SUMMARY.md`
- `STANDARDIZATION_STATUS.md`

Moved to `scripts/MRB/agent_generated/archive/`:
- `analysis_test_abundance.R`
- `figure_total_abundance_by_treatment.R`
- `treatment_effect_abundance.R`

---

## 6. Dependency Structure

```
1.libraries.R (packages)
    ↓
utils.R + mrb_figure_standards.R (shared functions/themes)
    ↓
2.taxonomic-coverage.R → 3.abundance.R → 4d.diversity.R → 5.fishes.R
    ↓
6.coral-growth.R (exports coral_growth.csv)
    ↓
7.coral-physiology.R (imports coral_growth.csv)
    ↓
8.coral-caffi.R (community-coral relationships)
    ↓
{9.null-models.R, 10.coral-physio.R, 12.nmds_permanova_cafi.R, 13.SLOSS.R}
    ↓
14.compile-manuscript-statistics.R (final compilation)
```

---

## 7. Agent-Generated Scripts

Located in `scripts/MRB/agent_generated/`:

| Script | Purpose | Status |
|--------|---------|--------|
| `comprehensive_sensitivity_analysis.R` | 10×10 rule sensitivity testing | Active |
| `figureS1_physiology_reformatted.R` | Supplemental Figure S1 | Active |
| `tableS4_sensitivity_analysis.R` | Table S4 formatting | Active |
| `verify_hochberg_pvalues.R` | Table S2 p-value verification | Validation |

**Archived** (moved to `agent_generated/archive/`):
- `treatment_effect_abundance.R` - Circular sourcing issue
- `figure_total_abundance_by_treatment.R` - Undefined functions
- `analysis_test_abundance.R` - Missing data sources

---

## 8. Documentation Files

### Active Documentation (Root)

| File | Purpose |
|------|---------|
| `README.md` | Main project overview |
| `CITATION.cff` | Citation metadata |
| `LICENSE` | CC-BY 4.0 license |
| `DATA_AVAILABILITY.md` | Data access statement |
| `DATA_SUBSETTING.md` | Filtering methodology |
| `REPRODUCIBILITY_GUIDE.md` | Reproduction instructions |
| `CODEBASE_AUDIT.md` | This file |

### Active Documentation (Output)

| File | Purpose |
|------|---------|
| `output/MRB/MANUSCRIPT_REVISION_SUMMARY.md` | Revision tracking |
| `output/MRB/KEY_STATISTICS_FOR_MANUSCRIPT.md` | Key results summary |
| `output/MRB/tables/sensitivity_analysis_summary.md` | Sensitivity results |

---

## 9. Quality Checks

### Before Running Pipeline

1. Ensure all packages installed (see `1.libraries.R`)
2. Verify data files exist in `data/MRB Amount/`
3. Check working directory is repository root

### After Running Pipeline

1. Verify `coral_growth.csv` exported to `data/processed/`
2. Check all figures generated in `output/MRB/figures/`
3. Review statistics in `output/MRB/tables/`

### Common Issues

| Issue | Solution |
|-------|----------|
| Coral ID mismatch | Ensure `strip_fe()` applied consistently |
| Missing packages | Run `install.packages()` from `1.libraries.R` |
| Theme not found | Source `mrb_figure_standards.R` at script top |
| Color mismatch | Use `TREATMENT_COLORS` from standards file |

---

## 10. Maintenance Notes

### Adding New Scripts

1. Use `HEADER_TEMPLATE.md` as header template
2. Source `mrb_figure_standards.R` at top
3. Use standard filtering constants
4. Apply `strip_fe()` when joining with CAFI data
5. Use `theme_publication()` for figures

### Updating Colors

**Only modify:** `mrb_figure_standards.R` lines 46-50

Do NOT define colors in:
- `mrb_config.R` (sources standards)
- Individual analysis scripts
- `generate_publication_figures.R`

### Archiving Files

Move outdated files to:
- `scripts/MRB/archive/` (scripts/backups)
- `scripts/MRB/archive/documentation/` (old docs)
- `output/MRB/archive/` (old outputs)

---

*Generated by CAFI Agent System - 2026-01-11*
