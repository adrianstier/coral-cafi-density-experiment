# Habitat Quantity Drives Community Assembly and Feedbacks to Coral Performance

**Authors:** Adrian C. Stier, Alexander Primo, Joseph S. Curtis, Craig W. Osenberg

**Repository:** Analysis code and data for coral density experiment examining CAFI community assembly and feedbacks to coral performance

**Status:** Ready for code review and reproducibility assessment

---

## Overview

This repository contains all data, code, and outputs for analyzing a field experiment that manipulated coral density (1, 3, or 6 *Pocillopora* colonies per reef) to test how habitat quantity affects:

1. Coral-associated fish and invertebrate (CAFI) community assembly
2. Feedbacks between CAFI communities and coral performance (growth and physiology)

**Key Findings:**
- CAFI abundance increased 5× and species richness doubled with increasing coral density
- CAFI community composition shifted with coral density
- Coral performance (integrated growth + physiology) declined with increasing coral density
- CAFI community composition predicted coral performance
- **Growth finding:** No treatment effect on coral growth after allometric size-correction (*p* = 0.267)

---

## Repository Structure

```
.
├── README.md                           # This file
├── REPRODUCIBILITY_GUIDE.md            # Step-by-step reproduction instructions
├── data/                               # Raw data files
│   └── MRB/
│       ├── coral/                      # Coral growth measurements (3D photogrammetry)
│       ├── cafi/                       # CAFI community data
│       └── physiology/                 # Coral tissue trait data
├── scripts/                            # Analysis scripts (run in order)
│   └── MRB/
│       ├── 1.data-organization.R       # Data import and organization
│       ├── 2.exploratory-figures.R     # Exploratory data visualization
│       ├── 3.cafi-community.R          # CAFI community analysis
│       ├── 4d.cafi-diversity.R         # CAFI diversity metrics
│       ├── 5.cafi-composition.R        # CAFI community composition (NMDS, PERMANOVA)
│       ├── 6.coral-growth.R            # Coral growth analysis (MAIN SCRIPT)
│       ├── 7.coral-physiology.R        # Coral physiological metrics
│       ├── 8.cafi-coral-community.R    # CAFI-coral performance relationships
│       ├── 12.publication-figures.R    # Generate publication-quality figures
│       └── 14.compile-manuscript-statistics.R  # Compile all stats for manuscript
├── output/                             # Analysis outputs
│   └── MRB/
│       ├── figures/                    # All generated figures
│       │   └── publication-figures/    # Final publication figures
│       ├── tables/                     # Statistical summary tables
│       ├── MANUSCRIPT_STATS_TABLE.csv  # Complete statistical results
│       ├── KEY_STATISTICS_FOR_MANUSCRIPT.md  # Quick reference for key results
│       ├── MANUSCRIPT_UPDATE_SUMMARY.md      # Summary of analysis updates
│       └── archive/                    # Archived development files
└── .gitignore                          # Git ignore rules
```

---

## Quick Start

### Prerequisites

**R version:** 4.3.x or higher

**Required R packages:**
```r
install.packages(c(
  "tidyverse",    # Data manipulation and visualization
  "here",         # File path management
  "lme4",         # Linear mixed-effects models
  "lmerTest",     # p-values for lme4
  "emmeans",      # Post-hoc comparisons
  "vegan",        # Community ecology analyses
  "gt",           # Table formatting
  "patchwork",    # Figure composition
  "cli"           # Console output formatting
))
```

### Running the Analysis

**Option 1: Run everything (recommended for full reproducibility)**
```bash
# From repository root
Rscript scripts/MRB/1.data-organization.R
Rscript scripts/MRB/2.exploratory-figures.R
Rscript scripts/MRB/3.cafi-community.R
Rscript scripts/MRB/4d.cafi-diversity.R
Rscript scripts/MRB/5.cafi-composition.R
Rscript scripts/MRB/6.coral-growth.R              # KEY SCRIPT
Rscript scripts/MRB/7.coral-physiology.R
Rscript scripts/MRB/8.cafi-coral-community.R
Rscript scripts/MRB/12.publication-figures.R
Rscript scripts/MRB/14.compile-manuscript-statistics.R
```

**Option 2: Run key manuscript analyses only**
```bash
Rscript scripts/MRB/6.coral-growth.R              # Coral growth (allometric models)
Rscript scripts/MRB/7.coral-physiology.R          # Physiology + integrated performance
Rscript scripts/MRB/8.cafi-coral-community.R      # CAFI-coral relationships
Rscript scripts/MRB/14.compile-manuscript-statistics.R  # Compile all stats
```

**Expected runtime:** ~5-10 minutes for all scripts

---

## Key Analyses and Files

### 1. Coral Growth Analysis (Script 6)

**File:** `scripts/MRB/6.coral-growth.R`

**What it does:**
- Tests for treatment-specific allometric growth relationships (interaction model)
- Estimates unified allometric exponent for size-correction (*b* = 0.6986)
- Tests treatment effect on size-corrected growth
- Generates 3-panel growth figure

**Key outputs:**
- `output/MRB/figures/coral/ANCOVA_Init_vs_Final_Volume_by_Treatment.png`
- `output/MRB/figures/coral/ANCOVA_TreatmentSpecific_Slopes.png`
- `output/MRB/figures/coral/SizeCorrected_Volume_Growth_by_Treatment.png`

**Key results:**
- Interaction: χ² = 6.48, *p* = 0.039 (treatment-specific slopes differ)
- Treatment-specific slopes: *b* = 0.78 (1 colony), -0.10 (3 colonies), 1.00 (6 colonies)
- Post-hoc comparisons: all *p* > 0.07 (no pairwise differences)
- **Main finding:** No treatment effect on size-corrected growth (χ² = 2.64, *p* = 0.267)

### 2. Coral Physiology (Script 7)

**File:** `scripts/MRB/7.coral-physiology.R`

**What it does:**
- Analyzes carbohydrate, protein, zooxanthellae, AFDW by treatment
- Creates integrated performance PC1 (physiology + growth)
- Multiple testing correction (Benjamini-Hochberg)

**Key results:**
- Carbohydrate: χ² = 10.0, *p* = 0.007 (BH-adjusted *p* = 0.039)
- PC1 (performance): χ² = 8.11, *p* = 0.017
- Individual metrics (protein, zooxanthellae, AFDW): NS

### 3. CAFI-Coral Relationships (Script 8)

**File:** `scripts/MRB/8.cafi-coral-community.R`

**What it does:**
- Tests whether CAFI community composition predicts coral performance
- Uses multiple data transformations for robustness
- Creates CAFI community PCA and coral performance PCA

**Key results:**
- SQRT_CS transformation: *p* = 0.00085
- HELLINGER transformation: *p* = 0.014
- SQRT transformation: *p* = 0.042
- **Robust, highly significant relationship**

### 4. Statistical Compilation (Script 14)

**File:** `scripts/MRB/14.compile-manuscript-statistics.R`

**What it does:**
- Compiles all statistical tests from Scripts 6-8
- Creates human-readable HTML table
- Generates CSV for easy import

**Key outputs:**
- `output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.csv`
- `output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.html`
- `output/MRB/MANUSCRIPT_STATS_TABLE.csv`

---

## Data Files

### Coral Growth Data
**Location:** `data/MRB/coral/`

**Files:**
- Individual coral volume/surface area measurements from 3D photogrammetry
- Initial (2019) and final (2021) measurements
- Tissue mortality assessments

**Sample size:** *n* = 44 corals (≥80% tissue alive)

### CAFI Community Data
**Location:** `data/MRB/cafi/`

**Files:**
- Species-level abundance data for all CAFI taxa
- Taxonomic identifications

**Sample size:** 23 experimental reefs (18 one-coral, 6 three-coral, 3 six-coral)

### Coral Physiology Data
**Location:** `data/MRB/physiology/`

**Files:**
- Carbohydrate content (mg/cm²)
- Protein content (mg/cm²)
- Zooxanthellae density (cells/cm²)
- Ash-free dry mass (AFDM, mg/cm²)

---

## Key Statistical Approaches

### Allometric Growth Model

We used a **unified allometric model** approach (Osenberg, pers. comm.) to handle size-dependent growth:

1. **Test for interaction:** `log(V_final) ~ log(V_initial) × treatment + (1|reef)`
   - Significant interaction (*p* = 0.039) indicated treatment-specific slopes
   - However, no pairwise differences after Tukey adjustment (all *p* > 0.07)

2. **Estimate unified exponent:** `log(V_final) ~ log(V_initial) + treatment + (1|reef)`
   - Unified *b* = 0.6986 (more precise than treatment-specific estimates)
   - Used for all size-corrections: `growth_vol_b = V_final / V_initial^0.6986`

3. **Test treatment effect:** Size-corrected growth ~ treatment + (1|reef)
   - Result: *p* = 0.267 (NS)

**Rationale:** Balances evidence for treatment-specific slopes against need for robust, comparable growth metric given high uncertainty in some slope estimates (SE = 0.38 for 3-colony treatment).

### Linear Mixed-Effects Models

All models used `lme4::lmer()` with:
- **Random effects:** Reef as random intercept (accounts for spatial structure)
- **Estimation:** REML for final models, ML for likelihood ratio tests
- **Inference:** Type III Wald χ² tests (lmerTest package)

### Multiple Testing Correction

Benjamini-Hochberg procedure for physiology metrics (Script 7)

---

## Reproducibility Notes

### Data Filtering

**Coral growth:** Only colonies with ≥80% tissue alive retained (*n* = 44 of 54)

**Rationale:** Partially dead colonies exhibit altered growth patterns

### Software Versions

- R 4.3.x
- lme4 1.1-35.x
- lmerTest 3.1-3
- emmeans 1.10.0
- vegan 2.6-x
- tidyverse 2.0.0

### Random Effects Structure

All models include `(1|reef)` random intercept to account for:
- Spatial clustering (reefs arranged in 9×3 grid)
- Shared environmental conditions within reef
- Non-independence of colonies on same reef

### Known Issues

None currently. All scripts run successfully with provided data.

---

## Output Files for Manuscript

### Main Statistics
- `output/MRB/MANUSCRIPT_STATS_TABLE.csv` - Complete statistical results table
- `output/MRB/KEY_STATISTICS_FOR_MANUSCRIPT.md` - Quick reference with copy-paste ready text

### Figures
- `output/MRB/figures/publication-figures/` - All publication-quality figures
  - Growth figures (3-panel: interaction, slopes, size-corrected)
  - Physiology figures
  - CAFI-coral relationship figures

### Tables
- `output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.csv` - All statistical tests
- `output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.html` - Formatted HTML table

---

## Citation

If you use this code or data, please cite:

Stier, A.C., Primo, A., Curtis, J.S., & Osenberg, C.W. (2025). Habitat quantity drives community assembly and feedbacks to coral performance in reef systems. *[Journal]*, *[Volume]*(*Issue*), *[Pages]*.

---

## Contact

**Corresponding Author:** Adrian Stier (astier@ucsb.edu)

**Code Questions:** Please open an issue on GitHub or email corresponding author

---

## License

This work is licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) (Creative Commons Attribution 4.0 International). See [LICENSE](LICENSE) for details.

---

## Acknowledgments

- Craig Osenberg for statistical consultation on unified allometric model approach
- Moorea Coral Reef LTER for field site access and support

---

**Last Updated:** November 15, 2025
**Status:** Ready for peer review and reproducibility assessment
