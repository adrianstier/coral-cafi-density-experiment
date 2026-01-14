# Contributing to CAFI Analysis Repository

Thank you for your interest in this research! This repository contains the analysis pipeline for:

**Stier, A.C., Primo, A., Curtis, J.S., & Osenberg, C.W.** Habitat quantity drives community assembly and feedbacks to coral performance.

## Ways to Contribute

### 1. Report Issues

If you find errors, bugs, or have questions about the analysis:

1. Check existing [Issues](https://github.com/adrianstier/coral-cafi-density-experiment/issues)
2. Create a new issue with:
   - Clear title describing the problem
   - Steps to reproduce (if applicable)
   - Expected vs. actual behavior
   - Your R version and operating system
   - Error messages (copy full error text)

### 2. Suggest Improvements

We welcome suggestions for:
- Code optimization or clarity
- Additional statistical tests or sensitivity analyses
- Figure improvements
- Documentation enhancements

Create an issue with the `enhancement` label.

### 3. Reproduce the Analysis

The best contribution is verification! Try to reproduce our results and report:
- ✓ Successful reproduction (with R/package versions)
- ✗ Failed reproduction (with error details)
- Any deviations from expected results

## Development Workflow

### Setting Up Your Environment

1. **Fork and Clone:**
   ```bash
   git clone https://github.com/YOUR_USERNAME/coral-cafi-density-experiment.git
   cd coral-cafi-density-experiment
   ```

2. **Install Dependencies:**
   ```r
   source("scripts/MRB/1.libraries.R")  # Installs required packages
   ```

3. **Verify Data Integrity:**
   ```bash
   md5sum -c data/checksums.txt  # Linux/WSL
   # or
   md5 -c data/checksums.txt     # macOS
   ```

### Code Style Guidelines

Follow these conventions when modifying scripts:

#### **1. Script Headers (Required)**

```r
# ==============================================================================
# Script: [NUMBER].[descriptive-name].R
# Purpose: [One-line description]
# Author: [Your Name]
# Date: [YYYY-MM-DD]
# ==============================================================================
#
# Inputs:
#   - data/MRB Amount/[specific files used]
#
# Outputs:
#   - output/MRB/figures/[specific figures]
#   - output/MRB/tables/[specific tables]
#
# Dependencies:
#   - Packages: tidyverse, lmerTest, vegan, etc.
#   - Scripts: Must run AFTER script X
#
# ==============================================================================
```

#### **2. Section Markers**

```r
# ---- 1. Load Data ----
# ---- 2. Data Processing ----
# ---- 3. Statistical Analysis ----
# ---- 4. Figures ----
# ---- 5. Save Outputs ----
```

#### **3. Inline Comments**

- Explain WHY, not WHAT (code shows what)
- Document statistical choices:
  ```r
  # Use Gower distance instead of Jaccard to handle mixed data types
  dist_matrix <- vegdist(comm_mat, method = "gower")
  ```

#### **4. Function Documentation**

Use roxygen2 style for functions in `utils.R`:

```r
#' Calculate treatment-specific color palette
#'
#' @param treatment Character vector of treatment levels
#' @return Named vector of hex colors
#' @export
#' @examples
#' get_treatment_colors(c("1", "3", "6"))
```

#### **5. Naming Conventions**

- **Variables:** snake_case (`coral_density`, `species_richness`)
- **Functions:** snake_case (`calculate_diversity`, `plot_abundance`)
- **Files:** Numbered with hyphens (`6.coral-growth.R`)
- **Figures:** Descriptive with underscores (`figure1_abundance_by_treatment.png`)

### Testing Your Changes

Before submitting:

1. **Run Modified Script:**
   ```bash
   Rscript scripts/MRB/YOUR_SCRIPT.R
   ```

2. **Run Full Pipeline:**
   ```bash
   ./run_all.sh
   ```

3. **Check Outputs:**
   - Figures render correctly
   - Tables contain expected columns
   - No errors or warnings (except expected ones)

4. **Save Session Info:**
   ```r
   sink("output/MRB/objects/sessionInfo_YOUR_SCRIPT.txt")
   sessionInfo()
   sink()
   ```

### Submitting Changes

1. **Create a branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make changes and commit:**
   ```bash
   git add your/changed/files
   git commit -m "Descriptive commit message"
   ```

3. **Push to your fork:**
   ```bash
   git push origin feature/your-feature-name
   ```

4. **Open a Pull Request:**
   - Describe what changed and why
   - Reference any related issues
   - Include output comparisons if results changed

## Adding New Analyses

### Creating a New Script

1. **Number it appropriately** (15, 16, etc.) based on execution order
2. **Use the header template** (see Code Style section)
3. **Source utilities:**
   ```r
   library(here)
   source(here("scripts/MRB/utils.R"))
   source(here("scripts/MRB/mrb_figure_standards.R"))
   ```

4. **Save outputs to organized directories:**
   ```r
   save_both(plot, here("output/MRB/figures/NEW_CATEGORY/figure_name"))
   write_csv(results, here("output/MRB/tables/analysis_results.csv"))
   ```

5. **Update Documentation:**
   - Add to `run_all.sh` if part of main pipeline
   - Document in `docs/REPRODUCIBILITY_GUIDE.md`
   - Update `README.md` if major addition

### Using the Python Agent System

To generate R code programmatically:

```python
from agents import StatsAgent

agent = StatsAgent()
script, path = agent.generate_lmm_analysis_script(
    analysis_name="new_analysis",
    response_var="your_variable",
    fixed_effects=["treatment", "reef_type"],
    random_effect="reef",
    run_script=True  # Execute immediately
)
```

Generated scripts go to `scripts/MRB/agent_generated/`

## Statistical Best Practices

When adding or modifying statistical analyses:

1. **Document Assumptions:**
   - Normality checks for parametric tests
   - Variance homogeneity
   - Independence of observations

2. **Use Mixed Models When Appropriate:**
   ```r
   # Account for reef-level clustering
   model <- lmer(response ~ treatment + (1|reef), data = data)
   ```

3. **Apply Multiple Testing Corrections:**
   ```r
   # Benjamini-Hochberg for family-wise error
   p_adjusted <- p.adjust(p_values, method = "BH")
   ```

4. **Report Effect Sizes:**
   - Cohen's d for differences
   - R² for model fit
   - Correlation coefficients

5. **Include Sensitivity Analyses** when making methodological choices

## Data Management

### Do NOT:
- Modify raw data files in `data/MRB Amount/`
- Commit large binary files (>.5 MB)
- Remove archive folders without justification

### Do:
- Save processed data to `data/processed/`
- Document any data transformations
- Use `here::here()` for all file paths
- Update `data/checksums.txt` if adding data

## Figure Standards

All publication figures should follow `mrb_figure_standards.R`:

```r
library(here)
source(here("scripts/MRB/mrb_figure_standards.R"))

plot <- ggplot(data, aes(x, y, color = treatment)) +
  geom_point() +
  scale_color_manual(values = treatment_colors) +  # Consistent colors
  theme_publication() +                             # Standard theme
  labs(title = "Your Title", x = "X Label", y = "Y Label")

save_both(plot, here("output/MRB/figures/category/figure_name"))  # Saves PNG + PDF
```

**Required figure specifications:**
- Resolution: 300 DPI (PNG), vector (PDF)
- Font: Arial or Helvetica, 10-12pt
- Colors: Use `treatment_colors` palette (orange/blue/green)
- Size: 6-8 inches wide for single-column, 12-14 inches for two-column

## Questions?

- **General questions:** Open an issue with `question` label
- **Code questions:** Check `docs/REPRODUCIBILITY_GUIDE.md` first
- **Statistical questions:** Email corresponding author (astier@ucsb.edu)

## Code of Conduct

This project follows standard open science principles:

- Be respectful and constructive in discussions
- Credit original authors when reusing code
- Report bugs and issues transparently
- Share improvements with the community

## License

By contributing, you agree that your contributions will be licensed under the same terms as the project:

- **Code:** MIT License
- **Data:** CC-BY 4.0

---

Thank you for helping improve open science in coral reef ecology!
