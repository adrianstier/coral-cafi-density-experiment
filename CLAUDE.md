# CLAUDE.md - Claude Code Configuration for CAFI MRB Analysis

## Project Overview

This repository contains the analysis pipeline for a coral reef field experiment examining how coral habitat quantity affects coral-associated fish and invertebrate (CAFI) community assembly in Moorea, French Polynesia.

**Key Scientific Question:** How does coral density (1, 3, or 6 colonies per reef patch) affect CAFI communities and coral performance over a 2-year period?

**Publication:** Stier, A.C., Primo, A., Curtis, J.S., & Osenberg, C.W. Habitat quantity drives community assembly and feedbacks to coral performance.

**Archive:** Zenodo DOI [10.5281/zenodo.18239647](https://doi.org/10.5281/zenodo.18239647)

## Repository Structure

```
├── data/MRB Amount/           # Raw data files
├── scripts/MRB/               # R analysis scripts (run in order 1-14)
├── output/MRB/                # Analysis outputs
├── docs/                      # Documentation
└── agents/                    # Python automation system (optional)
```

## Running the Analysis

### Option 1: Simple Bash Script (Recommended for Reviewers)

```bash
./run_all.sh
```

This runs the core analysis pipeline in order.

### Option 2: Individual Scripts (For Development)

```bash
Rscript scripts/MRB/1.libraries.R
Rscript scripts/MRB/6.coral-growth.R
Rscript scripts/MRB/7.coral-physiology.R
Rscript scripts/MRB/8.coral-caffi.R
Rscript scripts/MRB/14.compile-manuscript-statistics.R
```

### Option 3: Python Agent System (Advanced)

```python
import sys
sys.path.insert(0, 'agents')
from orchestrator import Orchestrator

orch = Orchestrator()
orch.run_full_analysis()  # Runs all scripts with validation
```

## Key Analysis Scripts

| Script | Purpose | Runtime | Key Outputs |
|--------|---------|---------|-------------|
| **1.libraries.R** | Load packages | <1 min | Package verification |
| **6.coral-growth.R** | Allometric growth analysis | ~2 min | Growth figures, LMM results |
| **7.coral-physiology.R** | Coral performance metrics | ~2 min | Physiology figures, PC1 |
| **8.coral-caffi.R** | CAFI-coral relationships | ~3 min | Community PCA, correlations |
| **14.compile-manuscript-statistics.R** | Compile all stats | ~1 min | Master statistics table |

**Expected Total Runtime:** 5-10 minutes for core analysis

## Key Findings (For Verification)

When reproducing the analysis, you should obtain these key results:

1. **CAFI Community Assembly:**
   - Abundance increases 5× with coral density (1 → 6 colonies)
   - Species richness doubles with coral density
   - Community composition shifts significantly (PERMANOVA p < 0.001)

2. **Coral Growth:**
   - No treatment effect on size-corrected growth (χ² = 2.64, p = 0.267)
   - Unified allometric exponent b = 0.6986
   - Interaction model shows treatment-specific slopes differ (p = 0.039)

3. **Coral Physiology:**
   - Carbohydrate content decreases with density (χ² = 10.0, p = 0.007)
   - Integrated performance (PC1) declines with density (χ² = 8.11, p = 0.017)

4. **CAFI-Coral Feedbacks:**
   - CAFI community composition predicts coral performance
   - Multiple transformations show robust relationship (p < 0.001 to p = 0.042)

## Data Filtering Decisions

**Important:** The analysis applies these filters documented in `docs/DATA_SUBSETTING.md`:

1. **Coral Growth:** Only corals with ≥80% tissue alive retained (n = 44 of 54)
   - Rationale: Partial mortality alters growth patterns

2. **CAFI Species:** 10×10 rule (≥10 corals AND ≥10 total individuals)
   - Analysis-level filtering: 38 species retained
   - Display-level filtering: Top 15-20 species for figures

## Documentation Files

| File | Purpose |
|------|---------|
| [README.md](README.md) | Project overview and quick start |
| [docs/REPRODUCIBILITY_GUIDE.md](docs/REPRODUCIBILITY_GUIDE.md) | Step-by-step reproduction instructions |
| [docs/DATA_AVAILABILITY.md](docs/DATA_AVAILABILITY.md) | Data access and citation |
| [docs/DATA_SUBSETTING.md](docs/DATA_SUBSETTING.md) | Filtering decisions and justification |
| [docs/MANUSCRIPT_TEXT_TEMPLATES.md](docs/MANUSCRIPT_TEXT_TEMPLATES.md) | Ready-to-use Methods/Results text |
| [docs/ZERO_INFLATION_PCA_ASSESSMENT.md](docs/ZERO_INFLATION_PCA_ASSESSMENT.md) | PCA validity with sparse data |

## R Package Dependencies

Core packages loaded by `scripts/MRB/1.libraries.R`:

- **Statistical Modeling:** lmerTest, emmeans, nlme
- **Community Ecology:** vegan, BiodiversityR
- **Data Manipulation:** tidyverse (dplyr, tidyr, ggplot2)
- **Utilities:** here, cli, patchwork

**R Version:** 4.3.x or higher recommended

## Common Tasks for Reviewers

### Verify Key Statistical Results

```bash
# Run growth analysis and check output
Rscript scripts/MRB/6.coral-growth.R

# Check for expected p-value (p = 0.267 for treatment effect)
grep "treatment effect" output/MRB/tables/growth_*.csv
```

### Regenerate Main Figures

```bash
Rscript scripts/MRB/generate_publication_figures.R
```

Outputs will be in `output/MRB/figures/publication-figures/`

### Check Data Integrity

```bash
# Verify checksums (if available)
md5sum -c data/checksums.txt

# Or use Python validation agent
python -m agents --validate
```

## File Naming Conventions

- **Scripts:** Numbered 1-14 in execution order
- **Figures:** `figure{N}_{descriptive_name}.{png|pdf}`
- **Tables:** `{analysis}_{transformation}.csv`
- **Data:** Original field names preserved (e.g., "MRB Amount")

## Notes for Claude Code Users

When asked to:

1. **"Run the analysis"** → Use `./run_all.sh` or Python orchestrator
2. **"Generate a figure"** → Use FigureAgent or relevant numbered script
3. **"Check statistical results"** → Read from `output/MRB/tables/`
4. **"Validate data"** → Use ValidationAgent or check `docs/REPRODUCIBILITY_GUIDE.md`

## Reproducibility Verification

Expected output file counts:

- **Figures:** ~80 total (14 in publication-figures/)
- **Tables:** ~90 CSV/HTML files
- **Session Info:** `output/MRB/objects/sessionInfo_*.txt`

Expected runtime on modern laptop:
- **Full analysis:** 5-10 minutes
- **Core scripts only:** 3-5 minutes

## Contact

- **Corresponding Author:** Adrian Stier (astier@ucsb.edu)
- **Code Issues:** https://github.com/adrianstier/coral-cafi-density-experiment/issues
- **Zenodo Archive:** https://doi.org/10.5281/zenodo.18239647

## License

- **Code:** MIT License
- **Data:** CC-BY 4.0

---

**Last Updated:** January 14, 2026
**Repository Version:** v1.0.0
