# Repository Publication Checklist

**Purpose:** Ensure repository is ready for public sharing with manuscript submission

**Date:** 2026-01-13

---

## ✅ Completed Items

### Documentation
- [x] **README.md** - Comprehensive overview with quick start guide
- [x] **LICENSE** - CC-BY 4.0 license file
- [x] **CITATION.cff** - Machine-readable citation information
- [x] **docs/DATA_AVAILABILITY.md** - Data availability statement for manuscript
- [x] **docs/DATA_SUBSETTING.md** - Transparent documentation of filtering decisions
- [x] **docs/CODEBASE_AUDIT.md** - Complete audit of code organization

### Analysis Documentation
- [x] **docs/ZERO_INFLATION_PCA_ASSESSMENT.md** - Assessment of PCA validity with zero-inflated data
- [x] **docs/FIGURE_GUIDE_FOR_PUBLICATION.md** - Guide to all publication figures

### Code Organization
- [x] All scripts have clear headers with purpose and dependencies
- [x] Scripts numbered in execution order (1, 2, 3, 4d, 5, 6, 7, 8, 9, 10, 12, 14)
- [x] Agent-generated scripts archived in `scripts/MRB/agent_generated/archive/`
- [x] Development documentation archived in `scripts/MRB/archive/documentation/`
- [x] Backup script files archived in `scripts/MRB/archive/`

### Data Files
- [x] Raw data in `data/MRB Amount/`
- [x] Data README explaining file structure
- [x] Processed data in `data/processed/`
- [x] No sensitive or personal information in data files
- [x] All data files documented with metadata

### Output Files
- [x] Publication figures archived with version history
- [x] Old figures in `output/MRB/archive/old-publication-figures/`
- [x] Development figures in `output/MRB/archive/publication-figures-raw-jan12/`
- [x] Statistical tables organized in `output/MRB/tables/`
- [x] Comprehensive sensitivity analyses included

### Reproducibility
- [x] All required R packages listed in README
- [x] `here::here()` used for all file paths (no hardcoded paths)
- [x] Session info saved in `output/MRB/objects/sessionInfo_*.txt`
- [x] Random seeds set where applicable
- [x] All scripts run successfully from project root

### Version Control
- [x] `.gitignore` properly configured
- [x] `.DS_Store` files removed
- [x] No large binary files tracked
- [x] Commit history preserved with descriptive messages
- [x] Repository pushed to GitHub

---

## 📋 Pre-Submission Verification

### Test Reproducibility (Run These Commands)

```bash
# From repository root
cd /Users/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT

# 1. Check R version
R --version

# 2. Install packages (if needed)
Rscript -e "install.packages(c('tidyverse', 'here', 'lme4', 'lmerTest', 'emmeans', 'vegan', 'gt', 'patchwork', 'cli'))"

# 3. Run key scripts to verify reproducibility
Rscript scripts/MRB/6.coral-growth.R
Rscript scripts/MRB/7.coral-physiology.R
Rscript scripts/MRB/8.coral-caffi.R

# 4. Generate publication figures
Rscript scripts/MRB/generate_publication_figures.R

# 5. Compile manuscript statistics
Rscript scripts/MRB/14.compile-manuscript-statistics.R

# 6. Run zero-inflation analysis
Rscript scripts/MRB/agent_generated/analyze_zero_inflation_pca.R

# 7. Verify all outputs generated
ls -lh output/MRB/tables/ | wc -l
ls -lh output/MRB/figures/ | wc -l
```

**Expected Results:**
- All scripts should run without errors
- Output files should be generated in `output/MRB/`
- Total runtime: ~5-10 minutes

---

## 🔍 Quality Assurance

### Code Quality
- [x] All scripts follow consistent naming conventions
- [x] Functions documented with clear comments
- [x] No hardcoded file paths (all use `here::here()`)
- [x] No personal/sensitive information in code
- [x] Console output informative with `cli` package

### Statistical Transparency
- [x] All filtering decisions documented
- [x] Sample sizes reported at each step
- [x] Multiple testing corrections applied where appropriate
- [x] Sensitivity analyses for key findings
- [x] P-values and effect sizes reported consistently

### Figure Quality
- [x] Publication figures at 300 DPI
- [x] Both PDF and PNG versions available
- [x] Consistent color scheme across figures
- [x] Font sizes readable (≥8pt)
- [x] Figure legends clear and complete

---

## 📊 Key Files for Reviewers

### Essential Files to Highlight
1. **README.md** - Start here for overview
2. **DATA_AVAILABILITY.md** - Data access and sharing policy
3. **DATA_SUBSETTING.md** - Filtering decisions and justification
4. **ZERO_INFLATION_PCA_ASSESSMENT.md** - PCA validity assessment

### Key Analysis Scripts
1. **scripts/MRB/6.coral-growth.R** - Coral growth analysis
2. **scripts/MRB/7.coral-physiology.R** - Physiology and integrated performance
3. **scripts/MRB/8.coral-caffi.R** - CAFI-coral relationships
4. **scripts/MRB/generate_publication_figures.R** - Publication figures

### Key Output Files
1. **output/MRB/tables/sensitivity_comprehensive_summary.md** - All sensitivity analyses
2. **output/MRB/tables/zero_inflation_species_sparsity.csv** - Species sparsity metrics
3. **output/MRB/figures/zero_inflation_diagnostic.png** - PCA diagnostic figure

---

## 🚀 Publication Actions

### Before Submission
- [ ] Update README.md with final manuscript citation placeholder
- [ ] Verify GitHub repository is public (or will be made public upon acceptance)
- [ ] Create Zenodo DOI for repository snapshot
- [ ] Add DOI badge to README.md
- [ ] Update "Last Updated" date in README.md

### In Manuscript Methods
Include the following text:

> **Data and Code Availability**
>
> All data and analysis code are available at: https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT
>
> A permanent archived version is available at Zenodo: [DOI will be inserted here]
>
> The repository includes:
> - Raw data files (coral growth, physiology, CAFI community surveys)
> - Complete analysis pipeline (R scripts)
> - Publication-quality figures
> - Comprehensive documentation of all filtering and analytical decisions
> - Sensitivity analyses for all major findings
>
> See docs/DATA_AVAILABILITY.md and docs/DATA_SUBSETTING.md for detailed information on data structure, filtering decisions, and reproducibility.

---

## 📝 Recommended Supplementary Materials

### Table S1: Data Filtering Summary
- Location: Already documented in DATA_SUBSETTING.md
- Consider converting to formatted table for supplement

### Table S2: R Package Versions
- Location: `output/MRB/objects/sessionInfo_*.txt`
- Consider creating formatted table with key packages

### Table S3: Sensitivity Analyses
- Location: `output/MRB/tables/sensitivity_comprehensive_summary.csv`
- Already comprehensive, can use as-is

### Table S4: PCA Sensitivity to Transformation
- Location: `output/MRB/tables/tableS4_pca_sensitivity.csv`
- Ready for inclusion

### Figure S1: Zero-Inflation Diagnostic
- Location: `output/MRB/figures/zero_inflation_diagnostic.pdf`
- Shows PCA validity assessment

### Figure S2: Species Prevalence Distribution
- Could extract Panel A from diagnostic figure

---

## 🔒 Final Checks Before Going Public

### Privacy and Ethics
- [x] No personal identifiable information in data
- [x] No sensitive location data (general "Moorea" is fine)
- [x] No email addresses or phone numbers (except in citation)
- [x] Institutional approval for data sharing confirmed

### Legal and Licensing
- [x] CC-BY 4.0 license clearly stated
- [x] All authors agree to open-source release
- [x] No copyright conflicts with journal policies
- [x] Funding agency requirements met for data sharing

### Technical
- [x] Repository size reasonable (<100 MB recommended)
- [x] No large binary files tracked in git
- [x] All links in README.md functional
- [x] GitHub Actions CI/CD not needed (optional)

---

## ✨ Optional Enhancements

### Could Add (but not required):
- [ ] GitHub Actions workflow to test reproducibility
- [ ] Docker container with R environment
- [ ] Binder/RStudio Cloud link for browser-based analysis
- [ ] Video walkthrough of analysis pipeline
- [ ] Interactive Shiny app for exploring results

### For Enhanced Visibility:
- [ ] Tweet thread summarizing findings with link to repo
- [ ] Blog post on lab website
- [ ] Add to relevant awesome-lists (e.g., awesome-ecology)
- [ ] Submit to Journal of Open Source Software (JOSS)

---

## 📧 Repository URL for Manuscript

**Include this in manuscript:**

```
Code and data available at: https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT
```

**Alternative (if using Zenodo):**

```
Code and data available at: https://doi.org/10.5281/zenodo.XXXXXXX
```

---

## ✅ Final Status

**Repository Status:** ✅ READY FOR PUBLIC RELEASE

**Last Checked:** 2026-01-13

**Checked By:** Adrian Stier

**Confidence Level:** HIGH - All essential documentation complete, code tested, outputs verified

---

## Notes for Future Updates

### If Reviewers Request Changes:
1. Create new branch for revisions: `git checkout -b revision-v1`
2. Make changes and commit with clear messages
3. Merge back to main after acceptance
4. Tag release version: `git tag -a v1.0 -m "Initial publication version"`

### After Acceptance:
1. Update README with final citation
2. Create GitHub release with DOI
3. Archive on Zenodo
4. Update manuscript with final DOI
5. Tweet/share widely!

---

**Ready to share! 🎉**
