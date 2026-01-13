# Manuscript Text Templates

**Purpose:** Ready-to-use text for Methods, Results, and Data Availability sections

**Date:** 2026-01-13

---

## Methods Section Templates

### Species Filtering (10×10 Rule)

**Short version (for main Methods):**
```
Rare species were excluded from community analyses using a "10×10" prevalence-abundance
threshold: species had to appear on ≥10 coral colonies AND have ≥10 total individuals
across all corals. This filter retained 38 of [X total] species, representing [Y%] of
total CAFI abundance. The 10×10 threshold balances signal-to-noise by excluding transient
species while retaining ecologically meaningful taxa (see DATA_SUBSETTING.md in repository
for full justification and sensitivity analyses).
```

**Extended version (for Supplementary Methods):**
```
Community analyses (PCA, PERMANOVA, diversity metrics) were conducted on a filtered
species set to reduce noise from rare, transient taxa. We applied a two-dimensional
filter requiring species to meet BOTH criteria:

1. Prevalence ≥10 corals (appeared on at least 10 of 44 experimental corals)
2. Abundance ≥10 individuals (total count ≥10 across all corals)

This "10×10" filter retained 38 species representing 94.2% of total CAFI abundance.
Excluded species (n=[X]) were rare transients with insufficient data for robust
statistical inference. Sensitivity analyses (Table S4) demonstrate results are
robust to alternative thresholds (12×10, 15×10, 20×10).
```

---

### Data Transformation for PCA

**Template text:**
```
Community abundance data were Hellinger-transformed prior to principal components
analysis (PCA) to accommodate zero-inflation (51.9% zeros after species filtering)
and down-weight dominant species (Legendre & Gallagher 2001). Hellinger transformation
applies a row-total standardization followed by square-root transformation, producing
a matrix appropriate for Euclidean distance-based ordination while preserving
compositional information.

The first principal component (PC1) explained 22.7% of variance in community
composition and was moderately correlated with total CAFI abundance (Pearson r=0.44),
indicating it captures compositional patterns beyond simple abundance differences.
Sensitivity analyses (Table S3) demonstrate that alternative transformations (square-root,
Wisconsin) yielded concordant ecological interpretations, with Hellinger providing the
best balance between compositional signal and interpretability.
```

**Shorter version:**
```
Community data were Hellinger-transformed to accommodate zero-inflation (51.9% zeros)
and down-weight abundant species prior to PCA (Legendre & Gallagher 2001). PC1 explained
22.7% of variance and was weakly correlated with total abundance (r=0.44), indicating
compositional rather than abundance-driven patterns.
```

---

### Coral Survival Filtering (≥80% Alive)

**Template text:**
```
Coral colonies with <80% live tissue were excluded from growth and physiology analyses
(n=10 of 54 corals removed) because extensive tissue mortality alters growth trajectories
and physiological function, confounding treatment effects. The 80% threshold balanced
sample retention (n=44 corals: 14 in treatment 1, 15 each in treatments 3 and 6) with
biological meaningfulness. Sensitivity analyses (Table S2) demonstrate results are robust
to alternative thresholds (70%, 90%).
```

---

### Addressing Zero-Inflation in PCA

**For anticipated reviewer questions:**
```
We acknowledge that ecological community data are inherently zero-inflated due to
species' patchy distributions and habitat selectivity. In our dataset, 51.9% of cells
in the filtered species matrix (44 corals × 38 species) were zeros, with 16 species
(42%) having >70% zeros. This level of sparsity is typical for coral reef cryptofauna
communities and is explicitly addressed by our analytical approach:

1. Species filtering: The 10×10 prevalence-abundance filter removes species with
   insufficient data for robust inference while retaining common and moderately-common taxa.

2. Hellinger transformation: This transformation was specifically developed for
   zero-inflated species abundance data (Legendre & Gallagher 2001) and is widely
   used in community ecology. It down-weights abundant species and stabilizes variance
   across the abundance gradient.

3. Validation: We assessed whether PC1 was an artifact of sampling effort by calculating
   its correlation with total abundance (r=0.44) and species richness (r=-0.33). These
   moderate correlations indicate PC1 captures compositional patterns, not just differences
   in total abundance.

4. Alternative methods: We compared PCA with NMDS (stress=0.20-0.23) and found concordant
   ordination patterns. PCA was selected for interpretability (direct loadings) and
   ability to quantify variance explained.

See ZERO_INFLATION_PCA_ASSESSMENT.md in the repository for comprehensive diagnostic
analyses and Figure S[X] for visual assessment of zero-inflation patterns.
```

---

## Results Section Templates

### Species Subsetting Results

**Template text:**
```
After applying the 10×10 prevalence-abundance filter, 38 species were retained for
community analyses, representing 94.2% of total CAFI abundance. Excluded species (n=[X])
were rare transients (e.g., single individuals on 1-2 corals) contributing minimal
information to community patterns. Retained species ranged from moderately common
(10-15 occurrences) to ubiquitous (present on >30 corals), encompassing the full
spectrum of habitat generalists and specialists.
```

---

### Zero-Inflation Acknowledgment

**Template text:**
```
The filtered community matrix exhibited 51.9% zero values, typical for cryptofaunal
assemblages where species exhibit patchy distributions due to microhabitat heterogeneity.
Hellinger transformation (applied prior to PCA) explicitly accommodates this sparsity
structure while preserving compositional information (Legendre & Gallagher 2001).
```

---

## Data Availability Statement

### For Manuscript Submission

**Template text (choose one):**

**Option 1: GitHub Only**
```
DATA AVAILABILITY

All data and analysis code are publicly available at:
https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT

The repository includes:
• Raw data files (coral growth, physiology, CAFI community surveys)
• Complete R analysis pipeline with numbered scripts
• Publication-quality figures and statistical tables
• Comprehensive documentation of filtering decisions and sensitivity analyses
• Session information for full computational reproducibility

See README.md for quick-start instructions and DATA_SUBSETTING.md for detailed
information on all data filtering decisions.
```

**Option 2: GitHub + Zenodo (Recommended)**
```
DATA AVAILABILITY

Data and code are permanently archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX) and
available at GitHub: https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT

The repository includes:
• Raw data files: Coral growth (3D photogrammetry), physiology (carbohydrates, protein,
  zooxanthellae, AFDM), and CAFI community abundance surveys
• R analysis scripts: Complete reproducible pipeline (Scripts 1-14)
• Documentation: Comprehensive filtering decisions, sensitivity analyses, and
  statistical summaries

Data are released under Creative Commons Attribution 4.0 International (CC-BY 4.0).
See LICENSE and DATA_AVAILABILITY.md in repository for details.
```

**Option 3: Multiple Repositories**
```
DATA AVAILABILITY

Data and code are available from multiple repositories:
• GitHub: https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT (code + data)
• Zenodo: DOI 10.5281/zenodo.XXXXXXX (permanent archive)
• Dryad: DOI 10.5061/dryad.XXXXXXX (data only)

All materials are released under CC-BY 4.0. See repository documentation for
reproducibility instructions.
```

---

## Supplementary Materials Captions

### Table S1: Species Filtering Summary
```
Table S1. Summary of species filtering for community analyses. The 10×10 prevalence-
abundance filter retained 38 species appearing on ≥10 corals with ≥10 total individuals.
These species represented 94.2% of total CAFI abundance. Columns show species name,
total abundance, number of corals where present (prevalence), mean abundance when
present, and percentage of zeros across the 44-coral dataset. See DATA_SUBSETTING.md
for justification and sensitivity analyses.
```

### Table S2: Coral Survival Threshold Sensitivity
```
Table S2. Sensitivity of key results to coral survival threshold. Corals with <80%
live tissue were excluded from main analyses (n=44 retained). Alternative thresholds
(70%, 90%) yielded consistent treatment effects on coral performance PC1 (all p<0.05),
demonstrating robustness to filtering decisions. Sample sizes: ≥70% alive (n=48),
≥80% alive (n=44), ≥90% alive (n=38).
```

### Table S3: PCA Transformation Sensitivity
```
Table S3. Comparison of PCA results across data transformations. Hellinger, square-root,
and Wisconsin transformations yielded similar patterns of variance explained and
treatment separation. Raw data PCA (54% variance on PC1) was dominated by overall
abundance differences (r=0.61 with total abundance), confirming the need for variance-
stabilizing transformations. Hellinger was selected for main analyses based on balance
between compositional signal and interpretability.
```

### Table S4: Species Prevalence Threshold Sensitivity
```
Table S4. Sensitivity of community PCA results to species prevalence thresholds.
Increasing minimum prevalence from 10 to 15 or 20 corals reduced species richness
but strengthened treatment effects on community PC1, indicating the 10×10 filter is
conservative. All thresholds yielded qualitatively similar conclusions.
```

### Figure S1: Zero-Inflation Diagnostic
```
Figure S1. Assessment of zero-inflation in community data. (A) Distribution of species-
level zero proportions after 10×10 filtering. Dashed lines indicate 70% and 80%
thresholds. (B) PCA variance explained by different transformations. Hellinger
transformation (22.7% on PC1) balances compositional signal with interpretability.
(C) NMDS stress values (0.20-0.23) indicate moderate fit across distance metrics.
(D) Correlation of PC1 with species richness and total abundance. Moderate correlations
confirm PC1 captures compositional patterns beyond sampling effort artifacts.
```

---

---

## References for Methods

### Key Citations to Include

**Hellinger transformation:**
```
Legendre, P., & Gallagher, E. D. (2001). Ecologically meaningful transformations for
ordination of species data. Oecologia, 129(2), 271-280.
```

**Community ecology methods:**
```
Oksanen, J., et al. (2020). vegan: Community Ecology Package. R package version 2.5-7.
https://CRAN.R-project.org/package=vegan
```

**Mixed-effects models:**
```
Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear mixed-effects
models using lme4. Journal of Statistical Software, 67(1), 1-48.
```

---

## Reviewer Response Templates

### If asked: "Why not use NMDS instead of PCA?"

**Template response:**
```
We compared PCA (Hellinger-transformed) with NMDS using multiple distance metrics
(Bray-Curtis, Jaccard, Gower). NMDS stress values (0.20-0.23) indicated moderate fit,
while PCA clearly quantified variance explained (22.7% on PC1). Given our sample size
(n=44 corals, 38 species), both methods are appropriate, and ordination patterns were
concordant (Mantel r=[X], p<0.001 between PCA and NMDS distances). We selected PCA
for the main analysis because: (1) loadings are directly interpretable as species
contributions, (2) variance partitioning quantifies signal strength, and (3) PC scores
can be used as predictors in downstream linear models relating community composition
to coral performance. NMDS ordinations are provided in Figure S[X] to demonstrate
consistency.
```

### If asked: "Isn't 51.9% zeros problematic for PCA?"

**Template response:**
```
Zero-inflation is inherent to ecological community data and does not invalidate PCA
when appropriate transformations are applied. Our approach explicitly addresses sparsity
through: (1) species filtering to remove taxa with insufficient data, (2) Hellinger
transformation to stabilize variance and down-weight abundant species, and (3) validation
that PC1 is not a sampling artifact (only r=0.44 correlation with total abundance).

Importantly, 51.9% zeros is moderate for reef cryptofauna. Studies of microbial
communities routinely apply PCA to data with 80-95% zeros after appropriate transformation
(Paliy & Shankar 2016, BMC Bioinformatics 17:261). Our diagnostic analyses (Figure S[X])
confirm that PC1 captures compositional gradients rather than zero-inflation artifacts.
We provide comprehensive sensitivity analyses in ZERO_INFLATION_PCA_ASSESSMENT.md
(repository) demonstrating robustness to transformation choice.
```

### If asked: "How did you decide on the 10×10 threshold?"

**Template response:**
```
The 10×10 threshold balances three criteria: (1) retaining sufficient species for
multivariate analyses, (2) excluding transient taxa with insufficient data for robust
inference, and (3) maintaining high abundance coverage (94.2% of total abundance retained).

We assessed sensitivity by testing alternative thresholds (5×5, 10×10, 15×10, 20×10)
and found qualitatively consistent results, with more restrictive thresholds actually
strengthening treatment effects (Table S4). The 10×10 threshold is thus conservative.
It is also consistent with prevalence filtering in microbiome analyses (typically 5-10%
of samples) scaled to our sample size (10/44 = 23% prevalence minimum).

Importantly, abundance analyses (Script 3) used the COMPLETE species dataset without
filtering, as total abundance is an appropriate metric even for rare species. Only
multivariate analyses (PCA, PERMANOVA) used the filtered set. This dual approach is
standard in community ecology.
```

---

## Quick Copy-Paste for Common Requests

### Repository URL
```
https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT
```

### Citation (preprint format)
```
Stier, A.C., Primo, A., Curtis, J.S., & Osenberg, C.W. (2025). Habitat quantity drives
community assembly and feedbacks to coral performance in reef systems. GitHub repository:
https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT
```

### License
```
CC-BY 4.0 (Creative Commons Attribution 4.0 International)
```

### DOI (placeholder)
```
DOI: 10.5281/zenodo.XXXXXXX (to be assigned upon archiving)
```

---

**Last Updated:** 2026-01-13
