# Comprehensive Sensitivity Analysis Summary

## Coverage of All Main Findings

| Finding | Test Type | % Robust | Conclusion |
|---------|-----------|----------|------------|
| Treatment → Community composition (PERMANOVA, Bray-Curtis) | Threshold grid (4×4×2) | 100% | ROBUST |
| Treatment → Community composition (PERMANOVA, Jaccard) | Threshold grid (4×4×2) | 100% | ROBUST |
| Treatment → CAFI PC1 (sqrt transform) | Threshold grid (4×4) | 100% | ROBUST |
| Treatment → CAFI PC1 (Hellinger transform) | Threshold grid (4×4) | 100% | ROBUST |
| CAFI PC1 → Coral condition PC1 (sqrt transform) | Threshold grid (4×4) | 100% | ROBUST |
| CAFI PC1 → Coral condition PC1 (Hellinger transform) | Threshold grid (4×4) | 100% | ROBUST |
| Treatment → Species richness | Threshold grid (4×4) | 0% | WEAK/VARIABLE |
| Treatment → Rarefied richness | Threshold grid (4×4) | 6% | WEAK/VARIABLE |
| Treatment → Shannon diversity | Threshold grid (4×4) | 12% | WEAK/VARIABLE |
| Treatment → Simpson diversity | Threshold grid (4×4) | 0% | WEAK/VARIABLE |
| Treatment → Evenness (Pielou's J) | Threshold grid (4×4) | 0% | WEAK/VARIABLE |
| Size × Treatment → Growth (interaction) [p=0.039]† | Bootstrap (n=100) | NOT TESTED† | NOT TESTED† |
| Treatment → Size-corrected growth [p=0.535, NS] | Bootstrap (n=100) | 0% | NOT SIG (as expected) |
| Species-coral relationships | Threshold grid (4×4) | 100% max | SOME ROBUST |
| Physiology trait correlations | Bootstrap (n=100) | 95% CI stable | ROBUST |


## Detailed Results

### Alpha Diversity
- Richness: 0/16 thresholds significant
- Rarefied richness: 1/16 thresholds significant
- Shannon: 2/16 thresholds significant
- Simpson: 0/16 thresholds significant
- Evenness (Pielou's J): 0/16 thresholds significant

### Beta Diversity (PERMANOVA)
- Bray-Curtis: 16/16 thresholds significant
- Jaccard: 16/16 thresholds significant

### CAFI PC1 (Treatment Effect)
- sqrt transform: 16/16 thresholds significant
- Hellinger transform: 16/16 thresholds significant

### CAFI-Coral Relationship
- sqrt transform: 16/16 thresholds significant
- Hellinger transform: 16/16 thresholds significant

### Growth Metrics
- Size-corrected growth (main effect, p=0.535): 0% of bootstrap iterations significant (as expected for NS finding)
- Size × Treatment interaction (p=0.039): NOT TESTED
  - †Requires raw volume data (vol_2019, vol_2021) from mesh processing
  - This significant but 'exploratory' finding should be interpreted cautiously

### Species-Coral Relationships
Robust species (sig in >=50% of thresholds):
- Dascyllus aruanus (100%)
- Fennera chacei (100%)

