# Data Availability Statement

## Overview

All data and code associated with this manuscript are publicly available.

## Data Repository

**Primary Repository:** GitHub
**URL:** https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT
**DOI:** [To be assigned upon Zenodo archival]

## Archived Version

Upon publication, a permanent archived version will be deposited at:
- **Zenodo:** [DOI to be added]
- **Dryad:** [DOI to be added, if applicable]

## Data Description

| Dataset | Description | Format | Size |
|---------|-------------|--------|------|
| CAFI community data | Fish and invertebrate observations | CSV | ~1 MB |
| Coral physiology | Tissue composition measurements | CSV | ~50 KB |
| 3D photogrammetry | Coral size and growth metrics | CSV | ~100 KB |
| Treatment assignments | Experimental design | CSV | ~5 KB |

## Access

- **License:** CC-BY 4.0 (Creative Commons Attribution)
- **Access:** Unrestricted public access
- **Format:** Plain-text CSV files (no proprietary formats)

## How to Access

1. **Clone the repository:**
   ```bash
   git clone https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT.git
   ```

2. **Download as ZIP:**
   - Visit the repository URL
   - Click "Code" → "Download ZIP"

3. **Archived version:**
   - Visit the Zenodo DOI link (available upon publication)

## Code Availability

All analysis code is included in the repository:

| Directory | Contents |
|-----------|----------|
| `scripts/MRB/` | R analysis scripts (numbered 1-14) |
| `scripts/MRB/utils.R` | Utility functions |
| `scripts/MRB/mrb_figure_standards.R` | Figure formatting |

## Reproducibility

To reproduce all analyses:

1. Install R (≥4.3.0) and required packages (see `scripts/MRB/1.libraries.R`)
2. Run scripts in numerical order (1 → 14)
3. See `REPRODUCIBILITY_GUIDE.md` for detailed instructions

## Citation

If you use this data or code, please cite:

```
Stier, A.C., Primo, A., Curtis, J.S., & Osenberg, C.W. (2025). Habitat quantity
drives community assembly and feedbacks to coral performance. Ecology Letters.
DOI: [to be added]
```

## Contact

For questions about data access:
- **Email:** adrian.stier@lifesci.ucsb.edu
- **Issues:** https://github.com/adrianstier/Stier-2025-CAFI136-MRB-AMOUNT/issues

## Ethical Statement

- Field work was conducted with appropriate permits from the French Polynesian government
- No endangered species were collected
- All coral manipulations followed established ethical guidelines for reef research
