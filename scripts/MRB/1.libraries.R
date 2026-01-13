# ==============================================================================
# MRB Analysis Script 1: Central Library Loading
# ==============================================================================
# Purpose: Loads all required R packages for the CAFI MRB analysis pipeline.
#          Source this file at the beginning of each analysis script to ensure
#          all necessary packages are available. Includes data manipulation,
#          visualization, statistical modeling, and ecological analysis packages.
#
# Inputs:  None (pure library loading)
# Outputs: None (loads packages into R session memory only)
#
# Depends:
#   R (>= 4.3), all packages listed below must be installed
#
# Run after:
#   None (this is always the first script to run)
#
# Author: CAFI Team
# Created: 2024-01-01 (estimated)
# Last updated: 2025-11-05
#
# Reproducibility notes:
#   - Explicit library() calls ensure package loading order
#   - All packages are from CRAN unless noted
#   - See sessionInfo() for exact package versions used
#   - Install missing packages with: install.packages("package_name")
# ==============================================================================

# Core tidyverse and data manipulation
library(tidyverse)  # Includes dplyr, ggplot2, tidyr, readr, purrr, stringr
library(dplyr)      # Explicitly load for clarity (included in tidyverse)
library(purrr)      # Functional programming (included in tidyverse)
library(stringr)    # String manipulation (included in tidyverse)
library(forcats)    # Factor handling

# File paths and I/O
library(here)       # Reproducible file paths
library(readxl)     # Excel file reading

# Visualization
library(ggplot2)    # Plotting (included in tidyverse)
library(patchwork)  # Combining plots
library(GGally)     # ggpairs and extensions
library(ggtext)     # Enhanced text rendering
library(ggalluvial) # Alluvial diagrams
library(ggrepel)    # Text label repelling
library(scales)     # Scale functions (comma_format, etc.)
library(gridExtra)  # Grid arrangements
library(htmlwidgets)# Interactive widgets
library(waffle)     # Waffle charts

# Tables and reporting
library(gt)         # Publication-quality tables
library(pagedown)   # HTML to PDF conversion

# Statistical modeling
library(lmerTest)    # Mixed effects models with p-values
library(broom)       # Tidy model outputs
library(broom.mixed) # Tidy mixed effects models
library(merTools)    # Mixed model diagnostics
library(emmeans)     # Estimated marginal means
library(MuMIn)       # Model selection

# Multivariate analysis
library(vegan)      # Community ecology analyses
library(cluster)    # Clustering algorithms
library(labdsv)     # Lab for data science in vegetation science
library(Hmisc)      # Harrell miscellaneous functions

# Additional specialized packages
library(indicspecies) # Indicator species analysis
library(ape)          # Phylogenetic analyses and PCoA
library(lme4)         # Linear mixed-effects models
library(reshape2)     # Data reshaping (for matrices)
library(viridis)      # Color scales for visualization
library(tibble)       # Enhanced data frames
library(readr)        # Fast CSV reading (included in tidyverse)
library(tidyr)        # Data reshaping (included in tidyverse)

# Additional packages for specific scripts
library(cli)          # Command line interface messaging
library(fitdistrplus) # Fitting distributions
library(performance)  # Model performance checks
library(car)          # Companion to Applied Regression
library(fs)           # File system operations
library(RColorBrewer) # Color palettes (used in script 12)

cat("✅ All libraries loaded successfully\n")
