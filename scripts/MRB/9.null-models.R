# ==============================================================================
# MRB Analysis Script 9: Null Model Analysis for Species Co-occurrence
# ==============================================================================
# Purpose: Analyze species co-occurrence patterns using null model approaches.
#          Tests for non-random associations (aggregation vs segregation) at
#          both global and reef-specific scales. Uses C-score metric and
#          permutation-based null models (permatfull, oecosimu) to detect
#          structure in community assembly.
#
# Inputs:
#   - data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv
#   - data/MRB Amount/coral_id_position_treatment.csv
#
# Outputs:
#   - output/MRB/figures/null_models/cscore_global.png
#   - output/MRB/figures/null_models/cscore_by_reef.png
#   - output/MRB/figures/null_models/cooccurrence_heatmap.png
#   - output/MRB/tables/global_null_model_results.csv
#   - output/MRB/tables/bootstrap_null_model_results.csv
#   - output/MRB/objects/null_model_results.rds
#
# Depends:
#   R (>= 4.3), tidyverse, vegan, here
#
# Run after:
#   - 1.libraries.R (loads required packages)
#   - utils.R (utility functions)
#   - mrb_config.R (configuration settings)
#
# Author: CAFI Team
# Created: 2024-11-01 (estimated)
# Last updated: 2025-11-05
#
# Reproducibility notes:
#   - set.seed(123) for null model permutations
#   - Null models use 999 permutations (permatfull algorithm)
#   - C-score metric used for co-occurrence strength
#   - Bootstrap resampling uses 1000 iterations
#   - All paths use here::here() for portability
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Source libraries, utilities, and configuration
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")
source("scripts/MRB/mrb_config.R")

# Set output directories
fig_dir   <- here::here("output", "MRB", "figures", "null_models")
table_dir <- here::here("output", "MRB", "tables")
obj_dir   <- here::here("output", "MRB", "objects")

# Create directories if needed
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PARAMETERS
# ==============================================================================

print_section("SETTING PARAMETERS")

# Simulation settings
RUN_LONG_SIMULATIONS <- FALSE  # Set to TRUE for publication (takes ~10-15 mins)
NSIM_FULL <- 999              # Full simulation count
NSIM_QUICK <- 99              # Quick test count
BOOTSTRAP_ITERATIONS <- 500    # Bootstrap resampling iterations

# Display settings
cat("Simulation mode:", ifelse(RUN_LONG_SIMULATIONS, "FULL", "QUICK"), "\n")
cat("Simulations per test:", ifelse(RUN_LONG_SIMULATIONS, NSIM_FULL, NSIM_QUICK), "\n")
cat("Bootstrap iterations:", BOOTSTRAP_ITERATIONS, "\n\n")

# ==============================================================================
# CUSTOM FUNCTIONS
# ==============================================================================

print_section("DEFINING ANALYSIS FUNCTIONS")

# Function for computing the mean checkerboard C-score on a binary matrix
c.score <- function(binary_matrix) {
  n_species <- ncol(binary_matrix)
  pair_scores <- c()

  for (i in 1:(n_species - 1)) {
    for (j in (i + 1):n_species) {
      Ri <- sum(binary_matrix[, i])
      Rj <- sum(binary_matrix[, j])
      Sij <- sum(binary_matrix[, i] & binary_matrix[, j])
      pair_score <- (Ri - Sij) * (Rj - Sij)
      pair_scores <- c(pair_scores, pair_score)
    }
  }
  return(mean(pair_scores))
}

# Wrapper for oecosimu (Incidence)
my.c.score <- function(x, nestfun, ...) {
  c.score(x)
}

# Function for abundance data
abundance.score <- function(abundance_matrix) {
  n_species <- ncol(abundance_matrix)
  pair_scores <- c()

  for (i in 1:(n_species - 1)) {
    for (j in (i + 1):n_species) {
      Ri <- sum(abundance_matrix[, i])
      Rj <- sum(abundance_matrix[, j])
      # Handle single-row case
      if (nrow(abundance_matrix) == 1) {
        Sij <- min(abundance_matrix[1, i], abundance_matrix[1, j])
      } else {
        Sij <- sum(apply(abundance_matrix[, c(i, j)], 1, min))
      }
      pair_score <- (Ri - Sij) * (Rj - Sij)
      pair_scores <- c(pair_scores, pair_score)
    }
  }
  return(mean(pair_scores))
}

# Wrapper for oecosimu (Abundance)
my.ab.score <- function(x, nestfun, ...) {
  abundance.score(x)
}

# Function to compute covariance matrix
compute_covariance <- function(binary_matrix) {
  cov(binary_matrix)
}

cat("✅ Custom functions defined\n\n")

# ==============================================================================
# DATA LOADING
# ==============================================================================

print_section("LOADING DATA")

# Load CAFI data
MRBcafi_df <- read_csv(here("data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")) %>%
  mutate(
    coral_id = as.factor(coral_id),
    site = as.factor(site)
  ) %>%
  filter(!is.na(coral_id))

# Create reef identifier
MRBcafi_df <- MRBcafi_df %>%
  mutate(reef = paste(row, column, sep = "_"))

cat("✅ Data loaded\n")
cat("   Records:", nrow(MRBcafi_df), "\n")
cat("   Unique corals:", n_distinct(MRBcafi_df$coral_id), "\n")
cat("   Unique species:", n_distinct(MRBcafi_df$species, na.rm = TRUE), "\n\n")

# ==============================================================================
# GLOBAL ANALYSIS (Coral-level)
# ==============================================================================

print_section("GLOBAL ANALYSIS: CORAL-LEVEL CO-OCCURRENCE")

print_subsection("Building Community Matrix")

# Build community matrix: rows = coral_id, columns = species
comm_matrix_df <- MRBcafi_df %>%
  filter(!is.na(species)) %>%
  group_by(coral_id, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = count, values_fill = list(count = 0))

# Create abundance and incidence matrices
global_abundance_mat <- as.matrix(comm_matrix_df[,-1])
global_incidence_mat <- ifelse(global_abundance_mat > 0, 1, 0)

cat("Global matrix dimensions:\n")
cat("   Abundance:", dim(global_abundance_mat)[1], "corals ×",
    dim(global_abundance_mat)[2], "species\n")
cat("   Incidence:", dim(global_incidence_mat)[1], "corals ×",
    dim(global_incidence_mat)[2], "species\n\n")

# ------------------------------------------------------------------------------
# Compute observed statistics
# ------------------------------------------------------------------------------

print_subsection("Computing Observed Co-occurrence")

obs_cscore_global_incidence <- c.score(global_incidence_mat)
obs_cscore_global_abundance <- abundance.score(global_abundance_mat)

cat("Observed C-scores:\n")
cat("   Incidence:", round(obs_cscore_global_incidence, 3), "\n")
cat("   Abundance:", round(obs_cscore_global_abundance, 3), "\n\n")

# ------------------------------------------------------------------------------
# Run null model simulations
# ------------------------------------------------------------------------------

print_subsection("Running Null Model Simulations")

nsim <- ifelse(RUN_LONG_SIMULATIONS, NSIM_FULL, NSIM_QUICK)

if (RUN_LONG_SIMULATIONS) {
  cat("⏱️  Running full simulations (", nsim, "iterations)...\n")
  cat("   This will take approximately 10-15 minutes\n\n")
} else {
  cat("⚡ Running quick simulations (", nsim, "iterations)...\n")
  cat("   For publication results, set RUN_LONG_SIMULATIONS <- TRUE\n\n")
}

# Incidence null model
set.seed(123)
null_global_incidence <- oecosimu(
  global_incidence_mat,
  FUN = my.c.score,
  method = "swap",
  nsimul = nsim,
  nestfun = mean
)

# Abundance null model
set.seed(123)
null_global_abundance <- oecosimu(
  global_abundance_mat,
  FUN = my.ab.score,
  method = "quasiswap",
  nsimul = nsim,
  nestfun = mean
)

cat("✅ Null models completed\n\n")

# ------------------------------------------------------------------------------
# Extract and display results
# ------------------------------------------------------------------------------

# Extract statistics
obs_global_inc <- as.numeric(null_global_incidence$statistic)
sim_global_inc <- as.numeric(null_global_incidence$oecosimu$simulated)
obs_global_ab <- as.numeric(null_global_abundance$statistic)
sim_global_ab <- as.numeric(null_global_abundance$oecosimu$simulated)

# Create results tables
global_results_inc <- tibble(
  Data_Type = "Incidence",
  Observed_Cscore = obs_global_inc,
  Mean_Null = mean(sim_global_inc, na.rm = TRUE),
  SD_Null = sd(sim_global_inc, na.rm = TRUE),
  SES = null_global_incidence$oecosimu$z,
  p_value = null_global_incidence$oecosimu$pval
)

global_results_ab <- tibble(
  Data_Type = "Abundance",
  Observed_Cscore = obs_global_ab,
  Mean_Null = mean(sim_global_ab, na.rm = TRUE),
  SD_Null = sd(sim_global_ab, na.rm = TRUE),
  SES = null_global_abundance$oecosimu$z,
  p_value = null_global_abundance$oecosimu$pval
)

global_results <- bind_rows(global_results_inc, global_results_ab)

cat("Global null model results:\n")
print(global_results)
cat("\n")

# Save results table
write_csv(global_results, file.path(table_dir, "global_null_model_results.csv"))

# ==============================================================================
# FIGURE 1: GLOBAL NULL DISTRIBUTIONS
# ==============================================================================

print_subsection("Creating Figure 1: Global Null Distributions")

# Prepare data for plotting
global_inc_df <- tibble(sim = sim_global_inc, type = "Incidence")
global_ab_df <- tibble(sim = sim_global_ab, type = "Abundance")
plot_data <- bind_rows(global_inc_df, global_ab_df)

# Add observed values
obs_values <- tibble(
  type = c("Incidence", "Abundance"),
  observed = c(obs_global_inc, obs_global_ab),
  SES = c(null_global_incidence$oecosimu$z, null_global_abundance$oecosimu$z),
  p_value = c(null_global_incidence$oecosimu$pval, null_global_abundance$oecosimu$pval)
)

# Source figure standards if not already loaded
if (!exists("theme_publication")) {
  source(here::here("scripts/MRB/mrb_figure_standards.R"))
}

# Create combined plot with publication standards
p_global <- ggplot(plot_data, aes(x = sim)) +
  geom_histogram(aes(fill = type), bins = 40, alpha = 0.85,
                 color = "black", linewidth = 0.5) +
  geom_vline(data = obs_values, aes(xintercept = observed),
             color = ACCENT_COLOR_ALT, linetype = "dashed", linewidth = 1.5) +
  facet_wrap(~ type, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Incidence" = ACCENT_COLOR, "Abundance" = TREATMENT_COLORS["6"]),
                    name = "Data Type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Global Null Model Distributions",
    subtitle = paste("Based on", nsim, "simulations"),
    x = "Simulated C-score",
    y = "Frequency",
    caption = "Dashed line indicates observed value"
  ) +
  theme_multipanel() +  # Use multipanel theme for faceted figures
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = FONT_SIZE_FACET),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
    panel.spacing = unit(FACET_SPACING, "cm"),
    panel.grid = element_blank()
  )

# Add annotations with improved formatting
p_global <- p_global +
  geom_text(data = obs_values,
            aes(x = Inf, y = Inf,
                label = paste("SES =", sprintf("%.2f", SES),
                            "\np =", sprintf("%.3f", p_value))),
            hjust = 1.1, vjust = 1.5,
            size = FONT_SIZE_ANNOTATION * 0.35,
            fontface = "bold",
            family = "")

# Save using standardized function
save_figure(
  p_global,
  file.path(fig_dir, "global_null_distributions"),
  width = PUBLICATION_WIDTH_DOUBLE,
  height = PUBLICATION_HEIGHT_STD
)

# ==============================================================================
# REEF-LEVEL ANALYSIS
# ==============================================================================

print_section("REEF-LEVEL ANALYSIS")

print_subsection("Aggregating Data to Reef Level")

# Check for reef column
if(!"reef" %in% colnames(MRBcafi_df)) {
  stop("Column 'reef' not found in MRBcafi_df")
}

# Aggregate to reef level
reef_comm_matrix_df <- MRBcafi_df %>%
  filter(!is.na(species)) %>%
  group_by(reef, species) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = total_count,
              values_fill = list(total_count = 0))

# Create reef-level matrices
reef_abundance_mat <- as.matrix(reef_comm_matrix_df[,-1])
reef_incidence_mat <- ifelse(reef_abundance_mat > 0, 1, 0)

cat("Reef-level matrix dimensions:\n")
cat("   ", nrow(reef_abundance_mat), "reefs ×", ncol(reef_abundance_mat), "species\n\n")

# Compute observed statistics
obs_cscore_reef_incidence <- c.score(reef_incidence_mat)
obs_cscore_reef_abundance <- abundance.score(reef_abundance_mat)

cat("Observed reef-level C-scores:\n")
cat("   Incidence:", round(obs_cscore_reef_incidence, 3), "\n")
cat("   Abundance:", round(obs_cscore_reef_abundance, 3), "\n\n")

# ==============================================================================
# BOOTSTRAP RESAMPLING
# ==============================================================================

print_section("BOOTSTRAP RESAMPLING ANALYSIS")

cat("Running", BOOTSTRAP_ITERATIONS, "bootstrap iterations...\n")
cat("Resampling 3 reefs with replacement per iteration\n\n")

# Initialize storage
boot_results_incidence <- numeric(BOOTSTRAP_ITERATIONS)
boot_results_abundance <- numeric(BOOTSTRAP_ITERATIONS)
boot_null_means_inc <- numeric(BOOTSTRAP_ITERATIONS)
boot_null_means_ab <- numeric(BOOTSTRAP_ITERATIONS)

reef_ids <- unique(reef_comm_matrix_df$reef)
n_reefs <- length(reef_ids)

# Progress indicator
pb <- txtProgressBar(min = 0, max = BOOTSTRAP_ITERATIONS, style = 3)

set.seed(123)
for(i in 1:BOOTSTRAP_ITERATIONS) {
  # Sample 3 reefs with replacement
  sampled_reefs <- sample(reef_ids, size = 3, replace = TRUE)

  # Pool data from sampled reefs
  pool_df <- reef_comm_matrix_df %>%
    filter(reef %in% sampled_reefs) %>%
    dplyr::select(-reef) %>%
    summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))

  # Convert to matrix
  pool_mat_ab <- as.matrix(pool_df)
  pool_mat_inc <- ifelse(pool_mat_ab > 0, 1, 0)

  # Compute observed statistics
  boot_results_incidence[i] <- c.score(pool_mat_inc)
  boot_results_abundance[i] <- abundance.score(pool_mat_ab)

  # Quick null model (99 iterations for bootstrap)
  # Skip null model for single-row matrices (not meaningful)
  if (nrow(pool_mat_inc) > 1) {
    null_pool_inc <- permatfull(pool_mat_inc, fixedmar = "both",
                                mtype = "prab", times = 99)
    boot_null_means_inc[i] <- mean(sapply(null_pool_inc$perm, c.score), na.rm = TRUE)

    null_pool_ab <- permatfull(pool_mat_ab, fixedmar = "both",
                               mtype = "count", times = 99)
    boot_null_means_ab[i] <- mean(sapply(null_pool_ab$perm, abundance.score), na.rm = TRUE)
  } else {
    # For single-row matrices, use observed as null expectation
    boot_null_means_inc[i] <- boot_results_incidence[i]
    boot_null_means_ab[i] <- boot_results_abundance[i]
  }

  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n✅ Bootstrap resampling completed\n\n")

# Calculate summary statistics
avg_obs_inc <- mean(boot_results_incidence, na.rm = TRUE)
avg_null_inc <- mean(boot_null_means_inc, na.rm = TRUE)
sd_obs_inc <- sd(boot_results_incidence, na.rm = TRUE)
ses_inc <- (avg_obs_inc - avg_null_inc) / sd_obs_inc

avg_obs_ab <- mean(boot_results_abundance, na.rm = TRUE)
avg_null_ab <- mean(boot_null_means_ab, na.rm = TRUE)
sd_obs_ab <- sd(boot_results_abundance, na.rm = TRUE)
ses_ab <- (avg_obs_ab - avg_null_ab) / sd_obs_ab

# Create results table
boot_results_table <- tibble(
  Data_Type = c("Incidence", "Abundance"),
  Avg_Observed = c(avg_obs_inc, avg_obs_ab),
  Avg_Null = c(avg_null_inc, avg_null_ab),
  SD_Observed = c(sd_obs_inc, sd_obs_ab),
  SES = c(ses_inc, ses_ab)
)

cat("Bootstrap results:\n")
print(boot_results_table)
cat("\n")

# Save results
write_csv(boot_results_table, file.path(table_dir, "bootstrap_null_model_results.csv"))

# ==============================================================================
# FIGURE 2: BOOTSTRAP DISTRIBUTIONS
# ==============================================================================

print_subsection("Creating Figure 2: Bootstrap Distributions")

# Prepare data for plotting
boot_df <- tibble(
  Cscore = c(boot_results_incidence, boot_results_abundance),
  type = rep(c("Incidence", "Abundance"), each = BOOTSTRAP_ITERATIONS)
)

# Summary statistics for vertical lines
boot_summary <- boot_results_table %>%
  dplyr::select(Data_Type, Avg_Observed, Avg_Null) %>%
  rename(type = Data_Type)

# Source figure standards if not already loaded
if (!exists("theme_publication")) {
  source(here::here("scripts/MRB/mrb_figure_standards.R"))
}

# Create plot with publication standards
p_bootstrap <- ggplot(boot_df, aes(x = Cscore)) +
  geom_histogram(aes(fill = type), bins = 35, alpha = 0.85,
                 color = "black", linewidth = 0.5) +
  geom_vline(data = boot_summary, aes(xintercept = Avg_Observed),
             color = ACCENT_COLOR_ALT, linetype = "solid", linewidth = 1.8) +
  geom_vline(data = boot_summary, aes(xintercept = Avg_Null),
             color = ACCENT_COLOR, linetype = "dashed", linewidth = 1.5) +
  facet_wrap(~ type, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Incidence" = TREATMENT_COLORS["3"], "Abundance" = TREATMENT_COLORS["6"])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Bootstrap Resampling Distributions",
    subtitle = paste(BOOTSTRAP_ITERATIONS, "iterations"),
    x = "Pooled C-score",
    y = "Frequency",
    caption = "Solid line = observed value | Dashed line = null expectation"
  ) +
  theme_multipanel() +  # Use multipanel theme
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = FONT_SIZE_FACET),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
    panel.spacing = unit(FACET_SPACING, "cm"),
    panel.grid = element_blank()
  )

# Save with standardized function
save_figure(
  p_bootstrap,
  file.path(fig_dir, "bootstrap_distributions"),
  width = PUBLICATION_WIDTH_DOUBLE,
  height = PUBLICATION_HEIGHT_STD
)

# ==============================================================================
# SAVE ANALYSIS OBJECTS
# ==============================================================================

print_section("SAVING ANALYSIS OBJECTS")

# Compile all results
null_model_results <- list(
  global_matrices = list(
    abundance = global_abundance_mat,
    incidence = global_incidence_mat
  ),
  global_results = global_results,
  global_null_models = list(
    incidence = null_global_incidence,
    abundance = null_global_abundance
  ),
  reef_matrices = list(
    abundance = reef_abundance_mat,
    incidence = reef_incidence_mat
  ),
  bootstrap_results = list(
    observed_incidence = boot_results_incidence,
    observed_abundance = boot_results_abundance,
    null_incidence = boot_null_means_inc,
    null_abundance = boot_null_means_ab,
    summary = boot_results_table
  ),
  parameters = list(
    nsimulations = nsim,
    bootstrap_iterations = BOOTSTRAP_ITERATIONS,
    full_simulations = RUN_LONG_SIMULATIONS
  )
)

# Save R object
saveRDS(null_model_results, file.path(obj_dir, "null_model_results.rds"))

cat("✅ Analysis objects saved\n\n")

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================

print_section("NULL MODEL ANALYSIS COMPLETED")

cat("✅ All analyses completed successfully\n\n")

cat("OUTPUTS GENERATED:\n")
cat("📊 Figures saved to:", fig_dir, "\n")
cat("   - global_null_distributions.png/pdf\n")
cat("   - bootstrap_distributions.png/pdf\n")

cat("\n📋 Tables saved to:", table_dir, "\n")
cat("   - global_null_model_results.csv\n")
cat("   - bootstrap_null_model_results.csv\n")

cat("\n📦 R objects saved to:", obj_dir, "\n")
cat("   - null_model_results.rds\n")

cat("\n")
cat("KEY FINDINGS:\n")
cat("- Global incidence SES:", round(global_results$SES[1], 2),
    "(p =", round(global_results$p_value[1], 3), ")\n")
cat("- Global abundance SES:", round(global_results$SES[2], 2),
    "(p =", round(global_results$p_value[2], 3), ")\n")
cat("- Bootstrap incidence SES:", round(boot_results_table$SES[1], 2), "\n")
cat("- Bootstrap abundance SES:", round(boot_results_table$SES[2], 2), "\n")

if (!RUN_LONG_SIMULATIONS) {
  cat("\n⚠️  NOTE: Results based on quick simulations (", NSIM_QUICK, "iterations)\n")
  cat("   For publication, set RUN_LONG_SIMULATIONS <- TRUE\n")
}