#!/usr/bin/env Rscript
# ==============================================================================
# Zero-Inflation Analysis for PCA Validity Assessment
# ==============================================================================
# Purpose: Investigate sparsity patterns in CAFI community data to assess
#          whether zero-inflation poses problems for PCA
#
# Key Questions:
# 1. What proportion of the data matrix is zeros?
# 2. How many zeros per species (after 10×10 filtering)?
# 3. Is zero-inflation causing PCA to extract spurious axes?
# 4. Are alternative approaches (NMDS, correspondence analysis) more appropriate?
#
# Generated: 2026-01-13
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(vegan)
  library(cli)
})

# --- Configuration ---
DATA_DIR <- here("data", "MRB Amount")
OUT_DIR  <- here("output", "MRB")
TAB_DIR  <- file.path(OUT_DIR, "tables")
FIG_DIR  <- file.path(OUT_DIR, "figures")

PREV_MIN  <- 10
ABUND_MIN <- 10

cli_h1("Zero-Inflation Analysis for PCA")

# =============================================================================
# 1. Load and filter data (matching Script 8 pipeline)
# =============================================================================

cli_h2("1. Loading CAFI community data")

# Helper function to strip "FE-" prefix
strip_fe <- function(x) gsub("^FE-", "", as.character(x))

# Load CAFI data
cafi_file <- file.path(DATA_DIR, "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")
cafi_raw <- read_csv(cafi_file, show_col_types = FALSE)

# Clean and pivot
cafi_neat <- cafi_raw %>%
  mutate(coral_id = strip_fe(coral_id)) %>%
  select(coral_id, species, count) %>%
  filter(!is.na(species), !is.na(count))

# Create wide matrix
comm_wide <- cafi_neat %>%
  group_by(coral_id, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0)

# Apply 10×10 filter
species_meta <- cafi_neat %>%
  group_by(species) %>%
  summarise(
    total_count = sum(count, na.rm = TRUE),
    n_corals    = n_distinct(coral_id[count > 0]),
    .groups     = "drop"
  ) %>%
  filter(total_count >= ABUND_MIN, n_corals >= PREV_MIN)

cli_alert_success("Species retained after 10×10 filter: {nrow(species_meta)}")

# Filter matrix to retained species
comm_filtered <- comm_wide %>%
  select(coral_id, all_of(species_meta$species))

# Load filtered coral list (≥80% alive)
growth_file <- here("data/processed/coral_growth.csv")
if (!file.exists(growth_file)) {
  cli_alert_warning("Growth file not found; using all corals from community data")
  keep_ids <- comm_filtered$coral_id
} else {
  growth <- read_csv(growth_file, show_col_types = FALSE) %>%
    mutate(coral_id = strip_fe(as.character(coral_id)))
  keep_ids <- growth$coral_id
  comm_filtered <- comm_filtered %>%
    filter(coral_id %in% keep_ids)
}

cli_alert_success("Corals retained: {nrow(comm_filtered)}")

# Convert to matrix
comm_mat <- comm_filtered %>%
  column_to_rownames("coral_id") %>%
  as.matrix()

n_corals  <- nrow(comm_mat)
n_species <- ncol(comm_mat)
n_cells   <- n_corals * n_species

cli_alert_info("Final matrix dimensions: {n_corals} corals × {n_species} species = {n_cells} cells")

# =============================================================================
# 2. Calculate sparsity metrics
# =============================================================================

cli_h2("2. Calculating zero-inflation metrics")

# Overall sparsity
n_zeros <- sum(comm_mat == 0)
prop_zeros <- n_zeros / n_cells

cli_alert_info("Total zeros: {n_zeros} / {n_cells} ({round(prop_zeros * 100, 1)}%)")

# Per-species sparsity
species_sparsity <- tibble(
  species = colnames(comm_mat),
  n_zeros = colSums(comm_mat == 0),
  prop_zeros = n_zeros / n_corals,
  n_presences = colSums(comm_mat > 0),
  prop_presences = n_presences / n_corals,
  total_count = colSums(comm_mat),
  mean_when_present = colSums(comm_mat) / pmax(n_presences, 1)
) %>%
  arrange(desc(prop_zeros))

# Per-coral sparsity
coral_sparsity <- tibble(
  coral_id = rownames(comm_mat),
  n_zeros = rowSums(comm_mat == 0),
  prop_zeros = n_zeros / n_species,
  n_species_present = rowSums(comm_mat > 0),
  prop_species_present = n_species_present / n_species,
  total_abundance = rowSums(comm_mat)
)

# Summary statistics
cli_alert_info("Species sparsity range: {round(min(species_sparsity$prop_zeros)*100, 1)}% to {round(max(species_sparsity$prop_zeros)*100, 1)}%")
cli_alert_info("Mean species sparsity: {round(mean(species_sparsity$prop_zeros)*100, 1)}%")
cli_alert_info("Median species sparsity: {round(median(species_sparsity$prop_zeros)*100, 1)}%")

# How many species have >70%, >80%, >90% zeros?
pct_gt70 <- mean(species_sparsity$prop_zeros > 0.70) * 100
pct_gt80 <- mean(species_sparsity$prop_zeros > 0.80) * 100
pct_gt90 <- mean(species_sparsity$prop_zeros > 0.90) * 100

cli_alert_warning("{round(pct_gt70, 1)}% of species have >70% zeros")
cli_alert_warning("{round(pct_gt80, 1)}% of species have >80% zeros")
cli_alert_warning("{round(pct_gt90, 1)}% of species have >90% zeros")

# =============================================================================
# 3. Compare transformation approaches
# =============================================================================

cli_h2("3. Comparing data transformations")

# Raw data
comm_raw <- comm_mat

# Square-root transform
comm_sqrt <- sqrt(comm_mat)

# Hellinger transform (used in Script 8)
comm_hellinger <- decostand(comm_mat, method = "hellinger")

# Presence-absence
comm_pa <- decostand(comm_mat, method = "pa")

# Log(x+1) transform
comm_log <- log1p(comm_mat)

# Wisconsin double standardization (species max = 1, then site total = 1)
comm_wisconsin <- wisconsin(comm_mat)

# =============================================================================
# 4. Run PCA on different transformations
# =============================================================================

cli_h2("4. Running PCA on different transformations")

transformations <- list(
  raw        = comm_raw,
  sqrt       = comm_sqrt,
  hellinger  = comm_hellinger,
  log1p      = comm_log,
  wisconsin  = comm_wisconsin
)

pca_results <- list()
pca_summaries <- tibble()

for (trans_name in names(transformations)) {
  cli_alert_info("Running PCA on {trans_name} transformation...")

  trans_mat <- transformations[[trans_name]]

  # Remove zero-variance columns
  trans_mat <- trans_mat[, apply(trans_mat, 2, var) > 0, drop = FALSE]

  # Run PCA
  pca <- prcomp(trans_mat, center = TRUE, scale. = FALSE)

  # Extract variance explained
  var_explained <- summary(pca)$importance[2, ]  # Proportion of Variance row

  # Store results
  pca_results[[trans_name]] <- pca

  pca_summaries <- bind_rows(
    pca_summaries,
    tibble(
      transformation = trans_name,
      n_species = ncol(trans_mat),
      pc1_var = var_explained[1] * 100,
      pc2_var = var_explained[2] * 100,
      pc3_var = var_explained[3] * 100,
      cumul_3pc = sum(var_explained[1:3]) * 100,
      total_variance = sum(apply(trans_mat, 2, var))
    )
  )

  cli_alert_success("  PC1: {round(var_explained[1]*100, 1)}%, PC2: {round(var_explained[2]*100, 1)}%, PC3: {round(var_explained[3]*100, 1)}%")
}

# =============================================================================
# 5. Assess PCA quality metrics
# =============================================================================

cli_h2("5. Assessing PCA quality")

# Kaiser-Meyer-Olkin (KMO) test for sampling adequacy
# Values > 0.5 suggest PCA is appropriate
kmo_results <- tibble()

for (trans_name in names(transformations)) {
  trans_mat <- transformations[[trans_name]]
  trans_mat <- trans_mat[, apply(trans_mat, 2, var) > 0, drop = FALSE]

  # Correlation matrix
  cor_mat <- cor(trans_mat)

  # Bartlett's test: tests if correlation matrix is identity matrix
  # Low p-value = correlations exist, PCA appropriate
  n <- nrow(trans_mat)
  p <- ncol(trans_mat)

  # Determinant of correlation matrix
  det_cor <- det(cor_mat)

  # KMO statistic (simplified version)
  # Full KMO calculation is complex; here we use determinant as proxy
  kmo_proxy <- ifelse(det_cor < 1e-10, "Poor",
                      ifelse(det_cor < 0.01, "Mediocre",
                             ifelse(det_cor < 0.1, "Middling", "Good")))

  kmo_results <- bind_rows(
    kmo_results,
    tibble(
      transformation = trans_name,
      det_correlation = det_cor,
      kmo_category = kmo_proxy
    )
  )
}

# =============================================================================
# 6. Compare with NMDS (non-metric multidimensional scaling)
# =============================================================================

cli_h2("6. Comparing PCA with NMDS")

# NMDS is designed for sparse, zero-inflated ecological data
# Run on multiple distance metrics

nmds_results <- tibble()

distance_methods <- c("bray", "jaccard", "gower")

for (dist_method in distance_methods) {
  cli_alert_info("Running NMDS with {dist_method} distance...")

  set.seed(123)
  nmds <- metaMDS(comm_mat, distance = dist_method, k = 2, try = 20,
                  autotransform = FALSE, trace = 0)

  nmds_results <- bind_rows(
    nmds_results,
    tibble(
      distance = dist_method,
      stress = nmds$stress,
      converged = nmds$converged,
      interpretation = case_when(
        nmds$stress < 0.05 ~ "Excellent",
        nmds$stress < 0.10 ~ "Good",
        nmds$stress < 0.20 ~ "Fair",
        nmds$stress < 0.30 ~ "Poor",
        TRUE ~ "Unacceptable"
      )
    )
  )

  cli_alert_success("  Stress = {round(nmds$stress, 3)} ({nmds_results$interpretation[nrow(nmds_results)]})")
}

# =============================================================================
# 7. Test for spurious correlations due to zeros
# =============================================================================

cli_h2("7. Testing for zero-inflation artifacts")

# Calculate correlation between species presence and PC1 scores
# If high, PC1 may be driven by presence/absence rather than abundance patterns

pca_raw <- pca_results$raw
pca_sqrt <- pca_results$sqrt
pca_hell <- pca_results$hellinger

# Species richness per coral
coral_richness <- rowSums(comm_mat > 0)

# Total abundance per coral
coral_abundance <- rowSums(comm_mat)

# Correlations with PC1
pc1_correlations <- tibble(
  transformation = c("raw", "sqrt", "hellinger"),
  cor_richness = c(
    cor(coral_richness, pca_raw$x[, 1]),
    cor(coral_richness, pca_sqrt$x[, 1]),
    cor(coral_richness, pca_hell$x[, 1])
  ),
  cor_abundance = c(
    cor(coral_abundance, pca_raw$x[, 1]),
    cor(coral_abundance, pca_sqrt$x[, 1]),
    cor(coral_abundance, pca_hell$x[, 1])
  )
)

cli_alert_info("Correlation between PC1 and species richness:")
print(pc1_correlations)

# Strong correlation (|r| > 0.7) suggests PC1 is a "sampling effort" axis
for (i in 1:nrow(pc1_correlations)) {
  r_rich <- pc1_correlations$cor_richness[i]
  r_abund <- pc1_correlations$cor_abundance[i]
  trans <- pc1_correlations$transformation[i]

  if (abs(r_rich) > 0.7) {
    cli_alert_warning("  {trans}: PC1 strongly correlated with richness (r={round(r_rich, 2)})")
  }
  if (abs(r_abund) > 0.7) {
    cli_alert_warning("  {trans}: PC1 strongly correlated with abundance (r={round(r_abund, 2)})")
  }
}

# =============================================================================
# 8. Export results
# =============================================================================

cli_h2("8. Exporting results")

# Species sparsity table
write_csv(species_sparsity, file.path(TAB_DIR, "zero_inflation_species_sparsity.csv"))
cli_alert_success("Saved: zero_inflation_species_sparsity.csv")

# Coral sparsity table
write_csv(coral_sparsity, file.path(TAB_DIR, "zero_inflation_coral_sparsity.csv"))
cli_alert_success("Saved: zero_inflation_coral_sparsity.csv")

# PCA comparison table
write_csv(pca_summaries, file.path(TAB_DIR, "zero_inflation_pca_comparison.csv"))
cli_alert_success("Saved: zero_inflation_pca_comparison.csv")

# KMO results
write_csv(kmo_results, file.path(TAB_DIR, "zero_inflation_kmo_results.csv"))
cli_alert_success("Saved: zero_inflation_kmo_results.csv")

# NMDS results
write_csv(nmds_results, file.path(TAB_DIR, "zero_inflation_nmds_comparison.csv"))
cli_alert_success("Saved: zero_inflation_nmds_comparison.csv")

# PC1 correlations
write_csv(pc1_correlations, file.path(TAB_DIR, "zero_inflation_pc1_correlations.csv"))
cli_alert_success("Saved: zero_inflation_pc1_correlations.csv")

# =============================================================================
# 9. Create diagnostic figure
# =============================================================================

cli_h2("9. Creating diagnostic figure")

library(patchwork)

# Panel A: Species sparsity distribution
p1 <- ggplot(species_sparsity, aes(x = prop_zeros * 100)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 70, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 80, linetype = "dashed", color = "darkred") +
  labs(
    title = "A. Species zero-inflation distribution",
    x = "Percent zeros per species",
    y = "Number of species"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 10))

# Panel B: PCA variance explained by transformation
p2 <- pca_summaries %>%
  select(transformation, PC1 = pc1_var, PC2 = pc2_var, PC3 = pc3_var) %>%
  pivot_longer(-transformation, names_to = "PC", values_to = "var_explained") %>%
  ggplot(aes(x = transformation, y = var_explained, fill = PC)) +
  geom_col(position = "dodge") +
  labs(
    title = "B. PCA variance by transformation",
    x = "Transformation",
    y = "Variance explained (%)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 10))

# Panel C: NMDS stress by distance
p3 <- ggplot(nmds_results, aes(x = distance, y = stress, fill = interpretation)) +
  geom_col() +
  geom_hline(yintercept = 0.20, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Excellent" = "darkgreen", "Good" = "green",
                                "Fair" = "yellow", "Poor" = "orange")) +
  labs(
    title = "C. NMDS stress by distance metric",
    x = "Distance metric",
    y = "Stress"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 10))

# Panel D: PC1 correlations
p4 <- pc1_correlations %>%
  pivot_longer(-transformation, names_to = "metric", values_to = "correlation") %>%
  mutate(metric = recode(metric,
                         cor_richness = "Species richness",
                         cor_abundance = "Total abundance")) %>%
  ggplot(aes(x = transformation, y = correlation, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = c(-0.7, 0.7), linetype = "dashed", color = "red") +
  labs(
    title = "D. PC1 correlation with richness/abundance",
    x = "Transformation",
    y = "Pearson correlation"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 10))

# Combine panels
combined <- (p1 | p2) / (p3 | p4)

ggsave(
  file.path(FIG_DIR, "zero_inflation_diagnostic.png"),
  combined,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(FIG_DIR, "zero_inflation_diagnostic.pdf"),
  combined,
  width = 12,
  height = 10
)

cli_alert_success("Saved: zero_inflation_diagnostic.png/.pdf")

# =============================================================================
# 10. Summary report
# =============================================================================

cli_h2("10. Summary Report")

cat("\n")
cat("===================================================================\n")
cat("ZERO-INFLATION ANALYSIS SUMMARY\n")
cat("===================================================================\n")
cat("\n")
cat("Matrix dimensions: ", n_corals, " corals × ", n_species, " species\n", sep = "")
cat("Overall sparsity:  ", round(prop_zeros * 100, 1), "% zeros\n", sep = "")
cat("\n")
cat("Species with >70% zeros: ", sum(species_sparsity$prop_zeros > 0.70),
    " (", round(pct_gt70, 1), "%)\n", sep = "")
cat("Species with >80% zeros: ", sum(species_sparsity$prop_zeros > 0.80),
    " (", round(pct_gt80, 1), "%)\n", sep = "")
cat("Species with >90% zeros: ", sum(species_sparsity$prop_zeros > 0.90),
    " (", round(pct_gt90, 1), "%)\n", sep = "")
cat("\n")
cat("PCA variance explained (PC1):\n")
for (i in 1:nrow(pca_summaries)) {
  cat("  ", pca_summaries$transformation[i], ": ",
      round(pca_summaries$pc1_var[i], 1), "%\n", sep = "")
}
cat("\n")
cat("NMDS stress (lower is better):\n")
for (i in 1:nrow(nmds_results)) {
  cat("  ", nmds_results$distance[i], ": ",
      round(nmds_results$stress, 3), " (", nmds_results$interpretation[i], ")\n", sep = "")
}
cat("\n")
cat("PC1 correlations with richness/abundance:\n")
for (i in 1:nrow(pc1_correlations)) {
  cat("  ", pc1_correlations$transformation[i],
      ": richness=", round(pc1_correlations$cor_richness[i], 2),
      ", abundance=", round(pc1_correlations$cor_abundance[i], 2), "\n", sep = "")
}
cat("\n")
cat("===================================================================\n")
cat("RECOMMENDATION:\n")
cat("===================================================================\n")
cat("\n")

# Decision logic
high_sparsity <- prop_zeros > 0.70
hell_pc1_good <- pca_summaries$pc1_var[pca_summaries$transformation == "hellinger"] > 15
pc1_confounded <- any(abs(pc1_correlations$cor_richness) > 0.7)

if (high_sparsity && pc1_confounded) {
  cat("⚠️  WARNING: High zero-inflation (", round(prop_zeros*100, 1), "%) detected.\n", sep = "")
  cat("⚠️  PC1 is confounded with richness/abundance.\n")
  cat("⚠️  RECOMMENDATION: Use NMDS instead of PCA, OR report both methods.\n")
} else if (high_sparsity) {
  cat("⚠️  Moderate zero-inflation (", round(prop_zeros*100, 1), "%) detected.\n", sep = "")
  cat("✅  Hellinger transformation appears appropriate for PCA.\n")
  cat("✅  RECOMMENDATION: Continue with PCA but mention sparsity in manuscript.\n")
} else {
  cat("✅  Zero-inflation is within acceptable range.\n")
  cat("✅  RECOMMENDATION: PCA is appropriate.\n")
}

cat("\n")
cat("===================================================================\n")
cat("\nAnalysis complete. Results saved to output/MRB/tables/\n\n")
