# Debug script to understand Trapezia serenei numbers
library(tidyverse)
library(here)

# Load data
treatment_df <- read_csv(here("data/MRB Amount/coral_id_position_treatment.csv"), show_col_types = FALSE)
MRBcafi_df <- read_csv(here("data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv"), show_col_types = FALSE)

# Process treatment data - same as main script
treatment_df <- treatment_df %>%
  mutate(
    row = str_extract(position, "^\\d+"),
    column = str_extract(position, "(?<=-)\\d+"),
    replicate = str_extract(position, "[A-Za-z]+")
  ) %>%
  mutate(across(c(row, column), as.integer)) %>%
  mutate(reef = paste0("Reef_", row, "-", column))

# Merge
MRBcafi_df <- MRBcafi_df %>% left_join(treatment_df, by = "coral_id")

# Create species matrix - each reef = each coral position
species_only_df <- MRBcafi_df %>%
  filter(!is.na(species) & species != "", !is.na(count))

species_abundance <- species_only_df %>%
  group_by(reef, species) %>%
  summarise(total_abundance = sum(count, na.rm = TRUE), .groups = "drop")

species_matrix <- species_abundance %>%
  pivot_wider(names_from = species, values_from = total_abundance, values_fill = 0)

# Metadata
metadata <- MRBcafi_df %>%
  distinct(reef, treatment) %>%
  right_join(species_matrix %>% select(reef), by = "reef") %>%
  mutate(reef = factor(reef), treatment = factor(treatment, levels = c("1","3","6")))

# Get Trapezia serenei data
cat("=== Trapezia serenei Analysis ===\n\n")

# Long format
species_long <- species_matrix %>%
  pivot_longer(-reef, names_to = "species", values_to = "abundance") %>%
  left_join(metadata, by = "reef")

# Filter for Trapezia serenei
trap_data <- species_long %>%
  filter(species == "Trapezia serenei")

cat("Sample sizes by treatment:\n")
trap_data %>% group_by(treatment) %>% summarise(n = n()) %>% print()

cat("\nRaw abundance per reef (sample unit) by treatment:\n")

cat("\nTreatment 1 reef abundances:\n")
t1_vals <- trap_data %>% filter(treatment == 1) %>% pull(abundance)
print(t1_vals)
cat("Mean:", mean(t1_vals), "  SD:", sd(t1_vals), "\n")

cat("\nTreatment 3 reef abundances:\n")
t3_vals <- trap_data %>% filter(treatment == 3) %>% pull(abundance)
print(t3_vals)
cat("Mean:", mean(t3_vals), "  SD:", sd(t3_vals), "\n")

cat("\nTreatment 6 reef abundances:\n")
t6_vals <- trap_data %>% filter(treatment == 6) %>% pull(abundance)
print(t6_vals)
cat("Mean:", mean(t6_vals), "  SD:", sd(t6_vals), "\n")

# Bootstrap expected - same function as main script
np_sum_ci <- function(x, k, B = 10000L, probs = c(0.025, 0.5, 0.975), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L || k < 1L || B < 1L) {
    return(setNames(rep(NA_real_, length(probs)), paste0("p", probs)))
  }
  draw <- sample.int(n, size = k * B, replace = TRUE)
  sm <- colSums(matrix(x[draw], nrow = k))
  stats::quantile(sm, probs = probs, names = TRUE, type = 7)
}

set.seed(1234)
q3 <- np_sum_ci(t1_vals, k = 3, B = 10000, seed = 1234)
q6 <- np_sum_ci(t1_vals, k = 6, B = 10000, seed = 1234)

cat("\n=== Bootstrap Expected Values ===\n")
cat("Expected at k=3 (sampling 3 reefs from T1):\n")
cat("  Lower (2.5%):", q3[[1]], "\n")
cat("  Median:", q3[[2]], "\n")
cat("  Upper (97.5%):", q3[[3]], "\n")

cat("\nExpected at k=6 (sampling 6 reefs from T1):\n")
cat("  Lower (2.5%):", q6[[1]], "\n")
cat("  Median:", q6[[2]], "\n")
cat("  Upper (97.5%):", q6[[3]], "\n")

cat("\n=== Comparison (what figure shows) ===\n")
cat("Observed MEAN at T=3:", mean(t3_vals), "\n")
cat("Expected CI at k=3: [", q3[[1]], ",", q3[[3]], "]\n")
cat("\nObserved MEAN at T=6:", mean(t6_vals), "\n")
cat("Expected CI at k=6: [", q6[[1]], ",", q6[[3]], "]\n")

cat("\n=== Also check Calcinus latens (Panel H) and Dascyllus flavicaudus (Panel J) ===\n\n")

# Calcinus latens
calc_data <- species_long %>% filter(species == "Calcinus latens")
t1_calc <- calc_data %>% filter(treatment == 1) %>% pull(abundance)
t6_calc <- calc_data %>% filter(treatment == 6) %>% pull(abundance)

q6_calc <- np_sum_ci(t1_calc, k = 6, B = 10000, seed = 1234)

cat("Calcinus latens:\n")
cat("  T1 values:", t1_calc, "\n")
cat("  T1 mean:", mean(t1_calc), "\n")
cat("  T6 mean:", mean(t6_calc), "\n")
cat("  Expected k=6 CI: [", q6_calc[[1]], ",", q6_calc[[3]], "]\n\n")

# Dascyllus flavicaudus
dasc_data <- species_long %>% filter(species == "Dascyllus flavicaudus")
t1_dasc <- dasc_data %>% filter(treatment == 1) %>% pull(abundance)
t6_dasc <- dasc_data %>% filter(treatment == 6) %>% pull(abundance)

q6_dasc <- np_sum_ci(t1_dasc, k = 6, B = 10000, seed = 1234)

cat("Dascyllus flavicaudus:\n")
cat("  T1 values:", t1_dasc, "\n")
cat("  T1 mean:", mean(t1_dasc), "\n")
cat("  T6 mean:", mean(t6_dasc), "\n")
cat("  Expected k=6 CI: [", q6_dasc[[1]], ",", q6_dasc[[3]], "]\n")
