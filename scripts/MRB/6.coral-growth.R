# ==============================================================================
# MRB Analysis Script 6: Coral Colony Growth Analysis
# ==============================================================================
# Purpose: Analyze coral colony growth from 3D photogrammetry data (2019-2021).
#          Processes mesh measurements, computes growth metrics, and fits mixed-
#          effects models with reef random effects to test treatment effects on
#          growth. Includes ANCOVA, size-corrected growth, and surface area-scaled
#          analyses. Exports growth data for use in script 7 (physiology integration).
#
# Inputs:
#   - data/MRB Amount/MRB_2019_200K_mesh_measure.csv
#   - data/MRB Amount/MRB_May_2021_200K_mesh_measure.csv
#   - data/MRB Amount/1. amount_manual_colony_measurements_dec2019_and_may2021.xlsx
#   - data/MRB Amount/coral_id_position_treatment.csv
#
# Outputs:
#   Figures (output/MRB/figures/coral/):
#     - percent_alive_hist.png                              # Survival histogram
#     - percent_alive_by_treatment.png                      # Survival by treatment
#     - ANCOVA_Init_vs_Final_Volume_by_Treatment.png        # Growth trajectories
#     - SizeCorrected_Volume_Growth_by_Treatment.png        # Size-corrected growth
#     - SA_Scaled_Volume_Growth_by_Treatment.png            # Surface area scaled
#     - DeltaVolume_vs_SA_ANCOVA.png                        # ANCOVA plot
#   Tables (output/MRB/tables/):
#     - (ANCOVA and growth model results printed to console)
#   Cached objects (data/processed/):
#     - coral_growth.rds                                    # Growth data object
#     - coral_growth.csv                                    # Growth data table
#   Data exports (data/MRB Amount/):
#     - coral_growth_surface_area_change_filtered.csv       # Filtered growth data
#
# Depends:
#   R (>= 4.3), tidyverse, lme4, lmerTest, car, here, readxl, scales
#
# Run after:
#   - 1.libraries.R (loads required packages)
#   - utils.R (utility functions)
#   - mrb_figure_standards.R (plotting standards)
#
# Author: Adrian C. Stier / CAFI Team
# Created: 2019-2021 (original analysis)
# Last updated: 2025-11-05
#
# Reproducibility notes:
#   - set.seed(42) for all stochastic operations
#   - Mixed models use REML estimation (default in lmer)
#   - Alive threshold: ≥80% live tissue for inclusion
#   - All paths use here::here() for portability
#   - All plots are both printed to console and saved to disk
#   - Growth data exported for use in script 7 (coral physiology integration)
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Source centralized libraries, utilities, and standards
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")
source("scripts/MRB/mrb_figure_standards.R")  # For colors, themes, save functions

# Backward-compatible aliases (show_and_save, save_both) provided by utils.R

# ---- 1. Setup ----------------------------------------------------------------
set.seed(42)

# Libraries already loaded by 1.libraries.R - no need to reload

# ---- 1.1 Standard Colors & Theme ---------------------------------------------
# Use TREATMENT_COLORS from mrb_figure_standards.R:
#   "1" = "#E69F00" (Orange)
#   "3" = "#56B4E9" (Sky Blue)
#   "6" = "#009E73" (Green)
# Use theme_publication() from mrb_figure_standards.R
# Use save_figure() or save_both() from utils.R

# Default geom tweaks
update_geom_defaults("violin",   list(linewidth = 0.4))
update_geom_defaults("boxplot",  list(linewidth = 0.5))
update_geom_defaults("point",    list(size = 2))
update_geom_defaults("line",     list(linewidth = 0.8))

# strip_fe() and ALIVE_THRESH defined in utils.R - sourced above

# Output dirs
out_dir_fig  <- here("output", "MRB", "figures", "coral")
out_dir_phys <- file.path(out_dir_fig, "physio")
dir.create(out_dir_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_phys, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load Raw Data --------------------------------------------------------
cli::cli_h2("Loading raw mesh, treatment, & manual % alive data")

data_2019 <- read_csv(here("data/MRB Amount/MRB_2019_200K_mesh_measure.csv"),
                      show_col_types = FALSE)
data_2021 <- read_csv(here("data/MRB Amount/MRB_May_2021_200K_mesh_measure.csv"),
                      show_col_types = FALSE)

treatment_df_raw <- read_csv(here("data/MRB Amount/coral_id_position_treatment.csv"),
                             show_col_types = FALSE)

manual_file <- here("data", "MRB Amount",
                    "1. amount_manual_colony_measurements_dec2019_and_may2021.xlsx")

manual_alive <- read_excel(manual_file) %>%
  mutate(
    coral_id      = strip_fe(as.character(coral_id)),
    percent_alive = percent_alive_may21 / 100
  ) %>%
  dplyr::select(coral_id, percent_alive)

# ---- 3. Clean Treatment / Reef Lookup ---------------------------------------
cli::cli_h2("Building treatment & reef lookup")

treatment_df <- treatment_df_raw %>%
  mutate(
    coral_id  = strip_fe(as.character(coral_id)),
    treatment = as.factor(treatment),
    row_num   = str_extract(position, "^\\d+"),
    col_num   = str_extract(position, "(?<=-)\\d+"),
    reef      = paste0("Reef_", row_num, "-", col_num)
  ) %>%
  dplyr::select(coral_id, treatment, reef) %>%
  distinct()

stopifnot(all(c("coral_id","treatment","reef") %in% names(treatment_df)))

# ---- 4. Process 2019 Mesh Data ----------------------------------------------
cli::cli_h2("Processing 2019 mesh data")

data_2019 <- data_2019 %>%
  mutate(
    reef            = str_extract(Chunk, "^\\d+-\\w*"),
    project_code    = str_extract(Chunk, "\\(([^)]+)\\)") %>% str_remove_all("[()]"),
    coral_id        = str_extract(Chunk, "POC\\d+") %>% strip_fe(),
    coral_replicate = str_extract(reef, "[A-Z]$") %>% replace_na("A"),
    open_closed     = str_extract(Model, "_(open|closed)_") %>% str_replace_all("_", "")
  ) %>%
  rename_with(str_trim)

data_2019 <- data_2019 %>%
  group_by(coral_id, open_closed) %>%
  mutate(version_count = n_distinct(Version)) %>%
  ungroup() %>%
  filter(!(Version == 1 & version_count > 1)) %>%
  dplyr::select(-version_count)

data_2019 <- data_2019 %>%
  mutate(
    reef_number = str_extract(reef, "^\\d+-\\d+"),
    reef_letter = str_remove(reef, "^\\d+-\\d") %>% replace_na("A")
  )

open_19 <- data_2019 %>%
  filter(open_closed == "open") %>%
  dplyr::select(coral_id, surface_area_open = `Surface Area (cm^2)`)

closed_19 <- data_2019 %>% filter(open_closed == "closed")

coral_2019 <- closed_19 %>%
  left_join(open_19, by = "coral_id") %>%
  mutate(`Surface Area (cm^2)` = surface_area_open) %>%
  dplyr::select(-surface_area_open)

# ---- 5. Process 2021 Mesh Data ----------------------------------------------
cli::cli_h2("Processing 2021 mesh data")

data_2021 <- data_2021 %>%
  mutate(
    reef            = str_extract(Chunk, "^\\d+-\\w*"),
    project_code    = str_extract(Chunk, "\\(([^)]+)\\)") %>% str_remove_all("[()]"),
    coral_id        = str_extract(Chunk, "POC\\d+") %>% strip_fe(),
    coral_replicate = str_extract(reef, "[A-Z]$") %>% replace_na("A"),
    open_closed     = str_extract(Model, "_(open|closed)_") %>% str_replace_all("_", "")
  ) %>%
  distinct() %>%
  mutate(
    reef_number = str_extract(reef, "^\\d+-\\d+"),
    reef_letter = str_remove(reef, "^\\d+-\\d") %>% replace_na("A")
  )

open_21 <- data_2021 %>%
  filter(open_closed == "open") %>%
  dplyr::select(coral_id, surface_area_open = `Surface Area (cm^2)`)

open_21u <- open_21 %>%
  group_by(coral_id) %>%
  summarise(surface_area_open = mean(surface_area_open, na.rm = TRUE), .groups = "drop")

closed_21 <- data_2021 %>% filter(open_closed == "closed")

coral_2021 <- closed_21 %>%
  left_join(open_21u, by = "coral_id") %>%
  mutate(`Surface Area (cm^2)` = surface_area_open) %>%
  dplyr::select(-surface_area_open)

# ---- 6. Combine Years & Derive Metrics --------------------------------------
cli::cli_h2("Combining years & computing growth metrics")

coral_2019 <- coral_2019 %>% mutate(dataset_year = "2019") %>%
  dplyr::select(-Version, -Issues_with_model)
coral_2021 <- coral_2021 %>% mutate(dataset_year = "2021")

common_cols <- intersect(names(coral_2019), names(coral_2021))
coral_2019  <- coral_2019 %>% dplyr::select(dplyr::all_of(common_cols), dataset_year)
coral_2021  <- coral_2021 %>% dplyr::select(dplyr::all_of(common_cols), dataset_year)

coral_all <- bind_rows(coral_2019, coral_2021)

names(coral_all) <- c(
  "chunk","model","model_type","surface_area_cm2","area_cm2","max_height_cm",
  "min_height_cm","height_range_cm","volume_cm3","extent_volume_cm3",
  "convex_hull_volume_cm3","reef","project_code","coral_id","coral_replicate",
  "open_closed","reef_number","reef_letter","dataset_year"
)

coral_all <- coral_all %>%
  mutate(interstitial_space_cm3 = convex_hull_volume_cm3 - volume_cm3)

coral_wide <- coral_all %>%
  dplyr::select(coral_id, dataset_year,
                surface_area_cm2, volume_cm3, max_height_cm, interstitial_space_cm3) %>%
  pivot_wider(names_from = dataset_year,
              values_from = c(surface_area_cm2, volume_cm3, max_height_cm, interstitial_space_cm3),
              names_glue = "{.value}_{dataset_year}")

coral_changes <- coral_wide %>%
  mutate(
    delta_surface_area  = surface_area_cm2_2021          - surface_area_cm2_2019,
    delta_volume        = volume_cm3_2021                - volume_cm3_2019,
    delta_max_height    = max_height_cm_2021             - max_height_cm_2019,
    delta_interstitial  = interstitial_space_cm3_2021    - interstitial_space_cm3_2019
  )

# ---- 7. Attach Treatment/Reef & % Alive; Filter -----------------------------
cli::cli_h2("Joining treatment/reef/% alive & filtering colonies")

vol_df <- coral_wide %>%
  dplyr::select(coral_id, volume_cm3_2019, volume_cm3_2021) %>%
  rename(vol_2019 = volume_cm3_2019,
         vol_2021 = volume_cm3_2021)

meta_df <- treatment_df %>%
  left_join(manual_alive, by = "coral_id")

coral_treatment <- coral_changes %>%
  left_join(meta_df, by = "coral_id") %>%
  left_join(vol_df,  by = "coral_id") %>%
  mutate(
    treatment     = droplevels(treatment),
    percent_alive = coalesce(percent_alive, 0)
  )

missing_tr <- coral_treatment %>% filter(is.na(treatment) | is.na(reef))
if (nrow(missing_tr) > 0) cli::cli_alert_warning("{nrow(missing_tr)} rows missing treatment/reef info.")

# --- 7A. Visualize % alive & justify threshold -------------------------------
cli::cli_h2("Visualizing percent alive to justify 80% cutoff")

p_alive_hist <- ggplot(coral_treatment, aes(percent_alive)) +
  geom_histogram(binwidth = 0.05, fill = "#3182BD", color = "white") +
  geom_vline(xintercept = ALIVE_THRESH, linetype = "dashed", linewidth = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "% Alive distribution (May 2021)",
       x = "Proportion alive", y = "Count") +
  theme_publication(base_size = 13)

show_and_save(p_alive_hist, file.path(out_dir_fig, "percent_alive_hist.png"), w = 6, h = 4)

meta_trt <- meta_df %>%
  mutate(percent_alive_scaled = percent_alive,
         treatment = droplevels(treatment))

cli::cli_h2("meta_trt (percent_alive + treatment) -- full print")
print(meta_trt %>% arrange(coral_id), n = Inf)

p_alive_trt <- meta_trt %>%
  ggplot(aes(x = treatment, y = percent_alive_scaled, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8, colour = "gray30") +
  geom_hline(yintercept = ALIVE_THRESH, linetype = 2) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_treatment(guide = "none") +
  labs(title = "Percent Alive by Treatment (dashed = 80% cutoff)",
       x = "Treatment", y = "Percent Alive") +
  theme_publication(base_size = 13)

show_and_save(p_alive_trt,
              file.path(out_dir_fig, "percent_alive_by_treatment.png"),
              w = 7, h = 5)

keep_ids <- coral_treatment %>%
  filter(percent_alive >= ALIVE_THRESH) %>%
  pull(coral_id)

cli::cli_alert_info("{length(keep_ids)} colonies retained (>= {ALIVE_THRESH*100}% alive).")

coral_treatment_f <- coral_treatment %>% filter(coral_id %in% keep_ids)

# ============================================================================
# 8. Growth Analyses (Random Reef Effects)  [Filtered] ------------------------
# ============================================================================
cli::cli_h2("Section 8: Volume ANCOVA & size-corrected growth (filtered)")

ancova_dat <- coral_treatment_f %>%
  dplyr::select(coral_id, reef, treatment, vol_2019, vol_2021) %>%
  drop_na(vol_2019, vol_2021, treatment, reef) %>%
  mutate(
    log_init  = log(vol_2019),
    log_final = log(vol_2021)
  )

# Fit both models first for comparison
m_noint  <- lmer(log_final ~ log_init + treatment + (1 | reef), data = ancova_dat)
m_full   <- lmer(log_final ~ log_init * treatment + (1 | reef), data = ancova_dat)

# Test interaction
ancova_full  <- car::Anova(m_full,  type = 3); print(ancova_full)

# ---- 8A.1 Post-hoc tests: Which treatment slopes differ? --------------------
# Since the interaction is significant (p = 0.039), we need to identify which
# treatment pairs have different allometric slopes. We use emmeans to extract
# and compare slopes (i.e., the log_init coefficient) for each treatment.

cli::cli_h3("Post-hoc: Pairwise comparisons of allometric slopes by treatment")

# Extract slopes (log_init coefficient) for each treatment using emtrends
# emtrends estimates the slope of the continuous predictor (log_init) at each level of treatment
slope_trends <- emmeans::emtrends(m_full,
                                  pairwise ~ treatment,  # Compare slopes across treatments
                                  var = "log_init",       # Slope of log_init (allometric exponent)
                                  adjust = "tukey")       # Tukey HSD for multiple comparisons

cat("\n=== Treatment-specific allometric slopes (b estimates) ===\n")
print(summary(slope_trends$emtrends))

cat("\n=== Pairwise comparisons of slopes (Tukey-adjusted) ===\n")
slope_contrasts <- summary(slope_trends$contrasts)
print(slope_contrasts)

# Interpret the results
sig_pairs <- slope_contrasts[slope_contrasts$p.value < 0.05, ]
if (nrow(sig_pairs) > 0) {
  cli::cli_alert_success("Significant pairwise differences in allometric slopes:")
  for (i in 1:nrow(sig_pairs)) {
    cat(sprintf("  • %s: diff = %.4f, p = %.4f\n",
                sig_pairs$contrast[i],
                sig_pairs$estimate[i],
                sig_pairs$p.value[i]))
  }
} else {
  cli::cli_alert_info("No significant pairwise differences in slopes (all p > 0.05)")
  cli::cli_alert_info("Interaction may be driven by overall heterogeneity, not specific pairs")
}

# Also extract treatment-specific intercepts for completeness
cat("\n=== Treatment-specific intercepts (at mean log_init) ===\n")
intercept_means <- emmeans::emmeans(m_full, ~ treatment, at = list(log_init = mean(ancova_dat$log_init)))
print(summary(intercept_means))

# Test if we should use interaction model vs parallel slopes model
cat("\n=== Model comparison: Interaction vs Parallel Slopes ===\n")
anova_comparison <- anova(m_noint, m_full)
print(anova_comparison)
if (anova_comparison$`Pr(>Chisq)`[2] < 0.05) {
  cli::cli_alert_warning("Interaction model is significantly better (p < 0.05)")
  cli::cli_alert_info("Recommendation: Report treatment-specific slopes, not a single b")
} else {
  cli::cli_alert_success("Parallel slopes model is adequate")
}

ancova_noint <- car::Anova(m_noint, type = 3); print(ancova_noint)

newdat <- ancova_dat %>%
  group_by(treatment) %>%
  summarise(
    xmin = min(log_init, na.rm = TRUE),
    xmax = max(log_init, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(log_init = list(seq(xmin, xmax, length.out = 100))) %>%
  unnest(log_init) %>%
  mutate(
    log_pred = predict(m_noint,
                       newdata = tibble(log_init = log_init, treatment = treatment),
                       re.form = NA),
    init_vol = exp(log_init),
    pred_vol = exp(log_pred)
  )

p_ancova <- ggplot(ancova_dat, aes(x = vol_2019, y = vol_2021, colour = treatment)) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_line(data = newdat,
            aes(x = init_vol, y = pred_vol, colour = treatment),
            size = 1.1, show.legend = FALSE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_color_treatment() +
  labs(
    title    = "Initial vs Final Volume by Treatment",
    subtitle = "Mixed-effects ANCOVA with random reef intercept (filtered ≥80% alive)",
    x        = "Initial Volume 2019 (cm³)",
    y        = "Final Volume 2021 (cm³)"
  ) +
  theme_publication(base_size = 14)

show_and_save(p_ancova, file.path(out_dir_fig, "ANCOVA_Init_vs_Final_Volume_by_Treatment.png"))

# ---- 8A.2 Visualization: Treatment-specific slopes --------------------------
# Create figure showing treatment-specific allometric relationships
# This visualizes the significant interaction by plotting separate regression lines

cli::cli_h3("Creating figure with treatment-specific allometric slopes")

# Generate predictions from the interaction model (separate slopes per treatment)
newdat_interaction <- ancova_dat %>%
  group_by(treatment) %>%
  summarise(
    xmin = min(log_init, na.rm = TRUE),
    xmax = max(log_init, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(log_init = list(seq(xmin, xmax, length.out = 100))) %>%
  unnest(log_init) %>%
  mutate(
    log_pred_interaction = predict(m_full,
                                    newdata = tibble(log_init = log_init, treatment = treatment),
                                    re.form = NA),
    init_vol = exp(log_init),
    pred_vol_interaction = exp(log_pred_interaction)
  )

# Extract slope values for plot annotation
slope_summary <- summary(slope_trends$emtrends)
slope_labels <- slope_summary %>%
  as.data.frame() %>%
  mutate(
    label = sprintf("b = %.3f", log_init.trend),
    # Position labels at high x values for visibility
    x_pos = max(ancova_dat$vol_2019) * 0.6,
    y_pos = case_when(
      treatment == "1" ~ max(ancova_dat$vol_2021) * 0.15,
      treatment == "3" ~ max(ancova_dat$vol_2021) * 0.35,
      treatment == "6" ~ max(ancova_dat$vol_2021) * 0.60
    )
  )

p_ancova_interaction <- ggplot(ancova_dat, aes(x = vol_2019, y = vol_2021, colour = treatment)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_line(data = newdat_interaction,
            aes(x = init_vol, y = pred_vol_interaction, colour = treatment),
            linewidth = 1.3, alpha = 0.9) +
  geom_text(data = slope_labels,
            aes(x = x_pos, y = y_pos, label = label, colour = treatment),
            size = 4, fontface = "bold", show.legend = FALSE) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_color_treatment() +
  labs(
    title = "Treatment-Specific Allometric Relationships",
    subtitle = "Interaction model: log(Vf) ~ log(Vi) × treatment + (1|reef) [p = 0.039]",
    x = "Initial Volume 2019 (cm³)",
    y = "Final Volume 2021 (cm³)",
    caption = "Different slopes suggest coral density affects size-growth scaling"
  ) +
  theme_publication(base_size = 14) +
  theme(plot.caption = element_text(hjust = 0, face = "italic", size = 11))

show_and_save(p_ancova_interaction,
              file.path(out_dir_fig, "ANCOVA_TreatmentSpecific_Slopes.png"),
              w = 9, h = 7)

# ---- 8B. Size-corrected growth ----------------------------------------------
# UPDATED 2025-11-14: Following Craig Osenberg's advice to fit unified model
# with treatment included when estimating allometric exponent b.
#
# RATIONALE FOR UNIFIED MODEL APPROACH:
# - Pools data across all treatments (N=44) for most precise b estimate
# - More statistically rigorous than separate b per treatment
# - Creates comparable growth metric across all corals
# - Conservative approach: doesn't assume treatment-specific slopes are correct
# - Post-hoc tests (above) show no significant pairwise slope differences
# - Model comparison shows parallel slopes model is adequate (LRT p=0.073)
#
# DECISION: Use pooled b from unified model for ALL growth calculations
# The interaction (p=0.039) is acknowledged as exploratory finding, but
# we use the unified model's b for primary analyses to maintain comparability
# and statistical rigor across the pipeline.
#
# Unified model: log(Vf) = constant + b*log(Vi) + treatment + (1|reef)
m_unified <- lmer(log(vol_2021) ~ log(vol_2019) + treatment + (1 | reef), data = ancova_dat)

# Extract allometric exponent b from unified model
b_vol <- fixef(m_unified)[["log(vol_2019)"]]

cat("\n=== Unified Model Allometric Exponent ===\n")
cat(sprintf("Pooled b estimate (from unified model): %.4f\n", b_vol))
cat("This b value is used for ALL growth_vol_b calculations\n")
cat("(Treatment-specific b values shown in post-hoc tests above are for exploration only)\n\n")

# Calculate size-corrected growth using unified-model-derived b
# NOTE: This uses the SAME b for all corals regardless of treatment
# This ensures growth_vol_b is comparable across treatments
coral_treat2 <- ancova_dat %>%
  mutate(growth_vol_b = vol_2021 / vol_2019^b_vol)

# Test treatment effect on size-corrected growth
# NOTE: This can also be assessed directly from m_unified via likelihood ratio test
growth_mod <- lmer(growth_vol_b ~ treatment + (1 | reef), data = coral_treat2)
print(cat("\n--- Size-corrected growth model (growth_vol_b ~ treatment + (1|reef)) ---\n"))
print(car::Anova(growth_mod, type = 3))

# Also print the unified model results for comparison
print(cat("\n--- Unified allometric model (log(Vf) ~ log(Vi) + treatment + (1|reef)) ---\n"))
print(car::Anova(m_unified, type = 3))

p_growth_b <- ggplot(coral_treat2, aes(x = treatment, y = growth_vol_b, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, colour = "gray30") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "gold") +
  scale_fill_treatment(guide = "none") +
  labs(
    title    = expression("Size-corrected Volume Growth: " * m[f] / m[i]^b),
    subtitle = paste0("b = ", round(b_vol, 3), "; random reef intercept; filtered ≥80% alive"),
    x        = "Treatment",
    y        = expression(m[f] / m[i]^b)
  ) +
  theme_publication(base_size = 14)

show_and_save(p_growth_b, file.path(out_dir_fig, "SizeCorrected_Volume_Growth_by_Treatment.png"))

# ============================================================================
# 9. SA-Scaled Growth & ΔV ~ SA ANCOVA (Random Reef) [Filtered] --------------
# ============================================================================
cli::cli_h2("Section 9: SA-scaled growth & ΔV ~ SA ANCOVA (filtered)")

sec9_df <- coral_changes %>%
  dplyr::select(coral_id, delta_volume) %>%
  left_join(coral_wide %>% dplyr::select(coral_id, surface_area_cm2_2019,
                                         volume_cm3_2019, volume_cm3_2021),
            by = "coral_id") %>%
  left_join(meta_df %>% dplyr::select(coral_id, treatment, reef, percent_alive),
            by = "coral_id") %>%
  filter(coral_id %in% keep_ids) %>%
  drop_na(delta_volume, surface_area_cm2_2019, treatment, reef) %>%
  mutate(growth_sa = delta_volume / surface_area_cm2_2019)

m_growth_sa <- lmer(growth_sa ~ treatment + (1 | reef), data = sec9_df)
print(car::Anova(m_growth_sa, type = 3))

p_growth_sa <- ggplot(sec9_df, aes(x = treatment, y = growth_sa, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8, colour = "gray30") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "gold") +
  scale_fill_treatment(guide = "none") +
  labs(
    title    = "Surface-Area–Scaled Volume Growth",
    subtitle = expression("(V"[2021] * " - V"[2019] * ") / SA"[2019] * " (filtered ≥80% alive)"),
    x        = "Treatment",
    y        = expression(Delta*V/SA[2019])
  ) +
  theme_publication(base_size = 16)

show_and_save(p_growth_sa, file.path(out_dir_fig, "SA_Scaled_Volume_Growth_by_Treatment.png"))

m_dv_full <- lmer(delta_volume ~ surface_area_cm2_2019 * treatment + (1 | reef), data = sec9_df)
print(car::Anova(m_dv_full, type = 3))

m_dv_par  <- lmer(delta_volume ~ surface_area_cm2_2019 + treatment + (1 | reef), data = sec9_df)
print(car::Anova(m_dv_par, type = 3))

newdat_sa <- sec9_df %>%
  group_by(treatment) %>%
  summarise(
    xmin = min(surface_area_cm2_2019, na.rm = TRUE),
    xmax = max(surface_area_cm2_2019, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(surface_area_cm2_2019 = list(seq(xmin, xmax, length.out = 100))) %>%
  unnest(surface_area_cm2_2019) %>%
  mutate(
    pred_deltaV = predict(m_dv_par,
                          newdata = tibble(surface_area_cm2_2019 = surface_area_cm2_2019,
                                           treatment = treatment),
                          re.form = NA)
  )

p_deltaV_sa <- ggplot(sec9_df,
                      aes(x = surface_area_cm2_2019, y = delta_volume, colour = treatment)) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_line(data = newdat_sa,
            aes(x = surface_area_cm2_2019, y = pred_deltaV, colour = treatment),
            size = 1.2, show.legend = FALSE) +
  scale_color_treatment() +
  labs(
    title    = expression(Delta*V~" vs. Initial Surface Area"),
    subtitle = "Mixed-effects ANCOVA: ΔV ~ SA_2019 + treatment + (1|reef); filtered ≥80% alive",
    x        = expression("Surface Area 2019 (cm"^2*")"),
    y        = expression(Delta*Volume~" (cm"^3*")")
  ) +
  theme_publication(base_size = 16)

show_and_save(p_deltaV_sa, file.path(out_dir_fig, "DeltaVolume_vs_SA_ANCOVA.png"))

# ---- 10. Export tidy growth data --------------------------------------------
cli::cli_h2("Exporting tidy growth data")

coral_growth_df <- coral_changes %>%
  left_join(meta_df, by = "coral_id") %>%
  mutate(size_corrected_volume_growth = delta_volume / surface_area_cm2_2019) %>%
  filter(coral_id %in% keep_ids) %>%
  dplyr::select(coral_id, treatment, reef, percent_alive,
                delta_surface_area, delta_volume, delta_max_height, delta_interstitial,
                size_corrected_volume_growth)

write_csv(coral_growth_df,
          here("data/MRB Amount/coral_growth_surface_area_change_filtered.csv"))

# ---- 11. FINAL EXPORT OBJECT (growth) ---------------------------------------
coral_growth <- coral_changes %>%
  left_join(meta_df, by = "coral_id") %>%
  left_join(
    coral_treat2 %>% dplyr::select(coral_id, growth_vol_b),
    by = "coral_id"
  ) %>%
  left_join(
    sec9_df %>% dplyr::select(coral_id, growth_sa),
    by = "coral_id"
  ) %>%
  mutate(size_corrected_volume_growth = growth_sa) %>%
  filter(coral_id %in% keep_ids) %>%
  dplyr::select(
    coral_id, treatment, reef, percent_alive,
    delta_surface_area, delta_volume, delta_max_height, delta_interstitial,
    size_corrected_volume_growth,
    growth_vol_b
  )

dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(coral_growth, here("data/processed/coral_growth.rds"))
write_csv(coral_growth, here("data/processed/coral_growth.csv"))

cli::cli_alert_success("Exported coral_growth (RDS & CSV).")

