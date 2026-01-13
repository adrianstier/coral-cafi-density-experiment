# ==============================================================================
# MRB Analysis Script 7: Coral Physiology and Growth Integration
# ==============================================================================
# Purpose: Integrate coral physiology data with growth metrics from script 6.
#          Performs PCA ordination of physiological metrics (protein, carbohydrates,
#          zooxanthellae, AFDW) combined with growth data. Fits mixed-effects models
#          to test treatment effects on physiology and examines relationships between
#          physiology and growth using multivariate approaches.
#
# Inputs:
#   - data/MRB Amount/1. amount_master_phys_data_v5.csv       # Physiology data
#   - data/processed/coral_growth.csv                         # From script 6
#
# Outputs:
#   Figures (output/MRB/figures/coral/physio/):
#     - physio_growth_pairs.png                             # Pairwise correlations
#     - physio_growth_pca_scree.png                         # PCA scree plot
#     - physio_growth_pca_loadings.png                      # PCA loadings
#     - physio_growth_pca_biplot.png                        # PCA biplot
#     - pc1_loadings_and_scores_paired.png                  # PC1 analysis
#     - physio_by_treatment.png                             # Physio metrics by trt
#     - pca_scores_by_treatment.png                         # PCA scores by trt
#     - univariate_metric_by_treatment.png                  # Individual metrics
#   Tables (output/MRB/figures/coral/physio/):
#     - merged_physio_growth_filtered.csv                   # Merged dataset
#     - physio_metrics_plus_growth_filtered.csv             # Analysis-ready data
#     - mixed_model_fixed_effects.csv                       # LMM results
#     - mixed_model_typeIII_tests.csv                       # Type III tests
#     - treatment_effect_only_anova.csv                     # Treatment effects
#     - univariate_anova_treatment_effects.csv              # Univariate ANOVAs
#     - univariate_anova_treatment_effects.html/.png        # Formatted tables
#     - section18_treatment_effect_summary_all.csv          # Combined summary
#     - section18_treatment_effect_physio_pca.csv           # Physio/PCA effects
#     - section18_treatment_effect_growth.csv               # Growth effects
#
# Depends:
#   R (>= 4.3), tidyverse, lme4, lmerTest, car, broom, broom.mixed,
#   GGally, patchwork, gt, here
#
# Run after:
#   - 1.libraries.R (loads required packages)
#   - utils.R (utility functions)
#   - mrb_figure_standards.R (plotting standards)
#   - 6.coral-growth.R (generates coral_growth.csv)
#
# Author: Adrian C. Stier / CAFI Team
# Created: 2019-2021 (original analysis)
# Last updated: 2025-11-14
#
# Reproducibility notes:
#   - No explicit seed set (ordinations are deterministic)
#   - PCA uses scaled data (scale. = TRUE)
#   - Mixed models use REML estimation (default in lmer)
#   - Alive threshold: ≥80% live tissue (filtering from script 6)
#   - Post-hoc tests use Benjamini-Hochberg adjustment
#   - All paths use here::here() for portability
#
# Growth metric notes (UPDATED 2025-11-14):
#   - growth_vol_b is calculated using UNIFIED allometric model (script 6)
#   - Uses pooled b estimate across all treatments (not treatment-specific)
#   - This ensures growth metric is comparable across treatments
#   - See script 6 section 8B for rationale and implementation details
# ==============================================================================

# ==============================================================================
# SETUP
# ==============================================================================

# Source centralized libraries, utilities, and standards
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")
source("scripts/MRB/mrb_figure_standards.R")  # For colors, themes, save functions

# Backward-compatible aliases (show_and_save, save_both) provided by utils.R

# ---- Setup from script 6 (needed for compatibility) --------------------------
# strip_fe() and ALIVE_THRESH defined in utils.R - sourced above

# Output dirs
out_dir_fig  <- here("output", "MRB", "figures", "coral")
out_dir_phys <- file.path(out_dir_fig, "physio")
dir.create(out_dir_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_phys, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Load growth data from script 6
# ============================================================================

cli::cli_h2("Loading growth data from script 6")

# Load growth data from script 6
coral_growth <- read_csv(here("data/processed/coral_growth.csv"), show_col_types = FALSE)

# Extract keep_ids from coral_growth (these are the ≥80% alive corals)
keep_ids <- coral_growth$coral_id

# Create meta_df from coral_growth
meta_df <- coral_growth %>%
  dplyr::select(coral_id, treatment, reef, percent_alive)

cli::cli_h2("Section 13: Load & clean physiology data")

physio_file <- here("data", "MRB Amount", "1. amount_master_phys_data_v5.csv")
physio_df <- read_csv(physio_file, show_col_types = FALSE) %>%
  mutate(coral_id = strip_fe(as.character(coral_id)))

# ---- 14. Merge Physiology with Growth & Select Metrics ----------------------
cli::cli_h2("Section 14: Merge physiology with growth metrics (filtered set)")

cg_metrics <- coral_growth %>% dplyr::select(-treatment, -reef, -percent_alive)

coral_physio_df <- physio_df %>%
  left_join(meta_df,    by = "coral_id") %>%
  left_join(cg_metrics, by = "coral_id") %>%
  filter(coral_id %in% keep_ids) %>%
  mutate(
    treatment = fct_drop(factor(treatment)),
    reef      = factor(reef)
  )

physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")
growth_vars <- c("growth_vol_b")

physio_metrics_df <- coral_physio_df %>%
  dplyr::select(coral_id, reef, treatment, percent_alive,
                dplyr::all_of(physio_vars), dplyr::all_of(growth_vars))

cli::cli_h2("physio_metrics_df (filtered ≥80% alive) -- full print")
print(physio_metrics_df %>% arrange(coral_id), n = Inf)

write_csv(coral_physio_df,
          file.path(out_dir_phys, "merged_physio_growth_filtered.csv"))
write_csv(physio_metrics_df,
          file.path(out_dir_phys, "physio_metrics_plus_growth_filtered.csv"))

# ---- 15. Correlations & Pairwise Plots --------------------------------------
cli::cli_h2("Section 15: Correlation / pairs plot")

pairs_df <- physio_metrics_df %>%
  dplyr::select(dplyr::all_of(physio_vars), dplyr::all_of(growth_vars)) %>%
  mutate(across(everything(), as.numeric)) %>%
  drop_na()

panel_cor <- function(data, mapping, digits = 2, ...) {
  x <- data[[rlang::as_name(mapping$x)]]
  y <- data[[rlang::as_name(mapping$y)]]
  r <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
  ggplot() +
    annotate("text", x = .5, y = .5,
             label = formatC(r, format = "f", digits = digits), size = 4) +
    theme_void()
}
panel_smooth <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    theme_publication()
}
panel_diag <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_density(linewidth = 0.6, fill = "grey85") +
    theme_publication()
}

p_pairs <- GGally::ggpairs(
  pairs_df,
  lower = list(continuous = panel_smooth),
  diag  = list(continuous = panel_diag),
  upper = list(continuous = panel_cor),
  title = "Physiology & Growth Metrics: Pairwise Relationships (Filtered ≥80% alive)",
  progress = FALSE
)

show_and_save(p_pairs,
              file.path(out_dir_phys, "physio_growth_pairs.png"),
              w = 10, h = 10)

# ---- 16. PCA on Physio + Growth ---------------------------------------------
cli::cli_h2("Section 16: PCA of physiology + growth metrics")

pca_input <- physio_metrics_df %>%
  dplyr::select(dplyr::all_of(physio_vars), dplyr::all_of(growth_vars)) %>%
  drop_na()

pca_res <- prcomp(pca_input, scale. = TRUE)

anchor_var <- if ("size_corrected_volume_growth" %in% names(pca_input))
  "size_corrected_volume_growth" else growth_vars[1]

flip_sign <- sign(cor(pca_res$x[, 1],
                      pca_input[[anchor_var]],
                      use = "pairwise.complete.obs"))

scores_fl   <- pca_res$x
rot_fl      <- pca_res$rotation
scores_fl[,1] <- scores_fl[,1] * flip_sign
rot_fl[,1]    <- rot_fl[,1]    * flip_sign

var_exp  <- pca_res$sdev^2 / sum(pca_res$sdev^2)
scree_df <- tibble(PC = paste0("PC", seq_along(var_exp)),
                   variance = var_exp)

p_scree <- ggplot(scree_df, aes(PC, variance)) +
  geom_col() +
  geom_text(aes(label = round(variance, 2)),
            vjust = -0.4, size = 3) +
  labs(title = "Scree Plot: Physio + Growth PCA",
       y = "Proportion Variance Explained") +
  theme_publication()

show_and_save(p_scree,
              file.path(out_dir_phys, "physio_growth_pca_scree.png"),
              w = 6, h = 4)

loadings_df <- as_tibble(rot_fl, rownames = "metric") %>%
  dplyr::select(metric, PC1, PC2)

p_loadings <- loadings_df %>%
  pivot_longer(cols = c(PC1, PC2),
               names_to = "PC", values_to = "loading") %>%
  ggplot(aes(x = reorder(metric, loading), y = loading, fill = PC)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "PCA Loadings (PC1 & PC2): Physio + Growth",
       x = NULL, y = "Loading") +
  theme_publication()

show_and_save(p_loadings,
              file.path(out_dir_phys, "physio_growth_pca_loadings.png"),
              w = 8, h = 4)

physio_metrics_df <- physio_metrics_df %>%
  mutate(
    PC1_physio_growth = scores_fl[, 1],
    PC2_physio_growth = scores_fl[, 2]
  )

# ---- 16b. PCA biplot (PC1 vs PC2) with convex hulls & clean labels ----------
cli::cli_h3("Section 16b: PCA biplot")

score_df <- tibble(
  coral_id  = physio_metrics_df$coral_id,
  treatment = physio_metrics_df$treatment,
  PC1 = scores_fl[, 1],
  PC2 = scores_fl[, 2]
)

load_df <- as_tibble(rot_fl[, 1:2], rownames = "metric") %>%
  dplyr::select(metric, PC1, PC2) %>%
  mutate(
    label = dplyr::recode(metric,
                          protein_mg_cm2 = "Protein~(mg~cm^{-2})",
                          carb_mg_cm2    = "Carbohydrate~(mg~cm^{-2})",
                          zoox_cells_cm2 = "Zooxanthellae~(cells~cm^{-2})",
                          afdw_mg_cm2    = "Tissue~biomass~(AFDW~mg~cm^{-2})",
                          growth_vol_b   = "Growth",
                          .default       = metric)
  )

# scale arrows
arrow_mult <- 0.9 * max(abs(c(score_df$PC1, score_df$PC2))) /
  max(abs(c(load_df$PC1, load_df$PC2)))
load_df <- load_df %>%
  mutate(PC1 = PC1 * arrow_mult,
         PC2 = PC2 * arrow_mult)

# convex hulls with guard
hull_df <- score_df %>%
  group_by(treatment) %>%
  filter(n() >= 3) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

pc1_axis_lab <- paste0("PC1 (", scales::percent(var_exp[1], accuracy = 0.1), ")")
pc2_axis_lab <- paste0("PC2 (", scales::percent(var_exp[2], accuracy = 0.1), ")")

p_biplot <- ggplot(score_df, aes(PC1, PC2, colour = treatment)) +
  geom_polygon(data = hull_df,
               aes(PC1, PC2, fill = treatment, group = treatment),
               alpha = 0.12, colour = NA, show.legend = FALSE) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.02, "npc")),
               linewidth = 0.6, colour = "grey30",
               inherit.aes = FALSE) +
  geom_text(data = load_df,
            aes(x = PC1, y = PC2, label = label),
            parse = TRUE,
            size = 3.2, colour = "grey15", vjust = -0.35,
            inherit.aes = FALSE) +
  scale_color_treatment() +
  scale_fill_treatment(guide = "none") +
  labs(title = "PCA Biplot: Physiology & Growth Metrics",
       x = pc1_axis_lab, y = pc2_axis_lab) +
  theme_publication()

show_and_save(p_biplot,
              file.path(out_dir_phys, "physio_growth_pca_biplot.png"),
              w = 7.5, h = 6.5)

# ---- 16c. PC1 loadings + PC1 by treatment (paired figure) -------------------
cli::cli_h3("Section 16c: PC1 loadings + scores by treatment")

pretty_labs <- c(
  protein_mg_cm2  = "Protein~(mg~cm^{-2})",
  carb_mg_cm2     = "Carbohydrate~(mg~cm^{-2})",
  zoox_cells_cm2  = "Zooxanthellae~(cells~cm^{-2})",
  afdw_mg_cm2     = "Tissue~biomass~(AFDW~mg~cm^{-2})",
  growth_vol_b    = "Growth"
)

load_pc1 <- rot_fl %>%
  as.data.frame() %>%
  tibble::rownames_to_column("metric") %>%
  dplyr::select(metric, loading = PC1) %>%
  mutate(metric_lab = dplyr::recode(metric, !!!pretty_labs),
         metric_lab = factor(metric_lab, levels = metric_lab[order(loading)]))

p_load_pc1 <- ggplot(load_pc1, aes(x = loading, y = metric_lab)) +
  geom_segment(aes(x = 0, xend = loading, y = metric_lab, yend = metric_lab),
               linewidth = 0.7, colour = "grey55") +
  geom_point(size = 3, colour = "black") +
  scale_y_discrete(labels = scales::label_parse()) +
  labs(x = expression(bold("Loading on coral condition (PC"[1][coral]*")")), y = NULL) +
  theme_publication() +
  theme(
    axis.text.y  = element_text(size = 13),
    axis.text.x  = element_text(size = 12),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.margin = margin(5, 15, 5, 5)
  )


pc1_score_lab <- paste0("PC1 score (", scales::percent(var_exp[1], accuracy = 0.1), ")")

#THIS IS THE A PUBLICATION FIGURE

 
# Build the data for the PC1-by-treatment plot
pc1_df <- score_df %>%
  dplyr::select(treatment, PC1) %>%
  dplyr::mutate(treatment = forcats::fct_drop(factor(treatment),
                                              only = c("1","3","6"))) %>%
  dplyr::arrange(treatment)

p_pc1_trt <- ggplot(pc1_df, aes(x = treatment, y = PC1, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.22, fill = "white", outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.14, size = 2, alpha = 0.65, colour = "grey30") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               aes(fill = treatment), colour = "black") +
  scale_fill_treatment(guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Number of corals", y = pc1_score_lab) +
  theme_publication() +
  theme(
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.margin = margin(5, 5, 5, 15)
  )

paired_pc1 <- p_load_pc1 + p_pc1_trt +
  patchwork::plot_layout(widths = c(0.45, 0.55)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )

show_and_save(paired_pc1,
              file.path(out_dir_phys, "pc1_loadings_and_scores_paired.png"),
              w = 10, h = 5)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(dirname(dirname(out_dir_phys)), "publication-figures", "pc1_loadings_and_scores_paired.pdf"),
  paired_pc1, width = 10, height = 5, dpi = 600, bg = "white"
)
ggsave(
  file.path(dirname(dirname(out_dir_phys)), "publication-figures", "pc1_loadings_and_scores_paired.png"),
  paired_pc1, width = 10, height = 5, dpi = 600, bg = "white"
)

# ---- 17. Mixed Models: Physio & PCA scores ~ treatment + (1|reef) -----------
cli::cli_h2("Section 17: Mixed models (physio & PCA scores ~ treatment + (1|reef))")

# helper for tidy fixed effects
safe_tidy <- function(tt) {
  need <- c("term","estimate","conf.low","conf.high","statistic","p.value")
  miss <- setdiff(need, names(tt))
  if (length(miss)) tt[miss] <- NA_real_
  tt %>%
    filter(term != "(Intercept)") %>%
    dplyr::select(dplyr::all_of(need))
}

# single std_anova defined ONCE (used again in Section 18)
std_anova <- function(x){
  df <- as.data.frame(x) %>% tibble::rownames_to_column("term")
  need <- c("Df","F value","Chisq","Pr(>F)","Pr(>Chisq)")
  miss <- setdiff(need, names(df))
  if (length(miss)) df[miss] <- NA_real_
  df %>%
    mutate(
      statistic = coalesce(`F value`, Chisq),
      p.value   = coalesce(`Pr(>F)`, `Pr(>Chisq)`)
    ) %>%
    dplyr::select(term, Df, statistic, p.value)
}

metrics_all <- c(physio_vars, "PC1_physio_growth", "PC2_physio_growth")

fit_mm <- function(response, df) {
  form <- as.formula(paste(response, "~ treatment + (1 | reef)"))
  m    <- lmer(form, data = df,
               control = lmerControl(check.conv.singular = "ignore"))
  tibble(
    metric = response,
    model  = list(m),
    anova  = list(car::Anova(m, type = 3)),
    tidy   = list(safe_tidy(broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE))),
    glance = list(broom.mixed::glance(m))
  )
}

mm_results <- map_dfr(metrics_all, fit_mm, df = physio_metrics_df)

mm_table <- mm_results %>%
  unnest(tidy) %>%
  dplyr::select(metric, term, estimate, conf.low, conf.high, statistic, p.value)

write_csv(mm_table,
          file.path(out_dir_phys, "mixed_model_fixed_effects.csv"))

anova_table <- mm_results %>%
  mutate(anova_tbl = map(anova, std_anova)) %>%
  unnest(anova_tbl)

write_csv(anova_table,
          file.path(out_dir_phys, "mixed_model_typeIII_tests.csv"))

treat_tests <- anova_table %>%
  filter(term == "treatment") %>%
  mutate(
    p.BH = p.adjust(p.value, method = "BH"),
    sig  = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

cli::cli_h3("Treatment effect (Type III)")
print(treat_tests)

write_csv(treat_tests,
          file.path(out_dir_phys, "treatment_effect_only_anova.csv"))


#THIS IS THE A PUBLICATION FIGURE


# Plots
# Create more intuitive facet labels
metric_labels <- c(
  "afdw_mg_cm2" = "AFDW (mg/cm²)",
  "carb_mg_cm2" = "Carbohydrate (mg/cm²)",
  "protein_mg_cm2" = "Protein (mg/cm²)",
  "zoox_cells_cm2" = "Zooxanthellae (cells/cm²)"
)

p_physio_trt <- physio_metrics_df %>%
  pivot_longer(cols = dplyr::all_of(physio_vars),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(metric_labels), labels = metric_labels)) %>%
  ggplot(aes(x = treatment, y = value, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, colour = "gray30") +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_treatment(guide = "none") +
  labs(title = "Physiology Metrics by Treatment (≥80% alive)",
       x = "Treatment", y = "Value") +
  theme_publication(base_size = 13)

show_and_save(p_physio_trt,
              file.path(out_dir_phys, "physio_by_treatment.png"),
              w = 10, h = 8)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(dirname(dirname(out_dir_phys)), "publication-figures", "physio_by_treatment.pdf"),
  p_physio_trt, width = 10, height = 8, dpi = 600, bg = "white"
)
ggsave(
  file.path(dirname(dirname(out_dir_phys)), "publication-figures", "physio_by_treatment.png"),
  p_physio_trt, width = 10, height = 8, dpi = 600, bg = "white"
)

p_pca_trt <- physio_metrics_df %>%
  pivot_longer(cols = c(PC1_physio_growth, PC2_physio_growth),
               names_to = "PC",
               values_to = "score") %>%
  ggplot(aes(x = treatment, y = score, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, alpha = 0.95) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, colour = "gray30") +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_treatment(guide = "none") +
  labs(title = "PCA Scores (PC1 & PC2) by Treatment",
       x = "Treatment", y = "Score") +
  theme_publication(base_size = 13)

show_and_save(p_pca_trt,
              file.path(out_dir_phys, "pca_scores_by_treatment.png"),
              w = 8, h = 6)

cli::cli_alert_success("Section 17 complete: mixed models + plots exported.")


# ============================================================================
# 17b. Univariate models for each physiological/growth metric ~ treatment ----
# ============================================================================

cli::cli_h2("Section 17b: Univariate treatment effects for individual coral metrics")

# Define clean, publication-ready labels for output and figures
pretty_labs <- c(
  protein_mg_cm2  = "Protein~(mg~cm^{-2})",
  carb_mg_cm2     = "Carbohydrate~(mg~cm^{-2})",
  zoox_cells_cm2  = "Zooxanthellae~(cells~cm^{-2})",
  afdw_mg_cm2     = "Tissue~biomass~(AFDW~mg~cm^{-2})",
  growth_vol_b    = "Growth"
)

pretty_labels_gt <- c(
  protein_mg_cm2  = "Protein (mg cm⁻²)",
  carb_mg_cm2     = "Carbohydrate (mg cm⁻²)",
  zoox_cells_cm2  = "Zooxanthellae (cells cm⁻²)",
  afdw_mg_cm2     = "Tissue biomass (AFDW mg cm⁻²)",
  growth_vol_b    = "Growth"
)

# Prepare long-form data
univar_df <- physio_metrics_df %>%
  dplyr::select(coral_id, reef, treatment, all_of(physio_vars), growth_vol_b) %>%
  pivot_longer(cols = -c(coral_id, reef, treatment),
               names_to = "metric", values_to = "value") %>%
  drop_na() %>%
  mutate(metric_label = pretty_labs[metric])

# Fit LMMs and extract ANOVA tables
univar_anova <- univar_df %>%
  group_by(metric) %>%
  group_modify(~ {
    m <- lmer(value ~ treatment + (1 | reef), data = .x,
              control = lmerControl(check.conv.singular = "ignore"))
    aov_tbl <- std_anova(car::Anova(m, type = 3))
    aov_tbl <- aov_tbl %>% filter(term == "treatment")
    aov_tbl$metric <- unique(.x$metric)
    aov_tbl
  }) %>%
  ungroup() %>%
  relocate(metric, .before = term) %>%
  mutate(
    p.BH = p.adjust(p.value, method = "BH"),
    sig  = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    metric_pretty = pretty_labels_gt[metric]
  )

# Save raw table
write_csv(univar_anova, file.path(out_dir_phys, "univariate_anova_treatment_effects.csv"))


# ----------------------------------------------------------------------------
# Create publication-quality gt table
# ----------------------------------------------------------------------------
univar_gt <- univar_anova %>%
  dplyr::select(metric_pretty, Df, statistic, p.value, sig) %>%  # <- removed p.BH
  arrange(p.value) %>%
  gt::gt() %>%
  gt::tab_header(
    title = gt::md("**Treatment effects on individual coral metrics**"),
    subtitle = "Univariate linear mixed-effects models with reef as random intercept"
  ) %>%
  gt::cols_label(
    metric_pretty = "Response Variable",
    Df            = "df",
    statistic     = "F",
    p.value       = "p",
    sig           = ""
  ) %>%
  gt::fmt_number(columns = c(statistic, p.value), decimals = 3) %>%  # <- removed p.BH
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_body(columns = sig, rows = sig != "")
  )

# Save
gt::gtsave(univar_gt, file.path(out_dir_phys, "univariate_anova_treatment_effects.html"))
gt::gtsave(univar_gt, file.path(out_dir_phys, "univariate_anova_treatment_effects.png"))

# ----------------------------------------------------------------------------
# Visualization: Violin plots
# ----------------------------------------------------------------------------
# --- Prepare long-form data (use plain, publication labels for facets) -------
pretty_facet <- c(
  protein_mg_cm2  = "Protein (mg cm⁻²)",
  carb_mg_cm2     = "Carbohydrate (mg cm⁻²)",
  zoox_cells_cm2  = "Zooxanthellae (cells cm⁻²)",
  afdw_mg_cm2     = "Tissue biomass (AFDW mg cm⁻²)",
  growth_vol_b    = "Growth"
)

# desired facet order
facet_levels <- pretty_facet[c(physio_vars, "growth_vol_b")]

univar_df <- physio_metrics_df %>%
  dplyr::select(coral_id, reef, treatment, dplyr::all_of(physio_vars), growth_vol_b) %>%
  tidyr::pivot_longer(cols = -c(coral_id, reef, treatment),
                      names_to = "metric", values_to = "value") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    metric_label = factor(pretty_facet[metric], levels = facet_levels)
  )

# --- Visualization: Violin plots (no parsing needed) -------------------------
p_univar_violin <- ggplot(univar_df, aes(x = treatment, y = value, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, colour = "gray30") +
  facet_wrap(~ metric_label, scales = "free_y") +   # <- plain labels, nice superscripts
  scale_fill_treatment(guide = "none") +
  labs(x = "Number of corals", y = "Performance metric") +
  theme_publication(base_size = 13)

show_and_save(
  p_univar_violin,
  file.path(out_dir_phys, "univariate_metric_by_treatment.png"),
  w = 10, h = 8
)

cli::cli_alert_success("Section 17b complete: univariate ANOVA results and figure exported.")


# ---- 18. Summary table of treatment effects (physio, PCA & growth) ----------
cli::cli_h2("Section 18: Treatment-effect p-values (physio, PCA, growth)")

fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            scales::pvalue(p, accuracy = 0.001))

if (!exists("anova_table")) stop("anova_table not found (run Section 17 first).")

treat_tests_physio <- anova_table %>%
  filter(term == "treatment") %>%
  mutate(source = "physio/PCA")

growth_models <- list(
  log_final_volume        = if (exists("m_noint"))    car::Anova(m_noint,  type = 3) else NULL,
  growth_vol_b            = if (exists("growth_mod")) car::Anova(growth_mod, type = 3) else NULL,
  growth_sa               = if (exists("m_growth_sa"))car::Anova(m_growth_sa, type = 3) else NULL,
  delta_volume_par_ANCOVA = if (exists("m_dv_par"))   car::Anova(m_dv_par, type = 3) else NULL
)

treat_tests_growth <- imap_dfr(growth_models, function(aov_obj, lab){
  if (is.null(aov_obj)) return(NULL)
  std_anova(aov_obj) %>%
    filter(term == "treatment") %>%
    mutate(metric = lab, source = "growth")
})

summary_tbl <- bind_rows(treat_tests_physio, treat_tests_growth) %>%
  mutate(
    p.BH   = p.adjust(p.value, method = "BH"),
    sig    = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    p_fmt   = fmt_p(p.value),
    pBH_fmt = fmt_p(p.BH)
  ) %>%
  dplyr::select(source, metric, Df, statistic, p.value, p_fmt, p.BH, pBH_fmt, sig) %>%
  arrange(p.value)

cli::cli_h3("Overall treatment tests (Type III) across metrics")
print(summary_tbl)

write_csv(summary_tbl,
          file.path(out_dir_phys, "section18_treatment_effect_summary_all.csv"))
write_csv(treat_tests_physio,
          file.path(out_dir_phys, "section18_treatment_effect_physio_pca.csv"))
write_csv(treat_tests_growth,
          file.path(out_dir_phys, "section18_treatment_effect_growth.csv"))

cli::cli_alert_success("Section 18 complete: combined treatment-effect table saved & printed.")



