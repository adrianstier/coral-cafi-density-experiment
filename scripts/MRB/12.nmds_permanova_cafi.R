# ==============================================================================
# File: 12.nmds_permanova_cafi.R
# Purpose: CAFI-136 — Compare community composition among 3 habitat treatments
#          using PERMDISP + PERMANOVA across Bray (sqrt abund), Jaccard (incidence),
#          and Gower (relative abundance), with NMDS plots, gt tables, and
#          asymmetric sample-size analysis. No sourcing of 8; reads files directly.
# ==============================================================================


# Source libraries and utilities
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")
source("scripts/MRB/mrb_figure_standards.R")



# ------------------------------------------------------------------------------
# 0) Setup
# ------------------------------------------------------------------------------
# All packages loaded via source("scripts/MRB/1.libraries.R") above

set.seed(1234)
options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)

# Output directories - use proper MRB output structure
fig_dir   <- here::here("output", "MRB", "figures", "diversity")
table_dir <- here::here("output", "MRB", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Save helpers (use centralized functions) -------------------------------
# Legacy aliases (save_both, show_and_save) provided by utils.R

# ------------------------------------------------------------------------------
# 1) Load and prepare data (no sourcing)
#    Files (adjust if your paths differ):
#      - CAFI counts: data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv
#      - Growth:      data/processed/coral_growth.csv      (expects coral_id, treatment [or in physio])
#      - Physio:      output/MRB/figures/coral/physio/physio_metrics_plus_growth_filtered.csv
# ------------------------------------------------------------------------------
# strip_fe() defined in utils.R - sourced above

cafi_path   <- here::here("data", "MRB Amount", "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")
growth_path <- here::here("data", "processed", "coral_growth.csv")
physio_path <- here::here("output", "MRB", "figures", "coral", "physio",
                          "physio_metrics_plus_growth_filtered.csv")

if (!file.exists(cafi_path))   stop("Missing CAFI file at: ", cafi_path)
if (!file.exists(growth_path)) stop("Missing growth file at: ", growth_path)
if (!file.exists(physio_path)) stop("Missing physio file at: ", physio_path)

cli::cli_alert_info("Reading CAFI community file…")
cafi_raw <- readr::read_csv(cafi_path, show_col_types = FALSE) %>%
  mutate(coral_id = strip_fe(coral_id),
         species  = trimws(as.character(species))) %>%
  filter(!is.na(coral_id), !is.na(species))

cli::cli_alert_info("Building coral × species abundance matrix…")
comm_abund_wide <- cafi_raw %>%
  group_by(coral_id, species) %>%
  summarise(abundance = sum(count, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

# Bring treatment from growth (preferred) or physio as fallback
cli::cli_alert_info("Reading growth & physio for metadata…")
growth <- readr::read_csv(growth_path, show_col_types = FALSE) %>%
  mutate(coral_id = strip_fe(coral_id))
physio <- readr::read_csv(physio_path, show_col_types = FALSE) %>%
  mutate(coral_id = strip_fe(coral_id))

meta_src <- if ("treatment" %in% names(growth)) growth else physio
if (!"treatment" %in% names(meta_src)) {
  stop("Could not find a 'treatment' column in growth or physio files.")
}
meta_df <- meta_src %>%
  dplyr::select(coral_id, treatment) %>%
  distinct()

# Align to shared coral_ids
ids <- intersect(comm_abund_wide$coral_id, meta_df$coral_id)
comm_abund_wide <- comm_abund_wide %>% filter(coral_id %in% ids)
meta_df <- meta_df %>% filter(coral_id %in% ids)
meta_df <- meta_df[match(comm_abund_wide$coral_id, meta_df$coral_id), , drop = FALSE]

# ------------------------------------------------------------------------------
# 2) Filtering + matrix construction (square-root & z-scored variants)
#    NOTE: We keep all species present; if you want thresholds, add them here.
# ------------------------------------------------------------------------------
# Raw abundance matrix
comm_raw_mat <- comm_abund_wide %>%
  tibble::column_to_rownames("coral_id") %>%
  as.matrix()

# Drop zero-variance columns that can break things
nonzero_cols <- apply(comm_raw_mat, 2, function(x) var(x, na.rm = TRUE) > 0)
comm_raw_mat <- comm_raw_mat[, nonzero_cols, drop = FALSE]

# Square-root transform
comm_sqrt_mat <- sqrt(comm_raw_mat)

# Z-scored (center/scale) version of sqrt abundance — available if you want it
comm_sqrt_z <- scale(comm_sqrt_mat)

# Presence/absence (for Jaccard)
comm_pa_mat <- vegan::decostand(comm_raw_mat, method = "pa")

# Relative abundance (row-wise proportions) for Gower
row_sums <- pmax(rowSums(comm_raw_mat), .Machine$double.eps)
comm_rel_mat <- comm_raw_mat / row_sums

# Treatment factor (use standardized colors from mrb_figure_standards.R)
meta_df$treatment <- factor(meta_df$treatment, levels = c("1", "3", "6"))

cli::cli_alert_success("Matrices ready: RAW={ncol(comm_raw_mat)} spp | SQRT={ncol(comm_sqrt_mat)} spp | REL={ncol(comm_rel_mat)} spp")

# ------------------------------------------------------------------------------
# 3) Distance builder using the *requested* inputs
#    - bray   -> Bray-Curtis on SQRT abundance
#    - jaccard-> Jaccard on presence/absence
#    - gower  -> Gower on relative abundance
# ------------------------------------------------------------------------------
make_distance <- function(metric) {
  metric <- match.arg(metric, c("bray","jaccard","gower"))
  if (metric == "bray") {
    vegdist(comm_sqrt_mat, method = "bray")
  } else if (metric == "jaccard") {
    vegdist(comm_pa_mat, method = "jaccard", binary = TRUE)
  } else { # gower
    as.dist(cluster::daisy(comm_rel_mat, metric = "gower"))
  }
}

# ------------------------------------------------------------------------------
# 4) PERMDISP + adonis2 for one metric
# ------------------------------------------------------------------------------
analyze_metric <- function(metric) {
  dist_mat <- make_distance(metric)
  
  # PERMDISP
  permdisp_mod <- betadisper(dist_mat, meta_df$treatment)
  permdisp_res <- permutest(permdisp_mod, permutations = 999)
  permdisp_df <- tibble(
    metric = metric,
    test   = "PERMDISP",
    F      = permdisp_res$tab[1, "F"],
    df     = paste(permdisp_res$tab[1, "Df"], permdisp_res$tab[2, "Df"], sep = ", "),
    p      = permdisp_res$tab[1, "Pr(>F)"]
  )
  
  # PERMANOVA
  adonis_res <- adonis2(dist_mat ~ treatment, data = meta_df, permutations = 999)
  adonis_df <- tibble(
    metric = metric,
    test   = "adonis2",
    F      = adonis_res$F[1],
    df     = paste(adonis_res$Df[1], adonis_res$Df[2], sep = ", "),
    p      = adonis_res$`Pr(>F)`[1]
  )
  
  bind_rows(permdisp_df, adonis_df)
}

# ------------------------------------------------------------------------------
# 5) Run for all metrics & save gt table
# ------------------------------------------------------------------------------
metrics <- c("bray","jaccard","gower")
results_all <- purrr::map_dfr(metrics, analyze_metric)

gt_results <- results_all %>%
  arrange(metric, desc(test)) %>%
  gt() %>%
  fmt_number(columns = c(F, p), decimals = 3) %>%
  tab_header(title = "Community Composition Tests by Distance Metric")

gtsave(gt_results, file = file.path(table_dir, "community_tests.html"))

# ------------------------------------------------------------------------------
# 6) NMDS plots for each metric
# ------------------------------------------------------------------------------
plot_nmds <- function(metric) {
  dist_mat <- make_distance(metric)
  nmds <- metaMDS(dist_mat, k = 2, trymax = 100)
  nmds_df <- as.data.frame(nmds$points) %>%
    mutate(treatment = meta_df$treatment)
  
  ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = treatment)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_treatment() +
    theme_publication() +
    theme(legend.position = "top") +
    labs(title = glue("NMDS ({metric})"),
         subtitle = glue("Stress = {round(nmds$stress, 3)}"),
         color = "Treatment")
}

# save with show_and_save (PNG+PDF)
purrr::walk(metrics, function(m) {
  p <- plot_nmds(m)
  show_and_save(
    p,
    file.path(fig_dir, glue::glue("NMDS_{m}.png")),
    width = 6, height = 5, dpi = 600
  )
})

# ------------------------------------------------------------------------------
# 7) Asymmetric sample-size effect (balanced subsampling by treatment) — VISUALS
# ------------------------------------------------------------------------------
min_n <- min(table(meta_df$treatment))
metrics <- c("bray","jaccard","gower")

subsample_once <- function(metric) {
  subs_ids <- meta_df %>%
    group_by(treatment) %>%
    slice_sample(n = min_n) %>%
    ungroup() %>%
    pull(coral_id)
  
  keep_idx <- rownames(comm_raw_mat) %in% subs_ids
  
  dist_mat <- switch(metric,
                     bray    = vegdist(comm_sqrt_mat[keep_idx, , drop=FALSE], method = "bray"),
                     jaccard = vegdist(comm_pa_mat  [keep_idx, , drop=FALSE], method = "jaccard", binary = TRUE),
                     gower   = as.dist(cluster::daisy(comm_rel_mat[keep_idx, , drop=FALSE], metric = "gower"))
  )
  
  subs_meta <- meta_df[match(rownames(as.matrix(dist_mat)), meta_df$coral_id), , drop = FALSE]
  
  permdisp_p <- permutest(betadisper(dist_mat, subs_meta$treatment),
                          permutations = 999)$tab[1, "Pr(>F)"]
  adonis_p <- adonis2(dist_mat ~ treatment, data = subs_meta,
                      permutations = 999)$`Pr(>F)`[1]
  
  tibble(metric = metric, permdisp_p = permdisp_p, adonis_p = adonis_p)
}

cli::cli_alert_info("Running balanced subsampling (500 reps × 3 metrics)…")
subsample_res <- purrr::map_dfr(metrics, function(m) {
  purrr::map_dfr(1:500, ~subsample_once(m))
})

# Long form for plotting
subsample_long <- subsample_res %>%
  tidyr::pivot_longer(cols = c(permdisp_p, adonis_p),
                      names_to = "test", values_to = "p") %>%
  dplyr::mutate(test = dplyr::recode(test,
                                     permdisp_p = "PERMDISP",
                                     adonis_p   = "adonis2"))

# Summary table (kept)
subsample_summary <- subsample_long %>%
  dplyr::group_by(metric, test) %>%
  dplyr::summarise(
    mean_p = mean(p),
    prop_sig_0.05 = mean(p < 0.05),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    se = sqrt((prop_sig_0.05 * (1 - prop_sig_0.05)) / n),
    ci_lo = pmax(0, prop_sig_0.05 - 1.96 * se),
    ci_hi = pmin(1, prop_sig_0.05 + 1.96 * se)
  )

gtsave(
  subsample_summary %>%
    gt() %>%
    fmt_number(columns = c(mean_p, prop_sig_0.05, se, ci_lo, ci_hi), decimals = 3) %>%
    tab_header(title = "Balanced subsampling: p-value stability across metrics"),
  file = file.path(table_dir, "subsample_summary.html")
)

# ---- 7a) P-value distributions (facet by test × metric) ----------------------
p_pdist <- ggplot2::ggplot(subsample_long, ggplot2::aes(x = p)) +
  ggplot2::geom_histogram(bins = 30, boundary = 0, closed = "left") +
  ggplot2::geom_vline(xintercept = 0.05, linetype = 2) +
  ggplot2::facet_grid(test ~ metric) +
  ggplot2::labs(
    title = "Subsampling p-value distributions",
    x = "p-value",
    y = "Count"
  ) +
  theme_publication() +
  ggplot2::theme(legend.position = "none")

show_and_save(
  p_pdist,
  file.path(fig_dir, "subsample_pvalue_distributions.png"),
  width = 8.5, height = 6.5, dpi = 600
)

# ---- 7b) Proportion significant with 95% CI ---------------------------------
p_propsig <- ggplot2::ggplot(subsample_summary,
                             ggplot2::aes(x = metric, y = prop_sig_0.05)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
  ggplot2::facet_wrap(~ test) +
  ggplot2::labs(
    title = "Proportion of subsamples with p < 0.05",
    x = "Distance metric",
    y = "Proportion significant"
  ) +
  ggplot2::ylim(0, 1) +
  theme_publication() +
  ggplot2::theme(legend.position = "none")

show_and_save(
  p_propsig,
  file.path(fig_dir, "subsample_prop_sig.png"),
  width = 7, height = 5, dpi = 600
)

cli::cli_alert_success("Subsampling visuals saved (PNG+PDF) in: {fig_dir}")
