# ==============================================================================
# File: 05_MRB_community_abundance_analysis.R
# Purpose: Build species/community matrices, quantify observed vs. expected totals,
#          richness metrics, and output publication-ready figures for MRB CAFI data.
# Inputs:  data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv
#          data/MRB Amount/coral_id_position_treatment.csv
# Outputs: output/MRB/figures/abundance/*.png  (plots)
#          output/MRB/data/top_species_treatment_data.csv (tidy top-species data)
#          output/MRB/objects/mrb_comm_summaries.rds       (key R objects)
# Depends: R (>= 4.3), tidyverse, here, stringr, vegan, patchwork, cli
# Run after: 01_download_raw.R, 02_clean_merge.R (or ensure raw files are present)
# Author: Adrian C. Stier
# Date: 2025-07-22
# Repro notes: set.seed() used for all bootstrapping
# ==============================================================================

# ==============================================================================
# 0. SETUP / HOUSEKEEPING
# ==============================================================================

# Source centralized figure standards (for theme and save functions)
source("scripts/MRB/mrb_figure_standards.R")

## 0.1 Reproducibility settings ----------------------------------------------
set.seed(1234)                           # for any stochastic procedures
options(stringsAsFactors = FALSE, scipen = 999)  # avoid factors and sci. notation

## 0.2 Required packages ------------------------------------------------------
required_pkgs <- c(
  "here","dplyr","tidyr","readr","ggplot2","patchwork",
  "vegan","gt","tibble","stringr","lme4","broom.mixed","car",
  "purrr","lmerTest","cli","fitdistrplus",
  "performance","emmeans","MASS","broom","fs","forcats"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Please install missing packages: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(required_pkgs, library, character.only = TRUE))

## 0.3 Paths ------------------------------------------------------------------
DATA_DIR   <- here("data", "MRB Amount")
OUTPUT_DIR <- here("output", "MRB")

# ---- 0.4 Parameters ---------------------------------------------------------
params <- list(
  bootstrap_B    = 10000,   # iterations for bootstrap
  sig_level      = 0.05,    # alpha for t-based CIs at species level
  top_n_species  = 20,      # how many species panels for the "top 20" plot
  target_orders  = c("Decapoda", "Perciformes", "Neogastropoda"), # focus orders
  fig_dpi        = 600,     # publication quality (updated from 300)
  fig_bg         = "white",
  fig_comm_size   = c(8, 5.5),   # single community plots
  fig_top20_size  = c(10, 7.5),  # top-20 species panel
  fig_threepanel  = c(6, 10.5),  # three-panel stacked figure
  fig_orders_size = c(11, 8),    # order-level figures (landscape fit)
  seed            = 1234L,       # unified seed for bootstraps
  col_obs         = "black",     # observed mean/line color
  col_exp         = "gray30",    # expected/bootstrap line color (dark gray)
  ribbon_alpha    = 0.15,        # alpha for gray ribbons
  ribbon_color    = "gray60"     # gray for confidence ribbons
)

# ---- 0.5 Paths --------------------------------------------------------------
paths <- list(
  data_dir       = here("data", "MRB Amount"),
  data_cafi      = here("data", "MRB Amount", "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv"),
  data_treatment = here("data", "MRB Amount", "coral_id_position_treatment.csv"),
  out_root       = here("output", "MRB"),
  out_dir_fig    = here("output", "MRB", "figures", "abundance"),
  out_dir_data   = here("output", "MRB", "data"),
  out_dir_obj    = here("output", "MRB", "objects"),
  out_dir_reports= here("output", "MRB", "reports")  # PDFs/gt tables, etc.
)

# Ensure output directories exist
invisible(lapply(paths[c("out_root","out_dir_fig","out_dir_data","out_dir_obj","out_dir_reports")],
                 dir.create, recursive = TRUE, showWarnings = FALSE))

# ==== 1. Load Data ------------------------------------------------------------

cli::cli_h2("Loading data")

# Quick helpers
require_cols <- function(df, cols, name) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop(sprintf("[%s] missing required columns: %s", name, paste(missing, collapse = ", ")))
  }
}

# Nonparametric k-of-n-with-replacement sum CI (fast, vectorized)
np_sum_ci <- function(x, k, B = 10000L, probs = c(0.025, 0.5, 0.975), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L || k < 1L || B < 1L) {
    return(setNames(rep(NA_real_, length(probs)), paste0("p", probs)))
  }
  draw <- sample.int(n, size = k * B, replace = TRUE)     # indices
  sm   <- colSums(matrix(x[draw], nrow = k))              # k x B → col sums
  stats::quantile(sm, probs = probs, names = TRUE, type = 7)
}

# CAFI observations (counts per coral_id/species)
MRBcafi_df <- readr::read_csv(paths$data_cafi, show_col_types = FALSE, progress = FALSE)

# Ensure expected columns exist
require_cols(MRBcafi_df, c("coral_id", "species", "count"), "CAFI data")

# Taxonomy guard: ensure taxonomic fields exist
for (nm in c("order","family","class")) {
  if (!nm %in% names(MRBcafi_df)) MRBcafi_df[[nm]] <- NA_character_
}

MRBcafi_df <- MRBcafi_df %>%
  mutate(
    coral_id = as.factor(coral_id),
    site     = if ("site" %in% names(.)) as.factor(site) else NULL,
    count    = as.numeric(count)  # coerce in case CSV typed it oddly
  ) %>%
  filter(!is.na(coral_id))  # defensive: drop incomplete keys

# Treatment metadata (placement, treatment assignment)
treatment_df <- readr::read_csv(paths$data_treatment, show_col_types = FALSE, progress = FALSE)
require_cols(treatment_df, c("coral_id", "position", "treatment"), "Treatment data")

# ==== 2. Process Treatment Data & Merge ---------------------------------------

cli::cli_h2("Processing treatment metadata")
# Extract layout info (row/column/replicate) from position code and build reef label.
# Example position format assumed: "12-4A" -> row=12, column=4, replicate=A
treatment_df <- treatment_df %>%
  mutate(
    row       = str_extract(position, "^\\d+"),
    column    = str_extract(position, "(?<=-)\\d+"),
    replicate = str_extract(position, "[A-Za-z]+")
  ) %>%
  mutate(across(c(row, column), as.integer)) %>%
  mutate(reef = paste0("Reef_", row, "-", column))

# Attach treatment info to observation rows
MRBcafi_df <- MRBcafi_df %>% left_join(treatment_df, by = "coral_id")

# Check if any rows lack treatment (useful to catch merging problems)
missing_treatment <- MRBcafi_df %>% filter(is.na(treatment)) %>% count()
cli::cli_alert_info("Rows missing treatment: {missing_treatment$n}")

# ==== 3. Species Matrices + Species×Treatment Summaries =======================

cli::cli_h2("Building species matrix & summarising species abundance")

# 3.0 Guardrails ---------------------------------------------------------------
MRBcafi_df <- MRBcafi_df %>%
  mutate(
    species = if ("species" %in% names(.)) trimws(species) else species,
    count   = as.numeric(count)
  )

# keep only identified species and non-missing counts
species_only_df <- MRBcafi_df %>%
  filter(!is.na(species) & species != "", !is.na(count)) %>%
  mutate(count = ifelse(count < 0, NA_real_, count))  # forbid negatives

neg_counts <- sum(is.na(species_only_df$count))
if (neg_counts > 0) {
  cli::cli_warn("{neg_counts} record(s) had negative counts; set to NA and dropped in sums.")
}

# 3.1 Aggregated long table: reef × species totals -----------------------------
species_abundance <- species_only_df %>%
  group_by(reef, species) %>%
  summarise(total_abundance = sum(count, na.rm = TRUE), .groups = "drop")

# 3.2 Wide species matrix (reef rows, species columns) -------------------------
species_matrix <- species_abundance %>%
  tidyr::pivot_wider(
    names_from  = species,
    values_from = total_abundance,
    values_fill = 0
  ) %>%
  relocate(reef) %>%
  distinct(reef, .keep_all = TRUE)

# 3.3 Metadata: one row per reef with treatment --------------------------------
metadata <- MRBcafi_df %>%
  distinct(reef, treatment) %>%
  right_join(species_matrix %>% dplyr::select(reef), by = "reef") %>%
  mutate(
    reef      = factor(reef),
    treatment = factor(treatment, levels = c("1","3","6"))
  ) %>%
  arrange(reef)

stopifnot(nrow(metadata) == nrow(species_matrix))
stopifnot(!any(duplicated(metadata$reef)))

# 3.4 Community matrix (numeric; rownames = reef) ------------------------------
community_matrix <- species_matrix %>%
  arrange(reef) %>%
  dplyr::select(-reef) %>%
  as.data.frame()
rownames(community_matrix) <- as.character(metadata$reef)

if (anyNA(community_matrix)) {
  cli::cli_warn("NAs found in community_matrix; converting to 0.")
  community_matrix[is.na(community_matrix)] <- 0
}

# 3.5 Diagnostics: zero-only rows/columns --------------------------------------
row_sums <- rowSums(community_matrix)
col_sums <- colSums(community_matrix)
n_zero_reefs   <- sum(row_sums == 0)
n_zero_species <- sum(col_sums == 0)
if (n_zero_reefs > 0)   cli::cli_warn("{n_zero_reefs} reef(s) have total abundance of 0.")
if (n_zero_species > 0) cli::cli_inform("{n_zero_species} species columns are all zeros.")

# 3.6 Species × Treatment summaries --------------------------------------------
species_long_by_reef <- species_matrix %>%
  arrange(reef) %>%
  tidyr::pivot_longer(-reef, names_to = "species", values_to = "abundance") %>%
  left_join(metadata, by = "reef")

species_abundance_treatment <- species_long_by_reef %>%
  group_by(treatment, species) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance   = sd(abundance,   na.rm = TRUE),
    n              = sum(!is.na(abundance)),
    .groups        = "drop"
  ) %>%
  mutate(se = sd_abundance / sqrt(n))

cli::cli_alert_success("Built species_matrix, community_matrix, metadata, and species_abundance_treatment.")

# ==== 4. Statistical tests: total community abundance =========================

cli::cli_h2("Model test: Treatment effect on total community abundance")

totals_df <- tibble::tibble(
  reef       = rownames(community_matrix),
  total_abun = rowSums(community_matrix)
) %>%
  dplyr::left_join(metadata, by = "reef") %>%
  dplyr::filter(!is.na(treatment)) %>%
  dplyr::mutate(
    treatment     = factor(treatment, levels = c("1","3","6")),
    treatment_num = as.numeric(as.character(treatment))
  )

m_pois <- glm(total_abun ~ treatment, family = poisson, data = totals_df)

disp_res <- performance::check_overdispersion(m_pois)
cli::cli_inform(sprintf("Overdispersion check: ratio = %.2f, p = %.3g",
                        disp_res$dispersion_ratio, disp_res$p))
use_nb   <- isTRUE(disp_res$dispersion_ratio > 1.2 && disp_res$p < 0.05)

if (use_nb) {
  cli::cli_inform(sprintf(
    "Overdispersion detected (ratio = %.2f, p = %.3g); using negative-binomial GLM.",
    disp_res$dispersion_ratio, disp_res$p
  ))
  m_final   <- MASS::glm.nb(total_abun ~ treatment, data = totals_df)
  fam_label <- "negative binomial"
} else {
  cli::cli_inform("No meaningful overdispersion; keeping Poisson GLM.")
  m_final   <- m_pois
  fam_label <- "poisson"
}

omni_tbl <- car::Anova(m_final, type = 3) %>% as.data.frame() %>% tibble::rownames_to_column("term")

coef_tbl <- broom::tidy(m_final, conf.int = TRUE) %>%
  dplyr::mutate(
    exp_est  = exp(estimate),
    exp_low  = exp(conf.low),
    exp_high = exp(conf.high)
  )

emm_trt  <- emmeans::emmeans(m_final, ~ treatment, type = "response")
emm_tbl  <- broom::tidy(emm_trt)

pairs_tbl <- emmeans::contrast(emm_trt, method = "pairwise") %>%
  summary(type = "response", infer = c(TRUE, TRUE)) %>%
  tibble::as_tibble()

cat("\n--- Community total abundance model (family:", fam_label, ") ---\n")
print(summary(m_final))
cat("\nOmnibus treatment effect:\n"); print(omni_tbl)
cat("\nEstimated means (count scale):\n"); print(emm_tbl)
cat("\nPairwise contrasts (count scale):\n"); print(pairs_tbl)

model_total_abundance <- list(
  family   = fam_label,
  poisson  = m_pois,
  final    = m_final,
  overdisp = disp_res,
  omnitab  = omni_tbl,
  coefs    = coef_tbl,
  emmeans  = emm_tbl,
  pairs    = pairs_tbl
)

# ==== 5. Total community: observed vs. nonparametric expected =================

cli::cli_h2("Total community: observed vs. bootstrap-expected")

# 5.0 single-reef totals at treatment 1
comm_reps <- species_matrix %>%
  left_join(metadata, by = "reef") %>%
  filter(treatment == 1) %>%
  transmute(replicate_total = rowSums(across(-c(reef, treatment))))

x1 <- comm_reps$replicate_total
B  <- params$bootstrap_B
set.seed(params$seed)

# 5.1 Observed (sum of species means; SEs combined by sqrt(sum(se^2)))
obs_comm_df <- species_abundance_treatment %>%
  group_by(treatment) %>%
  summarise(
    total_abundance = sum(mean_abundance),
    se_sum          = sqrt(sum(se^2)),
    .groups         = "drop"
  ) %>%
  mutate(
    treatment = as.integer(as.character(treatment)),
    lower_obs = total_abundance - 1.96 * se_sum,
    upper_obs = total_abundance + 1.96 * se_sum
  ) %>%
  filter(treatment %in% c(1,3,6))

# 5.2 NP bootstrap (k = 3,6) via single helper
q3 <- np_sum_ci(x1, 3, B = B, seed = params$seed)
q6 <- np_sum_ci(x1, 6, B = B, seed = params$seed)
np_ci_comm <- tibble(
  treatment = c(3,6),
  lower_np  = c(q3[[1]], q6[[1]]),
  median_np = c(q3[[2]], q6[[2]]),
  upper_np  = c(q3[[3]], q6[[3]])
)

# 5.3 Plot with raw data points
# Get raw data points for each treatment
raw_totals <- totals_df %>%
  mutate(treatment = as.integer(as.character(treatment))) %>%
  filter(treatment %in% c(1, 3, 6))

p_obs_np2 <- ggplot(obs_comm_df, aes(x = treatment, y = total_abundance)) +
  # Add raw data points with jitter and alpha
  geom_jitter(data = raw_totals,
              aes(x = treatment, y = total_abun, fill = factor(treatment)),
              width = 0.15, height = 0, alpha = 0.4, size = 2.5,
              shape = 21, color = "white", stroke = 0.5) +
  # Use treatment colors for raw points
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  # Black line and points for observed mean
  geom_line(color = params$col_obs, linewidth = 0.8) +
  geom_point(shape = 21, fill = params$col_obs, color = "white", size = 3) +
  # Expected values - gray ribbon, red dashed line
  geom_ribbon(data = np_ci_comm,
              aes(x = treatment, ymin = lower_np, ymax = upper_np),
              inherit.aes = FALSE,
              fill = params$ribbon_color, alpha = params$ribbon_alpha) +
  geom_line(data = np_ci_comm,
            aes(x = treatment, y = median_np, group = 1),
            inherit.aes = FALSE,
            color = params$col_exp, linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(breaks = c(1,3,6)) +
  labs(x = "Treatment (Number of Corals)", y = "Total Abundance") +
  theme_publication()

ggsave(
  file.path(paths$out_dir_fig, "community_abundance_observed_bootstrap.png"),
  p_obs_np2,
  width  = params$fig_comm_size[1],
  height = params$fig_comm_size[2],
  dpi    = params$fig_dpi,
  bg     = params$fig_bg
)

# ==== 6. Top species per order: observed vs. NP bootstrap =====================

cli::cli_h2("Focal orders: top species scaling (observed vs. bootstrap)")

# 6.1 Identify top-4 species per focal order (pooled over all treatments)
top_species <- MRBcafi_df %>%
  filter(order %in% params$target_orders, !is.na(species)) %>%
  group_by(order, species) %>%
  summarise(total_count = sum(count, na.rm=TRUE), .groups="drop") %>%
  group_by(order) %>%
  slice_max(total_count, n=4) %>%
  pull(species)

# 6.2 All treatment counts per species for raw data (includes zeros)
species_reps <- species_matrix %>%
  tidyr::pivot_longer(-reef, names_to="species", values_to="count") %>%
  left_join(metadata, by="reef") %>%
  filter(treatment %in% c(1, 3, 6), species %in% top_species)

# 6.2b t = 1 counts per species for bootstrap (includes zeros)
species_reps_t1 <- species_reps %>%
  filter(treatment == 1)

# 6.3 Observed summary by species & treatment (1, 3, 6)
obs_species_df <- species_abundance_treatment %>%
  filter(species %in% top_species, treatment %in% c(1,3,6)) %>%
  transmute(
    species,
    treatment   = as.integer(as.character(treatment)),
    mean_abund  = mean_abundance,
    lower_obs   = mean_abundance - 1.96 * se,
    upper_obs   = mean_abundance + 1.96 * se
  )

# 6.4 NP bootstrap for k = 3,6 (per species), using single helper ---------------
set.seed(params$seed)


boot_df <- tidyr::expand_grid(
  species   = top_species,
  treatment = c(3L, 6L)
) %>%
  mutate(
    qs = purrr::map2(
      species, treatment,
      ~ np_sum_ci(
        x = species_reps_t1$count[species_reps_t1$species == .x],
        k = .y, B = params$bootstrap_B, seed = params$seed
      )
    ),
    lower_np  = purrr::map_dbl(qs, ~ as.numeric(.x[[1]])),  # 2.5%
    median_np = purrr::map_dbl(qs, ~ as.numeric(.x[[2]])),  # 50%
    upper_np  = purrr::map_dbl(qs, ~ as.numeric(.x[[3]]))   # 97.5%
  ) %>%
  dplyr::select(-qs)

# 6.5 Plot: observed vs. NP bootstrap by species -------------------------------

# keep label_parsed happy
escape_for_parsed <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("'", "\\\\'", x, fixed = TRUE)
  x
}

# --- A) Build species → group map for the plotted species
species_in_plot <- union(obs_species_df$species, boot_df$species)

tax_map <- tibble(species = species_in_plot) %>%
  left_join(MRBcafi_df %>% distinct(species, class), by = "species") %>%
  mutate(
    group_order = case_when(
      class == "Actinopterygii" ~ "Fish",
      class == "Malacostraca"   ~ "Crustacean",
      class == "Gastropoda"     ~ "Snail",
      TRUE                      ~ NA_character_
    ),
    group_order = factor(group_order, levels = c("Fish","Crustacean","Snail"))
  ) %>%
  filter(!is.na(group_order)) %>%
  mutate(
    species_label = paste0("italic('", escape_for_parsed(species), "')"),
    panel_label   = paste0("atop('", as.character(group_order), "', ", species_label, ")")
  )

panel_levels <- tax_map %>% arrange(group_order, species) %>% pull(panel_label)

# --- B) Attach labels to plotting frames
obs_plot <- obs_species_df %>%
  inner_join(tax_map %>% dplyr::select(species, panel_label), by = "species") %>%
  mutate(panel = factor(panel_label, levels = panel_levels))

boot_plot <- boot_df %>%
  inner_join(tax_map %>% dplyr::select(species, panel_label), by = "species") %>%
  mutate(panel = factor(panel_label, levels = panel_levels),
         treatment = as.numeric(treatment)) %>%
  arrange(panel, treatment)


boot_plot <- boot_plot %>%
  # make sure x is numeric and rows are ordered for polygon assembly
  mutate(treatment = as.numeric(treatment)) %>%
  arrange(panel, treatment) %>%
  # guard against any accidental lower>upper (e.g., if a quantile call flips)
  mutate(
    .low = pmin(lower_np, upper_np, na.rm = TRUE),
    .up  = pmax(lower_np, upper_np, na.rm = TRUE)
  ) %>%
  mutate(lower_np = .low, upper_np = .up) %>%
  dplyr::select(-.low, -.up)

obs_plot <- obs_plot %>%
  mutate(treatment = as.numeric(treatment)) %>%
  arrange(panel, treatment)

# Ensure clean polygons for ribbons
boot_plot <- boot_plot %>%
  dplyr::mutate(treatment = as.numeric(treatment)) %>%
  dplyr::arrange(panel, treatment) %>%
  dplyr::mutate(
    .low = pmin(lower_np, upper_np, na.rm = TRUE),
    .up  = pmax(lower_np, upper_np, na.rm = TRUE)
  ) %>%
  dplyr::mutate(lower_np = .low, upper_np = .up) %>%
  dplyr::select(-c(.low, .up)) %>%
  dplyr::filter(!is.na(lower_np) & !is.na(upper_np))

obs_plot <- obs_plot %>%
  dplyr::mutate(treatment = as.numeric(treatment)) %>%
  dplyr::arrange(panel, treatment)

# --- C) Get raw data for each species/treatment combination
raw_focal_data <- species_reps %>%
  inner_join(tax_map %>% dplyr::select(species, panel_label), by = "species") %>%
  mutate(panel = factor(panel_label, levels = panel_levels),
         treatment = as.numeric(as.character(treatment))) %>%
  filter(!is.na(panel))

# --- D) Expected dashed line = k * mean at treatment 1
baseline_t1 <- obs_plot %>% filter(treatment == 1) %>%
  transmute(species, panel, baseline_mean = mean_abund)

exp_line_df <- baseline_t1 %>%
  tidyr::expand_grid(treatment = c(3L, 6L)) %>%
  mutate(expected = baseline_mean * treatment)

# --- D) Per-panel y ranges (include 0; small padding)
# Include observed CIs, bootstrap CIs, AND expected line values to ensure full visibility
rng_df <- dplyr::bind_rows(
  dplyr::transmute(obs_plot, panel, ymin = pmin(lower_obs, 0), ymax = upper_obs),
  dplyr::transmute(boot_plot, panel, ymin = pmin(lower_np,  0), ymax = upper_np),
  dplyr::transmute(exp_line_df, panel, ymin = 0, ymax = expected)  # Include expected line
) %>%
  group_by(panel) %>%
  summarise(ymin = min(ymin, na.rm = TRUE),
            ymax = max(ymax, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(pad = 0.05 * pmax(1e-9, ymax - ymin),  # Increased padding from 0.02 to 0.05
         ymin = ymin - pad,
         ymax = ymax + pad)

my_breaks <- function(lims) unique(sort(c(pretty(lims, n = 4), 0)))

theme_common <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.title.y       = element_text(size = 18, margin = margin(r = 8)),
    axis.title.x       = element_text(size = 18, margin = margin(t = 8)),
    axis.text          = element_text(size = 14, color = "black"),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.tag           = element_text(face = "bold", size = 18),
    plot.tag.position  = c(0.97, 0.98)
  )


# Helper: draw the ribbon bounds too (great for debugging visibility)
debug_bounds <- FALSE  # set TRUE to draw lower/upper as thin lines

#THIS IS THE A PUBLICATION FIGURE


p_focal_np <- ggplot() +
  geom_blank(data = rng_df, aes(x = 1, y = ymin)) +
  geom_blank(data = rng_df, aes(x = 6, y = ymax)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  # 1) RIBBON: explicit x, explicit per-panel group, gray fill
  geom_ribbon(
    data = boot_plot,
    aes(x = treatment, ymin = lower_np, ymax = upper_np, group = panel),
    fill  = params$ribbon_color,
    color = NA,   # no border
    alpha = params$ribbon_alpha,
    na.rm = TRUE
  ) +
  # (optional) draw the ribbon edges for sanity-checking
  { if (debug_bounds) geom_line(data = boot_plot, aes(x = treatment, y = lower_np, group = panel)) } +
  { if (debug_bounds) geom_line(data = boot_plot, aes(x = treatment, y = upper_np, group = panel)) } +
  # 2) Expected dashed line (k * mean at t=1)
  geom_line(
    data = exp_line_df,
    aes(x = treatment, y = expected, group = panel),
    color = params$col_exp, linetype = "dashed", linewidth = 0.9
  ) +
  # 3) Observed mean (black line and points, no error bars)
  geom_line(
    data = obs_plot,
    aes(x = treatment, y = mean_abund, group = panel),
    color = params$col_obs, linewidth = 0.9
  ) +
  geom_point(
    data = obs_plot,
    aes(x = treatment, y = mean_abund),
    color = "white", fill = params$col_obs,
    shape = 21, size = 2.7, stroke = 0.5
  ) +
  # 4) Raw data points with treatment colors (plotted last so they're on top)
  geom_jitter(
    data = raw_focal_data,
    aes(x = treatment, y = count, fill = factor(treatment)),
    width = 0.08, height = 0, alpha = 0.4, size = 1.8,
    shape = 21, color = "white", stroke = 0.4
  ) +
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  facet_wrap(~ panel, ncol = 4, scales = "free_y", labeller = ggplot2::label_parsed) +
  scale_x_continuous(breaks = c(1,3,6)) +
  scale_y_continuous(breaks = my_breaks, expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Number of corals", y = "Total abundance") +
  theme_common

outfile <- file.path(paths$out_dir_fig, "focal_order_species_np.png")
ggsave(outfile, p_focal_np, width = 11, height = 8,
       dpi = params$fig_dpi, bg = params$fig_bg)

p_focal_np_labeled <- p_focal_np + patchwork::plot_annotation(tag_levels = "A")

outfile_labeled <- file.path(paths$out_dir_fig, "focal_order_species_np_labeled.pdf")
ggsave(outfile_labeled, p_focal_np_labeled, width = 11, height = 8,
       dpi = params$fig_dpi, bg = params$fig_bg)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(dirname(paths$out_dir_fig), "publication-figures", "focal_order_species_np_labeled.pdf"),
  p_focal_np_labeled, width = 11, height = 8, dpi = params$fig_dpi, bg = params$fig_bg
)
ggsave(
  file.path(dirname(paths$out_dir_fig), "publication-figures", "focal_order_species_np_labeled.png"),
  p_focal_np_labeled, width = 11, height = 8, dpi = params$fig_dpi, bg = params$fig_bg
)

# 6.6 Top species by order: expected vs. observed table (saved to CSV) ---------
top_species_meta <- MRBcafi_df %>%
  dplyr::select(species, order, family) %>% distinct()

top_species_df <- species_abundance_treatment %>%
  filter(species %in% top_species) %>%
  left_join(top_species_meta, by = "species") %>%
  group_by(species) %>%
  mutate(
    baseline_mean = mean_abundance[treatment == 1][1],
    baseline_se   = se[treatment == 1][1],
    baseline_n    = n[treatment == 1][1]
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    t_crit = qt(0.975, df = baseline_n - 1),
    expected_abundance = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_mean * 3,
      treatment == 6 ~ baseline_mean * 6,
      TRUE ~ NA_real_
    ),
    expected_lower = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_mean * 3 - t_crit * sqrt(3) * baseline_se,
      treatment == 6 ~ baseline_mean * 6 - t_crit * sqrt(6) * baseline_se,
      TRUE ~ NA_real_
    ),
    expected_upper = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_mean * 3 + t_crit * sqrt(3) * baseline_se,
      treatment == 6 ~ baseline_mean * 6 + t_crit * sqrt(6) * baseline_se,
      TRUE ~ NA_real_
    ),
    lower_ci = mean_abundance - 1.96 * se,
    upper_ci = mean_abundance + 1.96 * se
  ) %>%
  ungroup() %>%
  mutate(
    treatment_num = as.numeric(as.character(treatment)),
    order_species = paste0(order, ": ", species)
  )

readr::write_csv(top_species_df, file.path(paths$out_dir_data, "top_species_treatment_data.csv"))
cli::cli_alert_success("Saved focal-order top-species table to {file.path(paths$out_dir_data, 'top_species_treatment_data.csv')}.")

# ==== 7. Richness & Rarefied Richness =========================================

cli::cli_h2("Richness metrics")

# 7.1 Compute reef-level richness and rarefied richness ------------------------
row_totals <- rowSums(community_matrix)
species_richness  <- rowSums(community_matrix > 0)                       # presence/absence
min_sample_size   <- max(1L, min(row_totals))                            # guard sample >= 1
rarefied_richness <- suppressWarnings(vegan::rarefy(community_matrix, sample = min_sample_size))

community_summary <- tibble(
  reef              = rownames(community_matrix),
  richness          = species_richness,
  rarefied_richness = rarefied_richness
) %>%
  left_join(metadata, by = "reef")

# 7.2 Summaries & CIs by treatment --------------------------------------------
richness_summary <- community_summary %>%
  group_by(treatment) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness   = sd(richness,   na.rm = TRUE),
    n             = dplyr::n(),
    .groups       = "drop"
  ) %>%
  mutate(
    se            = sd_richness / sqrt(n),
    lower_ci      = mean_richness - 1.96 * se,
    upper_ci      = mean_richness + 1.96 * se,
    treatment_num = as.numeric(as.character(treatment))
  )

baseline_rich_mean <- richness_summary %>% filter(treatment == 1) %>% pull(mean_richness) %>% .[1]
baseline_rich_se   <- richness_summary %>% filter(treatment == 1) %>% pull(se)            %>% .[1]

richness_summary <- richness_summary %>%
  mutate(
    expected_richness = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_rich_mean * 3,
      treatment == 6 ~ baseline_rich_mean * 6,
      TRUE ~ NA_real_
    ),
    expected_lower = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_rich_mean * 3 - 1.96 * (baseline_rich_se * 3),
      treatment == 6 ~ baseline_rich_mean * 6 - 1.96 * (baseline_rich_se * 6),
      TRUE ~ NA_real_
    ),
    expected_upper = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_rich_mean * 3 + 1.96 * (baseline_rich_se * 3),
      treatment == 6 ~ baseline_rich_mean * 6 + 1.96 * (baseline_rich_se * 6),
      TRUE ~ NA_real_
    )
  )

rarefied_summary <- community_summary %>%
  group_by(treatment) %>%
  summarise(
    mean_rarefied = mean(rarefied_richness, na.rm = TRUE),
    sd_rarefied   = sd(rarefied_richness,   na.rm = TRUE),
    n             = dplyr::n(),
    .groups       = "drop"
  ) %>%
  mutate(
    se            = sd_rarefied / sqrt(n),
    lower_ci      = mean_rarefied - 1.96 * se,
    upper_ci      = mean_rarefied + 1.96 * se,
    treatment_num = as.numeric(as.character(treatment))
  )

baseline_rare_mean <- rarefied_summary %>% filter(treatment == 1) %>% pull(mean_rarefied) %>% .[1]
baseline_rare_se   <- rarefied_summary %>% filter(treatment == 1) %>% pull(se)            %>% .[1]

rarefied_summary <- rarefied_summary %>%
  mutate(
    expected_rarefied = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_rare_mean * 3,
      treatment == 6 ~ baseline_rare_mean * 6,
      TRUE ~ NA_real_
    ),
    expected_lower = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_rare_mean * 3 - 1.96 * (baseline_rare_se * 3),
      treatment == 6 ~ baseline_rare_mean * 6 - 1.96 * (baseline_rare_se * 6),
      TRUE ~ NA_real_
    ),
    expected_upper = case_when(
      treatment == 1 ~ NA_real_,
      treatment == 3 ~ baseline_rare_mean * 3 + 1.96 * (baseline_rare_se * 3),
      treatment == 6 ~ baseline_rare_mean * 6 + 1.96 * (baseline_rare_se * 6),
      TRUE ~ NA_real_
    )
  )

# 7.3 Plots with raw data points -----------------------------------------------
# Prepare raw richness data (will be created in section 7.1)
raw_richness_7 <- community_summary %>%
  mutate(treatment_num = as.integer(as.character(treatment))) %>%
  filter(treatment_num %in% c(1, 3, 6))

p_richness <- ggplot(richness_summary, aes(x = treatment_num, y = mean_richness)) +
  # Add raw data points
  geom_jitter(data = raw_richness_7,
              aes(x = treatment_num, y = richness, fill = factor(treatment_num)),
              width = 0.15, height = 0, alpha = 0.4, size = 2.5,
              shape = 21, color = "white", stroke = 0.5) +
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  # Black line and points for observed mean
  geom_line(color = params$col_obs, linewidth = 1) +
  geom_point(shape = 21, fill = params$col_obs, color = "white", size = 3) +
  scale_x_continuous(breaks = richness_summary$treatment_num,
                     labels = richness_summary$treatment) +
  labs(
    title = "Observed Species Richness",
    x     = "Treatment (Number of Corals)",
    y     = "Species Richness"
  ) +
  theme_publication() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_rarefied <- ggplot(rarefied_summary, aes(x = treatment_num, y = mean_rarefied)) +
  # Add raw data points
  geom_jitter(data = raw_richness_7,
              aes(x = treatment_num, y = rarefied_richness, fill = factor(treatment_num)),
              width = 0.15, height = 0, alpha = 0.4, size = 2.5,
              shape = 21, color = "white", stroke = 0.5) +
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  # Black line and points for observed mean
  geom_line(color = params$col_obs, linewidth = 1) +
  geom_point(shape = 21, fill = params$col_obs, color = "white", size = 3) +
  scale_x_continuous(breaks = rarefied_summary$treatment_num,
                     labels = rarefied_summary$treatment) +
  labs(
    title = "Rarefied Species Richness",
    x     = "Treatment (Number of Corals)",
    y     = "Rarefied Richness"
  ) +
  theme_publication() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# 7.4–7.5 Modeling (unchanged except earlier logging) --------------------------
cli::cli_h2("Modeling raw species richness")
rich_model_df <- community_summary %>%
  mutate(treatment = factor(treatment, levels = c("1","3","6"))) %>%
  filter(!is.na(richness))

m_rich_pois <- glm(richness ~ treatment, family = poisson, data = rich_model_df)
disp_rich   <- performance::check_overdispersion(m_rich_pois)
use_nb_rich <- disp_rich$dispersion_ratio > 1.2 && disp_rich$p < 0.05

if (use_nb_rich) {
  message("Overdispersion detected in richness model; using negative binomial.")
  m_rich   <- MASS::glm.nb(richness ~ treatment, data = rich_model_df)
  fam_rich <- "negative binomial"
} else {
  message("No overdispersion in richness model; using Poisson.")
  m_rich   <- m_rich_pois
  fam_rich <- "poisson"
}

summary(m_rich)
car::Anova(m_rich, type = 3)
emm_rich <- emmeans::emmeans(m_rich, ~ treatment)
pairs_rich <- emmeans::contrast(emm_rich, method = "pairwise")
emm_rich_tbl  <- summary(emm_rich, type = "response") %>% tibble::as_tibble()
pairs_rich_tbl <- summary(pairs_rich, type = "response", infer = c(TRUE, TRUE)) %>% tibble::as_tibble()

cli::cli_h2("Modeling rarefied richness")
rarefied_model_df <- community_summary %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("1","3","6"))) %>%
  dplyr::filter(!is.na(rarefied_richness))

m_rare <- lm(rarefied_richness ~ treatment, data = rarefied_model_df)
summary(m_rare)
omni_rare_tbl <- car::Anova(m_rare, type = 3) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term")

emm_rare      <- emmeans::emmeans(m_rare, ~ treatment)
emm_rare_tbl  <- summary(emm_rare) %>% tibble::as_tibble()
pairs_rare    <- emmeans::contrast(emm_rare, method = "pairwise")
pairs_rare_tbl <- summary(pairs_rare, infer = c(TRUE, TRUE)) %>% tibble::as_tibble()

# ==== 8. Plotting Functions (publication-ready) ================================

cli::cli_h2("Generating figures")

# Build community_abundance once from §5 objects
community_abundance <- obs_comm_df %>%
  dplyr::select(treatment, total_abundance, lower_obs, upper_obs) %>%
  left_join(
    np_ci_comm %>% mutate(treatment = as.integer(treatment)),
    by = "treatment"
  ) %>%
  transmute(
    treatment_num   = treatment,
    total_abundance,
    lower_ci        = lower_obs,
    upper_ci        = upper_obs,
    expected_total  = median_np,
    expected_lower  = lower_np,
    expected_upper  = upper_np
  )

x_breaks <- sort(unique(c(
  community_abundance$treatment_num,
  unique(as.integer(as.character(richness_summary$treatment))),
  unique(as.integer(as.character(rarefied_summary$treatment)))
)))

x_axis_common <- list(scale_x_continuous(breaks = x_breaks, labels = as.character(x_breaks)))

theme_common <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", linewidth = 0.8),
    axis.title.y     = element_text(size = 18, face = "bold"),
    axis.title.x     = element_text(size = 18, face = "bold"),
    axis.text        = element_text(size = 16, color = "black")
  )

pt_sz <- 5.5

# 8.1 Community abundance panel with raw data points
p_comm <- ggplot(community_abundance, aes(x = treatment_num)) +
  # Add raw data points
  geom_jitter(data = raw_totals,
              aes(x = treatment, y = total_abun, fill = factor(treatment)),
              width = 0.15, height = 0, alpha = 0.5, size = 3,
              shape = 21, color = "white", stroke = 0.6) +
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  # Black line and points for observed mean
  geom_point(aes(y = total_abundance), size = pt_sz, shape = 21, fill = params$col_obs, color = "white", stroke = 0.8) +
  geom_line(aes(y = total_abundance, group = 1), color = params$col_obs, linewidth = 1.3) +
  # Expected values - gray ribbon, red dashed line
  geom_ribbon(aes(ymin = expected_lower, ymax = expected_upper), fill = params$ribbon_color, alpha = params$ribbon_alpha) +
  geom_line(aes(y = expected_total, group = 1), color = params$col_exp, linetype = "dashed", linewidth = 1.3) +
  labs(y = "Total abundance") +
  x_axis_common + theme_common +
  coord_cartesian(ylim = c(0, NA)) +
  theme(axis.title.x = element_blank())

# 8.2 Species richness panel with raw data points
# Get raw richness data
raw_richness <- community_summary %>%
  mutate(treatment_num = as.integer(as.character(treatment))) %>%
  filter(treatment_num %in% c(1, 3, 6))

p_rich <- ggplot(richness_summary, aes(x = treatment_num, y = mean_richness)) +
  # Add raw data points
  geom_jitter(data = raw_richness,
              aes(x = treatment_num, y = richness, fill = factor(treatment_num)),
              width = 0.15, height = 0, alpha = 0.5, size = 3,
              shape = 21, color = "white", stroke = 0.6) +
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  # Black line and points for observed mean
  geom_line(color = params$col_obs, linewidth = 1.3) +
  geom_point(shape = 21, fill = params$col_obs, color = "white", size = pt_sz, stroke = 0.8) +
  labs(y = "Species richness") +
  x_axis_common + theme_common +
  coord_cartesian(ylim = c(0, NA)) +
  theme(axis.title.x = element_blank())

# 8.3 Rarefied richness panel with raw data points
p_rare <- ggplot(rarefied_summary, aes(x = treatment_num, y = mean_rarefied)) +
  # Add raw data points
  geom_jitter(data = raw_richness,
              aes(x = treatment_num, y = rarefied_richness, fill = factor(treatment_num)),
              width = 0.15, height = 0, alpha = 0.5, size = 3,
              shape = 21, color = "white", stroke = 0.6) +
  scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
  # Black line and points for observed mean
  geom_line(color = params$col_obs, linewidth = 1.3) +
  geom_point(shape = 21, fill = params$col_obs, color = "white", size = pt_sz, stroke = 0.8) +
  labs(y = "Rarefied richness", x = "Number of corals") +
  x_axis_common + theme_common +
  coord_cartesian(ylim = c(0, NA))

# 8.4 Create publication-quality legend with all treatment colors
library(cowplot)

# Simpler approach: create legend with raw points plot
legend_data <- tibble(
  Element = factor(
    rep(c("Observed mean", "Expected mean", "95% CI (expected)", "1 coral", "3 corals", "6 corals"), each = 2),
    levels = c("Observed mean", "Expected mean", "95% CI (expected)", "1 coral", "3 corals", "6 corals")
  ),
  x = rep(1:2, 6),
  y = rep(1:6, each = 2)
)

p_legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Element, fill = Element,
                                          linetype = Element, shape = Element)) +
  geom_point() +
  scale_color_manual(
    values = c("Observed mean" = "black",
               "Expected mean" = params$col_exp,
               "95% CI (expected)" = "gray60",
               "1 coral" = "white",
               "3 corals" = "white",
               "6 corals" = "white")
  ) +
  scale_fill_manual(
    values = c("Observed mean" = "black",
               "Expected mean" = "white",
               "95% CI (expected)" = params$ribbon_color,
               "1 coral" = TREATMENT_COLORS["1"],
               "3 corals" = TREATMENT_COLORS["3"],
               "6 corals" = TREATMENT_COLORS["6"])
  ) +
  scale_linetype_manual(
    values = c("Observed mean" = "solid",
               "Expected mean" = "dashed",
               "95% CI (expected)" = "solid",
               "1 coral" = "blank",
               "3 corals" = "blank",
               "6 corals" = "blank")
  ) +
  scale_shape_manual(
    values = c("Observed mean" = 21,
               "Expected mean" = NA,
               "95% CI (expected)" = 22,
               "1 coral" = 21,
               "3 corals" = 21,
               "6 corals" = 21)
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = c(21, NA, 22, 21, 21, 21),
        size = c(2.8, NA, 3.5, 2.3, 2.3, 2.3),
        linewidth = c(0.9, 0.9, NA, NA, NA, NA),
        linetype = c("solid", "dashed", "blank", "blank", "blank", "blank"),
        stroke = c(0.5, NA, 0.7, 0.5, 0.5, 0.5),
        alpha = c(1, 1, 1, 0.5, 0.5, 0.5)
      ),
      label.theme = element_text(size = 10),
      keywidth = unit(0.7, "cm"),
      keyheight = unit(0.4, "cm")
    ),
    fill = "none",
    linetype = "none",
    shape = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(3, 0, 3, 0),
    legend.spacing.x = unit(0.25, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )

# Extract legend
legend_grob <- cowplot::get_legend(p_legend_plot)

# 8.5 Stack panels WITHOUT legend (publication version)
threepanel <- cowplot::plot_grid(
  p_comm, p_rich, p_rare,
  ncol = 1, align = "v", rel_heights = c(1, 1, 1),
  labels = c("a", "b", "c"),
  label_size = 18, label_fontface = "bold",
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1
)

print(threepanel)

#THIS IS THE A PUBLICATION FIGURE


ggsave(
  file.path(paths$out_dir_fig, "three_panel_nonparametric.pdf"),
  threepanel,
  width  = params$fig_threepanel[1],
  height = params$fig_threepanel[2],
  dpi    = params$fig_dpi,
  bg     = params$fig_bg
)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(dirname(paths$out_dir_fig), "publication-figures", "three_panel_nonparametric.pdf"),
  threepanel,
  width  = params$fig_threepanel[1],
  height = params$fig_threepanel[2],
  dpi    = params$fig_dpi,
  bg     = params$fig_bg
)
ggsave(
  file.path(dirname(paths$out_dir_fig), "publication-figures", "three_panel_nonparametric.png"),
  threepanel,
  width  = params$fig_threepanel[1],
  height = params$fig_threepanel[2],
  dpi    = params$fig_dpi,
  bg     = params$fig_bg
)

cli::cli_alert_success("Saved nonparametric three-panel figure (tagged a/b/c, shared x-axis).")

# =============================================================================
# 9. COMMUNITY SCALING: OBSERVED vs. EXPECTED (Nonparametric Bootstrap)
# =============================================================================

cli::cli_h2("Section 9: Observed vs. Expected Abundance by Treatment")

# 9.0: Long-format community data
comm_long <- species_matrix %>%
  tidyr::pivot_longer(-reef, names_to = "species", values_to = "abundance") %>%
  left_join(metadata, by = "reef") %>%
  filter(!is.na(abundance)) %>%
  mutate(treatment = as.integer(as.character(treatment)))

# 9.1: Identify Top 30 Most Abundant Species
top30_species <- comm_long %>%
  group_by(species) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(species)

# 9.2: Observed Means and 95% CIs by Treatment (3,6)
obs_summary_all <- comm_long %>%
  filter(species %in% top30_species, treatment %in% c(3, 6)) %>%
  group_by(species, treatment) %>%
  summarise(
    observed_mean = mean(abundance),
    observed_sd   = sd(abundance),
    n             = n(),
    lower_obs     = observed_mean - 1.96 * (observed_sd / sqrt(n)),
    upper_obs     = observed_mean + 1.96 * (observed_sd / sqrt(n)),
    .groups = "drop"
  )

# 9.3: Expected Means from Treatment 1 (Nonparametric Bootstrap, helper)
set.seed(params$seed)
boot_results_all <- comm_long %>%
  filter(species %in% top30_species, treatment == 1) %>%
  group_by(species) %>%
  summarise(
    q3 = list(np_sum_ci(abundance, 3, B = params$bootstrap_B, seed = params$seed)),
    q6 = list(np_sum_ci(abundance, 6, B = params$bootstrap_B, seed = params$seed)),
    .groups = "drop"
  ) %>%
  transmute(
    species,
    expected_median_3 = purrr::map_dbl(q3, ~ .x[2]),
    expected_lower_3  = purrr::map_dbl(q3, ~ .x[1]),
    expected_upper_3  = purrr::map_dbl(q3, ~ .x[3]),
    expected_median_6 = purrr::map_dbl(q6, ~ .x[2]),
    expected_lower_6  = purrr::map_dbl(q6, ~ .x[1]),
    expected_upper_6  = purrr::map_dbl(q6, ~ .x[3])
  )

# 9.4: Merge Observed and Expected
combined_results <- obs_summary_all %>%
  left_join(boot_results_all, by = "species") %>%
  mutate(
    expected_median = dplyr::if_else(treatment == 3, expected_median_3, expected_median_6),
    expected_lower  = dplyr::if_else(treatment == 3, expected_lower_3,  expected_lower_6),
    expected_upper  = dplyr::if_else(treatment == 3, expected_upper_3,  expected_upper_6),
    direction = dplyr::case_when(
      upper_obs < expected_lower ~ "Below Expected",
      lower_obs > expected_upper ~ "Above Expected",
      TRUE                       ~ "Proportional"
    ),
    diff = observed_mean - expected_median,
    sig  = direction != "Proportional"
  )

# 9.5: Tidy summary table (gt)
summary_table_all <- combined_results %>%
  arrange(desc(sig), desc(abs(diff))) %>%
  mutate(
    observed_mean   = round(observed_mean, 2),
    observed_sd     = round(observed_sd, 2),
    lower_obs       = round(lower_obs, 2),
    upper_obs       = round(upper_obs, 2),
    expected_median = round(expected_median, 2),
    expected_lower  = round(expected_lower, 2),
    expected_upper  = round(expected_upper, 2),
    diff            = round(diff, 2)
  ) %>%
  dplyr::select(
    species, treatment,
    observed_mean, observed_sd, lower_obs, upper_obs,
    expected_median, expected_lower, expected_upper,
    diff, direction, sig
  ) %>%
  gt::gt() %>%
  gt::cols_label(
    species         = "Species",
    treatment       = "Treatment",
    observed_mean   = "Observed Mean",
    observed_sd     = "Observed SD",
    lower_obs       = "Obs CI Lower",
    upper_obs       = "Obs CI Upper",
    expected_median = "Expected Median",
    expected_lower  = "Exp CI Lower",
    expected_upper  = "Exp CI Upper",
    diff            = "Δ (Obs - Exp)",
    direction       = "Direction",
    sig             = "Significant"
  ) %>%
  gt::data_color(
    columns = diff,
    fn = scales::col_bin(
      palette = c("steelblue", "grey90", "tomato"),
      domain = NULL,
      bins = c(-Inf, -1, 1, Inf)
    )
  ) %>%
  gt::sub_missing(columns = everything(), missing_text = "–")

# 9.6: Save CSV and PDF table (with PNG fallback)
readr::write_csv(combined_results, file.path(paths$out_dir_data, "top30_observed_vs_expected_all.csv"))

pdf_path <- file.path(paths$out_dir_fig, "top30_observed_vs_expected_all.pdf")
png_path <- file.path(paths$out_dir_fig, "top30_observed_vs_expected_all.png")
ok <- try(gt::gtsave(summary_table_all, filename = pdf_path), silent = TRUE)
if (inherits(ok, "try-error")) {
  cli::cli_warn("PDF save failed for gt table; falling back to PNG.")
  gt::gtsave(summary_table_all, filename = png_path)
}
summary_table_all

# 9.7: Visualization – Observed vs. Expected Abundance by Species (1, 3, 6)
cli::cli_h3("Section 9.7: Per-species plots with observed & bootstrap expected")

obs_summary_all_extended <- comm_long %>%
  filter(species %in% top30_species, treatment %in% c(1, 3, 6)) %>%
  group_by(species, treatment) %>%
  summarise(
    observed_mean = mean(abundance),
    observed_sd   = sd(abundance),
    n             = n(),
    .groups = "drop"
  )

plot_data <- combined_results %>%
  bind_rows(
    obs_summary_all_extended %>%
      filter(treatment == 1) %>%
      mutate(
        expected_median = NA_real_,
        expected_lower  = NA_real_,
        expected_upper  = NA_real_,
        direction       = "Observed only",
        diff            = NA_real_,
        sig             = FALSE
      )
  ) %>%
  arrange(species, treatment)

species_plot_list <- plot_data %>%
  mutate(species = forcats::fct_reorder(species, dplyr::desc(abs(diff)), .na_rm = TRUE)) %>%
  group_split(species) %>%
  purrr::map(~ {
    df <- .
    sp_name <- unique(df$species)

    # Get raw data for this species
    raw_species_data <- comm_long %>%
      filter(species == sp_name, treatment %in% c(1, 3, 6))

    expected_df <- df %>%
      filter(treatment %in% c(3, 6)) %>%
      dplyr::select(treatment, expected_median, expected_lower, expected_upper)

    ymax <- max(
      df$observed_mean + 1.96 * df$observed_sd / sqrt(df$n),
      expected_df$expected_upper,
      raw_species_data$abundance,
      na.rm = TRUE
    ) * 1.1

    ggplot(df, aes(x = treatment)) +
      # Expected ribbon with gray halo
      geom_ribbon(
        data = expected_df,
        aes(x = treatment, ymin = expected_lower, ymax = expected_upper),
        inherit.aes = FALSE,
        fill = params$ribbon_color, alpha = params$ribbon_alpha
      ) +
      geom_line(
        data = expected_df,
        aes(x = treatment, y = expected_median, group = 1),
        inherit.aes = FALSE,
        color = params$col_exp, linetype = "dashed", linewidth = 0.9
      ) +
      # Raw data points with treatment colors
      geom_jitter(
        data = raw_species_data,
        aes(x = treatment, y = abundance, fill = factor(treatment)),
        width = 0.12, height = 0, alpha = 0.4, size = 2,
        shape = 21, color = "white", stroke = 0.4
      ) +
      scale_fill_manual(values = TREATMENT_COLORS, guide = "none") +
      # Black line and points for observed mean
      geom_line(aes(y = observed_mean, group = 1),
                color = params$col_obs, linewidth = 0.9) +
      geom_point(aes(y = observed_mean),
                 shape = 21, fill = params$col_obs, color = "white", size = 3) +
      scale_x_continuous(breaks = c(1, 3, 6)) +
      labs(
        title = paste0("Species: ", sp_name),
        x     = "Treatment (Number of Corals)",
        y     = "Abundance"
      ) +
      theme_publication() +
      theme(plot.title = element_text(face = "bold")) +
      coord_cartesian(ylim = c(0, ymax))
  })

fs::dir_create(paths$out_dir_fig)
pdf(file.path(paths$out_dir_fig, "species_observed_vs_expected.pdf"), width = 7, height = 5)
purrr::walk(species_plot_list, print)
dev.off()

# ==== 10. Save Outputs ---------------------------------------------------------

cli::cli_h2("Saving figures & objects")

save_plot <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = params$fig_dpi,
    bg       = params$fig_bg
  )
}

if (exists("p_obs_np2")) {
  save_plot(
    p_obs_np2,
    file.path(paths$out_dir_fig, "community_total_abundance.png"),
    params$fig_comm_size[1], params$fig_comm_size[2]
  )
}

if (exists("p_richness")) {
  save_plot(
    p_richness,
    file.path(paths$out_dir_fig, "species_richness.png"),
    params$fig_comm_size[1], params$fig_comm_size[2]
  )
}

if (exists("p_rarefied")) {
  save_plot(
    p_rarefied,
    file.path(paths$out_dir_fig, "rarefied_richness.png"),
    params$fig_comm_size[1], params$fig_comm_size[2]
  )
}

# -- Data artifacts saved above:
# - Section 6: top_species_treatment_data.csv
# - Section 8: three_panel_nonparametric.pdf
# - Section 9: top30_observed_vs_expected_all.csv, top30_observed_vs_expected_all.pdf/PNG,
#              species_observed_vs_expected.pdf

# -- Save key R objects for reuse ---------------------------------------------
objs_to_save <- list()
if (exists("community_abundance")) objs_to_save$community_abundance <- community_abundance
if (exists("richness_summary"))    objs_to_save$richness_summary    <- richness_summary
if (exists("rarefied_summary"))    objs_to_save$rarefied_summary    <- rarefied_summary
if (exists("top_species_df"))      objs_to_save$top_species_df      <- top_species_df
if (exists("combined_results"))    objs_to_save$combined_results    <- combined_results
# Optional: model outputs from Section 4
if (exists("omni_tbl"))            objs_to_save$omni_tbl            <- omni_tbl
if (exists("coef_tbl"))            objs_to_save$coef_tbl            <- coef_tbl
if (exists("emm_tbl"))             objs_to_save$emm_tbl             <- emm_tbl
if (exists("pairs_tbl"))           objs_to_save$pairs_tbl           <- pairs_tbl

saveRDS(
  objs_to_save,
  file = file.path(paths$out_dir_obj, "mrb_comm_summaries.rds")
)

# ==== 10B. Session Info --------------------------------------------------------
sink(file.path(paths$out_dir_obj, "sessionInfo_05_MRB_comm.txt"))
print(sessionInfo())
sink()

cli::cli_alert_success("Saved figures, key objects (RDS), and session info.")