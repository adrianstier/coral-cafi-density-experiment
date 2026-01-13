# ==============================================================================
# File: 8.coral-caffi.R
# Purpose: Build coral × species community matrices; compute community PCAs
#          (RAW & Hellinger) and a coral *performance* PCA; align PC1 orientations;
#          relate Community PC1 to Performance PC1 with LMMs; fit per-species
#          LMMs for the top-20 RAW PC1 loadings; and export publication-ready
#          figures and tables.
#
# Inputs:
#   data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv   (CAFI counts)
#   data/processed/coral_growth.csv                                      (growth, reef, treatment)
#   output/MRB/figures/coral/physio/physio_metrics_plus_growth_filtered.csv
#       (performance metrics file: physiology + growth merged by coral_id)
#
# Outputs:
#   output/MRB/figures/combined/
#     - scree_COMM_{RAW,HELLINGER}.{png,pdf}                 # community PCA scree
#     - scree_COND_z.{png,pdf}                               # performance PCA scree
#     - PC1_loadings_COMM_{RAW,HELLINGER}{_aligned}.{png,pdf}
#     - PCA_LOADINGS_HELLINGER_3panel_clean.{png,pdf}
#     - PCA_LOADINGS_RAW_3panel_clean.{png,pdf}
#     - PCA_LOADINGS_RAW_2panel_clean.{png,pdf}
#     - LMM_COMM_vs_COND_comparison.{png,pdf}                # community PC1 → performance PC1
#     - LMM_species_top20_coefplot.{png,pdf}                 # per-species β ± 95% CI
#     - species_top20_LMM_panels.{png,pdf}                   # top-20 species panels (sqrt x)
#     - species_faceted_LMM_lines_rawX_sqrtAxis.{png,pdf}    # raw x w/ sqrt axis ticks
#
#   output/MRB/tables/
#     - PC1_loadings_COMM_{RAW,HELLINGER}.csv
#     - PC1_loadings_COMM_HELLINGER_aligned.csv
#     - LMM_{RAW,HELLINGER}_{fixed,anova}.csv
#     - LMM_COMM_vs_COND_summary.csv                         # fixed effects for RAW/HELL
#     - LMM_species_top20_RAWsqrt_vs_condPC1.csv             # per-species LMM results
#
# Depends: R (>= 4.3), here, glue, cli,
#          dplyr, tidyr, tibble, stringr, purrr, readr, forcats,
#          ggplot2, patchwork, scales,
#          vegan,
#          lme4, lmerTest, broom.mixed, car, gt
#
# Run after: required data files exist and the performance metrics file
#            (physiology + growth) has been generated.
#
# Author: Adrian C. Stier
# Date: 2025-08-07
# Last updated: 2025-11-14
#
# Repro notes:
#   - set.seed(1234) for any stochastic steps.
#   - All outputs are written to `output/MRB/figures/combined/` and `.../tables/`.
#   - Species filtering thresholds controlled by PREV_MIN_COMM and ABUND_MIN_COMM.
#
# Growth metric notes (UPDATED 2025-11-14):
#   - growth_vol_b (from data/processed/coral_growth.csv) uses UNIFIED allometric model
#   - Pooled b estimate across all treatments ensures comparable growth metric
#   - All CAFI-coral growth relationships use this unified growth metric
#   - See script 6 section 8B for details on unified model approach
# ==============================================================================

# =============================================================================
# 0) SETUP / HOUSEKEEPING
# =============================================================================

# Source centralized libraries, utilities, and figure standards
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")  # Provides ALIVE_THRESH, strip_fe(), load_cafi_data(), etc.
source("scripts/MRB/mrb_figure_standards.R")

## 0.1 Reproducibility & global options ----------------------------------------
set.seed(1234)  # used anywhere stochastic (e.g., bootstraps)
options(
  stringsAsFactors          = FALSE,  # never auto-convert to factors
  scipen                    = 999,    # avoid scientific notation in printing
  dplyr.summarise.inform    = FALSE,  # quiet dplyr's summarise() messages
  readr.show_col_types      = FALSE   # quiet readr column-type messages
)

## 0.2 Packages -----------------------------------------------------------------
# All packages loaded via source("scripts/MRB/1.libraries.R") above



## 0.3 Project directories ------------------------------------------------------
DATA_DIR <- here::here("data")
OUT_DIR  <- here::here("output", "MRB")
FIG_DIR  <- file.path(OUT_DIR, "figures", "cafi-coral")
TAB_DIR  <- file.path(OUT_DIR, "tables")

invisible(lapply(list(FIG_DIR, TAB_DIR), dir.create, recursive = TRUE, showWarnings = FALSE))
cli::cli_alert_info("Figure dir: {FIG_DIR}")
cli::cli_alert_info("Table  dir: {TAB_DIR}")

## 0.4 Analysis constants -------------------------------------------------------
# Filtering thresholds for community species
PREV_MIN_COMM  <- 10   # min number of corals a species must appear on
ABUND_MIN_COMM <- 10   # min total counts per species across all corals

# How many top (|loading|) features to show in loading plots
TOP_N_LOAD     <- 20

## 0.5 Utilities ---------------------------------------------------------------

# strip_fe() defined in utils.R - sourced above
# (Removes "FE-" prefix from coral IDs for consistent IDs across files)

# ---- Colors & Theme ----------------------------------------------------------
# Use TREATMENT_COLORS from mrb_figure_standards.R:
#   "1" = "#E69F00" (Orange)
#   "3" = "#56B4E9" (Sky Blue)
#   "6" = "#009E73" (Green)
# Use theme_publication() from mrb_figure_standards.R
# Use save_figure() from mrb_figure_standards.R
#
# Legacy aliases for backward compatibility (save_both, show_and_save provided by utils.R)
cols_trt <- TREATMENT_COLORS
theme_pub <- theme_publication


# === Bring in percent_alive (threshold defined in utils.R) ===
# ALIVE_THRESH is defined in utils.R (0.80) - sourced at top of script

# Option A (prefer): reuse the file you already export in the growth+physio pipeline
# data/MRB Amount/coral_growth_surface_area_change_filtered.csv includes percent_alive?
alive_src <- file.path(DATA_DIR, "MRB Amount", "coral_growth_surface_area_change_filtered.csv")
if (file.exists(alive_src)) {
  alive_meta <- readr::read_csv(alive_src, show_col_types = FALSE) %>%
    dplyr::mutate(coral_id = strip_fe(as.character(coral_id)))
} else {
  # Option B (fallback): read the original manual % alive file directly
  manual_file <- here::here("data","MRB Amount",
                            "1. amount_manual_colony_measurements_dec2019_and_may2021.xlsx")
  alive_meta <- readxl::read_excel(manual_file) %>%
    dplyr::transmute(
      coral_id      = strip_fe(as.character(coral_id)),
      percent_alive = percent_alive_may21 / 100
    )
}

stopifnot("coral_id" %in% names(alive_meta), "percent_alive" %in% names(alive_meta))

keep_ids_alive <- alive_meta %>%
  dplyr::filter(!is.na(percent_alive), percent_alive >= ALIVE_THRESH) %>%
  dplyr::pull(coral_id)

cli::cli_alert_info("{length(keep_ids_alive)} colonies retained by ≥{ALIVE_THRESH*100}% alive filter.")

# =============================================================================
# 1) LOAD, PROCESS, AND BUILD MATRICES
# =============================================================================

cli::cli_h2("1) Loading & processing data")

# --- 1.1 CAFI community abundances -------------------------------------------
cli::cli_alert_info("Reading CAFI community file…")

cafi_path <- here::here(DATA_DIR, "MRB Amount",
                        "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")
if (!file.exists(cafi_path)) stop("Missing CAFI file at: ", cafi_path)

cafi_raw <- readr::read_csv(cafi_path, col_types = readr::cols())

# Basic hygiene: keep needed cols, clean IDs, drop empties
cafi_neat <- cafi_raw %>%
  dplyr::mutate(
    coral_id = strip_fe(as.character(coral_id)),
    species  = trimws(as.character(species))
  ) %>%
  dplyr::filter(!is.na(coral_id), !is.na(species))

# Build coral x species wide matrix of total counts
cli::cli_alert_info("Building coral × species abundance matrix…")
comm_abund <- cafi_neat %>%
  dplyr::group_by(coral_id, species) %>%
  dplyr::summarise(abundance = sum(count, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

# Meta per species BEFORE any filtering (so we can report what was dropped)
cli::cli_alert_info("Summarising species prevalence & totals (pre-filter)…")
species_full_meta <- cafi_neat %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    total_count = sum(count, na.rm = TRUE),
    n_corals    = dplyr::n_distinct(coral_id[count > 0]),
    .groups     = "drop"
  )

# Apply filtering thresholds (prevalence + total abundance)
cli::cli_alert_info("Applying filters: total_count >= {ABUND_MIN_COMM}, n_corals >= {PREV_MIN_COMM}")
species_retained <- species_full_meta %>%
  dplyr::filter(total_count >= ABUND_MIN_COMM, n_corals >= PREV_MIN_COMM)

if (nrow(species_retained) == 0) {
  stop("No species passed the filters. Check ABUND_MIN_COMM and PREV_MIN_COMM thresholds.")
}

# Subset wide matrix to retained species only (species-level filter)
comm_abund <- comm_abund %>%
  dplyr::select(coral_id, dplyr::all_of(species_retained$species))

# Report retention (species)
n_total    <- nrow(species_full_meta)
n_retained <- nrow(species_retained)
n_dropped  <- n_total - n_retained
cli::cli_alert_success("Species retained: {n_retained}/{n_total} (dropped: {n_dropped})")

# Optional: helpers
retained_species <- species_retained %>% dplyr::arrange(dplyr::desc(total_count))
dropped_species  <- setdiff(species_full_meta$species, retained_species$species)

# --- 1.2 Growth + physiology (filtered set) -----------------------------------
cli::cli_alert_info("Reading growth & physiology files…")

growth_file <- here::here(DATA_DIR, "processed", "coral_growth.csv")
physio_file <- here::here(OUT_DIR, "figures", "coral", "physio",
                          "physio_metrics_plus_growth_filtered.csv")

if (!file.exists(growth_file)) stop("Missing growth file at: ", growth_file)
if (!file.exists(physio_file)) stop("Missing physio file at: ", physio_file)

growth <- readr::read_csv(growth_file, show_col_types = FALSE) %>%
  dplyr::mutate(coral_id = strip_fe(as.character(coral_id)))
physio <- readr::read_csv(physio_file, show_col_types = FALSE) %>%
  dplyr::mutate(coral_id = strip_fe(as.character(coral_id)))

req_growth_cols <- c("coral_id", "growth_vol_b")
req_physio_cols <- c("coral_id", "protein_mg_cm2", "carb_mg_cm2",
                     "zoox_cells_cm2", "afdw_mg_cm2")
if (!all(req_growth_cols %in% names(growth))) {
  stop("Growth file missing columns: ", paste(setdiff(req_growth_cols, names(growth)), collapse = ", "))
}
if (!all(req_physio_cols %in% names(physio))) {
  stop("Physio file missing columns: ", paste(setdiff(req_physio_cols, names(physio)), collapse = ", "))
}

# --- NEW: Apply the ≥80% alive coral filter to the community matrix ----------
# The physio CSV you loaded is already filtered to the ≥80% alive set in your other script.
# Use the intersection of corals present in the filtered growth/physio to define keep_ids.
keep_ids <- intersect(growth$coral_id, physio$coral_id)

n_comm_before <- nrow(comm_abund)
dropped_coral_ids <- setdiff(comm_abund$coral_id, keep_ids)

comm_abund <- comm_abund %>%
  dplyr::filter(coral_id %in% keep_ids)

cli::cli_alert_success(
  "Alive-threshold filter applied to community matrix: kept {nrow(comm_abund)}/{n_comm_before} corals; dropped {length(dropped_coral_ids)}."
)
if (length(dropped_coral_ids)) {
  cli::cli_alert_info("Example dropped coral IDs: {paste(utils::head(dropped_coral_ids, 10), collapse=', ')}")
}

cli::cli_alert_success("Data loaded: {nrow(comm_abund)} corals × {ncol(comm_abund)-1} species (post-filter)")
cli::cli_alert_success("Growth rows: {nrow(growth)} | Physio rows: {nrow(physio)}")

# --- 1.3 Build analysis matrices (condition + community) ----------------------
cli::cli_alert_info("Constructing condition (z-scored) and community matrices…")

## 1.3a Condition matrix (growth + physio), z-scored
cond_df <- growth %>%
  dplyr::select(coral_id, growth_vol_b) %>%
  dplyr::inner_join(
    physio %>%
      dplyr::select(coral_id, protein_mg_cm2, carb_mg_cm2, zoox_cells_cm2, afdw_mg_cm2),
    by = "coral_id"
  )

cond_mat <- cond_df %>%
  dplyr::select(-coral_id) %>%
  as.matrix()
rownames(cond_mat) <- cond_df$coral_id
cond_mat_z <- scale(cond_mat)

## 1.3b Community matrices (RAW, sqrt, Hellinger)
comm_raw_mat <- comm_abund %>%
  tibble::column_to_rownames("coral_id") %>%
  as.matrix()
# Drop zero-variance columns that can break PCA
comm_raw_mat <- comm_raw_mat[, apply(comm_raw_mat, 2, var) > 0, drop = FALSE]

# Square-root transform (post-filter)
comm_sqrt_mat <- sqrt(comm_raw_mat)
comm_sqrt_mat <- comm_sqrt_mat[, apply(comm_sqrt_mat, 2, var) > 0, drop = FALSE]

# Hellinger transform
comm_mat <- vegan::decostand(comm_raw_mat, method = "hellinger")
comm_mat <- comm_mat[, apply(comm_mat, 2, var) > 0, drop = FALSE]

cli::cli_alert_success("Matrices ready: RAW={ncol(comm_raw_mat)} spp, SQRT={ncol(comm_sqrt_mat)} spp, HELL={ncol(comm_mat)} spp")


# =============================================================================
# 2) PCA helper: run PCA on community + condition, export summaries/plots
# =============================================================================

cli::cli_h2("2) Defining PCA helper")

# run_pca()
# - comm_mat: numeric matrix [coral x species] to PCA
# - label:    short label used in column names and file outputs ("SQRT_CS","HELLINGER","SQRT")
# - make_plots: logical; if TRUE, save scree + loading plots
# Notes:
#   * Community PCA centering/scaling is controlled by comm_center / comm_scale.
#   * Condition PCA uses the already z-scored matrix (no center/scale here).
run_pca <- function(comm_mat, label, cond_mat_z, make_plots = TRUE,
                    comm_center = TRUE, comm_scale = TRUE) {
  # Align rows by shared coral_id
  rn_comm <- rownames(comm_mat)
  rn_cond <- rownames(cond_mat_z)
  common  <- intersect(rn_comm, rn_cond)
  if (length(common) == 0L)
    stop("No overlapping coral_id between community and condition matrices.")
  if (length(common) < length(rn_comm) || length(common) < length(rn_cond)) {
    cli::cli_alert_warning(
      "Aligning by coral_id: keeping {length(common)} shared rows; ",
      "dropping {length(setdiff(rn_comm, common))} from community and ",
      "{length(setdiff(rn_cond, common))} from condition."
    )
  }
  comm_aligned <- comm_mat [common, , drop = FALSE]
  cond_aligned <- cond_mat_z[common, , drop = FALSE]
  
  # PCAs
  pca_comm <- prcomp(comm_aligned,  center = comm_center,  scale. = comm_scale)
  pca_cond <- prcomp(cond_aligned,  center = FALSE,        scale. = FALSE)
  
  # Scores (build with temp names, then rename to avoid NSE issues)
  scores_comm <- tibble::tibble(
    coral_id = rownames(pca_comm$x),
    PC1_tmp  = pca_comm$x[, 1],
    PC2_tmp  = pca_comm$x[, 2]
  )
  names(scores_comm)[names(scores_comm) == "PC1_tmp"] <- paste0(label, "_PC1")
  names(scores_comm)[names(scores_comm) == "PC2_tmp"] <- paste0(label, "_PC2")
  
  scores_cond <- tibble::tibble(
    coral_id = rownames(pca_cond$x),
    cond_PC1 = pca_cond$x[, 1],
    cond_PC2 = pca_cond$x[, 2]
  )
  
  # Optional plots
  if (isTRUE(make_plots)) {
    var_comm <- pca_comm$sdev^2 / sum(pca_comm$sdev^2)
    p_scree_c <- tibble::tibble(PC = seq_along(var_comm), var = var_comm) |>
      ggplot2::ggplot(ggplot2::aes(PC, var)) +
      ggplot2::geom_col() +
      ggplot2::geom_text(ggplot2::aes(label = scales::percent(var, .1)),
                         vjust = -0.3, size = 3) +
      ggplot2::labs(title = paste(label, "Community PCA Scree"), y = "Prop var") +
      theme_pub()
    show_and_save(p_scree_c, file.path(FIG_DIR, paste0("scree_COMM_", label)), 5, 4)
    
    var_cond <- pca_cond$sdev^2 / sum(pca_cond$sdev^2)
    p_scree_d <- tibble::tibble(PC = seq_along(var_cond), var = var_cond) |>
      ggplot2::ggplot(ggplot2::aes(PC, var)) +
      ggplot2::geom_col() +
      ggplot2::geom_text(ggplot2::aes(label = scales::percent(var, .1)),
                         vjust = -0.3, size = 3) +
      ggplot2::labs(title = "Condition PCA Scree", y = "Prop var") +
      theme_pub()
    show_and_save(p_scree_d, file.path(FIG_DIR, "scree_COND_z"), 5, 4)
    
    comm_load <- tibble::tibble(
      feature = rownames(pca_comm$rotation),
      loading = pca_comm$rotation[, 1]
    ) |>
      dplyr::arrange(dplyr::desc(abs(loading)))
    readr::write_csv(comm_load, file.path(TAB_DIR, paste0("PC1_loadings_COMM_", label, ".csv")))
    
    p_load <- comm_load |>
      dplyr::slice_head(n = TOP_N_LOAD) |>
      ggplot2::ggplot(ggplot2::aes(stats::reorder(feature, loading), loading, fill = loading > 0)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste(label, "COMM PC1 Loadings"),
                    x = NULL, y = expression(Loading~on~PC[1])) +
      theme_pub()
    show_and_save(p_load, file.path(FIG_DIR, paste0("PC1_loadings_COMM_", label)), 6, 5)
  }
  
  list(
    pca_comm    = pca_comm,
    scores_comm = scores_comm,
    scores_cond = scores_cond
  )
}



# --- Section 2 (helpers) -- after run_pca() ---
.align_to_condition <- function(res_obj, label_for_files) {
  df <- dplyr::inner_join(res_obj$scores_comm, res_obj$scores_cond, by = "coral_id")
  comm_pc1_col <- grep("_PC1$", names(res_obj$scores_comm), value = TRUE)
  cc <- suppressWarnings(cor(df[[comm_pc1_col]], df$cond_PC1, use = "pairwise.complete.obs"))
  
  if (!is.na(cc) && cc < 0) {
    res_obj$pca_comm$rotation[, 1] <- -res_obj$pca_comm$rotation[, 1]
    res_obj$pca_comm$x[, 1]        <- -res_obj$pca_comm$x[, 1]
    res_obj$scores_comm[[comm_pc1_col]] <- -res_obj$scores_comm[[comm_pc1_col]]
    
    load_tbl <- tibble::tibble(
      feature = rownames(res_obj$pca_comm$rotation),
      loading = res_obj$pca_comm$rotation[, 1]
    ) |>
      dplyr::arrange(dplyr::desc(abs(loading)))
    readr::write_csv(load_tbl, file.path(TAB_DIR, paste0("PC1_loadings_COMM_", label_for_files, "_aligned.csv")))
    p_load <- load_tbl |>
      dplyr::slice_head(n = TOP_N_LOAD) |>
      ggplot2::ggplot(ggplot2::aes(stats::reorder(feature, loading), loading, fill = loading > 0)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste(label_for_files, "COMM PC1 Loadings (aligned to cond_PC1)"),
        x = NULL, y = expression(Loading~on~PC[1])
      ) + theme_pub()
    show_and_save(p_load, file.path(FIG_DIR, paste0("PC1_loadings_COMM_", label_for_files, "_aligned")), 6, 5)
  }
  res_obj
}


# =============================================================================
# 3) Community PCA variants & orientation alignment to Condition PC1
# =============================================================================
cli::cli_h2("3) Running PCAs & aligning PC1 to Condition PC1")

# 3.1 Run PCAs
res_sqrt_cs <- run_pca(comm_sqrt_mat, label = "SQRT_CS",
                       cond_mat_z = cond_mat_z, make_plots = TRUE,
                       comm_center = TRUE,  comm_scale = TRUE)

res_hell    <- run_pca(comm_mat,      label = "HELLINGER",
                       cond_mat_z = cond_mat_z, make_plots = TRUE,
                       comm_center = TRUE,  comm_scale = TRUE)

res_sqrt    <- run_pca(comm_sqrt_mat, label = "SQRT",
                       cond_mat_z = cond_mat_z, make_plots = TRUE,
                       comm_center = FALSE, comm_scale = FALSE)

res_raw     <- run_pca(comm_raw_mat,  label = "RAW",
                       cond_mat_z = cond_mat_z, make_plots = TRUE,
                       comm_center = TRUE,  comm_scale = TRUE)

# 3.2 Align PC1 signs so cor(Community PC1, cond_PC1) ≥ 0
res_sqrt_cs <- .align_to_condition(res_sqrt_cs, "SQRT_CS")
res_hell    <- .align_to_condition(res_hell,    "HELLINGER")
res_sqrt    <- .align_to_condition(res_sqrt,    "SQRT")
res_raw     <- .align_to_condition(res_raw,     "RAW")



# =============================================================================
# 4) Community PC1 → Condition PC1 (SQRT_CS, HELLINGER, SQRT): LMMs & plots
# =============================================================================

cli::cli_h2("4) Fitting LMMs and summarizing effects (SQRT_CS, HELLINGER, SQRT)")

# 4.1 Fit LMMs for all three
# Model: cond_PC1 ~ community_PC1 + (1 | reef_id)
# Notes:
# - 'reef' column from `growth` is used as reef_id (joined by coral_id)
# - Keep only rows present in both scores tables and growth

res_list <- list(SQRT_CS = res_sqrt_cs, HELLINGER = res_hell, SQRT = res_sqrt)

models <- purrr::imap(res_list, function(obj, lbl) {
  df <- obj$scores_comm %>%
    dplyr::inner_join(obj$scores_cond, by = "coral_id") %>%
    dplyr::left_join(growth %>% dplyr::select(coral_id, reef_id = reef), by = "coral_id") %>%
    dplyr::filter(!is.na(reef_id))
  
  # find the community PC1 column (SQRT_CS_PC1, HELLINGER_PC1, or SQRT_PC1)
  comm_col <- grep(paste0("^", lbl, "_PC1$"), names(df), value = TRUE)
  stopifnot(length(comm_col) == 1)
  
  # build formula: cond_PC1 ~ <comm_col> + (1 | reef_id)
  fm <- as.formula(paste0("cond_PC1 ~ ", comm_col, " + (1|reef_id)"))
  
  mm <- lmerTest::lmer(
    fm,
    data = df,
    control = lme4::lmerControl(check.conv.singular = "ignore")
  )
  
  # save per-model tables
  fix_tbl <- broom.mixed::tidy(mm, effects = "fixed", conf.int = TRUE)
  anova_tbl <- car::Anova(mm, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")
  
  readr::write_csv(fix_tbl,  file.path(TAB_DIR, paste0("LMM_", lbl, "_fixed.csv")))
  readr::write_csv(anova_tbl, file.path(TAB_DIR, paste0("LMM_", lbl, "_anova.csv")))
  
  list(model = mm, data = df, xcol = comm_col)
})

cli::cli_alert_success("Fitted LMMs for SQRT_CS, HELLINGER, and SQRT.")

# 4.2 Plot marginal fits side-by-side
cli::cli_h2("4.2) Plotting marginal fits")

plots <- purrr::imap(models, function(entry, lbl) {
  df   <- entry$data
  xcol <- entry$xcol[[1]]
  
  # fixed-effect line with observed points
  fe <- lme4::fixef(entry$model)
  df <- df %>% dplyr::mutate(pred = fe["(Intercept)"] + fe[xcol] * .data[[xcol]])
  
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[xcol]], y = cond_PC1)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = pred), linewidth = 1) +
    ggplot2::labs(
      title = paste(lbl, "community → condition"),
      x     = paste(lbl, "PC1"),
      y     = "Condition PC1"
    ) +
    theme_pub()
})

combo <- plots$SQRT_CS | plots$HELLINGER | plots$SQRT

show_and_save(
  combo,
  file.path(FIG_DIR, "LMM_COMM_vs_COND_comparison"),
  width = 18, height = 5
)

# 4.3 Tidy one-row-per-model summary and display
cli::cli_h2("4.3) Summarizing fixed effects")

lmm_summary <- purrr::imap_dfr(models, function(entry, lbl) {
  fix_tbl <- broom.mixed::tidy(entry$model, effects = "fixed", conf.int = TRUE)
  row <- dplyr::filter(fix_tbl, term == entry$xcol)  # SQRT_CS_PC1 / HELLINGER_PC1 / SQRT_PC1
  
  tibble::tibble(
    Matrix    = lbl,
    Estimate  = row$estimate,
    Std_Error = row$std.error,
    DF        = row$df,
    t_stat    = row$statistic,
    p_value   = row$p.value,
    CI_lower  = row$conf.low,
    CI_upper  = row$conf.high
  )
})

readr::write_csv(lmm_summary, file.path(TAB_DIR, "LMM_COMM_vs_COND_summary.csv"))

# Optional pretty table
lmm_summary %>%
  dplyr::mutate(
    Estimate  = signif(Estimate,  3),
    Std_Error = signif(Std_Error, 3),
    t_stat    = signif(t_stat,    3),
    p_value   = signif(p_value,   3),
    CI_lower  = signif(CI_lower,  3),
    CI_upper  = signif(CI_upper,  3)
  ) %>%
  gt::gt() %>%
  gt::tab_header(
    title    = gt::md("**LMM Summary: Community PC1 → Condition PC1**"),
    subtitle = "SQRT_CS (sqrt + center/scale), HELLINGER, SQRT (sqrt only)"
  ) %>%
  gt::cols_label(
    Matrix    = "Matrix",
    Estimate  = "Estimate",
    Std_Error = "Std. Error",
    DF        = "df",
    t_stat    = "t-stat",
    p_value   = "p-value",
    CI_lower  = "CI Lower",
    CI_upper  = "CI Upper"
  ) %>%
  gt::fmt_number(
    columns  = c("Estimate", "Std_Error", "t_stat", "CI_lower", "CI_upper"),
    decimals = 3
  ) %>%
  gt::fmt_number(columns = "p_value", decimals = 3) %>%
  print()



# --- 9-panel PCA overview (3 rows x 3 cols) ----------------------------------

library(ggplot2)
library(dplyr)
library(patchwork)

library(ggplot2)
library(dplyr)
library(patchwork)
library(forcats)

# Colors already defined at top: cols_trt = TREATMENT_COLORS

.make_triptych <- function(res_obj, label, top_n = TOP_N_LOAD, meta = NULL, cols_trt = cols_trt) {
  # Scree (community)
  # Scree (community)
  var_comm <- (res_obj$pca_comm$sdev ^ 2) / sum(res_obj$pca_comm$sdev ^ 2)
  ymax <- max(var_comm) + 0.10  # add 0.1 headroom
  
  p_scree <- tibble::tibble(PC = seq_along(var_comm), var = var_comm) |>
    ggplot(aes(PC, var)) +
    geom_col() +
    geom_text(aes(label = scales::percent(var, .1)), vjust = -0.3, size = 3) +
    labs(title = paste(label, "Community PCA Scree"), y = "Prop var", x = "PC") +
    coord_cartesian(ylim = c(0, ymax)) +  # ensure extra space
    theme_pub()
  
  # Loadings (PC1) — top N by |loading|
  load_tbl <- tibble::tibble(
    feature = rownames(res_obj$pca_comm$rotation),
    loading = res_obj$pca_comm$rotation[, 1]
  ) |>
    arrange(desc(abs(loading))) |>
    slice_head(n = top_n)
  
  p_load <- load_tbl |>
    ggplot(aes(x = stats::reorder(feature, loading), y = loading)) +
    geom_segment(aes(xend = feature, y = 0, yend = loading), colour = "grey55", linewidth = 0.5) +
    geom_point(size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
    labs(title = paste(label, "COMM PC1 Loadings"), x = NULL, y = expression(Loading~on~PC[1])) +
    theme_pub(base_size = 10)
  
  # Community PC1 vs Condition PC1 (scatter; color by treatment)
  comm_col <- grep(paste0("^", label, "_PC1$"), names(res_obj$scores_comm), value = TRUE)
  stopifnot(length(comm_col) == 1)
  
  df_sc <- res_obj$scores_comm |>
    inner_join(res_obj$scores_cond, by = "coral_id")
  
  # bring treatment if needed
  if (!"treatment" %in% names(df_sc)) {
    if (is.null(meta) || !"treatment" %in% names(meta)) {
      stop("treatment not found in scores; please pass meta = meta_df with coral_id + treatment.")
    }
    df_sc <- df_sc |>
      left_join(meta |> dplyr::select(coral_id, treatment), by = "coral_id")
  }
  
  # enforce factor + palette levels
  df_sc <- df_sc |>
    mutate(treatment = factor(as.character(treatment), levels = names(cols_trt)))
  
  p_scatter <- ggplot(df_sc, aes(x = .data[[comm_col]], y = cond_PC1, color = treatment)) +
    geom_point(alpha = 0.8, size = 1.9) +
    geom_smooth(method = "lm", se = TRUE, colour = "black") +
    scale_color_manual(values = cols_trt, drop = FALSE, name = "Treatment") +
    labs(title = paste(label, "Community PC1 vs Condition PC1"),
         x = paste(label, "PC1"), y = "Condition PC1") +
    theme_pub()
  
  list(scree = p_scree, load = p_load, scatter = p_scatter)
}

# After reading `growth` (and before calling .make_triptych)
meta_df <- growth %>%
  dplyr::transmute(
    coral_id,
    treatment = factor(as.character(treatment), levels = c("1","3","6"))
  )



# Build triptychs for the three community transforms
tri_sqrtcs   <- .make_triptych(res_sqrt_cs, "SQRT_CS", meta = meta_df, cols_trt = cols_trt)
tri_hell     <- .make_triptych(res_hell,    "HELLINGER", meta = meta_df, cols_trt = cols_trt)
tri_sqrtonly <- .make_triptych(res_sqrt,    "SQRT",      meta = meta_df, cols_trt = cols_trt)

# Arrange into a 3 x 3 grid
p9 <- (tri_sqrtcs$scree | tri_sqrtcs$load | tri_sqrtcs$scatter) /
  (tri_hell$scree   | tri_hell$load   | tri_hell$scatter)   /
  (tri_sqrtonly$scree | tri_sqrtonly$load | tri_sqrtonly$scatter) +
  plot_annotation(tag_levels = "A")

# Show & save
show_and_save(p9, file.path(FIG_DIR, "PCA_9panel_overview"), width = 15, height = 12, dpi = 600)


# Get triptych for SQRT_CS
tri_sqrtcs <- .make_triptych(res_sqrt_cs, "SQRT_CS", meta = meta_df, cols_trt = cols_trt)

# Panel A: COMM PC1 loadings
p_load <- tri_sqrtcs$load +
  labs(
    title = NULL,
    y = expression(Loading~on~PC[1][CAFI]),
    x = NULL
  )

# Panel B: Scatter
p_scatter <- tri_sqrtcs$scatter +
  labs(
    title = NULL,
    x = expression(Community~composition~(PC[1][CAFI])),
    y = expression(Coral~condition~(PC[1][coral]))
  )

# Combine with labels A (left) and B (right)
p_bc_sqrtcs <- p_load + p_scatter +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(legend.position = "top")

# Save
show_and_save(
  p_bc_sqrtcs,
  file.path(FIG_DIR, "SQRT_CS_PC1loadings_and_scatter"),
  width = 10, height = 5, dpi = 600
)

# =============================================================================
# 5) Publication-ready 3-panel figures (HELLINGER & RAW)
# =============================================================================

cli::cli_h2("5) Building publication-ready 3-panel figures")

library(ggplot2)
library(patchwork)
library(forcats)

# 5.0 Common styling & colors ---------------------------------------------------
# Colors already defined at top: cols_trt = TREATMENT_COLORS

# 5.1 Condition loadings once (shared across RAW/HELL) --------------------------
# Flip sign so positive loadings align with "better" condition (as discussed)
pca_cond_z <- prcomp(cond_mat_z, center = FALSE, scale. = FALSE)

df_cond_load <- tibble::tibble(
  feature = rownames(pca_cond_z$rotation),
  loading = -pca_cond_z$rotation[, 1]
) |>
  dplyr::arrange(dplyr::desc(loading)) |>   # strongest positive first
  dplyr::mutate(
    feature_label = forcats::fct_inorder(dplyr::recode(
      feature,
      carb_mg_cm2    = "Carbohydrate~(mg~cm^{-2})",
      protein_mg_cm2 = "Protein~(mg~cm^{-2})",
      afdw_mg_cm2    = "Tissue~biomass~(AFDW~mg~cm^{-2})",
      zoox_cells_cm2 = "Zooxanthellae~(cells~cm^{-2})",
      growth_vol_b   = "Growth"
    ))
  )

# ---- Small helpers ------------------------------------------------------------

# Build per-matrix data: scores + top-|loading| species with parse-ready labels
.make_row_data <- function(res, label) {
  df_pca <- res$scores_comm %>%
    dplyr::inner_join(res$scores_cond, by = "coral_id") %>%
    dplyr::left_join(growth %>% dplyr::select(coral_id, treatment), by = "coral_id") %>%
    dplyr::mutate(treatment = factor(treatment, levels = c("1","3","6")))
  
  df_comm_load <- tibble::tibble(
    feature = rownames(res$pca_comm$rotation),
    loading = res$pca_comm$rotation[, 1]
  ) %>%
    dplyr::arrange(dplyr::desc(loading)) %>%        # strongest positive first (top)
    dplyr::slice_head(n = TOP_N_LOAD) %>%
    dplyr::mutate(
      feature_label = paste0("italic(", gsub(" ", "~", feature), ")"),
      feature_label = forcats::fct_inorder(feature_label)
    )
  
  list(df_pca = df_pca, df_comm_load = df_comm_load)
}

# Scatter panel (no per-panel title; tags added via plot_annotation)
.make_scatter <- function(df_pca, xvar) {
  ggplot(df_pca, aes(x = .data[[xvar]], y = cond_PC1)) +
    geom_point(aes(fill = treatment),
               shape = 21, colour = "darkgray", size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, colour = "black") +
    scale_fill_manual(values = cols_trt) +
    labs(x = expression(Community~PC[1]), y = expression(Condition~PC[1]), fill = "Treatment") +
    theme_pub()
}

# Community loadings panel: positives at top (after coord_flip)
.make_comm_load <- function(df_comm_load) {
  df_comm_load <- df_comm_load %>%
    dplyr::arrange(dplyr::desc(loading)) %>%    # strongest positive first
    dplyr::mutate(feature_label = forcats::fct_inorder(feature_label))
  
  ggplot(df_comm_load, aes(x = feature_label, y = loading)) +
    geom_segment(aes(xend = feature_label, y = 0, yend = loading),
                 colour = "darkgray", linewidth = 0.4) +
    geom_point(colour = "black", size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05)),  # extra 5% space on both ends
      breaks = function(lims) unique(sort(c(scales::pretty_breaks(4)(lims), 0)))
    )+
    labs(x = NULL, y = expression(Loading~on~PC[1])) +
    theme_pub(base_size = 10)
}

# Condition loadings panel: already globally built & ordered
.make_cond_load <- function(df_cond_load) {
  df_cond_load <- df_cond_load %>%
    dplyr::arrange(dplyr::desc(loading)) %>%              # positives first
    dplyr::mutate(feature_label = forcats::fct_inorder(feature_label))
  
  ggplot(df_cond_load, aes(x = feature_label, y = loading)) +
    geom_segment(aes(xend = feature_label, y = 0, yend = loading),
                 colour = "darkgray", linewidth = 0.4) +
    geom_point(colour = "black", size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05)),           # extra space so dots don’t hug axes
      breaks = function(lims) unique(sort(c(scales::pretty_breaks(4)(lims), 0)))
    ) +
    labs(x = NULL, y = expression(Loading~on~PC[1])) +
    theme_pub(base_size = 10)
}


# Combine three panels and save with a single tag layer (A/B/C)
.combine_and_save <- function(pA, pB, pC, out_stub, width = 13, height = 5) {
  combo <- (pA | pB | pC) +
    patchwork::plot_layout(ncol = 3, widths = c(1.1, 1, 1)) +
    patchwork::plot_annotation(tag_levels = "A")
  show_and_save(combo, file.path(FIG_DIR, out_stub), width = width, height = height)
}

# 5.3 HELLINGER 3-panel ---------------------------------------------------------
cli::cli_h2("5.3) HELLINGER 3-panel")

hell <- .make_row_data(res_hell, "HELLINGER")
pA_hell <- .make_scatter(hell$df_pca, "HELLINGER_PC1")
pB_hell <- .make_comm_load(hell$df_comm_load)
pC_hell <- .make_cond_load(df_cond_load)

.combine_and_save(pA_hell, pB_hell, pC_hell, "PCA_LOADINGS_HELLINGER_3panel_clean")

# 5.4 RAW 3-panel ---------------------------------------------------------------
cli::cli_h2("5.4) RAW 3-panel")

raw <- .make_row_data(res_raw, "RAW")
pA_raw <- .make_scatter(raw$df_pca, "RAW_PC1") +
  labs(x = expression(Community~PC[1]~"(RAW)"))  # small label tweak
pB_raw <- .make_comm_load(raw$df_comm_load)
pC_raw <- .make_cond_load(df_cond_load)

.combine_and_save(pA_raw, pB_raw, pC_raw, "PCA_LOADINGS_RAW_3panel_clean")





#THIS IS THE A PUBLICATION FIGURE

# =============================================================================
# 5.5) RAW-only 2-panel: (A) Community PC1 vs Condition PC1,
#                        (B) Top species loadings on Community PC1
# =============================================================================

cli::cli_h2("5.5) RAW-only 2-panel (A/B)")

library(ggplot2)
library(patchwork)
library(forcats)

# Colors already defined at top: cols_trt = TREATMENT_COLORS

# Colors for taxonomic groups (distinct from treatment colors to avoid confusion)
cols_taxon <- c(
  "Fishes"        = "#2E86AB",  # blue
  "Shrimps/Crabs" = "#A23B72",  # magenta/pink
  "Snails"        = "#F18F01"   # orange
)

# Build data for RAW: scatter inputs + top-|loading| species with parse-ready labels
# Now includes taxonomic group assignment for coloring
.make_row_data_raw_2p <- function(res_raw) {
  df_pca <- res_raw$scores_comm %>%
    dplyr::inner_join(res_raw$scores_cond, by = "coral_id") %>%
    dplyr::left_join(growth %>% dplyr::select(coral_id, treatment), by = "coral_id") %>%
    dplyr::mutate(treatment = factor(treatment, levels = c("1","3","6")))

  df_comm_load <- tibble::tibble(
    feature = rownames(res_raw$pca_comm$rotation),
    loading = res_raw$pca_comm$rotation[, 1]
  ) %>%
    dplyr::arrange(dplyr::desc(loading)) %>%      # strongest positive first
    dplyr::slice_head(n = TOP_N_LOAD) %>%
    dplyr::mutate(
      feature_label = paste0("italic(", gsub(" ", "~", feature), ")"),
      # Assign taxonomic group based on known species
      taxon_group = dplyr::case_when(
        # Fishes (Actinopterygii)
        grepl("Dascyllus|Caracanthus|Halichoeres|Paragobiodon|Gobiodon|Plectroglyphidodon|Stegastes|Pomacentrus", feature) ~ "Fishes",
        # Shrimps and Crabs (Malacostraca)
        grepl("Trapezia|Tetralia|Alpheus|Synalpheus|Periclimenes|Harpiliopsis|Calcinus|Fennera|Galathea|Thor|Athanas|Cinetorhynchus|Saron|Urocaridella|Pagurixus|Luniella", feature) ~ "Shrimps/Crabs",
        # Snails (Gastropoda)
        grepl("Morula|Mitrella|Galeropsis|Macteola|Chlorodiella|Pascula|Cellana|Menaethius|Vexillum|Strigatella|Apatasia|Coralliophila|Drupella|Quoyula", feature) ~ "Snails",
        TRUE ~ "Fishes"  # Default fallback
      ),
      taxon_group = factor(taxon_group, levels = c("Fishes", "Shrimps/Crabs", "Snails"))
    )

  # Option 1: Keep ordering by loading (Craig's first suggestion)
  df_comm_load <- df_comm_load %>%
    dplyr::mutate(feature_label = forcats::fct_inorder(feature_label))

  # Option 2 (alternative): Order by group then by loading within group
  # Uncomment below to use Craig's second suggestion
  # df_comm_load <- df_comm_load %>%
  #   dplyr::arrange(taxon_group, dplyr::desc(loading)) %>%
  #   dplyr::mutate(feature_label = forcats::fct_inorder(feature_label))

  list(df_pca = df_pca, df_comm_load = df_comm_load)
}

# Panel A: scatter (RAW_PC1 vs Condition PC1)
# Panel A: scatter (RAW_PC1 vs Condition PC1)
.make_scatter_raw_2p <- function(df_pca) {
  ggplot(df_pca, aes(x = RAW_PC1, y = cond_PC1)) +
    geom_point(aes(fill = treatment), shape = 21, colour = "darkgray",
               size = 3, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, colour = "black") +
    scale_fill_manual(values = cols_trt, name = "Coral number") +
    labs(
      x = expression(PC1[CAFI]),   # PC1 with subscript "CAFI"
      y = expression(PC1[coral])   # PC1 with subscript "coral"
    ) +
    theme_pub() +
    theme(legend.position = "top")
}

# Panel B: top loadings on Community PC1 (positives at top)
# Now colored by taxonomic group (per Craig Osenberg's suggestion)
.make_comm_load_raw_2p <- function(df_comm_load) {
  df_comm_load <- df_comm_load %>%
    dplyr::arrange(dplyr::desc(loading)) %>%
    dplyr::mutate(feature_label = forcats::fct_inorder(feature_label))

  ggplot(df_comm_load, aes(x = feature_label, y = loading)) +
    geom_segment(aes(xend = feature_label, y = 0, yend = loading, colour = taxon_group),
                 linewidth = 0.5) +
    geom_point(aes(fill = taxon_group), shape = 21, colour = "black", size = 2.5, stroke = 0.3) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
    scale_x_discrete(labels = function(x) parse(text = x)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.08, 0.08)),   # extra space so dots don't sit on axes
      breaks = function(lims) unique(sort(c(scales::pretty_breaks(4)(lims), 0)))
    ) +
    scale_fill_manual(values = cols_taxon, name = "Taxon") +
    scale_colour_manual(values = cols_taxon, guide = "none") +  # match segment colors, hide duplicate legend
    labs(x = NULL, y = expression(Loading~on~PC1[CAFI])) +
    theme_pub(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8, margin = margin(r = 2)),
      plot.margin = margin(5, 5, 5, 10),
      legend.position = "top",
      legend.box.margin = margin(b = -5)
    )
}

# Build panels
raw2 <- .make_row_data_raw_2p(res_raw)
pA   <- .make_scatter_raw_2p(raw2$df_pca)
pB   <- .make_comm_load_raw_2p(raw2$df_comm_load)

# Combine + single set of panel tags
raw_two_panel <- (pB | pA) +
  patchwork::plot_layout(ncol = 2, widths = c(1.1, 1)) +
  patchwork::plot_annotation(tag_levels = "A")  # tags: A, B

# Show & save
show_and_save(
  raw_two_panel,
  file.path(FIG_DIR, "PCA_LOADINGS_RAW_2panel_clean"),
  width = 9.5, height = 6
)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(OUT_DIR, "figures", "publication-figures", "PCA_LOADINGS_RAW_2panel_clean.pdf"),
  raw_two_panel, width = 9.5, height = 6, dpi = 600, bg = "white"
)
ggsave(
  file.path(OUT_DIR, "figures", "publication-figures", "PCA_LOADINGS_RAW_2panel_clean.png"),
  raw_two_panel, width = 9.5, height = 6, dpi = 600, bg = "white"
)







# --- Bootstrap for §5.6 if §6.9 hasn't run yet --------------------------------
if (!exists("master_df") || !exists("species_name_map") || !exists("physio_vars")) {
  cli::cli_alert_warning("Bootstrapping master_df/species_name_map/physio_vars for §5.6")
  
  # 1) Required inputs from earlier sections must exist:
  stopifnot(exists("comm_raw_mat"), exists("growth"), exists("physio"))
  
  # 2) Clean species names (match §6.9 helpers)
  clean_species_name <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)   # replace spaces/punct with _
    x <- gsub("^_|_$", "", x)            # trim leading/trailing _
    x
  }
  
  # 3) Community counts as a data.frame with coral_id
  comm_df_raw <- comm_raw_mat |>
    as.data.frame() |>
    tibble::rownames_to_column("coral_id")
  
  orig_names  <- colnames(comm_raw_mat)
  clean_names <- clean_species_name(orig_names)
  colnames(comm_df_raw) <- c("coral_id", clean_names)
  
  species_name_map <- tibble::tibble(
    species_original = orig_names,
    species_clean    = clean_names
  )
  
  # 4) Traits/meta (match columns used elsewhere)
  #    - growth has coral_id, reef, treatment, growth_vol_b
  #    - physio has coral_id, protein_mg_cm2, carb_mg_cm2, zoox_cells_cm2, afdw_mg_cm2
  traits_df <- growth %>%
    dplyr::select(coral_id, reef_id = reef, treatment, growth_vol_b) %>%
    dplyr::inner_join(
      physio %>% dplyr::select(coral_id, protein_mg_cm2, carb_mg_cm2,
                               zoox_cells_cm2, afdw_mg_cm2),
      by = "coral_id"
    ) %>%
    dplyr::mutate(treatment = factor(as.character(treatment), levels = c("1","3","6")))
  
  master_df <- traits_df %>%
    dplyr::inner_join(comm_df_raw, by = "coral_id")
  
  # Optional: enforce the ≥80% alive filter if keep_ids_alive is available
  if (exists("keep_ids_alive")) {
    n0 <- nrow(master_df)
    master_df <- master_df %>% dplyr::filter(coral_id %in% keep_ids_alive)
    cli::cli_alert_info("Applied alive filter to master_df: kept {nrow(master_df)}/{n0} corals.")
  }
  
  # 5) Physio vars canon (no percent_alive in §5.6/§7)
  physio_vars <- c("growth_vol_b","protein_mg_cm2","carb_mg_cm2",
                   "zoox_cells_cm2","afdw_mg_cm2")
}

# Sanity check (original guard)
stopifnot(exists("master_df"), exists("physio_vars"), exists("species_name_map"),
          all(physio_vars %in% names(master_df)))




# ============================================================================
# 5.6) Species–performance correlations (raw & sqrt abundance)  [PC1 + sorting]
#       - Computes PC1 if missing (from physio_vars; centered & scaled)
#       - Coarse taxon groups: Fish / Shrimps+Crabs / Snails / Other
#       - Heatmap columns = metrics (incl. PC1); rows sorted within group by r(PC1) ↓
# ============================================================================

cli::cli_h2("5.6) Species–performance correlations (raw & sqrt) with PC1 & taxon-group ordering")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  library(stringr); library(readr)
})

# --- Required inputs ----------------------------------------------------------
stopifnot(
  exists("master_df"),
  exists("physio_vars"),
  exists("species_name_map"),
  exists("TAB_DIR"), exists("FIG_DIR")
)

# --- Ensure PC1(coral) exists -------------------------------------------------
pc1_coral_col <- dplyr::case_when(
  "PC1_physio_growth" %in% names(master_df) ~ "PC1_physio_growth",
  "PC1coral"          %in% names(master_df) ~ "PC1coral",
  "PC1_coral"         %in% names(master_df) ~ "PC1_coral",
  TRUE ~ NA_character_
)

if (is.na(pc1_coral_col)) {
  cli::cli_alert_info("Computing PC1(coral) from physio_vars (centered & scaled).")
  keep_cols <- intersect(physio_vars, names(master_df))
  stopifnot(length(keep_cols) >= 2)
  complete_idx <- stats::complete.cases(master_df[, keep_cols])
  
  pca_fit <- prcomp(master_df[complete_idx, keep_cols, drop = FALSE],
                    center = TRUE, scale. = TRUE)
  
  PC1_scores <- rep(NA_real_, nrow(master_df))
  PC1_scores[which(complete_idx)] <- pca_fit$x[, 1]
  
  # Orient PC1 to align positively with growth if available
  if ("growth_vol_b" %in% keep_cols) {
    sgn <- sign(stats::cor(PC1_scores[complete_idx],
                           master_df$growth_vol_b[complete_idx],
                           use = "complete.obs"))
    if (is.finite(sgn) && sgn < 0) PC1_scores <- -PC1_scores
  }
  
  master_df$PC1coral <- PC1_scores
  pc1_coral_col <- "PC1coral"
  
  # Save loadings for PC1
  load_tab <- tibble::tibble(
    variable = rownames(pca_fit$rotation),
    PC1_loading = as.numeric(pca_fit$rotation[,1])
  )
  out_load_csv <- file.path(TAB_DIR, "pc1_coral_loadings.csv")
  readr::write_csv(load_tab, out_load_csv)
  cli::cli_alert_success("Computed PC1(coral) and saved loadings → {out_load_csv}")
}

# --- Species columns ----------------------------------------------------------
species_cols_clean <- intersect(names(master_df), species_name_map$species_clean)
if (!length(species_cols_clean)) stop("No species columns found in master_df.")

# --- Metrics to correlate (PC1 + physio) -------------------------------------
metrics_use <- unique(c(pc1_coral_col, intersect(physio_vars, names(master_df))))

# --- Helper: compute Pearson correlations ------------------------------------
.compute_corr <- function(df, species_cols, metrics, transform = c("raw","sqrt")) {
  transform <- match.arg(transform)
  df_use <- df
  if (transform == "sqrt") for (sp in species_cols) df_use[[sp]] <- sqrt(df_use[[sp]])
  
  tidyr::crossing(species = species_cols, metric = metrics) %>%
    purrr::pmap_dfr(function(species, metric) {
      x <- df_use[[species]]; y <- df_use[[metric]]
      keep <- is.finite(x) & is.finite(y); n <- sum(keep)
      if (n >= 3) {
        ct <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "pearson"))
        tibble::tibble(species_clean = species, metric = metric, n = n,
                       r = unname(ct$estimate), p_value = ct$p.value)
      } else {
        tibble::tibble(species_clean = species, metric = metric, n = n,
                       r = NA_real_, p_value = NA_real_)
      }
    }) %>%
    dplyr::mutate(transform = transform)
}

corr_raw  <- .compute_corr(master_df, species_cols_clean, metrics_use, "raw")
corr_sqrt <- .compute_corr(master_df, species_cols_clean, metrics_use, "sqrt")

# --- Coarse taxon grouping from CAFI file ------------------------------------
tax_path <- here::here("data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")

# Robust column picker
cafi_head <- readr::read_csv(tax_path, n_max = 5, show_col_types = FALSE)
.pick_col <- function(df, candidates) {
  nx <- names(df); hit <- candidates[tolower(candidates) %in% tolower(nx)]
  if (length(hit)) nx[match(tolower(hit[1]), tolower(nx))] else NA_character_
}
genus_col   <- .pick_col(cafi_head,  c("genus","Genus"))
species_col <- .pick_col(cafi_head,  c("species","Species","epithet","species_epithet","sp","Sp"))
class_col   <- .pick_col(cafi_head,  c("class","Class"))
order_col   <- .pick_col(cafi_head,  c("order","Order"))
family_col  <- .pick_col(cafi_head,  c("family","Family"))

need_cols <- unique(na.omit(c(genus_col, species_col, class_col, order_col, family_col)))
cafi <- readr::read_csv(tax_path, show_col_types = FALSE) %>% dplyr::select(dplyr::all_of(need_cols))

cafi_tax_lookup <- cafi %>%
  dplyr::mutate(
    genus      = if (!is.na(genus_col))   .data[[genus_col]]   else NA_character_,
    species    = if (!is.na(species_col)) .data[[species_col]] else NA_character_,
    class_raw  = if (!is.na(class_col))   .data[[class_col]]   else NA_character_,
    order_raw  = if (!is.na(order_col))   .data[[order_col]]   else NA_character_,
    family_raw = if (!is.na(family_col))  .data[[family_col]]  else NA_character_
  ) %>%
  dplyr::mutate(
    species_original = stringr::str_squish(paste(genus, species)),
    norm_name        = tolower(species_original),
    genus_lc         = tolower(genus),
    family_lc        = tolower(family_raw),
    class_raw        = dplyr::na_if(class_raw, ""),
    order_raw        = dplyr::na_if(order_raw, "")
  ) %>%
  dplyr::filter(!is.na(norm_name), nzchar(norm_name)) %>%
  dplyr::mutate(
    class_inferred = dplyr::coalesce(
      class_raw,
      dplyr::case_when(
        order_raw %in% c("Perciformes","Scorpaeniformes","Eupercaria incertae sedis") ~ "Actinopterygii",
        order_raw %in% c("Decapoda","Tanaidacea") ~ "Malacostraca",
        order_raw %in% c("Neogastropoda","Caenogastropoda","[unassigned] Caenogastropoda","Littorinimorpha","Trochida") ~ "Gastropoda",
        TRUE ~ NA_character_
      )
    ),
    class_inferred = dplyr::case_when(
      class_inferred %in% c("Teleostei") ~ "Actinopterygii",
      TRUE ~ class_inferred
    ),
    taxon_group = dplyr::case_when(
      class_inferred == "Actinopterygii" ~ "Fishes",
      class_inferred == "Malacostraca"   ~ "Shrimps/Crabs",
      class_inferred == "Gastropoda"     ~ "Snails",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::select(norm_name, genus_lc, family_lc, taxon_group) %>%
  dplyr::distinct()

genus_to_group  <- cafi_tax_lookup %>% dplyr::distinct(genus_lc,  taxon_group) %>% dplyr::filter(!is.na(genus_lc),  nzchar(genus_lc))
family_to_group <- cafi_tax_lookup %>% dplyr::distinct(family_lc, taxon_group) %>% dplyr::filter(!is.na(family_lc), nzchar(family_lc))

# --- Combine correlations; attach names & taxon_group ------------------------
corr_all <- dplyr::bind_rows(corr_raw, corr_sqrt) %>%
  dplyr::left_join(species_name_map, by = "species_clean") %>%
  dplyr::mutate(
    species_original = dplyr::coalesce(species_original, gsub("_"," ", species_clean)),
    norm_name        = tolower(stringr::str_squish(species_original))
  ) %>%
  dplyr::left_join(cafi_tax_lookup, by = "norm_name") %>%
  dplyr::mutate(genus_from_name_lc = sub(" .*", "", norm_name)) %>%
  dplyr::left_join(genus_to_group  %>% dplyr::rename(taxon_from_genus  = taxon_group),
                   by = c("genus_from_name_lc" = "genus_lc")) %>%
  dplyr::mutate(taxon_group = dplyr::coalesce(taxon_group, taxon_from_genus)) %>%
  dplyr::select(-taxon_from_genus) %>%
  dplyr::left_join(family_to_group %>% dplyr::rename(taxon_from_family = taxon_group),
                   by = "family_lc") %>%
  dplyr::mutate(
    taxon_group = dplyr::coalesce(taxon_group, taxon_from_family),
    taxon_group = dplyr::case_when(is.na(taxon_group) | !nzchar(taxon_group) ~ "Other", TRUE ~ taxon_group),
    taxon_group = factor(taxon_group, levels = c("Fishes","Shrimps/Crabs","Snails","Other")),
    species_label = paste0("italic(", gsub("_","~", species_original), ")"),
    Effect       = dplyr::case_when(is.na(r) ~ NA_character_, r > 0 ~ "Positive", TRUE ~ "Negative")
  ) %>%
  dplyr::group_by(transform) %>%
  dplyr::mutate(p_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()

# --- Add a per-species r with PC1 (sqrt) for ordering & export ----------------
r_pc1_tbl <- corr_all %>%
  dplyr::filter(transform == "sqrt", metric == pc1_coral_col) %>%
  dplyr::select(species_clean, r_pc1_species = r)

corr_all <- corr_all %>%
  dplyr::left_join(r_pc1_tbl, by = "species_clean")

# Save full table
out_corr_csv <- file.path(TAB_DIR, "species_performance_correlations_raw_and_sqrt_w_groups.csv")
readr::write_csv(corr_all, out_corr_csv)
cli::cli_alert_success("Wrote correlation table (incl. taxon_group & r_pc1_species) → {out_corr_csv}")

# --- Heatmap (√ abundance only) ----------------------------------------------
TOP_SHOW <- 60

# Pretty labels for metrics with line breaks for long labels
pretty_metric <- function(x) dplyr::case_when(
  x == pc1_coral_col       ~ "PC1[coral]",
  x == "protein_mg_cm2"    ~ "atop(Protein, (mg~cm^{-2}))",
  x == "carb_mg_cm2"       ~ "atop(Carbohydrate, (mg~cm^{-2}))",
  x == "zoox_cells_cm2"    ~ "atop(Zooxanthellae, (cells~cm^{-2}))",
  x == "afdw_mg_cm2"       ~ "atop(Tissue~biomass, (AFDW~mg~cm^{-2}))",
  x == "growth_vol_b"      ~ "Growth",
  TRUE                     ~ x
)

heat_df <- corr_all %>%
  dplyr::filter(transform == "sqrt") %>%
  dplyr::mutate(metric_label = pretty_metric(metric))


# Make PC1 the leftmost column in the heatmap
metric_levels <- c("PC1[coral]", setdiff(unique(heat_df$metric_label), "PC1[coral]"))
heat_df <- heat_df %>% dplyr::mutate(metric_label = factor(metric_label, levels = metric_levels))

# Select top species (by max |r| across metrics) to keep figure readable
top_species <- heat_df %>%
  dplyr::group_by(species_original) %>%
  dplyr::summarise(max_abs_r = max(abs(r), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(max_abs_r)) %>%
  dplyr::slice_head(n = TOP_SHOW) %>%
  dplyr::pull(species_original)

heat_df <- heat_df %>%
  dplyr::filter(species_original %in% top_species)

# Order rows WITHIN EACH TAXON GROUP by r with PC1 (descending), fallback = max |r|
order_df <- heat_df %>%
  dplyr::group_by(taxon_group, species_original) %>%
  dplyr::summarise(
    r_pc1_ord = dplyr::first(na.omit(r_pc1_species)),
    fallback  = max(abs(r), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(order_key = dplyr::coalesce(r_pc1_ord, fallback)) %>%
  dplyr::arrange(taxon_group, dplyr::desc(order_key)) %>%
  dplyr::group_by(taxon_group) %>%
  dplyr::mutate(species_order = dplyr::row_number()) %>%
  dplyr::ungroup()

heat_df <- heat_df %>%
  dplyr::left_join(order_df %>% dplyr::select(taxon_group, species_original, species_order),
                   by = c("taxon_group","species_original")) %>%
  dplyr::mutate(species_display = stats::reorder(species_original, species_order))

#THIS IS THE A PUBLICATION FIGURE


# Plot - Publication quality heatmap
p_heat <- ggplot(heat_df, aes(x = metric_label, y = species_display, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.finite(r), sprintf("%.02f", r), "")),
            size = 3.2, color = "black", fontface = "plain") +
  scale_x_discrete(position = "top", labels = function(x) parse(text = x), expand = c(0, 0)) +
  scale_y_discrete(labels = function(y) parse(text = paste0("italic('", gsub("_"," ", y), "')")),
                   expand = c(0, 0)) +
  scale_fill_gradient2(
    name = "Correlation (r)",
    limits = c(-1, 1),
    low = "#B2182B",      # Dark red for negative correlations
    mid = "white",         # White for zero
    high = "#2166AC",      # Dark blue for positive correlations
    midpoint = 0,
    breaks = seq(-1, 1, 0.5),
    guide = guide_colorbar(
      direction = "vertical",
      barwidth = unit(1.2, "cm"),
      barheight = unit(12, "cm"),
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      ticks.colour = "black",
      ticks.linewidth = 0.5
    )
  ) +
  labs(x = NULL, y = NULL) +
  facet_grid(rows = vars(taxon_group), scales = "free_y", space = "free_y") +
  theme_pub(base_size = 14) +
  theme(
    # X-axis (top) - horizontal text with line breaks
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10, face = "plain", lineheight = 0.9),
    axis.ticks.x.top = element_line(linewidth = 0.5),
    axis.ticks.length.x.top = unit(0.2, "cm"),

    # Y-axis (species names)
    axis.text.y = element_text(size = 10, face = "italic", hjust = 1),
    axis.ticks.y = element_line(linewidth = 0.5),
    axis.ticks.length.y = unit(0.15, "cm"),

    # Facet strips (taxon groups) - on the right, VERTICAL text
    strip.text.y.right = element_text(size = 12, face = "bold", angle = 270),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),

    # Panel styling
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.y = unit(0.3, "cm"),
    panel.spacing.x = unit(0, "lines"),

    # Legend
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.margin = margin(l = 8, r = 5),
    legend.box.spacing = unit(0.4, "cm"),

    # Overall plot margins
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )

show_and_save(
  p_heat,
  file.path(FIG_DIR, "species_performance_corr_heatmap_sqrt_byGroup_orderedByPC1"),
  width = 10, height = 11
)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(OUT_DIR, "figures", "publication-figures", "species_performance_corr_heatmap_sqrt_byGroup_orderedByPC1.pdf"),
  p_heat, width = 10, height = 11, dpi = 600, bg = "white"
)
ggsave(
  file.path(OUT_DIR, "figures", "publication-figures", "species_performance_corr_heatmap_sqrt_byGroup_orderedByPC1.png"),
  p_heat, width = 10, height = 11, dpi = 600, bg = "white"
)

cli::cli_alert_success("§5.6 done: heatmap saved; rows ordered by r(species, PC1) within Fish / Shrimps+Crabs / Snails / Other.")








# =============================================================================
# 6) Species-level LMMs: top 20 RAW PC1 loadings → cond_PC1
# =============================================================================

cli::cli_h2("6) Species-level regressions (top-20 by |RAW PC1 loading|)")

# 6.1 Identify top-20 species by absolute RAW PC1 loading ----------------------
raw_loadings <- tibble::tibble(
  species = rownames(res_raw$pca_comm$rotation),
  loading = res_raw$pca_comm$rotation[, 1]
) |>
  dplyr::arrange(dplyr::desc(abs(loading))) |>
  dplyr::mutate(sign = dplyr::if_else(loading >= 0, "+", "-"))

top20 <- raw_loadings |> dplyr::slice_head(n = 20)
cli::cli_alert_success("Selected top {nrow(top20)} species by |RAW PC1 loading|.")

# 6.2 Build analysis frame: sqrt-abundance × cond_PC1 × reef -------------------
# long sqrt-abundance for all species
sqrt_long <- comm_sqrt_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("coral_id") |>
  tidyr::pivot_longer(-coral_id, names_to = "species", values_to = "sqrt_abund")

# keep only top-20 species
sqrt_long_top <- sqrt_long |> dplyr::semi_join(top20, by = "species")

# cond_PC1 (from Section 2 helper) + reef_id
cond_df <- res_raw$scores_cond |> dplyr::select(coral_id, cond_PC1)
reef_df <- growth            |> dplyr::select(coral_id, reef_id = reef)

anal_top <- sqrt_long_top |>
  dplyr::inner_join(cond_df, by = "coral_id") |>
  dplyr::left_join(reef_df, by = "coral_id") |>
  dplyr::filter(!is.na(reef_id))

# 6.3 Fit per-species LMM: cond_PC1 ~ sqrt_abund + (1|reef_id) ----------------
fit_one <- function(df) {
  tryCatch(
    lmerTest::lmer(cond_PC1 ~ sqrt_abund + (1 | reef_id),
                   data = df,
                   control = lme4::lmerControl(check.conv.singular = "ignore")),
    error = function(e) NULL
  )
}

lmm_top20 <- anal_top |>
  dplyr::group_by(species) |>
  dplyr::group_modify(~ {
    mod <- fit_one(.x)
    if (is.null(mod)) {
      out <- tibble::tibble(term = "sqrt_abund", estimate = NA_real_,
                            std.error = NA_real_, df = NA_real_,
                            statistic = NA_real_, p.value = NA_real_,
                            conf.low = NA_real_, conf.high = NA_real_)
    } else {
      out <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) |>
        dplyr::filter(term == "sqrt_abund") |>
        dplyr::select(term, estimate, std.error, df, statistic, p.value, conf.low, conf.high)
    }
    out
  }) |>
  dplyr::ungroup() |>
  dplyr::left_join(top20, by = "species") |>
  dplyr::mutate(
    p_adj_hochberg = p.adjust(p.value, method = "hochberg")) |>
  dplyr::arrange(dplyr::desc(loading))  # display strong positives first

# 6.4 Save table ---------------------------------------------------------------
out_csv <- file.path(TAB_DIR, "LMM_species_top20_RAWsqrt_vs_condPC1.csv")
readr::write_csv(lmm_top20, out_csv)
cli::cli_alert_success("Wrote species-level LMM table → {out_csv}")

# 6.5 Coefficient plot (estimate ± 95% CI), ordered by RAW loading ------------
# Explicit factor order so positively loading species appear at the TOP after coord_flip
order_levels <- lmm_top20 |>
  dplyr::arrange(loading) |>
  dplyr::mutate(lab = paste0("italic(", gsub(" ", "~", species), ")")) |>
  dplyr::pull(lab)

coef_plot <- lmm_top20 |>
  dplyr::mutate(
    species_label = paste0("italic(", gsub(" ", "~", species), ")"),
    species_label = factor(species_label, levels = order_levels, ordered = TRUE)
  ) |>
  ggplot2::ggplot(aes(x = species_label, y = estimate)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "gray60") +
  ggplot2::geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  ggplot2::geom_point() +
  ggplot2::coord_flip() +
  ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
  ggplot2::labs(
    title = "Top-20 species: effect of sqrt-abundance on Coral condition PC1",
    x = NULL,
    y = expression(widehat(beta)~" for "~sqrt(abundance))
  ) +
  theme_pub(10)

show_and_save(coef_plot,
              file.path(FIG_DIR, "LMM_species_top20_coefplot"),
              width = 7, height = 7)

# 6.6 (Optional) Pretty GT table ----------------------------------------------
lmm_top20 %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate(
    loading = signif(loading, 3),
    estimate = signif(estimate, 3),
    std.error = signif(std.error, 3),
    conf.low = signif(conf.low, 3),
    conf.high = signif(conf.high, 3),
    p.value = signif(p.value, 3),
    p_adj_hochberg = signif(p_adj_hochberg, 3)
  ) %>%
  gt::gt() %>%
  gt::tab_header(
    title = gt::md("**Top-20 species (by |RAW PC1 loading|): LMM on Coral condition PC1**"),
    subtitle = "Model: cond_PC1 ~ sqrt(abundance) + (1 | reef)"
  ) %>%
  gt::cols_label(
    species   = "Species",
    loading   = "RAW PC1 loading",
    estimate  = "β (sqrt abundance)",
    std.error = "SE",
    conf.low  = "CI low",
    conf.high = "CI high",
    p.value   = "p",
    p_adj_hochberg = "p (Hochberg)"
  ) %>%
  print()

# =============================================================================
# 6.6) Per-species regression panels (points + lmer fixed-effect line)
# =============================================================================

cli::cli_h2("6.6) Plotting species-level regressions with lmer fits")

suppressPackageStartupMessages({
  library(broom.mixed)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Sanity: make sure the predictor exists (rename here if your term differs)
stopifnot("sqrt_abund" %in% names(anal_top))

# Keep only top-20 species and complete cases
anal_top20 <- anal_top %>%
  semi_join(top20, by = "species") %>%
  filter(!is.na(sqrt_abund), !is.na(cond_PC1))

# Safe wrapper around your fit function
safe_fit <- possibly(~ fit_one(.x), otherwise = NULL, quiet = TRUE)

# One model per species (list-column)
mods_df <- anal_top20 %>%
  nest_by(species) %>%
  mutate(model = list(safe_fit(data))) %>%
  ungroup() %>%
  filter(!vapply(model, is.null, logical(1)))

# Unwrap to the actual lmer model if fit_one() returns a list
mods_df <- mods_df %>%
  mutate(
    model_obj = map(model, function(m) {
      if (inherits(m, "list")) {
        if (!is.null(m$model)) m <- m$model
        else if (!is.null(m$fit)) m <- m$fit
        else if (length(m) >= 1)  m <- m[[1]]
      }
      m
    })
  ) %>%
  filter(map_lgl(model_obj, ~ inherits(.x, c("lmerMod", "lmerModLmerTest"))))

# Extract fixed effects (intercept + slope for sqrt_abund)
coef_df <- mods_df %>%
  mutate(tidy = map(model_obj, ~ broom.mixed::tidy(.x, effects = "fixed"))) %>%
  unnest(tidy) %>%
  dplyr::select(species, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept = `(Intercept)`) %>%
  # if your predictor is named differently, change it here:
  rename(slope = `sqrt_abund`) %>%
  left_join(top20, by = "species")  # brings 'loading' for ordering

# Build prediction grid for fixed-effect line per species
pred_grid <- anal_top20 %>%
  group_by(species) %>%
  summarise(
    x_min = min(sqrt_abund, na.rm = TRUE),
    x_max = max(sqrt_abund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(sqrt_abund = list(seq(x_min, x_max, length.out = 100))) %>%
  unnest(sqrt_abund) %>%
  dplyr::select(species, sqrt_abund) %>%
  left_join(coef_df %>% dplyr::select(species, intercept, slope), by = "species") %>%
  mutate(pred_fixed = intercept + slope * sqrt_abund)

# Facet labels ordered by descending RAW PC1 loading (if available)
facet_levels <- coef_df %>%
  arrange(desc(loading)) %>%
  pull(species)

label_df <- tibble::tibble(
  species = facet_levels,
  species_label = paste0("italic(", gsub(" ", "~", facet_levels), ")")
)

plot_dat <- anal_top20 %>%
  left_join(label_df, by = "species") %>%
  mutate(species_label = factor(species_label,
                                levels = label_df$species_label,
                                ordered = TRUE))

pred_grid_lab <- pred_grid %>%
  left_join(label_df, by = "species") %>%
  mutate(species_label = factor(species_label,
                                levels = label_df$species_label,
                                ordered = TRUE))

# Plot: points + fixed-effect line
p_species_lmm <- ggplot(plot_dat, aes(x = sqrt_abund, y = cond_PC1)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_line(data = pred_grid_lab, aes(y = pred_fixed), linewidth = 0.8) +
  facet_wrap(~ species_label, ncol = 4, labeller = ggplot2::label_parsed,
             scales = "free_x") +
  labs(
    x = expression(sqrt(Abundance)),
    y = expression(Coral~condition~(PC[1][coral])),
    title = "Top-20 species (by |RAW PC1 loading|): cond_PC1 ~ sqrt(abundance) + (1|reef)"
  ) +
  theme_pub()  # your project theme

# Save
stub <- file.path(FIG_DIR, "species_top20_LMM_panels")
show_and_save(p_species_lmm, stub, width = 12, height = 14)
cli::cli_alert_success("Panel figure written: {stub}.png and {stub}.pdf")


# =============================================================================
# 6.7) Faceted species fits (2 cols × 3 rows): raw abundance on sqrt axis
#       – sqrt-spaced ticks, arithmetic labels, ≥3 ticks/facet, small left gap
# =============================================================================

cli::cli_h2("6.7) Faceted species fits (raw abundance, sqrt axis)")

plot_species <- c(
  "Caracanthus maculatus",
  "Luniella pugil",
  "Alpheus diadema",
  "Calcinus latens",
  "Harpiliopsis spinigera",
  "Dascyllus flavicaudus"
)
plot_species <- intersect(plot_species, unique(anal_top$species))

if (length(plot_species) == 0L) {
  cli::cli_alert_danger("None of the requested species are present. Skipping 6.7.")
} else {
  
  # --- Data used in the model (has sqrt_abund); add raw for plotting ----------
  dat_sub <- anal_top %>%
    dplyr::filter(species %in% plot_species) %>%
    dplyr::filter(!is.na(sqrt_abund), !is.na(cond_PC1)) %>%
    dplyr::mutate(raw_abund = sqrt_abund^2)
  
  # --- Per-species fits (reuse fit_one mechanics from §6.6) -------------------
  safe_fit <- purrr::possibly(~ fit_one(.x), otherwise = NULL, quiet = TRUE)
  
  mods_df <- dat_sub %>%
    dplyr::nest_by(species) %>%
    dplyr::mutate(model = list(safe_fit(data))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!vapply(model, is.null, logical(1))) %>%
    dplyr::mutate(
      model_obj = purrr::map(model, function(m) {
        if (inherits(m, "list")) {
          if (!is.null(m$model)) m <- m$model
          else if (!is.null(m$fit)) m <- m$fit
          else if (length(m) >= 1)  m <- m[[1]]
        }
        m
      })
    ) %>%
    dplyr::filter(purrr::map_lgl(model_obj, ~ inherits(.x, c("lmerMod","lmerModLmerTest"))))
  
  # Fixed-effect intercept & slope for sqrt_abund
  coef_df <- mods_df %>%
    dplyr::mutate(tidy = purrr::map(
      model_obj,
      ~ tryCatch(broom.mixed::tidy(.x, effects = "fixed"), error = function(e) NULL)
    )) %>%
    tidyr::unnest(tidy) %>%
    dplyr::select(species, term, estimate) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
    dplyr::rename(intercept = `(Intercept)`, slope = `sqrt_abund`)
  
  # Prediction grid in RAW units; predict via sqrt(raw_abund)
  rng_df <- dat_sub %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      x_min = max(0, min(raw_abund, na.rm = TRUE)),
      x_max = max(raw_abund, na.rm = TRUE),
      .groups = "drop"
    )
  
  pred_lines <- rng_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(raw_abund = list(seq(x_min, x_max, length.out = 100))) %>%
    tidyr::unnest(raw_abund) %>%
    dplyr::left_join(coef_df, by = "species") %>%
    dplyr::mutate(pred_fixed = intercept + slope * sqrt(raw_abund))
  
  # Parsed italic facet labels
  label_df <- tibble::tibble(
    species = plot_species,
    species_label = paste0("italic(", gsub(" ", "~", plot_species), ")")
  )
  
  dat_sub <- dat_sub %>%
    dplyr::left_join(label_df, by = "species") %>%
    dplyr::mutate(species_label = factor(species_label, levels = label_df$species_label))
  
  pred_lines <- pred_lines %>%
    dplyr::left_join(label_df, by = "species") %>%
    dplyr::mutate(species_label = factor(species_label, levels = label_df$species_label))
  
  # Palette per species using taxonomic colors (matching Figure 5 per Craig's suggestion)
  # Fishes = blue, Shrimps/Crabs = magenta/pink, Snails = orange
  taxon_cols <- c(
    "Fishes"        = "#2E86AB",  # blue
    "Shrimps/Crabs" = "#A23B72",  # magenta/pink
    "Snails"        = "#F18F01"   # orange
  )

  # Assign taxonomic group to each species
  species_taxon <- c(
    "Caracanthus maculatus"   = "Fishes",
    "Luniella pugil"          = "Shrimps/Crabs",
    "Alpheus diadema"         = "Shrimps/Crabs",
    "Calcinus latens"         = "Shrimps/Crabs",
    "Harpiliopsis spinigera"  = "Shrimps/Crabs",
    "Dascyllus flavicaudus"   = "Fishes"
  )

  # Map each species to its taxonomic color
  base_cols <- sapply(plot_species, function(sp) taxon_cols[species_taxon[sp]])
  pal_cols  <- stats::setNames(base_cols, plot_species)
  
  # ---- Breaks helper: rounder arithmetic labels, ≥3 ticks, sqrt spacing ------
  sqrt_breaks_smart <- function(lims) {
    lo <- max(0, lims[1]); hi <- max(lo, lims[2])
    if (!is.finite(hi) || hi <= 0) return(c(0, 1, 2))
    
    # Work in sqrt-space; get "pretty" positions there, then square back.
    k_lo <- sqrt(lo); k_hi <- sqrt(hi)
    k_seq <- scales::pretty_breaks(n = 4)(c(k_lo, k_hi))
    b <- (k_seq^2)
    
    # Keep within range and ensure at least 3 ticks
    b <- sort(unique(b[b >= lo & b <= hi]))
    if (length(b) < 3) {
      k_seq <- seq(k_lo, k_hi, length.out = 3)
      b <- sort(unique(c(b, k_seq^2)))
    }
    
    # Round labels to nice values (0 decimals for 0–20, else 1)
    if (max(b, na.rm = TRUE) <= 20) round(b, 0) else round(b, 1)
  }
  
  # Label helper: print integers cleanly; others with 1 decimal
  label_arith <- function(b) {
    sapply(b, function(x) {
      if (is.na(x)) return(NA_character_)
      if (abs(x - round(x)) < 1e-6) as.character(round(x)) else format(round(x, 1), nsmall = 1)
    })
  }
  
  # ---- Create individual plots with confidence intervals ----------------------
  # Generate predictions with confidence intervals using ggpredict
  pred_ci_list <- mods_df %>%
    dplyr::mutate(
      pred_ci = purrr::map(model_obj, function(m) {
        tryCatch({
          ggeffects::ggpredict(m, terms = "sqrt_abund [all]", type = "fixed")
        }, error = function(e) NULL)
      })
    ) %>%
    dplyr::select(species, pred_ci) %>%
    dplyr::filter(!purrr::map_lgl(pred_ci, is.null))

  # Create individual plots for each species
  plot_list <- list()

  for (i in seq_along(plot_species)) {
    sp <- plot_species[i]
    sp_data <- dat_sub %>% dplyr::filter(species == sp)

    # Get predictions with CI
    pred_data <- pred_ci_list %>%
      dplyr::filter(species == sp) %>%
      dplyr::pull(pred_ci) %>%
      .[[1]]

    if (is.null(pred_data)) next

    # Convert back to raw scale
    pred_df <- data.frame(
      raw_abund = pred_data$x^2,
      fit = pred_data$predicted,
      lower = pred_data$conf.low,
      upper = pred_data$conf.high
    )

    # Create plot
    p <- ggplot2::ggplot() +
      # Confidence ribbon
      ggplot2::geom_ribbon(
        data = pred_df,
        ggplot2::aes(x = raw_abund, ymin = lower, ymax = upper),
        fill = base_cols[i], alpha = 0.2
      ) +
      # Regression line
      ggplot2::geom_line(
        data = pred_df,
        ggplot2::aes(x = raw_abund, y = fit),
        color = base_cols[i], linewidth = 1.1
      ) +
      # Raw data points
      ggplot2::geom_point(
        data = sp_data,
        ggplot2::aes(x = raw_abund, y = cond_PC1),
        color = base_cols[i], alpha = 0.65, size = 1.8
      ) +
      # Add species name as subtitle
      ggplot2::labs(
        subtitle = bquote(italic(.(gsub("_", " ", sp)))),
        x = NULL,
        y = NULL
      ) +
      ggplot2::scale_x_sqrt(
        breaks = sqrt_breaks_smart,
        labels = label_arith,
        limits = c(0, NA),
        expand = ggplot2::expansion(mult = c(0, 0.03), add = c(0.3, 0))
      ) +
      theme_pub() +
      ggplot2::theme(
        axis.ticks.length.x = grid::unit(6, "pt"),
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0, face = "italic"),
        plot.margin = margin(5, 5, 5, 5)
      )

    plot_list[[i]] <- p
  }

  # Combine plots with cowplot and add letter labels
  p_combined <- cowplot::plot_grid(
    plotlist = plot_list,
    ncol = 2, nrow = 3,
    align = "hv",
    axis = "tblr",
    labels = LETTERS[1:6],
    label_size = 12,
    label_fontface = "bold",
    label_x = 0.95,
    label_y = 0.96,
    hjust = 1,
    vjust = 1
  )

  # Add shared axis labels
  y_label <- grid::textGrob(
    expression(PC1[coral]),
    gp = grid::gpar(fontsize = 13, fontface = "bold"),
    rot = 90
  )

  x_label <- grid::textGrob(
    "Abundance (no. / coral)",
    gp = grid::gpar(fontsize = 13)
  )

  # Create legend for taxonomic groups
  legend_data <- data.frame(
    taxon = factor(c("Fishes", "Shrimps/Crabs"), levels = c("Fishes", "Shrimps/Crabs")),
    x = 1:2, y = 1:2
  )

  p_legend <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y, fill = taxon)) +
    ggplot2::geom_point(shape = 21, size = 4, color = "black", stroke = 0.3) +
    ggplot2::scale_fill_manual(
      values = c("Fishes" = "#2E86AB", "Shrimps/Crabs" = "#A23B72"),
      name = NULL
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      direction = "horizontal",
      override.aes = list(size = 4)
    )) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      legend.key.size = grid::unit(0.8, "lines"),
      legend.spacing.x = grid::unit(0.5, "cm")
    )

  # Extract legend
  legend_grob <- cowplot::get_legend(p_legend)

  # Combine with axis labels and legend (legend on top)
  p_faceted <- cowplot::plot_grid(
    y_label,
    cowplot::plot_grid(
      legend_grob,
      p_combined,
      x_label,
      ncol = 1,
      rel_heights = c(0.06, 1, 0.05)
    ),
    nrow = 1,
    rel_widths = c(0.05, 1)
  )
  
  out_stub <- file.path(FIG_DIR, "species_faceted_LMM_lines_rawX_sqrtAxis")
  show_and_save(p_faceted, out_stub, width = 5, height = 6)
  cli::cli_alert_success(glue::glue("Faceted figure written: {out_stub}.png and {out_stub}.pdf"))

  # Also save to publication-figures folder (PDF and PNG)
  ggsave(
    file.path(OUT_DIR, "figures", "publication-figures", "species_faceted_LMM_lines_rawX_sqrtAxis.pdf"),
    p_faceted, width = 5, height = 6, dpi = 600, bg = "white"
  )
  ggsave(
    file.path(OUT_DIR, "figures", "publication-figures", "species_faceted_LMM_lines_rawX_sqrtAxis.png"),
    p_faceted, width = 5, height = 6, dpi = 600, bg = "white"
  )
}

#print out this figure



# =============================================================================
# 6.8) Faceted species fits (3 columns × 2 rows), PPT wide with title
# =============================================================================
cli::cli_h2("6.8) Faceted species fits (3x2, wide slide with title)")

if (length(plot_species) == 0L) {
  cli::cli_alert_danger("No species available for 6.8; skipping.")
} else {
  
  # If you didn't run 6.7 above, uncomment this block to (re)build the needed objects:
  # dat_sub <- anal_top %>%
  #   dplyr::filter(species %in% plot_species) %>%
  #   dplyr::filter(!is.na(sqrt_abund), !is.na(cond_PC1)) %>%
  #   dplyr::mutate(raw_abund = sqrt_abund^2)
  # safe_fit <- purrr::possibly(~ fit_one(.x), otherwise = NULL, quiet = TRUE)
  # mods_df <- dat_sub %>%
  #   dplyr::nest_by(species) %>% dplyr::mutate(model = list(safe_fit(data))) %>%
  #   dplyr::ungroup() %>% dplyr::filter(!vapply(model, is.null, logical(1))) %>%
  #   dplyr::mutate(model_obj = purrr::map(model, ~ if (inherits(.x, "list")) {
  #     if (!is.null(.x$model)) .x$model else if (!is.null(.x$fit)) .x$fit else .x[[1]]
  #   } else .x)) %>%
  #   dplyr::filter(purrr::map_lgl(model_obj, ~ inherits(.x, c("lmerMod","lmerModLmerTest"))))
  # coef_df <- mods_df %>%
  #   dplyr::mutate(tidy = purrr::map(model_obj, ~ broom.mixed::tidy(.x, effects = "fixed"))) %>%
  #   tidyr::unnest(tidy) %>% dplyr::select(species, term, estimate) %>%
  #   tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
  #   dplyr::rename(intercept = `(Intercept)`, slope = `sqrt_abund`)
  # rng_df <- dat_sub %>% dplyr::group_by(species) %>%
  #   dplyr::summarise(x_min = pmax(0, min(raw_abund, na.rm = TRUE)),
  #                    x_max = max(raw_abund, na.rm = TRUE), .groups = "drop")
  # pred_lines <- rng_df %>% dplyr::rowwise() %>%
  #   dplyr::mutate(raw_abund = list(seq(x_min, x_max, length.out = 100))) %>%
  #   tidyr::unnest(raw_abund) %>% dplyr::left_join(coef_df, by = "species") %>%
  #   dplyr::mutate(pred_fixed = intercept + slope * sqrt(raw_abund))
  # label_df <- tibble::tibble(
  #   species = plot_species,
  #   species_label = paste0("italic(", gsub(" ", "~", plot_species), ")")
  # )
  # dat_sub <- dat_sub %>% dplyr::left_join(label_df, by = "species") %>%
  #   dplyr::mutate(species_label = factor(species_label, levels = label_df$species_label))
  # pred_lines <- pred_lines %>% dplyr::left_join(label_df, by = "species") %>%
  #   dplyr::mutate(species_label = factor(species_label, levels = label_df$species_label))
  
  # Common title for the slide
  fig_title <- "Species-level links: Abundance (sqrt axis) vs Coral condition (PC1)"
  
  # 16:9 slide export size (PowerPoint default): 13.33 in × 7.5 in
  slide_w <- 13.33
  slide_h <- 7.5
  
  # Breaks for the sqrt x-axis (computed globally)
  max_x <- max(dat_sub$raw_abund, na.rm = TRUE)
  sqrt_breaks <- (0:floor(sqrt(max_x)))^2
  
  p_faceted_wide <- ggplot2::ggplot(
    dat_sub,
    ggplot2::aes(x = raw_abund, y = cond_PC1)
  ) +
    ggplot2::geom_point(ggplot2::aes(colour = species), alpha = 0.65, size = 1.8) +
    ggplot2::geom_line(
      data = pred_lines,
      ggplot2::aes(x = raw_abund, y = pred_fixed, colour = species),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(values = pal_cols, guide = "none") +
    # If your ggplot2 lacks scale_x_sqrt(), use:
    # ggplot2::scale_x_continuous(trans = "sqrt", breaks = sqrt_breaks, labels = sqrt_breaks)
    ggplot2::scale_x_sqrt(breaks = sqrt_breaks, labels = sqrt_breaks) +
    ggplot2::facet_wrap(
      ~ species_label,
      ncol = 3,                                # <-- 3 columns × 2 rows
      labeller = ggplot2::label_parsed,
      scales = "free_x"
    ) +
    ggplot2::labs(
      title = fig_title,
      x = "Abundance (square-root axis)",
      y = expression(Coral~condition~(PC[1][coral]))
    ) +
    theme_pub(base_size = 13) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 6)),
      strip.text      = ggplot2::element_text(face = "bold", size = 11),
      panel.spacing   = grid::unit(6, "pt"),
      plot.margin     = ggplot2::margin(t = 6, r = 8, b = 6, l = 8)
    )
  
  out_stub_wide <- file.path(FIG_DIR, "species_faceted_LMM_lines_rawX_sqrtAxis_wide_3x2")
  show_and_save(p_faceted_wide, out_stub_wide, width = slide_w, height = slide_h, dpi = 600)
  cli::cli_alert_success(glue::glue("Wrote PPT-wide figure → {out_stub_wide}.png / .pdf"))
}



# =============================================================================
# 6.9) Build master_df for Section 7 (traits + species + meta, by coral_id)
# =============================================================================
cli::cli_h2("6.9) Building master_df for Section 7")

# Helper: make clean/safe species column names that match code-friendly references
clean_species_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)  # replace spaces/punct with _
  x <- gsub("^_|_$", "", x)           # trim leading/trailing _
  x
}

# 1) Start from community abundance matrix (already filtered to core spp)
#    -> counts by coral_id × species (RAW)
comm_df_raw <- comm_raw_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("coral_id")

# clean species colnames and keep a mapping for reference
orig_names <- colnames(comm_raw_mat)
clean_names <- clean_species_name(orig_names)
colnames(comm_df_raw) <- c("coral_id", clean_names)

species_name_map <- tibble::tibble(
  species_original = orig_names,
  species_clean    = clean_names
)

# 2) Coral traits (use what you actually have in `growth`/`physio`)
#    From your glimpses: percent_alive, delta_surface_area, delta_volume,
#    delta_max_height, delta_interstitial, growth_vol_b + physio metrics.
req_growth_cols <- c("coral_id", "growth_vol_b")  # <- drop percent_alive here
req_physio_cols <- c("coral_id", "percent_alive", "protein_mg_cm2", "carb_mg_cm2",
                     "zoox_cells_cm2", "afdw_mg_cm2")

# Build traits_df (unchanged except growth now only contributes coral_id + growth_vol_b)
traits_df <- growth %>%
  dplyr::select(dplyr::all_of(req_growth_cols)) %>%
  dplyr::inner_join(
    physio %>% dplyr::select(dplyr::all_of(req_physio_cols), reef, treatment),
    by = "coral_id"
  ) %>%
  dplyr::rename(reef_id = reef) %>%
  dplyr::mutate(treatment = factor(as.character(treatment), levels = c("1","3","6")))

# 3) Merge traits with community counts
master_df <- traits_df |>
  dplyr::inner_join(comm_df_raw, by = "coral_id")

# 4) Create standardized versions of continuous predictors (z-scores)
z_cols <- intersect(
  c("percent_alive","delta_surface_area","delta_volume","delta_max_height",
    "delta_interstitial","growth_vol_b",
    "protein_mg_cm2","carb_mg_cm2","zoox_cells_cm2","afdw_mg_cm2"),
  names(master_df)
)

for (vv in z_cols) {
  master_df[[paste0(vv, "_z")]] <- as.numeric(scale(master_df[[vv]]))
}

# 5) Quick sanity + message on species columns
present_spp <- setdiff(names(master_df), c("coral_id","reef_id","treatment", z_cols,
                                           paste0(z_cols, "_z")))
present_spp <- intersect(present_spp, clean_names)

cli::cli_alert_success("master_df ready: {nrow(master_df)} corals × {length(present_spp)} focal species columns.")
cli::cli_alert_info("Example species columns: {paste(utils::head(present_spp, 6), collapse = ', ')}")

# OPTIONAL: write to disk for reuse
# readr::write_rds(master_df, file.path(OUT_DIR, "master_df_for_section7.rds"))


# =============================================================================
# 7. Network Analysis (reef-adjusted): All Species + Coral Traits (|r| ≥ 0.30)
# =============================================================================

cli::cli_h1("Section 7: Network analysis (reef-adjusted residual correlations)")

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(lme4)

# --- Key change A: enforce ≥80% alive and then drop 'percent_alive' from traits
# ALIVE_THRESH is defined in utils.R (0.80) - sourced at top of script

# If master_df has percent_alive, use it only to FILTER, then remove it from analysis
if ("percent_alive" %in% names(master_df)) {
  keep_ids_alive <- master_df %>%
    filter(!is.na(percent_alive), percent_alive >= ALIVE_THRESH) %>%
    pull(coral_id)
  n_before <- nrow(master_df)
  master_df <- master_df %>%
    filter(coral_id %in% keep_ids_alive)
  cli::cli_alert_success(
    "Alive filter applied (≥{ALIVE_THRESH*100}%): kept {nrow(master_df)}/{n_before} corals."
  )
} else if (exists("keep_ids_alive")) {
  # Fall back to any previously computed keep_ids_alive
  n_before <- nrow(master_df)
  master_df <- master_df %>% filter(coral_id %in% keep_ids_alive)
  cli::cli_alert_success(
    "Alive filter applied from keep_ids_alive: kept {nrow(master_df)}/{n_before} corals."
  )
} else {
  cli::cli_alert_warning(
    "No percent_alive column in master_df and no keep_ids_alive found; proceeding without extra filtering."
  )
}

# --- Key change B: REMOVE 'percent_alive' from the analysis variables ----------
physio_vars <- c(
  "growth_vol_b","protein_mg_cm2","carb_mg_cm2",
  "zoox_cells_cm2","afdw_mg_cm2"     # <-- percent_alive removed
)
meta_cols   <- c("coral_id","reef_id","treatment")

# Treat everything else numeric as "species" (drop any *_z leftovers)
species_vars <- names(master_df)
species_vars <- setdiff(
  species_vars,
  c(meta_cols, physio_vars, grep("_z$", species_vars, value = TRUE), "percent_alive")
)

# Build SEM input frame and drop incomplete cases
sem_df <- master_df %>%
  dplyr::select(dplyr::all_of(c(meta_cols, physio_vars, species_vars))) %>%
  tidyr::drop_na()

# 7.2 Remove reef-level signal with a random-intercept -------------------------
# Residualize EACH variable against reef_id: y ~ 1 + (1|reef_id)
reef_residualize <- function(y, reef) {
  try_fit <- tryCatch(lmer(y ~ 1 + (1 | reef), REML = TRUE),
                      error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(try_fit)) {
    res <- resid(try_fit)
  } else {
    # fallback: within-reef demeaning
    res <- y - ave(y, reef, FUN = function(x) mean(x, na.rm = TRUE))
  }
  as.numeric(res)
}

vars_use <- c(physio_vars, species_vars)

resid_df <- sem_df %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(vars_use),
                              ~ reef_residualize(.x, sem_df$reef_id),
                              .names = "{.col}_resid"))

# Standardize residuals (z-score) before correlation
sem_use <- resid_df %>%
  dplyr::mutate(dplyr::across(dplyr::ends_with("_resid"),
                              ~ as.numeric(scale(.x)),
                              .names = "{.col}_z")) %>%
  dplyr::select(dplyr::all_of(meta_cols), dplyr::ends_with("_resid_z")) %>%
  dplyr::rename_with(~ sub("_resid_z$", "", .x), dplyr::ends_with("_resid_z"))

# 7.3 Correlation matrix from reef-adjusted residuals --------------------------
num_mat <- sem_use %>%
  dplyr::select(-dplyr::all_of(meta_cols)) %>%
  as.matrix()

cor_mat <- cor(num_mat, use = "pairwise.complete.obs")

mat_upper <- cor_mat
mat_upper[lower.tri(mat_upper, diag = TRUE)] <- NA

cor_df <- as.data.frame(as.table(mat_upper), stringsAsFactors = FALSE) %>%
  dplyr::rename(from = Var1, to = Var2, weight = Freq) %>%
  dplyr::filter(!is.na(weight)) %>%
  dplyr::mutate(
    abs_weight = abs(weight),
    Effect = ifelse(weight > 0, "Positive", "Negative")
  )

# 7.4 Apply r cutoff, build helper tables -------------------------------------
r_cut <- 0.30
cor_df_filt <- cor_df %>% dplyr::filter(abs_weight >= r_cut)
cli::cli_alert_success("Edges kept at |r| ≥ {r_cut}: {nrow(cor_df_filt)}")

# For per-trait summaries, keep only trait–species pairs
edges_trait_species <- cor_df_filt %>%
  dplyr::mutate(
    trait = dplyr::case_when(
      from %in% physio_vars ~ from,
      to   %in% physio_vars ~ to,
      TRUE ~ NA_character_
    ),
    species = dplyr::case_when(
      from %in% physio_vars ~ to,
      to   %in% physio_vars ~ from,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(trait), !is.na(species)) %>%
  dplyr::arrange(trait, dplyr::desc(abs_weight))

K <- 10
trait_topk <- edges_trait_species %>%
  dplyr::group_by(trait) %>%
  dplyr::slice_head(n = K) %>%
  dplyr::ungroup() %>%
  dplyr::select(trait, species, r = weight, Effect, abs_weight)

# 7.5 Build network graph with pretty labels & shapes --------------------------
node_names <- sort(unique(c(cor_df_filt$from, cor_df_filt$to)))
is_trait   <- node_names %in% physio_vars

trait_label_map <- c(
  carb_mg_cm2    = "Carbohydrate",
  protein_mg_cm2 = "Protein",
  afdw_mg_cm2    = "Tissue~biomass",
  zoox_cells_cm2 = "Zooxanthellae",
  growth_vol_b   = "Growth"
  # percent_alive removed
)

label_parsed <- ifelse(
  is_trait,
  trait_label_map[node_names],
  paste0("italic(", gsub("_", "~", node_names), ")")
)

vertex_df <- tibble::tibble(
  name         = node_names,
  Type         = ifelse(is_trait, "Trait", "Species"),
  label_parsed = label_parsed
)

g <- igraph::graph_from_data_frame(
  d = cor_df_filt,
  directed = FALSE,
  vertices = vertex_df
)

set.seed(42)
lay <- igraph::layout_with_fr(g, weights = NA)

p_sem <- ggraph::ggraph(g, layout = "manual", x = lay[,1], y = lay[,2]) +
  # edges
  ggraph::geom_edge_link(
    ggplot2::aes(width = abs_weight, color = Effect),
    alpha = 0.65
  ) +
  ggraph::scale_edge_width_continuous(name = "|r|", range = c(0.3, 1.8)) +
  ggraph::scale_edge_color_manual(
    name = "Effect",
    values = c(Positive = "steelblue", Negative = "firebrick")
  ) +
  # nodes (triangles for traits, circles for species)
  ggraph::geom_node_point(
    ggplot2::aes(shape = Type, fill = Type),
    size = 3, stroke = 0.6, colour = "black"
  ) +
  ggplot2::scale_shape_manual(
    name = "Node",
    values = c(Trait = 24, Species = 21)
  ) +
  ggplot2::scale_fill_manual(
    name = "Node Fill",
    values = c(Trait = "white", Species = scales::alpha("gray50", 1))
  ) +
  # move labels LAST so they render above everything else
  ggraph::geom_node_text(
    ggplot2::aes(label = label_parsed),
    size = 3, parse = TRUE, repel = TRUE,
    bg.color = "white", bg.r = 0.15  # optional: halo background for clarity
  ) +
  ggplot2::labs(
    title    = "Correlation Network: All Species + Coral Traits",
    subtitle = paste0(
      "reef-adjusted residuals (random intercept for reef); edges shown for |r| ≥ ",
      r_cut
    )
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

print(p_sem)

# 7.6 Save outputs -------------------------------------------------------------
fs::dir_create(here::here("output","MRB","figures","network"))
fs::dir_create(here::here("output","MRB","objects"))
fs::dir_create(here::here("output","MRB","tables"))

ggplot2::ggsave(
  filename = here::here("output","MRB","figures","network","sem_network_all_species_r030_reefAdjusted.png"),
  plot = p_sem, width = 10, height = 8, dpi = 300
)

edges_for_graph <- cor_df_filt
saveRDS(
  list(
    sem_df          = sem_df,
    resid_df        = resid_df,
    sem_use         = sem_use,
    cor_mat         = cor_mat,
    cor_df          = cor_df,
    cor_df_filt     = cor_df_filt,
    edges_trait_sp  = edges_trait_species,
    trait_topk      = trait_topk,
    graph           = g
  ),
  file = here::here("output","MRB","objects","sem_network_all_species_r030_reefAdjusted.rds")
)

readr::write_csv(
  edges_for_graph,
  here::here("output","MRB","tables","sem_network_edges_r>=0.30_reefAdjusted.csv")
)


# =============================================================================
# 7.7 — Summaries for Methods & Results text
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(igraph)
  library(readr)
  library(cli)
  library(fs)
  library(here)
})

# Safety: defaults if not defined upstream
if (!exists("r_cut"))        r_cut        <- 0.30
if (!exists("physio_vars"))  physio_vars  <- c("growth_vol_b","protein_mg_cm2","carb_mg_cm2","zoox_cells_cm2","afdw_mg_cm2","percent_alive")
if (!exists("meta_cols"))    meta_cols    <- c("coral_id","reef_id","treatment")

# Output dir for tables
out_dir <- here::here("output","MRB","tables")
fs::dir_create(out_dir)

# --- A) What went into the analysis ------------------------------------------
n_corals_all <- if (exists("sem_df")) nrow(sem_df) else NA_integer_
n_traits     <- length(physio_vars)

# Species columns = everything not meta/traits (same rule used upstream)
if (exists("master_df")) {
  species_vars_all <- setdiff(names(master_df), c(meta_cols, physio_vars, grep("_z$", names(master_df), value = TRUE)))
  n_species_all    <- length(species_vars_all)
} else {
  species_vars_all <- character(0)
  n_species_all    <- NA_integer_
}

cli::cli_h2("Input overview")
cli::cli_alert_info("Corals (rows) used after drop_na: {n_corals_all}")
cli::cli_alert_info("Trait variables: {n_traits} → {paste(physio_vars, collapse=', ')}")
cli::cli_alert_info("Candidate species variables (pre-network): {n_species_all}")

# --- B) Edges/Nodes retained at threshold ------------------------------------
stopifnot(exists("cor_df_filt"), exists("g"))

n_edges_total   <- nrow(cor_df_filt)
n_edges_pos     <- sum(cor_df_filt$Effect == "Positive")
n_edges_neg     <- sum(cor_df_filt$Effect == "Negative")

# Vertex metadata
v_tbl <- tryCatch(
  igraph::as_data_frame(g, what = "vertices"),
  error = function(e) tibble::tibble(name = igraph::V(g)$name)
)

# Rebuild Type/label_parsed if missing
if (!"Type" %in% names(v_tbl)) {
  v_tbl <- v_tbl %>%
    dplyr::mutate(Type = ifelse(name %in% physio_vars, "Trait", "Species"))
}
if (!"label_parsed" %in% names(v_tbl)) {
  trait_label_map <- c(
    carb_mg_cm2    = "Carbohydrate~(mg~cm^{-2})",
    protein_mg_cm2 = "Protein~(mg~cm^{-2})",
    afdw_mg_cm2    = "Tissue~biomass~(AFDW~mg~cm^{-2})",
    zoox_cells_cm2 = "Zooxanthellae~(cells~cm^{-2})",
    growth_vol_b   = "Growth",
    percent_alive  = "Percent~alive"
  )
  v_tbl <- v_tbl %>%
    dplyr::mutate(label_parsed = ifelse(
      Type == "Trait",
      trait_label_map[name],
      paste0("italic(", gsub("_", "~", name), ")")
    ))
}

n_nodes_total   <- nrow(v_tbl)
n_nodes_traits  <- sum(v_tbl$Type == "Trait")
n_nodes_species <- sum(v_tbl$Type == "Species")

cli::cli_h2("Network overview (|r| >= {r_cut})")
cli::cli_alert_success("Edges retained: {n_edges_total} (Positive: {n_edges_pos}; Negative: {n_edges_neg})")
cli::cli_alert_success("Nodes in network: {n_nodes_total} (Traits: {n_nodes_traits}; Species: {n_nodes_species})")

# --- C) Strongest links (overall, +/-) ---------------------------------------
top_overall <- cor_df_filt %>%
  dplyr::arrange(dplyr::desc(abs_weight)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::select(rank, from, to, r = weight, Effect, abs_weight)

top_pos <- cor_df_filt %>%
  dplyr::filter(Effect == "Positive") %>%
  dplyr::arrange(dplyr::desc(abs_weight)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::select(rank, from, to, r = weight, abs_weight)

top_neg <- cor_df_filt %>%
  dplyr::filter(Effect == "Negative") %>%
  dplyr::arrange(dplyr::desc(abs_weight)) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::select(rank, from, to, r = weight, abs_weight)

cli::cli_h2("Top edges by |r|")
cli::cli_alert_info("Top 10 overall:\n{capture.output(print(utils::head(top_overall, 10), row.names = FALSE)) |> paste(collapse = '\n')}")
cli::cli_alert_info("Top 5 positive:\n{capture.output(print(utils::head(top_pos, 5), row.names = FALSE)) |> paste(collapse = '\n')}")
cli::cli_alert_info("Top 5 negative:\n{capture.output(print(utils::head(top_neg, 5), row.names = FALSE)) |> paste(collapse = '\n')}")

# --- D) For each TRAIT: top-k species links by |r| ---------------------------
K <- 10
edges_trait_species <- cor_df_filt %>%
  dplyr::mutate(
    trait = dplyr::case_when(
      from %in% physio_vars ~ from,
      to   %in% physio_vars ~ to,
      TRUE ~ NA_character_
    ),
    species = dplyr::case_when(
      from %in% physio_vars ~ to,
      to   %in% physio_vars ~ from,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(trait), !is.na(species)) %>%
  dplyr::arrange(trait, dplyr::desc(abs_weight))

trait_topk <- edges_trait_species %>%
  dplyr::group_by(trait) %>%
  dplyr::slice_head(n = K) %>%
  dplyr::ungroup() %>%
  dplyr::select(trait, species, r = weight, Effect, abs_weight)

# Per-trait one-liners without warnings or duplicates
if (exists("trait_topk") && nrow(trait_topk)) {
  trait_lines <- trait_topk %>%
    dplyr::group_by(trait) %>%
    dplyr::reframe(
      line = paste0(
        unique(trait), ": ",
        paste0(species, " (", sprintf("%.2f", r), ")", collapse = "; ")
      )
    ) %>%
    dplyr::ungroup()
  
  cat(paste0("  - ", trait_lines$line, collapse = "\n"), "\n")
}

cli::cli_h2("Per-trait top species links (|r| >= {r_cut}; top {K} per trait)")
for (tt in physio_vars) {
  tt_tbl <- trait_topk %>% dplyr::filter(trait == tt)
  if (nrow(tt_tbl)) {
    cli::cli_text("{.strong {tt}}:")
    print(tt_tbl, row.names = FALSE)
  } else {
    cli::cli_text("{.strong {tt}}: (no species at threshold)")
  }
}

# --- E) Centrality & degree summaries ----------------------------------------
# Degree (unweighted)
deg <- igraph::degree(g, mode = "all")

# Unweighted betweenness (explicitly ignore weights)
btw_unw <- igraph::betweenness(g, directed = FALSE, weights = NA, normalized = TRUE)

# Distance-weighted betweenness (recommended): define positive "cost"
E(g)$cost <- 1 / (E(g)$abs_weight + 1e-6)  # epsilon guards against 0
btw_w <- igraph::betweenness(g, directed = FALSE, weights = E(g)$cost, normalized = TRUE)

# Node strength (sum of absolute edge weights incident on node)
strength_abs <- igraph::strength(g, vids = igraph::V(g), mode = "all", weights = E(g)$abs_weight)

# Choose which betweenness to report
btw <- btw_w

cent_tbl <- v_tbl %>%
  dplyr::mutate(
    degree      = as.numeric(deg[name]),
    betweenness = as.numeric(btw[name]),
    strength    = as.numeric(strength_abs[name])
  ) %>%
  dplyr::arrange(dplyr::desc(degree), dplyr::desc(strength), dplyr::desc(betweenness))

top_nodes <- cent_tbl %>%
  dplyr::select(name, Type, degree, strength, betweenness) %>%
  dplyr::arrange(dplyr::desc(degree), dplyr::desc(strength), dplyr::desc(betweenness))

cli::cli_h2("Top nodes by degree (then strength, betweenness)")
print(utils::head(top_nodes, 15), row.names = FALSE)

# --- F) Export tidy tables for supplement ------------------------------------
edges_out <- cor_df_filt %>%
  dplyr::select(from, to, r = weight, abs_r = abs_weight, Effect) %>%
  dplyr::arrange(dplyr::desc(abs_r))

nodes_out <- v_tbl %>%
  dplyr::mutate(
    degree      = deg[name],
    betweenness = btw[name],
    strength    = strength_abs[name]
  ) %>%
  dplyr::select(name, Type, label_parsed, degree, strength, betweenness) %>%
  dplyr::arrange(dplyr::desc(degree), dplyr::desc(strength), dplyr::desc(betweenness))

trait_x_species_out <- edges_trait_species %>%
  dplyr::select(trait, species, r = weight, abs_r = abs_weight, Effect) %>%
  dplyr::arrange(trait, dplyr::desc(abs_r))

readr::write_csv(edges_out,           file.path(out_dir, "network_edges_r>=0.30.csv"))
readr::write_csv(nodes_out,           file.path(out_dir, "network_nodes_with_centrality.csv"))
readr::write_csv(trait_x_species_out, file.path(out_dir, "trait_species_edges_r>=0.30.csv"))
readr::write_csv(top_overall,         file.path(out_dir, "edges_top_overall_by_absr.csv"))
readr::write_csv(top_pos,             file.path(out_dir, "edges_top_positive_by_absr.csv"))
readr::write_csv(top_neg,             file.path(out_dir, "edges_top_negative_by_absr.csv"))

cli::cli_alert_success("Wrote tables → {out_dir}")

# --- G) One-liners you can paste into the manuscript --------------------------
cat("\n----- COPY/PASTE LINES -----\n")
cat(sprintf("We standardized all variables and retained Pearson correlations with |r| ≥ %.2f (n_edges = %d; positive = %d; negative = %d).\n",
            r_cut, n_edges_total, n_edges_pos, n_edges_neg))
cat(sprintf("The resulting network comprised %d nodes (%d trait variables; %d species variables).\n",
            n_nodes_total, n_nodes_traits, n_nodes_species))
if (nrow(top_overall) > 0) {
  cat(sprintf("Strongest association observed: %s – %s (r = %.2f).\n",
              top_overall$from[1], top_overall$to[1], top_overall$r[1]))
}
cat("Per-trait top species (|r|):\n")
if (exists("trait_topk") && nrow(trait_topk)) {
  trait_lines <- trait_topk %>%
    dplyr::group_by(trait) %>%
    dplyr::summarise(
      line = paste0(trait, ": ",
                    paste0(species, " (", sprintf("%.2f", r), ")", collapse = "; ")),
      .groups = "drop"
    )
  cat(paste0("  - ", trait_lines$line, collapse = "\n"), "\n")
}
cat("We computed unweighted degree, distance-weighted betweenness (edge cost = 1/|r|), and node strength (sum |r|).\n")
cat("----- END COPY/PASTE -----\n")


# --- 7.8 Save outputs ----------------------------------------------------------
fs::dir_create(here::here("output","MRB","figures","network"))
fs::dir_create(here::here("output","MRB","objects"))
fs::dir_create(here::here("output","MRB","tables"))

fig_path <- here::here("output","MRB","figures","network","sem_network_physio_panel.png")
ggplot2::ggsave(fig_path, p_sem, width = 8, height = 6, dpi = 300, bg = "white")
cli::cli_alert_success("Saved network figure → {fig_path}")

# Save underlying data (for reproducibility / re-plotting)
rds_path <- here::here("output","MRB","objects","sem_network_physio_panel.rds")
saveRDS(list(
  sem_df         = sem_df,
  sem_std        = sem_use,   # <- was sem_std
  cor_mat        = cor_mat,
  cor_df         = cor_df,
  cor_df_filt    = cor_df_filt,
  edges_for_graph= edges_for_graph,
  graph          = g
), rds_path)
# Optional: export edges table
tab_path <- here::here("output","MRB","tables","sem_network_edges_physio_panel.csv")
readr::write_csv(edges_for_graph, tab_path)
cli::cli_alert_success("Saved edges table → {tab_path}")

# Quick peek at strongest links
top_links <- cor_df_filt %>% dplyr::arrange(dplyr::desc(abs_weight)) %>% head(15)
cli::cli_alert_info("Top 15 |r| edges:\n{capture.output(print(top_links)) |> paste(collapse='\n')}")
