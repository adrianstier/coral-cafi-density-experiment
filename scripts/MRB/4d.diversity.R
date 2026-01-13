# ==============================================================================
# File:    4d.diversity.R
# Purpose: Community assembly analysis using a harmonized species set with
#          explicit 10×10 filtering (≥10 total individuals AND present on ≥10 corals).
#
# Distance flavors used throughout:
#   • Jaccard   — presence/absence
#   • Gower     — sqrt(counts) then z-score across species (no Hellinger)
#
# Style: publication-ready plots; treatment palette shared across scripts.
#
# --------------------------- QUESTIONS & APPROACH -----------------------------
# §0 Setup/Housekeeping
#   Q: How do we make the run reproducible and self-documenting?
#   A: Set seeds, load packages with checks, configure output dirs, define shared
#      themes/palettes/savers, and start robust console logging (tee).
#
# §1 Load & Join
#   Q: How are raw species counts joined to positions/treatments (reef, coral)?
#   A: Read CSVs, clean IDs, join, and construct a tidy base table.
#
# §2 Species filtering (the “10×10” rule)
#   Q: Which species are sufficiently common to analyze robustly?
#   A: Keep species with (prevalence ≥ 10 corals) AND (total abundance ≥ 10).
#      Propagate this filtered set consistently to all downstream objects.
#
# §3 Matrices & Distances
#   Q: What community matrices do we analyze?
#   A: Build reef- and coral-level matrices; create presence/absence, sqrt, and
#      sqrt→z flavors; compute Jaccard & Gower distances with aligned metadata.
#
# Later sections (α-diversity, rank-abundance, β-diversity PERMDISP/PERMANOVA,
# resampling, NMDS vectors, indicator species, sensitivity analyses) follow the
# same pattern: each starts with a brief “Q/A” header block explaining intent.
# ==============================================================================

# =============================================================================
# §0) SETUP / HOUSEKEEPING — reproducibility, packages, dirs, theme, helpers
# =============================================================================

# Source centralized libraries, utilities, and figure standards
source("scripts/MRB/1.libraries.R")
source("scripts/MRB/utils.R")  # Provides ALIVE_THRESH, strip_fe(), load functions
source("scripts/MRB/mrb_figure_standards.R")

set.seed(1234)
options(stringsAsFactors = FALSE,
        scipen = 999,
        dplyr.summarise.inform = FALSE,
        readr.show_col_types = FALSE)

# All packages loaded via source("scripts/MRB/1.libraries.R") above
# ALIVE_THRESH defined in utils.R (0.80)

# Paths
DATA_DIR <- here::here("data")
OUT_DIR  <- here::here("output", "MRB")
FIG_DIR  <- file.path(OUT_DIR, "figures", "diversity")
TAB_DIR  <- file.path(OUT_DIR, "tables")
invisible(lapply(list(FIG_DIR, TAB_DIR), dir.create, recursive = TRUE, showWarnings = FALSE))
cli::cli_alert_info("FIG_DIR: {FIG_DIR}")
cli::cli_alert_info("TAB_DIR: {TAB_DIR}")

# ---- Colors & Theme ----------------------------------------------------------
# Use TREATMENT_COLORS from mrb_figure_standards.R:
#   "1" = "#E69F00" (Orange)
#   "3" = "#56B4E9" (Sky Blue)
#   "6" = "#009E73" (Green)
# Use theme_publication() from mrb_figure_standards.R
# Use save_figure() from mrb_figure_standards.R or save_both() from utils.R
#
# Legacy aliases for backward compatibility (save_both, show_and_save provided by utils.R)
cols_trt <- TREATMENT_COLORS
theme_pub <- theme_publication

# ----------------------------- 10×10 FILTER CONFIG ---------------------------
# "10×10" means: keep species with prevalence ≥ PREV_MIN corals AND total abundance ≥ ABUND_MIN.
# Set USE_10x10_ONLY = TRUE to force the filtered species set across all downstream analyses.
FILTER_MODE    <- "threshold"
PREV_MIN       <- 10   # ≥10 corals with presence
ABUND_MIN      <- 10   # ≥10 total individuals
USE_10x10_ONLY <- TRUE # APPLY filtered set everywhere
FOCUS_LEVEL    <- "both"  # "reef", "coral", or "both"

# Resampling knobs (used in later sections)
N_ITER_PERM <- 1000
N_PER_TREAT <- 3

# ------------------------------ Utility helpers ------------------------------
# strip_fe() defined in utils.R - sourced above

# Guarded scaler: z-score columns; if a column has 0 variance, leave it as 0
safe_scale <- function(X) {
  X <- as.matrix(X)
  sds <- apply(X, 2, stats::sd, na.rm = TRUE)
  sds[sds == 0 | !is.finite(sds)] <- NA_real_
  Z <- sweep(X, 2, colMeans(X, na.rm = TRUE), "-")
  Z <- sweep(Z, 2, sds, "/")
  Z[, is.na(sds)] <- 0
  Z
}

# Build reef-level matrix from long data; presence=TRUE gives incidence
build_reef_matrix <- function(df, value_col = "count", presence = FALSE) {
  mat <- df %>%
    dplyr::group_by(reef, species) %>%
    dplyr::summarise(val = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = species, values_from = val, values_fill = 0) %>%
    tibble::column_to_rownames("reef") %>%
    as.matrix()
  if (presence) mat[] <- (mat > 0) * 1
  mat
}

# Align meta to a matrix’s rownames
meta_from_mat <- function(mat, df) {
  meta <- df %>% dplyr::select(reef, treatment) %>% dplyr::distinct()
  stopifnot(all(rownames(mat) %in% meta$reef))
  meta[match(rownames(mat), meta$reef), , drop = FALSE] %>%
    dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)))
}

# Centralized distance builders (reef-level)
distances_from_flavors <- function(mat_abund, mat_pa) {
  mat_sqrt   <- sqrt(mat_abund)
  mat_sqrt_z <- safe_scale(mat_sqrt)
  
  list(
    mat_abund   = mat_abund,
    mat_pa      = mat_pa,
    mat_sqrt    = mat_sqrt,
    mat_sqrt_z  = mat_sqrt_z,
    dist_jacc   = vegan::vegdist(mat_pa, method = "jaccard"),
    dist_gower  = as.dist(cluster::daisy(mat_sqrt_z, metric = "gower"))
  )
}

# Species filter (10×10) with a single source of truth
apply_species_filter <- function(spec_df, prev_min = PREV_MIN, abund_min = ABUND_MIN) {
  meta_tbl <- spec_df %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      total_abundance = sum(count, na.rm = TRUE),
      prevalence      = dplyr::n_distinct(coral_id[count > 0]),
      .groups = "drop"
    )
  keep <- meta_tbl %>%
    dplyr::filter(prevalence >= prev_min, total_abundance >= abund_min) %>%
    dplyr::pull(species) %>%
    unique() %>%
    sort()
  list(
    filtered_df = spec_df %>% dplyr::filter(species %in% keep),
    kept        = keep,
    dropped     = setdiff(meta_tbl$species, keep),
    meta_tbl    = meta_tbl
  )
}

# Handy labeler for treatment counts
label_treatment_counts <- function(meta_df) {
  tb <- table(meta_df$treatment)
  paste(sprintf("%s: n=%d", names(tb), as.integer(tb)), collapse = " | ")
}

# =============================================================================
# §1a) Coral inclusion filter: keep only corals with ≥80% live tissue
# =============================================================================
alive_src <- file.path(DATA_DIR, "MRB Amount", "coral_growth_surface_area_change_filtered.csv")

if (file.exists(alive_src)) {
  alive_meta <- readr::read_csv(alive_src, show_col_types = FALSE) %>%
    dplyr::mutate(coral_id = strip_fe(as.character(coral_id)))
  # expect a column named 'percent_alive' in 0–1 or 0–100; normalize if needed
  if (max(alive_meta$percent_alive, na.rm = TRUE) > 1.00001) {
    alive_meta <- dplyr::mutate(alive_meta, percent_alive = percent_alive / 100)
  }
} else {
  # Fallback: original manual file (adjust column names if yours differ)
  manual_file <- here::here("data","MRB Amount",
                            "1. amount_manual_colony_measurements_dec2019_and_may2021.xlsx")
  alive_meta <- readxl::read_excel(manual_file) %>%
    dplyr::transmute(
      coral_id      = strip_fe(as.character(coral_id)),
      percent_alive = percent_alive_may21 / 100
    )
}

stopifnot(all(c("coral_id","percent_alive") %in% names(alive_meta)))

keep_ids_alive <- alive_meta %>%
  dplyr::filter(!is.na(percent_alive), percent_alive >= ALIVE_THRESH) %>%
  dplyr::pull(coral_id) %>%
  unique()

cli::cli_alert_info("Alive filter (≥{ALIVE_THRESH*100}%): keeping {length(keep_ids_alive)} corals.")


# =============================================================================
# §1) LOAD & JOIN DATA
#   Q: How do we connect species counts to reef/coral/treatment metadata?
#   A: Read both CSVs, clean IDs, derive reef from position, and left-join.
# =============================================================================
cli::cli_h2("§1 Load & join data")

cafi_path  <- file.path(DATA_DIR, "MRB Amount", "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")
treat_path <- file.path(DATA_DIR, "MRB Amount", "coral_id_position_treatment.csv")
if (!file.exists(cafi_path))  stop("Missing: ", cafi_path)
if (!file.exists(treat_path)) stop("Missing: ", treat_path)

cafi_raw <- readr::read_csv(cafi_path) %>%
  dplyr::mutate(
    coral_id = strip_fe(coral_id),
    species  = trimws(species)
  ) %>%
  dplyr::filter(!is.na(coral_id), !is.na(species), !is.na(count))

treat_df <- readr::read_csv(treat_path) %>%
  dplyr::mutate(
    coral_id  = strip_fe(coral_id),
    row       = stringr::str_extract(position, "^\\d+") |> as.integer(),
    column    = stringr::str_extract(position, "(?<=-)\\d+") |> as.integer(),
    replicate = stringr::str_extract(position, "[A-Za-z]+"),
    reef      = glue::glue("Reef_{row}-{column}")
  ) %>%
  dplyr::select(coral_id, treatment, reef)

dat <- dplyr::left_join(cafi_raw, treat_df, by = "coral_id")

# >>> APPLY the ≥80% alive coral filter here <<<
dat <- dat %>% dplyr::filter(coral_id %in% keep_ids_alive)

if (any(is.na(dat$treatment)))
  cli::cli_alert_warning("Records missing treatment: {sum(is.na(dat$treatment))}")

species_only_df <- dat %>%
  dplyr::filter(!is.na(species), !is.na(reef), !is.na(treatment))

cli::cli_alert_success("Loaded {nrow(species_only_df)} species records with reef+treatment.")

# =============================================================================
# §2) SPECIES FILTERING — the 10×10 rule (global & consistent)
#   Q: Which species are sufficiently common to analyze robustly?
#   A: Keep species with (prevalence ≥ PREV_MIN corals) AND (total abundance ≥ ABUND_MIN).
#      If USE_10x10_ONLY = TRUE, all downstream objects use only this set.
# =============================================================================
cli::cli_h2("§2 10×10 species filtering")

flt <- apply_species_filter(species_only_df, prev_min = PREV_MIN, abund_min = ABUND_MIN)

# always keep these in the global env for downstream sections
filtered_df   <- flt$filtered_df
species_meta  <- flt$meta_tbl
.kept_species <- flt$kept

cli::cli_alert_success(
  "10×10 kept {length(.kept_species)} spp; dropped {length(flt$dropped)} spp. (prev ≥ {PREV_MIN}, total ≥ {ABUND_MIN})"
)

# Optionally enforce 10×10 everywhere
if (USE_10x10_ONLY) {
  species_only_df <- filtered_df
  cli::cli_alert_info("Using ONLY the 10×10 filtered species set downstream.")
} else {
  cli::cli_alert_info("Retaining all species downstream; some sections may still use the 10×10 subset explicitly.")
}

# Scope control (reef/coral/both)
if (FOCUS_LEVEL == "reef") {
  species_only_df <- species_only_df %>% dplyr::filter(!is.na(reef))
  cli::cli_alert_info("FOCUS_LEVEL=reef → keeping rows with reef IDs only.")
} else if (FOCUS_LEVEL == "coral") {
  species_only_df <- species_only_df %>% dplyr::filter(!is.na(coral_id))
  cli::cli_alert_info("FOCUS_LEVEL=coral → keeping rows with coral IDs only.")
} else {
  cli::cli_alert_info("FOCUS_LEVEL=both → reef and coral data remain available.")
}

# (Optional) write kept/dropped lists for records
readr::write_csv(
  species_meta %>% dplyr::mutate(kept = species %in% .kept_species),
  file.path(TAB_DIR, "02_species_filter_10x10_meta.csv")
)

# =============================================================================
# §3) REEF-LEVEL MATRICES & DISTANCES
#   Q: What are the core community matrices and distance flavors we analyze?
#   A: Build reef × species matrices for abundances and incidence; create sqrt
#      and sqrt→z versions; compute Jaccard (incidence) and Gower (sqrt→z).
# =============================================================================
cli::cli_h2("§3 Reef-level matrices & distances")

# reef-level matrices from the filtered (or full) working table
mat_abund_full <- build_reef_matrix(species_only_df, value_col = "count", presence = FALSE)
mat_pa_full    <- build_reef_matrix(species_only_df, value_col = "count", presence = TRUE)

# metadata aligned to those matrices
meta_abund <- meta_from_mat(mat_abund_full, species_only_df)
meta_pa    <- meta_from_mat(mat_pa_full,    species_only_df)

# bundle matrices + distances with centralized builder
D <- distances_from_flavors(mat_abund_full, mat_pa_full)

# expose these with the historical names you use later
mat_sqrt      <- D$mat_sqrt
mat_sqrt_z    <- D$mat_sqrt_z
dist_jaccard  <- D$dist_jacc
dist_gower    <- D$dist_gower

# quick snapshot to console
cli::cli_alert_info("Reef matrix dims — abund: {nrow(mat_abund_full)}×{ncol(mat_abund_full)}, PA: {nrow(mat_pa_full)}×{ncol(mat_pa_full)}")
cli::cli_alert_info("Treatment sizes (reef level): {label_treatment_counts(meta_abund)}")






# =============================================================================
# §4) α-DIVERSITY (REEF LEVEL)
#   Q: How does species richness and Shannon diversity vary across treatments?
#   A: Compute from reef-level matrices (filtered by 10×10 if USE_10x10_ONLY).
# =============================================================================
cli::cli_h2("§4 Reef α-diversity")

alpha_tbl <- tibble::tibble(
  reef       = rownames(mat_abund_full),
  richness   = vegan::specnumber(mat_abund_full),
  shannon    = vegan::diversity(mat_abund_full, index = "shannon"),
  treatment  = meta_abund$treatment
)

alpha_long <- tidyr::pivot_longer(alpha_tbl,
                                  cols = c(richness, shannon),
                                  names_to = "index",
                                  values_to = "value")

p_alpha <- ggplot(alpha_long, aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  facet_wrap(~ index, scales = "free_y") +
  labs(title = "§4 Reef α-diversity",
       x = "Treatment", y = "Diversity index") +
  theme_pub()
show_and_save(p_alpha, file.path(FIG_DIR, "04_alpha_diversity.png"),
              width = 7, height = 4, dpi = 600)

readr::write_csv(alpha_tbl, file.path(TAB_DIR, "04_alpha_diversity.csv"))
cli::cli_alert_success("§4 complete: α-diversity table + plot saved.")

# =============================================================================
# §5) RANK–ABUNDANCE CURVES
#   Q: Are assemblages dominated by a few taxa or evenly spread?
#   A: Plot log10 relative abundance by rank per treatment.
# =============================================================================
cli::cli_h2("§5 Rank–abundance curves")

rank_abund <- species_only_df %>%
  dplyr::group_by(species, treatment) %>%
  dplyr::summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(rank = dplyr::min_rank(desc(total)),
                rel  = total / sum(total)) %>%
  dplyr::ungroup()

p_rank <- ggplot(rank_abund, aes(rank, rel, color = treatment)) +
  geom_line() +
  scale_y_log10() +
  labs(title = "§5 Rank–abundance curves",
       x = "Rank", y = "Relative abundance (log10)") +
  theme_pub()
show_and_save(p_rank, file.path(FIG_DIR, "05_rank_abundance.png"),
              width = 7, height = 5, dpi = 600)

readr::write_csv(rank_abund, file.path(TAB_DIR, "05_rank_abundance.csv"))
cli::cli_alert_success("§5 complete: rank–abundance plot + table saved.")

# =============================================================================
# §6) β-DIVERSITY (REEF LEVEL)
#   Q: Do treatments differ in multivariate community composition?
#   A: Test homogeneity of dispersion (PERMDISP), main differences (PERMANOVA),
#      and visualize in NMDS ordinations.
# =============================================================================
cli::cli_h2("§6 Reef β-diversity")

# 6A) PERMDISP (homogeneity of multivariate dispersion)
betadisp_jac <- vegan::betadisper(dist_jaccard, group = meta_abund$treatment)
betadisp_gow <- vegan::betadisper(dist_gower,   group = meta_abund$treatment)

disp_test_jac <- vegan::permutest(betadisp_jac, permutations = 999)
disp_test_gow <- vegan::permutest(betadisp_gow, permutations = 999)

cli::cli_alert_info("PERMDISP Jaccard: F={round(disp_test_jac$tab[1,4],3)} p={disp_test_jac$tab[1,5]}")
cli::cli_alert_info("PERMDISP Gower:   F={round(disp_test_gow$tab[1,4],3)} p={disp_test_gow$tab[1,5]}")

# 6B) PERMANOVA (adonis2)
permanova_jac <- vegan::adonis2(dist_jaccard ~ treatment, data = meta_abund, permutations = 999)
permanova_gow <- vegan::adonis2(dist_gower   ~ treatment, data = meta_abund, permutations = 999)

cli::cli_alert_info("PERMANOVA Jaccard: R²={round(permanova_jac$R2[1],3)} p={permanova_jac$`Pr(>F)`[1]}")
cli::cli_alert_info("PERMANOVA Gower:   R²={round(permanova_gow$R2[1],3)} p={permanova_gow$`Pr(>F)`[1]}")

# 6C) NMDS ordinations
set.seed(123)
nmds_jac <- vegan::metaMDS(dist_jaccard, k = 2, trymax = 200)
nmds_gow <- vegan::metaMDS(dist_gower,   k = 2, trymax = 200)

nmds_df <- function(ord, meta, label) {
  as.data.frame(ord$points) %>%
    tibble::rownames_to_column("reef") %>%
    dplyr::left_join(meta, by = "reef") %>%
    dplyr::mutate(metric = label)
}

df_nmds <- dplyr::bind_rows(
  nmds_df(nmds_jac, meta_abund, "Jaccard"),
  nmds_df(nmds_gow, meta_abund, "Gower (sqrt→z)")
)

p_nmds <- ggplot(df_nmds, aes(MDS1, MDS2, color = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ metric) +
  labs(title = "§6 NMDS ordinations by treatment") +
  theme_pub()
show_and_save(p_nmds, file.path(FIG_DIR, "06_nmds.png"),
              width = 8, height = 4, dpi = 600)


# Helper: write adonis2 table to CSV
write_adonis <- function(fit, path) {
  df <- as.data.frame(fit)
  df <- tibble::rownames_to_column(df, "term")
  readr::write_csv(df, path)
}


# Helper: save PERMDISP (permutest on betadisper) results to CSV
save_disp_results <- function(permtest_obj, stub, out_dir) {
  # Global test table
  df <- as.data.frame(permtest_obj$tab)
  df <- tibble::rownames_to_column(df, "term")
  readr::write_csv(df, file.path(out_dir, paste0(stub, "_global.csv")))
}

# Save PERMANOVA + PERMDISP results to CSV
write_adonis(permanova_jac, file.path(TAB_DIR, "06_permanova_jaccard.csv"))
write_adonis(permanova_gow, file.path(TAB_DIR, "06_permanova_gower.csv"))
save_disp_results(disp_test_jac, "06_permdisp_jaccard", TAB_DIR)
save_disp_results(disp_test_gow, "06_permdisp_gower", TAB_DIR)

cli::cli_alert_success("§6 complete: PERMDISP, PERMANOVA, NMDS results saved.")




# =============================================================================
# §14) CORAL-LEVEL COMMUNITY STRUCTURE + PERMANOVAs
#   Qs:
#     (a) What does coral-level structure look like across treatments?
#     (b) Are treatment effects robust when permutations are constrained within reefs?
#     (c) Do conclusions match unit-correct reef-level models (corals aggregated)?
#   Notes:
#     - Uses the 10×10 species set if USE_10x10_ONLY = TRUE (from Block 2).
#     - Gower flavor = sqrt(counts) then z-score across species (no Hellinger).
# =============================================================================
cli::cli_h2("§14 Coral-level NMDS + PERMANOVAs (strata = reef)")

# ---- 14.0 Guards -------------------------------------------------------------
.required <- c("species_only_df", "flt", "FIG_DIR", "TAB_DIR", "cols_trt",
               "theme_pub", "show_and_save")
if (!all(vapply(.required, exists, logical(1)))) {
  stop("§14 requires objects from prior blocks: ", paste(.required, collapse = ", "))
}
.kept_species <- flt$kept
if (!length(.kept_species)) stop("§14: no species passed the 10×10 filter.")

# ---- 14.1 Coral-level matrices (filtered species) ---------------------------
meta_coral <- species_only_df %>%
  dplyr::select(coral_id, reef, treatment) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    treatment = factor(treatment, levels = c(1,3,6)),
    reef      = factor(reef)
  )

coral_long <- species_only_df %>%
  dplyr::filter(species %in% .kept_species) %>%
  dplyr::group_by(coral_id, species) %>%
  dplyr::summarise(val = sum(count, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(
    coral_id = meta_coral$coral_id,
    species  = .kept_species,
    fill     = list(val = 0)
  )

mat_coral_abund <- coral_long %>%
  tidyr::pivot_wider(names_from = species, values_from = val, values_fill = 0) %>%
  tibble::column_to_rownames("coral_id") %>%
  as.matrix()

# Align metadata row-for-row with matrix
stopifnot(all(rownames(mat_coral_abund) %in% meta_coral$coral_id))
meta_coral <- meta_coral[match(rownames(mat_coral_abund), meta_coral$coral_id), , drop = FALSE]

# Derived matrices
mat_coral_pa     <- (mat_coral_abund > 0) * 1
mat_coral_sqrt_z <- scale(sqrt(mat_coral_abund))  # Gower flavor

# Distances (coral level)
dist_coral_jac <- vegan::vegdist(mat_coral_pa, method = "jaccard")
dist_coral_gow <- as.dist(cluster::daisy(mat_coral_sqrt_z, metric = "gower"))



# ---- Build 'scal' (coral NMDS scores for both metrics) ----------------------

# Ordinations
nmds_coral_jac <- vegan::metaMDS(dist_coral_jac, k = 2, trymax = 200, trace = FALSE)
nmds_coral_gow <- vegan::metaMDS(dist_coral_gow, k = 2, trymax = 200, trace = FALSE)

.scores_coral <- function(ord, metric_label) {
  df <- vegan::scores(ord, display = "sites") %>%
    as.data.frame()
  # Standardize column names to NMDS1/NMDS2 regardless of what vegan gives us
  if (all(c("MDS1","MDS2") %in% names(df))) {
    df <- dplyr::rename(df, NMDS1 = MDS1, NMDS2 = MDS2)
  } else if (all(c("NMDS1","NMDS2") %in% names(df))) {
    # already good
  } else {
    # fall back: take first two numeric columns as axes
    ax <- names(df)[vapply(df, is.numeric, logical(1))]
    stopifnot(length(ax) >= 2)
    df <- dplyr::rename(df, NMDS1 = !!ax[1], NMDS2 = !!ax[2])
  }
  df %>%
    tibble::rownames_to_column("coral_id") %>%
    dplyr::left_join(meta_coral, by = "coral_id") %>%
    dplyr::mutate(metric = metric_label)
}

# Combined scores for Jaccard and Gower
scal <- dplyr::bind_rows(
  .scores_coral(nmds_coral_jac, "Jaccard"),
  .scores_coral(nmds_coral_gow, "Gower (sqrt→z)")
)

# ---- 14.2 NMDS (coral level) — add reef centroids + spiders ------------------
# 'scal' already exists and contains coral scores for both metrics (Jaccard/Gower)
# with columns: coral_id, NMDS1, NMDS2, reef, treatment, metric

# 1) Compute reef centroids within each metric (mean of coral scores per reef)
centroids_coral <- scal %>%
  dplyr::group_by(metric, reef) %>%
  dplyr::summarise(
    cx = mean(NMDS1, na.rm = TRUE),
    cy = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    scal %>% dplyr::distinct(reef, treatment), by = "reef"
  )
# (centroids get the reef's treatment for coloring/legend consistency)

# 2) Make a segment table that joins each coral to its reef centroid (per metric)
segments_coral <- scal %>%
  dplyr::left_join(centroids_coral, by = c("metric","reef")) %>%
  dplyr::filter(is.finite(NMDS1), is.finite(NMDS2), is.finite(cx), is.finite(cy))

# 3) Plot with spiders + centroids + points
p_coral_nmds <- ggplot(scal, aes(NMDS1, NMDS2, color = treatment)) +
  # spiders: coral → reef centroid (light gray so they don’t overpower)
  geom_segment(
    data = segments_coral,
    aes(x = NMDS1, y = NMDS2, xend = cx, yend = cy, group = reef),
    inherit.aes = FALSE,
    color = "grey70", linewidth = 0.3, alpha = 0.6
  ) +
  # reef centroids: larger points with black outline
  geom_point(
    data = centroids_coral,
    aes(x = cx, y = cy, fill = treatment),
    inherit.aes = FALSE,
    shape = 21, size = 3.5, color = "black", stroke = 0.4, alpha = 0.95
  ) +
  # individual corals
  geom_point(size = 1.6, alpha = 0.9) +
  facet_wrap(~ metric) +
  scale_color_manual(values = cols_trt, name = "Treatment") +
  scale_fill_manual(values = cols_trt, guide = "none") +
  coord_equal(expand = TRUE) +
  labs(
    title = "§14 Coral-level NMDS by treatment (with reef centroids & spiders)",
    subtitle = paste0("stress: Jaccard=", round(nmds_coral_jac$stress, 3),
                      " ; Gower=", round(nmds_coral_gow$stress, 3)),
    x = "NMDS1", y = "NMDS2"
  ) +
  theme_pub() +
  theme(legend.position = "top")

show_and_save(p_coral_nmds, file.path(FIG_DIR, "14_coral_nmds_centroids_spiders.png"),
              width = 9, height = 4.8, dpi = 600)

# ---- 14.3 PERMANOVAs (coral level, restricted perms within reefs) -----------
# IMPORTANT: strata = reef to avoid pseudo-replication of corals within reefs
permanova_coral_jac <- vegan::adonis2(
  dist_coral_jac ~ treatment,
  data = meta_coral,
  permutations = 999,
  strata = meta_coral$reef
)
permanova_coral_gow <- vegan::adonis2(
  dist_coral_gow ~ treatment,
  data = meta_coral,
  permutations = 999,
  strata = meta_coral$reef
)

# ---- 14.4 Unit-correct reef-level models (corals aggregated to reefs) -------
agg_to_reef <- function(mat_coral, meta) {
  mm <- as.data.frame(mat_coral) %>%
    tibble::rownames_to_column("coral_id") %>%
    dplyr::left_join(meta, by = "coral_id") %>%
    dplyr::select(-coral_id) %>%
    dplyr::group_by(reef, treatment) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE),
                     .groups = "drop")
  X <- mm %>% dplyr::select(-reef, -treatment) %>% as.matrix()
  rownames(X) <- mm$reef
  list(X = X, meta = mm %>% dplyr::select(reef, treatment))
}

reef_abund    <- agg_to_reef(mat_coral_abund, meta_coral)       # abundance
reef_pa       <- (reef_abund$X > 0) * 1
reef_sqrt_z   <- scale(sqrt(reef_abund$X))

dist_reef_jac <- vegan::vegdist(reef_pa, method = "jaccard")
dist_reef_gow <- as.dist(cluster::daisy(reef_sqrt_z, metric = "gower"))
meta_reef     <- reef_abund$meta %>% dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)))

permanova_reef_jac <- vegan::adonis2(dist_reef_jac ~ treatment, data = meta_reef, permutations = 999)
permanova_reef_gow <- vegan::adonis2(dist_reef_gow ~ treatment, data = meta_reef, permutations = 999)

# ---- 14.5 Pretty tables + disk outputs --------------------------------------
# helper: format adonis2 row 1 (treatment)
.adonis_row1 <- function(fit) {
  tab <- as.data.frame(fit)
  tibble::tibble(
    Term = rownames(tab)[1],
    Df   = suppressWarnings(as.numeric(tab$Df[1])),
    F    = suppressWarnings(as.numeric(tab$F[1])),
    R2   = suppressWarnings(as.numeric(tab$R2[1])),
    p    = suppressWarnings(as.numeric(tab[1, "Pr(>F)"]))
  )
}

sum_tbl <- dplyr::bind_rows(
  .adonis_row1(permanova_coral_jac) %>% dplyr::mutate(Level = "Coral (strata=reef)", Metric = "Jaccard"),
  .adonis_row1(permanova_coral_gow) %>% dplyr::mutate(Level = "Coral (strata=reef)", Metric = "Gower (sqrt→z)"),
  .adonis_row1(permanova_reef_jac)  %>% dplyr::mutate(Level = "Reef (unit-correct)", Metric = "Jaccard"),
  .adonis_row1(permanova_reef_gow)  %>% dplyr::mutate(Level = "Reef (unit-correct)", Metric = "Gower (sqrt→z)")
) %>%
  dplyr::relocate(Level, Metric, Term, Df, F, R2, p) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ ""
  ))

readr::write_csv(sum_tbl, file.path(TAB_DIR, "14_permanova_summary.csv"))

gt_tbl <- gt::gt(sum_tbl, groupname_col = "Level") %>%
  gt::cols_label(Metric = "Distance", Term = "Term", Df = "df", F = "F",
                 R2 = "R²", p = "p", sig = "Sig.") %>%
  gt::fmt_number(columns = c(Df, F, R2, p), decimals = 3) %>%
  gt::tab_header(title = "§14 — Treatment effects: coral (strata=reef) vs reef (unit-correct)") %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold"),
    locations = gt::cells_body(columns = sig, rows = sig %in% c("***","**","*"))
  )
gt::gtsave(gt_tbl, file.path(TAB_DIR, "14_permanova_summary_gt.html"))

cli::cli_alert_success("§14 complete: coral NMDS, coral PERMANOVAs (strata=reef), reef PERMANOVAs, and GT table saved.")





# =============================================================================
# §15) TREATMENT EFFECTS AT THE REEF LEVEL
#   Qs:
#     (a) Is there an overall treatment effect among reefs (PERMANOVA)?
#     (b) Which treatment pairs differ (pairwise PERMANOVA, BH-adjusted)?
#     (c) Are group dispersions comparable (PERMDISP)? If not, flag caution.
#     (d) How big are pairwise R² and how do between-reef distances distribute?
#   Notes:
#     - Uses reef-level distance objects from §14 (unit-correct).
#     - Distances: Jaccard (incidence) and Gower (sqrt→z on abundance).
# =============================================================================
cli::cli_h2("§15 Treatment effects among reefs (PERMANOVA, PERMDISP, visuals)")

# ---- 15.0 Guards -------------------------------------------------------------
.required15 <- c("meta_reef","dist_reef_jac","dist_reef_gow","FIG_DIR","TAB_DIR",
                 "theme_pub","show_and_save","cols_trt")
if (!all(vapply(.required15, exists, logical(1)))) {
  stop("§15 requires: ", paste(.required15, collapse = ", "),
       "\nRun Blocks 1–4 first.")
}
meta_reef <- meta_reef %>% dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)))

# ---- 15.1 Utilities (scoped to §15 to avoid name clashes) -------------------
align_meta_to_dist15 <- function(dist_obj, meta, reef_col = "reef") {
  labs <- attr(dist_obj, "Labels")
  if (is.null(labs)) stop("Distance object has no 'Labels'.")
  meta2 <- meta %>% dplyr::mutate(reef = as.character(.data[[reef_col]]))
  idx   <- match(labs, meta2$reef)
  if (any(is.na(idx))) {
    stop("Reefs in distance not found in meta: ", paste(setdiff(labs, meta2$reef), collapse = ", "))
  }
  meta2[idx, , drop = FALSE]
}
pval_to_stars15 <- function(p) dplyr::case_when(
  is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ ""
)
adonis_row1_15 <- function(fit) {
  tab <- as.data.frame(fit)
  tibble::tibble(
    df = suppressWarnings(as.numeric(tab$Df[1])),
    F  = suppressWarnings(as.numeric(tab$F[1])),
    R2 = suppressWarnings(as.numeric(tab$R2[1])),
    p  = suppressWarnings(as.numeric(tab[1,"Pr(>F)"]))
  )
}

pairwise_permanova15 <- function(dist_obj, meta, group_var = "treatment", nperm = 999) {
  metaA <- align_meta_to_dist15(dist_obj, meta, reef_col = "reef")
  g     <- droplevels(factor(metaA[[group_var]]))
  lev   <- levels(g); L <- length(lev)
  Dmat  <- as.matrix(dist_obj)
  out   <- vector("list", choose(L, 2)); k <- 1
  for (i in 1:(L - 1)) for (j in (i + 1):L) {
    keep <- metaA[[group_var]] %in% c(lev[i], lev[j])
    dsub <- as.dist(Dmat[keep, keep, drop = FALSE])
    msub <- droplevels(metaA[keep, , drop = FALSE])
    msub[[group_var]] <- droplevels(factor(msub[[group_var]]))
    fmla <- stats::as.formula(sprintf("dsub ~ %s", group_var))
    fit  <- vegan::adonis2(fmla, data = msub, permutations = nperm)
    r1   <- adonis_row1_15(fit)
    out[[k]] <- tibble::tibble(
      group1 = as.character(lev[i]),
      group2 = as.character(lev[j]),
      df = r1$df, F = r1$F, R2 = r1$R2, p = r1$p
    )
    k <- k + 1
  }
  dplyr::bind_rows(out)
}

pairwise_permdisp15 <- function(dist_obj, meta, group_var = "treatment", nperm = 999) {
  metaA <- align_meta_to_dist15(dist_obj, meta, reef_col = "reef")
  bd <- vegan::betadisper(dist_obj, droplevels(factor(metaA[[group_var]])))
  pt <- vegan::permutest(bd, permutations = nperm)  # global F and p
  tk <- tryCatch(TukeyHSD(bd), error = function(e) NULL)
  if (!is.null(tk) && length(tk)) {
    df <- as.data.frame(tk[[1]])
    df$contrast <- rownames(df)
    pw <- df %>%
      tibble::as_tibble() %>%
      tidyr::separate(contrast, into = c("group1","group2"), sep = "-") %>%
      dplyr::transmute(group1, group2,
                       diff = .data$diff, lwr = .data$lwr, upr = .data$upr, p_adj = .data$`p adj`)
  } else {
    pw <- tibble::tibble(group1=character(), group2=character(),
                         diff=numeric(), lwr=numeric(), upr=numeric(), p_adj=numeric())
  }
  list(global = pt, pairwise = pw, bd = bd)
}

plot_pw_heatmap15 <- function(df_metric, title_) {
  ggplot(df_metric, aes(x = pair, y = "R²", fill = R2)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(sprintf("%.2f", R2), "\n", stars)), size = 3) +
    scale_fill_viridis_c(option = "D") +
    labs(title = title_, x = NULL, y = NULL, fill = "R² (pairwise)") +
    theme_pub() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

tidy_dist_by_treatment15 <- function(dist_obj, meta) {
  metaA <- align_meta_to_dist15(dist_obj, meta, reef_col = "reef")
  labs  <- attr(dist_obj, "Labels"); M <- as.matrix(dist_obj)
  comb  <- utils::combn(labs, 2)
  dvals <- mapply(function(a,b) M[a,b], comb[1,], comb[2,])
  tibble::tibble(reef1 = comb[1,], reef2 = comb[2,], d = as.numeric(dvals)) %>%
    dplyr::left_join(metaA, by = c("reef1" = "reef")) %>% dplyr::rename(trt1 = treatment) %>%
    dplyr::left_join(metaA, by = c("reef2" = "reef")) %>% dplyr::rename(trt2 = treatment) %>%
    dplyr::mutate(contrast = paste(pmin(trt1,trt2), "vs", pmax(trt1,trt2)))
}

plot_dist_box15 <- function(df, title_) {
  ggplot(df, aes(x = contrast, y = d)) +
    geom_boxplot(outlier.shape = NA, fill = "grey85", color = "grey30") +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
    labs(title = title_, x = "Treatment contrast (reef pairs)", y = "Dissimilarity") +
    theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# ---- 15.2 Main PERMANOVA (reef level) ---------------------------------------
mainmod_jac <- vegan::adonis2(dist_reef_jac ~ treatment, data = meta_reef, permutations = 999)
mainmod_gow <- vegan::adonis2(dist_reef_gow ~ treatment, data = meta_reef, permutations = 999)

jac1 <- adonis_row1_15(mainmod_jac)
gow1 <- adonis_row1_15(mainmod_gow)

mainmod_df <- tibble::tibble(
  metric = c("Jaccard (incidence)", "Gower (sqrt→z)"),
  df = c(jac1$df, gow1$df),
  F  = c(jac1$F,  gow1$F),
  R2 = c(jac1$R2, gow1$R2),
  p  = c(jac1$p,  gow1$p)
) %>% dplyr::mutate(sig = pval_to_stars15(p))

print(mainmod_df)

gt_main <- mainmod_df %>%
  gt::gt() %>%
  gt::tab_header(title = "§15 — Main PERMANOVA (reef-level treatment effect)") %>%
  gt::cols_label(metric="Distance metric", df="df", F="F", R2="R²", p="p", sig="Sig.") %>%
  gt::fmt_number(columns = c(df, F, R2, p), decimals = 3) %>%
  gt::tab_style(gt::cell_text(weight="bold"),
                locations = gt::cells_body(columns = sig, rows = sig %in% c("***","**","*")))
gt::gtsave(gt_main, file.path(TAB_DIR, "15_main_permanova_gt.html"))

# ---- 15.3 Pairwise PERMANOVA -------------------------------------------------
pw_jac <- pairwise_permanova15(dist_reef_jac, meta_reef, "treatment", nperm = 999) %>%
  dplyr::mutate(metric = "Jaccard (incidence)")
pw_gow <- pairwise_permanova15(dist_reef_gow, meta_reef, "treatment", nperm = 999) %>%
  dplyr::mutate(metric = "Gower (sqrt→z)")

pw_all <- dplyr::bind_rows(pw_jac, pw_gow) %>%
  dplyr::mutate(p_adj = p.adjust(p, method = "BH"),
                stars = pval_to_stars15(p_adj))

print(pw_all)
readr::write_csv(pw_all, file.path(TAB_DIR, "15_pairwise_permanova_reef_treatment.csv"))

gt_pw <- pw_all %>%
  dplyr::mutate(pair = paste(group1, "vs", group2)) %>%
  dplyr::select(metric, pair, df, F, R2, p, p_adj, stars) %>%
  gt::gt(groupname_col = "metric") %>%
  gt::tab_header(title = "§15 — Pairwise PERMANOVA among reefs (BH-adjusted p)") %>%
  gt::fmt_number(columns = c(df, F, R2, p, p_adj), decimals = 3) %>%
  gt::cols_label(pair="Contrast", p_adj="p (BH)", stars="Sig.") %>%
  gt::tab_style(gt::cell_text(weight="bold"),
                locations = gt::cells_body(columns = stars, rows = stars %in% c("***","**","*")))
gt::gtsave(gt_pw, file.path(TAB_DIR, "15_pairwise_permanova_gt.html"))

# ---- 15.4 PERMDISP (global + Tukey) -----------------------------------------
disp_jac <- pairwise_permdisp15(dist_reef_jac, meta_reef, "treatment", nperm = 999)
disp_gow <- pairwise_permdisp15(dist_reef_gow, meta_reef, "treatment", nperm = 999)

get_permF15 <- function(pt) {
  tab <- pt$tab
  F  <- suppressWarnings(as.numeric(tab[1, "F"]))
  p  <- suppressWarnings(as.numeric(tab[1, "Pr(>F)"]))
  list(F=F, p=p)
}
gp_j <- get_permF15(disp_jac$global); gp_g <- get_permF15(disp_gow$global)

permdisp_global_df <- tibble::tibble(
  metric = c("Jaccard (incidence)", "Gower (sqrt→z)"),
  F      = c(gp_j$F, gp_g$F),
  p      = c(gp_j$p, gp_g$p)
) %>% dplyr::mutate(sig = pval_to_stars15(p))

gt_permdisp <- permdisp_global_df %>%
  gt::gt() %>%
  gt::tab_header(title = "§15 — PERMDISP: homogeneity of multivariate dispersion (reef level)") %>%
  gt::cols_label(metric="Distance metric", F="F", p="p", sig="Sig.") %>%
  gt::fmt_number(columns = c(F, p), decimals = 3) %>%
  gt::tab_style(gt::cell_text(weight="bold"),
                locations = gt::cells_body(columns = sig, rows = sig %in% c("***","**","*")))
gt::gtsave(gt_permdisp, file.path(TAB_DIR, "15_permdisp_overall_gt.html"))

# Save full PERMDISP printouts (global + pairwise Tukey) as text
sink(file.path(TAB_DIR, "15_permdisp_pairwise_details.txt"))
cat("PERMDISP global tests (reef level):\n\n")
cat("Jaccard:\n"); print(disp_jac$global); cat("\n---\n")
cat("Gower:\n");   print(disp_gow$global); cat("\n\n")
cat("Pairwise Tukey on distances to centroid (Jaccard):\n"); print(disp_jac$pairwise); cat("\n---\n")
cat("Pairwise Tukey on distances to centroid (Gower):\n");   print(disp_gow$pairwise); cat("\n")
sink()

# ---- 15.5 Visual summaries ---------------------------------------------------
heatmap_df <- pw_all %>%
  dplyr::mutate(pair = paste(group1, "vs", group2)) %>%
  dplyr::select(metric, pair, R2, p_adj, stars)

p_heat_jac <- plot_pw_heatmap15(dplyr::filter(heatmap_df, metric == "Jaccard (incidence)"),
                                "Pairwise PERMANOVA R² — Jaccard (reef level)")
p_heat_gow <- plot_pw_heatmap15(dplyr::filter(heatmap_df, metric == "Gower (sqrt→z)"),
                                "Pairwise PERMANOVA R² — Gower (reef level)")

p_heat <- p_heat_jac / p_heat_gow +
  patchwork::plot_annotation(title = "§15 Pairwise treatment contrasts among reefs (R² with BH-adjusted significance)")
show_and_save(p_heat, file.path(FIG_DIR, "15_pairwise_permanova_heatmaps.png"), width = 12, height = 6)

td_jac <- tidy_dist_by_treatment15(dist_reef_jac, meta_reef)
td_gow <- tidy_dist_by_treatment15(dist_reef_gow, meta_reef)

p_box_jac <- plot_dist_box15(td_jac, "Between-reef dissimilarities by treatment contrast — Jaccard")
p_box_gow <- plot_dist_box15(td_gow, "Between-reef dissimilarities by treatment contrast — Gower (sqrt→z)")
p_boxes   <- p_box_jac / p_box_gow + patchwork::plot_annotation(
  title = "§15 Distribution of between-reef distances across treatment contrasts"
)
show_and_save(p_boxes, file.path(FIG_DIR, "15_between_reef_distance_boxplots.png"), width = 11, height = 8)

# ---- 15.6 Notes file for reporting ------------------------------------------
sink(file.path(TAB_DIR, "15_reporting_notes.txt"))
cat("§15 — Reporting notes\n\n")
cat("* Main PERMANOVA table: 15_main_permanova_gt.html\n")
cat("* Pairwise PERMANOVA CSV: 15_pairwise_permanova_reef_treatment.csv\n")
cat("* Pairwise PERMANOVA GT:  15_pairwise_permanova_gt.html\n")
cat("* PERMDISP global GT:     15_permdisp_overall_gt.html\n")
cat("* Always interpret PERMANOVA alongside PERMDISP (dispersion differences can inflate Type I error).\n")
sink()

cli::cli_alert_success("§15 complete: main + pairwise PERMANOVA, PERMDISP, heatmaps and boxplots saved.")





# =============================================================================
# §16) Gower contrasts: (Abundances vs Proportions) × (Coral vs Reef) — hull NMDS
# -----------------------------------------------------------------------------
# Questions:
#   Q1. Does using abundances vs row-wise proportions change separation by treatment?
#   Q2. Are conclusions consistent at coral and reef scales?
#   Q3. Are treatment effects significant in each combination (PERMANOVA)?
# Approach:
#   - Build four matrices: Coral×Abund, Coral×Prop, Reef×Abund, Reef×Prop.
#   - Distances: Gower on z-scored columns; for Abund use sqrt(counts)→z; for Prop use row-proportions→z.
#   - Ordinate (NMDS, k=2), overlay convex hulls by treatment (no ellipses).
#   - PERMANOVA: coral-level uses strata = reef; reef-level uses reefs as units.
#   - Save 2×2 NMDS panel + PERMANOVA summary (CSV + GT).
# =============================================================================
cli::cli_h2("§16 Gower contrasts: Abundance vs Proportion × Coral vs Reef")

# ---- 16.0 Guards -------------------------------------------------------------
.required16 <- c("FIG_DIR","TAB_DIR","cols_trt",
                 "mat_coral_abund","meta_coral",
                 "reef_abund","meta_reef","theme_pub","show_and_save")
if (!all(vapply(.required16, exists, logical(1)))) {
  stop("§16 requires objects from §14/§15: mat_coral_abund, meta_coral, reef_abund, meta_reef, FIG_DIR, TAB_DIR, cols_trt.")
}
meta_coral <- meta_coral %>% dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)))
meta_reef  <- meta_reef  %>% dplyr::mutate(treatment = factor(treatment,  levels = c(1,3,6)))

# ---- 16.1 Helpers (scoped to §16) -------------------------------------------
row_props16 <- function(X) {                       # row-wise proportions with safe zero guard
  X <- as.matrix(X)
  rs <- rowSums(X, na.rm = TRUE); rs[rs == 0] <- 1
  sweep(X, 1, rs, "/")
}

# Robust Gower builders: sqrt->z (abund) or prop->z (prop), drop zero-variance cols
gower_from_abund16 <- function(X) {
  Z <- scale(sqrt(as.matrix(X)))
  Z[!is.finite(Z)] <- 0
  sds <- apply(Z, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  if (!any(keep)) stop("gower_from_abund16: no varying species columns remain.")
  as.dist(cluster::daisy(Z[, keep, drop = FALSE], metric = "gower"))
}
gower_from_prop16  <- function(X) {
  Z <- scale(as.matrix(X))
  Z[!is.finite(Z)] <- 0
  sds <- apply(Z, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  if (!any(keep)) stop("gower_from_prop16: no varying species columns remain.")
  as.dist(cluster::daisy(Z[, keep, drop = FALSE], metric = "gower"))
}

# convex hull per group (requires n ≥ 3)
hull_by_group16 <- function(df, grp = "treatment", x = "NMDS1", y = "NMDS2") {
  df |>
    dplyr::group_by(.data[[grp]]) |>
    dplyr::filter(dplyr::n() >= 3, is.finite(.data[[x]]), is.finite(.data[[y]])) |>
    dplyr::slice(chull(.data[[x]], .data[[y]])) |>
    dplyr::ungroup()
}

# unified NMDS plot (hulls + points) for coral or reef meta
plot_gower_nmds_hulls16 <- function(dist_obj, meta_df, title_stub, palette = cols_trt) {
  set.seed(42)
  fit <- vegan::metaMDS(dist_obj, k = 2, trymax = 200, trace = FALSE)
  sc <- vegan::scores(fit, display = "sites") |>
    as.data.frame() |>
    tibble::rownames_to_column("unit")
  # unify join key to 'unit'
  if ("coral_id" %in% names(meta_df)) {
    meta_use <- meta_df |> dplyr::select(unit = coral_id, treatment, reef)
  } else {
    meta_use <- meta_df |> dplyr::select(unit = reef, treatment)
  }
  sc <- dplyr::left_join(sc, meta_use, by = "unit") |>
    dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)))
  hulls <- hull_by_group16(sc, grp = "treatment", x = "NMDS1", y = "NMDS2")
  ggplot(sc, aes(NMDS1, NMDS2, color = treatment)) +
    (if (nrow(hulls) > 0)
      geom_polygon(
        data = hulls,
        aes(x = NMDS1, y = NMDS2, group = treatment, fill = treatment),
        inherit.aes = FALSE, alpha = 0.12, color = NA
      ) else NULL) +
    geom_point(size = if ("reef" %in% names(meta_use)) 3 else 2, alpha = 0.95) +
    scale_color_manual(values = palette, name = "Treatment") +
    scale_fill_manual(values = palette, guide = "none") +
    coord_equal() +
    labs(title = sprintf("%s (stress = %.3f)", title_stub, fit$stress),
         x = "NMDS1", y = "NMDS2") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

# tidy extractor for adonis2 (row 1 = treatment)
adonis_row1_16 <- function(fit) {
  tab <- as.data.frame(fit)
  tibble::tibble(
    Df = suppressWarnings(as.numeric(tab$Df[1])),
    F  = suppressWarnings(as.numeric(tab$F[1])),
    R2 = suppressWarnings(as.numeric(tab$R2[1])),
    p  = suppressWarnings(as.numeric(tab[1, "Pr(>F)"]))
  )
}

# ---- 16.2 Matrices (Abundance vs Proportion; Coral vs Reef) -----------------
# Coral-level matrices (already 10×10-filtered species earlier)
X_coral_abund <- mat_coral_abund
X_coral_prop  <- row_props16(mat_coral_abund)

# Reef-level matrices (already aggregated in reef_abund$X)
X_reef_abund  <- reef_abund$X
X_reef_prop   <- row_props16(reef_abund$X)

# ---- 16.2b Reef-level "density" (per-colony) matrix -------------------------
# Divide each reef's species counts by its treatment (1, 3, or 6) to get per-colony density.
reef_trt_vec <- meta_reef %>%
  dplyr::select(reef, treatment) %>%
  dplyr::mutate(treatment = as.numeric(as.character(treatment))) %>%
  dplyr::slice(match(rownames(X_reef_abund), reef)) %>%
  dplyr::pull(treatment)
stopifnot(all(is.finite(reef_trt_vec)), length(reef_trt_vec) == nrow(X_reef_abund))
X_reef_density <- sweep(X_reef_abund, 1, reef_trt_vec, "/")

# ---- 16.3 Distances (Gower) -------------------------------------------------
dist_coral_gow_abund <- gower_from_abund16(X_coral_abund)
dist_coral_gow_prop  <- gower_from_prop16 (X_coral_prop)
dist_reef_gow_abund  <- gower_from_abund16(X_reef_abund)
dist_reef_gow_prop   <- gower_from_prop16 (X_reef_prop)
# Density uses the abundance flavor transform (sqrt -> z -> Gower)
dist_reef_gow_density <- gower_from_abund16(X_reef_density)

# ---- 16.4 Ordinations & panels (create all plots before any comparisons) ----
p_coral_abund <- plot_gower_nmds_hulls16(dist_coral_gow_abund, meta_coral, "Coral × Abundance (Gower)")
p_coral_prop  <- plot_gower_nmds_hulls16(dist_coral_gow_prop,  meta_coral, "Coral × Proportion (Gower)")
p_reef_abund  <- plot_gower_nmds_hulls16(dist_reef_gow_abund,  meta_reef,  "Reef × Abundance (Gower)")
p_reef_prop   <- plot_gower_nmds_hulls16(dist_reef_gow_prop,   meta_reef,  "Reef × Proportion (Gower)")
p_reef_density<- plot_gower_nmds_hulls16(dist_reef_gow_density,meta_reef,  "Reef × Density (Abund/1,3,6; Gower)")

# 2×2 and 2×3 panels
p_2x2_gower <- (p_coral_abund | p_coral_prop) / (p_reef_abund | p_reef_prop)
p_2x3_gower <- (p_coral_abund | p_coral_prop | p_reef_abund) /
  (p_reef_prop | p_reef_density | patchwork::plot_spacer()) +
  patchwork::plot_annotation(
    title = "Gower distance — Abundance vs Proportion × Coral vs Reef (+ Density)",
    subtitle = "Density = reef counts ÷ treatment size (1,3,6). Abundance: sqrt(counts)→z; Proportion: row-prop→z"
  )

show_and_save(p_2x2_gower, file.path(FIG_DIR, "16_gower_2x2_coral_reef_abund_prop.png"), width = 12, height = 10)
show_and_save(p_2x3_gower, file.path(FIG_DIR, "16_gower_2x3_with_density.png"),       width = 15, height = 10)

# Optional quick comparison figure: Abundance vs Density (reef only)
p_reef_abund_vs_density <- (p_reef_abund | p_reef_density) +
  patchwork::plot_annotation(
    title = "Reef-level: Abundance vs Per-colony Density (Gower, sqrt→z)",
    subtitle = "Density = reef counts divided by treatment size (1, 3, 6)"
  )
show_and_save(p_reef_abund_vs_density, file.path(FIG_DIR, "16_reef_abundance_vs_density.png"), width = 10, height = 5)

# ---- 16.5 PERMANOVAs for each combination -----------------------------------
# Coral-level: restrict permutations WITHIN reefs (reef as strata)
fit_coral_abund <- vegan::adonis2(dist_coral_gow_abund ~ treatment, data = meta_coral, permutations = 999, strata = meta_coral$reef)
fit_coral_prop  <- vegan::adonis2(dist_coral_gow_prop  ~ treatment, data = meta_coral, permutations = 999, strata = meta_coral$reef)

# Reef-level: reefs are the units (no strata)
fit_reef_abund  <- vegan::adonis2(dist_reef_gow_abund  ~ treatment, data = meta_reef,  permutations = 999)
fit_reef_prop   <- vegan::adonis2(dist_reef_gow_prop   ~ treatment, data = meta_reef,  permutations = 999)
fit_reef_density<- vegan::adonis2(dist_reef_gow_density~ treatment, data = meta_reef,  permutations = 999)

sum16 <- dplyr::bind_rows(
  adonis_row1_16(fit_coral_abund)  |> dplyr::mutate(scale = "Coral (strata=reef)", flavor = "Abundance"),
  adonis_row1_16(fit_coral_prop )  |> dplyr::mutate(scale = "Coral (strata=reef)", flavor = "Proportion"),
  adonis_row1_16(fit_reef_abund )  |> dplyr::mutate(scale = "Reef",                 flavor = "Abundance"),
  adonis_row1_16(fit_reef_prop  )  |> dplyr::mutate(scale = "Reef",                 flavor = "Proportion"),
  adonis_row1_16(fit_reef_density)|> dplyr::mutate(scale = "Reef",                 flavor = "Per-colony density")
) |>
  dplyr::relocate(scale, flavor, .before = Df) |>
  dplyr::mutate(sig = dplyr::case_when(
    p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ ""
  ))

print(sum16)
readr::write_csv(sum16, file.path(TAB_DIR, "16_permanova_summary.csv"))

# nice GT table
gt16 <- sum16 |>
  dplyr::rename(df = Df) |>
  gt::gt(groupname_col = "scale") |>
  gt::tab_header(title = "§16 — PERMANOVA summary: Abundance vs Proportion × Scale (+ Density)") |>
  gt::fmt_number(columns = c(df, F, R2, p), decimals = 3) |>
  gt::cols_label(flavor = "Flavor", df = "df", F = "F", R2 = "R²", p = "p", sig = "Sig.") |>
  gt::tab_style(gt::cell_text(weight = "bold"),
                locations = gt::cells_body(columns = sig, rows = sig %in% c("***","**","*")))
gt::gtsave(gt16, file.path(TAB_DIR, "16_permanova_summary.html"))

cli::cli_alert_success("§16 complete: panels saved; PERMANOVAs run and summaries written.")


# ---- Build objects used by §16 vertical 2-panel ----
# (scores_prop, hulls_prop, centroids_prop, top15)

# 1) Reef × proportion NMDS scores (Gower)
set.seed(42)
fit_reef_prop <- vegan::metaMDS(dist_reef_gow_prop, k = 2, trymax = 200, trace = FALSE)

scores_prop <- vegan::scores(fit_reef_prop, display = "sites") |>
  as.data.frame()

# Standardize column names to NMDS1/NMDS2
if (all(c("MDS1","MDS2") %in% names(scores_prop))) {
  scores_prop <- dplyr::rename(scores_prop, NMDS1 = MDS1, NMDS2 = MDS2)
} else if (!all(c("NMDS1","NMDS2") %in% names(scores_prop))) {
  ax <- names(scores_prop)[vapply(scores_prop, is.numeric, logical(1))]
  scores_prop <- dplyr::rename(scores_prop, NMDS1 = !!ax[1], NMDS2 = !!ax[2])
}

scores_prop <- scores_prop |>
  tibble::rownames_to_column("reef") |>
  dplyr::left_join(meta_reef |> dplyr::select(reef, treatment), by = "reef") |>
  dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)))

# Convex hulls and centroids for Panel A
hulls_prop <- scores_prop |>
  dplyr::group_by(treatment) |>
  dplyr::filter(dplyr::n() >= 3, is.finite(NMDS1), is.finite(NMDS2)) |>
  dplyr::slice(chull(NMDS1, NMDS2)) |>
  dplyr::ungroup()

centroids_prop <- scores_prop |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(cx = mean(NMDS1, na.rm = TRUE),
                   cy = mean(NMDS2, na.rm = TRUE),
                   .groups = "drop")

# 2) Top-15 Δ z-scored per-colony density (6 − 1) for Panel B
# NOTE: If you prefer sqrt→z for density, change the next line to: scale(sqrt(X_reef_density))
Z_den <- scale(X_reef_density)
Z_den[!is.finite(Z_den)] <- 0

den_long <- as.data.frame(Z_den) |>
  tibble::rownames_to_column("reef") |>
  dplyr::left_join(meta_reef |> dplyr::select(reef, treatment), by = "reef") |>
  tidyr::pivot_longer(-c(reef, treatment), names_to = "species", values_to = "z")

top15 <- den_long |>
  dplyr::filter(treatment %in% c(1,6)) |>
  dplyr::group_by(species, treatment) |>
  dplyr::summarise(m = mean(z, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(treatment = as.character(treatment)) |>
  tidyr::pivot_wider(names_from = treatment, values_from = m, names_prefix = "t") |>
  dplyr::mutate(diff = t6 - t1,
                adiff = abs(diff),
                dir = ifelse(diff > 0, "↑ 6 > 1", "↓ 6 < 1")) |>
  dplyr::arrange(dplyr::desc(adiff)) |>
  dplyr::slice(1:15)

# =========================
# §16.6/16.8 — Vertical 2-panel figure
# Panels:
#   A) Reef NMDS (proportion, Gower) with hulls + centroid stars
#   B) Top-15 Δ z-scored per-colony density (6 − 1) dumbbell
# =========================
cli::cli_h3("§16 — Vertical 2-panel: NMDS (prop) + Top-15 density change")

# ---- Safety: build centroids if missing ----
if (!exists("centroids_prop", inherits = FALSE)) {
  centroids_prop <- scores_prop %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarise(
      cx = mean(NMDS1, na.rm = TRUE),
      cy = mean(NMDS2, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------------------------------
# Panel A — NMDS (reef × proportion, Gower)
# -------------------------------------------
pA_nmds <- ggplot(scores_prop, aes(NMDS1, NMDS2, color = treatment)) +
  # convex hulls (skip if <3 reefs per level)
  { if (exists("hulls_prop", inherits=FALSE) && nrow(hulls_prop))
    geom_polygon(
      data = hulls_prop,
      aes(x = NMDS1, y = NMDS2, group = treatment, fill = treatment),
      inherit.aes = FALSE, alpha = 0.15, color = NA
    )
  } +
  # reef points (circles) + centroid stars
  geom_point(shape = 16, size = 3.5, alpha = 0.8, stroke = 0.5) +
  geom_point(
    data = centroids_prop,
    aes(x = cx, y = cy),
    inherit.aes = TRUE,
    shape = 8, size = 5, stroke = 1.2
  ) +
  scale_color_manual(
    values = cols_trt,
    name = "Coral number",
    labels = c("1", "3", "6")
  ) +
  scale_fill_manual(values = cols_trt, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(x = "NMDS1", y = "NMDS2") +
  # Add stress annotation in upper-left corner
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Stress = %.3f", fit_reef_prop$stress),
           hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "italic") +
  theme_pub() +
  theme(
    legend.position = c(0.38, 0.08),  # Bottom-center-left
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.5),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.55, "cm"),
    legend.margin = margin(4, 6, 4, 6),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    plot.title   = element_blank(),
    plot.subtitle= element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# -------------------------------------------
# Panel B — Top-15 Δ density (z; 6 − 1) dumbbell
# -------------------------------------------
# Symmetric x-limits (balanced around 0)
xmax <- max(abs(c(top15$t1, top15$t6)), na.rm = TRUE)
xlim_sym <- c(-xmax, xmax) * 1.05

# ---- Reorder species for Panel B (↓ 6<1 on top, then ↑ 6>1 on bottom) ----
# Build an ordering key:
#   dir_flag = 1 for "↓ 6 < 1" (negative, on top), 0 for "↑ 6 > 1" (positive, on bottom)
#   Negative group (↓): sort by magnitude ascending (least negative first)
#   Positive group (↑): sort by magnitude descending (largest first)
ord_tbl <- top15 %>%
  dplyr::mutate(
    dir_flag = ifelse(dir == "↓ 6 < 1", 1L, 0L),  # SWAPPED: negative now gets 1
    sort_key = ifelse(dir == "↑ 6 > 1", -adiff, adiff)  # negative for desc, positive for asc
  ) %>%
  dplyr::arrange(dplyr::desc(dir_flag), sort_key)

# Factor levels determine vertical order in a horizontal dumbbell:
# last level is at the TOP, so reverse the vector to put ↓ block on top.
levs <- rev(ord_tbl$species)

top15 <- top15 %>%
  dplyr::mutate(species = factor(species, levels = levs))

# ---- Dumbbell plot (Panel B) -------------------------------------------------
# ---- Dumbbell plot (Panel B) as arrows --------------------------------------
pB_db <- ggplot(top15) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.5, linetype = "dashed") +
  # segments become single-headed arrows from t1 to t6
  geom_segment(
    aes(x = t1, xend = t6, y = species, yend = species, color = dir),
    linewidth = 1.2,
    arrow = arrow(type = "closed", length = unit(3.5, "mm"))
  ) +
  # endpoints - treatment 1 (orange) and treatment 6 (green)
  geom_point(aes(x = t1, y = species), shape = 21, size = 3.5,
             fill = cols_trt[["1"]], color = "black", stroke = 0.6) +
  geom_point(aes(x = t6, y = species), shape = 21, size = 3.5,
             fill = cols_trt[["6"]], color = "black", stroke = 0.6) +
  scale_color_manual(
    values = c("↑ 6 > 1" = "#CC79A7", "↓ 6 < 1" = "#0072B2"),  # Purple for positive, Teal for negative
    name   = NULL,
    guide  = "none"
  ) +
  labs(x = "Standardized per-colony density (z-score)", y = NULL) +
  coord_cartesian(xlim = xlim_sym) +
  theme_pub() +
  theme(
    legend.position   = "none",
    axis.title.x      = element_text(size = 14, face = "bold"),
    axis.text.x       = element_text(size = 11),
    axis.text.y       = element_text(face = "italic", size = 11),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.9),
    plot.title        = element_blank(),
    plot.subtitle     = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# -------------------------------------------
# Arrange side-by-side, tag panels A/B, save
# -------------------------------------------

#THIS IS THE A PUBLICATION FIGURE

p_horizontal <- pA_nmds + pB_db +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(widths = c(1, 1.3)) &
  theme(
    plot.tag = element_text(face = "bold", size = 18, hjust = 1, vjust = 1),
    plot.tag.position = c(0.98, 0.98)
  )

show_and_save(
  p_horizontal,
  file.path(FIG_DIR, "16_horizontal_2panel_nmds_prop_plus_top15_density.png"),
  width = 12, height = 5.5, dpi = 600
)

# Also save to publication-figures folder (PDF and PNG)
ggsave(
  file.path(OUT_DIR, "figures", "publication-figures", "16_horizontal_2panel_nmds_prop_plus_top15_density.pdf"),
  p_horizontal, width = 12, height = 5.5, dpi = 600, bg = "white"
)
ggsave(
  file.path(OUT_DIR, "figures", "publication-figures", "16_horizontal_2panel_nmds_prop_plus_top15_density.png"),
  p_horizontal, width = 12, height = 5.5, dpi = 600, bg = "white"
)

cli::cli_alert_success("Saved: 16_horizontal_2panel_nmds_prop_plus_top15_density.png")










# --- NMDS (Gower on row-wise proportions) — standalone -----------------------
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(vegan); library(cluster)

# If you start from a long table `df` with columns reef, species, count:
# build_reef_matrix <- function(df) {
#   df %>% group_by(reef, species) %>%
#     summarise(val = sum(count, na.rm = TRUE), .groups = "drop") %>%
#     pivot_wider(names_from = species, values_from = val, values_fill = 0) %>%
#     column_to_rownames("reef") %>% as.matrix()
# }

# --- NMDS (Gower on per-colony density) — stars colored by treatment ---------
library(dplyr); library(tidyr); library(tibble)
library(ggplot2); library(vegan); library(cluster)

safe_z <- function(X) {
  X <- as.matrix(X); sds <- apply(X, 2, sd, na.rm = TRUE)
  Z <- scale(X); Z[, !is.finite(sds) | sds == 0] <- 0; Z
}
hulls_by <- function(df, grp = "treatment") {
  df %>% group_by(.data[[grp]]) %>%
    filter(n() >= 3, is.finite(NMDS1), is.finite(NMDS2)) %>%
    slice(chull(NMDS1, NMDS2)) %>% ungroup()
}

plot_reef_nmds_gower_density <- function(
    X_reef_abund, meta_reef,
    cols_trt = TREATMENT_COLORS
) {
  stopifnot(is.matrix(X_reef_abund), all(rownames(X_reef_abund) %in% meta_reef$reef))
  
  # align meta and get numeric treatment sizes (1,3,6)
  meta <- meta_reef %>% mutate(treatment = factor(treatment, levels = c(1,3,6)))
  trt_num <- meta %>%
    slice(match(rownames(X_reef_abund), reef)) %>%
    pull(treatment) %>% as.character() %>% as.numeric()
  
  # per-colony density and Gower on sqrt→z
  X_density <- sweep(X_reef_abund, 1, trt_num, "/")
  Z <- safe_z(sqrt(X_density))
  D <- as.dist(cluster::daisy(Z, metric = "gower"))
  
  set.seed(42)
  ord <- vegan::metaMDS(D, k = 2, trymax = 200, trace = FALSE)
  
  sc <- vegan::scores(ord, display = "sites") %>% as.data.frame() %>%
    tibble::rownames_to_column("reef") %>% left_join(meta, by = "reef")
  names(sc)[names(sc) %in% c("MDS1","NMDS1")] <- "NMDS1"
  names(sc)[names(sc) %in% c("MDS2","NMDS2")] <- "NMDS2"
  sc <- sc %>% mutate(treatment = factor(treatment, levels = c(1,3,6)))
  
  cent  <- sc %>% group_by(treatment) %>%
    summarise(cx = mean(NMDS1), cy = mean(NMDS2), .groups = "drop")
  hulls <- hulls_by(sc, "treatment")
  
  p <- ggplot(sc, aes(NMDS1, NMDS2, color = treatment)) +
    (if (nrow(hulls)) geom_polygon(
      data = hulls,
      aes(x = NMDS1, y = NMDS2, group = treatment, fill = treatment),
      alpha = 0.12, color = NA, inherit.aes = FALSE)) +
    geom_point(size = 3, alpha = 0.95) +
    # ★ Stars explicitly mapped to treatment color
    geom_point(
      data = cent,
      aes(x = cx, y = cy, color = treatment),
      inherit.aes = FALSE,
      shape = 8, size = 5, stroke = 1.05
    ) +
    scale_color_manual(values = cols_trt, name = "Coral number") +
    scale_fill_manual(values = cols_trt, guide = "none") +
    coord_equal() +
    labs(x = "NMDS1", y = "NMDS2")+
         # subtitle = sprintf("Gower on sqrt(per-colony density) → z (stress = %.3f)", ord$stress)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
  
  print(p)  # ensure it draws now
  invisible(list(plot = p, stress = ord$stress, scores = sc, X_density = X_density))
}

# Example:
# res <- plot_reef_nmds_gower_density(X_reef_abund, meta_reef)
# ggsave(file.path(FIG_DIR, "nmds_gower_density_stars.png"),
#        res$plot, width = 6, height = 6, dpi = 600, bg = "white")

res <- plot_reef_nmds_gower_density(X_reef_abund, meta_reef)
ggsave(file.path(FIG_DIR, "nmds_gower_density_square.png"),
res$plot, width = 6, height = 6, dpi = 600, bg = "white")


# =============================================================================
# §17) Resampling strategies: Rarefaction, Bootstrap, N-sensitivity
# -----------------------------------------------------------------------------
# Purpose:
#   - Test robustness of PERMANOVA treatment effects to sample-size imbalance.
#   - Three strategies:
#       A) Balanced rarefaction (downsample to n_min per treatment)
#       B) Balanced bootstrap (upsample to n_max per treatment)
#       C) N-sensitivity (run across N = 2..n_min per treatment)
#   - Compare each strategy to main unbalanced §15 PERMANOVA.
# Outputs:
#   - Iteration tables:   res_rare, res_boot, res_sens
#   - Summaries:          sum_rare, sum_boot, sens_summary, sum_sens_overall
#   - Concordance:        concordance vs main (§15) + GT table
#   - Figures:            p-value hists, F/R² curves & bars
# =============================================================================

cli::cli_h1("§17 — Resampling strategies")

# ---- 17.0 Guards & defaults --------------------------------------------------
.required17 <- c("meta_reef","TAB_DIR","FIG_DIR","theme_pub","show_and_save")
if (!all(vapply(.required17, exists, logical(1)))) {
  stop("§17 requires: meta_reef, TAB_DIR, FIG_DIR, theme_pub(), show_and_save(). Run earlier blocks first.")
}

# Distance set from §15; build if missing
if (!exists("dist_full", inherits = FALSE)) {
  if (!exists("dist_reef_jac", inherits = FALSE) && !exists("dist_reef_gow", inherits = FALSE))
    stop("Need dist_full (preferred) or dist_reef_jac/dist_reef_gow present from §15.")
  dist_full <- list()
  if (exists("dist_reef_jac", inherits = FALSE)) dist_full$jac <- dist_reef_jac
  if (exists("dist_reef_gow", inherits = FALSE)) dist_full$sz  <- dist_reef_gow
  if (exists("dist_reef_gow_prop", inherits = FALSE)) dist_full$pz <- dist_reef_gow_prop
}

# Tidy treatment factor & ids_by_trt
meta_reef <- meta_reef %>% dplyr::mutate(treatment = factor(treatment, levels = c(1,3,6)),
                                         reef = as.character(reef))
ids_by_trt <- split(meta_reef$reef, meta_reef$treatment, drop = TRUE)

# Permutations
if (!exists("N_PERM", inherits = FALSE)) N_PERM <- 999L

# Iterations
if (!exists("B_RARE", inherits = FALSE)) B_RARE <- 199L
if (!exists("B_BOOT", inherits = FALSE)) B_BOOT <- 199L
if (!exists("B_SENS", inherits = FALSE)) B_SENS <- 199L
set.seed(42)

# Small helpers
`%||%` <- function(a,b) if (!is.null(a)) a else b
stars_  <- function(p) dplyr::case_when(is.na(p) ~ "", p < .001 ~ "***", p < .01 ~ "**",
                                        p < .05 ~ "*", p < .1 ~ ".", TRUE ~ "")
.canon_metric <- function(x) {
  lx <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    grepl("^jac|jaccard|incid", lx)                                  ~ "Jaccard (incidence)",
    grepl("^sz$|sqrt.?->?z|sqrt.?z|abund|gower.*sqrt", lx)            ~ "Gower (sqrt→z)",
    grepl("^gow$|^gower$", lx) & !grepl("prop|proportion", lx)        ~ "Gower (sqrt→z)",
    grepl("^pz$|prop.?->?z|proportion|row-?prop|gower.*prop", lx)     ~ "Gower (prop→z)",
    grepl("gower", lx) & grepl("prop", lx)                            ~ "Gower (prop→z)",
    grepl("gower", lx)                                                ~ "Gower (sqrt→z)",
    TRUE ~ x
  )
}

# Normalise names on dist_full (important for joins later)
if (is.null(names(dist_full))) names(dist_full) <- paste0("metric_", seq_along(dist_full))
names(dist_full) <- .canon_metric(names(dist_full))

# Extract row 1 from adonis2
.adonis_row1 <- function(fit) {
  tab <- as.data.frame(fit)
  tibble::tibble(
    df = suppressWarnings(as.numeric(tab$Df[1])),
    F  = suppressWarnings(as.numeric(tab$F[1])),
    R2 = suppressWarnings(as.numeric(tab$R2[1])),
    p  = suppressWarnings(as.numeric(tab[1, "Pr(>F)"]))
  )
}

# Quiet alignment + adonis2
.adonis_aligned_quiet <- function(D, meta, nperm = 999) {
  stopifnot(inherits(D, "dist"))
  labs_all <- attr(D, "Labels"); if (is.null(labs_all)) stop("Distance object has no 'Labels'.")
  meta2 <- meta %>% dplyr::mutate(reef = as.character(reef),
                                  treatment = factor(treatment, levels = c(1,3,6)))
  keep <- labs_all %in% meta2$reef
  labs <- labs_all[keep]
  if (length(labs) < 3L) return(NULL)
  M <- as.matrix(D)
  Dsub  <- stats::as.dist(M[keep, keep, drop = FALSE])
  metaA <- meta2[match(labs, meta2$reef), , drop = FALSE]
  if (nlevels(droplevels(metaA$treatment)) < 2) return(NULL)
  suppressMessages(vegan::adonis2(Dsub ~ treatment, data = metaA, permutations = nperm))
}

# Align meta to a distance (optionally subset to a pick)
.align_meta_to_dist17 <- function(dist_obj, meta, pick = NULL) {
  labs <- attr(dist_obj, "Labels"); if (is.null(labs)) stop("Distance object has no 'Labels'.")
  if (!is.null(pick)) labs <- intersect(labs, pick)
  if (length(labs) < 3L) return(NULL)
  meta2 <- meta %>% dplyr::mutate(reef = as.character(.data$reef))
  idx   <- match(labs, meta2$reef)
  if (any(is.na(idx))) return(NULL)
  meta2[idx, , drop = FALSE]
}

# Worker: run PERMANOVA (across all metrics) for one draw of reef ids
.run_perms_one <- function(pick, label = "Draw", iter_idx = NA_integer_, verbose = FALSE) {
  outs <- lapply(names(dist_full), function(metric_name) {
    dist_obj <- dist_full[[metric_name]]
    labs_all <- attr(dist_obj, "Labels")
    if (is.null(labs_all)) return(NULL)
    keep <- intersect(labs_all, pick)
    if (length(keep) < 3L) return(NULL)
    Dm   <- as.matrix(dist_obj)
    idx  <- match(keep, labs_all)
    Dsub <- as.dist(Dm[idx, idx, drop = FALSE])
    meta_sub <- .align_meta_to_dist17(dist_obj, meta_reef, pick = keep)
    if (is.null(meta_sub)) return(NULL)
    g <- droplevels(meta_sub$treatment)
    if (nlevels(g) < 2L) return(NULL)
    fit <- tryCatch(
      suppressMessages(vegan::adonis2(Dsub ~ treatment, data = meta_sub, permutations = N_PERM)),
      error = function(e) { if (verbose) message("adonis2 error: ", conditionMessage(e)); NULL }
    )
    if (is.null(fit)) return(NULL)
    r1 <- .adonis_row1(fit)
    tibble::tibble(
      iter      = iter_idx,
      label     = label,
      metric    = .canon_metric(metric_name),
      df        = r1$df, F = r1$F, R2 = r1$R2, p = r1$p,
      n_units   = length(keep),
      n_trt_lvls= nlevels(g)
    )
  })
  dplyr::bind_rows(outs)
}

# ---- 17.1 Design info --------------------------------------------------------
n_per_trt <- vapply(ids_by_trt, length, integer(1))
n_min <- min(n_per_trt); n_max <- max(n_per_trt)
cli::cli_alert_info("Treatments: min N = {n_min}, max N = {n_max}")

# ---- 17.2 Rarefaction (balanced downsample) ---------------------------------
cli::cli_h3("17A) Balanced rarefaction")
if (n_min < 2) {
  cli::cli_alert_warning("Skipping rarefaction: n_min < 2")
  res_rare <- tibble::tibble()
} else {
  rare_list <- vector("list", B_RARE); n_ok <- 0L
  for (i in seq_len(B_RARE)) {
    pick <- tryCatch(unlist(Map(\(v) sample(v, n_min, replace = FALSE), ids_by_trt),
                            use.names = FALSE),
                     error = function(e) { cli::cli_alert_warning("Rare iter {i}: {conditionMessage(e)}"); NULL })
    if (is.null(pick)) next
    out_i <- tryCatch(.run_perms_one(pick, label = "Rarefied", iter_idx = i),
                      error = function(e) { cli::cli_alert_warning(".run_perms_one() failed (rare {i}): {conditionMessage(e)}"); NULL })
    if (!is.null(out_i) && nrow(out_i)) { n_ok <- n_ok + 1L; rare_list[[n_ok]] <- out_i }
    if (i %% 25 == 0) cli::cli_alert_info("… {i}/{B_RARE} rarefaction iterations")
  }
  res_rare <- dplyr::bind_rows(rare_list[seq_len(n_ok)])
  if (nrow(res_rare)) {
    readr::write_csv(res_rare, file.path(TAB_DIR, "section17_results_rarefaction.csv"))
    cli::cli_alert_success("Rarefaction: {nrow(res_rare)} rows across {n_ok} iterations.")
  } else cli::cli_alert_warning("Rarefaction produced no valid iterations.")
}

# ---- 17.3 Bootstrap (balanced upsample) -------------------------------------
cli::cli_h3("17B) Balanced bootstrap")
.draw_boot_ids <- function(ids_by_trt, n_max, max_tries = 50) {
  for (t in seq_len(max_tries)) {
    pick_by_trt <- lapply(ids_by_trt, function(v) sample(v, size = n_max, replace = TRUE))
    ok <- all(vapply(pick_by_trt, function(v) length(unique(v)) >= 2, logical(1)))
    if (ok) return(unlist(pick_by_trt, use.names = FALSE))
  }
  stop("Could not satisfy >=2 unique reefs per treatment within max_tries.")
}
if (n_max < 2) {
  cli::cli_alert_warning("Skipping bootstrap: n_max < 2")
  res_boot <- tibble::tibble()
} else {
  boot_list <- vector("list", B_BOOT); n_ok <- 0L
  for (i in seq_len(B_BOOT)) {
    pick <- tryCatch(.draw_boot_ids(ids_by_trt, n_max),
                     error = function(e) { cli::cli_alert_warning("Boot draw {i} failed: {conditionMessage(e)}"); NULL })
    if (is.null(pick)) next
    out_i <- tryCatch(.run_perms_one(pick, label = "Bootstrap", iter_idx = i),
                      error = function(e) { cli::cli_alert_warning(".run_perms_one() failed (boot {i}): {conditionMessage(e)}"); NULL })
    if (!is.null(out_i) && nrow(out_i)) { n_ok <- n_ok + 1L; boot_list[[n_ok]] <- out_i }
    if (i %% 25 == 0) cli::cli_alert_info("… {i}/{B_BOOT} bootstrap iterations")
  }
  res_boot <- dplyr::bind_rows(boot_list[seq_len(n_ok)])
  if (nrow(res_boot)) {
    readr::write_csv(res_boot, file.path(TAB_DIR, "section17_results_bootstrap.csv"))
    cli::cli_alert_success("Bootstrap: {nrow(res_boot)} rows across {n_ok} iterations.")
  } else cli::cli_alert_warning("Bootstrap produced no valid iterations.")
}

# ---- 17.4 N-sensitivity (N = 2..n_min) --------------------------------------
cli::cli_h3("17C) N-sensitivity")
if (n_min < 2) {
  cli::cli_alert_warning("Skipping N-sensitivity: n_min < 2")
  res_sens <- tibble::tibble(); sens_summary <- tibble::tibble()
} else {
  sens_list <- list(); row_ix <- 0L
  for (n_each in 2:n_min) {
    cli::cli_alert_info("§17C running N = {n_each} (B_SENS = {B_SENS})")
    for (i in seq_len(B_SENS)) {
      pick <- tryCatch(unlist(lapply(ids_by_trt, function(v) sample(v, n_each, replace = FALSE)), use.names = FALSE),
                       error = function(e) { cli::cli_alert_warning("N={n_each}, iter {i} draw failed: {conditionMessage(e)}"); NULL })
      if (is.null(pick)) next
      out_i <- tryCatch(.run_perms_one(pick, label = sprintf("N-sens (n=%d/trt)", n_each), iter_idx = i),
                        error = function(e) { cli::cli_alert_warning(".run_perms_one() failed (N={n_each}, iter {i}): {conditionMessage(e)}"); NULL })
      if (!is.null(out_i) && nrow(out_i)) { row_ix <- row_ix + 1L; sens_list[[row_ix]] <- dplyr::mutate(out_i, n_each = n_each) }
    }
  }
  res_sens <- dplyr::bind_rows(sens_list)
  if (nrow(res_sens)) {
    readr::write_csv(res_sens, file.path(TAB_DIR, "section17_results_sensitivity_raw.csv"))
    cli::cli_alert_success("N-sensitivity: {nrow(res_sens)} rows total.")
    sens_summary <- res_sens %>%
      dplyr::group_by(metric, n_each) %>%
      dplyr::summarise(
        mean_F  = mean(F,  na.rm = TRUE),
        q25_F   = stats::quantile(F,  0.25, na.rm = TRUE),
        q75_F   = stats::quantile(F,  0.75, na.rm = TRUE),
        mean_R2 = mean(R2, na.rm = TRUE),
        q25_R2  = stats::quantile(R2, 0.25, na.rm = TRUE),
        q75_R2  = stats::quantile(R2, 0.75, na.rm = TRUE),
        med_p   = stats::median(p, na.rm = TRUE),
        prop_sig= mean(p < 0.05, na.rm = TRUE),
        .groups = "drop"
      )
    readr::write_csv(sens_summary, file.path(TAB_DIR, "section17_results_sensitivity_summary.csv"))
  } else {
    cli::cli_alert_warning("N-sensitivity produced no valid iterations.")
    sens_summary <- tibble::tibble()
  }
}

# ---- 17.5 Summaries & concordance vs main -----------------------------------
cli::cli_h3("17D) Concordance vs §15 main PERMANOVA")

# Summarise res_rare / res_boot
.summarise_iter_tbl <- function(tbl, label) {
  if (is.null(tbl) || !is.data.frame(tbl) || !nrow(tbl)) return(tibble::tibble())
  stopifnot(all(c("metric","F","R2","p") %in% names(tbl)))
  tbl %>%
    dplyr::mutate(metric = .canon_metric(metric)) %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(
      technique = label,
      F_med     = stats::median(F,  na.rm = TRUE),
      R2_med    = stats::median(R2, na.rm = TRUE),
      p_med     = stats::median(p,  na.rm = TRUE),
      prop_p_lt_0.05 = mean(p < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
}
sum_rare <- .summarise_iter_tbl(res_rare, "Rarefaction")
sum_boot <- .summarise_iter_tbl(res_boot, "Bootstrap")

# N-sensitivity: per-N and across-N rollup
if (nrow(sens_summary)) {
  sens_by_n <- sens_summary %>%
    dplyr::transmute(
      metric = .canon_metric(metric), n_each,
      F_med = mean_F, R2_med = mean_R2, p_med = med_p, prop_p_lt_0.05 = prop_sig
    )
  sum_sens_overall <- sens_by_n %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(
      technique = "N-sensitivity (across N)",
      F_med     = stats::median(F_med,  na.rm = TRUE),
      R2_med    = stats::median(R2_med, na.rm = TRUE),
      p_med     = stats::median(p_med,  na.rm = TRUE),
      prop_p_lt_0.05 = stats::median(prop_p_lt_0.05, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(sens_by_n,        file.path(TAB_DIR, "section17_sensitivity_byN.csv"))
  readr::write_csv(sum_sens_overall, file.path(TAB_DIR, "section17_summary_sensitivity_overall.csv"))
} else {
  sens_by_n <- tibble::tibble()
  sum_sens_overall <- tibble::tibble()
}

# Main unbalanced models (rebuild from dist_full)
main_full <- purrr::imap_dfr(dist_full, function(D, nm) {
  fit <- .adonis_aligned_quiet(D, meta_reef, nperm = N_PERM)
  if (is.null(fit)) return(tibble::tibble())
  .adonis_row1(fit) %>%
    dplyr::mutate(metric = .canon_metric(nm))
}) %>%
  dplyr::rename(F_main = F, R2_main = R2, p_main = p)

if (!nrow(main_full)) cli::cli_alert_warning("Main PERMANOVAs did not run (too few units/levels?).")

# Concordance table
sum_all <- dplyr::bind_rows(
  sum_rare %>% dplyr::select(metric, technique, F_med, R2_med, p_med, prop_p_lt_0.05),
  sum_boot %>% dplyr::select(metric, technique, F_med, R2_med, p_med, prop_p_lt_0.05),
  sum_sens_overall %>% dplyr::select(metric, technique, F_med, R2_med, p_med, prop_p_lt_0.05)
) %>% dplyr::mutate(metric = .canon_metric(metric))

concordance <- sum_all %>%
  dplyr::left_join(main_full, by = "metric") %>%
  dplyr::mutate(
    stars_med  = stars_(p_med),
    stars_main = stars_(p_main)
  ) %>%
  dplyr::relocate(metric, technique, F_med, R2_med, p_med, stars_med,
                  F_main, R2_main, p_main, stars_main, prop_p_lt_0.05)

readr::write_csv(concordance, file.path(TAB_DIR, "section17_concordance_vs_main.csv"))
cli::cli_alert_success("Concordance written: {nrow(concordance)} rows.")

# GT table
if (nrow(concordance)) {
  gt_conc <- concordance %>%
    gt::gt(groupname_col = "metric") %>%
    gt::tab_header(
      title   = "§17 — Concordance of resampling medians vs §15 main PERMANOVA",
      subtitle= "Medians across iterations; proportion significant at α = 0.05"
    ) %>%
    gt::fmt_number(columns = c(F_med, R2_med, p_med, F_main, R2_main, p_main, prop_p_lt_0.05), decimals = 3) %>%
    gt::cols_label(
      technique = "Method",
      F_med = "F (med)", R2_med = "R² (med)", p_med = "p (med)", stars_med = "Sig.",
      F_main = "F (main)", R2_main = "R² (main)", p_main = "p (main)", stars_main = "Sig. (main)",
      prop_p_lt_0.05 = "Prop. p<0.05"
    ) %>%
    gt::tab_style(gt::cell_text(weight = "bold"),
                  locations = list(
                    gt::cells_body(columns = stars_med,  rows = stars_med  %in% c("***","**","*")),
                    gt::cells_body(columns = stars_main, rows = stars_main %in% c("***","**","*"))
                  ))
  gt::gtsave(gt_conc, file.path(TAB_DIR, "section17_concordance_vs_main.html"))
}

# ---- 17.6 Figures ------------------------------------------------------------
cli::cli_h3("17.6 Figures")

# A) p-value histograms across techniques × metric
dist_df <- dplyr::bind_rows(
  if (exists("res_rare") && nrow(res_rare)) res_rare %>% dplyr::mutate(tech = "Rarefied"),
  if (exists("res_boot") && nrow(res_boot)) res_boot %>% dplyr::mutate(tech = "Bootstrap"),
  if (exists("res_sens") && nrow(res_sens)) res_sens %>% dplyr::mutate(tech = paste0("N-sens (n=", n_each, ")"))
)
if (nrow(dist_df)) {
  p_hist <- ggplot2::ggplot(dist_df, ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_grid(metric ~ tech, scales = "free_y") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed") +
    ggplot2::labs(title = "§17 — PERMANOVA p-value distributions", x = "p-value", y = "Count") +
    theme_pub()
  show_and_save(p_hist, file.path(FIG_DIR, "section17_hist_pvalues.png"), 13, 8)
}

# B) N-sensitivity curves (F with IQR ribbons)
if (nrow(sens_summary)) {
  p_sens <- ggplot2::ggplot(sens_summary, ggplot2::aes(n_each, mean_F)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q25_F, ymax = q75_F), alpha = 0.2, fill = "grey60") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(title = "§17C — N-sensitivity: PERMANOVA F vs N (IQR ribbons)",
                  x = "N per treatment", y = "PERMANOVA F") +
    theme_pub()
  show_and_save(p_sens, file.path(FIG_DIR, "section17_sensitivity_curves.png"), 9, 7)
}

# C) Bars: median p and prop. sig by technique × metric
bars_df <- dplyr::bind_rows(
  if (exists("sum_rare") && nrow(sum_rare)) sum_rare  %>% dplyr::transmute(metric, technique, p_med, prop_p_lt_0.05),
  if (exists("sum_boot") && nrow(sum_boot)) sum_boot  %>% dplyr::transmute(metric, technique, p_med, prop_p_lt_0.05),
  if (exists("sum_sens_overall") && nrow(sum_sens_overall)) sum_sens_overall %>% dplyr::transmute(metric, technique, p_med, prop_p_lt_0.05)
)
if (nrow(bars_df)) {
  p_bars <- bars_df %>%
    tidyr::pivot_longer(c(p_med, prop_p_lt_0.05), names_to = "what", values_to = "val") %>%
    dplyr::mutate(what = dplyr::recode(what, p_med = "Median p-value", prop_p_lt_0.05 = "Prop. p<0.05")) %>%
    ggplot2::ggplot(ggplot2::aes(technique, val, fill = metric)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::facet_wrap(~ what, scales = "free_y") +
    ggplot2::labs(title = "§17 — Resampling summaries by technique and distance",
                  x = NULL, y = NULL, fill = "Distance") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  show_and_save(p_bars, file.path(FIG_DIR, "section17_bars_technique_metric.png"), 11, 7)
}

cli::cli_alert_success("§17 complete: resampling results, summaries, concordance, and figures written.")