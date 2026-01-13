# ==============================================================================
# TITLE: Generate Publication Figures for Stier et al. CAFI 2025
# PURPOSE: Create all main text and supplemental figures from manuscript
#
# This script generates:
# - Figure 1: Community abundance and richness by coral number
# - Figure 2: Species-specific responses to coral number
# - Figure 3: Community composition (NMDS) and density changes
# - Figure 4: PCA of coral performance metrics
# - Figure 5: CAFI community composition vs coral performance
# - Figure 6: Individual species effects on coral performance
# - Supplemental Figure S1: Individual coral performance metrics
# - Supplemental Table S1: Treatment effects on performance
# - Supplemental Table S2: Species-performance correlations
# ==============================================================================

# Setup ------------------------------------------------------------------------
set.seed(1234)
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
  library(vegan)
  library(scales)
  library(gt)
})

# Source centralized figure standards (provides TREATMENT_COLORS, theme_publication, etc.)
source("scripts/MRB/mrb_figure_standards.R")

# Publication directories (follows standard output/MRB/figures/ structure)
PUB_DIR <- here("output", "MRB", "figures", "publication")
SUP_DIR <- here("output", "MRB", "figures", "publication", "supplemental")
dir.create(PUB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

# Use standard treatment colors from mrb_figure_standards.R
cols_trt <- TREATMENT_COLORS

# Theme for publication
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA),
      axis.ticks.length = unit(2.5, "mm"),
      axis.text = element_text(colour = "black"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

cat("\n=======================================================================\n")
cat("  GENERATING PUBLICATION FIGURES FOR STIER ET AL. CAFI 2025\n")
cat("=======================================================================\n\n")

# Load required data -----------------------------------------------------------
cat("Loading data...\n")

# Load CAFI data
cafi_data <- read_csv(here("data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv"),
                      show_col_types = FALSE) %>%
  filter(!is.na(coral_id), !is.na(species), !is.na(count))

# Load treatment/metadata
metadata <- read_csv(here("data/MRB Amount/coral_id_position_treatment.csv"),
                     show_col_types = FALSE) %>%
  rename(reef = coral_id)

# Load coral performance data (from mesh measurements and physiology)
coral_2019 <- read_csv(here("data/MRB Amount/MRB_2019_200K_mesh_measure.csv"),
                       show_col_types = FALSE)
coral_2021 <- read_csv(here("data/MRB Amount/MRB_May_2021_200K_mesh_measure.csv"),
                       show_col_types = FALSE)

cat("✓ Data loaded successfully\n\n")

# FIGURE 1: Community abundance and richness ----------------------------------
cat("Generating Figure 1: Community responses to coral number...\n")

# Calculate community metrics by reef
community_metrics <- cafi_data %>%
  group_by(reef = coral_id) %>%
  summarise(
    total_abundance = sum(count, na.rm = TRUE),
    richness = n_distinct(species[count > 0]),
    .groups = "drop"
  ) %>%
  left_join(metadata, by = "reef") %>%
  mutate(treatment = factor(treatment, levels = c(1, 3, 6)))

# Calculate rarefied richness (using vegan)
# This requires a community matrix
comm_matrix <- cafi_data %>%
  filter(!is.na(coral_id), !is.na(species)) %>%
  group_by(coral_id, species) %>%
  summarise(abundance = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = abundance,
              values_fill = 0, id_cols = coral_id) %>%
  column_to_rownames("coral_id") %>%
  as.matrix()

rarefied <- rarefy(comm_matrix, sample = min(rowSums(comm_matrix)))

community_metrics$rarefied_richness <- rarefied[match(community_metrics$reef, names(rarefied))]

# Bootstrap for expected values (Field of Dreams hypothesis)
set.seed(42)
treatment_1 <- community_metrics %>% filter(treatment == 1)

bootstrap_expected <- function(n_corals, n_boot = 10000) {
  replicate(n_boot, {
    sample(treatment_1$total_abundance, size = n_corals, replace = TRUE) %>% sum()
  })
}

expected_3 <- bootstrap_expected(3)
expected_6 <- bootstrap_expected(6)

# Panel A: Total abundance
fig1a <- community_metrics %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(total_abundance),
    se = sd(total_abundance) / sqrt(n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = as.numeric(as.character(treatment)), y = mean)) +
  # Expected ribbon for Field of Dreams
  geom_ribbon(data = tibble(
    treatment = c(3, 6),
    lower = c(quantile(expected_3, 0.025), quantile(expected_6, 0.025)),
    upper = c(quantile(expected_3, 0.975), quantile(expected_6, 0.975))
  ), aes(ymin = lower, ymax = upper, x = treatment),
  fill = "salmon", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "steelblue") +
  geom_line(data = tibble(treatment = c(3, 6),
                          mean = c(median(expected_3), median(expected_6))),
            linetype = "dashed", color = "salmon") +
  scale_x_continuous(breaks = c(1, 3, 6)) +
  labs(x = NULL, y = "Total abundance", tag = "a") +
  theme_pub() +
  theme(plot.tag = element_text(face = "bold", size = 14))

# Panel B: Species richness
fig1b <- community_metrics %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(richness),
    se = sd(richness) / sqrt(n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = as.numeric(as.character(treatment)), y = mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "steelblue") +
  scale_x_continuous(breaks = c(1, 3, 6)) +
  labs(x = NULL, y = "Species richness", tag = "b") +
  theme_pub() +
  theme(plot.tag = element_text(face = "bold", size = 14))

# Panel C: Rarefied richness
fig1c <- community_metrics %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(rarefied_richness, na.rm = TRUE),
    se = sd(rarefied_richness, na.rm = TRUE) / sqrt(n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = as.numeric(as.character(treatment)), y = mean)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "steelblue") +
  scale_x_continuous(breaks = c(1, 3, 6)) +
  labs(x = "Number of corals", y = "Rarefied richness", tag = "c") +
  theme_pub() +
  theme(plot.tag = element_text(face = "bold", size = 14))

# Combine Figure 1
figure1 <- fig1a / fig1b / fig1c
ggsave(file.path(PUB_DIR, "Figure1.png"), figure1,
       width = 6, height = 10, dpi = 600, bg = "white")
ggsave(file.path(PUB_DIR, "Figure1.pdf"), figure1,
       width = 6, height = 10, device = cairo_pdf, bg = "white")

cat("✓ Figure 1 saved\n\n")

# FIGURE 2: Species-specific responses -----------------------------------------
cat("Generating Figure 2: Species-specific abundance responses...\n")

# Prepare species-level data with treatment
species_by_reef <- cafi_data %>%
  group_by(coral_id, species, class) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  left_join(metadata %>% select(reef = reef, treatment),
            by = c("coral_id" = "reef")) %>%
  filter(!is.na(treatment)) %>%
  mutate(treatment = factor(treatment, levels = c(1, 3, 6)))

# Get top species by taxonomic group
get_top_species <- function(data, tax_class, n_species = 4) {
  data %>%
    filter(class == tax_class) %>%
    group_by(species) %>%
    summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>%
    slice_head(n = n_species) %>%
    pull(species)
}

top_fishes <- get_top_species(species_by_reef, "Actinopterygii", 4)
top_crustaceans <- get_top_species(species_by_reef, "Malacostraca", 4)
top_molluscs <- get_top_species(species_by_reef, "Gastropoda", 4)

all_focal_species <- c(top_fishes, top_crustaceans, top_molluscs)

# Bootstrap expected values for each species
treatment_1_species <- species_by_reef %>% filter(treatment == 1)

bootstrap_species <- function(species_name, n_corals, n_boot = 10000) {
  species_abund_trt1 <- treatment_1_species %>%
    filter(species == species_name) %>%
    pull(count)

  if(length(species_abund_trt1) == 0) return(rep(0, n_boot))

  replicate(n_boot, {
    sample(species_abund_trt1, size = n_corals, replace = TRUE) %>% sum()
  })
}

# Create panels for each species
fig2_panels <- lapply(all_focal_species, function(sp) {
  # Observed data
  obs_data <- species_by_reef %>%
    filter(species == sp) %>%
    group_by(treatment) %>%
    summarise(
      mean = mean(count, na.rm = TRUE),
      se = sd(count, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  # Expected values
  exp_3 <- bootstrap_species(sp, 3)
  exp_6 <- bootstrap_species(sp, 6)

  expected_data <- tibble(
    treatment = c(3, 6),
    exp_mean = c(median(exp_3), median(exp_6)),
    exp_lower = c(quantile(exp_3, 0.025), quantile(exp_6, 0.025)),
    exp_upper = c(quantile(exp_3, 0.975), quantile(exp_6, 0.975))
  )

  ggplot(obs_data, aes(x = as.numeric(as.character(treatment)), y = mean)) +
    geom_ribbon(data = expected_data,
                aes(ymin = exp_lower, ymax = exp_upper, x = treatment),
                fill = "salmon", alpha = 0.3, inherit.aes = FALSE) +
    geom_line(data = expected_data, aes(y = exp_mean, x = treatment),
              linetype = "dashed", color = "salmon") +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 2) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                  width = 0.2, color = "steelblue") +
    scale_x_continuous(breaks = c(1, 3, 6)) +
    labs(title = gsub("(.{1,30})(\\s|$)", "\\1\n", sp),
         x = NULL, y = "Abundance") +
    theme_pub() +
    theme(plot.title = element_text(size = 9, face = "italic"))
})

# Combine into 4x3 grid
figure2 <- wrap_plots(fig2_panels, ncol = 3)
ggsave(file.path(PUB_DIR, "Figure2.png"), figure2,
       width = 10, height = 13, dpi = 600, bg = "white")
ggsave(file.path(PUB_DIR, "Figure2.pdf"), figure2,
       width = 10, height = 13, device = cairo_pdf, bg = "white")

cat("✓ Figure 2 saved\n\n")

# FIGURE 3: NMDS and species density changes -----------------------------------
cat("Generating Figure 3: Community composition and density changes...\n")

# Source publication standards if not already loaded
if (!exists("theme_publication")) {
  source(here::here("scripts/MRB/mrb_figure_standards.R"))
}

# Load NMDS results if available
nmds_file <- here("output/objects_mrb/nmds_gower.rds")
if (file.exists(nmds_file)) {
  nmds_result <- readRDS(nmds_file)

  # Extract scores
  nmds_scores <- as.data.frame(scores(nmds_result)) %>%
    rownames_to_column("reef") %>%
    left_join(metadata, by = "reef") %>%
    mutate(treatment = factor(treatment, levels = c(1, 3, 6)))

  # Calculate centroids for each treatment
  centroids <- nmds_scores %>%
    group_by(treatment) %>%
    summarise(
      cent_NMDS1 = mean(NMDS1),
      cent_NMDS2 = mean(NMDS2),
      .groups = "drop"
    )

  # Panel A: NMDS ordination with publication standards
  fig3a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment)) +
    stat_ellipse(geom = "polygon", alpha = 0.25, level = 0.95, linewidth = 1.2) +
    geom_point(size = 5, alpha = 0.9, shape = 21, stroke = 1.5, color = "black") +
    # Add centroid stars
    geom_point(data = centroids,
               aes(x = cent_NMDS1, y = cent_NMDS2, color = treatment),
               size = 10, shape = 8, stroke = 2, show.legend = FALSE) +
    scale_color_manual(
      values = c("1" = "#E69F00", "3" = "#56B4E9", "6" = "#009E73"),
      labels = c("1", "3", "6"),
      name = "# Corals"
    ) +
    scale_fill_manual(
      values = c("1" = "#E69F00", "3" = "#56B4E9", "6" = "#009E73"),
      labels = c("1", "3", "6"),
      name = "# Corals"
    ) +
    labs(
      x = "NMDS1",
      y = "NMDS2",
      subtitle = paste0("Stress: ", round(nmds_result$stress, 3))
    ) +
    theme_ordination() +
    theme(
      legend.position = c(0.15, 0.9),
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = FONT_SIZE_LEGEND_TITLE),
      legend.text = element_text(size = FONT_SIZE_LEGEND_TEXT),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      plot.subtitle = element_text(hjust = 0.95, vjust = 0, size = FONT_SIZE_SUBTITLE),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.25)
    ) +
    coord_equal()
} else {
  # Create NMDS on the fly
  comm_matrix <- species_by_reef %>%
    group_by(coral_id, species) %>%
    summarise(abundance = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = species, values_from = abundance,
                values_fill = 0, id_cols = coral_id) %>%
    column_to_rownames("coral_id") %>%
    as.matrix()

  # Transform and standardize for Gower-like approach
  comm_trans <- sqrt(comm_matrix)
  comm_std <- decostand(comm_trans, method = "standardize")

  set.seed(1234)
  nmds_result <- metaMDS(comm_std, distance = "euclidean", k = 2, trymax = 200)

  nmds_scores <- as.data.frame(scores(nmds_result)) %>%
    rownames_to_column("reef") %>%
    left_join(metadata, by = "reef") %>%
    mutate(treatment = factor(treatment, levels = c(1, 3, 6)))

  # Calculate centroids for each treatment
  centroids <- nmds_scores %>%
    group_by(treatment) %>%
    summarise(
      cent_NMDS1 = mean(NMDS1),
      cent_NMDS2 = mean(NMDS2),
      .groups = "drop"
    )

  # Panel A: NMDS ordination with publication standards
  fig3a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment)) +
    stat_ellipse(geom = "polygon", alpha = 0.25, level = 0.95, linewidth = 1.2) +
    geom_point(size = 5, alpha = 0.9, shape = 21, stroke = 1.5, color = "black") +
    # Add centroid stars
    geom_point(data = centroids,
               aes(x = cent_NMDS1, y = cent_NMDS2, color = treatment),
               size = 10, shape = 8, stroke = 2, show.legend = FALSE) +
    scale_color_manual(
      values = c("1" = "#E69F00", "3" = "#56B4E9", "6" = "#009E73"),
      labels = c("1", "3", "6"),
      name = "# Corals"
    ) +
    scale_fill_manual(
      values = c("1" = "#E69F00", "3" = "#56B4E9", "6" = "#009E73"),
      labels = c("1", "3", "6"),
      name = "# Corals"
    ) +
    labs(
      x = "NMDS1",
      y = "NMDS2",
      subtitle = paste0("Stress: ", round(nmds_result$stress, 3))
    ) +
    theme_ordination() +
    theme(
      legend.position = c(0.15, 0.9),
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = FONT_SIZE_LEGEND_TITLE),
      legend.text = element_text(size = FONT_SIZE_LEGEND_TEXT),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      plot.subtitle = element_text(hjust = 0.95, vjust = 0, size = FONT_SIZE_SUBTITLE),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.25)
    ) +
    coord_equal()
}

# Panel B: Top 15 species density changes with publication standards
top_15_species <- species_by_reef %>%
  group_by(species) %>%
  summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(species)

density_changes <- species_by_reef %>%
  filter(species %in% top_15_species) %>%
  group_by(species, treatment) %>%
  summarise(mean_density = mean(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = mean_density, values_fill = 0) %>%
  mutate(
    change_6vs1 = `6` - `1`,
    species_short = gsub("(.{1,30})(\\s|$)", "\\1", species),  # Slightly longer names
    direction = ifelse(change_6vs1 > 0, "↑ 6 > 1", "↓ 6 < 1")  # Add direction indicators
  ) %>%
  arrange(change_6vs1)

# Create bar chart with publication formatting
fig3b <- ggplot(density_changes, aes(x = change_6vs1, y = reorder(species_short, change_6vs1))) +
  # Add vertical line at zero
  geom_vline(xintercept = 0, linewidth = 0.8, color = "grey30") +
  # Bars with publication colors
  geom_col(aes(fill = change_6vs1 > 0), alpha = 0.9, color = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),  # Blue for increase, orange for decrease
    guide = "none"
  ) +
  # Add small silhouettes or symbols for taxa if needed
  labs(
    x = "Standardized per-colony density (z)",
    y = NULL
  ) +
  theme_publication() +
  theme(
    axis.text.y = element_text(size = FONT_SIZE_AXIS_TEXT - 1, face = "italic"),
    axis.text.x = element_text(size = FONT_SIZE_AXIS_TEXT),
    axis.title.x = element_text(size = FONT_SIZE_AXIS_TITLE),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 0.2, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )

# Combine panels with patchwork using publication dimensions
figure3 <- fig3a / fig3b +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = FONT_SIZE_TITLE, face = "bold"))

# Save with standardized function
ggsave(file.path(PUB_DIR, "Figure3.png"), figure3,
       width = PUBLICATION_WIDTH_SINGLE,
       height = PUBLICATION_HEIGHT_TALL,
       dpi = PUBLICATION_DPI, bg = "white")
ggsave(file.path(PUB_DIR, "Figure3.pdf"), figure3,
       width = PUBLICATION_WIDTH_SINGLE,
       height = PUBLICATION_HEIGHT_TALL,
       device = cairo_pdf, bg = "white")

cat("✓ Figure 3 saved\n\n")

# FIGURE 4: PCA of coral performance -------------------------------------------
cat("Generating Figure 4: Coral performance PCA...\n")

# Process coral performance data
# Extract coral ID from Chunk column and use closed models
coral_perf <- coral_2019 %>%
  filter(`Model type` == "Mesh", grepl("closed", Model)) %>%
  mutate(reef = gsub(".*\\((.*)\\).*", "\\1", Chunk)) %>%
  select(reef, sa_2019 = `Surface Area (cm^2)`, vol_2019 = `Volume (cm^3)`) %>%
  left_join(
    coral_2021 %>%
      filter(`Model type` == "Mesh", grepl("closed", Model)) %>%
      mutate(reef = gsub(".*\\((.*)\\).*", "\\1", Chunk)) %>%
      select(reef, sa_2021 = `Surface Area (cm^2)`, vol_2021 = `Volume (cm^3)`),
    by = "reef"
  ) %>%
  mutate(
    delta_sa = sa_2021 - sa_2019,
    delta_vol = vol_2021 - vol_2019,
    growth_rate = delta_vol / sa_2019
  ) %>%
  left_join(metadata, by = "reef") %>%
  filter(!is.na(treatment), !is.na(growth_rate)) %>%
  mutate(treatment = factor(treatment, levels = c(1, 3, 6)))

# Perform PCA (placeholder - would need full performance metrics)
# For now, create a simple version
fig4a <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "PCA loadings\n(requires full coral metrics)",
           size = 5) +
  theme_pub() +
  labs(tag = "a") +
  theme(plot.tag = element_text(face = "bold", size = 14))

fig4b <- ggplot(coral_perf, aes(x = treatment, y = growth_rate, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_color_manual(values = cols_trt, guide = "none") +
  labs(x = "Number of corals", y = "Growth rate", tag = "b") +
  theme_pub() +
  theme(plot.tag = element_text(face = "bold", size = 14))

figure4 <- fig4a | fig4b
ggsave(file.path(PUB_DIR, "Figure4.png"), figure4,
       width = 10, height = 5, dpi = 600, bg = "white")
ggsave(file.path(PUB_DIR, "Figure4.pdf"), figure4,
       width = 10, height = 5, device = cairo_pdf, bg = "white")

cat("✓ Figure 4 saved\n\n")

# FIGURE 5 & 6: Placeholder (require additional data processing) ---------------
cat("Figures 5-6 require additional coral physiology data integration...\n")
cat("  (These will be generated once full performance metrics are processed)\n\n")

# SUPPLEMENTAL FIGURE S1: Individual coral metrics -----------------------------
cat("Generating Supplemental Figure S1: Individual coral metrics...\n")

# This would show violin plots of protein, carbohydrate, AFDM, zooxanthellae, growth
# Placeholder for now
figS1 <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "Supplemental Figure S1\nCoral physiology metrics\n(requires full physiology data)",
           size = 5) +
  theme_pub()

ggsave(file.path(SUP_DIR, "FigureS1.png"), figS1,
       width = 10, height = 8, dpi = 600, bg = "white")
ggsave(file.path(SUP_DIR, "FigureS1.pdf"), figS1,
       width = 10, height = 8, device = cairo_pdf, bg = "white")

cat("✓ Supplemental Figure S1 saved\n\n")

# Create a summary document ----------------------------------------------------
cat("Creating publication figure summary document...\n")

summary_md <- paste0(
  "# Publication Figures - Stier et al. CAFI 2025\n\n",
  "**Generated:** ", Sys.time(), "\n\n",
  "## Main Text Figures\n\n",
  "- **Figure 1**: Effect of coral number on CAFI abundance and richness\n",
  "  - Panel A: Total abundance with Field of Dreams expectations\n",
  "  - Panel B: Species richness\n",
  "  - Panel C: Rarefied richness\n",
  "  - Files: `Figure1.png`, `Figure1.pdf`\n\n",
  "- **Figure 2**: Species-specific responses to coral number\n",
  "  - 12 panels (4 fishes, 4 crustaceans, 4 molluscs)\n",
  "  - Observed vs expected abundance patterns\n",
  "  - Files: `Figure2.png`, `Figure2.pdf`\n\n",
  "- **Figure 3**: Community composition and density changes\n",
  "  - Panel A: NMDS ordination of CAFI communities\n",
  "  - Panel B: Top 15 species density changes (6 vs 1 coral)\n",
  "  - Files: `Figure3.png`, `Figure3.pdf`\n\n",
  "- **Figure 4**: Coral performance metrics\n",
  "  - Panel A: PCA loadings (placeholder)\n",
  "  - Panel B: Growth rate by treatment\n",
  "  - Files: `Figure4.png`, `Figure4.pdf`\n\n",
  "- **Figures 5-6**: In development (require coral physiology integration)\n\n",
  "## Supplemental Figures\n\n",
  "- **Figure S1**: Individual coral metrics (placeholder)\n",
  "  - Files: `supplemental/FigureS1.png`, `supplemental/FigureS1.pdf`\n\n",
  "## Notes\n\n",
  "- All figures generated at 600 dpi\n",
  "- Color palette: Treatment 1 (#FFD92F), Treatment 3 (#8DA0CB), Treatment 6 (#66C2A5)\n",
  "- Theme: Publication-ready with black borders and axis text\n",
  "- Bootstrap iterations: 10,000 for expected values\n",
  "- NMDS: k=2, trymax=200, Gower-like distance\n\n",
  "## To Complete\n\n",
  "1. Figure 5: CAFI PC vs coral performance PC (requires full PCA)\n",
  "2. Figure 6: Species-coral performance relationships (6 panels)\n",
  "3. Supplemental Table S1: Treatment effects (from ANOVA/ANCOVA results)\n",
  "4. Supplemental Table S2: Species-performance correlations (heatmap)\n"
)

writeLines(summary_md, file.path(PUB_DIR, "README.md"))

cat("\n=======================================================================\n")
cat("  FIGURE GENERATION COMPLETE!\n")
cat("=======================================================================\n\n")
cat("Output locations:\n")
cat("  Main figures:         ", PUB_DIR, "\n")
cat("  Supplemental figures: ", SUP_DIR, "\n\n")
cat("To generate remaining figures, ensure all analysis scripts have run.\n")
cat("See STATISTICAL_RESULTS_SUMMARY.md for verification of values.\n\n")
