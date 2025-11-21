# ==============================================================================
# NMDS with Species Loading Arrows - Draft for Joe's Suggestion
# ==============================================================================
# Purpose: Test adding species loading arrows to NMDS plot
# Author: CAFI Team
# Date: November 20, 2025
# ==============================================================================

library(here)
library(tidyverse)
library(vegan)
library(ggrepel)

# Source figure standards
source(here("scripts/MRB/mrb_figure_standards.R"))

# ==============================================================================
# 1. Load Data
# ==============================================================================

# Load treatment info
treatment_info <- read_csv(here("data/MRB Amount/coral_id_position_treatment.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    coral_id = gsub("^FE-", "", coral_id),
    # Extract reef number from position (e.g., "2-7A" -> "2")
    reef = as.numeric(gsub("-.*", "", position)),
    treatment = factor(treatment, levels = c(1, 3, 6))
  )

# Load CAFI raw data
cafi_raw <- read_csv(here("data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv"),
                     show_col_types = FALSE)

# Create species by coral data (each coral is a sample unit)
species_by_coral <- cafi_raw %>%
  mutate(
    species = paste(genus, species),
    coral_id = gsub("^FE-", "", coral_id)
  ) %>%
  left_join(treatment_info %>% select(coral_id, reef, treatment), by = "coral_id") %>%
  group_by(coral_id, species, reef, treatment) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

# Create metadata at coral level
metadata <- species_by_coral %>%
  filter(!is.na(reef)) %>%
  mutate(reef = as.character(reef)) %>%
  distinct(coral_id, reef, treatment)

# Keep species_by_reef for top 15 calculation
species_by_reef <- species_by_coral

# ==============================================================================
# 2. Create Community Matrix
# ==============================================================================

# Build community matrix at coral level
comm_matrix <- species_by_coral %>%
  filter(!is.na(reef)) %>%
  group_by(coral_id, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = count,
              values_fill = list(count = 0)) %>%
  column_to_rownames("coral_id") %>%
  as.matrix()

# Transform (square root) - standard for abundance data
comm_trans <- sqrt(comm_matrix)

# Remove species with zero total (can't contribute to ordination)
comm_trans <- comm_trans[, colSums(comm_trans) > 0]

# ==============================================================================
# 3. Run NMDS
# ==============================================================================

set.seed(1234)
# Use Bray-Curtis distance (standard for community data)
nmds_result <- metaMDS(comm_trans, distance = "bray", k = 2, trymax = 200)

cat("NMDS Stress:", round(nmds_result$stress, 3), "\n")

# ==============================================================================
# 4. Extract Scores
# ==============================================================================

# Site scores
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites")) %>%
  rownames_to_column("coral_id") %>%
  left_join(metadata, by = "coral_id") %>%
  mutate(treatment = factor(treatment, levels = c(1, 3, 6)))

# Species scores
species_scores <- as.data.frame(scores(nmds_result, display = "species")) %>%
  rownames_to_column("species")

# Calculate centroids
centroids <- nmds_scores %>%
  group_by(treatment) %>%
  summarise(
    cent_NMDS1 = mean(NMDS1),
    cent_NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

# ==============================================================================
# 5. Identify Top 15 Species (same as Panel B)
# ==============================================================================

top_15_species <- species_by_reef %>%
  group_by(species) %>%
  summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(species)

# Filter species scores to top 15
top_species_scores <- species_scores %>%
  filter(species %in% top_15_species)

cat("\nTop 15 species for arrows:\n")
print(top_species_scores)

# ==============================================================================
# 6. Create Plot with Arrows
# ==============================================================================

# Scale factor for arrows (adjust to fit plot)
arrow_scale <- 1.5

p_nmds_arrows <- ggplot() +
  # Ellipses
  stat_ellipse(data = nmds_scores,
               aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment),
               geom = "polygon", alpha = 0.25, level = 0.95, linewidth = 1.2) +

  # Species arrows
  geom_segment(data = top_species_scores,
               aes(x = 0, y = 0,
                   xend = NMDS1 * arrow_scale,
                   yend = NMDS2 * arrow_scale),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "grey30", alpha = 0.7, linewidth = 0.5) +

  # Species labels
  geom_text_repel(data = top_species_scores,
                  aes(x = NMDS1 * arrow_scale,
                      y = NMDS2 * arrow_scale,
                      label = species),
                  size = 2.5, fontface = "italic",
                  color = "grey20",
                  max.overlaps = 20,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color = "grey60",
                  segment.size = 0.3) +

  # Site points
  geom_point(data = nmds_scores,
             aes(x = NMDS1, y = NMDS2, fill = treatment),
             size = 4, alpha = 0.9, shape = 21, stroke = 1, color = "black") +

  # Centroid stars
  geom_point(data = centroids,
             aes(x = cent_NMDS1, y = cent_NMDS2, color = treatment),
             size = 8, shape = 8, stroke = 2, show.legend = FALSE) +

  # Colors
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

  # Labels
  labs(
    title = "NMDS with Species Loading Arrows (Top 15)",
    x = "NMDS1",
    y = "NMDS2",
    subtitle = paste0("Stress = ", round(nmds_result$stress, 3))
  ) +

  # Theme
  theme_ordination() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.95, size = 10)
  ) +
  coord_equal()

# ==============================================================================
# 7. Save Plot
# ==============================================================================

out_dir <- here("output/MRB/figures/nmds_draft")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "nmds_with_species_arrows.png"), p_nmds_arrows,
       width = 10, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "nmds_with_species_arrows.pdf"), p_nmds_arrows,
       width = 10, height = 10, device = cairo_pdf, bg = "white")

cat("\n✅ Saved: nmds_with_species_arrows.png/pdf\n")
cat("   Location:", out_dir, "\n")

# ==============================================================================
# 8. Alternative: Cleaner Version with Abbreviated Names
# ==============================================================================

# Create abbreviated species names (genus + first 3 letters of species)
top_species_scores_abbrev <- top_species_scores %>%
  mutate(
    species_abbrev = sapply(species, function(s) {
      parts <- strsplit(s, " ")[[1]]
      if (length(parts) >= 2) {
        paste0(substr(parts[1], 1, 1), ". ", substr(parts[2], 1, 4), ".")
      } else {
        s
      }
    })
  )

p_nmds_arrows_clean <- ggplot() +
  # Ellipses
  stat_ellipse(data = nmds_scores,
               aes(x = NMDS1, y = NMDS2, color = treatment, fill = treatment),
               geom = "polygon", alpha = 0.25, level = 0.95, linewidth = 1.2) +

  # Species arrows
  geom_segment(data = top_species_scores_abbrev,
               aes(x = 0, y = 0,
                   xend = NMDS1 * arrow_scale,
                   yend = NMDS2 * arrow_scale),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "grey40", alpha = 0.8, linewidth = 0.4) +

  # Species labels (abbreviated)
  geom_text_repel(data = top_species_scores_abbrev,
                  aes(x = NMDS1 * arrow_scale,
                      y = NMDS2 * arrow_scale,
                      label = species_abbrev),
                  size = 2.2, fontface = "italic",
                  color = "grey30",
                  max.overlaps = 20,
                  box.padding = 0.25,
                  segment.color = "grey70",
                  segment.size = 0.2) +

  # Site points
  geom_point(data = nmds_scores,
             aes(x = NMDS1, y = NMDS2, fill = treatment),
             size = 4, alpha = 0.9, shape = 21, stroke = 1, color = "black") +

  # Centroid stars
  geom_point(data = centroids,
             aes(x = cent_NMDS1, y = cent_NMDS2, color = treatment),
             size = 8, shape = 8, stroke = 2, show.legend = FALSE) +

  # Colors
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

  # Labels
  labs(
    title = "NMDS with Species Arrows (Abbreviated)",
    x = "NMDS1",
    y = "NMDS2",
    subtitle = paste0("Stress = ", round(nmds_result$stress, 3))
  ) +

  # Theme
  theme_ordination() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.95, size = 10)
  ) +
  coord_equal()

ggsave(file.path(out_dir, "nmds_with_species_arrows_abbreviated.png"),
       p_nmds_arrows_clean,
       width = 10, height = 10, dpi = 300, bg = "white")

cat("✅ Saved: nmds_with_species_arrows_abbreviated.png\n")

cat("\n========================================\n")
cat("Review both versions to see which works better.\n")
cat("The abbreviated version may be cleaner.\n")
cat("========================================\n")
