# ==============================================================================
# Snail Corallivore Analysis
# ==============================================================================
# Purpose: Analyze corallivorous snail patterns in response to coral number
#          and their relationships with coral performance
#
# Addresses reviewer comment: "Seems odd not to mention the snails, especially
# given their propensity to eat coral and that of the three broad groups they
# showed the strongest and most consistent hyperaggregation on 6-coral reefs"
#
# Outputs:
#   - output/MRB/snail_analysis/snail_corallivore_analysis.md
#   - output/MRB/snail_analysis/snail_corallivore_analysis.html
#   - Associated figures
# ==============================================================================

# Load libraries
library(here)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(knitr)
library(rmarkdown)

# Load figure standards
source(here("scripts/MRB/mrb_figure_standards.R"))

# Set output directory
out_dir <- here("output/MRB/snail_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load Data
# ==============================================================================

# Load CAFI community data
cafi_raw <- read_csv(here("data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv"),
                     show_col_types = FALSE)

# Load physiology data (has treatment info, percent_alive, and performance)
physio_file <- here("output/MRB/figures/coral/physio/physio_metrics_plus_growth_filtered.csv")
physio_df <- read_csv(physio_file, show_col_types = FALSE)

# Use physio_df as growth_df (already filtered to >=80% alive)
growth_df <- physio_df

# ==============================================================================
# 2. Identify Corallivorous Snails
# ==============================================================================

# Known corallivorous/coral-associated gastropods
# Based on literature and taxonomic knowledge
corallivore_genera <- c(

  "Drupella",        # Well-known corallivore

  "Coralliophila",   # Corallivore
  "Quoyula",         # Corallivore
  "Galeropsis",      # Found in coral analysis
  "Morula",          # Predatory snail
  "Strigatella",     # Found in coral habitats
  "Macteola",        # Found in coral analysis
  "Vexillum",        # Found in coral habitats
  "Mitrella"         # Found in coral analysis
)

# Create species_name column and standardize coral_id
cafi_raw <- cafi_raw %>%
  mutate(
    species_name = paste(genus, species),
    # Strip "FE-" prefix from coral_id to match physio data
    coral_id = gsub("^FE-", "", coral_id)
  )

# Get all unique species from CAFI data
all_species <- unique(cafi_raw$species_name)

# Identify snail species (gastropods)
snail_species <- all_species[grepl(paste(corallivore_genera, collapse = "|"),
                                    all_species, ignore.case = TRUE)]

cat("Identified snail species:\n")
print(snail_species)

# ==============================================================================
# 3. Build Snail Abundance Matrix
# ==============================================================================

# Filter to snails and summarize by coral
snail_data <- cafi_raw %>%
  filter(grepl(paste(corallivore_genera, collapse = "|"),
               species_name, ignore.case = TRUE)) %>%
  group_by(coral_id, species_name) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species_name,
              values_from = count, values_fill = 0)

# Calculate total snail abundance per coral
snail_totals <- snail_data %>%
  mutate(total_snails = rowSums(select(., -coral_id))) %>%
  select(coral_id, total_snails)

# Add species richness
snail_richness <- cafi_raw %>%
  filter(grepl(paste(corallivore_genera, collapse = "|"),
               species_name, ignore.case = TRUE)) %>%
  group_by(coral_id) %>%
  summarise(
    snail_richness = n_distinct(species_name[count > 0]),
    .groups = "drop"
  )

# Merge with treatment data (growth_df already filtered to >=80% alive)
snail_summary <- growth_df %>%
  select(coral_id, reef, treatment, percent_alive) %>%
  left_join(snail_totals, by = "coral_id") %>%
  left_join(snail_richness, by = "coral_id") %>%
  mutate(
    total_snails = replace_na(total_snails, 0),
    snail_richness = replace_na(snail_richness, 0),
    treatment = factor(treatment, levels = c("1", "3", "6")),
    reef = as.factor(reef)
  )

# ==============================================================================
# 4. Compute Coral Performance (PC1)
# ==============================================================================

# Get physiology variables for PCA
if ("protein_mg_cm2" %in% names(physio_df)) {
  physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

  # Check which vars exist
  physio_vars <- physio_vars[physio_vars %in% names(physio_df)]

  if (length(physio_vars) >= 2) {
    # Compute PC1 of coral performance
    physio_complete <- physio_df %>%
      filter(complete.cases(select(., all_of(physio_vars))))

    pca_result <- prcomp(physio_complete[, physio_vars], center = TRUE, scale. = TRUE)

    physio_complete$PC1_coral <- pca_result$x[, 1]

    # Merge with snail data
    snail_summary <- snail_summary %>%
      left_join(physio_complete %>% select(coral_id, PC1_coral, all_of(physio_vars)),
                by = "coral_id")

    pc1_var_explained <- summary(pca_result)$importance[2, 1] * 100
  }
} else {
  # Use growth as performance proxy
  snail_summary <- snail_summary %>%
    left_join(growth_df %>% select(coral_id, growth_vol_b), by = "coral_id") %>%
    mutate(PC1_coral = scale(growth_vol_b)[,1])

  pc1_var_explained <- 100
}

# ==============================================================================
# 5. Statistical Analyses
# ==============================================================================

results <- list()

# 5.1 Snail abundance by treatment
tryCatch({
  model_abundance <- lmer(total_snails ~ treatment + (1|reef),
                          data = snail_summary)
  results$abundance_anova <- anova(model_abundance)
  results$abundance_emmeans <- emmeans(model_abundance, pairwise ~ treatment)
}, error = function(e) {
  cat("Warning: Abundance model failed:", e$message, "\n")
  # Use simpler model without random effects
  model_abundance <- lm(total_snails ~ treatment, data = snail_summary)
  results$abundance_anova <<- anova(model_abundance)
})

# 5.2 Snail richness by treatment
tryCatch({
  model_richness <- lmer(snail_richness ~ treatment + (1|reef),
                         data = snail_summary)
  results$richness_anova <- anova(model_richness)
  results$richness_emmeans <- emmeans(model_richness, pairwise ~ treatment)
}, error = function(e) {
  cat("Warning: Richness model failed:", e$message, "\n")
  model_richness <- lm(snail_richness ~ treatment, data = snail_summary)
  results$richness_anova <<- anova(model_richness)
})

# 5.3 Snail-performance correlation
if ("PC1_coral" %in% names(snail_summary)) {
  snail_perf <- snail_summary %>% filter(!is.na(PC1_coral))

  # Overall correlation
  results$snail_pc1_cor <- cor.test(snail_perf$total_snails,
                                    snail_perf$PC1_coral)

  # Mixed model
  tryCatch({
    model_perf <- lmer(PC1_coral ~ total_snails + (1|reef),
                       data = snail_perf)
    results$perf_model <- summary(model_perf)
    results$perf_anova <- anova(model_perf)
  }, error = function(e) {
    cat("Warning: Performance model failed:", e$message, "\n")
  })
}

# 5.4 Per-species patterns
species_cols <- setdiff(names(snail_data), "coral_id")
species_results <- list()

for (sp in species_cols) {
  sp_data <- snail_summary %>%
    left_join(snail_data %>% select(coral_id, all_of(sp)), by = "coral_id") %>%
    rename(sp_count = all_of(sp)) %>%
    mutate(sp_count = replace_na(sp_count, 0))

  # Treatment effect
  if (sum(sp_data$sp_count > 0) >= 5) {
    model_sp <- tryCatch(
      lmer(sp_count ~ treatment + (1|reef), data = sp_data),
      error = function(e) NULL
    )

    if (!is.null(model_sp)) {
      sp_anova <- anova(model_sp)

      # Means by treatment
      sp_means <- sp_data %>%
        group_by(treatment) %>%
        summarise(
          mean = mean(sp_count),
          se = sd(sp_count) / sqrt(n()),
          .groups = "drop"
        )

      # Correlation with PC1
      if ("PC1_coral" %in% names(sp_data)) {
        sp_cor <- cor.test(sp_data$sp_count, sp_data$PC1_coral, use = "complete.obs")
      } else {
        sp_cor <- NULL
      }

      species_results[[sp]] <- list(
        anova = sp_anova,
        means = sp_means,
        correlation = sp_cor
      )
    }
  }
}

results$species <- species_results

# ==============================================================================
# 6. Create Figures
# ==============================================================================

# Treatment colors
cols_trt <- c("1" = "#E69F00", "3" = "#56B4E9", "6" = "#009E73")

# 6.1 Snail abundance by treatment
p_abundance <- ggplot(snail_summary, aes(x = treatment, y = total_snails, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  scale_fill_manual(values = cols_trt, guide = "none") +
  labs(
    title = "Snail Abundance by Coral Number",
    x = "Number of Corals",
    y = "Total Snail Count per Coral"
  ) +
  theme_publication()

ggsave(file.path(out_dir, "snail_abundance_by_treatment.png"), p_abundance,
       width = 6, height = 5, dpi = 600)
ggsave(file.path(out_dir, "snail_abundance_by_treatment.pdf"), p_abundance,
       width = 6, height = 5)

# 6.2 Snail richness by treatment
p_richness <- ggplot(snail_summary, aes(x = treatment, y = snail_richness, fill = treatment)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  scale_fill_manual(values = cols_trt, guide = "none") +
  labs(
    title = "Snail Species Richness by Coral Number",
    x = "Number of Corals",
    y = "Number of Snail Species per Coral"
  ) +
  theme_publication()

ggsave(file.path(out_dir, "snail_richness_by_treatment.png"), p_richness,
       width = 6, height = 5, dpi = 600)
ggsave(file.path(out_dir, "snail_richness_by_treatment.pdf"), p_richness,
       width = 6, height = 5)

# 6.3 Snail abundance vs coral performance
if ("PC1_coral" %in% names(snail_summary)) {
  p_performance <- ggplot(snail_summary %>% filter(!is.na(PC1_coral)),
                          aes(x = total_snails, y = PC1_coral)) +
    geom_point(aes(fill = treatment), shape = 21, size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_fill_manual(values = cols_trt, name = "Coral number") +
    labs(
      title = "Snail Abundance vs Coral Performance",
      x = "Total Snail Count",
      y = expression(PC1[coral]~(performance))
    ) +
    theme_publication()

  ggsave(file.path(out_dir, "snail_vs_performance.png"), p_performance,
         width = 7, height = 5, dpi = 600)
  ggsave(file.path(out_dir, "snail_vs_performance.pdf"), p_performance,
         width = 7, height = 5)
}

# 6.4 Species-level heatmap of treatment means
if (length(species_results) > 0) {
  species_means_df <- map_dfr(names(species_results), function(sp) {
    species_results[[sp]]$means %>%
      mutate(species = sp)
  })

  # Calculate fold change from 1-coral to 6-coral
  fold_changes <- species_means_df %>%
    select(species, treatment, mean) %>%
    pivot_wider(names_from = treatment, values_from = mean) %>%
    mutate(
      fold_change = ifelse(`1` > 0, `6` / `1`, `6`),
      direction = ifelse(`6` > `1`, "Higher on 6-coral", "Higher on 1-coral")
    ) %>%
    arrange(desc(fold_change))

  p_species <- species_means_df %>%
    mutate(species = factor(species, levels = fold_changes$species)) %>%
    ggplot(aes(x = treatment, y = species, fill = mean)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#009E73", name = "Mean\ncount") +
    labs(
      title = "Snail Species Abundance by Treatment",
      x = "Number of Corals",
      y = "Species"
    ) +
    theme_publication() +
    theme(axis.text.y = element_text(face = "italic", size = 8))

  ggsave(file.path(out_dir, "snail_species_heatmap.png"), p_species,
         width = 7, height = max(4, length(species_results) * 0.4), dpi = 600)
  ggsave(file.path(out_dir, "snail_species_heatmap.pdf"), p_species,
         width = 7, height = max(4, length(species_results) * 0.4))
}

# ==============================================================================
# 7. Generate Markdown Report
# ==============================================================================

# Summary statistics
summary_stats <- snail_summary %>%
  group_by(treatment) %>%
  summarise(
    n = n(),
    mean_abundance = mean(total_snails),
    se_abundance = sd(total_snails) / sqrt(n()),
    mean_richness = mean(snail_richness),
    se_richness = sd(snail_richness) / sqrt(n()),
    .groups = "drop"
  )

# Extract key statistics
abund_F <- results$abundance_anova$`F value`[1]
abund_p <- results$abundance_anova$`Pr(>F)`[1]
rich_F <- results$richness_anova$`F value`[1]
rich_p <- results$richness_anova$`Pr(>F)`[1]

# Performance correlation
if (!is.null(results$snail_pc1_cor)) {
  perf_r <- results$snail_pc1_cor$estimate
  perf_p <- results$snail_pc1_cor$p.value
} else {
  perf_r <- NA
  perf_p <- NA
}

# Create markdown content
md_content <- sprintf('---
title: "Corallivorous Snail Analysis"
author: "MRB CAFI Analysis"
date: "%s"
output:
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
    highlight: tango
---

# Overview

This analysis examines corallivorous and coral-associated snail (gastropod) patterns in response to experimental coral density manipulation. We address the reviewer comment regarding the notable hyperaggregation of snails on multi-coral reefs.
## Key Questions

1. **How did snail abundance shift with coral number?**
2. **Were snail distributions clustered (hyperaggregated) on 6-coral reefs?**
3. **How were snails correlated with coral performance?**

---

# Data Summary

## Snail Species Identified

The following gastropod genera were identified as coral-associated snails:

%s

**Total snail species in dataset:** %d

---

# Results

## 1. Snail Abundance by Coral Number

![Snail Abundance](snail_abundance_by_treatment.png)

### Summary Statistics

| Treatment | n | Mean Abundance | SE |
|-----------|---|----------------|-----|
| 1 coral | %d | %.2f | %.2f |
| 3 corals | %d | %.2f | %.2f |
| 6 corals | %d | %.2f | %.2f |

### Statistical Test

**Linear Mixed Model:** `total_snails ~ treatment + (1|reef)`

- **F-statistic:** %.2f
- **p-value:** %.4f

%s

---

## 2. Snail Species Richness

![Snail Richness](snail_richness_by_treatment.png)

### Summary Statistics

| Treatment | Mean Richness | SE |
|-----------|---------------|-----|
| 1 coral | %.2f | %.2f |
| 3 corals | %.2f | %.2f |
| 6 corals | %.2f | %.2f |

### Statistical Test

**Linear Mixed Model:** `snail_richness ~ treatment + (1|reef)`

- **F-statistic:** %.2f
- **p-value:** %.4f

%s

---

## 3. Hyperaggregation Pattern

To assess whether snails showed hyperaggregation (stronger-than-expected increase with coral number), we compare the fold-increase from 1-coral to 6-coral treatments.

**Expected under proportional scaling:** 6× increase (6 corals vs 1 coral)

**Observed:**
- Mean snails on 1-coral reefs: %.2f
- Mean snails on 6-coral reefs: %.2f
- **Fold increase:** %.2f×

%s

---

## 4. Snail Abundance vs Coral Performance

![Snail vs Performance](snail_vs_performance.png)

### Correlation Analysis

%s

---

## 5. Species-Level Patterns

![Species Heatmap](snail_species_heatmap.png)

### Individual Species Results

%s

---
# Summary and Interpretation

## Key Findings

1. **Snail abundance %s with coral number** (F = %.2f, p = %.4f)

2. **Snail richness %s with coral number** (F = %.2f, p = %.4f)

3. **Hyperaggregation:** Snails showed a %.1f-fold increase from 1-coral to 6-coral treatments, %s the 6× expected under proportional scaling.

4. **Performance relationship:** %s

## Biological Interpretation

%s

---

# Conclusions

This analysis provides quantitative support for the observation that snails showed %s on multi-coral reefs. The %s relationship with coral performance suggests that %s.

These patterns are consistent with the hypothesis that %s.

---

*Analysis generated: %s*
',
  format(Sys.Date(), "%%B %%d, %%Y"),
  paste("- *", snail_species, "*", collapse = "\n"),
  length(snail_species),
  summary_stats$n[1], summary_stats$mean_abundance[1], summary_stats$se_abundance[1],
  summary_stats$n[2], summary_stats$mean_abundance[2], summary_stats$se_abundance[2],
  summary_stats$n[3], summary_stats$mean_abundance[3], summary_stats$se_abundance[3],
  abund_F, abund_p,
  ifelse(abund_p < 0.05, "**Significant treatment effect detected.**", "Treatment effect not statistically significant."),
  summary_stats$mean_richness[1], summary_stats$se_richness[1],
  summary_stats$mean_richness[2], summary_stats$se_richness[2],
  summary_stats$mean_richness[3], summary_stats$se_richness[3],
  rich_F, rich_p,
  ifelse(rich_p < 0.05, "**Significant treatment effect detected.**", "Treatment effect not statistically significant."),
  summary_stats$mean_abundance[1],
  summary_stats$mean_abundance[3],
  ifelse(summary_stats$mean_abundance[1] > 0,
         summary_stats$mean_abundance[3] / summary_stats$mean_abundance[1],
         summary_stats$mean_abundance[3]),
  ifelse(summary_stats$mean_abundance[3] / max(summary_stats$mean_abundance[1], 0.1) > 6,
         "This represents **hyperaggregation** - snails increased more than proportionally with coral number.",
         ifelse(summary_stats$mean_abundance[3] / max(summary_stats$mean_abundance[1], 0.1) < 6,
                "This represents **hypoaggregation** - snails increased less than proportionally.",
                "This represents approximately proportional scaling.")),
  ifelse(!is.na(perf_r),
         sprintf("- **Pearson correlation (r):** %.3f\n- **p-value:** %.4f\n\n%s",
                 perf_r, perf_p,
                 ifelse(perf_p < 0.05,
                        ifelse(perf_r < 0,
                               "**Significant negative correlation:** Higher snail abundance associated with lower coral performance.",
                               "**Significant positive correlation:** Higher snail abundance associated with better coral performance."),
                        "No significant correlation between snail abundance and coral performance.")),
         "Correlation analysis not available due to missing performance data."),
  # Species results section
  ifelse(length(species_results) > 0,
         paste(sapply(names(species_results), function(sp) {
           res <- species_results[[sp]]
           sprintf("**%s**\n- Treatment effect: F = %.2f, p = %.4f\n- Correlation with PC1: r = %.3f, p = %.4f\n",
                   sp,
                   res$anova$`F value`[1],
                   res$anova$`Pr(>F)`[1],
                   ifelse(!is.null(res$correlation), res$correlation$estimate, NA),
                   ifelse(!is.null(res$correlation), res$correlation$p.value, NA))
         }), collapse = "\n"),
         "No species had sufficient data for individual analysis."),
  # Summary section
  ifelse(abund_p < 0.05, "increased significantly", "did not change significantly"),
  abund_F, abund_p,
  ifelse(rich_p < 0.05, "increased significantly", "did not change significantly"),
  rich_F, rich_p,
  ifelse(summary_stats$mean_abundance[1] > 0,
         summary_stats$mean_abundance[3] / summary_stats$mean_abundance[1],
         summary_stats$mean_abundance[3]),
  ifelse(summary_stats$mean_abundance[3] / max(summary_stats$mean_abundance[1], 0.1) > 6,
         "exceeding", "below"),
  ifelse(!is.na(perf_r) && perf_p < 0.05,
         sprintf("Snail abundance was %s correlated with coral performance (r = %.3f, p = %.4f)",
                 ifelse(perf_r < 0, "negatively", "positively"), perf_r, perf_p),
         "Snail abundance was not significantly correlated with coral performance"),
  # Biological interpretation
  ifelse(!is.na(perf_r) && perf_r < 0 && perf_p < 0.05,
         "The negative relationship between snail abundance and coral performance could indicate: (1) corallivorous snails directly harm corals through tissue consumption, (2) snails preferentially colonize already-stressed corals, or (3) both processes operate simultaneously. Distinguishing between these hypotheses requires experimental manipulation of snail presence.",
         ifelse(!is.na(perf_r) && perf_r > 0,
                "The positive relationship between snail abundance and coral performance suggests these gastropods may not be primarily corallivorous, or that healthy corals support larger snail populations through increased habitat complexity or food availability.",
                "The lack of relationship between snail abundance and coral performance suggests these gastropods may have neutral effects on coral hosts, or that positive and negative effects cancel out across species.")),
  # Conclusions
  ifelse(summary_stats$mean_abundance[3] / max(summary_stats$mean_abundance[1], 0.1) > 6,
         "strong hyperaggregation",
         "aggregation patterns"),
  ifelse(!is.na(perf_r) && perf_p < 0.05,
         ifelse(perf_r < 0, "negative", "positive"),
         "non-significant"),
  ifelse(!is.na(perf_r) && perf_r < 0 && perf_p < 0.05,
         "increased snail density on multi-coral reefs may contribute to the observed decline in coral performance, though causality cannot be established from correlational data",
         "snail abundance patterns are largely independent of coral performance metrics"),
  ifelse(summary_stats$mean_abundance[3] / max(summary_stats$mean_abundance[1], 0.1) > 6 &&
           !is.na(perf_r) && perf_r < 0,
         "the hyperaggregation of potentially harmful gastropods on multi-coral reefs represents one mechanism through which increased coral density leads to decreased per-capita coral performance",
         "snail community assembly on experimental reefs follows density-dependent patterns that warrant further investigation"),
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

# Write markdown file
md_file <- file.path(out_dir, "snail_corallivore_analysis.md")
writeLines(md_content, md_file)

cat("\nMarkdown report written to:", md_file, "\n")

# Render to HTML
html_file <- file.path(out_dir, "snail_corallivore_analysis.html")
tryCatch({
  rmarkdown::render(md_file, output_file = html_file, quiet = TRUE)
  cat("HTML report written to:", html_file, "\n")
}, error = function(e) {
  cat("Note: Could not render HTML. Markdown file available at:", md_file, "\n")
})

# Update todos
cat("\n✅ Snail corallivore analysis complete!\n")
cat("   Outputs saved to:", out_dir, "\n")

# Save results object for reference
saveRDS(results, file.path(out_dir, "snail_analysis_results.rds"))
