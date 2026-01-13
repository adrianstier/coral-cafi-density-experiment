# ==============================================================================
# MRB Analysis Configuration File
# ==============================================================================
# Purpose: Centralized configuration for all MRB analysis scripts
#          Ensures consistency in colors, themes, and formatting
# Author: CAFI Team
# Created: 2025-11-04
# ==============================================================================

# Load figure standards
source("scripts/MRB/mrb_figure_standards.R")

# ==============================================================================
# TREATMENT LABELS AND ORDER
# ==============================================================================
# NOTE: TREATMENT_COLORS are defined in mrb_figure_standards.R (sourced above)
#       Do NOT redefine them here to avoid divergence

# Alternative names for treatments
TREATMENT_LABELS <- c(
  "1" = "1 Colony",
  "3" = "3 Colonies",
  "6" = "6 Colonies"
)

# Order for treatment display
TREATMENT_ORDER <- c("1", "3", "6")

# ==============================================================================
# PUBLICATION THEME SETTINGS
# ==============================================================================

# Base theme function for all plots
theme_mrb <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # Panel
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),

      # Axes
      axis.text = element_text(color = "black", size = base_size * 0.9),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),

      # Title and labels
      plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      plot.caption = element_text(size = base_size * 0.8, hjust = 1),

      # Legend
      legend.position = "bottom",
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size * 0.9),
      legend.key.size = unit(1.2, "lines"),
      legend.background = element_rect(fill = "white", color = NA),

      # Facets
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
      strip.text = element_text(size = base_size * 0.9, face = "bold"),

      # Spacing
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# FIGURE SPECIFICATIONS
# ==============================================================================

# Standard figure sizes (inches)
FIG_SIZES <- list(
  single = c(width = 8, height = 6),      # Single panel figures
  double = c(width = 12, height = 6),     # Two-panel figures
  triple = c(width = 16, height = 6),     # Three-panel figures
  square = c(width = 8, height = 8),      # Square figures (e.g., ordination)
  tall = c(width = 8, height = 10),       # Tall figures
  wide = c(width = 12, height = 8)        # Wide figures
)

# Publication quality settings
PUB_DPI <- 600  # DPI for PNG output
PUB_DEVICE <- cairo_pdf  # Device for PDF output

# ==============================================================================
# STATISTICAL SETTINGS
# ==============================================================================

# Significance levels
ALPHA <- 0.05
CI_LEVEL <- 0.95

# Bootstrap settings
BOOT_ITERATIONS <- 10000  # For publication
BOOT_ITERATIONS_QUICK <- 1000  # For testing

# Permutation settings
PERM_ITERATIONS <- 999  # Standard for permutation tests
PERM_ITERATIONS_QUICK <- 99  # For testing

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Save figure in both PNG and PDF formats with consistent settings
#'
#' @param plot ggplot object
#' @param filename Base filename (without extension)
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi DPI for PNG output
save_publication_figure <- function(plot,
                                   filename,
                                   width = FIG_SIZES$single["width"],
                                   height = FIG_SIZES$single["height"],
                                   dpi = PUB_DPI) {

  # Ensure directory exists
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  # Save PNG
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  # Save PDF
  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = PUB_DEVICE,
    bg = "white"
  )

  message("✅ Saved: ", basename(filename), " (.png and .pdf)")
}

#' Apply consistent treatment colors to a ggplot
#'
#' @param plot ggplot object
#' @return ggplot object with treatment colors applied
apply_treatment_colors <- function(plot) {
  plot +
    scale_color_manual(
      values = TREATMENT_COLORS,
      labels = TREATMENT_LABELS,
      name = "Treatment"
    ) +
    scale_fill_manual(
      values = TREATMENT_COLORS,
      labels = TREATMENT_LABELS,
      name = "Treatment"
    )
}

#' Format p-values consistently
#'
#' @param p p-value
#' @return Formatted string
format_p_value <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  } else if (p < 0.01) {
    return(sprintf("p = %.3f", p))
  } else if (p < 0.05) {
    return(sprintf("p = %.3f", p))
  } else {
    return(sprintf("p = %.2f", p))
  }
}

#' Create consistent axis labels
#'
#' @param label Text for label
#' @param unit Unit in parentheses (optional)
#' @return Formatted label
make_axis_label <- function(label, unit = NULL) {
  if (!is.null(unit)) {
    return(bquote(bold(.(label)) ~ "(" * .(unit) * ")"))
  } else {
    return(bquote(bold(.(label))))
  }
}

# ==============================================================================
# DATA PROCESSING HELPERS
# ==============================================================================

#' Ensure treatment is ordered factor
#'
#' @param data Data frame with treatment column
#' @return Data frame with treatment as ordered factor
order_treatments <- function(data) {
  if ("treatment" %in% names(data)) {
    data %>%
      mutate(treatment = factor(treatment,
                                levels = TREATMENT_ORDER,
                                labels = TREATMENT_LABELS,
                                ordered = TRUE))
  } else {
    data
  }
}

#' Add significance stars
#'
#' @param p p-value
#' @return Significance stars
get_significance_stars <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ "ns"
  )
}

# ==============================================================================
# TABLE FORMATTING
# ==============================================================================

#' Create consistent gt table theme
#'
#' @param gt_object gt table object
#' @return Styled gt table
style_gt_table <- function(gt_object) {
  gt_object %>%
    tab_options(
      table.font.size = px(12),
      heading.title.font.size = px(16),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = px(14),
      column_labels.font.weight = "bold",
      column_labels.font.size = px(12),
      table.border.top.color = "black",
      table.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      heading.border.bottom.color = "black",
      column_labels.border.bottom.color = "black",
      table.width = pct(100)
    )
}

# ==============================================================================
# LOGGING
# ==============================================================================

#' Print section header
#'
#' @param text Section title
print_section <- function(text) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(text, "\n")
  cat(strrep("=", 70), "\n\n", sep = "")
}

#' Print subsection header
#'
#' @param text Subsection title
print_subsection <- function(text) {
  cat("\n", strrep("-", 50), "\n", sep = "")
  cat(text, "\n")
  cat(strrep("-", 50), "\n\n", sep = "")
}

# ==============================================================================
# EXPORT MESSAGE
# ==============================================================================

cat("✅ MRB configuration loaded\n")
cat("   • Treatment colors defined\n")
cat("   • Publication theme available: theme_mrb()\n")
cat("   • Helper functions loaded\n")
cat("   • Figure specifications set (DPI:", PUB_DPI, ")\n\n")