# ==============================================================================
# MRB Figure Standardization System
# ==============================================================================
# Purpose: Ensure ALL figures are publication quality with consistent formatting
# Author: CAFI Team
# Created: 2025-11-04
# ==============================================================================

# ==============================================================================
# CORE STANDARDS
# ==============================================================================

# Publication settings
PUBLICATION_DPI <- 600
PUBLICATION_WIDTH_SINGLE <- 8      # Single column figure
PUBLICATION_WIDTH_DOUBLE <- 12     # Two-column figure
PUBLICATION_WIDTH_WIDE <- 16       # Full page width
PUBLICATION_HEIGHT_STD <- 6        # Standard height
PUBLICATION_HEIGHT_TALL <- 10      # Tall figures
PUBLICATION_HEIGHT_SQUARE <- 8     # Square figures

# Font sizes (all in points)
FONT_SIZE_BASE <- 14               # Base font size
FONT_SIZE_TITLE <- 18              # Plot titles
FONT_SIZE_SUBTITLE <- 14          # Subtitles
FONT_SIZE_AXIS_TITLE <- 14        # Axis titles
FONT_SIZE_AXIS_TEXT <- 12         # Axis text
FONT_SIZE_LEGEND_TITLE <- 14      # Legend title
FONT_SIZE_LEGEND_TEXT <- 12       # Legend text
FONT_SIZE_FACET <- 12              # Facet labels
FONT_SIZE_ANNOTATION <- 10        # Annotations

# Spacing (in lines or cm)
MARGIN_TOP <- 0.5
MARGIN_RIGHT <- 0.5
MARGIN_BOTTOM <- 0.5
MARGIN_LEFT <- 0.5
LEGEND_KEY_SIZE <- 1.5             # Size of legend keys
FACET_SPACING <- 0.5               # Space between facets

# ==============================================================================
# COLOR PALETTES
# ==============================================================================

# Treatment colors (colorblind-friendly)
TREATMENT_COLORS <- c(
  "1" = "#E69F00",   # Orange
  "3" = "#56B4E9",   # Sky blue
  "6" = "#009E73"    # Green
)

# Extended color palette for multiple categories
CATEGORY_COLORS <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky blue
  "#009E73",  # Green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermillion
  "#CC79A7",  # Reddish purple
  "#999999"   # Grey
)

# Gradient colors
GRADIENT_LOW <- "#0571B0"
GRADIENT_MID <- "#F7F7F7"
GRADIENT_HIGH <- "#CA0020"

# Single-color accent for non-treatment bar charts
# Uses CATEGORY_COLORS[5] (Blue) to distinguish from treatment colors
ACCENT_COLOR <- "#0072B2"         # Primary accent (Blue)
ACCENT_COLOR_ALT <- "#D55E00"     # Secondary accent (Vermillion)

# ==============================================================================
# MASTER THEME FUNCTION
# ==============================================================================

theme_publication <- function(base_size = FONT_SIZE_BASE,
                            base_family = "",
                            grid = FALSE) {

  # Start with minimal theme
  theme_bw(base_size = base_size, base_family = base_family) +
  theme(
    # Panel
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = if(grid) element_line(color = "grey90", linewidth = 0.25) else element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),

    # Axes
    axis.title.x = element_text(size = FONT_SIZE_AXIS_TITLE,
                               face = "bold",
                               margin = margin(t = 10)),
    axis.title.y = element_text(size = FONT_SIZE_AXIS_TITLE,
                               face = "bold",
                               margin = margin(r = 10)),
    axis.text.x = element_text(size = FONT_SIZE_AXIS_TEXT,
                               color = "black",
                               margin = margin(t = 5)),
    axis.text.y = element_text(size = FONT_SIZE_AXIS_TEXT,
                               color = "black",
                               margin = margin(r = 5)),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),

    # Title and labels
    plot.title = element_text(size = FONT_SIZE_TITLE,
                             face = "bold",
                             hjust = 0.5,
                             margin = margin(b = 10)),
    plot.subtitle = element_text(size = FONT_SIZE_SUBTITLE,
                                 hjust = 0.5,
                                 margin = margin(b = 10)),
    plot.caption = element_text(size = FONT_SIZE_ANNOTATION,
                                hjust = 1,
                                color = "grey40"),

    # Legend
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = FONT_SIZE_LEGEND_TITLE,
                               face = "bold"),
    legend.text = element_text(size = FONT_SIZE_LEGEND_TEXT),
    legend.key.size = unit(LEGEND_KEY_SIZE, "lines"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(t = 10),

    # Facets
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 1),
    strip.text = element_text(size = FONT_SIZE_FACET,
                             face = "bold",
                             margin = margin(3, 3, 3, 3)),
    panel.spacing = unit(FACET_SPACING, "cm"),

    # Overall margins
    plot.margin = margin(t = MARGIN_TOP,
                        r = MARGIN_RIGHT,
                        b = MARGIN_BOTTOM,
                        l = MARGIN_LEFT,
                        unit = "cm"),

    # Background
    plot.background = element_rect(fill = "white", color = NA)
  )
}

# ==============================================================================
# SPECIALIZED THEMES
# ==============================================================================

# Theme for multi-panel figures
theme_multipanel <- function() {
  theme_publication(base_size = 12) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "horizontal"
  )
}

# Theme for ordination plots
theme_ordination <- function() {
  theme_publication(grid = TRUE) +
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
    legend.position = "right"
  )
}

# Theme for heatmaps
theme_heatmap <- function() {
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )
}

# ==============================================================================
# COLOR SCALE FUNCTIONS
# ==============================================================================

# Apply treatment colors consistently
scale_color_treatment <- function(...) {
  scale_color_manual(
    values = TREATMENT_COLORS,
    labels = c("1" = "1 Colony", "3" = "3 Colonies", "6" = "6 Colonies"),
    name = "Treatment",
    ...
  )
}

scale_fill_treatment <- function(...) {
  scale_fill_manual(
    values = TREATMENT_COLORS,
    labels = c("1" = "1 Colony", "3" = "3 Colonies", "6" = "6 Colonies"),
    name = "Treatment",
    ...
  )
}

# Categorical colors for non-treatment factors
scale_color_category <- function(n = NULL, ...) {
  if (is.null(n)) {
    scale_color_manual(values = CATEGORY_COLORS, ...)
  } else {
    scale_color_manual(values = CATEGORY_COLORS[1:n], ...)
  }
}

scale_fill_category <- function(n = NULL, ...) {
  if (is.null(n)) {
    scale_fill_manual(values = CATEGORY_COLORS, ...)
  } else {
    scale_fill_manual(values = CATEGORY_COLORS[1:n], ...)
  }
}

# Gradient scales
scale_fill_gradient_publication <- function(...) {
  scale_fill_gradient2(
    low = GRADIENT_LOW,
    mid = GRADIENT_MID,
    high = GRADIENT_HIGH,
    midpoint = 0,
    ...
  )
}

scale_color_gradient_publication <- function(...) {
  scale_color_gradient2(
    low = GRADIENT_LOW,
    mid = GRADIENT_MID,
    high = GRADIENT_HIGH,
    midpoint = 0,
    ...
  )
}

# ==============================================================================
# FIGURE SAVING FUNCTION
# ==============================================================================

save_figure <- function(plot,
                       filename,
                       width = PUBLICATION_WIDTH_SINGLE,
                       height = PUBLICATION_HEIGHT_STD,
                       dpi = PUBLICATION_DPI,
                       units = "in") {

  # Ensure directory exists
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  # Save as PNG
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = units,
    bg = "white"
  )

  # Save as PDF
  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    units = units,
    bg = "white"
  )

  message("✅ Saved: ", basename(filename), " (.png and .pdf)")
  message("   Dimensions: ", width, " x ", height, " inches at ", dpi, " DPI")
}

# ==============================================================================
# PLOT ENHANCEMENT FUNCTIONS
# ==============================================================================

# Add significance stars to plot
add_significance <- function(plot, p_value, x, y, size = 5) {
  sig_label <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ ".",
    TRUE ~ "ns"
  )

  plot +
    annotate("text", x = x, y = y, label = sig_label,
             size = size, fontface = "bold")
}

# Format axis labels
format_axis_scientific <- function(x) {
  parse(text = gsub("e\\+", " %*% 10^", scales::scientific(x)))
}

format_axis_percent <- function(x) {
  scales::percent(x, accuracy = 1)
}

format_axis_comma <- function(x) {
  scales::comma(x, accuracy = 1)
}

# ==============================================================================
# STANDARD PLOT TYPES
# ==============================================================================

# Standard boxplot with points
create_boxplot <- function(data, x, y, fill = NULL,
                          title = "", subtitle = "",
                          xlab = "", ylab = "") {

  p <- ggplot(data, aes_string(x = x, y = y, fill = fill)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2),
               alpha = 0.5, size = 2) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_publication()

  if (!is.null(fill) && fill == "treatment") {
    p <- p + scale_fill_treatment()
  }

  return(p)
}

# Standard bar plot with error bars
create_barplot <- function(data, x, y, error = NULL, fill = NULL,
                          title = "", subtitle = "",
                          xlab = "", ylab = "") {

  p <- ggplot(data, aes_string(x = x, y = y, fill = fill))

  if (!is.null(error)) {
    p <- p +
      geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.5) +
      geom_errorbar(aes_string(ymin = paste0(y, " - ", error),
                               ymax = paste0(y, " + ", error)),
                    width = 0.2, linewidth = 0.5)
  } else {
    p <- p + geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.5)
  }

  p <- p +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab) +
    theme_publication()

  if (!is.null(fill) && fill == "treatment") {
    p <- p + scale_fill_treatment()
  }

  return(p)
}

# ==============================================================================
# EXPORT MESSAGE
# ==============================================================================

cat("✅ MRB Figure Standards loaded\n")
cat("   • Base font size: ", FONT_SIZE_BASE, "\n")
cat("   • Publication DPI: ", PUBLICATION_DPI, "\n")
cat("   • Treatment colors defined\n")
cat("   • theme_publication() available\n")
cat("   • Specialized themes loaded\n")
cat("   • Helper functions ready\n\n")