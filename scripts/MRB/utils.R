# ==============================================================================
# Utility Functions for MRB CAFI Analysis
# ==============================================================================
# Purpose: Centralized utility functions to reduce code duplication
# Author: CAFI Team
# Created: 2025-10-24
# License: CC-BY 4.0
# ==============================================================================
# Functions included:
#   DATA LOADING:
#     - load_cafi_data()
#     - load_treatment_data()
#     - load_coral_growth_data()
#
#   PLOTTING THEMES:
#     - theme_cafi()
#     - theme_cafi_pub()
#
#   PLOTTING HELPERS:
#     - save_both()
#     - save_figure()
#
#   COMMUNITY ANALYSIS:
#     - create_community_matrix()
#     - filter_rare_taxa()
#
#   STATISTICAL HELPERS:
#     - fit_reef_lmm()
#
#   UTILITIES:
#     - get_taxonomic_resolution()
#     - safe_view()
#     - strip_fe()
#
#   CONSTANTS:
#     - ALIVE_THRESH (coral survival threshold: 0.80)
# ==============================================================================

# ==============================================================================
# SHARED CONSTANTS
# ==============================================================================

#' Coral survival threshold
#'
#' Corals with tissue survival >= this threshold are included in analyses.
#' Standard value: 80% alive (0.80)
ALIVE_THRESH <- 0.80

# ==============================================================================
# DATA TRANSFORMATION HELPERS
# ==============================================================================

#' Strip FE- prefix from coral IDs
#'
#' Standardizes coral ID format by removing the "FE-" prefix.
#' Raw data uses "FE-POC61" format, processed data uses "POC61".
#'
#' @param x Character vector of coral IDs
#' @return Character vector with FE- prefix removed
#' @examples
#' strip_fe("FE-POC61")  # Returns "POC61"
#' strip_fe(c("FE-POC61", "FE-POC62"))  # Returns c("POC61", "POC62")
strip_fe <- function(x) {
  stringr::str_remove(as.character(x), "^FE-")
}

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

#' Load CAFI community data
#'
#' Loads and preprocesses the main CAFI community dataset with standard
#' transformations and filtering.
#'
#' @param use_cache Logical. If TRUE and cache exists, load from RDS cache
#' @param show_message Logical. Show loading message
#' @return Tibble with CAFI community data
#' @examples
#' cafi_data <- load_cafi_data()
load_cafi_data <- function(use_cache = TRUE, show_message = TRUE) {
  cache_file <- here::here("data", "processed", "cafi_data.rds")
  data_file <- here::here("data", "MRB Amount", "1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv")

  # Check if data file exists
  if (!file.exists(data_file)) {
    stop("CAFI data file not found at: ", data_file, "\n",
         "Please check that the data file exists in the expected location.")
  }

  # Load from cache if available
  if (use_cache && file.exists(cache_file)) {
    if (show_message) message("📦 Loading CAFI data from cache...")
    return(readRDS(cache_file))
  }

  # Load from CSV
  if (show_message) message("📂 Loading CAFI data from CSV...")

  data <- readr::read_csv(data_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      coral_id = as.factor(coral_id),
      site = as.factor(site)
    ) %>%
    dplyr::filter(!is.na(coral_id))

  # Create cache directory if needed
  cache_dir <- dirname(cache_file)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Save cache
  saveRDS(data, cache_file)
  if (show_message) message("✅ CAFI data loaded (", nrow(data), " rows)")

  return(data)
}

#' Load treatment data
#'
#' Loads coral treatment assignments and spatial positions
#'
#' @param show_message Logical. Show loading message
#' @return Tibble with treatment data including reef identifiers
#' @examples
#' treatment_data <- load_treatment_data()
load_treatment_data <- function(show_message = TRUE) {
  data_file <- here::here("data", "MRB Amount", "coral_id_position_treatment.csv")

  if (!file.exists(data_file)) {
    stop("Treatment data file not found at: ", data_file)
  }

  if (show_message) message("📂 Loading treatment data...")

  data <- readr::read_csv(data_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      coral_id = as.factor(coral_id),
      treatment = as.factor(treatment)
    ) %>%
    # Extract reef information from position
    dplyr::mutate(
      row_num = stringr::str_extract(position, "^\\d+"),
      col_num = stringr::str_extract(position, "(?<=-)\\d+"),
      replicate = stringr::str_extract(position, "[A-Za-z]+$"),
      reef = paste0("Reef_", row_num, "-", col_num)
    )

  if (show_message) message("✅ Treatment data loaded (", nrow(data), " rows)")

  return(data)
}

#' Load coral growth data
#'
#' Loads processed coral growth metrics
#'
#' @param show_message Logical. Show loading message
#' @return Tibble with coral growth metrics
#' @examples
#' growth_data <- load_coral_growth_data()
load_coral_growth_data <- function(show_message = TRUE) {
  data_file <- here::here("data", "MRB Amount", "coral_growth_surface_area_change.csv")

  if (!file.exists(data_file)) {
    stop("Coral growth data file not found at: ", data_file)
  }

  if (show_message) message("📂 Loading coral growth data...")

  data <- readr::read_csv(data_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      coral_id = paste0("FE-", coral_id)
    )

  if (show_message) message("✅ Coral growth data loaded (", nrow(data), " rows)")

  return(data)
}

# ==============================================================================
# PLOTTING THEMES
# ==============================================================================

#' Standard CAFI plotting theme
#'
#' Clean, minimal theme for everyday plotting
#'
#' @param base_size Base font size
#' @return A ggplot2 theme object
#' @examples
#' ggplot(data, aes(x, y)) + geom_point() + theme_cafi()
theme_cafi <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}

#' Publication-quality CAFI plotting theme
#'
#' Clean theme suitable for publication with borders and clear axes
#'
#' @param base_size Base font size
#' @return A ggplot2 theme object
#' @examples
#' ggplot(data, aes(x, y)) + geom_point() + theme_cafi_pub()
theme_cafi_pub <- function(base_size = 16) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "black", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# ==============================================================================
# PLOTTING HELPERS
# ==============================================================================

#' Save plot in PNG and PDF formats
#'
#' Convenience function to save a plot in both PNG (high res) and PDF formats
#'
#' @param plot ggplot object to save
#' @param file_stub File path without extension (e.g., "output/figure1")
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution for PNG
#' @param bg Background color
#' @examples
#' p <- ggplot(mtcars, aes(mpg, hp)) + geom_point()
#' save_both(p, "output/scatter", w = 8, h = 6)
save_both <- function(plot, file_stub, width = 8, height = 6, dpi = 600, bg = "white") {
  # Ensure output directory exists
  output_dir <- dirname(file_stub)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save PNG
  ggplot2::ggsave(
    paste0(file_stub, ".png"),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    units = "in"
  )

  # Save PDF
  ggplot2::ggsave(
    paste0(file_stub, ".pdf"),
    plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = bg,
    units = "in"
  )

  message("✅ Saved: ", file_stub, ".png and .pdf")
}

#' Save figure with standard naming
#'
#' Wrapper around ggsave with consistent defaults
#'
#' @param filename Full file path with extension
#' @param plot ggplot object (default: last plot)
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
#' @param bg Background color
#' @examples
#' ggplot(mtcars, aes(mpg, hp)) + geom_point()
#' save_figure("output/scatter.png")
save_figure <- function(filename, plot = ggplot2::last_plot(),
                       width = 8, height = 6, dpi = 300, bg = "white") {
  # Ensure output directory exists
  output_dir <- dirname(filename)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    units = "in"
  )

  message("✅ Saved: ", filename)
}

#' Show and save a plot (backward compatibility alias)
#'
#' Prints plot to screen and saves to file. Wrapper for backward compatibility.
#' Supports both (w, h) and (width, height) argument naming for flexibility.
#'
#' @param plot ggplot object
#' @param filename File path (extension stripped automatically)
#' @param w Width in inches (alias for width)
#' @param h Height in inches (alias for height)
#' @param width Width in inches (overrides w if both provided)
#' @param height Height in inches (overrides h if both provided)
#' @param ... Additional arguments passed to save_figure (excluding width/height)
show_and_save <- function(plot, filename, w = 8, h = 6, width = NULL, height = NULL, ...) {
  # Support both (w, h) and (width, height) naming conventions
  final_width <- if (!is.null(width)) width else w
  final_height <- if (!is.null(height)) height else h
  print(plot)
  save_figure(plot = plot, filename = tools::file_path_sans_ext(filename),
              width = final_width, height = final_height, ...)
}

# ==============================================================================
# COMMUNITY ANALYSIS FUNCTIONS
# ==============================================================================

#' Create community matrix from long-format data
#'
#' Converts long-format species abundance data to a wide community matrix
#' suitable for vegan analyses
#'
#' @param data Data frame with species abundance data
#' @param id_col Column name for site/sample identifier
#' @param species_col Column name for species
#' @param count_col Column name for abundance counts
#' @param include_metadata Return metadata alongside matrix
#' @return List with community matrix and optionally metadata
#' @examples
#' result <- create_community_matrix(cafi_data, "coral_id", "species", "count")
create_community_matrix <- function(data, id_col = "coral_id",
                                   species_col = "species", count_col = "count",
                                   include_metadata = TRUE) {

  # Create matrix
  matrix_data <- data %>%
    dplyr::group_by(.data[[id_col]], .data[[species_col]]) %>%
    dplyr::summarise(total_count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = .data[[species_col]],
      values_from = total_count,
      values_fill = 0
    )

  # Extract metadata if requested
  if (include_metadata) {
    metadata <- matrix_data %>%
      dplyr::select(.data[[id_col]])

    community_matrix <- matrix_data %>%
      dplyr::select(-.data[[id_col]]) %>%
      as.data.frame()

    rownames(community_matrix) <- pull(metadata, .data[[id_col]])

    return(list(
      matrix = community_matrix,
      metadata = metadata,
      species_names = colnames(community_matrix),
      n_sites = nrow(community_matrix),
      n_species = ncol(community_matrix)
    ))
  } else {
    community_matrix <- matrix_data %>%
      dplyr::select(-.data[[id_col]]) %>%
      as.data.frame()

    return(community_matrix)
  }
}

#' Filter rare taxa from community matrix
#'
#' Removes rare species based on minimum abundance and occurrence thresholds
#'
#' @param comm_matrix Community matrix (sites × species)
#' @param min_total Minimum total abundance across all sites
#' @param min_sites Minimum number of sites where species occurs
#' @return Filtered community matrix
#' @examples
#' filtered_comm <- filter_rare_taxa(comm_matrix, min_total = 20, min_sites = 5)
filter_rare_taxa <- function(comm_matrix, min_total = 10, min_sites = 5) {
  col_sums <- colSums(comm_matrix, na.rm = TRUE)
  col_nonzero <- colSums(comm_matrix > 0, na.rm = TRUE)
  keep <- (col_sums >= min_total) & (col_nonzero >= min_sites)

  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message("🔧 Removed ", n_removed, " rare species (",
            round(100 * n_removed / length(keep), 1), "%)")
  }

  return(comm_matrix[, keep, drop = FALSE])
}

# ==============================================================================
# STATISTICAL HELPERS
# ==============================================================================

#' Fit mixed-effects model with reef as random effect
#'
#' Convenience wrapper for fitting linear mixed models with reef-level
#' random effects, a common pattern in the CAFI analyses
#'
#' @param response Character name of response variable
#' @param predictor Character name of predictor variable
#' @param data Data frame
#' @param random_effect Character name of random effect grouping variable
#' @return List with model, ANOVA, tidy summary, and emmeans
#' @examples
#' results <- fit_reef_lmm("growth", "treatment", data, random_effect = "reef")
fit_reef_lmm <- function(response, predictor, data, random_effect = "reef") {
  # Build formula
  formula_str <- paste(response, "~", predictor, "+ (1 |", random_effect, ")")
  model_formula <- as.formula(formula_str)

  # Fit model
  model <- lmerTest::lmer(model_formula, data = data)

  # Return comprehensive results
  list(
    model = model,
    anova = car::Anova(model, type = 3),
    tidy = broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE),
    emmeans = emmeans::emmeans(model, specs = predictor),
    formula = model_formula
  )
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Get finest taxonomic resolution for an observation
#'
#' Checks taxonomic columns in hierarchy from species to phylum and returns
#' the most specific level that contains non-NA data
#'
#' @param row Named list or data frame row with taxonomic columns
#' @return Character: "Species", "Genus", "Family", "Order", "Class", "Phylum", or "Unknown"
#' @examples
#' get_taxonomic_resolution(list(species = "Dascyllus aruanus", genus = "Dascyllus"))
#' # Returns: "Species"
get_taxonomic_resolution <- function(row) {
  resolution_levels <- c("species", "genus", "family", "order", "class", "phylum")

  for (level in resolution_levels) {
    if (level %in% names(row)) {
      val <- row[[level]]
      # Check that it's not NA and not empty string
      if (!is.na(val) && nzchar(as.character(val))) {
        return(stringr::str_to_title(level))
      }
    }
  }

  return("Unknown")
}

#' Safe View wrapper for non-interactive sessions
#'
#' Only opens View() in interactive sessions, otherwise prints head()
#'
#' @param x Object to view
#' @param title Optional title for View window
#' @param n Number of rows to print if non-interactive
#' @examples
#' safe_view(my_data)
safe_view <- function(x, title = deparse(substitute(x)), n = 10) {
  if (interactive()) {
    View(x, title = title)
  } else {
    cat("📋", title, "( first", n, "rows ):\n")
    print(head(x, n))
  }
}

# ==============================================================================
# PACKAGE CHECK HELPER
# ==============================================================================

#' Check if required packages are installed
#'
#' @param packages Character vector of package names
#' @return Invisible TRUE if all installed, stops with message otherwise
check_required_packages <- function(packages) {
  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "), "\n",
         "Install with: install.packages(c(",
         paste(paste0('"', missing, '"'), collapse = ", "), "))")
  }

  invisible(TRUE)
}

# ==============================================================================
# END OF UTILS
# ==============================================================================

message("✅ CAFI utils loaded successfully")
message("📦 Functions available:")
message("  Data: load_cafi_data(), load_treatment_data(), load_coral_growth_data()")
message("  Themes: theme_cafi(), theme_cafi_pub()")
message("  Helpers: save_both(), save_figure(), show_and_save(), create_community_matrix()")
message("  Constants: ALIVE_THRESH, strip_fe()")
