#!/usr/bin/env Rscript
# ==============================================================================
# Script: validate_pipeline.R
# Purpose: Validate data integrity and verify pipeline outputs
# Author: Repository maintenance
# Date: 2026-01-14
# ==============================================================================
#
# Usage:
#   Rscript scripts/MRB/validate_pipeline.R
#
# This script replaces the Python validation agent with pure R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cli)
  library(digest)
})

# ---- Functions ----

#' Print colored status messages
print_status <- function(msg, type = "info") {
  if (type == "success") {
    cli_alert_success(msg)
  } else if (type == "error") {
    cli_alert_danger(msg)
  } else if (type == "warning") {
    cli_alert_warning(msg)
  } else {
    cli_alert_info(msg)
  }
}

#' Verify checksums for data files
verify_checksums <- function(checksum_file = here("data/checksums.txt")) {
  if (!file.exists(checksum_file)) {
    print_status("Checksum file not found - skipping verification", "warning")
    return(invisible(FALSE))
  }

  print_status("Verifying data file checksums...")

  # Read checksums
  checksums <- read_delim(checksum_file, delim = "  ",
                          col_names = c("hash", "file"),
                          show_col_types = FALSE)

  data_dir <- here("data/MRB Amount")
  all_valid <- TRUE

  for (i in 1:nrow(checksums)) {
    file_path <- file.path(data_dir, checksums$file[i])
    expected_hash <- checksums$hash[i]

    if (!file.exists(file_path)) {
      print_status(paste("Missing:", checksums$file[i]), "error")
      all_valid <- FALSE
      next
    }

    actual_hash <- digest(file_path, algo = "md5", file = TRUE)

    if (actual_hash == expected_hash) {
      print_status(paste("✓", checksums$file[i]), "success")
    } else {
      print_status(paste("✗", checksums$file[i], "(checksum mismatch)"), "error")
      all_valid <- FALSE
    }
  }

  return(all_valid)
}

#' Check required data files exist
check_data_files <- function() {
  print_status("Checking required data files...")

  required_files <- c(
    "data/MRB Amount/1. mrb_fe_cafi_summer_2021_v4_AP_updated_2024.csv",
    "data/MRB Amount/coral_id_position_treatment.csv",
    "data/MRB Amount/1. amount_master_phys_data_v5.csv",
    "data/MRB Amount/MRB_2019_200K_mesh_measure.csv",
    "data/MRB Amount/MRB_May_2021_200K_mesh_measure.csv"
  )

  all_present <- TRUE
  for (file in required_files) {
    file_path <- here(file)
    if (file.exists(file_path)) {
      print_status(paste("✓", basename(file)), "success")
    } else {
      print_status(paste("✗ Missing:", basename(file)), "error")
      all_present <- FALSE
    }
  }

  return(all_present)
}

#' Verify key output files exist
check_outputs <- function() {
  print_status("Checking pipeline outputs...")

  expected_outputs <- c(
    "output/MRB/tables/MANUSCRIPT_STATISTICAL_TESTS.csv",
    "output/MRB/figures/coral/growth_treatment_comparison.png",
    "output/MRB/objects/sessionInfo_05_MRB_comm.txt"
  )

  all_present <- TRUE
  for (file in expected_outputs) {
    file_path <- here(file)
    if (file.exists(file_path)) {
      file_info <- file.info(file_path)
      size_kb <- round(file_info$size / 1024, 1)
      print_status(paste("✓", basename(file), paste0("(", size_kb, " KB)")), "success")
    } else {
      print_status(paste("✗ Missing:", basename(file)), "warning")
      all_present <- FALSE
    }
  }

  return(all_present)
}

#' Get repository status
get_status <- function() {
  cli_h1("Repository Status")

  # Check data
  cli_h2("Data Files")
  data_ok <- check_data_files()

  # Check checksums
  cli_h2("Data Integrity")
  checksums_ok <- verify_checksums()

  # Check outputs
  cli_h2("Pipeline Outputs")
  outputs_ok <- check_outputs()

  # Summary
  cli_h2("Summary")
  if (data_ok && checksums_ok && outputs_ok) {
    print_status("All checks passed! Repository is ready.", "success")
    return(0)
  } else {
    if (!data_ok) print_status("Some required data files are missing", "warning")
    if (!checksums_ok) print_status("Some checksums don't match", "warning")
    if (!outputs_ok) print_status("Some outputs missing (run pipeline)", "warning")
    return(1)
  }
}

#' Count files in output directories
count_outputs <- function() {
  cli_h1("Output Statistics")

  dirs <- c(
    "output/MRB/figures" = "Figures",
    "output/MRB/tables" = "Tables",
    "output/MRB/objects" = "R Objects"
  )

  for (dir_path in names(dirs)) {
    full_path <- here(dir_path)
    if (dir.exists(full_path)) {
      files <- list.files(full_path, recursive = TRUE, pattern = "\\.(png|pdf|csv|html|txt|rds)$")
      print_status(paste(dirs[dir_path], ":", length(files), "files"), "info")
    }
  }
}

# ---- Main Execution ----

main <- function() {
  cli_rule("CAFI Pipeline Validation")
  cat("\n")

  exit_code <- get_status()
  cat("\n")
  count_outputs()

  cat("\n")
  cli_rule("Validation Complete")

  quit(status = exit_code)
}

# Run if called as script
if (!interactive()) {
  main()
}
