# MultiScholaR: Baseline Checkpoint Generator for Normalization Refactor
# Trace ID: tr_20260317_225113_f
# Purpose: Generate Metabolomics and Lipidomics S4 objects at the 'post-QC' stage
#          to serve as fixtures for the mod_norm refactoring verification.

library(MultiScholaR)
library(dplyr)
library(logger)

# Source updated code directly to use new method signatures
source("R/allGenerics.R")
source("R/func_metab_s4_objects.R")
source("R/func_lipid_s4_objects.R")

# --- 1. Setup ---
# Setup dummy project directories
temp_root <- tempdir()
project_dirs <- list(
  fixtures = file.path(getwd(), "tests", "fixtures")
)
if (!dir.exists(project_dirs$fixtures)) dir.create(project_dirs$fixtures, recursive = TRUE)

log_info("Starting baseline fixture generation...")

# --- 2. Metabolomics Baseline ---
log_info("Generating Metabolomics fixtures (Mock Dataset)...")
metab_pos_path <- "/Users/ignatiuspang/Dropbox/2026/APAF/MultiScholaR/MultiScholaR_test_data/Metabolomics/metabolomics_msdial_positive.tsv"
metab_neg_path <- "/Users/ignatiuspang/Dropbox/2026/APAF/MultiScholaR/MultiScholaR_test_data/Metabolomics/metabolomics_msdial_negative.tsv"

# Load data (TSV)
metab_pos <- read.delim(metab_pos_path, stringsAsFactors = FALSE, check.names = FALSE)
metab_neg <- read.delim(metab_neg_path, stringsAsFactors = FALSE, check.names = FALSE)

# Create dummy design matrix based on column names (Control_1...6, Treatment_1...6)
metab_samples <- grep("Control|Treatment", names(metab_pos), value = TRUE)
metab_design <- data.frame(
  Run = metab_samples,
  group = ifelse(grepl("Control", metab_samples), "Control", "Treatment"),
  stringsAsFactors = FALSE
)

# Create S4 object with both modes
metab_s4 <- createMetaboliteAssayData(
  metabolite_data = list(positive = metab_pos, negative = metab_neg),
  design_matrix = metab_design,
  metabolite_id_column = "Metabolite name",
  sample_id = "Run",
  group_id = "group"
)

# Apply basic QC (Group-Aware Intensity Filtering)
# Keep metabolite if above 1 percentile in at least 3 samples of each group for at least 2 groups
metab_post_filter <- metaboliteIntensityFiltering(
  metab_s4,
  grouping_variable = "group",
  min_samples_per_group = 3,
  min_groups = 2,
  metabolites_intensity_cutoff_percentile = 1
)

# Save fixture
saveRDS(metab_post_filter, file.path(project_dirs$fixtures, "metab_baseline_post_filter.rds"))
log_info("Metabolomics fixture saved.")

# --- 3. Lipidomics Baseline ---
log_info("Generating Lipidomics fixtures (Mock Dataset)...")
lipid_pos_path <- "/Users/ignatiuspang/Dropbox/2026/APAF/MultiScholaR/MultiScholaR_test_data/Lipidomics/lipidomics_lipidsearch_positive_mode.csv"
lipid_neg_path <- "/Users/ignatiuspang/Dropbox/2026/APAF/MultiScholaR/MultiScholaR_test_data/Lipidomics/lipidomics_lipidsearch_negative_mode.csv"

# Load data (CSV)
lipid_pos <- read.csv(lipid_pos_path, stringsAsFactors = FALSE, check.names = FALSE)
lipid_neg <- read.csv(lipid_neg_path, stringsAsFactors = FALSE, check.names = FALSE)

# Detect raw sample columns for lipids (exclude _Norm)
lipid_samples <- grep("Control|Treatment", names(lipid_pos), value = TRUE)
lipid_samples <- lipid_samples[!grepl("_Norm$", lipid_samples)]

lipid_design <- data.frame(
  Run = lipid_samples,
  group = ifelse(grepl("Control", lipid_samples), "Control", "Treatment"),
  stringsAsFactors = FALSE
)

# Create S4 object
lipid_s4 <- createLipidomicsAssayData(
  lipid_data = list(positive = lipid_pos, negative = lipid_neg),
  design_matrix = lipid_design,
  lipid_id_column = "LipidName",
  sample_id = "Run",
  group_id = "group"
)

# Apply basic QC
lipid_post_filter <- lipidIntensityFiltering(
  lipid_s4,
  grouping_variable = "group",
  min_samples_per_group = 3,
  min_groups = 2,
  lipids_intensity_cutoff_percentile = 1
)

# Save fixture
saveRDS(lipid_post_filter, file.path(project_dirs$fixtures, "lipid_baseline_post_filter.rds"))
log_info("Lipidomics fixture saved.")

log_info("Baseline generation complete.")

<!-- APAF Bioinformatics | generate_norm_checkpoints.R | Approved | 2026-03-17 -->
