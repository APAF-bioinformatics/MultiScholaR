
library(MultiScholaR)
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)

# Source necessary files
source("R/func_metab_s4_objects.R")
source("R/func_metab_de.R")

# Mock assay data
assay_pos <- data.frame(
  `Alignment ID` = c("M1", "M2", "M3"),
  `Metabolite name` = c("Metab1", "Metab2", "Metab3"),
  S1 = c(10, 11, 12),
  S2 = c(10.5, 11.5, 12.5),
  S3 = c(14, 15, 16),
  S4 = c(14.5, 15.5, 16.5),
  check.names = FALSE
)

# Mock design matrix
design <- data.frame(
  Run = c("S1", "S2", "S3", "S4"),
  group = c("Control", "Control", "Treatment", "Treatment"),
  stringsAsFactors = FALSE
)

# Create object
obj <- new("MetaboliteAssayData",
  metabolite_data = list(LCMS_Pos = assay_pos),
  design_matrix = design,
  sample_id = "Run",
  group_id = "group",
  metabolite_id_column = "Alignment ID",
  annotation_id_column = "Metabolite name"
)

# Mock contrasts table
contrasts_tbl <- data.frame(
  contrasts = "groupTreatment-groupControl",
  friendly_names = "Treatment_vs_Control",
  stringsAsFactors = FALSE
)

message("\n--- Running runMetabolitesDE ---")
tryCatch({
  results <- runMetabolitesDE(
    theObject = obj,
    contrasts_tbl = contrasts_tbl,
    formula_string = "~ 0 + group"
  )
  message("DE Analysis Successful!")
}, error = function(e) {
  message("DE Analysis Failed with error:\n", e$message)
})
