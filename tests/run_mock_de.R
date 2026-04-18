if (requireNamespace("MultiScholaR", quietly = TRUE)) {
  library(MultiScholaR)
} else {
  tryCatch(
    devtools::load_all("."),
    error = function(e) {
      stop(
        paste(
          "Unable to load MultiScholaR from the local checkout.",
          "Install the package or install the missing package imports required by devtools::load_all().",
          paste("Original error:", conditionMessage(e))
        ),
        call. = FALSE
      )
    }
  )
}

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

# Create object using the package API rather than hard-coded source files
obj <- createMetaboliteAssayData(
  metabolite_data = list(LCMS_Pos = assay_pos),
  design_matrix = design,
  sample_id = "Run",
  group_id = "group",
  metabolite_id_column = "Alignment ID",
  annotation_id_column = "Metabolite name",
  database_identifier_type = "InternalName"
)

# Mock contrasts table
contrasts_tbl <- data.frame(
  contrasts = "groupTreatment-groupControl",
  friendly_names = "Treatment_vs_Control",
  stringsAsFactors = FALSE
)

message("\n--- Running runMetabolitesDA ---")
tryCatch({
  results <- runMetabolitesDA(
    theObject = obj,
    contrasts_tbl = contrasts_tbl,
    formula_string = "~ 0 + group"
  )
  message("DA Analysis Successful!")
  message("Result entries: ", paste(names(results), collapse = ", "))
}, error = function(e) {
  message("DA Analysis Failed with error:\n", e$message)
})
