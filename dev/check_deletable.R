# Check which source files can be safely deleted (no S4 content)

source_files <- c(
  "annotation.R",
  "enrichment_functions.R", 
  "file_management.R",
  "multiomics_functions_MOFA.R",
  "de_proteins_functions.R",
  "string_enrichment_functions_refactored.R",
  "de_analysis_function_wrapper.R",
  "multiomics_enrichment_functions.R",
  "protein_de_analysis_wrapper.R",
  "shiny_applets.R",
  "get_best_accession_helper.R",
  "phosphoproteomics_helper.R",
  "helper_functions.R",
  "qc_and_rollup.R",
  "QC_visualisation.R",
  "metabolite_de_analysis_wrapper.R",
  "metabolite_normalization.R",
  "metabolite_qc.R",
  "metabolites_de_analysis_wrapper.R"
)

deletable <- character(0)
has_s4 <- character(0)

for (f in source_files) {
  path <- file.path("R", f)
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    has_setClass <- any(grepl("^setClass|<- setClass", lines))
    has_setMethod <- any(grepl("^setMethod", lines))
    
    if (!has_setClass && !has_setMethod) {
      cat("DELETABLE:", f, "\n")
      deletable <- c(deletable, f)
    } else {
      cat("HAS S4:", f, "(setClass:", has_setClass, ", setMethod:", has_setMethod, ")\n")
      has_s4 <- c(has_s4, f)
    }
  } else {
    cat("NOT FOUND:", f, "\n")
  }
}

cat("\n=== SUMMARY ===\n")
cat("\nFiles safe to delete (", length(deletable), "):\n")
for (f in deletable) cat("  - R/", f, "\n", sep = "")

cat("\nFiles with S4 content (KEEP):\n")
for (f in has_s4) cat("  - R/", f, "\n", sep = "")

