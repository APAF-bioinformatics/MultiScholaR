# Phase 3: Verification Script
# Purpose: Compare functions in func_*.R files against source files to identify gaps

library(stringr)

cat("=== Phase 3: Function Extraction Verification ===\n\n")

# Get all func_*.R files
func_files <- list.files("R", pattern = "^func_.*\\.R$", full.names = TRUE)
cat("Found", length(func_files), "func_*.R files\n")

# Get all mod_*.R files
mod_files <- list.files("R", pattern = "^mod_.*\\.R$", full.names = TRUE)
cat("Found", length(mod_files), "mod_*.R files\n\n")

# Extract all function names from func_*.R and mod_*.R files
extracted_funcs <- character(0)

for (f in c(func_files, mod_files)) {
  content <- readLines(f, warn = FALSE)
  # Match: funcName <- function( or funcName = function(
  matches <- str_match(content, "^([a-zA-Z][a-zA-Z0-9_.]*?)\\s*(<-|=)\\s*function\\s*\\(")
  funcs <- na.omit(matches[,2])
  extracted_funcs <- c(extracted_funcs, funcs)
}

cat("Total functions in func_*.R and mod_*.R files:", length(extracted_funcs), "\n\n")

# Source files to audit
source_files <- c(
  "R/annotation.R",
  "R/de_proteins_functions.R",
  "R/enrichment_functions.R",
  "R/file_management.R",
  "R/helper_functions.R",
  "R/qc_and_rollup.R",
  "R/QC_visualisation.R",
  "R/string_enrichment_functions_refactored.R",
  "R/de_analysis_function_wrapper.R",
  "R/metabolite_de_analysis_wrapper.R",
  "R/metabolite_normalization.R",
  "R/metabolite_qc.R",
  "R/metabolites_de_analysis_wrapper.R",
  "R/multiomics_enrichment_functions.R",
  "R/multiomics_functions_MOFA.R",
  "R/protein_de_analysis_wrapper.R",
  "R/shiny_applets.R",
  "R/get_best_accession_helper.R",
  "R/phosphoproteomics_helper.R",
  "R/allGenerics.R",
  "R/peptideVsSamplesS4Objects.R",
  "R/proteinVsSamplesS4Objects.R",
  "R/metaboliteVsSamplesS4Objects.R",
  "R/functional_enrichment.R",
  "R/limpa_functions.R"
)

# Function to check a source file
check_source <- function(source_file, extracted_funcs) {
  if (!file.exists(source_file)) {
    return(list(file = basename(source_file), exists = FALSE))
  }
  
  content <- readLines(source_file, warn = FALSE)
  
  # Find standalone function definitions
  matches <- str_match(content, "^([a-zA-Z][a-zA-Z0-9_.]*?)\\s*(<-|=)\\s*function\\s*\\(")
  source_funcs <- na.omit(matches[,2])
  
  # Find S4 content
  has_setClass <- any(grepl("^setClass\\(", content))
  has_setMethod <- any(grepl("^setMethod\\(", content))
  has_setGeneric <- any(grepl("^setGeneric\\(", content))
  
  # Count S4 definitions
  n_setClass <- sum(grepl("^setClass\\(", content))
  n_setMethod <- sum(grepl("^setMethod\\(", content))
  n_setGeneric <- sum(grepl("^setGeneric\\(", content))
  
  # Find functions NOT in func_*.R files
  missing <- source_funcs[!source_funcs %in% extracted_funcs]
  already_extracted <- source_funcs[source_funcs %in% extracted_funcs]
  
  list(
    file = basename(source_file),
    exists = TRUE,
    total_funcs = length(source_funcs),
    extracted = length(already_extracted),
    missing = missing,
    n_missing = length(missing),
    has_s4 = has_setClass || has_setMethod || has_setGeneric,
    n_setClass = n_setClass,
    n_setMethod = n_setMethod,
    n_setGeneric = n_setGeneric
  )
}

# Run verification on all source files
results <- list()
for (f in source_files) {
  results[[f]] <- check_source(f, extracted_funcs)
}

# Print results
cat("=== VERIFICATION RESULTS ===\n\n")

# Categorize files
deletable <- character(0)
needs_extraction <- character(0)
s4_files <- character(0)
all_missing_funcs <- list()

for (r in results) {
  if (!r$exists) {
    cat("MISSING:", r$file, "\n")
    next
  }
  
  cat("\n---", r$file, "---\n")
  cat("  Total standalone functions:", r$total_funcs, "\n")
  cat("  Already extracted:", r$extracted, "\n")
  cat("  Missing:", r$n_missing, "\n")
  
  if (r$has_s4) {
    cat("  S4 content: setClass=", r$n_setClass, ", setMethod=", r$n_setMethod, 
        ", setGeneric=", r$n_setGeneric, "\n")
    s4_files <- c(s4_files, r$file)
  }
  
  if (r$n_missing > 0) {
    cat("  Functions to extract:\n")
    for (fn in r$missing) {
      cat("    -", fn, "\n")
    }
    needs_extraction <- c(needs_extraction, r$file)
    all_missing_funcs[[r$file]] <- r$missing
  }
  
  if (r$n_missing == 0 && !r$has_s4 && r$total_funcs > 0) {
    deletable <- c(deletable, r$file)
  }
}

cat("\n\n=== SUMMARY ===\n")
cat("\nFiles ready for deletion (all functions extracted, no S4):\n")
if (length(deletable) > 0) {
  for (f in deletable) cat("  -", f, "\n")
} else {
  cat("  (none)\n")
}

cat("\nFiles needing function extraction:\n")
if (length(needs_extraction) > 0) {
  for (f in needs_extraction) cat("  -", f, "\n")
} else {
  cat("  (none)\n")
}

cat("\nFiles with S4 content (need S4 extraction):\n")
if (length(s4_files) > 0) {
  for (f in s4_files) cat("  -", f, "\n")
} else {
  cat("  (none)\n")
}

# Save detailed results to file
sink("dev/verification_results.txt")
cat("=== PHASE 3 VERIFICATION RESULTS ===\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("TOTAL FUNCTIONS IN func_*.R + mod_*.R:", length(extracted_funcs), "\n\n")

cat("=== MISSING FUNCTIONS BY FILE ===\n\n")
for (name in names(all_missing_funcs)) {
  cat(name, ":\n")
  for (fn in all_missing_funcs[[name]]) {
    cat("  -", fn, "\n")
  }
  cat("\n")
}

cat("=== FILES READY FOR DELETION ===\n")
for (f in deletable) cat(f, "\n")

cat("\n=== FILES WITH S4 CONTENT ===\n")
for (f in s4_files) cat(f, "\n")

sink()

cat("\n\nDetailed results saved to: dev/verification_results.txt\n")

