# ============================================================================
# validate_dependencies.R
# ============================================================================
# Purpose: Validate that all function dependencies in mod_*.R and app files
#          are available in func_*.R files before deleting source files.
#
# Usage:
#   source("dev/validate_dependencies.R")
#   validate_app_dependencies()
#
# ============================================================================

library(stringr)

#' Extract all function calls from R code
#' 
#' @param code Character vector of R code lines
#' @return Character vector of function names called
extract_function_calls <- function(code) {
  # Pattern for function calls: word followed by (
  # Exclude common R keywords and control structures
  all_matches <- str_extract_all(code, "[a-zA-Z][a-zA-Z0-9_.]*(?=\\s*\\()")
  calls <- unique(unlist(all_matches))
  
  # Exclude R keywords and common non-function patterns
  keywords <- c("if", "else", "for", "while", "repeat", "function", "in",
                "next", "break", "TRUE", "FALSE", "NULL", "NA", "NaN", "Inf",
                "return", "library", "require", "c", "list", "data.frame",
                "matrix", "vector", "array", "factor", "names", "class",
                "length", "nrow", "ncol", "dim", "is", "as", "print", "cat",
                "paste", "paste0", "sprintf", "message", "warning", "stop",
                "tryCatch", "try", "switch", "match", "which", "any", "all",
                "sum", "mean", "max", "min", "abs", "sqrt", "log", "exp",
                "round", "floor", "ceiling", "seq", "rep", "unique", "sort",
                "order", "rev", "head", "tail", "subset", "merge", "rbind",
                "cbind", "lapply", "sapply", "vapply", "mapply", "apply",
                "do.call", "Reduce", "Filter", "Map", "with", "within",
                "exists", "get", "assign", "rm", "ls", "new", "setClass",
                "setMethod", "setGeneric", "setValidity", "validObject",
                "ns", "NS", "session", "input", "output", "reactive",
                "reactiveVal", "reactiveValues", "observe", "observeEvent",
                "eventReactive", "isolate", "req", "validate", "need",
                "renderUI", "renderText", "renderPrint", "renderPlot",
                "renderTable", "renderDataTable", "renderDT", "uiOutput",
                "textOutput", "plotOutput", "tableOutput", "dataTableOutput",
                "DTOutput", "actionButton", "actionLink", "textInput",
                "numericInput", "selectInput", "checkboxInput", "radioButtons",
                "fileInput", "downloadButton", "downloadHandler", "tags",
                "tagList", "div", "span", "p", "h1", "h2", "h3", "h4", "h5",
                "br", "hr", "icon", "HTML", "includeHTML", "includeCSS",
                "fluidPage", "fluidRow", "column", "wellPanel", "tabPanel",
                "navbarPage", "sidebarLayout", "sidebarPanel", "mainPanel",
                "conditionalPanel", "showModal", "modalDialog", "removeModal",
                "showNotification", "removeNotification", "updateSelectInput",
                "updateNumericInput", "updateTextInput", "updateCheckboxInput",
                "updateRadioButtons", "updateTabsetPanel", "insertUI",
                "removeUI", "runApp", "shinyApp", "moduleServer", "callModule"
                )
  
  calls <- calls[!calls %in% keywords]
  calls <- calls[!grepl("^(shiny|shinydashboard|DT|plotly|ggplot2|dplyr|tidyr|purrr|stringr|logger|golem)::", calls)]
  
  return(calls)
}

#' Get all defined functions in a file
#'
#' @param file_path Path to R file
#' @return Character vector of function names defined
get_defined_functions <- function(file_path) {
  if (!file.exists(file_path)) return(character(0))
  
  content <- readLines(file_path, warn = FALSE)
  
  # Match function definitions
  matches <- str_match(content, "^([a-zA-Z][a-zA-Z0-9_.]*?)\\s*(<-|=)\\s*function\\s*\\(")
  funcs <- na.omit(matches[, 2])
  
  return(as.character(funcs))
}

#' Validate all module dependencies
#'
#' @return List with validation results
#' @export
validate_app_dependencies <- function() {
  cat("=== MultiScholaR App Dependency Validation ===\n\n")
  
  # Get all func_*.R defined functions
  func_files <- list.files("R", pattern = "^func_.*\\.R$", full.names = TRUE)
  func_defined <- character(0)
  for (f in func_files) {
    func_defined <- c(func_defined, get_defined_functions(f))
  }
  cat("Functions defined in func_*.R files:", length(unique(func_defined)), "\n")
  
  # Get all mod_*.R defined functions
  mod_files <- list.files("R", pattern = "^mod_.*\\.R$", full.names = TRUE)
  mod_defined <- character(0)
  for (f in mod_files) {
    mod_defined <- c(mod_defined, get_defined_functions(f))
  }
  cat("Functions defined in mod_*.R files:", length(unique(mod_defined)), "\n")
  
  # Get app infrastructure functions
  app_files <- c("R/app_server.R", "R/app_ui.R", "R/run_app.R", "R/launch_app.R")
  app_defined <- character(0)
  for (f in app_files) {
    if (file.exists(f)) {
      app_defined <- c(app_defined, get_defined_functions(f))
    }
  }
  cat("Functions defined in app infrastructure:", length(unique(app_defined)), "\n")
  
  # Get utils functions
  util_files <- list.files("R", pattern = "^utils_.*\\.R$", full.names = TRUE)
  util_defined <- character(0)
  for (f in util_files) {
    util_defined <- c(util_defined, get_defined_functions(f))
  }
  cat("Functions defined in utils_*.R files:", length(unique(util_defined)), "\n\n")
  
  # ALL available functions (in new architecture)
  all_available <- unique(c(func_defined, mod_defined, app_defined, util_defined))
  
  # Source files that will be deleted
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
    "R/phosphoproteomics_helper.R"
  )
  
  # Get functions ONLY in source files (would be lost)
  source_only <- character(0)
  for (f in source_files) {
    if (file.exists(f)) {
      source_only <- c(source_only, get_defined_functions(f))
    }
  }
  source_only <- unique(source_only)
  
  # Functions that would be lost (in source but not in func_*)
  would_be_lost <- source_only[!source_only %in% all_available]
  
  cat("=== FUNCTIONS THAT WOULD BE LOST IF SOURCE FILES DELETED ===\n")
  if (length(would_be_lost) > 0) {
    for (fn in sort(would_be_lost)) {
      cat("  -", fn, "\n")
    }
    cat("\nTotal:", length(would_be_lost), "functions need extraction before deletion!\n")
  } else {
    cat("  NONE - All source file functions are available in func_*.R\n")
  }
  
  # Now check what mod_*.R files actually CALL
  cat("\n=== CHECKING MOD_*.R FUNCTION DEPENDENCIES ===\n\n")
  
  all_mod_calls <- character(0)
  mod_dependencies <- list()
  
  for (mod_file in c(mod_files, app_files)) {
    if (!file.exists(mod_file)) next
    
    content <- readLines(mod_file, warn = FALSE)
    calls <- extract_function_calls(content)
    
    # Filter to non-base, non-package prefixed calls
    relevant_calls <- calls[!grepl("::", calls)]
    relevant_calls <- relevant_calls[nchar(relevant_calls) > 2]  # Skip very short names
    
    mod_dependencies[[basename(mod_file)]] <- relevant_calls
    all_mod_calls <- c(all_mod_calls, relevant_calls)
  }
  
  all_mod_calls <- unique(all_mod_calls)
  
  # Check which called functions are NOT in func_*.R or mod_*.R
  missing_calls <- all_mod_calls[!all_mod_calls %in% all_available]
  
  # Further filter - check if they're in source files
  all_source_funcs <- character(0)
  for (f in source_files) {
    if (file.exists(f)) {
      all_source_funcs <- c(all_source_funcs, get_defined_functions(f))
    }
  }
  
  # S4 files (will keep)
  s4_files <- c(
    "R/allGenerics.R",
    "R/peptideVsSamplesS4Objects.R",
    "R/proteinVsSamplesS4Objects.R",
    "R/metaboliteVsSamplesS4Objects.R",
    "R/functional_enrichment.R",
    "R/limpa_functions.R"
  )
  
  s4_funcs <- character(0)
  for (f in s4_files) {
    if (file.exists(f)) {
      s4_funcs <- c(s4_funcs, get_defined_functions(f))
    }
  }
  
  # Functions called by modules that are ONLY in to-be-deleted source files
  critical_missing <- missing_calls[missing_calls %in% source_only]
  critical_missing <- critical_missing[!critical_missing %in% s4_funcs]
  
  cat("Functions called by mod_*.R/app files:\n")
  cat("  Total unique calls:", length(all_mod_calls), "\n")
  cat("  Available in func_*/mod_*/utils_*:", sum(all_mod_calls %in% all_available), "\n")
  cat("  Available in S4 files (keeping):", sum(all_mod_calls %in% s4_funcs), "\n")
  cat("  CRITICAL - In source files only:", length(critical_missing), "\n")
  
  if (length(critical_missing) > 0) {
    cat("\n=== CRITICAL: THESE FUNCTIONS MUST BE EXTRACTED ===\n")
    for (fn in sort(critical_missing)) {
      # Find which source file contains it
      for (f in source_files) {
        if (file.exists(f)) {
          if (fn %in% get_defined_functions(f)) {
            cat("  -", fn, "  (from", basename(f), ")\n")
            break
          }
        }
      }
    }
  }
  
  # Save detailed report
  sink("dev/dependency_validation_report.txt")
  cat("=== MULTISCHOLAR APP DEPENDENCY VALIDATION REPORT ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("=== SUMMARY ===\n")
  cat("Functions in func_*.R:", length(unique(func_defined)), "\n")
  cat("Functions in mod_*.R:", length(unique(mod_defined)), "\n")
  cat("Functions in utils_*.R:", length(unique(util_defined)), "\n")
  cat("Functions in app infrastructure:", length(unique(app_defined)), "\n")
  cat("Functions in S4 files (keeping):", length(unique(s4_funcs)), "\n")
  cat("Functions in source files (to delete):", length(unique(source_only)), "\n")
  cat("\n")
  
  cat("=== FUNCTIONS THAT WOULD BE LOST ===\n")
  for (fn in sort(would_be_lost)) cat(fn, "\n")
  
  cat("\n=== CRITICAL DEPENDENCIES (called by mods, only in source) ===\n")
  for (fn in sort(critical_missing)) cat(fn, "\n")
  
  cat("\n=== PER-MODULE DEPENDENCIES ===\n")
  for (mod_name in names(mod_dependencies)) {
    calls <- mod_dependencies[[mod_name]]
    critical <- calls[calls %in% critical_missing]
    if (length(critical) > 0) {
      cat("\n", mod_name, ":\n")
      for (fn in critical) cat("  -", fn, "\n")
    }
  }
  
  sink()
  
  cat("\n\nDetailed report saved to: dev/dependency_validation_report.txt\n")
  
  # Return validation status
  return(list(
    valid = length(critical_missing) == 0 && length(would_be_lost) == 0,
    would_be_lost = would_be_lost,
    critical_missing = critical_missing,
    all_available = all_available
  ))
}

# Run validation
cat("\n")
result <- validate_app_dependencies()

if (result$valid) {
  cat("\n*** VALIDATION PASSED - Safe to delete source files ***\n")
} else {
  cat("\n*** VALIDATION FAILED - Extract missing functions first ***\n")
}

