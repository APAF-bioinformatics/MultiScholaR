formatProtEnrichStatusText <- function(analysisComplete,
                                       methodInfo = NULL,
                                       selectedContrast = NULL,
                                       upCutoff = NULL,
                                       downCutoff = NULL,
                                       qCutoff = NULL,
                                       gprofilerResults = NULL,
                                       clusterprofilerResults = NULL,
                                       stringdbResults = NULL) {
  if (isTRUE(analysisComplete)) {
    gprofilerCount <- if (!is.null(gprofilerResults)) nrow(gprofilerResults) else 0
    clusterprofilerCount <- if (!is.null(clusterprofilerResults)) nrow(clusterprofilerResults) else 0
    stringdbCount <- if (!is.null(stringdbResults)) nrow(stringdbResults) else 0

    return(paste(
      "[OK] Analysis Complete\n",
      sprintf("Method: %s\n", methodInfo$method),
      sprintf("Contrast: %s\n", selectedContrast),
      sprintf("Up log2FC cutoff: %.1f\n", upCutoff),
      sprintf("Down log2FC cutoff: %.1f\n", downCutoff),
      sprintf("Q-value cutoff: %.3f\n", qCutoff),
      sprintf("Organism: %s\n", methodInfo$species_name),
      "",
      "Results Available:",
      sprintf("* gprofiler2: %d terms", gprofilerCount),
      sprintf("* clusterProfileR: %d terms", clusterprofilerCount),
      sprintf("* STRING-DB: %d networks", stringdbCount),
      "",
      "[OK] Results saved to workflow state",
      sep = "\n"
    ))
  }

  paste(
    "[WAITING] Ready for analysis\n",
    "",
    "Steps:",
    "1. Select contrast from DE results",
    "2. Set log fold change cutoffs",
    "3. Set Q-value cutoff (significance threshold)",
    "4. Click 'Run Enrichment Analysis'",
    "",
    "Method automatically determined by organism.",
    sep = "\n"
  )
}

formatProtEnrichAnalysisMethodText <- function(methodInfo) {
  if (isTRUE(methodInfo$supported)) {
    return(paste(
      "[OK] SUPPORTED ORGANISM\n",
      sprintf("Method: %s\n", methodInfo$method),
      sprintf("Species: %s\n", methodInfo$species_name),
      "All enrichment methods available"
    ))
  }

  paste(
    "[WARNING] CUSTOM ORGANISM\n",
    sprintf("Method: %s\n", methodInfo$method),
    sprintf("Organism: %s\n", methodInfo$species_name),
    "Using UniProt GO annotations"
  )
}

formatProtEnrichContrastsText <- function(contrastsAvailable) {
  if (!is.null(contrastsAvailable)) {
    return(paste(contrastsAvailable, collapse = "\n"))
  }

  "No contrasts available.\nComplete differential expression\nanalysis first."
}

