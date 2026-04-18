formatProtEnrichGprofilerSummaryText <- function(gprofilerResults,
                                                 directionFilter = "all") {
  if (is.null(gprofilerResults) || nrow(gprofilerResults) == 0) {
    return("No gprofiler2 results available.")
  }

  tryCatch({
    if (!identical(directionFilter, "all") &&
        "directionality" %in% names(gprofilerResults)) {
      directionValue <- if (identical(directionFilter, "up")) "positive" else "negative"
      filteredResults <- gprofilerResults |> dplyr::filter(directionality == directionValue)
      displayedCount <- nrow(filteredResults)

      if (identical(directionFilter, "up")) {
        messageText <- sprintf("Showing %d up-regulated pathways", displayedCount)
      } else {
        messageText <- sprintf("Showing %d down-regulated pathways", displayedCount)
      }
    } else {
      totalTerms <- nrow(gprofilerResults)
      positiveTerms <- sum(gprofilerResults$directionality == "positive", na.rm = TRUE)
      negativeTerms <- sum(gprofilerResults$directionality == "negative", na.rm = TRUE)

      messageText <- paste(
        sprintf("Total enrichment terms: %d", totalTerms),
        sprintf("Up-regulated pathways: %d", positiveTerms),
        sprintf("Down-regulated pathways: %d", negativeTerms),
        sep = "\n"
      )
    }

    paste(
      messageText,
      "",
      "Results displayed in table below.",
      sep = "\n"
    )
  }, error = function(e) {
    paste("Error calculating statistics:", e$message)
  })
}

formatProtEnrichClusterProfilerSummaryText <- function(clusterprofilerResults,
                                                       directionFilter = "all") {
  if (is.null(clusterprofilerResults) || nrow(clusterprofilerResults) == 0) {
    return("No clusterProfileR results available.")
  }

  tryCatch({
    if (!identical(directionFilter, "all") &&
        "directionality" %in% names(clusterprofilerResults)) {
      filteredResults <- clusterprofilerResults |>
        dplyr::filter(directionality == directionFilter)
      displayedCount <- nrow(filteredResults)

      if (identical(directionFilter, "up")) {
        messageText <- sprintf("Showing %d up-regulated GO terms", displayedCount)
      } else {
        messageText <- sprintf("Showing %d down-regulated GO terms", displayedCount)
      }
    } else {
      totalTerms <- nrow(clusterprofilerResults)
      upTerms <- sum(clusterprofilerResults$directionality == "up", na.rm = TRUE)
      downTerms <- sum(clusterprofilerResults$directionality == "down", na.rm = TRUE)

      messageText <- paste(
        sprintf("Total GO terms: %d", totalTerms),
        sprintf("Up-regulated: %d", upTerms),
        sprintf("Down-regulated: %d", downTerms),
        sep = "\n"
      )
    }

    paste(
      messageText,
      "",
      "Results displayed in table below.",
      sep = "\n"
    )
  }, error = function(e) {
    paste("Error calculating statistics:", e$message)
  })
}

formatProtEnrichStringDbSummaryText <- function(stringdbResults) {
  if (is.null(stringdbResults) || nrow(stringdbResults) == 0) {
    return("STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics.")
  }

  paste(
    "STRING-DB Network Analysis",
    "Status: Implementation pending",
    "",
    "Features planned:",
    "* Protein-protein interaction networks",
    "* Functional cluster identification",
    "* Network topology analysis",
    "* Interactive network visualization",
    sep = "\n"
  )
}

