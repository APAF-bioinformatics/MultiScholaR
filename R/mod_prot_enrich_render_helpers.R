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

renderProtEnrichGprofilerResultsTable <- function(gprofilerResults,
                                                  directionFilter = "all",
                                                  renderDtFn = DT::renderDT,
                                                  datatableFn = DT::datatable,
                                                  formatRoundFn = DT::formatRound,
                                                  catFn = cat) {
  renderDtFn({
    if (is.null(gprofilerResults) || nrow(gprofilerResults) == 0) {
      return(datatableFn(data.frame(
        Message = "No gprofiler2 results available. Run analysis first."
      )))
    }

    tryCatch({
      currentResults <- gprofilerResults

      if (!identical(directionFilter, "all") &&
          "directionality" %in% names(currentResults)) {
        directionValue <- if (identical(directionFilter, "up")) "positive" else "negative"
        currentResults <- currentResults |>
          dplyr::filter(directionality == directionValue)
      }

      datatableFn(
        currentResults,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons"
      ) |>
        formatRoundFn(
          columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(currentResults)),
          digits = 4
        )
    }, error = function(e) {
      catFn(paste("*** ERROR in gprofiler2 results table:", e$message, "\n"))
      datatableFn(data.frame(Message = paste("Error:", e$message)))
    })
  })
}

renderProtEnrichGprofilerPlot <- function(analysisComplete,
                                          enrichmentResultsFull,
                                          rawContrast,
                                          directionFilter = "all",
                                          renderPlotlyFn = plotly::renderPlotly,
                                          plotLyFn = plotly::plot_ly,
                                          addTextFn = plotly::add_text) {
  renderPlotlyFn({
    buildPlaceholderPlot <- function(messageText) {
      plotLyFn() |>
        addTextFn(
          x = 0.5,
          y = 0.5,
          text = messageText,
          showlegend = FALSE
        )
    }

    if (!analysisComplete) {
      return(buildPlaceholderPlot("Run enrichment analysis first"))
    }

    if (is.null(enrichmentResultsFull)) {
      return(buildPlaceholderPlot("No enrichment results"))
    }

    tryCatch({
      plots <- protEnrichInteractivePlots(
        enrichmentResultsFull,
        rawContrast
      )

      if (identical(directionFilter, "up") && !is.null(plots$up)) {
        return(plots$up)
      }

      if (identical(directionFilter, "down") && !is.null(plots$down)) {
        return(plots$down)
      }

      if (identical(directionFilter, "all")) {
        if (!is.null(plots$up)) {
          return(plots$up)
        }

        if (!is.null(plots$down)) {
          return(plots$down)
        }
      }

      if (identical(directionFilter, "up")) {
        return(buildPlaceholderPlot("No up-regulated enrichment data"))
      }

      if (identical(directionFilter, "down")) {
        return(buildPlaceholderPlot("No down-regulated enrichment data"))
      }

      buildPlaceholderPlot("No plot data for this contrast")
    }, error = function(e) {
      buildPlaceholderPlot(paste("Plot error:", e$message))
    })
  })
}

renderProtEnrichClusterProfilerPlot <- function(analysisComplete,
                                                enrichmentResultsFull,
                                                rawContrast,
                                                directionFilter = "all",
                                                renderPlotlyFn = plotly::renderPlotly,
                                                plotLyFn = plotly::plot_ly,
                                                addTextFn = plotly::add_text) {
  renderPlotlyFn({
    buildPlaceholderPlot <- function(messageText) {
      plotLyFn() |>
        addTextFn(
          x = 0.5,
          y = 0.5,
          text = messageText,
          showlegend = FALSE
        )
    }

    if (!analysisComplete) {
      return(buildPlaceholderPlot("Run enrichment analysis first"))
    }

    if (is.null(enrichmentResultsFull)) {
      return(buildPlaceholderPlot("No enrichment results"))
    }

    tryCatch({
      plots <- protEnrichInteractivePlots(
        enrichmentResultsFull,
        rawContrast
      )

      if (identical(directionFilter, "up") && !is.null(plots$up)) {
        return(plots$up)
      }

      if (identical(directionFilter, "down") && !is.null(plots$down)) {
        return(plots$down)
      }

      if (identical(directionFilter, "all")) {
        if (!is.null(plots$up)) {
          return(plots$up)
        }

        if (!is.null(plots$down)) {
          return(plots$down)
        }
      }

      if (identical(directionFilter, "up")) {
        return(buildPlaceholderPlot("No up-regulated enrichment data"))
      }

      if (identical(directionFilter, "down")) {
        return(buildPlaceholderPlot("No down-regulated enrichment data"))
      }

      buildPlaceholderPlot("No plot data for this contrast")
    }, error = function(e) {
      buildPlaceholderPlot(paste("Plot error:", e$message))
    })
  })
}

renderProtEnrichClusterProfilerResultsTable <- function(clusterprofilerResults,
                                                        directionFilter = "all",
                                                        renderDtFn = DT::renderDT,
                                                        datatableFn = DT::datatable,
                                                        formatRoundFn = DT::formatRound,
                                                        catFn = cat) {
  renderDtFn({
    if (is.null(clusterprofilerResults) || nrow(clusterprofilerResults) == 0) {
      return(datatableFn(data.frame(
        Message = "No clusterProfileR results available. Run analysis first."
      )))
    }

    tryCatch({
      currentResults <- clusterprofilerResults

      if (!identical(directionFilter, "all") &&
          "directionality" %in% names(currentResults)) {
        currentResults <- currentResults |>
          dplyr::filter(directionality == directionFilter)
      }

      datatableFn(
        currentResults,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons"
      ) |>
        formatRoundFn(
          columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(currentResults)),
          digits = 4
        )
    }, error = function(e) {
      catFn(sprintf("*** ERROR in clusterProfileR results table: %s ***\n", e$message))
      datatableFn(data.frame(Message = paste("Error:", e$message)))
    })
  })
}

renderProtEnrichStringDbResultsTable <- function(stringdbResults,
                                                 filterSignificant = FALSE,
                                                 enrichmentPValThresh = NULL,
                                                 maxResults = Inf,
                                                 renderDtFn = DT::renderDT,
                                                 datatableFn = DT::datatable) {
  renderDtFn({
    if (is.null(stringdbResults) || nrow(stringdbResults) == 0) {
      return(datatableFn(data.frame(
        Message = "STRING-DB analysis not yet implemented.",
        Note = "This tab will show protein-protein interaction network enrichment results.",
        Status = "Coming soon..."
      )))
    }

    tryCatch({
      currentResults <- stringdbResults

      if (isTRUE(filterSignificant)) {
        currentResults <- currentResults |>
          dplyr::filter(p_value < enrichmentPValThresh)
      }

      if (nrow(currentResults) > maxResults) {
        currentResults <- currentResults |>
          dplyr::slice_head(n = maxResults)
      }

      datatableFn(
        currentResults,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons"
      )
    }, error = function(e) {
      datatableFn(data.frame(Message = paste("STRING-DB Error:", e$message)))
    })
  })
}

renderProtEnrichStringDbPlot <- function(renderPlotlyFn = plotly::renderPlotly,
                                         plotLyFn = plotly::plot_ly,
                                         addTextFn = plotly::add_text) {
  renderPlotlyFn({
    plotLyFn() |>
      addTextFn(
        x = 0.5,
        y = 0.5,
        text = "STRING-DB coming soon",
        showlegend = FALSE
      )
  })
}

buildProtEnrichResultsDownloadFilename <- function(selectedContrast,
                                                   date = Sys.Date(),
                                                   prefix = "Enrichment_results") {
  contrastSafe <- gsub("[^A-Za-z0-9_.-]", "_", selectedContrast)
  paste0(prefix, "_", contrastSafe, "_", date, ".zip")
}

writeProtEnrichResultsDownloadArchive <- function(file,
                                                  selectedContrast,
                                                  methodInfo,
                                                  organismTaxid,
                                                  upCutoff,
                                                  downCutoff,
                                                  qCutoff,
                                                  gprofilerResults = NULL,
                                                  clusterprofilerResults = NULL,
                                                  stringdbResults = NULL,
                                                  tempDir = tempdir(),
                                                  writeTsvFn = readr::write_tsv,
                                                  writeLinesFn = writeLines,
                                                  zipFn = utils::zip,
                                                  sysTimeFn = Sys.time) {
  tryCatch({
    filesToZip <- character()

    if (!is.null(gprofilerResults) && nrow(gprofilerResults) > 0) {
      gprofilerFile <- file.path(tempDir, "gprofiler2_results.tsv")
      writeTsvFn(gprofilerResults, gprofilerFile)
      filesToZip <- c(filesToZip, gprofilerFile)
    }

    if (!is.null(clusterprofilerResults) && nrow(clusterprofilerResults) > 0) {
      clusterprofilerFile <- file.path(tempDir, "clusterProfileR_results.tsv")
      writeTsvFn(clusterprofilerResults, clusterprofilerFile)
      filesToZip <- c(filesToZip, clusterprofilerFile)
    }

    if (!is.null(stringdbResults) && nrow(stringdbResults) > 0) {
      stringdbFile <- file.path(tempDir, "stringdb_results.tsv")
      writeTsvFn(stringdbResults, stringdbFile)
      filesToZip <- c(filesToZip, stringdbFile)
    }

    summaryFile <- file.path(tempDir, "enrichment_analysis_summary.txt")
    summaryContent <- paste(
      "# MultiScholaR Enrichment Analysis Results",
      paste("Date:", sysTimeFn()),
      paste("Contrast:", selectedContrast),
      paste("Analysis Method:", methodInfo$method),
      paste("Organism:", methodInfo$species_name),
      paste("Taxonomy ID:", organismTaxid),
      "",
      "## Analysis Parameters:",
      paste("Up log2FC Cutoff:", upCutoff),
      paste("Down log2FC Cutoff:", downCutoff),
      paste("Q-value Cutoff:", qCutoff),
      "",
      "## Results Summary:",
      paste("gprofiler2 terms:", if (!is.null(gprofilerResults)) nrow(gprofilerResults) else 0),
      paste("clusterProfileR terms:", if (!is.null(clusterprofilerResults)) nrow(clusterprofilerResults) else 0),
      paste("STRING-DB networks:", if (!is.null(stringdbResults)) nrow(stringdbResults) else 0),
      "",
      "## Files Included:",
      if (length(filesToZip) > 0) paste("*", basename(filesToZip), collapse = "\n") else "* No result files (analysis may have failed)",
      "",
      "Generated by MultiScholaR Enrichment Analysis Module",
      sep = "\n"
    )

    writeLinesFn(summaryContent, summaryFile)
    filesToZip <- c(filesToZip, summaryFile)

    if (length(filesToZip) > 0) {
      zipFn(zipfile = file, files = filesToZip, flags = "-j")
    } else {
      noteFile <- file.path(tempDir, "no_results.txt")
      writeLinesFn("No enrichment results available to download.", noteFile)
      zipFn(zipfile = file, files = noteFile, flags = "-j")
      filesToZip <- noteFile
    }

    invisible(list(
      status = "ok",
      files = filesToZip,
      zipfile = file
    ))
  }, error = function(e) {
    errorFile <- file.path(tempDir, "download_error.txt")
    writeLinesFn(paste("Error creating download:", e$message), errorFile)
    zipFn(zipfile = file, files = errorFile, flags = "-j")

    invisible(list(
      status = "error",
      files = errorFile,
      zipfile = file,
      error = e$message
    ))
  })
}
