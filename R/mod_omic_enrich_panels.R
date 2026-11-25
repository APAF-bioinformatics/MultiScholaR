# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-107 W3 (D2, D3): the enrichment tab's views, moved unchanged with their names from
# mod_prot_enrich_render_helpers.R (the `format*` and `render*` functions, `protEnrichRenderValue`, the download's name
# and archive), and the two reactives and the restore that stay in the GUI: createProtEnrichCurrentAnalysisMethodReactive
# (mod_prot_enrich_execution_helpers.R), createProtEnrichRawContrastNameReactive (mod_prot_enrich_input_helpers.R) and
# restoreProtEnrichResultsFromState (mod_prot_enrich_state_helpers.R). W3 deletes the eight helper files. Names are kept,
# not made neutral: MCI-010, its helper and test-workflow-dia-products.R call them by name (recorded divergence: a shared
# file holds *ProtEnrich* names). The restore's body calls the core (`restoreWorkflowEnrichment`,
# utils_workflow_enrichment.R), and its guard treats an empty holder as nothing restored (F107-02).

#' The label of a significance rule of the clusterProfiler route (KIT-203); NULL (a run made before the rule) is the
#' q-value's.
#' @noRd
enrichmentSignificanceLabel <- function(significance_rule) {
  rules <- ENRICHMENT_SIGNIFICANCE_RULES
  rules$label[match(significance_rule %||% "qvalue", rules$rule)]
}

formatProtEnrichStatusText <- function(analysisComplete,
                                       methodInfo = NULL,
                                       selectedContrast = NULL,
                                       upCutoff = NULL,
                                       downCutoff = NULL,
                                       qCutoff = NULL,
                                       gprofilerResults = NULL,
                                       clusterprofilerResults = NULL,
                                       stringdbResults = NULL,
                                       significanceRule = NULL) {
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
      # KIT-203: the run's significance rule, on the clusterProfiler route only. It is joined to the cutoff's line
      # (`paste` would turn a NULL argument into an empty line), so the g:Profiler text is unchanged.
      paste0(sprintf("Q-value cutoff: %.3f\n", qCutoff), if (identical(methodInfo$method, "clusterprofiler"))
        sprintf("\nSignificance: %s\n", enrichmentSignificanceLabel(significanceRule))),
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

#' The reason a list was not tested (KIT-215, VAL-170): each route has one.
#' @param route `"gprofiler2"` or `"clusterprofiler"`.
#' @param proteins The list's size.
#' @noRd
enrichmentNotTestedReason <- function(route, proteins) {
  ifelse(route == "clusterprofiler",
    sprintf("none of its %d proteins has a Gene Ontology annotation that this analysis tests", proteins),
    sprintf("none of its %d proteins has a gene name, and g:Profiler tests gene names", proteins))
}

#' The summary's lines for the selected comparison's directions that were not tested or had no list (KIT-215): one
#' line per such direction, in the form of the summary's direction lines; none when there are none.
#' @param notTested,noList The tab's `enrichment_not_tested`, `enrichment_no_list` (NULL when absent).
#' @param selectedContrast The selected comparison (display name).
#' @param label `"Up-regulated"` and `"Down-regulated"`, or the g:Profiler summary's `"... pathways"` labels.
#' @return Named character vector (`up`, `down`), only the directions concerned.
#' @noRd
protEnrichDirectionNotes <- function(notTested, noList, selectedContrast, label) {
  notes <- character()
  for (direction in c("up", "down")) {
    hit <- if (NROW(notTested)) notTested[notTested$comparison %in% selectedContrast &
      notTested$direction == direction, , drop = FALSE]
    empty <- NROW(noList) && any(noList$comparison %in% selectedContrast & noList$direction == direction)
    if (NROW(hit)) {
      notes[[direction]] <- sprintf("%s: not tested (%s)", label[[direction]],
        enrichmentNotTestedReason(hit$route[[1L]], hit$proteins[[1L]]))
    } else if (isTRUE(empty)) {
      notes[[direction]] <- sprintf("%s: no list (no protein passed the cutoffs)", label[[direction]])
    }
  }
  notes
}

formatProtEnrichGprofilerSummaryText <- function(gprofilerResults,
                                                 directionFilter = "all",
                                                 notTested = NULL,
                                                 noList = NULL,
                                                 selectedContrast = NULL) {
  # KIT-215: the selected comparison's directions not tested or with no list replace their counts.
  notes <- protEnrichDirectionNotes(notTested, noList, selectedContrast,
    c(up = "Up-regulated pathways", down = "Down-regulated pathways"))
  if (is.null(gprofilerResults) || nrow(gprofilerResults) == 0) {
    return(paste(c("No gprofiler2 results available.", notes), collapse = "\n"))
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
        # KIT-215 (VAL-170 (5)): a filtered direction that was not tested or had no list says so, not a count of 0.
        if (!is.na(notes[directionFilter])) messageText <- notes[[directionFilter]]
    } else {
      totalTerms <- nrow(gprofilerResults)
      positiveTerms <- sum(gprofilerResults$directionality == "positive", na.rm = TRUE)
      negativeTerms <- sum(gprofilerResults$directionality == "negative", na.rm = TRUE)

      messageText <- paste(
        sprintf("Total enrichment terms: %d", totalTerms),
        if (is.na(notes["up"])) sprintf("Up-regulated pathways: %d", positiveTerms) else notes[["up"]],
        if (is.na(notes["down"])) sprintf("Down-regulated pathways: %d", negativeTerms) else notes[["down"]],
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
                                                       directionFilter = "all",
                                                       notTested = NULL,
                                                       noList = NULL,
                                                       selectedContrast = NULL) {
  # KIT-215: the selected comparison's directions not tested or with no list replace their counts.
  notes <- protEnrichDirectionNotes(notTested, noList, selectedContrast,
    c(up = "Up-regulated", down = "Down-regulated"))
  if (is.null(clusterprofilerResults) || nrow(clusterprofilerResults) == 0) {
    return(paste(c("No clusterProfileR results available.", notes), collapse = "\n"))
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
        # KIT-215 (VAL-170 (5)): a filtered direction that was not tested or had no list says so, not a count of 0.
        if (!is.na(notes[directionFilter])) messageText <- notes[[directionFilter]]
    } else {
      totalTerms <- nrow(clusterprofilerResults)
      upTerms <- sum(clusterprofilerResults$directionality == "up", na.rm = TRUE)
      downTerms <- sum(clusterprofilerResults$directionality == "down", na.rm = TRUE)

      messageText <- paste(
        sprintf("Total GO terms: %d", totalTerms),
        if (is.na(notes["up"])) sprintf("Up-regulated: %d", upTerms) else notes[["up"]],
        if (is.na(notes["down"])) sprintf("Down-regulated: %d", downTerms) else notes[["down"]],
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

#' Resolve a render-helper argument that may be a zero-argument accessor (KIT-064, roar BUG-08).
#' @param value A value, or a function of no arguments returning it (a reactive read).
#' @return The value.
#' @noRd
protEnrichRenderValue <- function(value) {
  if (is.function(value)) value() else value
}

renderProtEnrichGprofilerResultsTable <- function(gprofilerResults,
                                                  directionFilter = "all",
                                                  renderDtFn = DT::renderDT,
                                                  datatableFn = DT::datatable,
                                                  formatRoundFn = DT::formatRound,
                                                  catFn = cat,
                                                  gprofilerResultsTable = NULL) {
  pager <- new.env(parent = emptyenv())
  renderDtFn({
    # KIT-064 (roar BUG-08): the arguments may be zero-argument accessors of reactive values; they are read here,
    # inside the render expression, so every render re-evaluates them (a value read in the enclosing frame would be
    # memoised by R's promise on the first render and never update).
    gprofilerResultsNow <- protEnrichRenderValue(gprofilerResults)
    directionFilterNow <- protEnrichRenderValue(directionFilter)
    gprofilerResultsTableNow <- protEnrichRenderValue(gprofilerResultsTable)
    pager$read <- NULL
    if (is.null(gprofilerResultsNow) || nrow(gprofilerResultsNow) == 0) {
      return(datatableFn(data.frame(
        Message = "No gprofiler2 results available. Run analysis first."
      )))
    }

    tryCatch({
      if (isWorkflowTableReference(gprofilerResultsTableNow)) {
        # KIT-024: the saved display table is paged server-side from the store (the direction filter applied in
        # the store read); the widget shows the row count.
        filter <- if (identical(directionFilterNow, "all")) NULL else
          list(column = "directionality", value = if (identical(directionFilterNow, "up")) "positive" else "negative")
        manager <- getDIAActiveManager()
        first <- readDIADisplayTablePage(gprofilerResultsTableNow, filter = filter, limit = 25L, manager = manager)
        pager$rows <- first$rows
        pager$read <- function(offset, limit) readDIADisplayTablePage(gprofilerResultsTableNow, filter = filter,
          offset = offset, limit = limit, manager = manager)$data
        currentResults <- first$data
      } else {
      currentResults <- gprofilerResultsNow

      if (!identical(directionFilterNow, "all") &&
          "directionality" %in% names(currentResults)) {
        directionValue <- if (identical(directionFilterNow, "up")) "positive" else "negative"
        currentResults <- currentResults |>
          dplyr::filter(directionality == directionValue)
      }
      }

      datatableFn(
        currentResults,
        options = if (is.null(pager$read)) list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ) else diaDisplayTableOptions(25L, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
        extensions = "Buttons",
        rownames = FALSE
      ) |>
        formatRoundFn(
          columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(currentResults)),
          digits = 4
        )
    }, error = function(e) {
      catFn(paste("*** ERROR in gprofiler2 results table:", e$message, "\n"))
      datatableFn(data.frame(Message = paste("Error:", e$message)))
    })
  }, server = TRUE, funcFilter = diaDisplayTableFilter(pager))
}

renderProtEnrichGprofilerPlot <- function(analysisComplete,
                                          enrichmentResultsFull,
                                          rawContrast,
                                          directionFilter = "all",
                                          renderPlotlyFn = plotly::renderPlotly,
                                          plotLyFn = plotly::plot_ly,
                                          addTextFn = plotly::add_text,
                                          notTested = NULL) {
  renderPlotlyFn({
    # KIT-064 (roar BUG-08): the arguments may be zero-argument accessors of reactive values; they are read here,
    # inside the render expression, so every render re-evaluates them (a value read in the enclosing frame would be
    # memoised by R's promise on the first render and never update).
    analysisCompleteNow <- protEnrichRenderValue(analysisComplete)
    enrichmentResultsFullNow <- protEnrichRenderValue(enrichmentResultsFull)
    rawContrastNow <- protEnrichRenderValue(rawContrast)
    directionFilterNow <- protEnrichRenderValue(directionFilter)
    buildPlaceholderPlot <- function(messageText) {
      plotLyFn() |>
        addTextFn(
          x = 0.5,
          y = 0.5,
          text = messageText,
          showlegend = FALSE
        )
    }

    if (!analysisCompleteNow) {
      return(buildPlaceholderPlot("Run enrichment analysis first"))
    }

    if (is.null(enrichmentResultsFullNow)) {
      return(buildPlaceholderPlot("No enrichment results"))
    }

    tryCatch({
      # KIT-119: a plot the results hold is shown as stored; otherwise it is built here from the stored result, by
      # the builder the run wrote its saved file with (`resolveEnrichmentPlotly`). "all" shows up, else down.
      directions <- if (identical(directionFilterNow, "all")) c("up", "down") else
        intersect(directionFilterNow, c("up", "down"))
      for (direction in directions) {
        plot <- resolveEnrichmentPlotly(enrichmentResultsFullNow, rawContrastNow, direction)
        if (!is.null(plot)) {
          return(plot)
        }
      }

      # KIT-215 (VAL-170): a direction that was not tested says so, with its reason.
      notTestedNow <- protEnrichRenderValue(notTested)
      if (NROW(notTestedNow) && isTRUE(directionFilterNow %in% c("up", "down"))) {
        hit <- notTestedNow[notTestedNow$contrast %in% rawContrastNow & notTestedNow$direction == directionFilterNow, ,
          drop = FALSE]
        if (NROW(hit)) {
          return(buildPlaceholderPlot(paste0("Not tested: ", enrichmentNotTestedReason(hit$route[[1L]],
            hit$proteins[[1L]]), ". An empty result here does not mean that no term is enriched.")))
        }
      }

      if (identical(directionFilterNow, "up")) {
        return(buildPlaceholderPlot("No up-regulated enrichment data"))
      }

      if (identical(directionFilterNow, "down")) {
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
                                                addTextFn = plotly::add_text,
                                                notTested = NULL) {
  renderPlotlyFn({
    # KIT-064 (roar BUG-08): the arguments may be zero-argument accessors of reactive values; they are read here,
    # inside the render expression, so every render re-evaluates them (a value read in the enclosing frame would be
    # memoised by R's promise on the first render and never update).
    analysisCompleteNow <- protEnrichRenderValue(analysisComplete)
    enrichmentResultsFullNow <- protEnrichRenderValue(enrichmentResultsFull)
    rawContrastNow <- protEnrichRenderValue(rawContrast)
    directionFilterNow <- protEnrichRenderValue(directionFilter)
    buildPlaceholderPlot <- function(messageText) {
      plotLyFn() |>
        addTextFn(
          x = 0.5,
          y = 0.5,
          text = messageText,
          showlegend = FALSE
        )
    }

    if (!analysisCompleteNow) {
      return(buildPlaceholderPlot("Run enrichment analysis first"))
    }

    if (is.null(enrichmentResultsFullNow)) {
      return(buildPlaceholderPlot("No enrichment results"))
    }

    tryCatch({
      # KIT-119: a plot the results hold is shown as stored; otherwise it is built here from the stored result, by
      # the builder the run wrote its saved file with (`resolveEnrichmentPlotly`). "all" shows up, else down.
      directions <- if (identical(directionFilterNow, "all")) c("up", "down") else
        intersect(directionFilterNow, c("up", "down"))
      for (direction in directions) {
        plot <- resolveEnrichmentPlotly(enrichmentResultsFullNow, rawContrastNow, direction)
        if (!is.null(plot)) {
          return(plot)
        }
      }

      # KIT-215 (VAL-170): a direction that was not tested says so, with its reason.
      notTestedNow <- protEnrichRenderValue(notTested)
      if (NROW(notTestedNow) && isTRUE(directionFilterNow %in% c("up", "down"))) {
        hit <- notTestedNow[notTestedNow$contrast %in% rawContrastNow & notTestedNow$direction == directionFilterNow, ,
          drop = FALSE]
        if (NROW(hit)) {
          return(buildPlaceholderPlot(paste0("Not tested: ", enrichmentNotTestedReason(hit$route[[1L]],
            hit$proteins[[1L]]), ". An empty result here does not mean that no term is enriched.")))
        }
      }

      if (identical(directionFilterNow, "up")) {
        return(buildPlaceholderPlot("No up-regulated enrichment data"))
      }

      if (identical(directionFilterNow, "down")) {
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
                                                        catFn = cat,
                                                        clusterprofilerResultsTable = NULL) {
  pager <- new.env(parent = emptyenv())
  renderDtFn({
    # KIT-064 (roar BUG-08): the arguments may be zero-argument accessors of reactive values; they are read here,
    # inside the render expression, so every render re-evaluates them (a value read in the enclosing frame would be
    # memoised by R's promise on the first render and never update).
    clusterprofilerResultsNow <- protEnrichRenderValue(clusterprofilerResults)
    directionFilterNow <- protEnrichRenderValue(directionFilter)
    clusterprofilerResultsTableNow <- protEnrichRenderValue(clusterprofilerResultsTable)
    pager$read <- NULL
    if (is.null(clusterprofilerResultsNow) || nrow(clusterprofilerResultsNow) == 0) {
      return(datatableFn(data.frame(
        Message = "No clusterProfileR results available. Run analysis first."
      )))
    }

    tryCatch({
      if (isWorkflowTableReference(clusterprofilerResultsTableNow)) {
        # KIT-024: server-side pages from the store, the direction filter applied in the read.
        filter <- if (identical(directionFilterNow, "all")) NULL else list(column = "directionality", value = directionFilterNow)
        manager <- getDIAActiveManager()
        first <- readDIADisplayTablePage(clusterprofilerResultsTableNow, filter = filter, limit = 25L, manager = manager)
        pager$rows <- first$rows
        pager$read <- function(offset, limit) readDIADisplayTablePage(clusterprofilerResultsTableNow, filter = filter,
          offset = offset, limit = limit, manager = manager)$data
        currentResults <- first$data
      } else {
      currentResults <- clusterprofilerResultsNow

      if (!identical(directionFilterNow, "all") &&
          "directionality" %in% names(currentResults)) {
        currentResults <- currentResults |>
          dplyr::filter(directionality == directionFilterNow)
      }
      }

      datatableFn(
        currentResults,
        options = if (is.null(pager$read)) list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ) else diaDisplayTableOptions(25L, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
        extensions = "Buttons",
        rownames = FALSE
      ) |>
        formatRoundFn(
          columns = intersect(c("pvalue", "p.adjust", "qvalue", "p_adjust_family"), names(currentResults)),
          digits = 4
        )
    }, error = function(e) {
      catFn(sprintf("*** ERROR in clusterProfileR results table: %s ***\n", e$message))
      datatableFn(data.frame(Message = paste("Error:", e$message)))
    })
  }, server = TRUE, funcFilter = diaDisplayTableFilter(pager))
}

renderProtEnrichStringDbResultsTable <- function(stringdbResults,
                                                 filterSignificant = FALSE,
                                                 enrichmentPValThresh = NULL,
                                                 maxResults = Inf,
                                                 renderDtFn = DT::renderDT,
                                                 datatableFn = DT::datatable) {
  renderDtFn({
    # KIT-064 (roar BUG-08): the arguments may be zero-argument accessors of reactive values; they are read here,
    # inside the render expression, so every render re-evaluates them (a value read in the enclosing frame would be
    # memoised by R's promise on the first render and never update).
    stringdbResultsNow <- protEnrichRenderValue(stringdbResults)
    filterSignificantNow <- protEnrichRenderValue(filterSignificant)
    enrichmentPValThreshNow <- protEnrichRenderValue(enrichmentPValThresh)
    maxResultsNow <- protEnrichRenderValue(maxResults)
    if (is.null(stringdbResultsNow) || nrow(stringdbResultsNow) == 0) {
      return(datatableFn(data.frame(
        Message = "STRING-DB analysis not yet implemented.",
        Note = "This tab will show protein-protein interaction network enrichment results.",
        Status = "Coming soon..."
      )))
    }

    tryCatch({
      currentResults <- stringdbResultsNow

      if (isTRUE(filterSignificantNow)) {
        currentResults <- currentResults |>
          dplyr::filter(p_value < enrichmentPValThreshNow)
      }

      if (nrow(currentResults) > maxResultsNow) {
        currentResults <- currentResults |>
          dplyr::slice_head(n = maxResultsNow)
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

#' The download summary's lines for the selected comparison's directions not tested or with no list (KIT-215): the
#' no-list lines, then a section "## Not tested:"; NULL when there are none (the summary is then as before).
#' @noRd
protEnrichDownloadNotes <- function(notTested, noList, selectedContrast) {
  shown <- c(up = "up-regulated", down = "down-regulated")
  empty <- if (NROW(noList)) noList[noList$comparison %in% selectedContrast, , drop = FALSE]
  hit <- if (NROW(notTested)) notTested[notTested$comparison %in% selectedContrast, , drop = FALSE]
  lines <- character()
  if (NROW(empty)) {
    lines <- sprintf("%s, %s: no list (no protein passed the cutoffs)", empty$comparison, shown[empty$direction])
  }
  if (NROW(hit)) {
    lines <- c(lines, "", "## Not tested:", sprintf("%s, %s: %s", hit$comparison, shown[hit$direction],
      enrichmentNotTestedReason(hit$route, hit$proteins)))
  }
  if (length(lines)) paste(lines, collapse = "\n")
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
                                                  sysTimeFn = Sys.time,
                                                  significanceRule = NULL,
                                                  notTested = NULL,
                                                  noList = NULL) {
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
    stringdbLine <- paste("STRING-DB networks:", if (!is.null(stringdbResults)) nrow(stringdbResults) else 0)
    # KIT-215: the selected comparison's directions with no list, then the lists not tested, after the term counts.
    downloadNotes <- protEnrichDownloadNotes(notTested, noList, selectedContrast)
    if (!is.null(downloadNotes)) stringdbLine <- paste(stringdbLine, downloadNotes, sep = "\n")
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
      # KIT-203 (OD9): on the clusterProfiler method the cutoff is named for both its uses, and the run's rule is
      # given; the g:Profiler summary is unchanged.
      if (identical(methodInfo$method, "clusterprofiler")) paste0("Significance cutoff: ", qCutoff, "\nSignificance: ",
        enrichmentSignificanceLabel(significanceRule)) else paste("Q-value Cutoff:", qCutoff),
      "",
      "## Results Summary:",
      paste("gprofiler2 terms:", if (!is.null(gprofilerResults)) nrow(gprofilerResults) else 0),
      paste("clusterProfileR terms:", if (!is.null(clusterprofilerResults)) nrow(clusterprofilerResults) else 0),
      stringdbLine,
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

createProtEnrichCurrentAnalysisMethodReactive <- function(input,
                                                          enrichmentData,
                                                          supportedOrganismsFn,
                                                          reactiveFn = shiny::reactive,
                                                          reqFn = shiny::req,
                                                          resolveAnalysisMethodFn = resolveProtEnrichAnalysisMethod) {
  reactiveFn({
    reqFn(input$organism_taxid)

    methodState <- resolveAnalysisMethodFn(
      organismTaxid = input$organism_taxid,
      supportedOrganisms = supportedOrganismsFn()
    )

    enrichmentData$analysis_method <- methodState$analysisMethod
    enrichmentData$organism_supported <- methodState$organismSupported

    methodState$methodInfo
  })
}

createProtEnrichRawContrastNameReactive <- function(input,
                                                    workflowData = NULL,
                                                    globalEnv = NULL,
                                                    reactiveFn = shiny::reactive,
                                                    reqFn = shiny::req,
                                                    existsFn = exists,
                                                    getFn = get,
                                                    resolveRawContrastNameFn = resolveProtEnrichRawContrastName) {
  reactiveFn({
    reqFn(input$selected_contrast)

    contrastsTbl <- resolveProtEnrichContrastsTable(workflowData, globalEnv, existsFn, getFn)

    resolveRawContrastNameFn(input$selected_contrast, contrastsTbl)$rawContrastName
  })
}

#' Restore the enrichment results from the completed state (KIT-024): after an Open the tab holds no results although
#' the state carries the KIT-003 `EnrichmentResults`; on first use they are rebuilt into the per-contrast tables (and
#' saved display tables) exactly as the run builds them, by the core (`restoreWorkflowEnrichment`, KIT-107 D3).
#' F107-02: the tab starts with an empty list, so an empty holder counts as nothing restored (the start commit tested
#' for NULL, and the restore never ran after an Open).
#' @param enrichmentData The tab's reactive values.
#' @param workflowData The bus (its `state_manager`).
#' @param config The omic's `enrichment` group (`state_name`); the frame passes its own.
#' @return TRUE when restored, FALSE when nothing was needed or nothing was available.
#' @noRd
restoreProtEnrichResultsFromState <- function(enrichmentData, workflowData,
                                              config = omicGuiConfig("proteomics")$enrichment) {
  if (length(shiny::isolate(enrichmentData$all_enrichment_results))) return(FALSE)
  manager <- if (is.null(workflowData)) NULL else shiny::isolate(workflowData$state_manager)
  restored <- restoreWorkflowEnrichment(manager, config)
  if (is.null(restored)) return(FALSE)
  enrichmentData$all_enrichment_results <- restored$all_enrichment_results
  enrichmentData$enrichment_results_full <- restored$enrichment_results_full
  holdOmicEnrichOptionalFields(enrichmentData, restored) # KIT-215
  enrichmentData$analysis_complete <- TRUE
  TRUE
}
