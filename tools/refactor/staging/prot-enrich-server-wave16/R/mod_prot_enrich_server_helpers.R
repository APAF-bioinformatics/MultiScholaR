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
      plots <- enrichmentResultsFull@enrichment_plotly[[rawContrast]]

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
      plots <- enrichmentResultsFull@enrichment_plotly[[rawContrast]]

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

