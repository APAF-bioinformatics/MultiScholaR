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

