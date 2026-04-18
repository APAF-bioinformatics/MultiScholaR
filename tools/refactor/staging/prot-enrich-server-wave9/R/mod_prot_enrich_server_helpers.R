registerProtEnrichGprofilerResultsTableOutput <- function(output,
                                                          input,
                                                          enrichmentData,
                                                          renderResultsTableFn = renderProtEnrichGprofilerResultsTable) {
  gprofilerResultsTable <- renderResultsTableFn(
    gprofilerResults = enrichmentData$gprofiler_results,
    directionFilter = input$gprofiler_direction_filter
  )

  output$gprofiler_results_table <- gprofilerResultsTable
  gprofilerResultsTable
}

registerProtEnrichGprofilerSummaryOutput <- function(output,
                                                     input,
                                                     enrichmentData,
                                                     renderTextFn = shiny::renderText,
                                                     formatSummaryTextFn = formatProtEnrichGprofilerSummaryText) {
  gprofilerSummaryStats <- renderTextFn({
    formatSummaryTextFn(
      gprofilerResults = enrichmentData$gprofiler_results,
      directionFilter = input$gprofiler_direction_filter
    )
  })

  output$gprofiler_summary_stats <- gprofilerSummaryStats
  gprofilerSummaryStats
}

registerProtEnrichClusterProfilerResultsTableOutput <- function(output,
                                                                input,
                                                                enrichmentData,
                                                                renderResultsTableFn = renderProtEnrichClusterProfilerResultsTable) {
  clusterprofilerResultsTable <- renderResultsTableFn(
    clusterprofilerResults = enrichmentData$clusterprofiler_results,
    directionFilter = input$clusterprofiler_direction_filter
  )

  output$clusterprofiler_results_table <- clusterprofilerResultsTable
  clusterprofilerResultsTable
}

registerProtEnrichClusterProfilerSummaryOutput <- function(output,
                                                           input,
                                                           enrichmentData,
                                                           renderTextFn = shiny::renderText,
                                                           formatSummaryTextFn = formatProtEnrichClusterProfilerSummaryText) {
  clusterprofilerSummaryStats <- renderTextFn({
    formatSummaryTextFn(
      clusterprofilerResults = enrichmentData$clusterprofiler_results,
      directionFilter = input$clusterprofiler_direction_filter
    )
  })

  output$clusterprofiler_summary_stats <- clusterprofilerSummaryStats
  clusterprofilerSummaryStats
}

registerProtEnrichStringDbResultsTableOutput <- function(output,
                                                         input,
                                                         enrichmentData,
                                                         renderResultsTableFn = renderProtEnrichStringDbResultsTable) {
  stringdbResultsTable <- renderResultsTableFn(
    stringdbResults = enrichmentData$stringdb_results,
    filterSignificant = input$stringdb_filter_significant,
    enrichmentPValThresh = input$enrichment_p_val_thresh,
    maxResults = input$stringdb_max_results
  )

  output$stringdb_results_table <- stringdbResultsTable
  stringdbResultsTable
}

registerProtEnrichStringDbSummaryOutput <- function(output,
                                                    enrichmentData,
                                                    renderTextFn = shiny::renderText,
                                                    formatSummaryTextFn = formatProtEnrichStringDbSummaryText) {
  stringdbSummaryStats <- renderTextFn({
    formatSummaryTextFn(
      stringdbResults = enrichmentData$stringdb_results
    )
  })

  output$stringdb_summary_stats <- stringdbSummaryStats
  stringdbSummaryStats
}

registerProtEnrichStatusOutput <- function(output,
                                           input,
                                           enrichmentData,
                                           currentAnalysisMethodFn,
                                           renderTextFn = shiny::renderText,
                                           formatStatusTextFn = formatProtEnrichStatusText) {
  enrichmentStatus <- renderTextFn({
    methodInfo <- NULL
    if (enrichmentData$analysis_complete) {
      methodInfo <- currentAnalysisMethodFn()
    }

    formatStatusTextFn(
      analysisComplete = enrichmentData$analysis_complete,
      methodInfo = methodInfo,
      selectedContrast = input$selected_contrast,
      upCutoff = input$up_cutoff,
      downCutoff = input$down_cutoff,
      qCutoff = input$q_cutoff,
      gprofilerResults = enrichmentData$gprofiler_results,
      clusterprofilerResults = enrichmentData$clusterprofiler_results,
      stringdbResults = enrichmentData$stringdb_results
    )
  })

  output$enrichment_status <- enrichmentStatus
  enrichmentStatus
}

