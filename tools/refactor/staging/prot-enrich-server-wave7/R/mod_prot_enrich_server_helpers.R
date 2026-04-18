registerProtEnrichPlotOutputs <- function(output,
                                          input,
                                          enrichmentData,
                                          rawContrastNameFn,
                                          renderGprofilerPlotFn = renderProtEnrichGprofilerPlot,
                                          renderClusterProfilerPlotFn = renderProtEnrichClusterProfilerPlot) {
  gprofilerPlot <- renderGprofilerPlotFn(
    analysisComplete = enrichmentData$analysis_complete,
    enrichmentResultsFull = enrichmentData$enrichment_results_full,
    rawContrast = rawContrastNameFn(),
    directionFilter = input$gprofiler_direction_filter
  )

  clusterprofilerPlot <- renderClusterProfilerPlotFn(
    analysisComplete = enrichmentData$analysis_complete,
    enrichmentResultsFull = enrichmentData$enrichment_results_full,
    rawContrast = rawContrastNameFn(),
    directionFilter = input$clusterprofiler_direction_filter
  )

  output$gprofiler_plot <- gprofilerPlot
  output$clusterprofiler_plot <- clusterprofilerPlot

  list(
    gprofilerPlot = gprofilerPlot,
    clusterprofilerPlot = clusterprofilerPlot
  )
}

registerProtEnrichStringDbPlotOutput <- function(output,
                                                 renderStringDbPlotFn = renderProtEnrichStringDbPlot) {
  stringdbPlot <- renderStringDbPlotFn()
  output$stringdb_plot <- stringdbPlot
  stringdbPlot
}

registerProtEnrichResultsDownloadHandler <- function(output,
                                                     input,
                                                     enrichmentData,
                                                     currentAnalysisMethodFn,
                                                     downloadHandlerFn = shiny::downloadHandler,
                                                     buildFilenameFn = buildProtEnrichResultsDownloadFilename,
                                                     writeArchiveFn = writeProtEnrichResultsDownloadArchive) {
  downloadHandler <- downloadHandlerFn(
    filename = function() {
      buildFilenameFn(input$selected_contrast)
    },
    content = function(file) {
      writeArchiveFn(
        file = file,
        selectedContrast = input$selected_contrast,
        methodInfo = currentAnalysisMethodFn(),
        organismTaxid = input$organism_taxid,
        upCutoff = input$up_cutoff,
        downCutoff = input$down_cutoff,
        qCutoff = input$q_cutoff,
        gprofilerResults = enrichmentData$gprofiler_results,
        clusterprofilerResults = enrichmentData$clusterprofiler_results,
        stringdbResults = enrichmentData$stringdb_results
      )
    }
  )

  output$download_enrichment_results <- downloadHandler
  downloadHandler
}

