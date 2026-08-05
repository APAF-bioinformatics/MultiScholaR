renderProtNormQcImage <- function(
  filename,
  altText,
  normData,
  proteinQcDir,
  renderImageFn = shiny::renderImage,
  buildQcImagePayloadFn = buildProtNormQcImagePayload
) {
  renderImageFn({
    normData$plot_refresh_trigger

    buildQcImagePayloadFn(proteinQcDir, filename, altText)
  }, deleteFile = FALSE)
}

registerProtNormQcImageOutputs <- function(
  output,
  normData,
  proteinQcDir,
  renderQcImageFn = renderProtNormQcImage
) {
  output$pca_post_filtering <- renderQcImageFn(
    filename = "pre_norm_pca.png",
    altText = "PCA Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$density_post_filtering <- renderQcImageFn(
    filename = "pre_norm_density.png",
    altText = "Density Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$rle_post_filtering <- renderQcImageFn(
    filename = "pre_norm_rle.png",
    altText = "RLE Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$correlation_post_filtering <- renderQcImageFn(
    filename = "pre_norm_correlation.png",
    altText = "Correlation Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$pca_post_normalization <- renderQcImageFn(
    filename = "post_norm_pca.png",
    altText = "PCA Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$density_post_normalization <- renderQcImageFn(
    filename = "post_norm_density.png",
    altText = "Density Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$rle_post_normalization <- renderQcImageFn(
    filename = "post_norm_rle.png",
    altText = "RLE Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$correlation_post_normalization <- renderQcImageFn(
    filename = "post_norm_correlation.png",
    altText = "Correlation Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$pca_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_pca.png",
    altText = "PCA RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$density_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_density.png",
    altText = "Density RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$rle_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_rle.png",
    altText = "RLE RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$correlation_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_correlation.png",
    altText = "Correlation RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )

  invisible(output)
}

renderProtNormPostNormFilteringSummary <- function(
  normData,
  renderPlotFn = shiny::renderPlot,
  gridDrawFn = grid::grid.draw,
  plotNewFn = graphics::plot.new,
  textFn = graphics::text
) {
  renderPlotFn({
    if (!is.null(normData$post_norm_filtering_plot)) {
      gridDrawFn(normData$post_norm_filtering_plot)
    } else {
      plotNewFn()
      textFn(
        0.5,
        0.5,
        "Complete normalization and RUV correction\nto generate filtering summary",
        cex = 1.2
      )
    }
  })
}

checkProtNormMemoryUsage <- function(
  threshold_gb = 8,
  context = "",
  gcFn = gc,
  warningFn = warning,
  messageFn = message
) {
  mem_info <- gcFn()
  mem_used_mb <- sum(mem_info[, 2])
  mem_used_gb <- mem_used_mb / 1024

  if (mem_used_gb > threshold_gb) {
    warningFn(sprintf(
      "*** HIGH MEMORY WARNING [%s]: %.1f GB used (threshold: %.1f GB) ***",
      context,
      mem_used_gb,
      threshold_gb
    ))
    messageFn(sprintf(
      "*** HIGH MEMORY WARNING [%s]: %.1f GB used ***",
      context,
      mem_used_gb
    ))
  } else {
    messageFn(sprintf(
      "*** MEMORY CHECK [%s]: %.1f GB used ***",
      context,
      mem_used_gb
    ))
  }

  invisible(mem_used_gb)
}

createProtNormReactiveState <- function(
  reactiveValuesFn = shiny::reactiveValues
) {
  reactiveValuesFn(
    pre_norm_qc_generated = FALSE,
    normalization_complete = FALSE,
    ruv_complete = FALSE,
    correlation_filtering_complete = FALSE,
    QC_composite_figure = NULL,
    qc_plots = list(
      post_filtering = list(
        pca = NULL,
        density = NULL,
        rle = NULL,
        correlation = NULL
      ),
      post_normalization = list(
        pca = NULL,
        density = NULL,
        rle = NULL,
        correlation = NULL
      ),
      ruv_corrected = list(
        pca = NULL,
        density = NULL,
        rle = NULL,
        correlation = NULL
      )
    ),
    normalized_protein_obj = NULL,
    ruv_normalized_obj = NULL,
    correlation_filtered_obj = NULL,
    best_k = NULL,
    control_genes_index = NULL,
    correlation_vector = NULL,
    correlation_threshold = NULL,
    final_qc_plot = NULL,
    final_filtering_plot = NULL,
    post_norm_filtering_plot = NULL,
    filtering_summary_text = NULL,
    plot_refresh_trigger = 0,
    ruv_optimization_result = NULL
  )
}

renderProtNormFilteringSummaryText <- function(
  normData,
  renderTextFn = shiny::renderText,
  getSummaryTextFn = getProtNormFilteringSummaryText
) {
  renderTextFn({
    getSummaryTextFn(normData$filtering_summary_text)
  })
}

renderProtNormFinalQcPlot <- function(
  normData,
  renderPlotFn = shiny::renderPlot,
  resolveFinalQcRenderStateFn = resolveProtNormFinalQcRenderState,
  gridNewpageFn = grid::grid.newpage,
  pushViewportFn = grid::pushViewport,
  popViewportFn = grid::popViewport,
  viewportFn = grid::viewport,
  gridLayoutFn = grid::grid.layout,
  gridDrawFn = grid::grid.draw,
  ggplotGrobFn = ggplotGrob,
  plotNewFn = graphics::plot.new,
  textFn = graphics::text,
  messageFn = message
) {
  renderPlotFn({
    messageFn("--- DEBUG66 [final_qc_plot]: Rendering ---")

    final_qc_state <- resolveFinalQcRenderStateFn(
      finalQcPlot = normData$final_qc_plot,
      finalFilteringPlot = normData$final_filtering_plot
    )

    if (final_qc_state$mode == "combined") {
      messageFn("   [final_qc_plot] Drawing PCA + Filtering Progression...")
      gridNewpageFn()
      pushViewportFn(
        viewportFn(layout = gridLayoutFn(2, 1, heights = c(0.4, 0.6)))
      )

      pushViewportFn(viewportFn(layout.pos.row = 1))
      gridDrawFn(ggplotGrobFn(final_qc_state$finalQcPlot))
      popViewportFn()

      pushViewportFn(viewportFn(layout.pos.row = 2))
      gridDrawFn(final_qc_state$finalFilteringPlot)
      popViewportFn()

      popViewportFn()
      messageFn("   [final_qc_plot] Render complete.")
    } else if (final_qc_state$mode == "filtering_only") {
      messageFn("   [final_qc_plot] Drawing Filtering Progression ONLY...")
      gridDrawFn(final_qc_state$finalFilteringPlot)
      messageFn("   [final_qc_plot] Render complete.")
    } else if (final_qc_state$mode == "pca_only") {
      messageFn("   [final_qc_plot] Drawing PCA ONLY...")
      final_qc_state$finalQcPlot
    } else {
      plotNewFn()
      textFn(0.5, 0.5, "Apply correlation filter to generate final QC plot", cex = 1.2)
    }
  })
}

renderProtNormRuvCanonicalCorrelationPlot <- function(
  normData,
  renderPlotFn = shiny::renderPlot,
  getCanonicalCorrelationPlotFn = getProtNormRuvCanonicalCorrelationPlot,
  plotNewFn = graphics::plot.new,
  textFn = graphics::text
) {
  renderPlotFn({
    cancor_plot <- getCanonicalCorrelationPlotFn(normData$ruv_optimization_result)

    if (!is.null(cancor_plot)) {
      cancor_plot
    } else {
      plotNewFn()
      textFn(0.5, 0.5, "Run normalization to generate RUV canonical correlation plot", cex = 1.2)
    }
  })
}

renderProtNormRuvOptimizationSummary <- function(
  normData,
  ruvMode,
  groupingVariable,
  renderTextFn = shiny::renderText,
  buildSummaryFn = buildProtNormRuvOptimizationSummary
) {
  renderTextFn({
    buildSummaryFn(
      ruvOptimizationResult = normData$ruv_optimization_result,
      ruvMode = ruvMode,
      groupingVariable = groupingVariable
    )
  })
}

renderProtNormRuvOptimizationTable <- function(
  normData,
  renderDataTableFn = DT::renderDataTable,
  prepareOptimizationResultsTableFn = prepareProtNormOptimizationResultsTable,
  datatableFn = DT::datatable,
  formatStyleFn = DT::formatStyle,
  styleEqualFn = DT::styleEqual
) {
  renderDataTableFn({
    optimization_table <- prepareOptimizationResultsTableFn(normData$ruv_optimization_result)

    datatable_output <- datatableFn(
      optimization_table$results,
      options = list(
        pageLength = 10,
        scrollY = "300px",
        scrollCollapse = TRUE,
        dom = "t"
      ),
      rownames = FALSE
    )

    if (!optimization_table$hasResults) {
      return(datatable_output)
    }

    formatStyleFn(
      datatable_output,
      "percentage",
      target = "row",
      backgroundColor = styleEqualFn(optimization_table$bestPercentage, "#e6f3ff")
    )
  })
}

registerProtNormRenderOutputs <- function(
  output,
  normData,
  ruvMode,
  groupingVariable,
  renderPostNormFilteringSummaryFn = renderProtNormPostNormFilteringSummary,
  renderFilteringSummaryTextFn = renderProtNormFilteringSummaryText,
  renderFinalQcPlotFn = renderProtNormFinalQcPlot,
  renderRuvCanonicalCorrelationPlotFn = renderProtNormRuvCanonicalCorrelationPlot,
  renderRuvOptimizationSummaryFn = renderProtNormRuvOptimizationSummary,
  renderRuvOptimizationTableFn = renderProtNormRuvOptimizationTable
) {
  output$post_norm_filtering_summary <- renderPostNormFilteringSummaryFn(
    normData = normData
  )
  output$filtering_summary_text <- renderFilteringSummaryTextFn(
    normData = normData
  )
  output$final_qc_plot <- renderFinalQcPlotFn(
    normData = normData
  )
  output$ruv_canonical_correlation_plot <- renderRuvCanonicalCorrelationPlotFn(
    normData = normData
  )
  output$ruv_optimization_summary <- renderRuvOptimizationSummaryFn(
    normData = normData,
    ruvMode = ruvMode,
    groupingVariable = groupingVariable
  )
  output$ruv_optimization_table <- renderRuvOptimizationTableFn(
    normData = normData
  )

  invisible(output)
}

