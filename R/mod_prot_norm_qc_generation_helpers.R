resolveProtNormQcStateObject <- function(stateManager, reqFn = shiny::req, messageFn = message) {
  reqFn(stateManager)

  current_state <- stateManager$current_state
  messageFn(sprintf("Current state: '%s'", current_state))

  snapshot <- stateManager$.context$snapshot
  if (current_state == "protein_replicate_filtered" || is.null(snapshot$qc_sources$pre_normalization)) {
    snapshot$qc_sources$pre_normalization <- getDIAStateReference(stateManager, current_state)
    publishDIAStateHead(stateManager$.context, snapshot)
  }
  current_s4 <- snapshot$qc_sources$pre_normalization
  messageFn(sprintf("S4 object is NULL: %s", is.null(current_s4)))

  if (!is.null(current_s4)) {
    messageFn(sprintf("S4 object class: %s", class(current_s4)))
  }

  if (is.null(current_s4)) {
    stop("No S4 object available for QC plot generation")
  }

  messageFn("Step 2: S4 object retrieved successfully")

  list(currentState = current_state, currentS4 = current_s4)
}

#' Preserve established injected plot seams while using products for real S4 data.
#' @param s4Object,aesthetics,plotType Product inputs.
#' @param plotPcaFn,plotRleFn,buildDensityFn,buildCorrelationFn Existing producers.
#' @param buildProductFn Strict real-S4 product producer.
#' @return Existing plot value.
#' @noRd
resolveProtNormQcPlotValue <- function(
  s4Object,
  aesthetics,
  plotType,
  plotPcaFn,
  plotRleFn,
  buildDensityFn,
  buildCorrelationFn,
  buildProductFn
) {
  if (isSavedWorkflowReference(s4Object)) {
    return(resolveDIAReferenceQcPlot(s4Object, aesthetics, plotType, buildProductFn,
      plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn))
  }
  if (methods::is(s4Object, "ProteinQuantitativeData")) {
    return(buildProductFn(s4Object, aesthetics, plotType,
      plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn)$plot)
  }
  switch(plotType,
    pca = plotPcaFn(s4Object, grouping_variable = aesthetics$color_var,
      label_column = "", shape_variable = aesthetics$shape_var,
      title = "", font_size = 8),
    rle = plotRleFn(s4Object, group = aesthetics$color_var,
      yaxis_limit = c(-6, 6)),
    density = buildDensityFn(s4Object, aesthetics),
    correlation = buildCorrelationFn(s4Object, aesthetics$color_var)
  )
}

generateProtNormPreNormalizationQcArtifacts <- function(
  stateManager,
  qcDir,
  aesthetics,
  qcPlotPaths = NULL,
  reqFn = shiny::req,
  messageFn = message,
  gcFn = gc,
  timeFn = Sys.time,
  plotPcaFn = plotPca,
  plotRleFn = plotRle,
  buildDensityFn = buildProtNormDensityPlot,
  buildCorrelationFn = buildProtNormCorrelationPlot,
  saveArtifactFn = saveProtNormQcPlotArtifact,
  recordPathFn = recordProtNormQcPlotPath,
  initPathsFn = initializeProtNormQcPlotPaths,
  buildProductFn = buildProtNormQcPlotProduct
) {
  messageFn("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")
  messageFn("Step 1: Getting current S4 object...")

  qc_state <- resolveProtNormQcStateObject(
    stateManager = stateManager,
    reqFn = reqFn,
    messageFn = messageFn
  )
  current_s4 <- qc_state$currentS4

  qcPlotPaths <- initPathsFn(qcPlotPaths)

  messageFn("*** PRE-NORM QC: Generating PCA plot (Step 1/4) ***")
  pca_start_time <- timeFn()
  pca_plot <- resolveProtNormQcPlotValue(current_s4, aesthetics, "pca",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)

  pca_path <- saveArtifactFn(qcDir, "pre_norm_pca.png", pca_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "pca", pca_path)

  rm(pca_plot)
  gcFn()

  pca_elapsed <- as.numeric(difftime(timeFn(), pca_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: PCA plot completed in %.1f seconds ***", pca_elapsed))

  messageFn("*** PRE-NORM QC: Generating RLE plot (Step 2/4) ***")
  rle_start_time <- timeFn()
  rle_plot <- resolveProtNormQcPlotValue(current_s4, aesthetics, "rle",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)

  rle_path <- saveArtifactFn(qcDir, "pre_norm_rle.png", rle_plot, width = 10, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "rle", rle_path)

  rm(rle_plot)
  gcFn()

  rle_elapsed <- as.numeric(difftime(timeFn(), rle_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: RLE plot completed in %.1f seconds ***", rle_elapsed))

  messageFn("*** PRE-NORM QC: Generating Density plot (Step 3/4) ***")
  density_start_time <- timeFn()

  density_plot <- resolveProtNormQcPlotValue(current_s4, aesthetics, "density",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)
  density_path <- saveArtifactFn(qcDir, "pre_norm_density.png", density_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "density", density_path)

  rm(density_plot)
  gcFn()

  density_elapsed <- as.numeric(difftime(timeFn(), density_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: Density plot completed in %.1f seconds ***", density_elapsed))

  messageFn("*** PRE-NORM QC: Generating Pearson correlation plot (Step 4/4) ***")
  num_samples <- if (isSavedWorkflowReference(current_s4)) {
    meta <- readDIAStateMetadata(stateManager, current_s4)$data
    length(setdiff(meta$dimensions$columns, meta$slots$protein_id_column))
  } else length(setdiff(colnames(current_s4@protein_quant_table), current_s4@protein_id_column))
  estimated_pairs <- choose(num_samples, 2)
  messageFn(sprintf("*** PRE-NORM QC: Sample count = %d, Expected pairs ~= %d ***", num_samples, estimated_pairs))

  pearson_start_time <- timeFn()

  correlation_plot <- resolveProtNormQcPlotValue(current_s4, aesthetics, "correlation",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)
  corr_path <- saveArtifactFn(qcDir, "pre_norm_correlation.png", correlation_plot, width = 10, height = 8)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "correlation", corr_path)

  rm(correlation_plot)
  gcFn()

  pearson_elapsed <- as.numeric(difftime(timeFn(), pearson_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: Pearson correlation completed in %.1f seconds ***", pearson_elapsed))

  messageFn("*** PRE-NORM QC: Running garbage collection ***")
  gcFn()
  messageFn("Pre-normalization QC plots generated successfully")

  qcPlotPaths
}

generateProtNormPostNormalizationQcArtifacts <- function(
  normalizedS4,
  qcDir,
  aesthetics,
  qcPlotPaths = NULL,
  messageFn = message,
  gcFn = gc,
  plotPcaFn = plotPca,
  plotRleFn = plotRle,
  buildDensityFn = buildProtNormDensityPlot,
  buildCorrelationFn = buildProtNormCorrelationPlot,
  saveArtifactFn = saveProtNormQcPlotArtifact,
  recordPathFn = recordProtNormQcPlotPath,
  initPathsFn = initializeProtNormQcPlotPaths,
  buildProductFn = buildProtNormQcPlotProduct
) {
  messageFn("=== GENERATING POST-NORMALIZATION QC PLOTS ===")

  qcPlotPaths <- initPathsFn(qcPlotPaths)

  messageFn("*** POST-NORM QC: Generating PCA plot ***")
  pca_plot <- resolveProtNormQcPlotValue(normalizedS4, aesthetics, "pca",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)

  pca_path <- saveArtifactFn(qcDir, "post_norm_pca.png", pca_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "pca", pca_path)

  rm(pca_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Generating RLE plot ***")
  rle_plot <- resolveProtNormQcPlotValue(normalizedS4, aesthetics, "rle",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)

  rle_path <- saveArtifactFn(qcDir, "post_norm_rle.png", rle_plot, width = 10, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "rle", rle_path)

  rm(rle_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Generating Density plot ***")
  density_plot <- resolveProtNormQcPlotValue(normalizedS4, aesthetics, "density",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)
  density_path <- saveArtifactFn(qcDir, "post_norm_density.png", density_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "density", density_path)

  rm(density_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Generating Pearson correlation plot ***")
  correlation_plot <- resolveProtNormQcPlotValue(normalizedS4, aesthetics, "correlation",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)
  corr_path <- saveArtifactFn(qcDir, "post_norm_correlation.png", correlation_plot, width = 10, height = 8)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "correlation", corr_path)

  rm(correlation_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Running garbage collection ***")
  gcFn()
  messageFn("Post-normalization QC plots generated successfully")

  qcPlotPaths
}

generateProtNormRuvCorrectedQcArtifacts <- function(
  ruvCorrectedS4,
  qcDir,
  aesthetics,
  qcPlotPaths = NULL,
  messageFn = message,
  gcFn = gc,
  plotPcaFn = plotPca,
  plotRleFn = plotRle,
  buildDensityFn = buildProtNormDensityPlot,
  buildCorrelationFn = buildProtNormCorrelationPlot,
  saveArtifactFn = saveProtNormQcPlotArtifact,
  recordPathFn = recordProtNormQcPlotPath,
  initPathsFn = initializeProtNormQcPlotPaths,
  buildProductFn = buildProtNormQcPlotProduct
) {
  messageFn("=== GENERATING RUV-CORRECTED QC PLOTS ===")

  qcPlotPaths <- initPathsFn(qcPlotPaths)

  messageFn("*** RUV QC: Generating PCA plot ***")
  pca_plot <- resolveProtNormQcPlotValue(ruvCorrectedS4, aesthetics, "pca",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)

  pca_path <- saveArtifactFn(qcDir, "ruv_corrected_pca.png", pca_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "pca", pca_path)

  rm(pca_plot)
  gcFn()

  messageFn("*** RUV QC: Generating RLE plot ***")
  rle_plot <- resolveProtNormQcPlotValue(ruvCorrectedS4, aesthetics, "rle",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)

  rle_path <- saveArtifactFn(qcDir, "ruv_corrected_rle.png", rle_plot, width = 10, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "rle", rle_path)

  rm(rle_plot)
  gcFn()

  messageFn("*** RUV QC: Generating Density plot ***")
  density_plot <- resolveProtNormQcPlotValue(ruvCorrectedS4, aesthetics, "density",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)
  density_path <- saveArtifactFn(qcDir, "ruv_corrected_density.png", density_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "density", density_path)

  rm(density_plot)
  gcFn()

  messageFn("*** RUV QC: Generating Pearson correlation plot ***")
  correlation_plot <- resolveProtNormQcPlotValue(ruvCorrectedS4, aesthetics, "correlation",
    plotPcaFn, plotRleFn, buildDensityFn, buildCorrelationFn, buildProductFn)
  corr_path <- saveArtifactFn(qcDir, "ruv_corrected_correlation.png", correlation_plot, width = 10, height = 8)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "correlation", corr_path)

  rm(correlation_plot)
  gcFn()

  messageFn("*** RUV QC: Running garbage collection ***")
  gcFn()
  messageFn("RUV-corrected QC plots generated successfully")

  qcPlotPaths
}
