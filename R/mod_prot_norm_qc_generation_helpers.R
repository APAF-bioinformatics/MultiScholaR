resolveProtNormQcStateObject <- function(stateManager, reqFn = shiny::req, messageFn = message) {
  reqFn(stateManager)

  current_state <- stateManager$current_state
  messageFn(sprintf("Current state: '%s'", current_state))

  current_s4 <- stateManager$getState(current_state)
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
  initPathsFn = initializeProtNormQcPlotPaths
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
  pca_plot <- plotPcaFn(
    current_s4,
    grouping_variable = aesthetics$color_var,
    label_column = "",
    shape_variable = aesthetics$shape_var,
    title = "",
    font_size = 8
  )

  pca_path <- saveArtifactFn(qcDir, "pre_norm_pca.png", pca_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "pca", pca_path)

  rm(pca_plot)
  gcFn()

  pca_elapsed <- as.numeric(difftime(timeFn(), pca_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: PCA plot completed in %.1f seconds ***", pca_elapsed))

  messageFn("*** PRE-NORM QC: Generating RLE plot (Step 2/4) ***")
  rle_start_time <- timeFn()
  rle_plot <- plotRleFn(
    current_s4,
    group = aesthetics$color_var,
    yaxis_limit = c(-6, 6)
  )

  rle_path <- saveArtifactFn(qcDir, "pre_norm_rle.png", rle_plot, width = 10, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "rle", rle_path)

  rm(rle_plot)
  gcFn()

  rle_elapsed <- as.numeric(difftime(timeFn(), rle_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: RLE plot completed in %.1f seconds ***", rle_elapsed))

  messageFn("*** PRE-NORM QC: Generating Density plot (Step 3/4) ***")
  density_start_time <- timeFn()

  density_plot <- buildDensityFn(current_s4, aesthetics)
  density_path <- saveArtifactFn(qcDir, "pre_norm_density.png", density_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_filtering", "density", density_path)

  rm(density_plot)
  gcFn()

  density_elapsed <- as.numeric(difftime(timeFn(), density_start_time, units = "secs"))
  messageFn(sprintf("*** PRE-NORM QC: Density plot completed in %.1f seconds ***", density_elapsed))

  messageFn("*** PRE-NORM QC: Generating Pearson correlation plot (Step 4/4) ***")
  num_samples <- length(setdiff(colnames(current_s4@protein_quant_table), current_s4@protein_id_column))
  estimated_pairs <- choose(num_samples, 2)
  messageFn(sprintf("*** PRE-NORM QC: Sample count = %d, Expected pairs ~= %d ***", num_samples, estimated_pairs))

  pearson_start_time <- timeFn()

  correlation_plot <- buildCorrelationFn(current_s4, aesthetics$color_var)
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
  initPathsFn = initializeProtNormQcPlotPaths
) {
  messageFn("=== GENERATING POST-NORMALIZATION QC PLOTS ===")

  qcPlotPaths <- initPathsFn(qcPlotPaths)

  messageFn("*** POST-NORM QC: Generating PCA plot ***")
  pca_plot <- plotPcaFn(
    normalizedS4,
    grouping_variable = aesthetics$color_var,
    label_column = "",
    shape_variable = aesthetics$shape_var,
    title = "",
    font_size = 8
  )

  pca_path <- saveArtifactFn(qcDir, "post_norm_pca.png", pca_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "pca", pca_path)

  rm(pca_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Generating RLE plot ***")
  rle_plot <- plotRleFn(
    normalizedS4,
    group = aesthetics$color_var,
    yaxis_limit = c(-6, 6)
  )

  rle_path <- saveArtifactFn(qcDir, "post_norm_rle.png", rle_plot, width = 10, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "rle", rle_path)

  rm(rle_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Generating Density plot ***")
  density_plot <- buildDensityFn(normalizedS4, aesthetics)
  density_path <- saveArtifactFn(qcDir, "post_norm_density.png", density_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "post_normalization", "density", density_path)

  rm(density_plot)
  gcFn()

  messageFn("*** POST-NORM QC: Generating Pearson correlation plot ***")
  correlation_plot <- buildCorrelationFn(normalizedS4, aesthetics$color_var)
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
  initPathsFn = initializeProtNormQcPlotPaths
) {
  messageFn("=== GENERATING RUV-CORRECTED QC PLOTS ===")

  qcPlotPaths <- initPathsFn(qcPlotPaths)

  messageFn("*** RUV QC: Generating PCA plot ***")
  pca_plot <- plotPcaFn(
    ruvCorrectedS4,
    grouping_variable = aesthetics$color_var,
    label_column = "",
    shape_variable = aesthetics$shape_var,
    title = "",
    font_size = 8
  )

  pca_path <- saveArtifactFn(qcDir, "ruv_corrected_pca.png", pca_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "pca", pca_path)

  rm(pca_plot)
  gcFn()

  messageFn("*** RUV QC: Generating RLE plot ***")
  rle_plot <- plotRleFn(
    ruvCorrectedS4,
    group = aesthetics$color_var,
    yaxis_limit = c(-6, 6)
  )

  rle_path <- saveArtifactFn(qcDir, "ruv_corrected_rle.png", rle_plot, width = 10, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "rle", rle_path)

  rm(rle_plot)
  gcFn()

  messageFn("*** RUV QC: Generating Density plot ***")
  density_plot <- buildDensityFn(ruvCorrectedS4, aesthetics)
  density_path <- saveArtifactFn(qcDir, "ruv_corrected_density.png", density_plot, width = 8, height = 6)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "density", density_path)

  rm(density_plot)
  gcFn()

  messageFn("*** RUV QC: Generating Pearson correlation plot ***")
  correlation_plot <- buildCorrelationFn(ruvCorrectedS4, aesthetics$color_var)
  corr_path <- saveArtifactFn(qcDir, "ruv_corrected_correlation.png", correlation_plot, width = 10, height = 8)
  qcPlotPaths <- recordPathFn(qcPlotPaths, "ruv_corrected", "correlation", corr_path)

  rm(correlation_plot)
  gcFn()

  messageFn("*** RUV QC: Running garbage collection ***")
  gcFn()
  messageFn("RUV-corrected QC plots generated successfully")

  qcPlotPaths
}

