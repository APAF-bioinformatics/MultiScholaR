# KIT-105 W3: proteomics' normalisation core. Back-references #10 to #22 (the functions the DIA reference path, the
# DIA sessions and the lipid/metab RUV declaration call), their module-file callees and the DIA_LIMPA pass-through,
# moved unchanged and with their names kept from R/mod_prot_norm_{workflow,ruv,qc_generation,qc_support,correlation,
# session}_helpers.R (design 2.2, 4.7). The tab's adapters, renderers and observers stay in its module files; KIT-111
# moves the run out of the GUI process.

# --- RUV-III parameters and results (#11, #12, #10; the manual result #11 calls) ---

buildProtNormManualRuvResult <- function(
  percentageAsNegCtrl,
  ruvK,
  controlGenesIndex,
  cancorPlot
) {
  list(
    best_percentage = percentageAsNegCtrl,
    best_k = ruvK,
    best_control_genes_index = controlGenesIndex,
    best_cancor_plot = cancorPlot,
    optimization_results = data.frame(
      percentage = percentageAsNegCtrl,
      separation_score = NA,
      best_k = ruvK,
      composite_score = NA,
      num_controls = sum(controlGenesIndex, na.rm = TRUE),
      valid_plot = TRUE
    ),
    separation_metric_used = "manual",
    k_penalty_weight = NA,
    adaptive_k_penalty_used = FALSE
  )
}

validateProtNormRuvK <- function(ruvK, paramLabel = "ruv_number_k") {
  if (is.null(ruvK) ||
      length(ruvK) != 1L ||
      is.na(ruvK) ||
      !is.numeric(ruvK) ||
      ruvK < 1) {
    stop(sprintf(
      "The %s = %s value is invalid.",
      paramLabel,
      paste(ruvK, collapse = ",")
    ), call. = FALSE)
  }

  as.integer(ruvK)
}

# Kernel-only RUV parameter determination (KIT-004): the same control selection, k validation and
# canonical-correlation plot as resolveProtNormRuvParameters, without progress, messages, holders or
# the config_list audit trail. The module keeps the adapter below; the operation runner calls this.
resolveProtRuvParameters <- function(
  normalizedS4,
  settings,
  ruvGroupingVariable,
  findBestNegCtrlPercentageFn = findBestNegCtrlPercentage,
  getNegCtrlProtAnovaFn = getNegCtrlProtAnova,
  ruvCancorFn = ruvCancor,
  buildManualRuvResultFn = buildProtNormManualRuvResult
) {
  if (identical(settings$ruv_mode, "automatic")) {
    percentage_min <- if (is.null(settings$auto_percentage_min) || is.na(settings$auto_percentage_min)) 1 else settings$auto_percentage_min
    percentage_max <- if (is.null(settings$auto_percentage_max) || is.na(settings$auto_percentage_max)) 20 else settings$auto_percentage_max
    if (percentage_min >= percentage_max) {
      stop("Minimum percentage must be less than maximum percentage")
    }
    result <- findBestNegCtrlPercentageFn(
      normalised_protein_matrix_obj = normalizedS4,
      percentage_range = seq(percentage_min, percentage_max, by = 1),
      separation_metric = settings$separation_metric,
      k_penalty_weight = settings$k_penalty_weight,
      adaptive_k_penalty = settings$adaptive_k_penalty,
      ruv_grouping_variable = ruvGroupingVariable,
      verbose = TRUE
    )
    ruv_k <- validateProtNormRuvK(result$best_k, "ruv_k")
    return(list(result = result, k = ruv_k, controls = result$best_control_genes_index,
      percentage = result$best_percentage))
  }
  percentage <- settings$ruv_percentage
  ruv_k <- if (is.null(settings$ruv_k) || is.na(settings$ruv_k)) 3 else settings$ruv_k
  ruv_k <- validateProtNormRuvK(ruv_k, "ruv_k")
  controls <- getNegCtrlProtAnovaFn(
    normalizedS4,
    ruv_grouping_variable = ruvGroupingVariable,
    percentage_as_neg_ctrl = percentage,
    ruv_qval_cutoff = 0.05,
    ruv_fdr_method = "BH"
  )
  # The adapter tolerates a failed canonical-correlation plot and leaves the result unset.
  result <- tryCatch({
    cancor_plot <- ruvCancorFn(normalizedS4, ctrl = controls, num_components_to_impute = 2,
      ruv_grouping_variable = ruvGroupingVariable)
    buildManualRuvResultFn(percentageAsNegCtrl = percentage, ruvK = ruv_k,
      controlGenesIndex = controls, cancorPlot = cancor_plot)
  }, error = function(e) {
    message(paste("Warning: Could not generate canonical correlation plot for manual mode:", e$message))
    NULL
  })
  list(result = result, k = ruv_k, controls = controls, percentage = result$best_percentage)
}

buildProtNormSkippedRuvResult <- function() {
  list(
    best_percentage = NA,
    best_k = NA,
    best_control_genes_index = NA,
    best_cancor_plot = NULL,
    ruv_skipped = TRUE,
    skip_reason = "User selected skip due to dataset constraints"
  )
}

# --- QC plot products (#13 to #16) ---

#' Remove captured plotting environments while preserving plot data and layers.
#' @param plotObject Existing ggplot or patchwork object.
#' @return Plot without ownership of producer-local scientific inputs.
#' @noRd
stripProtNormQcPlotOwners <- function(plotObject) {
  if (inherits(plotObject, c("gg", "ggplot"))) {
    plotObject$plot_env <- rlang::base_env()
  }
  if (!is.null(plotObject$patches) && is.list(plotObject$patches$plots)) {
    plotObject$patches$plots <- lapply(
      plotObject$patches$plots,
      stripProtNormQcPlotOwners
    )
  }
  plotObject
}

#' Produce one existing normalization QC plot as a bounded reusable product.
#' @param s4Object Existing ProteinQuantitativeData.
#' @param aesthetics Existing color/shape selection.
#' @param plotType One of pca, rle, density or correlation.
#' @param plotPcaFn,plotRleFn,buildDensityFn,buildCorrelationFn Existing producers.
#' @return Product containing ordered population metadata, plot data and baseline plot.
#' @noRd
buildProtNormQcPlotProduct <- function(
  s4Object,
  aesthetics,
  plotType,
  plotPcaFn = plotPca,
  plotRleFn = plotRle,
  buildDensityFn = buildProtNormDensityPlot,
  buildCorrelationFn = buildProtNormCorrelationPlot
) {
  if (!methods::is(s4Object, "ProteinQuantitativeData") ||
      !is.list(aesthetics) ||
      !identical(sort(names(aesthetics)), c("color_var", "shape_var")) ||
      !is.character(plotType) || length(plotType) != 1L || is.na(plotType) ||
      !plotType %in% c("pca", "rle", "density", "correlation")) {
    cli::cli_abort("invalid protein normalization QC product request",
      class = "saved_entity_invalid")
  }
  previous_plot <- ggplot2::get_last_plot()
  on.exit(ggplot2::set_last_plot(previous_plot), add = TRUE)
  custom_producer <- switch(plotType,
    pca = !identical(plotPcaFn, plotPca),
    rle = !identical(plotRleFn, plotRle),
    density = !identical(buildDensityFn, buildProtNormDensityPlot),
    correlation = !identical(buildCorrelationFn, buildProtNormCorrelationPlot))
  quant <- methods::slot(s4Object, "protein_quant_table")
  key <- methods::slot(s4Object, "protein_id_column")
  design <- methods::slot(s4Object, "design_matrix")
  sample_id <- methods::slot(s4Object, "sample_id")
  samples <- as.character(design[[sample_id]])
  annotations <- design[match(samples, design[[sample_id]]), , drop = FALSE]
  plot <- switch(plotType,
    pca = plotPcaFn(s4Object, grouping_variable = aesthetics$color_var,
      label_column = "", shape_variable = aesthetics$shape_var,
      title = "", font_size = 8),
    rle = plotRleFn(s4Object, group = aesthetics$color_var,
      yaxis_limit = c(-6, 6)),
    density = buildDensityFn(s4Object, aesthetics),
    correlation = buildCorrelationFn(s4Object, aesthetics$color_var)
  )
  data <- if (!is.null(plot$data) && is.data.frame(plot$data)) {
    as.data.frame(plot$data)
  } else {
    data.frame()
  }
  population_role <- if (identical(plotType, "density")) {
    "density_independent_pca_call"
  } else {
    paste0(plotType, "_existing_producer")
  }
  coordinate_limits <- if (!is.null(plot$coordinates$limits)) {
    plot$coordinates$limits
  } else {
    list(x = NULL, y = NULL)
  }
  product <- list(
    schema = "prot_norm_qc_plot_product_v1",
    plot_type = plotType,
    population_role = population_role,
    feature_ids = as.character(quant[[key]]),
    sample_ids = samples,
    annotations = annotations,
    scientific_parameters = list(
      color_variable = aesthetics$color_var,
      shape_variable = aesthetics$shape_var,
      pca_cv_percentile = 0.90,
      rle_yaxis_limit = c(-6, 6),
      correlation_exclude_pool_samples = TRUE
    ),
    render_metadata = list(
      x_label = if (is.null(plot$labels$x)) "" else as.character(plot$labels$x),
      y_label = if (is.null(plot$labels$y)) "" else as.character(plot$labels$y),
      x_limits = coordinate_limits$x,
      y_limits = coordinate_limits$y
    ),
    data = data,
    plot = NULL
  )
  # The reconstruction frame receives calculated display data, never this S4 or
  # the native plot's captured matrices/fits. Preserve explicitly injected seams.
  if (custom_producer) {
    product$plot <- stripProtNormQcPlotOwners(plot)
  } else {
    plot <- NULL
    product$plot <- stripProtNormQcPlotOwners(rebuildWorkflowDiaQcPlot(product))
  }
  product
}

buildProtNormDensityPlot <- function(
  s4Object,
  aesthetics,
  plotPcaFn = plotPca,
  plotPcaBoxFn = plotPcaBox
) {
  pcaPlot <- plotPcaFn(
    s4Object,
    grouping_variable = aesthetics$color_var,
    label_column = "",
    shape_variable = aesthetics$shape_var,
    title = "",
    font_size = 8
  )

  plotPcaBoxFn(
    pcaPlot,
    grouping_variable = aesthetics$color_var,
    show_legend = TRUE
  )
}

buildProtNormCorrelationPlot <- function(
  s4Object,
  colorVar,
  isRunningFn = shiny::isRunning,
  withProgressFn = shiny::withProgress,
  incProgressFn = shiny::incProgress,
  plotPearsonFn = plotPearson
) {
  if (isRunningFn()) {
    return(withProgressFn(
      message = "Generating Pearson correlation plot...",
      detail = "Calculating pairwise correlations...",
      value = 0.5,
      {
        result <- plotPearsonFn(
          s4Object,
          correlation_group = colorVar,
          exclude_pool_samples = TRUE
        )
        incProgressFn(0.5)
        result
      }
    ))
  }

  plotPearsonFn(
    s4Object,
    correlation_group = colorVar,
    exclude_pool_samples = TRUE
  )
}

# KIT-122: moved unchanged, with its name kept, from R/mod_prot_norm_qc_support_helpers.R: the worker's QC figures
# action (computeWorkflowDiaQcFiguresAction, utils_workflow_dia_products.R) renders through it, so that the image it
# saves is the one the tab's generators would have drawn.
saveProtNormQcPlotArtifact <- function(
  qcDir,
  filename,
  plotObject,
  width,
  height,
  dpi = 150,
  savePlotFn = ggplot2::ggsave
) {
  if (is.null(qcDir)) {
    return(NULL)
  }

  artifactPath <- file.path(qcDir, filename)
  manager <- getDIAActiveManager()
  product_reference <- attr(plotObject, "dia_product_reference", exact = TRUE)
  if (isSavedWorkflowReference(product_reference)) {
    dimensions <- as.integer(c(width, height) * dpi)
    aesthetics <- attr(plotObject, "dia_render_aesthetics", exact = TRUE)
    aesthetics <- aesthetics[!vapply(aesthetics, is.null, logical(1))]
    identity <- newSavedWorkflowRenderIdentity(product_reference,
      c(aesthetics, list(dpi = dpi)),
      "png", dimensions, manager$.context$software)
    hit <- resolveSavedWorkflowAsset(manager$.context$store, product_reference, identity)
    if (identical(hit$status, "hit")) {
      writeBin(hit$value$bytes, artifactPath)
      return(artifactPath)
    }
    savePlotFn(artifactPath, plotObject, width = width, height = height, dpi = dpi)
    # KIT-131: the PNG device on Windows writes nothing, and says nothing, from 260 characters.
    checkWorkflowFigureFiles(artifactPath, "a figure was drawn but its file was not written")
    asset <- newSavedWorkflowRenderedAsset(readBin(artifactPath, "raw", n = file.size(artifactPath)),
      "image/png", dimensions)
    reference <- saveWorkflowDiaQcAsset(manager$.context$store, asset,
      basename(tempfile("dia-qc-image-")), identity)
    rootDIAProduct(manager, reference, list(product_reference))
    return(artifactPath)
  }
  savePlotFn(artifactPath, plotObject, width = width, height = height, dpi = dpi)
  artifactPath
}

# --- KIT-259 (B2): a normalisation QC figure that cannot be drawn ---
# Before, a figure that failed to draw after the normalisation was committed left the step unfinished and every later
# step refused until Reset; the project kept the earlier run's file under that name. Now the figure is replaced, under
# its own name and size, by a placeholder image naming it and the reason, and the step completes with one warning.

#' A normalisation QC figure's name for the user, from its file name.
#' @param filename The figure's file name (`<stage>_<kind>.png`).
#' @param adjusted Whether the RUV figures show the RUV-III-adjusted matrix drawn for quality control only (the DIA
#'   with limpa RUV covariates step; KIT-259 third pass, VAL-307, W11), not a corrected one.
#' @param limpa Whether the figure belongs to the DIA with limpa route, whose first column of Figure 1 is DPC-Quant's
#'   estimates, not the table before the between-sample normalisation (KIT-266, VAL-316 (1)).
#' @return The name, or the file name when it is not a stage figure.
#' @noRd
protNormFigureTitle <- function(filename, adjusted = FALSE, limpa = FALSE) {
  stages <- c(pre_norm = "Pre-normalisation", post_norm = "Post-normalisation", ruv_corrected = "RUV-corrected",
    ruv_covariates = "RUV-III covariates")
  kinds <- c(pca = "PCA", rle = "RLE", density = "PCA density", correlation = "Pearson correlation",
    cancor = "canonical correlation")
  parts <- regmatches(filename, regexec("^(pre_norm|post_norm|ruv_corrected|ruv_covariates)_(pca|rle|density|correlation|cancor)[.]png$",
    filename))[[1L]]
  if (length(parts) != 3L) return(filename)
  if (adjusted && identical(parts[[2L]], "ruv_corrected")) stages[["ruv_corrected"]] <- "RUV-III-adjusted (for quality control)"
  if (limpa && identical(parts[[2L]], "pre_norm")) stages[["pre_norm"]] <- "DPC-Quant estimates"
  paste(stages[[parts[[2L]]]], kinds[[parts[[3L]]]])
}

#' Write the placeholder image of a QC figure that could not be drawn, at the figure's size, and its marker
#' `<path>.placeholder` (the reason, one line; KIT-259 second correction, VAL-287 (ii): the report reads it to name the
#' Figure 1 panels that were replaced). No partial file is left when either cannot be written.
#' @param path The figure's path.
#' @param figure The figure's file name.
#' @param reason Why it could not be drawn.
#' @param width,height,dpi The figure's size, as it is saved.
#' @param limpa Whether the figure names DPC-Quant's estimates as its first column (`protNormFigureTitle`).
#' @return `path`.
#' @noRd
writeProtNormFigurePlaceholder <- function(path, figure, reason, width, height, dpi = 150, limpa = FALSE) {
  label <- paste(c(paste0(protNormFigureTitle(figure, limpa = limpa), ": could not be drawn."), strwrap(reason, 70L)), collapse = "\n")
  plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = 5) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  tryCatch({
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
    writeLines(gsub("[\r\n\t]+", " ", reason), paste0(path, ".placeholder"), useBytes = TRUE)
    checkWorkflowFigureFiles(c(path, paste0(path, ".placeholder")),
      "a placeholder for a figure that could not be drawn was not written")
  }, error = function(e) {
    unlink(c(path, paste0(path, ".placeholder")))
    stop(e)
  })
  path
}

#' Draw one normalisation QC figure; when it cannot be drawn, write its placeholder under its name and signal
#' `prot_norm_figure_unavailable` (fields `figure`, `reason`). Diverges from conventions/R.md (catch specific
#' conditions): ggplot2's drawing errors carry no stable class, so a drawing error is caught generically, as the
#' lipidomics and metabolomics generators do, and logged with its message. A figure that could not be written
#' (`workflow_figure_unwritten`), a condition of the store (`saved_entity_*`), the head-ownership conditions
#' (`workflow_execution_stale`, `workflow_execution_busy`), R's memory-allocation error and an interrupt still stop
#' (KIT-259 second correction, VAL-287; AUDIT-KIT-259-impl I05). A figure drawn removes its earlier placeholder's
#' marker.
#' @param qcDir The figure's directory; NULL writes no placeholder and signals nothing: a figure that cannot be drawn
#'   is then only logged and NULL is returned.
#' @param filename The figure's file name.
#' @param width,height,dpi The figure's size.
#' @param draw `function()` that draws and writes the figure, returning its path.
#' @param placeholder Whether a placeholder replaces the figure (FALSE: the figure is skipped).
#' @param signal Whether `prot_norm_figure_unavailable` is signalled (the worker only logs).
#' @return The path written, the figure's or its placeholder's; NULL when nothing was written.
#' @noRd
drawProtNormQcFigure <- function(qcDir, filename, width, height, draw, dpi = 150, placeholder = TRUE, signal = TRUE) {
  drawn <- tryCatch(draw(), error = function(e) {
    if (inherits(e, c("workflow_figure_unwritten", "workflow_execution_stale", "workflow_execution_busy")) ||
        any(startsWith(class(e), "saved_entity")) ||
        grepl("cannot allocate (vector|memory)", conditionMessage(e))) stop(e)
    reason <- cli::ansi_strip(conditionMessage(e))
    logger::log_warn(paste("The QC figure", filename, "could not be drawn:", reason))
    if (!placeholder || is.null(qcDir)) return(NULL)
    path <- writeProtNormFigurePlaceholder(file.path(qcDir, filename), filename, reason, width, height, dpi)
    if (signal) rlang::warn(paste0("The QC figure '", filename, "' could not be drawn: ", reason),
      class = "prot_norm_figure_unavailable", figure = filename, reason = reason)
    structure(path, placeholder = TRUE)
  })
  if (isTRUE(attr(drawn, "placeholder"))) return(as.vector(drawn))
  if (!is.null(qcDir)) unlink(file.path(qcDir, paste0(filename, ".placeholder")))
  drawn
}

#' KIT-259 third pass (VAL-307; closure check C01): the panels of the combined QC figure that are placeholders, written
#' beside the figure when the publication builds it (`<omic>_composite_QC_figure.placeholders`, one line per panel: its
#' letter and its figure, tab-separated), so that the report names exactly the panels its Figure 1 shows as notes. The
#' panels are read from the markers beside the figure's inputs in `qcDir` (the publication's staging, where the
#' pre-normalisation markers are copied with their images). No placeholder, or no combined figure: no file.
#' @param ruvMode The run's RUV mode (the figure's layout).
#' @param qcDir The folder holding the combined figure and its inputs.
#' @param omicType The omic (the figure's file name).
#' @param resolveInputsFn The layout's resolver, `function(ruvMode, qcDir)`: the figure's own (KIT-266: the DIA with
#'   limpa route's, whose canonical-correlation image has another name).
#' @return The file written, or NULL, invisibly.
#' @noRd
writeProtNormCompositePlaceholders <- function(ruvMode, qcDir, omicType, resolveInputsFn = resolveProtNormCompositeFigureInputs) {
  if (!file.exists(file.path(qcDir, paste0(omicType, "_composite_QC_figure.png")))) return(invisible(NULL))
  inputs <- resolveInputsFn(ruvMode, qcDir)
  panels <- sub("[)]$", "", unlist(inputs$rowLabels, use.names = FALSE))
  figures <- basename(inputs$plotFiles)
  marked <- !is.na(inputs$plotFiles) & file.exists(paste0(inputs$plotFiles, ".placeholder"))
  if (!any(marked)) return(invisible(NULL))
  path <- file.path(qcDir, paste0(omicType, "_composite_QC_figure.placeholders"))
  writeLines(paste(panels[marked], sub("[.]png$", "", figures[marked]), sep = "\t"), path, useBytes = TRUE)
  invisible(path)
}

#' The warning notice naming the figures that could not be drawn.
#' @param figures A list of `list(figure, reason)`.
#' @param step The step that finished (`"normalisation"`, `"RUV-III covariates step"`), or NULL when the figures were
#'   drawn before it (a pre-normalisation figure redrawn on a managed session).
#' @return The notice's text.
#' @noRd
buildProtNormFigureUnavailableNotice <- function(figures, step = "normalisation") {
  adjusted <- identical(step, "RUV-III covariates step")
  named <- vapply(figures, function(item) paste0("'", protNormFigureTitle(item$figure, adjusted), "'"), "")
  reasons <- vapply(figures, function(item) item$reason, "")
  lead <- if (is.null(step)) "The " else paste0("The ", step, " finished, but the ")
  if (length(figures) == 1L) {
    return(paste0(lead, "figure ", named, " could not be drawn: ", reasons,
      ". The normalised data are not affected. In this tab, the composite QC figure and the report the figure is ",
      "replaced by a note giving the reason."))
  }
  paste0(if (is.null(step)) "These " else paste0("The ", step, " finished, but these "), "figures could not be drawn: ",
    paste0(named, " (", reasons, ")", collapse = "; "),
    ". The normalised data are not affected. In this tab, the composite QC figure and the report each figure is ",
    "replaced by a note giving the reason.")
}

#' Evaluate a normalisation run's completion; afterwards, one warning notice names the figures it could not draw.
#' @param code The run and its completion.
#' @param showNotificationFn The notice.
#' @param step The step that finished, for the notice (`buildProtNormFigureUnavailableNotice`).
#' @return The value of `code`.
#' @noRd
withProtNormFigureNotice <- function(code, showNotificationFn = shiny::showNotification, step = "normalisation") {
  figures <- list()
  value <- withCallingHandlers(code, prot_norm_figure_unavailable = function(w) {
    figures[[length(figures) + 1L]] <<- list(figure = w$figure, reason = w$reason)
    invokeRestart("muffleWarning")
  })
  if (length(figures)) {
    showNotificationFn(buildProtNormFigureUnavailableNotice(figures, step), type = "warning", duration = NULL)
  }
  value
}

# --- The composite QC figure and the workflow's completion (#21, #22; the composite's callees) ---

buildProtNormLabelPlot <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
    ggplot2::xlim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      panel.background = ggplot2::element_blank()
    )
}

buildProtNormTitlePlot <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
    ggplot2::xlim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 10, 5),
      panel.background = ggplot2::element_blank()
    )
}

loadProtNormImageAsPlot <- function(filePath, messageFn = message) {
  if (is.na(filePath) || !file.exists(filePath)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  tryCatch({
    img <- png::readPNG(filePath)
    g <- grid::rasterGrob(img, interpolate = TRUE)
    stripProtNormQcPlotOwners(ggplot2::ggplot() +
      ggplot2::annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::theme_void())
  }, error = function(e) {
    messageFn(sprintf("   [generateCompositeFromFiles] Could not load image: %s", filePath))
    ggplot2::ggplot() + ggplot2::theme_void()
  })
}

generateProtNormCompositeFromFiles <- function(
  plotFiles,
  ncol = 3,
  rowLabels = NULL,
  columnLabels = NULL,
  messageFn = message,
  warningFn = warning
) {
  messageFn(sprintf("   [generateCompositeFromFiles] Generating composite from %d files...", length(plotFiles)))

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warningFn("patchwork package required for composite generation")
    return(NULL)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warningFn("ggplot2 package required for composite generation")
    return(NULL)
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    warningFn("png package required for composite generation")
    return(NULL)
  }

  tryCatch({
    n_files <- length(plotFiles)
    n_plot_types <- n_files / ncol

    if (is.null(rowLabels)) {
      all_labels <- letters[seq_len(n_files)] # KIT-064 (roar BUG-16): as the lipid copy already does
      rowLabels <- split(paste0(all_labels, ")"), rep(seq_len(n_plot_types), each = ncol))
      names(rowLabels) <- paste0("row", seq_len(n_plot_types))
    }

    plot_sections <- list()
    height_values <- c()

    if (!is.null(columnLabels)) {
      if (length(columnLabels) != ncol) {
        warningFn("Number of column labels does not match ncol")
      } else {
        title_plots <- lapply(columnLabels, buildProtNormTitlePlot)
        plot_sections <- append(plot_sections, list(
          patchwork::wrap_plots(title_plots, ncol = ncol)
        ))
        height_values <- c(height_values, 0.2)
        messageFn("   [generateCompositeFromFiles] Added column titles")
      }
    }

    row_names <- names(rowLabels)

    for (i in seq_along(row_names)) {
      row_name <- row_names[i]
      labels <- rowLabels[[row_name]]

      start_idx <- (i - 1) * ncol + 1
      end_idx <- min(i * ncol, n_files)
      row_files <- plotFiles[start_idx:end_idx]

      has_files <- any(!is.na(row_files) & sapply(row_files, function(f) !is.na(f) && file.exists(f)))

      if (has_files || row_name == "cancor") {
        label_plots <- lapply(labels, buildProtNormLabelPlot)
        image_plots <- lapply(row_files, loadProtNormImageAsPlot, messageFn = messageFn)

        plot_sections <- append(plot_sections, list(
          patchwork::wrap_plots(label_plots, ncol = ncol),
          patchwork::wrap_plots(image_plots, ncol = ncol)
        ))
        # The single RUV control plot takes a shorter row than the three-panel rows.
        height_values <- c(height_values, 0.1, if (row_name == "cancor") 0.75 else 1)

        messageFn(sprintf("   [generateCompositeFromFiles] Added row: %s", row_name))
      } else {
        messageFn(sprintf("   [generateCompositeFromFiles] Skipping empty row: %s", row_name))
      }
    }

    if (length(plot_sections) == 0) {
      warningFn("No valid plot sections to combine")
      return(NULL)
    }

    messageFn("   [generateCompositeFromFiles] Combining plot sections...")
    combined_plot <- patchwork::wrap_plots(plot_sections, ncol = 1) +
      patchwork::plot_layout(heights = height_values)

    plot_width <- 4 + (ncol * 3)
    plot_height <- 4 + (length(height_values) * 2)

    rm(plot_sections)
    gc()

    list(plot = combined_plot, width = plot_width, height = plot_height)
  }, error = function(e) {
    messageFn(paste("   [generateCompositeFromFiles] Error:", e$message))
    NULL
  })
}

resolveProtNormCompositeFigureInputs <- function(ruvMode, qcDir) {
  if (ruvMode == "skip") {
    ncol_composite <- 2
    file_names <- c(
      "pre_norm_pca.png", "post_norm_pca.png",
      "pre_norm_density.png", "post_norm_density.png",
      "pre_norm_rle.png", "post_norm_rle.png",
      "pre_norm_correlation.png", "post_norm_correlation.png"
    )
    row_labels <- list(
      pca = c("a)", "b)"),
      density = c("c)", "d)"),
      rle = c("e)", "f)"),
      correlation = c("g)", "h)")
    )
    column_labels <- c("Pre-Normalisation", "Post-Normalisation")
  } else {
    ncol_composite <- 3
    file_names <- c(
      "pre_norm_pca.png", "post_norm_pca.png", "ruv_corrected_pca.png",
      "pre_norm_density.png", "post_norm_density.png", "ruv_corrected_density.png",
      "pre_norm_rle.png", "post_norm_rle.png", "ruv_corrected_rle.png",
      "pre_norm_correlation.png", "post_norm_correlation.png", "ruv_corrected_correlation.png",
      # The RUV control plot stands alone on its row: centred (the middle column), not under the RUV column.
      NA, "ruv_corrected_cancor.png", NA
    )
    row_labels <- list(
      pca = c("a)", "b)", "c)"),
      density = c("d)", "e)", "f)"),
      rle = c("g)", "h)", "i)"),
      correlation = c("j)", "k)", "l)"),
      cancor = c("", "m)", "")
    )
    column_labels <- c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
  }

  plot_files <- if (!is.null(qcDir)) {
    vapply(
      file_names,
      function(fn) {
        if (is.na(fn)) NA_character_ else file.path(qcDir, fn)
      },
      character(1)
    )
  } else {
    NULL
  }

  list(
    ncol = ncol_composite,
    plotFiles = plot_files,
    rowLabels = row_labels,
    columnLabels = column_labels
  )
}

generateProtNormCompositeQcFigure <- function(
  ruvMode,
  qcDir,
  omicType,
  resolveInputsFn = resolveProtNormCompositeFigureInputs,
  generateCompositeFn = generateProtNormCompositeFromFiles,
  savePlotFn = savePlot,
  messageFn = message
) {
  messageFn("*** STEP 6B: Generating composite QC figure from saved images ***")

  if (is.null(qcDir)) {
    return(invisible(NULL))
  }

  # ggplot2 retains the last constructed plot independently of our GUI holders.
  # Keep the caller's plot while releasing this export's raster composite on exit.
  previous_plot <- ggplot2::get_last_plot()
  on.exit(ggplot2::set_last_plot(previous_plot), add = TRUE)
  composite_inputs <- resolveInputsFn(ruvMode, qcDir)

  tryCatch({
    composite_res <- generateCompositeFn(
      plotFiles = composite_inputs$plotFiles,
      ncol = composite_inputs$ncol,
      rowLabels = composite_inputs$rowLabels,
      columnLabels = composite_inputs$columnLabels
    )

    if (!is.null(composite_res)) {
      savePlotFn(
        composite_res$plot,
        base_path = qcDir,
        plot_name = paste0(omicType, "_composite_QC_figure"),
        width = composite_res$width,
        height = composite_res$height,
        dpi = 150,
        limitsize = FALSE
      )
      messageFn("*** STEP 6B: Composite QC figure saved to protein_qc_dir ***")
    }
  }, error = function(e) {
    messageFn(paste("Warning: Could not generate composite QC figure:", e$message))
  })

  invisible(NULL)
}

finalizeProtNormWorkflowState <- function(
  normData,
  gcFn = gc,
  messageFn = message
) {
  normData$QC_composite_figure <- NULL
  messageFn("*** STEP 6: Clearing redundant plot objects from memory ***")
  gcFn()

  normData$ruv_complete <- TRUE
  messageFn("*** STEP 6: RUV-corrected workflow completed ***")
  messageFn("*** STEP 7: Enabling correlation filtering step ***")
  messageFn("*** STEP 7: Normalization and RUV workflow completed - ready for correlation filtering ***")

  invisible(NULL)
}

# --- The DIA_LIMPA route's normalisation pass-through (beside runDIAReferenceNormalization's callees) ---

#' The DIA_LIMPA route's normalisation action (KIT-041): the between-sample normalisation and the peptide-count
#' filter ran in the protein QC stage, so the tab's normalisation button designates the current state as the
#' correlation step's input, draws the post-normalisation QC from it and completes; the RUV covariates run after
#' the correlation filter.
#' @noRd
runDIAReferenceLimpaNormalisationPassThrough <- function(settings, workflow, norm_data, post_qc) {
  manager <- workflow$state_manager
  withDIAStateScope(manager, function() {
    assertDIAExportReady(manager, allow_action = "normalization")
    source <- getDIAStateReference(manager)
    # GUI audit P-35, KIT-244 (VAL-220 (1)): pressed again on the tab's own RUV state or correlation result, the button
    # keeps the input it designated before; the RUV state as the correlation's input would let the covariates that
    # follow compound, and the correlation result as the input would let the filter compound (a re-applied filter
    # acting at the higher of the two thresholds, a Skip keeping the first filter's removals).
    if (manager$current_state %in% c("ruv_covariates", "ruv_corrected", "correlation_filtered") &&
      isSavedWorkflowReference(norm_data$normalization_input_reference)) source <- norm_data$normalization_input_reference
    norm_data$normalization_input_reference <- source
    norm_data$normalized_protein_obj <- source
    norm_data$ruv_normalized_obj <- source
    norm_data$ruv_optimization_result <- NULL
    norm_data$best_k <- NA
    norm_data$control_genes_index <- NA
    post_qc(source)
    norm_data$normalization_complete <- TRUE
    norm_data$filtering_summary_text <- paste0(
      "DIA with limpa: between-sample normalisation (", diaStateHistory(manager) |> (\(h) if ("protein_normalised" %in% h) "applied in the protein QC stage" else "not run")(),
      ") and the peptide count filter are protein QC steps; RUV-III factors are computed as covariates after the correlation filter.")
    TRUE
  })
}

#' The source the DIA_LIMPA pass-through will designate, by its own rule (GUI audit P-35), for the managed tab (KIT-122):
#' the worker saves that source's QC figures before the pass-through draws them. The pass-through stays the authority;
#' were the two to differ, its figures would be drawn in process as before.
#' @param manager,norm_data Session owner and the tab's state.
#' @return The source reference.
#' @noRd
diaLimpaPassThroughSource <- function(manager, norm_data) {
  if (manager$current_state %in% c("ruv_covariates", "ruv_corrected", "correlation_filtered") &&
    isSavedWorkflowReference(norm_data$normalization_input_reference)) return(norm_data$normalization_input_reference)
  getDIAStateReference(manager)
}

# --- KIT-266 (VAL-308, VAL-316): the DIA_LIMPA route's combined QC figure (Figure 1) ---
# The route's normalisation commits no state, so the classic publication's builder never ran on it and the report described
# a Figure 1 it did not contain. The figure is built from the stage images in the project's `protein_qc` folder by the
# classic builder (`generateProtNormCompositeQcFigure`), with the layout of `resolveDiaLimpaCompositeFigureInputs`: two
# columns (DPC-Quant's estimates; the table after the peptide count filter) after the pass-through, three after the RUV
# covariates or correction step. Display files only: no state, number or table changes.

#' The between-sample method the route's protein QC stage recorded, for Figure 1's centre caption: the `norm_method` of the
#' `protein_normalised` state's configuration (the template reads `limpa_norm_method` of the analysis record; both take
#' an absent method as `"none"`).
#' @param manager The session owner.
#' @return The method name; `"none"` when the stage was not run or recorded none.
#' @noRd
diaLimpaNormMethod <- function(manager) {
  if (!inherits(manager, "DIAStateAdapter") || !"protein_normalised" %in% diaStateHistory(manager)) return("none")
  method <- resolveDIAStateValue(manager$.context, getDIAStateReference(manager, "protein_normalised"))$config$norm_method
  if (is.character(method) && length(method) == 1L && !is.na(method) && nzchar(method)) method else "none"
}

#' The saved table Figure 1's left column shows on this route: the DPC-Quant state's reference.
#' @param manager The session owner.
#' @return The reference, or NULL when the history has no `protein_dpc_quant` state.
#' @noRd
diaLimpaPreNormalisationSource <- function(manager) {
  if (!"protein_dpc_quant" %in% diaStateHistory(manager)) return(NULL)
  getDIAStateReference(manager, "protein_dpc_quant")
}

#' Make the tab's pre-normalisation figures' source the DPC-Quant state (as `diaNormalizationRequest` does for the
#' classic route): when the persisted source is absent or another table, set it and publish the head.
#' @param manager The session owner.
#' @return TRUE when the source was set, so the figures are to be drawn; FALSE when they are current or there is no
#'   DPC-Quant state.
#' @noRd
setDiaLimpaPreNormalisationSource <- function(manager) {
  reference <- diaLimpaPreNormalisationSource(manager)
  if (is.null(reference)) return(FALSE)
  snapshot <- manager$.context$snapshot
  drawn <- snapshot$qc_sources$pre_normalization
  if (isSavedWorkflowReference(drawn) && identical(drawn$entity_id, reference$entity_id)) return(FALSE)
  snapshot$qc_sources$pre_normalization <- reference
  publishDIAStateHead(manager$.context, snapshot)
  TRUE
}

#' Draw the tab's pre-normalisation figures from the DPC-Quant state, in process, on a press that designates the head:
#' the source is set, the figures are drawn with the tab's own generator and the tab does not draw them again on entry.
#' @param manager The session owner.
#' @param normData The tab's state.
#' @param preQc `function()`, the tab's pre-normalisation figure generator.
#' @return TRUE when the figures were drawn, invisibly.
#' @noRd
drawDiaLimpaPreNormalisationFigures <- function(manager, normData, preQc) {
  if (!inherits(manager, "DIAStateAdapter")) return(invisible(FALSE))
  withDIAStateScope(manager, function() {
    if (!setDiaLimpaPreNormalisationSource(manager)) return(invisible(FALSE))
    preQc()
    normData$pre_norm_qc_generated <- TRUE
    invisible(TRUE)
  })
}

#' Whether the persisted source of the left column's images is the DPC-Quant state.
#' @param manager The session owner.
#' @return TRUE or FALSE.
#' @noRd
diaLimpaLeftColumnCurrent <- function(manager) {
  expected <- diaLimpaPreNormalisationSource(manager)
  drawn <- manager$.context$snapshot$qc_sources$pre_normalization
  !is.null(expected) && isSavedWorkflowReference(drawn) && identical(drawn$entity_id, expected$entity_id)
}

#' The source check of every limpa build (VAL-316 (1)): when the left column's images are not drawn from the DPC-Quant
#' state (a session saved by an earlier build, or no DPC-Quant state), the four left cells are a labelled note instead.
#' The notes are written straight into the project, not through `drawProtNormQcFigure`, so no step notice is signalled
#' (the report's sentence names the panels).
#' @param manager The session owner.
#' @param qc_dir The project's `protein_qc` folder.
#' @return TRUE when the notes were written, invisibly.
#' @noRd
writeDiaLimpaPreNormalisationNotes <- function(manager, qc_dir) {
  if (is.null(qc_dir) || diaLimpaLeftColumnCurrent(manager)) return(invisible(FALSE))
  reason <- paste("the figures saved in this session for this column come from a later table, or no saved table from",
    "before the between-sample normalisation step is available.")
  for (kind in names(diaQcStageFigureSizes)) {
    figure <- paste0("pre_norm_", kind, ".png")
    size <- diaQcStageFigureSizes[[kind]]
    writeProtNormFigurePlaceholder(file.path(qc_dir, figure), figure, reason, size[[1L]], size[[2L]], limpa = TRUE)
  }
  invisible(TRUE)
}

#' The between-sample method's name in Figure 1's centre caption.
#' @noRd
diaLimpaNormMethodName <- function(method) {
  switch(method, cyclicloess = "cyclic loess", quantile = "quantile", scale = "scale", method)
}

#' The limpa route's column captions of Figure 1 (drafts; the owner sets the final wording, VAL-41), at most three lines
#' of thirty characters, broken explicitly.
#' @param ruvMode `skip` (two columns), `covariates` or `correction`.
#' @param normMethod The recorded between-sample method.
#' @return The captions, one per column.
#' @noRd
diaLimpaCompositeColumnLabels <- function(ruvMode, normMethod) {
  centre <- if (identical(normMethod, "none")) "After the peptide count filter\n(no between-sample\nnormalisation applied)" else
    paste0("After ", diaLimpaNormMethodName(normMethod), " normalisation\nand the peptide count filter")
  right <- if (identical(ruvMode, "skip")) NULL else if (identical(ruvMode, "correction")) "RUV-III-corrected\n(the data analysed)" else
    "RUV-III-adjusted\n(for quality control only;\nnot the data analysed)"
  c("DPC-Quant estimates\n(log2)", centre, right)
}

#' The layout of the limpa route's combined figure: the classic resolver's, with the limpa captions and the covariates
#' step's own canonical-correlation image (panel m, only when it exists).
#' @param ruvMode `skip`, `covariates` or `correction`.
#' @param qcDir The folder holding the stage images.
#' @param norm_method The recorded between-sample method.
#' @return `list(ncol, plotFiles, rowLabels, columnLabels)`, as `resolveProtNormCompositeFigureInputs`.
#' @noRd
resolveDiaLimpaCompositeFigureInputs <- function(ruvMode, qcDir, norm_method = "none") {
  inputs <- resolveProtNormCompositeFigureInputs(ruvMode, qcDir)
  inputs$columnLabels <- diaLimpaCompositeColumnLabels(ruvMode, norm_method)
  if (identical(ruvMode, "skip") || is.null(inputs$plotFiles)) return(inputs)
  cancor <- file.path(qcDir, "ruv_covariates_cancor.png")
  if (file.exists(cancor)) {
    inputs$plotFiles[[14L]] <- cancor
  } else {
    inputs$plotFiles <- inputs$plotFiles[1:12]
    inputs$rowLabels$cancor <- NULL
  }
  inputs
}

#' The layout of the limpa figure as a string, for the identities of its files in the store: the mode, the method and (for
#' the three-column layouts) the canonical-correlation image's content, so that an identity follows the layout and its
#' captions. It cannot collide with a classic RUV mode (the `limpa` prefix).
#' @noRd
diaLimpaCompositeLayoutString <- function(mode, norm_method, qc_dir) {
  cancor <- file.path(qc_dir, "ruv_covariates_cancor.png")
  paste("limpa", mode, norm_method, if (!identical(mode, "skip") && file.exists(cancor))
    digest::digest(file = cancor, algo = "sha256") else "none", sep = "|")
}

#' The list of replaced panels beside the limpa figure (KIT-259's helper through the limpa layout), written afresh: an
#' earlier build's list is removed first, since the helper writes only when a panel is replaced.
#' @noRd
writeDiaLimpaCompositePlaceholders <- function(mode, qc_dir, omic_type, resolveInputsFn) {
  unlink(file.path(qc_dir, paste0(omic_type, "_composite_QC_figure.placeholders")))
  writeProtNormCompositePlaceholders(mode, qc_dir, omic_type, resolveInputsFn = resolveInputsFn)
}

#' Build the limpa route's combined QC figure in the project's `protein_qc` folder, in the application (the in-process
#' route, and the managed route's fallback): the source check of the left column, the classic builder through the limpa
#' layout, then KIT-259's list of replaced panels. Nothing without a project folder or an owner.
#' @param paths The omic's project paths (`protein_qc_dir`).
#' @param omic_type The omic, which names the files.
#' @param ruv_mode `skip` (two columns), `covariates` or `correction`.
#' @param manager The session owner.
#' @return NULL, invisibly.
#' @noRd
buildDiaLimpaCombinedQcFigure <- function(paths, omic_type, ruv_mode, manager) {
  qc_dir <- paths$protein_qc_dir
  if (is.null(qc_dir) || !inherits(manager, "DIAStateAdapter")) return(invisible(NULL))
  withDIAStateScope(manager, function() {
    norm_method <- diaLimpaNormMethod(manager)
    writeDiaLimpaPreNormalisationNotes(manager, qc_dir)
    resolver <- function(ruvMode, qcDir) resolveDiaLimpaCompositeFigureInputs(ruvMode, qcDir, norm_method)
    generateProtNormCompositeQcFigure(ruv_mode, qc_dir, omic_type, resolveInputsFn = resolver)
    writeDiaLimpaCompositePlaceholders(ruv_mode, qc_dir, omic_type, resolver)
  })
  invisible(NULL)
}

#' The managed route's delivery of the limpa figure: the files the worker saved are written into the project (their
#' identity names the layout and the stage images), then KIT-259's list. When they are not found saved the figure is built
#' here, as the classic publication does, and the reason is logged: it re-creates the application's memory event the
#' worker branch removes.
#' @param manager The session owner.
#' @param state `list(mode, method, source)`: the layout, the method and the state the action named.
#' @param paths,omic_type As `buildDiaLimpaCombinedQcFigure`.
#' @return NULL, invisibly.
#' @noRd
deliverDiaLimpaCompositeFigure <- function(manager, state, paths, omic_type) {
  qc_dir <- paths$protein_qc_dir
  if (is.null(qc_dir)) return(invisible(NULL))
  layout <- diaLimpaCompositeLayoutString(state$mode, state$method, qc_dir)
  if (!writeDIACompositeFigureFiles(manager, state$source, layout, qc_dir, omic_type)) {
    logger::log_warn(paste("The combined QC figure of the DIA with limpa route was not found saved by the worker",
      "(identities not found saved); it is built in the application:", layout))
    return(buildDiaLimpaCombinedQcFigure(paths, omic_type, state$mode, manager))
  }
  resolver <- function(ruvMode, qcDir) resolveDiaLimpaCompositeFigureInputs(ruvMode, qcDir, state$method)
  writeDiaLimpaCompositePlaceholders(state$mode, qc_dir, omic_type, resolver)
  invisible(NULL)
}

#' The layout the head has now (the managed composite action reads it when it starts, so a press that waited behind another
#' reads what the action before it left): three columns after a RUV step that applied RUV (`ruv_covariates`: covariates,
#' `ruv_corrected`: correction), two otherwise.
#' @param manager The session owner (in its scope).
#' @param normData The tab's state (`best_k`).
#' @return `skip`, `covariates` or `correction`.
#' @noRd
diaLimpaCompositeLayout <- function(manager, normData) {
  applied <- isTRUE(shiny::isolate(normData$best_k) > 0)
  state <- manager$current_state
  if (applied && identical(state, "ruv_covariates")) "covariates" else if (applied && identical(state, "ruv_corrected")) "correction" else "skip"
}

#' Whether a correlation filter about to run on a DIA with limpa project leaves a RUV result behind it: the head, read
#' before the correlation runs, is a RUV state, and the filter commits `correlation_filtered` after it, so the three-column
#' figure is rebuilt with two (a first correlation after the pass-through, or after a Reset, rebuilds nothing).
#' @param manager The session owner.
#' @return TRUE or FALSE.
#' @noRd
protNormLimpaRebuildAtCorrelation <- function(manager) {
  inherits(manager, "DIAStateAdapter") && identical(manager$workflow_type, "DIA_LIMPA") &&
    manager$current_state %in% c("ruv_covariates", "ruv_corrected")
}

#' The markers of the RUV figures' placeholders that a RUV covariates or correction run replaces, removed by exact name
#' before the run, so a failed earlier run cannot make the list of replaced panels name a drawn panel. KIT-276 (note 1.5,
#' AUDIT-KIT-266-impl I03): also the canonical-correlation image itself, so that a declared-rule run with no sweep does not
#' show an earlier run's plot as panel m.
#' @param paths The omic's project paths (`protein_qc_dir`).
#' @return The removed paths, invisibly; a condition of class `workflow_figure_unwritten` when one cannot be removed.
#' @noRd
removeDiaLimpaRuvMarkers <- function(paths) {
  names <- paste0(c(paste0("ruv_corrected_", c("pca", "rle", "density", "correlation")), "ruv_covariates_cancor"), ".png.placeholder")
  removeWorkflowOwnFiles(file.path(paths$protein_qc_dir, c(names, "ruv_covariates_cancor.png")), omicNormGuardTexts()$norm_unremoved)
}

#' KIT-276 (S02): the RUV run's own project file `ruv_optimization_results.RDS` (a session export writes it from the RUV
#' result), removed by exact name when the table no longer carries the RUV result (a correlation applied again after a RUV
#' step, or "Skip RUV"), so that the study parameters do not print that run's sweep beside a table without RUV.
#' @param paths The omic's project paths (`source_dir`).
#' @return The removed paths, invisibly; a condition of class `workflow_figure_unwritten` when one cannot be removed.
#' @noRd
removeProtNormLeftoverRuvResults <- function(paths) {
  if (is.null(paths$source_dir)) return(invisible(character()))
  removeWorkflowOwnFiles(file.path(paths$source_dir, "ruv_optimization_results.RDS"), omicNormGuardTexts()$norm_unremoved)
}

# --- Correlation result files (#17, #18) ---

resolveProtNormCorrelationResultFilenames <- function(
  normRuvOptimizationResult = NULL,
  workflowRuvOptimizationResult = NULL,
  messagePrefix = "*** CORRELATION",
  messageFn = message
) {
  ruv_was_skipped <- isTRUE(normRuvOptimizationResult$ruv_skipped) ||
    isTRUE(workflowRuvOptimizationResult$ruv_skipped)

  if (ruv_was_skipped) {
    messageFn(sprintf("%s: RUV was skipped - using 'normalised' file names ***", messagePrefix))
    return(list(
      ruvWasSkipped = TRUE,
      tsvFilename = "normalised_results_cln_with_replicates.tsv",
      rdsFilename = "normalised_results_cln_with_replicates.RDS"
    ))
  }

  messageFn(sprintf("%s: RUV was applied - using 'ruv_normalised' file names ***", messagePrefix))

  list(
    ruvWasSkipped = FALSE,
    tsvFilename = "ruv_normalised_results_cln_with_replicates.tsv",
    rdsFilename = "ruv_normalised_results_cln_with_replicates.RDS"
  )
}

#' KIT-271 (VAL-317): the file names of the correlation filter's results table on the route that writes it. On DIA
#' with limpa the filter's table is the protein table before any RUV step in every RUV mode (the covariates are computed
#' after the filter and the RUV matrices are other files), so it takes the `normalised` names whatever RUV result the
#' tab holds, which on that route is absent or applied by design. The pinned resolver stays the one place that holds the
#' names, called once on every route; the other routes get its value unchanged.
#' @param manager The session owner (only its `workflow_type` is read), or NULL.
#' @param normRuvOptimizationResult,workflowRuvOptimizationResult The tab's and the workflow's RUV results.
#' @param messagePrefix,messageFn The console line's prefix and writer.
#' @return The resolver's list (`ruvWasSkipped`, `tsvFilename`, `rdsFilename`).
#' @noRd
resolveProtNormCorrelationTableNames <- function(
  manager,
  normRuvOptimizationResult = NULL,
  workflowRuvOptimizationResult = NULL,
  messagePrefix = "*** CORRELATION",
  messageFn = message
) {
  workflow_type <- tryCatch(manager$workflow_type, error = function(e) NULL)
  if (!identical(workflow_type, "DIA_LIMPA")) {
    return(resolveProtNormCorrelationResultFilenames(normRuvOptimizationResult, workflowRuvOptimizationResult,
      messagePrefix, messageFn))
  }
  table_names <- resolveProtNormCorrelationResultFilenames(list(ruv_skipped = TRUE), NULL, messagePrefix,
    messageFn = function(...) invisible(NULL))
  messageFn(sprintf(
    "%s: DIA with limpa - the correlation filter's table precedes any RUV step; using 'normalised' file names ***",
    messagePrefix))
  table_names
}

#' KIT-257 (VAL-250): before a proteomics normalisation or correlation run, remove its own earlier files from
#' `protein_qc_dir`, by exact name: both names of the correlation's results table (TSV and the S4 compatibility RDS),
#' and at a normalisation also the RUV-corrected figures and the cancor figure (written with RUV only). Copy prefers the
#' RUV table whenever it exists, so an earlier run's table left there was delivered as the current one.
#' @param paths The omic's project paths (`protein_qc_dir`).
#' @param stage `normalization` or `correlation`.
#' @return The removed paths, invisibly; a condition of class `workflow_figure_unwritten` when one cannot be removed.
#' @noRd
removeProtNormOwnFiles <- function(paths, stage = c("normalization", "correlation")) {
  stage <- match.arg(stage)
  tables <- paste0(rep(c("ruv_normalised", "normalised"), 2L), "_results_cln_with_replicates.", rep(c("tsv", "RDS"), each = 2L))
  # KIT-259 second correction (VAL-287 (ii)): and the markers of the post-normalisation and RUV-corrected figures'
  # placeholders, which the report reads.
  names <- if (identical(stage, "normalization")) c(paste0("ruv_corrected_", c("pca", "rle", "density", "correlation", "cancor"),
    ".png"), tables, paste0(rep(c("post_norm_", "ruv_corrected_"), each = 5L),
      c("pca", "rle", "density", "correlation", "cancor"), ".png.placeholder"),
    # KIT-259 third pass (VAL-307, C01): the combined figure's list of replaced panels.
    "proteomics_composite_QC_figure.placeholders") else tables
  texts <- omicNormGuardTexts()
  removeWorkflowOwnFiles(file.path(paths$protein_qc_dir, names),
    if (identical(stage, "normalization")) texts$norm_unremoved else texts$correlation_unremoved)
}

#' AUDIT-KIT-257-impl I05: a failure after this run removed its earlier files carries them (class
#' `prot_norm_own_files_removed`, field `removed_own_files`), so the tab's notice can say so; the message is unchanged.
#' @param error The condition; `removed` the paths the run removed.
#' @return The condition.
#' @noRd
protNormTagRemovedFiles <- function(error, removed) {
  if (!length(removed) || inherits(error, "prot_norm_own_files_removed")) return(error)
  error$removed_own_files <- removed
  class(error) <- c("prot_norm_own_files_removed", class(error))
  error
}

#' Run `code`; an error after the removal is raised again carrying the removed paths (`protNormTagRemovedFiles`).
#' @noRd
withProtNormRemovedFiles <- function(removed, code) {
  tryCatch(code, error = function(e) stop(protNormTagRemovedFiles(e, removed)))
}

#' The tab's error notice of a run that stopped after removing its earlier files: the reason's first line gains
#' `removed_then_stopped` (the first line, so the in-process and worker notices keep the same first line).
#' @noRd
protNormRemovedFilesNotice <- function(error) {
  if (!inherits(error, "prot_norm_own_files_removed")) return(error)
  text <- omicNormGuardTexts()$removed_then_stopped
  message <- error$message %||% conditionMessage(error)
  cut <- regexpr("\n", message, fixed = TRUE)
  error$message <- if (cut < 0) paste(message, text) else paste0(substr(message, 1L, cut - 1L), " ", text, substring(message, cut))
  error
}

saveProtNormCorrelationResults <- function(
  finalS4ForDe,
  experimentPaths,
  normRuvOptimizationResult = NULL,
  workflowRuvOptimizationResult = NULL,
  messagePrefix = "*** CORRELATION STEP 4",
  writeTsvFn = readr::write_tsv,
  saveRdsFn = saveRDS,
  messageFn = message,
  compatibilityExport = !isDIANativeSession(getDIAActiveManager())
) {
  file_names <- resolveProtNormCorrelationResultFilenames(
    normRuvOptimizationResult = normRuvOptimizationResult,
    workflowRuvOptimizationResult = workflowRuvOptimizationResult,
    messagePrefix = messagePrefix,
    messageFn = messageFn
  )

  if (!is.null(experimentPaths) && "protein_qc_dir" %in% names(experimentPaths)) {
    writeTsvFn(
      finalS4ForDe@protein_quant_table,
      file.path(experimentPaths$protein_qc_dir, file_names$tsvFilename)
    )

    if (isTRUE(compatibilityExport)) {
      saveRdsFn(
        finalS4ForDe,
        file.path(experimentPaths$protein_qc_dir, file_names$rdsFilename)
      )
    } else {
      file_names$rdsFilename <- NULL
    }

    messageFn(sprintf(
      "%s: Saved files: %s ***",
      messagePrefix,
      paste(c(file_names$tsvFilename, file_names$rdsFilename), collapse = ", ")
    ))
  }

  invisible(file_names)
}

# --- Session export metadata and summary (#19, #20) ---

buildProtNormExportSummaryContent <- function(
  sessionData,
  sessionFilename,
  timeFn = Sys.time,
  formatTimeFn = format
) {
  sprintf(
    paste(
      "Filtered Session Data Export Summary",
      "=====================================",
      "",
      "Export Timestamp: %s",
      "Session File: %s",
      "",
      "Data Summary:",
      "- Proteins: %d",
      "- Samples: %d",
      "- Contrasts: %d",
      "- Normalization: %s",
      "- RUV Mode: %s",
      "- RUV K: %s",
      "- Correlation Threshold: %s",
      "- Workflow Type: %s",
      "- Report Template: %s",
      "- limpa DPC-Quant: %s",
      "- limpa Quantification Method: %s",
      "",
      "Contrasts:",
      "%s",
      "",
      "This data is ready for differential expression analysis.",
      "Use 'Load Filtered Session' in the DE tab to import.",
      sep = "\n"
    ),
    formatTimeFn(timeFn(), "%Y-%m-%d %H:%M:%S"),
    sessionFilename,
    sessionData$final_protein_count,
    sessionData$final_sample_count,
    ifelse(is.null(sessionData$contrasts_tbl), 0, nrow(sessionData$contrasts_tbl)),
    sessionData$normalization_method,
    sessionData$ruv_mode,
    ifelse(is.null(sessionData$ruv_k), "NA", sessionData$ruv_k),
    ifelse(is.null(sessionData$correlation_threshold), "NA", sessionData$correlation_threshold),
    ifelse(is.null(sessionData$workflow_type), "NA", sessionData$workflow_type),
    ifelse(is.null(sessionData$report_template), "NA", sessionData$report_template),
    ifelse(isTRUE(sessionData$limpa_applied), "yes", "no"),
    ifelse(
      is.null(sessionData$limpa_dpc_quant_results$quantification_method),
      "NA",
      sessionData$limpa_dpc_quant_results$quantification_method
    ),
    if (!is.null(sessionData$contrasts_tbl)) paste(sessionData$contrasts_tbl$friendly_names, collapse = "\n") else "None"
  )
}

saveProtNormExportMetadataFiles <- function(
  sessionData,
  sourceDir,
  saveRdsFn = saveRDS,
  fileExistsFn = file.exists,
  messageFn = message,
  strict = FALSE
) {
  tryCatch({
    if (!is.null(sessionData$accession_cleanup_results)) {
      saveRdsFn(sessionData$accession_cleanup_results, file.path(sourceDir, "accession_cleanup_results.RDS"))
      messageFn("*** EXPORT: Saved accession_cleanup_results.RDS ***")
    }

    if (!is.null(sessionData$qc_params)) {
      saveRdsFn(sessionData$qc_params, file.path(sourceDir, "qc_params.RDS"))
      messageFn("*** EXPORT: Saved qc_params.RDS ***")
    }

    if (!is.null(sessionData$fasta_metadata)) {
      fasta_metadata_file <- file.path(sourceDir, "fasta_metadata.RDS")
      if (!fileExistsFn(fasta_metadata_file)) {
        saveRdsFn(sessionData$fasta_metadata, fasta_metadata_file)
        messageFn("*** EXPORT: Saved fasta_metadata.RDS ***")
      }
    }

    if (!is.null(sessionData$ruv_optimization_result)) {
      ruv_file <- file.path(sourceDir, "ruv_optimization_results.RDS")
      if (!fileExistsFn(ruv_file)) {
        saveRdsFn(sessionData$ruv_optimization_result, ruv_file)
        messageFn("*** EXPORT: Saved ruv_optimization_results.RDS ***")
      }
    }

    if (!is.null(sessionData$mixed_species_analysis)) {
      mixed_species_file <- file.path(sourceDir, "mixed_species_analysis.RDS")
      saveRdsFn(sessionData$mixed_species_analysis, mixed_species_file)
      messageFn(sprintf(
        "*** EXPORT: Saved mixed_species_analysis.RDS (enabled: %s) ***",
        isTRUE(sessionData$mixed_species_analysis$enabled)
      ))
    }

    if (!is.null(sessionData$limpa_parameters)) {
      saveRdsFn(sessionData$limpa_parameters, file.path(sourceDir, "limpa_parameters.RDS"))
      messageFn("*** EXPORT: Saved limpa_parameters.RDS ***")
    }

    if (!is.null(sessionData$limpa_dpc_quant_results)) {
      saveRdsFn(sessionData$limpa_dpc_quant_results, file.path(sourceDir, "limpa_dpc_quant_results.RDS"))
      messageFn("*** EXPORT: Saved limpa_dpc_quant_results.RDS ***")
    }
  }, error = function(e) {
    messageFn(sprintf("*** WARNING: Some metadata files could not be saved: %s ***", e$message))
    if (strict) stop(e)
  })

  invisible(NULL)
}
