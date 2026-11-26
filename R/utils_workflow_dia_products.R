#' Build the complete identity for one existing DIA normalization QC producer.
#' @param s4Object Existing ProteinQuantitativeData input.
#' @param stage Existing QC stage name.
#' @param plotType pca, rle, density or correlation.
#' @param aesthetics Existing color/shape selections.
#' @param software Installed software/build identity.
#' @return Scientific identity for the shared saved-entity store.
#' @noRd
newWorkflowDiaQcProductIdentity <- function(s4Object, stage, plotType,
                                            aesthetics, software) {
  if (!methods::is(s4Object, "ProteinQuantitativeData") ||
      !plotType %in% c("pca", "rle", "density", "correlation")) {
    cli::cli_abort("invalid DIA QC identity input", class = "saved_entity_invalid")
  }
  checkSavedEntityString(stage, "stage")
  checkSavedEntityString(software, "software")
  quant <- methods::slot(s4Object, "protein_quant_table")
  design <- methods::slot(s4Object, "design_matrix")
  key <- methods::slot(s4Object, "protein_id_column")
  sample_id <- methods::slot(s4Object, "sample_id")
  input <- digest::digest(list(key = key, feature_ids = quant[[key]],
    sample_columns = setdiff(names(quant), key), values = quant), algo = "sha256")
  design_id <- digest::digest(list(sample_id = sample_id, design = design),
    algo = "sha256")
  parameters <- digest::digest(list(plot_type = plotType,
    color_variable = aesthetics$color_var,
    shape_variable = aesthetics$shape_var,
    pca_cv_percentile = 0.90,
    rle_yaxis_limit = c(-6, 6),
    correlation_exclude_pool_samples = TRUE), algo = "sha256")
  population <- if (identical(plotType, "density"))
    "density_independent_pca_call" else paste0(plotType, "_existing_producer")
  selectors <- digest::digest(list(stage = stage, population = population),
    algo = "sha256")
  newSavedWorkflowScientificIdentity(input, design_id, parameters, selectors, software)
}

#' Save one bounded DIA QC product through the shared product contract.
#' @param store Shared store.
#' @param product Product returned by buildProtNormQcPlotProduct.
#' @param request_id Idempotent request identity.
#' @param identity Complete scientific identity.
#' @param stage Existing workflow stage.
#' @return Accepted scientific-product descriptor.
#' @noRd
saveWorkflowDiaQcProduct <- function(store, product, request_id, identity, stage) {
  if (!is.list(product) || !identical(product$schema,
      "prot_norm_qc_plot_product_v1")) {
    cli::cli_abort("invalid DIA QC product", class = "saved_entity_invalid")
  }
  persisted_product <- product
  persisted_product["plot"] <- list(NULL)
  envelope <- newSavedWorkflowScientificProduct(persisted_product,
    paste0("dia_normalization_", product$plot_type),
    provenance = list(route = "iq", stage = stage,
      population_role = product$population_role))
  saveWorkflowConsumerEntity(store, envelope, request_id, identity,
    "scientific_product", factors = list(route = "iq", stage = stage,
      plot_type = product$plot_type))
}

#' Rebuild an existing QC display from bounded calculated product data.
#' @param product Saved or newly produced QC product.
#' @return Baseline-equivalent plot without access to the source assay.
#' @noRd
rebuildWorkflowDiaQcPlot <- function(product) {
  data <- product$data
  parameters <- product$scientific_parameters
  metadata <- product$render_metadata
  type <- product$plot_type
  if (!is.data.frame(data) || !type %in% c("pca", "rle", "density", "correlation")) {
    cli::cli_abort("invalid saved DIA QC plot product", class = "saved_entity_invalid")
  }
  if (identical(type, "density")) {
    seed <- ggplot2::ggplot(data, ggplot2::aes(.data[["PC1"]], .data[["PC2"]]))
    return(plotPcaBox(seed, grouping_variable = parameters$color_variable,
      show_legend = TRUE))
  }
  if (identical(type, "correlation")) {
    # KIT-259 (A2): no within-group pair with a correlation gives a labelled figure, not an error.
    if (!any(is.finite(data$pearson_correlation))) {
      return(buildProtPearsonUnavailablePlot(data, parameters$color_variable))
    }
    breaks <- seq(min(round(data$pearson_correlation - 0.5, 2), na.rm = TRUE),
      1, 0.001)
    return(ggplot2::ggplot(data, ggplot2::aes(.data[["pearson_correlation"]])) +
      ggplot2::geom_histogram(breaks = breaks) +
      protPearsonCountScale() +
      ggplot2::xlab("Pearson Correlation") + ggplot2::ylab("Counts") +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()))
  }
  if (identical(type, "rle")) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["rle.x.factor"]])) +
      ggplot2::geom_boxplot(ggplot2::aes(lower = .data[["lower"]],
        middle = .data[["middle"]], upper = .data[["upper"]],
        max = .data[["max"]], min = .data[["min"]]), stat = "identity") +
      ggplot2::theme_bw() + ggplot2::geom_hline(yintercept = 0) +
      ggplot2::coord_cartesian(ylim = parameters$rle_yaxis_limit) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.5)))
    if ("rowinfo" %in% names(data)) {
      plot <- plot + ggplot2::aes(fill = .data[["rowinfo"]]) + ggplot2::labs(fill = "")
    }
    return(plot)
  }
  mapping <- if (is.null(parameters$shape_variable) ||
      !nzchar(parameters$shape_variable)) {
    ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
      color = .data[[parameters$color_variable]])
  } else {
    ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
      color = .data[[parameters$color_variable]],
      shape = .data[[parameters$shape_variable]])
  }
  plot <- ggplot2::ggplot(data, mapping) + ggplot2::geom_point(size = 3) +
    ggplot2::xlab(metadata$x_label) + ggplot2::ylab(metadata$y_label) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(labels = function(x) format(x, scientific = FALSE, digits = 3)) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE, digits = 3)) +
    ggplot2::scale_color_manual(values = getCategoricalColourPalette()) +
    ggplot2::coord_cartesian(xlim = metadata$x_limits, ylim = metadata$y_limits)
  if (!is.null(parameters$shape_variable) && nzchar(parameters$shape_variable)) {
    plot <- plot + ggplot2::scale_shape_manual(
      values = c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4, 5, 6, 9, 10, 11))
  }
  plot
}

#' Resolve one exact saved DIA QC product without invoking its producer.
#' @param store Shared store.
#' @param identity Scientific identity.
#' @param stage,plotType Exact typed factors.
#' @return Explicit hit or miss.
#' @noRd
resolveWorkflowDiaQcProduct <- function(store, identity, stage, plotType) {
  resolveSavedWorkflowEntity(store, "scientific_product", identity,
    factors = list(route = "iq", stage = stage, plot_type = plotType))
}

#' Apply bounded render-only style to a product plot.
#' @param product Existing QC product.
#' @param style Named title and font_size.
#' @return Styled plot without rerunning scientific production.
#' @noRd
styleWorkflowDiaQcPlot <- function(product, style) {
  if (!is.list(product) || !identical(product$schema,
      "prot_norm_qc_plot_product_v1") ||
      !is.list(style) || !identical(sort(names(style)), c("font_size", "title")) ||
      !is.character(style$title) || length(style$title) != 1L || is.na(style$title) ||
      !is.numeric(style$font_size) || length(style$font_size) != 1L ||
      !is.finite(style$font_size) || style$font_size <= 0) {
    cli::cli_abort("invalid DIA QC render style", class = "saved_entity_invalid")
  }
  product_plot <- product[["plot", exact = TRUE]]
  plot <- if (is.null(product_plot)) rebuildWorkflowDiaQcPlot(product) else product_plot
  if (inherits(plot, "patchwork")) {
    plot + patchwork::plot_annotation(title = style$title,
      theme = ggplot2::theme(text = ggplot2::element_text(size = style$font_size)))
  } else {
    plot + ggplot2::labs(title = style$title) +
      ggplot2::theme(text = ggplot2::element_text(size = style$font_size))
  }
}

#' Render one existing QC product to self-contained bytes.
#' @param product Existing QC product, never the source assay.
#' @param style Render-only style.
#' @param format png or pdf.
#' @param dimensions Integer pixel width and height.
#' @param dpi Raster resolution.
#' @param savePlotFn Existing ggsave-compatible renderer.
#' @return Explicit rendered-asset envelope.
#' @noRd
renderWorkflowDiaQcProduct <- function(product, style, format, dimensions,
    dpi = 150, savePlotFn = ggplot2::ggsave) {
  if (!format %in% c("png", "pdf") || !is.integer(dimensions) ||
      length(dimensions) != 2L || anyNA(dimensions) || any(dimensions < 1L) ||
      !is.numeric(dpi) || length(dpi) != 1L || !is.finite(dpi) || dpi <= 0) {
    cli::cli_abort("invalid DIA QC render request", class = "saved_entity_invalid")
  }
  plot <- styleWorkflowDiaQcPlot(product, style)
  path <- tempfile(fileext = paste0(".", format))
  on.exit(unlink(path), add = TRUE)
  savePlotFn(path, plot = plot, width = dimensions[[1L]] / dpi,
    height = dimensions[[2L]] / dpi, dpi = dpi, units = "in")
  if (!file.exists(path) || is.na(file.size(path)) || file.size(path) <= 0) {
    cli::cli_abort("DIA QC renderer produced no asset", class = "saved_entity_write")
  }
  bytes <- readBin(path, "raw", n = file.size(path))
  newSavedWorkflowRenderedAsset(bytes,
    if (identical(format, "png")) "image/png" else "application/pdf",
    dimensions)
}

#' Save one rendered QC asset independently of its source product.
#' @param store Shared store.
#' @param asset Rendered envelope.
#' @param request_id Idempotent request identity.
#' @param render_identity Identity returned by newSavedWorkflowRenderIdentity.
#' @return Accepted rendered-asset descriptor.
#' @noRd
saveWorkflowDiaQcAsset <- function(store, asset, request_id, render_identity) {
  saveWorkflowConsumerEntity(store, asset, request_id, render_identity,
    "rendered_asset")
}

#' One product's root in a snapshot, not yet published.
#' @param manager Current component owner.
#' @param snapshot The head snapshot to add to.
#' @param reference Accepted existing-store product descriptor.
#' @param dependencies What the product depends on.
#' @return The snapshot.
#' @noRd
addWorkflowProductRoot <- function(manager, snapshot, reference, dependencies = list()) {
    snapshot$products[[reference$entity_id]] <- reference
    roots <- collectDIANativeReferences(dependencies)
    for (root in roots) {
        if (identical(root$kind, "checkpoint") && identical(root$payload$codec, "rds")) {
            record <- readDIAStateMetadata(manager, root)
            roots <- collectDIANativeReferences(list(record$data$components, record$data$support,
                record$data$products, record$data$dependencies), roots)
        }
    }
    snapshot$product_dependencies[[reference$entity_id]] <- roots
    snapshot
}

#' Root a reusable scientific or display product in the durable session head. A worker holds no head (KIT-078): an
#' owner whose context carries `deferred_roots` collects the root instead, and the controller's completion publishes it.
#' @param manager Current component owner.
#' @param reference Accepted existing-store product descriptor.
#' @return Descriptor invisibly, without hydrating its source components.
#' @noRd
rootWorkflowProduct <- function(manager, reference, dependencies = list()) {
    if (is.list(manager$.context$deferred_roots)) {
        manager$.context$deferred_roots[[length(manager$.context$deferred_roots) + 1L]] <- list(reference = reference,
            dependencies = collectDIANativeReferences(dependencies))
        return(invisible(reference))
    }
    publishDIAStateHead(manager$.context, addWorkflowProductRoot(manager, manager$.context$snapshot, reference, dependencies))
    invisible(reference)
}

#' DIA name of `rootWorkflowProduct` (KIT-045): the alias the proteomics modules, the harness and the tests keep.
#' @param ... The arguments of `rootWorkflowProduct`.
#' @return What `rootWorkflowProduct` returns.
#' @noRd
rootDIAProduct <- rootWorkflowProduct

#' Identify science using immutable component identities and ordered input roles.
#' @param manager,reference Current owner and source checkpoint.
#' @param operation,parameters Explicit producer and population-affecting parameters.
#' @return Existing store scientific identity, without any assay read.
#' @noRd
workflowComponentScientificIdentity <- function(manager, reference, operation, parameters) {
    spec <- readDIAStateMetadata(manager, reference)$data
    # KIT-031: the class declaration names the assay component (the proteomics identities unchanged: their assay is
    # one table); a `lists` assay (lipidomics, metabolomics) holds one component per named table and is identified by
    # the digest of its sorted element ids (the KIT-070 harness rule), its roles the declaration's.
    declaration <- workflowClassDeclaration(spec$class)
    proteomics <- spec$class %in% c("PeptideQuantitativeData", "ProteinQuantitativeData")
    assay <- if (proteomics) spec$components$protein_quant_table %||% spec$components$peptide_data else
        spec$components[[declaration$assay]]
    input <- if (is.character(assay$entity_id)) assay$entity_id else
        digest::digest(sort(unname(vapply(assay, function(element) element$entity_id, character(1)))), algo = "sha256")
    roles <- if (proteomics) c("protein_id_column", "sample_id", "group_id", "technical_replicate_id") else declaration$roles
    scientific_slots <- spec$slots
    # The DA step's argument records are not identity: the UI block and the publication settings (KIT-004), and the
    # DA defaults the step commits (`differentialAbundanceAnalysis`, its helper), which duplicate the contrast
    # operation's explicit parameters (KIT-022: the commit precedes the contrasts, so the head carries them).
    scientific_slots$args[c("daAnalysisUI", "outputDaResultsAllContrasts", "differentialAbundanceAnalysis",
        "differentialAbundanceAnalysisHelper")] <- NULL
    newSavedWorkflowScientificIdentity(input,
        spec$components$design_matrix$entity_id,
        digest::digest(list(operation, parameters, slots = scientific_slots), algo = "sha256"),
        digest::digest(list(columns = spec$dimensions$columns, roles = spec$slots[roles]), algo = "sha256"),
        manager$.context$software)
}

#' DIA name of `workflowComponentScientificIdentity` (KIT-045): the alias the proteomics modules, the harness and the tests keep.
#' @param ... The arguments of `workflowComponentScientificIdentity`.
#' @return What `workflowComponentScientificIdentity` returns.
#' @noRd
diaComponentScientificIdentity <- workflowComponentScientificIdentity

#' Resolve one normalization display product before assembling any matrix.
#' @param reference Input checkpoint descriptor.
#' @param aesthetics,plot_type Current render choices and existing scientific producer.
#' @param build_product,plot_pca,plot_rle,build_density,build_correlation Existing kernels.
#' @return Plot with a product reference and bounded render identity only.
#' @noRd
resolveDIAReferenceQcPlot <- function(reference, aesthetics, plot_type, build_product,
    plot_pca, plot_rle, build_density, build_correlation) {
    manager <- getDIAActiveManager()
    spec <- readDIAStateMetadata(manager, reference)$data
    # Correlation grouping selects pairs; PCA/RLE colour and shape only render them.
    parameters <- list(plot_type = plot_type, pca_cv_percentile = 0.90,
        correlation_group = if (plot_type == "correlation") aesthetics$color_var else NULL,
        exclude_pool_samples = TRUE, population = if (plot_type == "density")
            "density_independent_pca_call" else paste0(plot_type, "_existing_producer"))
    identity <- diaComponentScientificIdentity(manager, reference, "normalization-qc", parameters)
    hit <- resolveWorkflowDiaQcProduct(manager$.context$store, identity, "component-qc", plot_type)
    if (identical(hit$status, "hit")) {
        product <- hit$value$value
        product_reference <- hit$reference
    } else {
        object <- readDIAExecutionComponents(manager, reference)
        product <- build_product(object, aesthetics, plot_type, plot_pca, plot_rle,
            build_density, build_correlation)
        product$sample_id <- object@sample_id
        product$plot <- NULL
        object <- NULL
        product_reference <- saveWorkflowDiaQcProduct(manager$.context$store, product,
            basename(tempfile("dia-qc-product-")), identity, "component-qc")
        rootDIAProduct(manager, product_reference, list(reference))
    }
    for (field in c(aesthetics$color_var, aesthetics$shape_var)) {
        if (!field %in% names(product$annotations)) cli::cli_abort("unknown QC annotation")
    }
    product$scientific_parameters$color_variable <- aesthetics$color_var
    product$scientific_parameters$shape_variable <- aesthetics$shape_var
    if (plot_type == "rle") {
        at <- match(as.character(product$data$rle.x.factor),
            as.character(product$annotations[[product$sample_id]]))
        product$data$rowinfo <- product$annotations[[aesthetics$color_var]][at]
    }
    plot <- rebuildWorkflowDiaQcPlot(product)
    attr(plot, "dia_product_reference") <- product_reference
    attr(plot, "dia_render_aesthetics") <- aesthetics
    plot
}

#' The normalisation + RUV operation (KIT-004): one declaration producing the existing
#' `dia-normalization-ruv` product. The normalised object is checkpointed before RUV so the corrected
#' object inherits it as its source, exactly as before.
#' @param ruv_design KIT-126. `"legacy"` (the default): today's declaration and kernel path, the corrected matrix
#'   tested by the DA; its ten parameters and identities are unchanged. `"declared"`: the declaration gains
#'   the parameters `ruv_design` (`"groups_all_features"`, `"technical_replicates"` or `"first_pass_nulls"`, KIT-129
#'   section 1e; the replicate sets and controls by `protRuvModelReplicateSets` and `computeProtRuvModelFactors`),
#'   `ruv_k_rule` and `ruv_k_seed`, and the kernel runs KIT-124's arm 7 (groups as replicate sets, every protein a
#'   control), arm 16 (the declared units as replicate sets) or arm 23 (the first pass's controls): k by the shared resolver
#'   (KIT-129 item 3) under `ruv_k_rule`: `num_sv_be` (the estimated count, seeded by `ruv_k_seed`, NULL for the
#'   default seed) or `sweep_bounded` (the sweep's estimate; undefined with every protein a control, so refused) with
#'   `ruv_mode = "automatic"`, `user` (the manual k) with `ruv_mode = "manual"`. The final state's table is the
#'   corrected matrix and it carries the factors and the normalised matrix for the DA (`attachProtRuvModelFactors`);
#'   a count of 0 applies no RUV (the final state is the normalised one, `ruv_result$ruv_applied` FALSE with the
#'   reason). The count's record is `ruv_result$k_estimate` (`plotRuvFactorCount`).
#' @param supplied KIT-129 fifth pass item 3: the declared declaration of the design `supplied_controls`, whose identity
#'   adds the list's hash and count (`RUV_SUPPLIED_LIST_PARAMETERS`); the list itself is the request's `ruv_control_ids`
#'   (`ruvDeclaredSuppliedList`), and the controls are the listed proteins of the table (arm 8 of simulator S1: the
#'   groups as replicate sets). Listed ids that are not proteins of the table are named in a warning
#'   (`ruv_control_ids_unmatched`) and in `ruv_result`; a list that matches none is refused.
#' @noRd
diaNormalizationOperation <- function(ruv_design = c("legacy", "declared"), supplied = FALSE) {
    ruv_design <- match.arg(ruv_design)
    factors_in_model <- identical(ruv_design, "declared")
    supplied <- factors_in_model && isTRUE(supplied)
    newWorkflowOperation("normalization-ruv", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        output_class = "ProteinQuantitativeData", product_type = "dia-normalization-ruv",
        parameters = c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
            "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k", "grouping",
            if (factors_in_model) c("ruv_design", "ruv_k_rule", "ruv_k_seed"), if (supplied) RUV_SUPPLIED_LIST_PARAMETERS),
        kernel = function(object, parameters, store) {
            # The kernels record their arguments under the calling name; the legacy path calls them through
            # these aliases, so the aliases keep the recorded args (and every checkpoint hash) identical.
            normaliseBetweenSamplesFn <- normaliseBetweenSamples
            # KIT-126 audit F-06: today's declaration does not know the design; handed it, it refuses rather than
            # run today's design under a record that names the other.
            if (!factors_in_model) assertProtNormLegacyRuvParameters(parameters)
            # Audit F-08: under the new design a manual k is one whole number, refused otherwise (not truncated, no default).
            if (factors_in_model && identical(parameters$ruv_mode, "manual")) validateProtRuvModelK(parameters$ruv_k)
            if (factors_in_model) assertProtRuvModelKRule(parameters)
            # KIT-129 fifth pass item 3: the supplied list, checked against the request's identity of it.
            listed <- if (factors_in_model) ruvDeclaredSuppliedList(parameters)
            if (!is.null(listed) && !supplied) {
                cli::cli_abort("the design supplied_controls runs under its own declaration", class = "workflow_execution_invalid")
            }
            # KIT-129 integration item 2: the option set's None is the design `none`, with the step's skip mode only.
            none <- factors_in_model && identical(parameters$ruv_design, "none")
            if (none && !identical(parameters$ruv_mode, "skip"))
                cli::cli_abort("the design none runs with the skip mode only", class = "workflow_execution_invalid")
            normalized <- normaliseBetweenSamplesFn(object, normalisation_method = parameters$norm_method)
            normalized_reference <- store$checkpoint(normalized, "normalized-science")
            object <- NULL
            if (parameters$ruv_mode == "skip") {
                ruv_result <- buildProtNormSkippedRuvResult()
                final_reference <- normalized_reference
                # The new option set's None declares itself, so that its DA is the new design's (decision M30).
                if (none) final_reference <- store$checkpoint(declareProtRuvNoFactors(normalized, list(ruv_design = "none",
                    reason = "the design none: no correction for unwanted variation")), "ruv-none-science", from = normalized_reference)
                best_k <- NA
                controls <- NA
                best_percentage <- NA
            } else {
                resolved <- resolveProtRuvParameters(normalized, parameters, parameters$grouping)
                ruv_result <- resolved$result
                if (!is.null(ruv_result$best_cancor_plot)) {
                    ruv_result$best_cancor_plot <- ggplot2::ggplotGrob(
                        stripProtNormQcPlotOwners(ruv_result$best_cancor_plot))
                }
                ruvIII_C_VaryingFn <- ruvIII_C_Varying
                removeRowsWithMissingValuesPercentFn <- removeRowsWithMissingValuesPercent
                if (factors_in_model) {
                    # KIT-126: the sweep (or the manual k) gives k only; every protein is a control, the groups are
                    # the replicate sets, and the corrected matrix uses the same controls, k and sets as the factors.
                    # KIT-129 third pass item 3: the technical-replicate design (the declared units as the replicate
                    # sets) and first-pass nulls (its controls; the replicate sets the design resolves) as well.
                    replicate_sets <- protRuvModelReplicateSets(normalized, parameters$ruv_design, parameters$grouping)
                    # KIT-129 fifth pass item 3: the listed proteins of the table are the controls (arm 8).
                    list_controls <- NULL
                    unmatched <- character()
                    if (!is.null(listed)) {
                        protein_ids <- as.character(normalized@protein_quant_table[[normalized@protein_id_column]])
                        if (!any(listed %in% protein_ids)) {
                            cli::cli_abort(sprintf("none of the %d listed control identifiers is a protein of this data set", length(listed)),
                                class = "workflow_execution_invalid", call = NULL)
                        }
                        unmatched <- warnRuvSuppliedListUnmatched(listed, protein_ids, "proteins")
                        list_controls <- protein_ids %in% listed
                    }
                    # KIT-129 fifth pass item 4: with a supplied list the separation rule reads that list (arm 22 of
                    # simulator S1), not the sweep's own controls; its plot is the one shown.
                    if (!is.null(listed) && identical(parameters$ruv_k_rule, "sweep_bounded")) {
                        separated <- protRuvSeparationRuleK(normalized, list_controls, parameters$grouping)
                        resolved$k <- separated$k
                        ruv_result$best_cancor_plot <- ggplot2::ggplotGrob(stripProtNormQcPlotOwners(separated$plot))
                    }
                    ruv_model <- computeProtRuvModelFactors(normalized,
                        ruv_k = if (!identical(parameters$ruv_k_rule, "num_sv_be")) resolved$k, ruv_k_rule = parameters$ruv_k_rule,
                        controls = list_controls,
                        ruv_grouping_variable = replicate_sets, ruv_k_seed = parameters$ruv_k_seed %||% RUV_FACTOR_COUNT_SEED,
                        ruv_design = switch(parameters$ruv_design, first_pass_nulls = "first_pass_nulls", supplied_controls = "supplied_controls",
                            "groups_all_features"))
                    applied <- identical(ruv_model$status, "ok")
                    resolved$k <- ruv_model$k
                    resolved$controls <- if (applied) ruv_model$controls else NA
                    resolved$percentage <- if (applied) 100 * mean(ruv_model$controls) else NA
                    # Audit F-05: the result every reader takes (the panel text, the study parameters, the optimisation
                    # RDS, the reports) names what was used: every protein a control, the k of the factors. The
                    # control selection's own values (the sweep's, or the manual percentage's) stay as a diagnostic;
                    # the canonical-correlation plot and the optimisation table are that selection's too.
                    if (!is.null(ruv_result)) ruv_result$sweep_diagnostic <-
                        ruv_result[c("best_percentage", "best_k", "best_control_genes_index")]
                    ruv_result$best_percentage <- resolved$percentage
                    ruv_result$best_k <- ruv_model$k
                    ruv_result$best_control_genes_index <- if (applied) ruv_model$controls
                    ruv_result$ruv_design <- protRuvModelDesign(ruv_model)
                    # KIT-129 item 3: the rule that gave k, the count's record (for `plotRuvFactorCount`) and whether RUV
                    # was applied, by the names the other omics' products use.
                    ruv_result$ruv_k_rule <- ruv_model$ruv_k_rule
                    ruv_result$ruv_k_estimate <- ruv_model$ruv_k_estimate
                    ruv_result$k_estimate <- ruv_model$k_estimate
                    ruv_result$ruv_applied <- applied
                    ruv_result$ruv_not_applied_reason <- ruv_model$reason
                    # KIT-129 fifth pass item 3: the supplied list by section 1e's names (its hash, the ids listed, those not
                    # found); the controls used are `ruv_model$control_count`.
                    if (!is.null(listed)) {
                        ruv_result$ruv_control_rule <- "supplied_list"
                        ruv_result[c("ruv_control_list_sha256", "ruv_control_list_count")] <- list(parameters$ruv_control_list_sha256,
                            length(listed))
                        ruv_result$ruv_control_list_unmatched <- unmatched
                    }
                    # KIT-129 integration item 3: the count beside a k the user chose (section 1e), once, here; and the
                    # seed the count was (or would be) made with.
                    ruv_result$ruv_k_seed <- as.integer(parameters$ruv_k_seed %||% RUV_FACTOR_COUNT_SEED)
                    if (identical(parameters$ruv_k_rule, "user")) {
                        samples <- as.character(normalized@design_matrix[[normalized@sample_id]])
                        technical <- resolveProtTechnicalReplicateBlock(normalized@design_matrix, normalized@group_id,
                            normalized@technical_replicate_id)
                        beside <- recordRuvFactorCountBesideUserK(as.matrix(normalized@protein_quant_table[, samples, drop = FALSE]),
                            normalized@design_matrix[[parameters$grouping]], ruv_result$ruv_k_seed,
                            units = if (isTRUE(technical$technical)) technical$block)
                        ruv_result$ruv_k_estimate <- beside$k_estimate
                        ruv_result$ruv_k_estimate_reason <- beside$k_estimate_reason
                        ruv_result$k_estimate <- beside$estimate
                    }
                    if (applied) {
                        corrected <- ruvIII_C_VaryingFn(protRuvGroupReplicateObject(normalized),
                            ruv_grouping_variable = replicate_sets, ruv_number_k = resolved$k, ctrl = resolved$controls)
                        corrected <- removeRowsWithMissingValuesPercentFn(theObject = corrected)
                        corrected <- attachProtRuvModelFactors(corrected, normalized, ruv_model)
                    }
                    # KIT-129 fourth pass item 4: the factor-tied-to-groups check on this run's W (decision M27).
                    ruv_result$ruv_group_association <- ruvGroupAssociationCheck(ruv_model$W, ruvFactorGroups(ruv_model$W,
                        normalized@design_matrix, normalized@sample_id, normalized@group_id), "proteomics", parameters$ruv_design,
                        ruv_model$ruv_k_rule, ruv_model$k, ruv_result$ruv_k_estimate %||% NA_integer_)
                    ruv_model$W <- ruv_model$controls <- ruv_model$k_estimate <- NULL
                } else {
                corrected <- ruvIII_C_VaryingFn(normalized,
                    ruv_grouping_variable = parameters$grouping, ruv_number_k = resolved$k, ctrl = resolved$controls)
                corrected <- removeRowsWithMissingValuesPercentFn(theObject = corrected)
                }
                # A count of 0 under the new design applies no RUV: the final state is the normalised one, declaring
                # the design without factors (KIT-129 integration item 2), so that its DA is the new design's.
                final_reference <- if (factors_in_model && !applied) {
                    store$checkpoint(declareProtRuvNoFactors(normalized, c(list(ruv_design = parameters$ruv_design),
                        ruv_model[c("ruv_k_rule", "ruv_k_estimate", "reason")])), "ruv-none-science", from = normalized_reference)
                } else store$checkpoint(corrected, "ruv-science", from = normalized_reference)
                corrected <- NULL
                best_k <- resolved$k
                controls <- resolved$controls
                best_percentage <- resolved$percentage
            }
            normalized <- NULL
            value <- list(normalized = normalized_reference, final = final_reference, ruv_result = ruv_result,
                best_k = best_k, controls = controls, percentage = best_percentage,
                parameters = parameters[c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
                    "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k", "grouping")])
            # Audit F-08: a skipped RUV step applied no design, so none is recorded on its product; the design none is
            # recorded as itself (KIT-129 integration item 2).
            if (none) value$parameters[c("ruv_design", "ruv_k_rule", "ruv_k_seed")] <- parameters[c("ruv_design", "ruv_k_rule", "ruv_k_seed")]
            if (factors_in_model && parameters$ruv_mode != "skip") {
                value$parameters[c("ruv_design", "ruv_k_rule", "ruv_k_seed")] <- parameters[c("ruv_design", "ruv_k_rule", "ruv_k_seed")]
                if (supplied) value$parameters[RUV_SUPPLIED_LIST_PARAMETERS] <- parameters[RUV_SUPPLIED_LIST_PARAMETERS]
                value$ruv_model <- c(list(ruv_design = protRuvModelDesign(ruv_model)), ruv_model)
            }
            value
        },
        # KIT-187 (VAL-71, D2): the clean-up after RUV decides by counts; a recorded step without this marker replays
        # with the earlier rounded-percentage comparison (`workflowReplayRuleMarkers`).
        fixed = list(missing_value_comparison = "counts"))
}

#' Today's normalisation declaration takes none of the new design's parameters (KIT-126 audit F-06, KIT-129 B1): its
#' identity reads its ten parameters only, so `runSavedOperation` calls this before the identity, and a request that
#' names the new design is refused rather than served today's saved product. The condition is the legacy refusal of
#' every omic (`workflow_execution_invalid`).
#' @param parameters A request's parameters.
#' @return `parameters`, invisibly.
#' @noRd
assertProtNormLegacyRuvParameters <- function(parameters) {
    if (!is.null(parameters[["ruv_design", exact = TRUE]])) {
        cli::cli_abort(c("the normalisation declaration of today's design was handed ruv_design = {.val {parameters$ruv_design}}",
            i = "declare it with diaNormalizationOperation(\"declared\")"), class = "workflow_execution_invalid", call = NULL)
    }
    given <- Filter(\(name) !is.null(parameters[[name, exact = TRUE]]), c(PROT_NORM_RUV_K_PARAMETERS, "ruv_control_ids",
        RUV_SUPPLIED_LIST_PARAMETERS))
    if (length(given)) {
        cli::cli_abort(c("the normalisation declaration of today's design takes no k rule; the request carries {.val {given}}",
            i = "declare it with diaNormalizationOperation(\"declared\")"), class = "workflow_execution_invalid", call = NULL)
    }
    invisible(parameters)
}

#' The new design's k parameters of the classic normalisation (KIT-129 item 3), after `ruv_design`.
#' @noRd
PROT_NORM_RUV_K_PARAMETERS <- c("ruv_k_rule", "ruv_k_seed")

#' The k rule of the classic new design agrees with the step's mode (KIT-129 item 3): the user's k under `manual`, the
#' estimated count or the separation rule under `automatic`; a seed is a whole number or NULL (the default seed).
#' @param parameters The declaration's parameters.
#' @return `parameters`, invisibly; an error of class `workflow_execution_invalid` otherwise.
#' @noRd
assertProtRuvModelKRule <- function(parameters) {
    invalid <- function(message) cli::cli_abort(message, class = "workflow_execution_invalid", call = NULL)
    rule <- parameters$ruv_k_rule
    if (identical(parameters$ruv_mode, "skip")) return(invisible(parameters))
    if (!checkmate::test_choice(rule, RUV_K_RULES)) invalid(paste("ruv_k_rule must be one of", paste(RUV_K_RULES, collapse = ", ")))
    if (identical(parameters$ruv_mode, "manual") != identical(rule, "user")) {
        invalid("ruv_k_rule is user with ruv_mode manual, and num_sv_be or sweep_bounded with ruv_mode automatic")
    }
    if (!is.null(parameters$ruv_k_seed) && !checkmate::test_int(parameters$ruv_k_seed)) invalid("ruv_k_seed must be a whole number")
    invisible(parameters)
}

#' Run or reuse the normalisation + RUV declaration and publish its outputs.
#' @param settings Plain list of the normalisation settings (`norm_method`, `ruv_mode`, the RUV fields,
#'   `color_variable`, `shape_variable`); the module edge builds it from its inputs.
#' @param workflow,norm_data,paths Existing GUI owners and publication paths.
#' @param omic_type,experiment_label Existing publication labels.
#' @param post_qc,ruv_qc,grouping Existing callbacks.
#' @return TRUE after science, required tables and plots have completed.
#' @noRd
#' Save a display table (KIT-018): a bounded data frame the module renders, saved once as a result table
#' (content-addressed, so a repeated save is a hit) and paged on demand by `readDIADisplayTablePage`.
#' @param manager Session owner.
#' @param table Data frame; NULL or zero rows save nothing.
#' @return Table reference, or NULL.
#' @noRd
saveWorkflowDisplayTable <- function(manager, table) {
    if (!is.data.frame(table) || !nrow(table) || !ncol(table)) return(NULL)
    # Only atomic columns are displayable and saveable (a g:Profiler result carries list columns such as
    # `parents`); the module's download handlers keep the full frame.
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    table <- table[vapply(table, function(column) is.atomic(column) && is.null(dim(column)), logical(1))]
    if (!ncol(table)) return(NULL)
    saveDIATableComponent(manager$.context, table)
}

#' DIA name of `saveWorkflowDisplayTable` (KIT-045): the alias the proteomics modules, the harness and the tests keep.
#' @param ... The arguments of `saveWorkflowDisplayTable`.
#' @return What `saveWorkflowDisplayTable` returns.
#' @noRd
saveDIADisplayTable <- saveWorkflowDisplayTable

#' Read one page of a saved display table (KIT-018), optionally restricted to the rows whose `filter$column`
#' equals `filter$value` (the column is scanned in 10,000-row reads; the selected rows are read by ordinal).
#' @param reference Table reference from `saveDIADisplayTable`.
#' @param columns Columns to read; NULL for every column the table has.
#' @param filter NULL or `list(column, value)`.
#' @param offset,limit Page position and size (at most 10,000 rows per read).
#' @return `list(data, rows)`: the page and the matching row count.
#' @noRd
readWorkflowDisplayTablePage <- function(reference, columns = NULL, filter = NULL, offset = 0L, limit = 1000L,
    manager = getDIAActiveManager()) {
    # KIT-017 screen run: the widget's page requests are served by Shiny's request handler, which carries no
    # reactive domain and no open scope, so a render resolves the owner once and the pager passes it here.
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("no state owner for the display table page", class = "dia_owner_missing")
    store <- manager$.context$store
    if (!isWorkflowTableReference(reference)) cli::cli_abort("not a saved display table", class = "saved_entity_invalid")
    scope <- newWorkflowTableReadScope(store, reference, "duckdb")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    if (is.null(columns)) columns <- scope$metadata$columns
    rows <- reference$metadata$shape[[1L]]
    read <- function(columns, offset, limit) readSavedWorkflowTable(store, reference, columns, offset = offset,
        limit = limit, scope = scope, reader = "duckdb")$data
    if (is.null(filter)) return(list(data = read(columns, offset, limit), rows = rows))
    selected <- integer()
    at <- 0L
    while (at < rows) {
        count <- min(10000L, rows - at)
        part <- read(filter$column, at, count)
        selected <- c(selected, at + which(part[[filter$column]] == filter$value))
        at <- at + count
    }
    chosen <- selected[seq.int(offset + 1L, length.out = max(0L, min(limit, length(selected) - offset)))]
    page <- if (length(chosen)) readWorkflowSelectedOrdinals(store, reference, columns, chosen, scope = scope) else
        read(columns, 0L, 1L)[0L, , drop = FALSE]
    list(data = page, rows = length(selected))
}

#' DIA name of `readWorkflowDisplayTablePage` (KIT-045): the alias the proteomics modules, the harness and the tests keep.
#' @param ... The arguments of `readWorkflowDisplayTablePage`.
#' @return What `readWorkflowDisplayTablePage` returns.
#' @noRd
readDIADisplayTablePage <- readWorkflowDisplayTablePage

#' A server-side page filter for the table widget over a saved table (KIT-024): the widget asks for a page
#' (`start`, `length`, `draw`) and the store serves it through `read_page(offset, limit)`; ordering and searching
#' are not applied (the renders disable them), so the row count is the table's. Without an installed pager
#' (`pager$read` NULL: an in-memory table) the widget's own filter runs.
#' @param pager Environment with `read` (`function(offset, limit)` returning a data frame) and `rows`.
#' @return `function(data, params)` for `DT::renderDT(funcFilter = )`.
#' @noRd
workflowDisplayTableFilter <- function(pager) {
    default_filter <- utils::getFromNamespace("dataTablesFilter", "DT")
    # The widget reads each row as an array, so the page is unnamed the way the widget's own filter unnames it
    # (KIT-017 screen run: a named page serialised as objects and the widget's error dialog blocked the browser).
    clean <- utils::getFromNamespace("cleanDataFrame", "DT")
    function(data, params) {
        if (is.null(pager$read)) return(default_filter(data, params))
        rows <- as.integer(pager$rows)
        start <- max(0L, as.integer(params$start %||% 0L))
        length <- as.integer(params$length %||% 10L)
        if (is.na(length) || length < 0L) length <- rows
        length <- min(max(length, 1L), 10000L)
        page <- if (start < rows) pager$read(start, length) else pager$read(0L, 1L)[0L, , drop = FALSE]
        page <- as.data.frame(lapply(page, function(column)
            if (is.factor(column) || inherits(column, c("Date", "POSIXt"))) as.character(column) else column),
            stringsAsFactors = FALSE, check.names = FALSE)
        rownames(page) <- NULL
        if (!identical(params$escape, "false")) for (j in seq_len(ncol(page)))
            if (is.character(page[[j]])) page[[j]] <- htmltools::htmlEscape(page[[j]])
        list(draw = as.integer(params$draw %||% 0L), recordsTotal = rows, recordsFiltered = rows, data = clean(page),
            DT_rows_all = seq_len(rows), DT_rows_current = seq.int(start + 1L, length.out = nrow(page)))
    }
}

#' DIA name of `workflowDisplayTableFilter` (KIT-045): the alias the proteomics modules, the harness and the tests keep.
#' @param ... The arguments of `workflowDisplayTableFilter`.
#' @return What `workflowDisplayTableFilter` returns.
#' @noRd
diaDisplayTableFilter <- workflowDisplayTableFilter

#' The table widget options the paged renders share (KIT-024): server-side pages, no client ordering or search
#' (the store serves rows in source order), a visible row count.
#' @noRd
diaDisplayTableOptions <- function(page_length = 25L, ...) list(pageLength = page_length, ordering = FALSE,
    searching = FALSE, serverSide = TRUE, processing = TRUE, ...)

#' The enrichment method an `EnrichmentResults` object was produced with (KIT-024): `clusterprofiler` when its
#' entries are `enrichResult` objects, else `gprofiler2`.
#' @noRd
diaEnrichmentMethodOf <- function(results) {
    for (contrast in names(results@enrichment_data)) for (direction in names(results@enrichment_data[[contrast]])) {
        entry <- results@enrichment_data[[contrast]][[direction]]
        if (isS4(entry) && methods::is(entry, "enrichResult")) return("clusterprofiler")
        if (is.list(entry)) return("gprofiler2")
    }
    "gprofiler2"
}

#' The RUV optimisation result of the accepted normalisation (KIT-018): read from the saved normalisation
#' product the current normalised state references, never from a module copy or an RDS file.
#' @param manager Session owner; NULL or a non-reference owner gives NULL.
#' @return The `ruv_result` list the normalisation declaration produced, or NULL before normalisation.
#' @noRd
resolveDIARuvOptimizationResult <- function(manager) {
    if (!usesDIAReferenceStates(manager)) return(NULL)
    history <- manager$getHistory()
    for (state in intersect(c("correlation_filtered", "ruv_corrected", "normalized"), history)) {
        record <- tryCatch(readDIAStateMetadata(manager, getDIAStateReference(manager, state))$data, error = function(e) NULL)
        product <- record$products$normalization_science
        if (!is.null(product)) return(hydrateSavedWorkflowEntity(manager$.context$store, product)$value$ruv_result)
    }
    NULL
}

#' The RUV mode and grouping of the run that produced the accepted RUV result (KIT-122, audit F2): the declared
#' parameters its saved product carries (`grouping` for the normalisation, `ruv_grouping_variable` for the DIA_LIMPA
#' covariates), read from the head state's product, else from the newest state that references one (after the DIA
#' route's correlation filter, or an Open). The tab's RUV summary reads this, never the tab's inputs.
#' @param manager Session owner; NULL or a non-reference owner gives NULL.
#' @return `list(ruv_mode, grouping)`, or NULL before normalisation.
#' @noRd
resolveDIARuvRunDescription <- function(manager) {
    if (!usesDIAReferenceStates(manager)) return(NULL)
    states <- unique(c(manager$current_state, rev(manager$getHistory())))
    for (state in intersect(states, c("ruv_covariates", "correlation_filtered", "ruv_corrected", "normalized"))) {
        record <- tryCatch(readDIAStateMetadata(manager, getDIAStateReference(manager, state))$data, error = function(e) NULL)
        product <- record$products$ruv_covariates_science %||% record$products$normalization_science
        if (is.null(product)) next
        value <- hydrateSavedWorkflowEntity(manager$.context$store, product)$value
        parameters <- value$parameters
        # KIT-129 (audit F2): a declared-rule run that applied no RUV (a count of 0, or a refused automatic k) is
        # described as such, whatever mode its parameters asked for.
        ruv_mode <- if (identical(value$ruv_applied, FALSE)) "none" else parameters$ruv_mode
        return(list(ruv_mode = ruv_mode, grouping = parameters$grouping %||% parameters$ruv_grouping_variable))
    }
    NULL
}

#' The in-process normalisation: the gate, the request, the compute half, the commit and the publication half.
#' @param pre_qc KIT-182 (F-e): the tab's pre-normalisation figure generator, called when the run's source differs from
#'   the table the figures were drawn from (`diaNormalizationRequest`'s `redraw_pre`), before the publication copies
#'   them; NULL draws nothing, as before.
#' @noRd
runDIAReferenceNormalization <- function(settings, workflow, norm_data, paths, omic_type,
    experiment_label, post_qc, ruv_qc, grouping, pre_qc = NULL) {
    prior_plot <- ggplot2::get_last_plot()
    on.exit(ggplot2::set_last_plot(prior_plot), add = TRUE)
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "normalization")
        request <- diaNormalizationRequest(manager, settings, norm_data, grouping())
        # KIT-257 (VAL-250): the earlier run's own files go before anything is drawn, computed or committed.
        removed <- removeProtNormOwnFiles(paths, "normalization")
        # AUDIT-KIT-257-impl I05: a failure after the removal carries the removed files, for the tab's notice.
        withProtNormRemovedFiles(removed, {
            if (isTRUE(request$redraw_pre) && is.function(pre_qc)) pre_qc()
            run <- computeDIAReferenceNormalization(manager, request$source, request$parameters)
            commit <- diaNormalizationCommit(manager, settings, run)
            if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            publishDIAReferenceNormalization(manager, settings, run, workflow, norm_data, paths, omic_type, experiment_label,
                post_qc, ruv_qc)
        })
    })
}

# --- KIT-182 (VAL-48, VAL-91, VAL-99): the normalisation reads the current quality-control result ---

#' The normalisation tab's own states (KIT-182): pressed there, "Run Normalization & RUV" re-runs from the table the
#' previous run read (the restore treats the same three names as "normalisation active").
#' @noRd
DIA_NORMALISATION_TAB_STATES <- c("normalized", "ruv_corrected", "correlation_filtered")

#' The save labels that mark a table as normalised (KIT-182): the tab's own states and the normalisation's
#' between-sample checkpoint, saved by the operation's store as `owner:normalized-science`.
#' @noRd
DIA_NORMALISED_LABELS <- c(DIA_NORMALISATION_TAB_STATES, "owner:normalized-science")

#' The heads on which a protein quality-control step is refused by name (KIT-182, VAL-91 D1): the normalisation's and the
#' DIA_LIMPA RUV covariates' states, the correlation filter's, and the result table enrichment commits after DA.
#' @noRd
DIA_QC_REFUSED_HEADS <- c("normalized", "ruv_corrected", "ruv_covariates", "correlation_filtered", "enrichment_completed")

#' The quality-control steps as the Quality Control tab names them (KIT-182), for the refusals' "Revert until" step.
#' @noRd
DIA_QC_STEP_TITLES <- c(protein_s4_initial = "the design step (before quality control)",
    protein_s4_created = "IQ Protein Rollup", protein_accession_cleaned = "Accession Cleanup",
    protein_intensity_filtered = "Protein Intensity Filter", duplicates_removed = "Duplicate Removal",
    protein_replicate_filtered = "Protein Replicate Filter", protein_dpc_quant = "DPC Quantification",
    protein_normalised = "Between-Sample Normalisation", peptide_count_filtered = "Peptide Count Filter")

#' A saved state's bounded record for the lineage checks (KIT-182), or NULL when it is not a readable component
#' checkpoint (a legacy payload, or no component metadata): logged, and the caller keeps today's outcome.
#' @param manager Session owner.
#' @param reference Saved state descriptor, or NULL.
#' @return The record's `data` (the component checkpoint), or NULL.
#' @noRd
readDIALineageRecord <- function(manager, reference) {
    if (!isSavedWorkflowReference(reference)) return(NULL)
    unread <- function(condition) {
        logger::log_info(paste("KIT-182: a state record was not read for the normalisation lineage:",
            conditionMessage(condition)))
        NULL
    }
    tryCatch(readComponentMetadata(manager, reference)$data,
        saved_entity_unsupported = unread, saved_entity_corrupt = unread)
}

#' Whether a protein state has normalisation in its lineage (KIT-182, the ticket's 7.1): (i) its record carries the
#' normalisation product; or (ii) a reference among its record's cumulative dependencies was saved under a
#' normalisation label (`DIA_NORMALISED_LABELS`: the label every descriptor keeps unaltered); or (iii), for the head
#' only, the head was adopted from a legacy session (`adopted` provenance) with `normalized` or `ruv_corrected` earlier
#' in the active history.
#' @param manager Session owner.
#' @param spec The state's record (`readDIALineageRecord`).
#' @param head TRUE when `spec` is the head's record, which enables clause (iii).
#' @return TRUE or FALSE.
#' @noRd
diaNormalisedLineage <- function(manager, spec, head = FALSE) {
    if (!is.null(spec$products$normalization_science)) return(TRUE)
    labelled <- vapply(spec$dependencies, function(dependency)
        is.list(dependency) && isTRUE(dependency$identity$parameters %in% DIA_NORMALISED_LABELS), logical(1))
    if (any(labelled)) return(TRUE)
    if (!head || !isTRUE(manager$.context$snapshot$step_provenance[[manager$current_state]]$adopted)) return(FALSE)
    history <- diaStateHistory(manager)
    any(match(c("normalized", "ruv_corrected"), history) < match(manager$current_state, history), na.rm = TRUE)
}

#' The last quality-control step before normalisation, as the Quality Control tab names it (KIT-182, the ticket's 7.3):
#' the active history walked back from the position before its first normalisation-tab state (from the head's position
#' when there is none), to the first state whose record satisfies neither lineage clause (i) nor (ii).
#' @param manager Session owner.
#' @return The step's title (`DIA_QC_STEP_TITLES`), or "the peptide quality-control steps" for any other state.
#' @noRd
diaPreNormalisationStepName <- function(manager) {
    history <- diaStateHistory(manager)
    tab <- which(history %in% c("normalized", "ruv_corrected", "ruv_covariates", "correlation_filtered"))
    start <- if (length(tab)) tab[[1L]] - 1L else match(manager$current_state, history, nomatch = length(history))
    title <- function(state) if (isTRUE(state %in% names(DIA_QC_STEP_TITLES))) DIA_QC_STEP_TITLES[[state]] else
        "the peptide quality-control steps"
    for (state in rev(history[seq_len(start)])) {
        spec <- readDIALineageRecord(manager, getDIAStateReference(manager, state))
        if (is.null(spec) || !diaNormalisedLineage(manager, spec)) return(title(state))
    }
    title(NA_character_)
}

#' The refusal of "Run Normalization & RUV" on a normalised input (KIT-182, VAL-91 D1; wording proposed, VAL-41).
#' @param manager Session owner.
#' @noRd
abortDIANormalisedInput <- function(manager) {
    rlang::abort(paste0("this normalisation would start from a table that has already been normalised ",
        "(quality-control steps were applied after normalisation, or the project recorded a normalised table as its ",
        "starting table), so it cannot be normalised again. In the Quality Control tab press Revert until the current ",
        "step is \"", diaPreNormalisationStepName(manager), "\", the last quality-control result before normalisation, ",
        "apply again from there any quality-control steps you want, then press \"Run Normalization & RUV\"."),
        class = "prot_norm_input_after_normalisation")
}

#' Today's chain for a press at the tab's own state (KIT-182 `guard`): a source saved under a normalisation-tab state's
#' name, or whose record satisfies lineage clause (i) or (ii), is refused; a NULL source, a value that is not a saved
#' reference or an unreadable record passes unchanged (today's outcome).
#' @param manager Session owner.
#' @param source The candidate source.
#' @return `source`.
#' @noRd
guardDIANormalizationSource <- function(manager, source) {
    if (!isSavedWorkflowReference(source)) return(source)
    normalised <- isTRUE(source$identity$parameters %in% DIA_NORMALISATION_TAB_STATES)
    if (!normalised) {
        spec <- readDIALineageRecord(manager, source)
        normalised <- !is.null(spec) && diaNormalisedLineage(manager, spec)
    }
    if (normalised) abortDIANormalisedInput(manager)
    source
}

#' The normalisation's source (KIT-182, VAL-48; the ticket's 7.1). At the tab's own state (`normalized`,
#' `ruv_corrected`, `correlation_filtered`; checked first) the table the previous run read, as before (R3), so a
#' settings change never normalises a normalised table. Otherwise the head when it is a protein or peptide state with no
#' normalisation in its lineage (R1: the current quality-control result; a peptide head meets the operation's existing
#' refusal), refused when it has (R2: quality control applied after normalisation, in a session saved before the fix or
#' by a script). Any other head (a result table, an unreadable record): today's chain (R3). Every R3 result is guarded.
#' A NULL or non-reference owner keeps today's chain unguarded.
#' @param manager,norm_data Session owner and the tab's state.
#' @return The source reference.
#' @noRd
diaNormalizationSource <- function(manager, norm_data) {
    today <- function() norm_data$normalization_input_reference %||%
        manager$.context$snapshot$qc_sources$pre_normalization %||%
        getDIAStateReference(manager, "protein_replicate_filtered") %||% getDIAStateReference(manager)
    if (!usesDIAReferenceStates(manager)) return(today())
    if (manager$current_state %in% DIA_NORMALISATION_TAB_STATES) return(guardDIANormalizationSource(manager, today()))
    head <- getDIAStateReference(manager)
    spec <- readDIALineageRecord(manager, head)
    if (!is.null(spec) && isTRUE(spec$class %in% c("ProteinQuantitativeData", "PeptideQuantitativeData"))) {
        if (!diaNormalisedLineage(manager, spec, head = TRUE)) return(head)
        abortDIANormalisedInput(manager)
    }
    guardDIANormalizationSource(manager, today())
}

#' The refusal of a protein quality-control Apply on a normalised table (KIT-182, VAL-91 D1). It decides by name
#' first (`DIA_QC_REFUSED_HEADS`); only then does it read the head's record (`readComponentMetadata` through
#' `readDIALineageRecord`, before any read of the history or provenance) and refuse a protein state with normalisation
#' in its lineage (a session saved before the fix). A no-op on a NULL or non-reference owner. Called first by the six
#' protein QC sub-tabs' steps in process and before their managed `apply`; not by the kernels, the registry handlers or
#' replay.
#' @param workflowData The workflow bus (`state_manager`).
#' @return TRUE invisibly, or an error of class `prot_qc_input_normalised`.
#' @noRd
assertProtQcHeadNotNormalised <- function(workflowData) {
    manager <- workflowData$state_manager
    if (is.null(manager) || !usesDIAReferenceStates(manager)) return(invisible(TRUE))
    refuse <- function() rlang::abort(paste0("quality-control steps cannot be applied to a table that has already ",
        "been normalised. To change quality control, press Revert in this tab until the current step is \"",
        diaPreNormalisationStepName(manager), "\", the last quality-control result before normalisation, then apply ",
        "the steps again and press \"Run Normalization & RUV\" again."), class = "prot_qc_input_normalised")
    if (isTRUE(manager$current_state %in% DIA_QC_REFUSED_HEADS)) refuse()
    spec <- readDIALineageRecord(manager, getDIAStateReference(manager))
    if (identical(spec$class, "ProteinQuantitativeData") && diaNormalisedLineage(manager, spec, head = TRUE)) refuse()
    invisible(TRUE)
}

#' The managed protein QC sub-tab's check before `apply` (KIT-182, VAL-91 D1): the refusal, or any error reading the
#' head, is shown as the managed step's own failures are ("Error applying <label>: ...") and the step is not started.
#' @param workflowData The workflow bus.
#' @param label The label the sub-tab gives `installWorkflowManagedStep`.
#' @return TRUE when `apply` may run.
#' @noRd
protQcManagedApplyAllowed <- function(workflowData, label) {
    tryCatch({
        withDIAStateScope(workflowData$state_manager, function() assertProtQcHeadNotNormalised(workflowData))
        TRUE
    }, error = function(e) {
        message <- paste0("Error applying ", label, ": ", conditionMessage(e))
        logger::log_error(message)
        shiny::showNotification(message, type = "error", duration = 15)
        FALSE
    })
}

#' The DIA_LIMPA pass-through's refusal of a normalised table (KIT-182, VAL-91 D1, second part; the third re-check's
#' C02): the table it would designate (`diaLimpaPassThroughSource`) is a protein state with normalisation in its lineage,
#' as in a session saved before the fix in which Between-Sample Normalisation or Peptide Count Filter ran on
#' `ruv_corrected`, `ruv_covariates` or `correlation_filtered`. KIT-244: a peptide table (the press before DPC
#' Quantification) is refused first, in plain words, before anything is designated, drawn or queued. A no-op on a
#' non-reference owner.
#' @param manager,normData Session owner and the tab's state.
#' @return TRUE invisibly, or an error of class `prot_norm_limpa_peptide_input` or
#'   `prot_norm_input_after_normalisation`.
#' @noRd
assertDIALimpaDesignationNotNormalised <- function(manager, normData) {
    if (!usesDIAReferenceStates(manager)) return(invisible(TRUE))
    source <- diaLimpaPassThroughSource(manager, normData)
    spec <- readDIALineageRecord(manager, source)
    # KIT-244 (draft wording, VAL-41): names no state, so it holds at every peptide head.
    if (identical(spec$class, "PeptideQuantitativeData")) rlang::abort(paste0("on DIA with limpa, \"Run Normalization & ",
        "RUV\" needs the protein table that DPC Quantification makes, and the current table is still a peptide table. ",
        "In the Quality Control tab apply the remaining peptide steps and DPC Quantification (normally followed by ",
        "Between-Sample Normalisation and the Peptide Count Filter), then press \"Run Normalization & RUV\" again."),
        class = "prot_norm_limpa_peptide_input")
    head <- identical(source$entity_id, getDIAStateReference(manager)$entity_id)
    if (identical(spec$class, "ProteinQuantitativeData") && diaNormalisedLineage(manager, spec, head = head))
        abortDIANormalisedInput(manager)
    invisible(TRUE)
}

#' The normalisation's source and declared parameters, as the tab sends them in process and to the worker (KIT-111 D3):
#' the source by `diaNormalizationSource` (KIT-182), recorded as the tab's normalisation input. KIT-182: after a completed
#' run, a protein source other than the previous run's signals the classed message `prot_norm_input_current`, which the
#' tab shows as a notice; a source other than the table the pre-normalisation figures were drawn from becomes the
#' figures' source (`qc_sources$pre_normalization`) and `redraw_pre` is TRUE.
#' @param manager,norm_data Session owner and the tab's state.
#' @param settings The normalisation settings (`runDIAReferenceNormalization`).
#' @param grouping The RUV grouping variable.
#' @return `list(source, parameters, redraw_pre)`.
#' @noRd
diaNormalizationRequest <- function(manager, settings, norm_data, grouping) {
    previous <- norm_data$normalization_input_reference
    source <- diaNormalizationSource(manager, norm_data)
    norm_data$normalization_input_reference <- source
    redraw_pre <- FALSE
    if (usesDIAReferenceStates(manager) && isSavedWorkflowReference(source)) {
        if (isTRUE(norm_data$normalization_complete) && isSavedWorkflowReference(previous) &&
            !identical(previous$entity_id, source$entity_id) &&
            identical(readDIALineageRecord(manager, source)$class, "ProteinQuantitativeData")) {
            message(structure(class = c("prot_norm_input_current", "message", "condition"), list(message = paste(
                "This normalisation uses your current quality-control result. The previous normalisation used an",
                "earlier table, from before your latest quality-control changes."), call = NULL)))
        }
        snapshot <- manager$.context$snapshot
        drawn <- snapshot$qc_sources$pre_normalization
        if (!is.null(drawn) && !identical(drawn$entity_id, source$entity_id)) {
            snapshot$qc_sources$pre_normalization <- source
            publishDIAStateHead(manager$.context, snapshot)
            redraw_pre <- TRUE
        }
    }
    fields <- c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
        "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k")
    parameters <- stats::setNames(lapply(fields, \(field) settings[[field]]), fields)
    parameters$grouping <- grouping
    # KIT-129 integration item 1: with the gate on, a new design the tab names travels with its k rule and seed.
    if (omicGuiRuvDesignEnabled() && is.character(settings$ruv_design) && !identical(settings$ruv_design, "legacy")) {
        parameters[c("ruv_design", PROT_NORM_RUV_K_PARAMETERS)] <- list(settings$ruv_design, settings$ruv_k_rule, settings$ruv_k_seed)
        # KIT-129 fifth pass item 3: the supplied list beside them (its identity is set by `workflowRequestOperation`).
        if (identical(settings$ruv_design, "supplied_controls")) parameters$ruv_control_ids <- settings$ruv_control_ids
    }
    list(source = source, parameters = parameters, redraw_pre = redraw_pre)
}

#' The normalisation's compute half (KIT-111 D3): the declaration on the source, reused on a hit. Worker-safe: no head
#' write, no Shiny, no tab state.
#' @return `runSavedOperation`'s run.
#' @noRd
computeDIAReferenceNormalization <- function(manager, source, parameters) {
    # KIT-129 integration item 1: the new design's declaration when the request names it (`ruv_design`), on both paths.
    chosen <- workflowRequestOperation(manager, workflowOperationEntry("normalization-ruv"), source, parameters)
    runSavedOperation(manager, chosen$declaration, source, chosen$parameters)
}

#' The normalisation's commit description (KIT-111 D3), shared by the in-process runner and the worker's completion.
#' @param settings The settings or the request's parameters (`norm_method`, `ruv_mode`).
#' @param run The compute half's run.
#' @return `list(state_name, checkpoint, config, description)`, or NULL when the head already holds the state with
#'   this product.
#' @noRd
diaNormalizationCommit <- function(manager, settings, run) {
    saved <- run$value
    skipped <- identical(settings$ruv_mode, "skip")
    final <- readDIAStateMetadata(manager, saved$final)$data
    config <- list(norm_method = settings$norm_method, ruv_mode = settings$ruv_mode,
        ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage)
    if (skipped) config$ruv_applied <- FALSE
    # KIT-129 seventh pass item 2: a new design whose k rule gave no factor (a count of 0) keeps the state's name and its
    # normalised matrix, and records that no RUV was applied (the integration input then says not corrected).
    if (!is.null(saved$parameters$ruv_design) && identical(saved$ruv_result$ruv_applied, FALSE)) config$ruv_applied <- FALSE
    state_name <- if (skipped) "normalized" else "ruv_corrected"
    final$products$normalization_science <- run$reference
    previous <- if (manager$current_state == state_name) readDIAStateMetadata(manager) else NULL
    if (!is.null(previous) && identical(previous$data$products$normalization_science, run$reference)) return(NULL)
    list(state_name = state_name, checkpoint = final, config = config,
        description = if (skipped) "Post-normalization complete: RUV-III skipped by user" else
            "Post-normalization complete: RUV-III correction and missing value cleanup completed")
}

#' The normalisation's publication half (KIT-111 D3), in the controller on both paths: the tab's state, the RUV
#' optimisation table and the step's files and display, on the accepted head.
#' @param run The compute half's run (`value`, `reference`); on the worker path the delivered product.
#' @inheritParams runDIAReferenceNormalization
#' @return TRUE after the required tables and plots have completed.
#' @noRd
publishDIAReferenceNormalization <- function(manager, settings, run, workflow, norm_data, paths, omic_type,
    experiment_label, post_qc, ruv_qc) {
    saved <- run$value
    product <- run$reference
    norm_data$normalized_protein_obj <- saved$normalized
    norm_data$ruv_normalized_obj <- saved$final
    norm_data$best_k <- saved$best_k
    norm_data$control_genes_index <- saved$controls
    norm_data$ruv_optimization_result <- saved$ruv_result
    workflow$ruv_optimization_result <- saved$ruv_result
    # KIT-129 fourth pass item 3: what a new design's run did, for the RUV design section (nothing with the gate off).
    omicNormSetRuvDesignResult(norm_data, protNormRuvDesignResult(saved))
    # KIT-018: the optimisation grid is a saved display table the module pages on demand.
    norm_data$ruv_optimization_table <- saveDIADisplayTable(manager, saved$ruv_result$optimization_results)
    skipped <- identical(settings$ruv_mode, "skip")
    final <- readDIAStateMetadata(manager, saved$final)$data
    accepted <- getDIAStateReference(manager)
    # KIT-052: drawn in a short folder and copied into the action staging folder (`withWorkflowShortStaging`); drawn
    # there, a QC figure, the cancor figure and the composite were lost without a word from 260 characters on Windows.
    display <- publishDIAActionOutputs(manager, accepted, "normalization",
        paths[c("source_dir", "protein_qc_dir", "time_dir")], function(staging) withWorkflowShortStaging(staging, function(staged) {
            norm_data$publication_paths <- staged
            on.exit({ norm_data$publication_paths <- NULL }, add = TRUE)
            # KIT-259 third pass (VAL-307, C01): with their placeholders' markers, which the combined figure's list reads.
            pre_files <- list.files(paths$protein_qc_dir, pattern = "^pre_norm_.*[.]png([.]placeholder)?$", full.names = TRUE)
            if (length(pre_files) && !all(file.copy(pre_files, staged$protein_qc_dir))) {
                cli::cli_abort("cannot stage pre-normalization display inputs")
            }
            norm_spec <- readDIAStateMetadata(manager, saved$normalized)$data
            writeDIADesignTable(manager, norm_spec$components$protein_quant_table,
                file.path(staged$protein_qc_dir, "normalized_protein_matrix_pre_ruv.tsv"), vroom::vroom_write)
            saveRDS(saved$ruv_result, file.path(staged$source_dir, "ruv_optimization_results.RDS"))
            post_qc(saved$normalized)
            if (!skipped) {
                ruv_qc(saved$final)
                summary <- readDIAProteinProgress(manager, saved$final)
                directories <- stats::setNames(list(list(time_dir = staged$time_dir,
                    publication_graphs_dir = paths$publication_graphs_dir)),
                    paste(omic_type, experiment_label, sep = "_"))
                norm_data$post_norm_filtering_plot <- updateProteinFiltering(summary,
                    "11_RUV_filtered", omic_type, experiment_label, overwrite = TRUE,
                    return_grid = TRUE, project_dirs = directories)
                workflow$protein_counts$after_ruv_filtering <- final$metadata$proteins
                norm_data$filtering_summary_text <- sprintf(
                    "Filtering Summary through Normalization & RUV:\n\n* Pre-normalization proteins: [Previous step]\n* Post-RUV filtering: %d proteins\n* Proteins removed by RUV: [Calculated from difference]\n\nRUV Parameters:\n* Normalization method: %s\n* RUV mode: %s\n* RUV k value: %d\n* Negative control %%: %.1f",
                    final$metadata$proteins, settings$norm_method, settings$ruv_mode, saved$best_k, saved$percentage)
                # KIT-259 (B2): a canonical-correlation figure that cannot be drawn is replaced by its placeholder.
                if (!is.null(saved$ruv_result$best_cancor_plot)) drawProtNormQcFigure(staged$protein_qc_dir,
                    "ruv_corrected_cancor.png", 8, 6, function() ggplot2::ggsave(
                    file.path(staged$protein_qc_dir, "ruv_corrected_cancor.png"),
                    saved$ruv_result$best_cancor_plot, width = 8, height = 6, dpi = 150))
            }
            # KIT-122: on a managed session the worker built the composite from these staged images and saved its files
            # (`saveWorkflowDiaCompositeFigure`); their bytes are written here. Not found saved (no runner, or other
            # images than the worker read), the figure is built here, as before.
            if (!writeDIACompositeFigureFiles(manager, product, settings$ruv_mode, staged$protein_qc_dir, omic_type))
                generateProtNormCompositeQcFigure(settings$ruv_mode, staged$protein_qc_dir, omic_type)
            # KIT-259 third pass (VAL-307, C01): the combined figure's replaced panels, written with it.
            writeProtNormCompositePlaceholders(settings$ruv_mode, staged$protein_qc_dir, omic_type)
            list(qc_paths = norm_data$qc_plot_paths, filtering_plot = norm_data$post_norm_filtering_plot,
                filtering_text = norm_data$filtering_summary_text, staged_qc_dir = staged$protein_qc_dir)
        }, "the normalisation's files could not be copied into the project's staging folder"), identity = list(product = product$entity_id, color = settings$color_variable,
            shape = settings$shape_variable, output_schema = 2L))
    norm_data$qc_plot_paths <- lapply(display$qc_paths, function(group) lapply(group, function(path) {
        if (is.character(path)) sub(display$staged_qc_dir, paths$protein_qc_dir, path, fixed = TRUE) else path
    }))
    norm_data$post_norm_filtering_plot <- display$filtering_plot
    norm_data$filtering_summary_text <- display$filtering_text
    norm_data$normalization_complete <- TRUE
    finalizeProtNormWorkflowState(norm_data)
    TRUE
}

#' Count the accepted protein display population with projected columns.
#' @param manager,reference Current owner and accepted checkpoint.
#' @return Bounded FilteringProgress input, with no retained assay.
#' @noRd
readDIAProteinProgress <- function(manager, reference, summary = NULL) {
    spec <- readDIAStateMetadata(manager, reference)$data
    context <- manager$.context
    scope <- newWorkflowTableReadScope(context$store, spec$components$protein_quant_table,
        reader = "arrow", metadata_only = TRUE)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    identity <- newSavedWorkflowScientificIdentity(workflowTablePopulationId(scope),
        spec$components$design_matrix$entity_id, digest::digest(spec$slots[c("sample_id", "protein_id_column")],
            algo = "sha256"), "protein-filtering-counts-v1", context$software)
    hit <- resolveSavedWorkflowEntity(context$store, "scientific_product", identity)
    if (identical(hit$status, "hit")) return(hit$value$value$summary)
    if (is.null(summary)) {
        design <- readDIAExecutionComponents(manager, reference, fields = "design_matrix")@design_matrix
        columns <- as.character(design[[spec$slots$sample_id]])
        if (length(columns) > 10000L || anyNA(columns) || !all(columns %in% scope$metadata$columns))
            cli::cli_abort("declared progress samples are invalid", class = "saved_entity_invalid")
        releaseWorkflowTableReadScope(scope)
        scope <- newWorkflowTableReadScope(context$store, spec$components$protein_quant_table)
        assertWorkflowPhysicalTableScope(scope)
        unique_columns <- unique(columns)
        counts <- stats::setNames(integer(length(unique_columns)), unique_columns)
        for (block in split(unique_columns, ceiling(seq_along(unique_columns) / 64L))) {
            quoted <- as.character(DBI::dbQuoteIdentifier(scope$connection, block))
            expressions <- paste0("count(*) FILTER (WHERE ", quoted,
                " IS NOT NULL AND NOT isnan(", quoted, ")) AS ", quoted)
            value <- DBI::dbGetQuery(scope$connection, paste0("SELECT ",
                paste(expressions, collapse = ", "), " FROM read_parquet(?)"), params = list(scope$path))
            counts[block] <- as.integer(value[1L, ])
        }
        per_run <- data.frame(Run = columns, n_proteins = counts[columns])
        summary <- structure(list(is_protein_quant = TRUE, proteins = spec$dimensions$rows,
            proteins_per_run = dplyr::arrange(per_run, .data$Run), total_peptides = NA_integer_,
            peptides_per_protein = data.frame(), peptides_per_run = data.frame()), class = "DIAFilteringSummary")
    }
    value <- list(summary = summary, dependencies = list(scope$physical_reference,
        spec$components$design_matrix))
    product <- saveWorkflowConsumerEntity(context$store,
        newSavedWorkflowScientificProduct(value, "dia-protein-filtering-counts"),
        basename(tempfile("dia-protein-counts-")), identity, "scientific_product")
    snapshot <- context$snapshot
    snapshot$output_products[[paste0("counts-", reference$entity_id)]] <- product
    publishDIAStateHead(context, snapshot)
    summary
}

#' The two correlation operations (KIT-004): within-group Pearson pairs (a product) and the sample
#' selection at a threshold (a checkpointed object whose reference is the product).
#' @noRd
diaPearsonPairsOperation <- function() {
    newWorkflowOperation("pearson-pairs", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        product_type = "dia-pearson-pairs", parameters = c("group", "tech_rep_remove_regex"),
        kernel = function(object, parameters, store) {
            pearsonCorForSamplePairs(object, tech_rep_remove_regex = parameters$tech_rep_remove_regex,
                correlation_group = parameters$group)
        })
}

diaPearsonSelectionOperation <- function() {
    newWorkflowOperation("pearson-selection", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        output_class = "ProteinQuantitativeData", product_type = "dia-pearson-selection",
        parameters = c("pairs", "threshold"),
        kernel = function(object, parameters, store) {
            filterSamplesFn <- filterSamplesByProteinCorrelationThreshold
            filtered <- filterSamplesFn(object, pearson_correlation_per_pair = parameters$correlations,
                min_pearson_correlation_threshold = parameters$threshold)
            # KIT-041: on the DIA_LIMPA route the Limpa results travel with the object; a dropped sample leaves the
            # quantified EList (and its standard errors) as well as the matrix, so the EList the DA branch fits
            # matches the design. A DIA route object carries no such entry and is untouched.
            filtered <- subsetLimpaResultsToSamples(filtered)
            store$checkpoint(filtered, "pearson-filtered-science")
        })
}

#' The correlation filter's Skip as a step (KIT-248, VAL-238): no parameters, the fixed constant `skipped = TRUE`, no
#' kernel work. The Skip commits the normalised source as `correlation_filtered` in the controller
#' (`runDIAReferenceCorrelation`); the declaration names that step for the analysis record, the study parameters and
#' the replay.
#' @noRd
diaPearsonSkipOperation <- function() {
    newWorkflowOperation("pearson-skip", "ProteinQuantitativeData", fields = "design_matrix", unit = "store",
        parameters = character(), product_type = "dia-pearson-skip", fixed = list(skipped = TRUE),
        kernel = function(object, parameters, store) list())
}

#' Keep the Limpa results' sample columns in step with the design (KIT-041): the quantified EList, its standard
#' errors and observation counts are subset to the design's samples when a sample was dropped.
#' @noRd
subsetLimpaResultsToSamples <- function(object) {
    results <- object@args$limpa_dpc_quant_results
    if (!is.list(results) || is.null(results$quantified_elist)) return(object)
    samples <- as.character(object@design_matrix[[object@sample_id]])
    present <- colnames(results$quantified_elist$E)
    if (is.null(present) || setequal(present, samples)) return(object)
    keep <- intersect(present, samples)
    results$quantified_elist <- results$quantified_elist[, keep]
    for (field in c("standard_errors", "n_observations", "y_peptide_for_dpc")) {
        value <- results[[field]]
        if (is.matrix(value) && !is.null(colnames(value)) && all(keep %in% colnames(value)))
            results[[field]] <- value[, keep, drop = FALSE]
    }
    object@args$limpa_dpc_quant_results <- results
    object
}

#' RUV-III factors as covariates (KIT-041, the limpa script's chunk 19b) as a whole-matrix declaration run in
#' process: the script's negative-control sweep (`findBestNegCtrlPercentageCapped`), then `computeRuvCovariates`
#' (RUV-III with the design's replicate structure, the W factors); the output state carries the quantified matrix
#' unchanged and the design with `W1..Wk` appended; the optimisation table, the cancor plot and the RUV-adjusted
#' matrix (a component checkpoint, never a state) are the product.
#' @param ruv_design KIT-129 W3. `"legacy"` (the default): today's declaration and kernel path, the sweep gives the
#'   controls and k; its ten parameters and identities are unchanged. `"declared"`: the declaration gains the
#'   parameters of `diaRuvCovariatesRuleParameters` and the kernel is `diaRuvCovariatesByRuleKernel` (the control
#'   rule and the k rule as declared; the sweep runs only when the k rule is today's estimate).
#' @noRd
diaRuvCovariatesOperation <- function(ruv_design = c("legacy", "declared")) {
    declared_rules <- identical(match.arg(ruv_design), "declared")
    newWorkflowOperation("protein_ruv_covariates", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        output_class = "ProteinQuantitativeData", product_type = "dia-ruv-covariates",
        parameters = c("percentage_min", "percentage_max", "num_components_to_impute", "ruv_grouping_variable",
            "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_mode",
            if (declared_rules) diaRuvCovariatesRuleParameters),
        kernel = function(object, parameters, store) {
            if (declared_rules) return(diaRuvCovariatesByRuleKernel(object, parameters, store))
            assertDiaRuvSweepParameters(parameters)
            findBestNegCtrlPercentageCappedFn <- findBestNegCtrlPercentageCapped
            sweep <- suppressMessages(findBestNegCtrlPercentageCappedFn(object,
                percentage_range = seq(parameters$percentage_min, parameters$percentage_max, by = 1),
                num_components_to_impute = parameters$num_components_to_impute,
                ruv_grouping_variable = parameters$ruv_grouping_variable, ruv_qval_cutoff = parameters$ruv_qval_cutoff,
                ruv_fdr_method = parameters$ruv_fdr_method, separation_metric = parameters$separation_metric,
                k_penalty_weight = parameters$k_penalty_weight, adaptive_k_penalty = isTRUE(parameters$adaptive_k_penalty),
                verbose = FALSE))
            computeRuvCovariatesFn <- computeRuvCovariates
            factors <- computeRuvCovariatesFn(object, sweep$best_k, sweep$best_control_genes_index, parameters$ruv_grouping_variable)
            with_w <- object
            samples <- as.character(with_w@design_matrix[[with_w@sample_id]])
            for (i in seq_len(factors$k)) with_w@design_matrix[[paste0("W", i)]] <- unname(factors$W[samples, i])
            control_count <- sum(sweep$best_control_genes_index, na.rm = TRUE)
            with_w@args$ruvCovariates <- list(k = factors$k, percentage_as_neg_ctrl = sweep$best_percentage,
                control_count = control_count, replicate_structure = factors$replicate_structure,
                block_reason = factors$block_reason, blocks = factors$blocks,
                ruv_grouping_variable = parameters$ruv_grouping_variable, separation_score = sweep$best_separation_score,
                composite_score = sweep$best_composite_score, max_acceptable_k = sweep$max_acceptable_k,
                ruv_mode = parameters$ruv_mode %||% "covariates")
            covariates_reference <- store$checkpoint(with_w, "ruv-covariates-science")
            adjusted <- object
            adjusted_table <- as.data.frame(factors$adjusted)
            adjusted_table[[adjusted@protein_id_column]] <- rownames(factors$adjusted)
            rownames(adjusted_table) <- NULL
            adjusted@protein_quant_table <- adjusted_table[, colnames(object@protein_quant_table), drop = FALSE]
            adjusted_reference <- store$checkpoint(adjusted, "ruv-adjusted-science", from = covariates_reference)
            final_reference <- covariates_reference
            if (identical(parameters$ruv_mode, "correction")) {
                # KIT-056 (user decision 2026-10-22): the DIA route's form on this route. The RUV-III adjusted matrix
                # (Y - W alpha, RUV-III's newY, the same fit, k and controls as the covariates mode) becomes the state's
                # protein table and the owned DPC EList's E; the DPC standard errors and observation counts are kept as
                # quantified (a linear shift of E leaves both valid); the design carries no W terms; dpcDE then fits the
                # plain formula on the corrected EList.
                corrected <- adjusted
                results <- corrected@args$limpa_dpc_quant_results
                if (!is.list(results) || is.null(results$quantified_elist))
                    cli::cli_abort("the RUV-III correction needs the DPC quantification's owned EList", class = "workflow_execution_invalid")
                E <- results$quantified_elist$E
                if (!all(rownames(E) %in% rownames(factors$adjusted)) || !all(colnames(E) %in% colnames(factors$adjusted)))
                    cli::cli_abort("the RUV-III adjusted matrix does not cover the quantified EList", class = "workflow_execution_invalid")
                results$quantified_elist$E <- factors$adjusted[rownames(E), colnames(E), drop = FALSE]
                corrected@args$limpa_dpc_quant_results <- results
                corrected@args$ruvCovariates <- with_w@args$ruvCovariates
                final_reference <- store$checkpoint(corrected, "ruv-corrected-science", from = covariates_reference)
                corrected <- results <- E <- NULL
            }
            cancor_plot <- if (!is.null(sweep$best_cancor_plot)) ggplot2::ggplotGrob(stripProtNormQcPlotOwners(sweep$best_cancor_plot)) else NULL
            ruv_result <- list(best_percentage = sweep$best_percentage, best_k = factors$k,
                best_control_genes_index = sweep$best_control_genes_index, best_separation_score = sweep$best_separation_score,
                best_composite_score = sweep$best_composite_score, optimization_results = sweep$optimization_results,
                best_cancor_plot = cancor_plot, max_acceptable_k = sweep$max_acceptable_k, sample_size = sweep$sample_size,
                separation_metric_used = parameters$separation_metric, k_penalty_weight = parameters$k_penalty_weight,
                adaptive_k_penalty_used = isTRUE(parameters$adaptive_k_penalty),
                ruv_mode = parameters$ruv_mode %||% "covariates", replicate_structure = factors$replicate_structure, block_reason = factors$block_reason)
            object <- with_w <- adjusted <- NULL
            list(covariates = covariates_reference, adjusted = adjusted_reference, final = final_reference,
                ruv_mode = parameters$ruv_mode %||% "covariates", ruv_result = ruv_result,
                W = factors$W, best_k = factors$k, percentage = sweep$best_percentage, controls = sweep$best_control_genes_index,
                control_count = control_count, replicate_structure = factors$replicate_structure,
                block_reason = factors$block_reason, blocks = factors$blocks,
                parameters = parameters[c("percentage_min", "percentage_max", "num_components_to_impute", "ruv_grouping_variable",
                    "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_mode")])
        })
}

#' The `[ruvCovariates]` parameters of a project (KIT-041): the project's configuration section when present, else
#' the shipped defaults; `ruv_grouping_variable` from the tab's selection when given; `ruv_mode` (KIT-056) the tab's
#' choice, `covariates` (the script's method, the default) or `correction` (the DIA route's form).
#' @param ruv_control_rule,ruv_k_rule,ruv_k,ruv_control_ids,ruv_k_seed KIT-129 W3. With the defaults the value is
#'   today's ten parameters, unchanged. A control rule other than today's sweep appends the five rule parameters
#'   (`resolveDiaRuvRuleParameters`), for `diaRuvCovariatesOperation("declared")`.
#' @noRd
resolveDIARuvCovariatesParameters <- function(config_list, grouping = NULL, ruv_mode = "covariates",
    ruv_control_rule = "sweep_least_significant", ruv_k_rule = NULL, ruv_k = NULL, ruv_control_ids = NULL, ruv_k_seed = NULL) {
    if (!is.character(ruv_mode) || length(ruv_mode) != 1L || !ruv_mode %in% c("covariates", "correction"))
        cli::cli_abort("ruv_mode must be covariates or correction", class = "workflow_execution_invalid")
    section <- config_list$ruvCovariates %||% list()
    shipped <- readWorkflowShippedDefaults()$ruvCovariates %||% list()
    pick <- function(name, default) {
        value <- section[[name]] %||% shipped[[name]] %||% default
        value
    }
    numeric <- function(name, default) as.numeric(pick(name, default))
    parameters <- list(percentage_min = numeric("percentage_min", 1), percentage_max = numeric("percentage_max", 30),
        num_components_to_impute = numeric("num_components_to_impute", 5),
        ruv_grouping_variable = grouping %||% as.character(pick("ruv_grouping_variable", "group")),
        ruv_qval_cutoff = numeric("ruv_qval_cutoff", 0.05), ruv_fdr_method = as.character(pick("ruv_fdr_method", "BH")),
        separation_metric = as.character(pick("separation_metric", "max_difference")),
        k_penalty_weight = numeric("k_penalty_weight", 0.5),
        adaptive_k_penalty = tolower(as.character(pick("adaptive_k_penalty", "true"))) %in% c("true", "1", "yes"),
        ruv_mode = ruv_mode)
    todays_sweep <- identical(ruv_control_rule, "sweep_least_significant") && is.null(ruv_k_rule) && is.null(ruv_k) &&
        is.null(ruv_control_ids) && is.null(ruv_k_seed)
    if (todays_sweep) return(parameters)
    c(parameters, resolveDiaRuvRuleParameters(ruv_mode, ruv_control_rule, ruv_k_rule, ruv_k, ruv_control_ids, ruv_k_seed))
}

#' The declared parameters of the RUV covariates' control rule and k rule (KIT-129 W3), after the sweep's ten.
#' @noRd
diaRuvCovariatesRuleParameters <- c("ruv_control_rule", "ruv_k_rule", "ruv_k", "ruv_control_ids", "ruv_k_seed")

#' Today's declaration takes none of the rule parameters (KIT-129 audit F3): a request that carries one is refused,
#' never run as the sweep. Its identity reads the sweep's ten parameters only, so `runSavedOperation` returns a saved
#' sweep product for such a request before any kernel runs; the route's compute half
#' (`computeDIAReferenceRuvCovariates`) therefore checks before the runner, and the kernel checks on a computed run.
#' @param parameters A request's parameters.
#' @return `parameters`, invisibly; an error of class `workflow_execution_invalid` (the legacy refusal of every omic,
#'   KIT-129 item 1) when a rule parameter is present with a value.
#' @noRd
assertDiaRuvSweepParameters <- function(parameters) {
    given <- Filter(\(name) !is.null(parameters[[name]]), diaRuvCovariatesRuleParameters)
    if (length(given)) {
        cli::cli_abort(c("today's RUV covariates declaration (the sweep) takes no rule parameter; the request carries {.val {given}}",
            i = "a control rule or a k rule runs under {.code diaRuvCovariatesOperation(\"declared\")}"),
            class = "workflow_execution_invalid", call = NULL)
    }
    invisible(parameters)
}

#' The control rule and the k rule of the RUV covariates as declared parameters (KIT-129 W3). The control rules are
#' the designs of KIT-129 section 1b that leave today's sweep: `all_features` (arms 7 and 16 of simulator S2),
#' `supplied_list` (arm 8) and `first_pass_drop_top_25` (arm 5). The k rule, when none is given, is the one the
#' design was measured with: the estimated count, or today's estimate for the first-pass nulls. Offered in the
#' covariates mode only.
#' @param ruv_mode,ruv_control_rule,ruv_k_rule,ruv_k,ruv_control_ids,ruv_k_seed The tab's choices; `ruv_k` for the
#'   k rule `user` only, `ruv_control_ids` for `supplied_list` only, `ruv_k_seed` the seed of the estimated count.
#' @return The five parameters of `diaRuvCovariatesRuleParameters`.
#' @noRd
resolveDiaRuvRuleParameters <- function(ruv_mode, ruv_control_rule, ruv_k_rule = NULL, ruv_k = NULL,
                                        ruv_control_ids = NULL, ruv_k_seed = NULL) {
    invalid <- function(message) cli::cli_abort(message, class = "workflow_execution_invalid", call = NULL)
    declared <- setdiff(RUV_CONTROL_RULES, "sweep_least_significant")
    if (!checkmate::test_choice(ruv_control_rule, declared)) {
        invalid(paste0("a k rule, a k, a control list or a seed needs ruv_control_rule to be one of ",
            paste(declared, collapse = ", "), "; today's sweep takes none of them"))
    }
    if (!identical(ruv_mode, "covariates")) invalid("a declared control rule is offered in the covariates mode only")
    ruv_k_rule <- ruv_k_rule %||% if (identical(ruv_control_rule, "first_pass_drop_top_25")) "sweep_bounded" else "num_sv_be"
    if (!checkmate::test_choice(ruv_k_rule, RUV_K_RULES)) {
        invalid(paste("ruv_k_rule must be one of", paste(RUV_K_RULES, collapse = ", ")))
    }
    if (identical(ruv_k_rule, "user") != checkmate::test_int(ruv_k, lower = 1)) {
        invalid("ruv_k is a whole number of 1 or more for the k rule user, and is given for no other rule")
    }
    if (identical(ruv_control_rule, "supplied_list") != checkmate::test_character(ruv_control_ids, min.len = 1L, any.missing = FALSE)) {
        invalid("ruv_control_ids names the control proteins for the control rule supplied_list, and is given for no other rule")
    }
    ruv_k_seed <- ruv_k_seed %||% RUV_FACTOR_COUNT_SEED
    if (!checkmate::test_int(ruv_k_seed)) invalid("ruv_k_seed must be a whole number")
    list(ruv_control_rule = ruv_control_rule, ruv_k_rule = ruv_k_rule, ruv_k = if (!is.null(ruv_k)) as.integer(ruv_k),
        ruv_control_ids = if (!is.null(ruv_control_ids)) sort(unique(ruv_control_ids)), ruv_k_seed = as.integer(ruv_k_seed))
}

#' The first pass of the first-pass-nulls control rule on the DIA_LIMPA route (KIT-129 W3; arm 5 of simulator S2):
#' the route's own test without RUV, `limpa::dpcDE` on the quantified EList with the groups as the model and the
#' design's technical-replicate block, then the moderated F test of the differences between the groups (trend and
#' robust, the shipped settings). The contrasts of the DA are not known at the RUV step; with two groups the F test
#' is the test of their difference, which is what the arm ranked.
#' @param object The quantified `ProteinQuantitativeData` object, carrying the DPC quantification's EList.
#' @param grouping The RUV grouping variable.
#' @return The p-values, named by protein.
#' @noRd
diaRuvFirstPassPValues <- function(object, grouping) {
    results <- object@args$limpa_dpc_quant_results
    if (!is.list(results) || is.null(results$quantified_elist) || !requireNamespace("limpa", quietly = TRUE)) {
        cli::cli_abort("the first-pass test needs limpa and the DPC quantification's owned EList", class = "workflow_execution_invalid")
    }
    design <- as.data.frame(object@design_matrix)
    rownames(design) <- as.character(design[[object@sample_id]])
    model <- stats::model.matrix(stats::as.formula(paste("~ 0 +", grouping)), design)
    samples <- rownames(model)
    elist <- results$quantified_elist
    if (ncol(model) < 2L || !all(samples %in% colnames(elist$E))) {
        cli::cli_abort("the first-pass test needs two groups or more and every sample of the design in the quantified EList",
            class = "workflow_execution_invalid")
    }
    replicate <- resolveProtTechnicalReplicateBlock(design[samples, , drop = FALSE], object@group_id, object@technical_replicate_id)
    fit <- limpa::dpcDE(elist[, samples], model, plot = FALSE, block = if (isTRUE(replicate$technical)) replicate$block else NULL)
    differences <- rbind(-1, diag(ncol(model) - 1L))
    rownames(differences) <- colnames(model)
    tested <- limma::eBayes(limma::contrasts.fit(fit, differences), trend = TRUE, robust = TRUE)
    stats::setNames(tested$F.p.value, rownames(tested$coefficients))
}

#' The separation rule for k on a supplied control list on the DIA with limpa route (KIT-129 sixth pass item 3), as
#' simulator S2's `separationK` computes it for arm 22: the canonical-correlation curve of the controls against all
#' proteins (`prepareRuvCancorInput` with the declaration's components, `ruvCancorFromInput` with the grouping), up to
#' `min(calculateAdaptiveMaxK(samples), controls - 1)`, k at the first largest all-minus-control difference.
#' @param object The quantified `ProteinQuantitativeData` object.
#' @param controls A logical vector over the protein rows.
#' @param grouping The grouping column. @param components `num_components_to_impute`.
#' @return `list(k, plot)`; refuses with `ruv_k_not_estimated` (and `ruv_k_refused`) when the curve gives no k.
#' @noRd
diaRuvSeparationRuleK <- function(object, controls, grouping, components) {
    plot <- ruvCancorFromInput(prepareRuvCancorInput(object, components), object, unname(controls), grouping)
    curve <- plot$data[plot$data$K <= min(calculateAdaptiveMaxK(nrow(object@design_matrix)), sum(controls) - 1), , drop = FALSE]
    of_controls <- curve$featureset == "Control"
    difference <- curve$cc[curve$featureset == "All"] - curve$cc[of_controls]
    best <- which(difference == max(difference, na.rm = TRUE))
    if (!length(best)) {
        cli::cli_abort("the number of unwanted factors (k) is refused: the separation rule returned no k for the supplied list",
            class = c("ruv_k_not_estimated", "ruv_k_refused"), call = NULL)
    }
    list(k = as.integer(curve$K[of_controls][[best[[1L]]]]), plot = plot)
}

#' The controls, k and factors of the RUV covariates by declared rules (KIT-129 W3), with no store: the science of
#' `diaRuvCovariatesByRuleKernel`. The controls follow the control rule; k comes from `resolveRuvFactorCount` (the
#' estimated count on the quantified matrix, today's sweep estimate, or the user's number) under the limit the
#' design carries; the factors are `computeRuvCovariates` with the design's replicate structure, as today (the
#' declared technical replicates where the design has them, else the groups).
#'
#' A user's k above the limit, the separation rule with all proteins as controls and a missing `sva` are errors
#' for the user to act on. Under an automatic rule, what the data cannot give (too few complete proteins or
#' samples for the estimate, a k above the limit) leaves the run without RUV, with a warning of class
#' `dia_ruv_not_applied`; an estimate of 0 does the same with the resolver's message (`ruv_k_zero`).
#' @param object The quantified `ProteinQuantitativeData` object.
#' @param parameters The declaration's parameters (`resolveDIARuvCovariatesParameters` with a control rule).
#'
#' Supplied ids that are not in the protein table raise a warning of class `dia_ruv_control_ids_unmatched` with
#' the fields `matched`, `listed` and `unmatched` (audit F6); a list that matches none is refused, as before.
#' @return `list(status, reason, k, factors, controls, controls_unmatched, sweep, replicate_structure, block_reason,
#'   blocks)`: `status` is `ok`, `none` (k = 0) or `refused`; `k` the resolver's record; `factors` the value of
#'   `computeRuvCovariates` or NULL; `controls` a logical vector named by protein; `sweep` today's sweep when it ran.
#' @noRd
computeDiaRuvCovariatesByRule <- function(object, parameters) {
    invalid <- function(message) cli::cli_abort(message, class = "workflow_execution_invalid", call = NULL)
    control_rule <- parameters$ruv_control_rule
    k_rule <- parameters$ruv_k_rule
    grouping <- parameters$ruv_grouping_variable
    design <- object@design_matrix
    if (!checkmate::test_choice(control_rule, setdiff(RUV_CONTROL_RULES, "sweep_least_significant")) ||
        !checkmate::test_choice(k_rule, RUV_K_RULES) || !identical(parameters$ruv_mode, "covariates")) {
        invalid("the declared RUV covariates need a control rule, a k rule and the covariates mode")
    }
    if (!checkmate::test_choice(grouping, names(design))) invalid("the RUV grouping variable is not a column of the design")
    ids <- as.character(object@protein_quant_table[[object@protein_id_column]])
    samples <- as.character(design[[object@sample_id]])

    # --- Controls ---
    unmatched <- character()
    controls <- if (identical(control_rule, "all_features")) {
        rep(TRUE, length(ids))
    } else if (identical(control_rule, "supplied_list")) {
        unmatched <- setdiff(parameters$ruv_control_ids, ids)
        matched <- ids %in% parameters$ruv_control_ids
        if (!any(matched)) invalid("none of the supplied control ids is in the protein table")
        if (length(unmatched)) {
            listed <- length(unique(parameters$ruv_control_ids))
            cli::cli_warn(c("{sum(matched)} of {listed} supplied control ids are in the protein table; the controls are those {sum(matched)}",
                i = "not in the table: {.val {utils::head(unmatched, 5L)}}{if (length(unmatched) > 5L) ' and others'}"),
                class = "dia_ruv_control_ids_unmatched", matched = sum(matched), listed = listed, unmatched = unmatched)
        }
        matched
    } else {
        first_pass <- diaRuvFirstPassPValues(object, grouping)[ids]
        if (anyNA(first_pass)) invalid("the first-pass test gave no p-value for some proteins of the table")
        # arm 5 of simulator S2: every protein outside the 25% most significant of the first pass
        unname(rank(first_pass, ties.method = "first") > 0.25 * length(ids))
    }

    # --- k ---
    replicates <- resolveProtRuvReplicateBlocks(design, object@sample_id, object@group_id, object@technical_replicate_id, grouping)
    block <- resolveProtTechnicalReplicateBlock(design, object@group_id, object@technical_replicate_id)
    # Decision M28: on a design with declared technical replicates the biological units carry the test, so k is limited
    # by them and the count removes them.
    units <- if (isTRUE(block$technical)) block$block
    limit <- ruvFactorCountLimit(length(samples), replicates$blocks, length(unique(as.character(design[[grouping]]))), sum(controls),
        units = if (!is.null(units)) length(unique(units)))
    # Today's sweep runs only where the resolver reads it: the separation rule, which it refuses with all proteins
    # as controls before asking for a k. KIT-129 sixth pass item 3: with a supplied list the separation rule reads that
    # list (simulator S2's arm 22, `diaRuvSeparationRuleK`), as in the other omics; the sweep is not run.
    # A refusal of the rule is raised where the resolver reads it, so that an automatic rule leaves the run without RUV.
    separation <- if (identical(k_rule, "sweep_bounded") && identical(control_rule, "supplied_list")) {
        tryCatch(diaRuvSeparationRuleK(object, controls, grouping, parameters$num_components_to_impute),
            ruv_k_not_estimated = \(condition) condition)
    }
    sweep <- if (identical(k_rule, "sweep_bounded") && !identical(control_rule, "all_features") && is.null(separation)) {
        suppressMessages(findBestNegCtrlPercentageCapped(object,
            percentage_range = seq(parameters$percentage_min, parameters$percentage_max, by = 1),
            num_components_to_impute = parameters$num_components_to_impute,
            ruv_grouping_variable = grouping, ruv_qval_cutoff = parameters$ruv_qval_cutoff,
            ruv_fdr_method = parameters$ruv_fdr_method, separation_metric = parameters$separation_metric,
            k_penalty_weight = parameters$k_penalty_weight, adaptive_k_penalty = isTRUE(parameters$adaptive_k_penalty),
            verbose = FALSE))
    }
    resolveK <- function() {
        resolveRuvFactorCount(k_rule, control_rule, limit, k_user = parameters$ruv_k,
            # The count removes the groups, or the biological units on a design with technical replicates (M28).
            estimate_count = \() estimateRuvFactorCount(as.matrix(object@protein_quant_table[, samples, drop = FALSE]),
                design[[grouping]], seed = parameters$ruv_k_seed, technical_replicates = isTRUE(block$technical), units = units),
            separation_rule = \() if (inherits(separation, "condition")) stop(separation) else if (!is.null(separation)) separation$k else
                sweep$best_k)
    }
    keep <- \(condition) condition
    k <- if (identical(k_rule, "user")) {
        resolveK()
    } else {
        tryCatch(resolveK(), ruv_factor_count_too_few_features = keep, ruv_factor_count_too_few_samples = keep,
            ruv_k_above_limit = keep, ruv_k_not_estimated = keep)
    }
    if (inherits(k, "condition")) {
        # What the data cannot give under an automatic rule: no RUV for this run, said and recorded.
        reason <- conditionMessage(k)
        cli::cli_warn("RUV is not applied on this run: {reason}", class = "dia_ruv_not_applied")
        k <- list(k = 0L, k_rule = k_rule, k_estimate = NA_integer_, k_estimate_reason = reason, k_max = limit$k_max,
            control_rule = control_rule, status = "refused", reason = reason, estimate = NULL)
    }

    # --- Factors ---
    factors <- if (identical(k$status, "ok")) computeRuvCovariates(object, k$k, controls, grouping)
    list(status = k$status, reason = k$reason, k = k, factors = factors, controls = stats::setNames(controls, ids),
        controls_unmatched = unmatched, sweep = sweep, separation = if (!inherits(separation, "condition")) separation,
        replicate_structure = replicates$structure,
        block_reason = block$reason, blocks = replicates$blocks)
}

#' The kernel of the RUV covariates by declared rules (KIT-129 W3): `computeDiaRuvCovariatesByRule`, then the same
#' products as today's kernel. The output state carries the quantified matrix unchanged, `W1..Wk` in the design and
#' the declaration `args$ruvCovariates` with the rules and k's source; the RUV-adjusted matrix is a component
#' checkpoint for the QC plots.
#'
#' The product describes what this run used (audit F2 of 2026-11-18). A run without RUV (an estimate of 0, or a
#' refusal) returns the input as the state and no adjusted matrix of its own; its mode is `"none"` in the product,
#' in `ruv_result` and in the state's `args`, with `ruv_applied = FALSE` and the reason, and it names no controls
#' and no control count as used (`controls` NULL, `control_count` NA). The sweep's percentage is not a quantity of a
#' declared run (`percentage`, `best_percentage` NA); the count the control rule selected is
#' `ruv_controls_selected`. When the sweep gave k, its own record (its k, percentage, control count, scores, table
#' and plot) is `ruv_result$k_sweep`, never the run's own fields (`optimization_results` and `best_cancor_plot`
#' absent, the scores NA). The estimated count's record (`k_estimate`, for `plotRuvFactorCount`) is kept whenever
#' it was computed.
#' @param object,parameters,store As a declaration's kernel.
#' @return Today's product (`covariates`, `adjusted`, `final`, `ruv_mode`, `ruv_result`, `W`, `best_k`,
#'   `percentage`, `controls`, `control_count`, `replicate_structure`, `block_reason`, `blocks`, `parameters`) and
#'   `ruv_control_rule`, `ruv_k_rule`, `ruv_k_estimate`, `ruv_applied`, `ruv_not_applied_reason`,
#'   `ruv_controls_selected`.
#' @noRd
diaRuvCovariatesByRuleKernel <- function(object, parameters, store) {
    result <- computeDiaRuvCovariatesByRule(object, parameters)
    applied <- identical(result$status, "ok")
    factors <- result$factors
    sweep <- result$sweep
    best_k <- if (applied) factors$k else 0L
    ruv_mode <- if (applied) "covariates" else "none"
    controls <- if (applied) result$controls
    control_count <- if (applied) sum(result$controls) else NA_integer_
    rules <- list(ruv_control_rule = result$k$control_rule, ruv_k_rule = result$k$k_rule, ruv_k_estimate = result$k$k_estimate,
        ruv_applied = applied, ruv_not_applied_reason = result$reason, ruv_controls_selected = sum(result$controls))
    with_w <- object
    samples <- as.character(with_w@design_matrix[[with_w@sample_id]])
    for (i in seq_len(best_k)) with_w@design_matrix[[paste0("W", i)]] <- unname(factors$W[samples, i])
    with_w@args$ruvCovariates <- c(list(k = best_k, percentage_as_neg_ctrl = NA_real_, control_count = control_count,
        replicate_structure = result$replicate_structure, block_reason = result$block_reason, blocks = result$blocks,
        ruv_grouping_variable = parameters$ruv_grouping_variable, separation_score = NA_real_, composite_score = NA_real_,
        max_acceptable_k = NA_integer_, ruv_mode = ruv_mode), rules)
    covariates_reference <- store$checkpoint(with_w, "ruv-covariates-science")
    adjusted_reference <- covariates_reference
    if (applied) {
        adjusted <- object
        adjusted_table <- as.data.frame(factors$adjusted)
        adjusted_table[[adjusted@protein_id_column]] <- rownames(factors$adjusted)
        rownames(adjusted_table) <- NULL
        adjusted@protein_quant_table <- adjusted_table[, colnames(object@protein_quant_table), drop = FALSE]
        adjusted_reference <- store$checkpoint(adjusted, "ruv-adjusted-science", from = covariates_reference)
    }
    # KIT-129 sixth pass item 3: the separation rule on a supplied list keeps its k and plot where the sweep's are kept.
    k_sweep <- if (!is.null(result$separation)) {
        list(best_k = result$separation$k, control_count = sum(result$controls), rule = "separation_on_supplied_list",
            cancor_plot = ggplot2::ggplotGrob(stripProtNormQcPlotOwners(result$separation$plot)))
    } else if (!is.null(sweep)) {
        list(best_k = sweep$best_k, best_percentage = sweep$best_percentage,
            control_count = sum(sweep$best_control_genes_index, na.rm = TRUE), separation_score = sweep$best_separation_score,
            composite_score = sweep$best_composite_score, max_acceptable_k = sweep$max_acceptable_k,
            separation_metric = parameters$separation_metric, k_penalty_weight = parameters$k_penalty_weight,
            adaptive_k_penalty = isTRUE(parameters$adaptive_k_penalty), optimization_results = sweep$optimization_results,
            cancor_plot = if (!is.null(sweep$best_cancor_plot)) ggplot2::ggplotGrob(stripProtNormQcPlotOwners(sweep$best_cancor_plot)))
    }
    ruv_result <- c(list(best_percentage = NA_real_, best_k = best_k, best_control_genes_index = controls,
        best_separation_score = NA_real_, best_composite_score = NA_real_, max_acceptable_k = NA_integer_,
        sample_size = length(samples), ruv_mode = ruv_mode, replicate_structure = result$replicate_structure,
        block_reason = result$block_reason), rules,
        list(ruv_k_estimate_reason = result$k$k_estimate_reason, k_max = result$k$k_max, k_estimate = result$k$estimate,
            k_sweep = k_sweep, controls_unmatched = result$controls_unmatched,
            # KIT-129 fourth pass item 4: the factor-tied-to-groups check on this run's W (no threshold validated on this
            # route, decision M27: the value without a verdict).
            ruv_group_association = ruvGroupAssociationCheck(factors$W, ruvFactorGroups(factors$W, object@design_matrix,
                object@sample_id, object@group_id), "dia_limpa", NA_character_, result$k$k_rule)))
    c(list(covariates = covariates_reference, adjusted = adjusted_reference, final = covariates_reference,
        ruv_mode = ruv_mode, ruv_result = ruv_result, W = factors$W, best_k = best_k, percentage = NA_real_,
        controls = controls, control_count = control_count, replicate_structure = result$replicate_structure,
        block_reason = result$block_reason, blocks = result$blocks,
        parameters = parameters[c("percentage_min", "percentage_max", "num_components_to_impute", "ruv_grouping_variable",
            "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_mode",
            diaRuvCovariatesRuleParameters)]), rules)
}

#' Run or reuse the RUV-covariates declaration on the DIA_LIMPA route and publish its outputs (KIT-041): after
#' the correlation filter, the state `ruv_covariates` (label `13_ruv_covariates`) commits the quantified matrix
#' unchanged with the W design; the adjusted matrix's QC plots, the cancor plot, the optimisation table and the
#' W design are published as the step's files.
#' @param settings `ruv_mode` (`covariates`, or `correction` since KIT-056: the state `ruv_corrected` commits the
#'   adjusted matrix as the protein table and the owned EList's E, no W terms), `color_variable`, `shape_variable`.
#' @param workflow,norm_data,paths Existing GUI owners and publication paths.
#' @param ruv_qc The module's RUV-corrected QC callback (takes a reference).
#' @param grouping Function returning the RUV grouping variable.
#' @return The accepted `ruv_covariates` reference.
#' @noRd
runDIAReferenceRuvCovariates <- function(settings, workflow, norm_data, paths, omic_type, experiment_label, ruv_qc, grouping) {
    prior_plot <- ggplot2::get_last_plot()
    on.exit(ggplot2::set_last_plot(prior_plot), add = TRUE)
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "ruv-covariates")
        source <- diaRuvCovariatesSource(manager, norm_data)
        parameters <- diaRuvCovariatesRequestParameters(workflow$config_list, grouping(), settings)
        run <- computeDIAReferenceRuvCovariates(manager, source, parameters)
        commit <- diaRuvCovariatesCommit(manager, parameters, run)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        publishDIAReferenceRuvCovariates(manager, settings, run, workflow, norm_data, paths, ruv_qc)
    })
}

#' The RUV covariates' source, as the tab sends it in process and to the worker (GUI audit P-35, F111-02): the recorded
#' input, else the tab's correlation result, else the head. The source is recorded as the step's input, as the
#' normalisation records its own (`diaNormalizationRequest`), because the publication half points the tab's correlation
#' result at the committed RUV state: a repeated click or a mode switch then starts from the correlation result again,
#' never from a state that already holds the correction or the W terms. The correlation's publication clears the record.
#' @param manager,norm_data Session owner and the tab's state.
#' @return The source reference.
#' @noRd
diaRuvCovariatesSource <- function(manager, norm_data) {
    source <- norm_data$ruv_covariates_input_reference %||% norm_data$correlation_filtered_obj %||% getDIAStateReference(manager)
    if (!isSavedWorkflowReference(source)) cli::cli_abort("the correlation step must complete before the RUV covariates")
    norm_data$ruv_covariates_input_reference <- source
    source
}

#' The RUV covariates' compute half (KIT-111 D3): the declaration on the source, reused on a hit. Worker-safe.
#' @return `runSavedOperation`'s run.
#' @noRd
computeDIAReferenceRuvCovariates <- function(manager, source, parameters) {
    # KIT-129 integration item 1: the declared rules' declaration when the request names a control rule, on both paths.
    chosen <- workflowRequestOperation(manager, workflowOperationEntry("protein_ruv_covariates"), source, parameters)
    # KIT-129 audit F3: before the runner, whose identity lookup would return today's saved product for a request
    # that carries rule parameters.
    if (identical(chosen$declaration$parameters, workflowOperationEntry("protein_ruv_covariates")$parameters)) assertDiaRuvSweepParameters(parameters)
    runSavedOperation(manager, chosen$declaration, source, chosen$parameters)
}

#' The RUV covariates' parameters the tab sends, in process and to the worker (KIT-129 integration item 1): today's
#' (`resolveDIARuvCovariatesParameters`), and with the gate on the rules of a declared design the settings name
#' (`ruv_control_rule`, `ruv_k_rule`, `ruv_k`, `ruv_control_ids`, `ruv_k_seed`).
#' @param config_list,grouping As `resolveDIARuvCovariatesParameters`.
#' @param settings The tab's settings (`ruv_mode`, and the rules).
#' @return The declaration's parameters.
#' @noRd
diaRuvCovariatesRequestParameters <- function(config_list, grouping, settings) {
    mode <- settings$ruv_mode %||% "covariates"
    if (!omicGuiRuvDesignEnabled() || is.null(settings$ruv_control_rule)) return(resolveDIARuvCovariatesParameters(config_list, grouping, mode))
    resolveDIARuvCovariatesParameters(config_list, grouping, mode, ruv_control_rule = settings$ruv_control_rule,
        ruv_k_rule = settings$ruv_k_rule, ruv_k = settings$ruv_k, ruv_control_ids = settings$ruv_control_ids, ruv_k_seed = settings$ruv_k_seed)
}

#' The RUV covariates' commit description (KIT-111 D3), shared by the in-process runner and the worker's completion.
#' Today's product commits today's fields. A declared-rule product (KIT-129) adds its rules to the configuration
#' (`ruv_control_rule`, `ruv_k_rule`, `ruv_k_estimate`, `ruv_not_applied_reason`, `ruv_controls_selected`) and to
#' the outcomes; a run without RUV commits the state `ruv_covariates` with `ruv_applied = FALSE` and mode `none`.
#' @param parameters The declaration's parameters (`resolveDIARuvCovariatesParameters`).
#' @param run The compute half's run.
#' @return `list(state_name, checkpoint, config, description)`, or NULL when the head already holds the state with
#'   this product.
#' @noRd
diaRuvCovariatesCommit <- function(manager, parameters, run) {
    saved <- run$value
    ruv_mode <- saved$ruv_mode %||% "covariates"
    correction <- identical(ruv_mode, "correction")
    # KIT-129 (audit F2): whether RUV was applied is the product's own statement. Today's product (the sweep) makes
    # none and always applies; a declared-rule product says, with its reason and its rules.
    declared <- !is.null(saved$ruv_control_rule)
    applied <- !identical(saved$ruv_applied, FALSE)
    state_name <- if (correction) "ruv_corrected" else "ruv_covariates"
    final <- readDIAStateMetadata(manager, saved$final %||% saved$covariates)$data
    final$products$ruv_covariates_science <- run$reference
    config <- list(ruv_mode = ruv_mode, ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage,
        control_count = saved$control_count, separation_score = saved$ruv_result$best_separation_score,
        composite_score = saved$ruv_result$best_composite_score, max_acceptable_k = saved$ruv_result$max_acceptable_k,
        replicate_structure = saved$replicate_structure, block_reason = saved$block_reason, replicate_blocks = saved$blocks,
        ruv_grouping_variable = parameters$ruv_grouping_variable, ruv_applied = applied,
        outcomes = list(k = saved$best_k, percentage_as_neg_ctrl = saved$percentage, control_count = saved$control_count,
            separation_score = saved$ruv_result$best_separation_score, composite_score = saved$ruv_result$best_composite_score,
            replicate_structure = saved$replicate_structure, block_reason = saved$block_reason, ruv_mode = ruv_mode))
    description <- if (correction) {
        sprintf("RUV-III correction applied to the quantified matrix and the DPC EList (k = %d; Y - W alpha); no W terms in the design", saved$best_k)
    } else {
        sprintf("RUV-III factors appended to the design as covariates (W1..W%d); the quantified matrix unchanged", saved$best_k)
    }
    if (declared) {
        # What this run used, from its product: the control rule, k and the rule that gave it, and the reason when
        # no RUV was applied.
        config <- c(config, saved[c("ruv_control_rule", "ruv_k_rule", "ruv_k_estimate", "ruv_not_applied_reason", "ruv_controls_selected")])
        config$outcomes <- c(config$outcomes, list(ruv_applied = applied, control_rule = saved$ruv_control_rule,
            k_rule = saved$ruv_k_rule, k_estimate = saved$ruv_k_estimate, not_applied_reason = saved$ruv_not_applied_reason,
            controls_selected = saved$ruv_controls_selected,
            # KIT-129 integration item 3: the count's seed, with the rule and the estimate.
            k_seed = saved$parameters$ruv_k_seed %||% parameters$ruv_k_seed),
            # KIT-129 fourth pass item 4: the factor-tied-to-groups check (statistic, threshold, outcome).
            saved$ruv_result$ruv_group_association[ruvGroupAssociationRecordFields])
        description <- if (applied) {
            sprintf("RUV-III factors appended to the design as covariates (W1..W%d; controls by the rule %s, %d proteins; k by the rule %s); the quantified matrix unchanged",
                saved$best_k, saved$ruv_control_rule, saved$control_count, saved$ruv_k_rule)
        } else {
            sprintf("RUV not applied (%s); the quantified matrix and the design unchanged", saved$ruv_not_applied_reason)
        }
    }
    previous <- if (identical(manager$current_state, state_name)) readDIAStateMetadata(manager) else NULL
    if (!is.null(previous) && identical(previous$data$products$ruv_covariates_science, run$reference)) return(NULL)
    list(state_name = state_name, checkpoint = final, config = config, description = description)
}

#' The RUV covariates' publication half (KIT-111 D3), in the controller on both paths: the tab's state and the step's
#' files and display, on the accepted head.
#' @param run The compute half's run (`value`, `reference`); on the worker path the delivered product.
#' @inheritParams runDIAReferenceRuvCovariates
#' @return The accepted reference.
#' @noRd
publishDIAReferenceRuvCovariates <- function(manager, settings, run, workflow, norm_data, paths, ruv_qc) {
    saved <- run$value
    product <- run$reference
    norm_data$ruv_optimization_result <- saved$ruv_result
    workflow$ruv_optimization_result <- saved$ruv_result
    norm_data$ruv_optimization_table <- saveDIADisplayTable(manager, saved$ruv_result$optimization_results)
    # KIT-129 fourth pass item 3: what a declared rule's run did, for the RUV design section (nothing with the gate off).
    omicNormSetRuvDesignResult(norm_data, diaRuvDesignResult(saved, settings$ruv_design))
    norm_data$best_k <- saved$best_k
    norm_data$control_genes_index <- saved$controls
    accepted <- getDIAStateReference(manager)
    norm_data$ruv_normalized_obj <- accepted
    norm_data$correlation_filtered_obj <- accepted
    workflow$ruv_normalised_for_da_analysis_obj <- accepted
    # KIT-052: drawn in a short folder and copied into the action staging folder, as the normalisation.
    display <- publishDIAActionOutputs(manager, accepted, "ruv-covariates", paths[c("source_dir", "protein_qc_dir")],
        function(staging) withWorkflowShortStaging(staging, function(staged) {
            norm_data$publication_paths <- staged
            on.exit({ norm_data$publication_paths <- NULL }, add = TRUE)
            writeDIADesignTable(manager, readDIAStateMetadata(manager, saved$covariates)$data$components$design_matrix,
                file.path(staged$protein_qc_dir, "ruv_covariates_design.tsv"), vroom::vroom_write)
            adjusted_spec <- readDIAStateMetadata(manager, saved$adjusted)$data
            writeDIADesignTable(manager, adjusted_spec$components$protein_quant_table,
                file.path(staged$protein_qc_dir, "ruv_adjusted_protein_matrix.tsv"), vroom::vroom_write)
            # KIT-129 fourth pass item 3: a declared rule's run has no sweep table unless the sweep gave k.
            if (is.data.frame(saved$ruv_result$optimization_results)) {
                vroom::vroom_write(saved$ruv_result$optimization_results, file.path(staged$protein_qc_dir, "ruv_covariates_optimization.tsv"))
            }
            saveRDS(saved$ruv_result[setdiff(names(saved$ruv_result), "best_cancor_plot")],
                file.path(staged$source_dir, "ruv_covariates_optimization_results.RDS"))
            # KIT-259 second correction (VAL-287; AUDIT-KIT-259-impl I04): guarded as the normalisation's cancor figure.
            if (!is.null(saved$ruv_result$best_cancor_plot)) drawProtNormQcFigure(staged$protein_qc_dir,
                "ruv_covariates_cancor.png", 8, 6, function() ggplot2::ggsave(
                file.path(staged$protein_qc_dir, "ruv_covariates_cancor.png"),
                saved$ruv_result$best_cancor_plot, width = 8, height = 6, dpi = 150))
            ruv_qc(saved$adjusted)
            list(completed = TRUE, qc_paths = norm_data$qc_plot_paths, staged_qc_dir = staged$protein_qc_dir)
        }, "the RUV covariates' files could not be copied into the project's staging folder"), identity = list(product = product$entity_id, color = settings$color_variable,
            shape = settings$shape_variable, output_schema = 1L))
    # GUI audit A111-08 (F111-05): the tab's QC images are the published files, not the staging copy, as after the
    # normalisation. A display saved before this fix names no staging directory and is left as it was.
    if (!is.null(display$staged_qc_dir) && length(display$qc_paths)) {
        norm_data$qc_plot_paths <- lapply(display$qc_paths, function(group) lapply(group, function(path) {
            if (is.character(path)) sub(display$staged_qc_dir, paths$protein_qc_dir, path, fixed = TRUE) else path
        }))
    }
    norm_data$ruv_covariates_complete <- TRUE
    accepted
}

#' Apply, skip or reapply Pearson filtering from the original saved population.
#' @param workflow,norm_data,paths Existing module owners and publication paths.
#' @param threshold,skipped User filtering decision.
#' @param omic_type,experiment_label Existing output labels.
#' @param grouping,aesthetics Current scientific grouping and render selections.
#' @return Final reference, after required output completion; the module edge completes the GUI state.
#' @noRd
runDIAReferenceCorrelation <- function(workflow, norm_data, paths, threshold,
    skipped, omic_type, experiment_label, grouping, aesthetics) {
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- diaCorrelationSource(manager, norm_data)
        # KIT-257 (VAL-250): the earlier results table (either name) goes before anything is computed or committed.
        removed <- removeProtNormOwnFiles(paths, "correlation")
        withProtNormRemovedFiles(removed, {
            original <- readDIAStateMetadata(manager, source)$data
            run <- if (skipped) list(final = source) else computeDIAReferenceCorrelation(manager, source, grouping(), threshold)
            commit <- diaCorrelationCommit(manager, run$final, threshold, skipped)
            if (!is.null(commit)) {
                # KIT-248 (VAL-238): the Skip records itself, only when it commits (a second Skip commits nothing and leaves
                # nothing pending); a stale pending Apply is dropped so it is not chained into the Skip's entry.
                if (skipped) {
                    clearDIAStepProvenance(manager)
                    setDIAStepProvenance(manager, "pearson-skip", "whole", character(), fixed = list(skipped = TRUE))
                }
                commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            }
            publishDIAReferenceCorrelation(manager, run, workflow, norm_data, paths, threshold, skipped, omic_type,
                experiment_label, aesthetics)
        })
    })
}

#' The correlation filter's source, as the tab sends it in process and to the worker (GUI audit P-35): the normalisation's
#' final. On the DIA_LIMPA route that is the pass-through's recorded input (`normalization_input_reference`), not the
#' tab's `ruv_normalized_obj`, which the RUV covariates' publication points at the committed RUV state: a correlation
#' re-apply after the covariates then starts where the first correlation did, never from a state that already holds the
#' correction or the W terms. On a first click the two are the same reference. Every other route reads
#' `ruv_normalized_obj`, as before. KIT-182 (VAL-91 D2, VAL-99 (4)): refused when the Quality Control tab's table has
#' changed since normalisation (`assertDIACorrelationHeadCurrent`).
#' @param manager,norm_data Session owner and the tab's state.
#' @return The source reference.
#' @noRd
diaCorrelationSource <- function(manager, norm_data) {
    source <- if (identical(manager$workflow_type, "DIA_LIMPA"))
        norm_data$normalization_input_reference %||% norm_data$ruv_normalized_obj else norm_data$ruv_normalized_obj
    if (!isSavedWorkflowReference(source)) cli::cli_abort("normalization must complete before correlation filtering")
    if (usesDIAReferenceStates(manager)) assertDIACorrelationHeadCurrent(manager, source)
    source
}

#' The correlation filter's refusal after a Quality Control change (KIT-182, VAL-91 D2; the ticket's 7.5), for Apply and
#' Skip on both paths. DIA, LFQ, TMT: allowed at the tab's own states, at `normalized` or `ruv_corrected` only when the
#' head's table is the source's (VAL-99 (4): a Revert between two normalisations leaves an earlier run's state current
#' while the tab holds the later run's result); refused at any other protein or peptide head (a step applied or
#' reverted since normalisation). DIA_LIMPA (its normalisation commits no state): allowed at its tab states
#' (`ruv_covariates`, `ruv_corrected`, `correlation_filtered`) or when the head is the table the pass-through designated;
#' refused at any other protein or peptide head. A result-class head or an unreadable record: allowed, as today.
#' @param manager Reference-backed session owner.
#' @param source The correlation's source (`diaCorrelationSource`).
#' @return TRUE invisibly, or an error of class `prot_norm_correlation_input_stale`.
#' @noRd
assertDIACorrelationHeadCurrent <- function(manager, source) {
    refuse <- function() rlang::abort(paste0("the table in the Quality Control tab has changed since normalisation was ",
        "last run (a step was applied or reverted), so the correlation filter would use a table that does not match it. ",
        "Press \"Run Normalization & RUV\" again, then apply or skip the correlation filter."),
        class = "prot_norm_correlation_input_stale")
    limpa <- identical(manager$workflow_type, "DIA_LIMPA")
    head <- getDIAStateReference(manager)
    if (!limpa && manager$current_state %in% c("normalized", "ruv_corrected")) {
        head_spec <- readDIALineageRecord(manager, head)
        source_spec <- readDIALineageRecord(manager, source)
        if (!is.null(head_spec) && !is.null(source_spec) && !identical(head_spec$components, source_spec$components)) refuse()
        return(invisible(TRUE))
    }
    own <- if (limpa) c("ruv_covariates", "ruv_corrected", "correlation_filtered") else "correlation_filtered"
    if (manager$current_state %in% own) return(invisible(TRUE))
    if (limpa && identical(head$entity_id, source$entity_id)) return(invisible(TRUE))
    spec <- readDIALineageRecord(manager, head)
    if (is.null(spec) || !isTRUE(spec$class %in% c("ProteinQuantitativeData", "PeptideQuantitativeData")))
        return(invisible(TRUE))
    refuse()
}

#' The within-group Pearson pairs' parameters (KIT-111 D3): the grouping and the pool-sample pattern the tab removes.
#' @noRd
diaPearsonPairsParameters <- function(grouping) list(group = grouping, tech_rep_remove_regex = "pool")

#' The correlation filter's compute half (KIT-111 D3): the pairs, then the selection on them at the threshold, each
#' reused on a hit. Worker-safe; the worker runs the two as two actions (`computeWorkflowDiaPearsonPairsAction`,
#' `computeWorkflowDiaPearsonSelectionAction`).
#' @return `list(pairs, final)`: the pairs product and the selection's checkpoint.
#' @noRd
computeDIAReferenceCorrelation <- function(manager, source, grouping, threshold) {
    pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), source, diaPearsonPairsParameters(grouping))
    selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), source,
        list(pairs = pairs$reference$entity_id, threshold = threshold, correlations = pairs$value),
        dependencies = list(pairs$reference))
    list(pairs = pairs$reference, final = selection$value)
}

#' The correlation filter's commit description (KIT-111 D3), shared by the in-process runner and the worker's completion.
#' @param final The selection's checkpoint (the source when skipped).
#' @param threshold,skipped The user's filtering decision.
#' @return `list(state_name, checkpoint, config, description)`, or NULL when the head is `correlation_filtered` with
#'   the same components and config.
#' @noRd
diaCorrelationCommit <- function(manager, final, threshold, skipped) {
    spec <- readDIAStateMetadata(manager, final)$data
    config <- if (skipped) list(min_pearson_correlation_threshold = 0, skipped = TRUE) else
        list(min_pearson_correlation_threshold = threshold)
    previous <- if (manager$current_state == "correlation_filtered") readDIAStateMetadata(manager) else NULL
    if (!is.null(previous) && identical(previous$data$components, spec$components) &&
        identical(previous$config, config)) return(NULL)
    list(state_name = "correlation_filtered", checkpoint = spec, config = config,
        description = if (skipped) "Skipped final sample correlation filter" else "Applied final sample correlation filter (chunk 28)")
}

#' The correlation filter's publication half (KIT-111 D3), in the controller on both paths: the tab's state and the
#' step's files and display, on the accepted head.
#' @param run `list(pairs, final)` (`computeDIAReferenceCorrelation`; `final` alone when skipped).
#' @inheritParams runDIAReferenceCorrelation
#' @return The selection's checkpoint (the source when skipped).
#' @noRd
publishDIAReferenceCorrelation <- function(manager, run, workflow, norm_data, paths, threshold, skipped, omic_type,
    experiment_label, aesthetics) {
    final <- run$final
    if (!skipped) norm_data$correlation_vector <- run$pairs
    norm_data$correlation_filtered_obj <- final
    # GUI audit P-35: a new correlation result is the RUV covariates' new input.
    if (!is.null(norm_data$ruv_covariates_input_reference)) norm_data$ruv_covariates_input_reference <- NULL
    norm_data$correlation_threshold <- if (skipped) 0 else threshold
    spec <- readDIAStateMetadata(manager, final)$data
    accepted <- getDIAStateReference(manager)
    display <- publishDIAActionOutputs(manager, accepted, "correlation",
        paths[c("protein_qc_dir", "time_dir")], function(staged) {
            summary <- readDIAProteinProgress(manager, final)
            directories <- stats::setNames(list(list(time_dir = staged$time_dir,
                publication_graphs_dir = paths$publication_graphs_dir)),
                paste(omic_type, experiment_label, sep = "_"))
            # KIT-041: the DIA_LIMPA route labels the correlation step as the limpa script does.
            step_label <- if (identical(manager$workflow_type, "DIA_LIMPA")) "12_samplecorr_filtered" else "12_correlation_filtered"
            filtering_plot <- updateProteinFiltering(summary, step_label,
                omic_type, experiment_label, overwrite = TRUE, return_grid = TRUE, project_dirs = directories)
            plot <- resolveDIAReferenceQcPlot(final, aesthetics(), "pca", buildProtNormQcPlotProduct,
                plotPca, plotRle, buildProtNormDensityPlot, buildProtNormCorrelationPlot)
            qc_plot <- stripProtNormQcPlotOwners(plot + ggplot2::labs(title = if (skipped)
                "Final Data (Correlation Filter Skipped)" else "Final Correlation-Filtered Data"))
            names <- resolveProtNormCorrelationResultFilenames(norm_data$ruv_optimization_result,
                workflow$ruv_optimization_result)
            writeDIADesignTable(manager, spec$components$protein_quant_table,
                file.path(staged$protein_qc_dir, names$tsvFilename), readr::write_tsv)
            list(filtering_plot = filtering_plot, qc_plot = qc_plot)
        }, identity = list(result = final$entity_id, skipped = skipped, threshold = threshold,
            aesthetics = aesthetics(), output_schema = 2L))
    norm_data$final_filtering_plot <- display$filtering_plot
    norm_data$final_qc_plot <- display$qc_plot
    final
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-111 D3: the normalisation tab's declarations on the managed worker's `omics` role. Each handler runs its compute
# half on the request's state and returns references and the commit description only; the controller's completion
# (`completeWorkflowOmicsAction`) commits, and the tab's delivery runs the publication half. A hit is rooted with the
# saved references the run names, never with the product's value (F111-01).
# ---------------------------------------------------------------------------------------------------------------

#' The worker's value for a state-producing DIA action: the committed record (the head's when the commit description is
#' NULL), the product and the commit description without its checkpoint.
#' @noRd
diaWorkerCommitValue <- function(manager, run, commit) {
    list(checkpoint = if (is.null(commit)) readWorkflowStateMetadata(manager)$data else commit$checkpoint,
        reference = run$reference, commit = if (!is.null(commit)) commit[c("state_name", "config", "description")])
}

#' The normalisation on the worker (`normalization-ruv`).
#' @noRd
computeWorkflowDiaNormalizationAction <- function(manager, request, variant) {
    source <- request$input_references$state
    run <- computeDIAReferenceNormalization(manager, source, request$parameters)
    if (identical(run$status, "hit")) rootDIAProduct(manager, run$reference, list(source, run$value$normalized, run$value$final))
    diaWorkerCommitValue(manager, run, diaNormalizationCommit(manager, request$parameters, run))
}

#' The within-group Pearson pairs on the worker (`pearson-pairs`, a product): the source's record is returned.
#' @noRd
computeWorkflowDiaPearsonPairsAction <- function(manager, request, variant) {
    source <- request$input_references$state
    run <- runSavedOperation(manager, diaPearsonPairsOperation(), source, request$parameters)
    if (identical(run$status, "hit")) rootDIAProduct(manager, run$reference, list(source))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = run$reference)
}

#' The sample selection on the worker (`pearson-selection`): the correlations come from the pairs product its `pairs`
#' parameter names in the session; the request carries no vector.
#' @noRd
computeWorkflowDiaPearsonSelectionAction <- function(manager, request, variant) {
    source <- request$input_references$state
    pairs <- manager$.context$snapshot$products[[request$parameters$pairs]]
    if (is.null(pairs)) cli::cli_abort("the Pearson pairs product is not in the session", class = "workflow_operation_parameters")
    parameters <- c(request$parameters, list(correlations = hydrateSavedWorkflowEntity(manager$.context$store, pairs)$value))
    run <- runSavedOperation(manager, diaPearsonSelectionOperation(), source, parameters, dependencies = list(pairs))
    if (identical(run$status, "hit")) rootDIAProduct(manager, run$reference, list(source, pairs, run$value))
    diaWorkerCommitValue(manager, run, diaCorrelationCommit(manager, run$value, request$parameters$threshold, FALSE))
}

#' The RUV covariates on the worker (`protein_ruv_covariates`, DIA_LIMPA).
#' @noRd
computeWorkflowDiaRuvCovariatesAction <- function(manager, request, variant) {
    source <- request$input_references$state
    run <- computeDIAReferenceRuvCovariates(manager, source, request$parameters)
    if (identical(run$status, "hit"))
        rootDIAProduct(manager, run$reference, c(list(source), run$value[c("covariates", "adjusted", "final")]))
    diaWorkerCommitValue(manager, run, diaRuvCovariatesCommit(manager, request$parameters, run))
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-122 (audit E1, F3): the normalisation tab's QC figures as a worker product, as the lipid and metab tabs' are
# (`lipid_qc_figures`, `lipid_final_pca`). The action commits nothing. It saves, for the request's named state, what the
# tab's generators and publication halves would otherwise compute in the application: the QC plot products
# (`resolveDIAReferenceQcPlot`, which hydrates the state on a miss) and, for the stage figures, their rendered images
# (`saveProtNormQcPlotArtifact`). The controller's completion roots them. The tab's code is unchanged: run after the
# delivery, it finds the product saved (no state is hydrated, nothing is computed) and writes the saved image's bytes
# (nothing is drawn). With no runner, or if a figure is not found saved, it computes as before.
# ---------------------------------------------------------------------------------------------------------------

#' The stage figures' kinds, in the generators' order, with the width and height in inches each is saved at
#' (`generateProtNorm*QcArtifacts`, mod_prot_norm_qc_generation_helpers.R).
#' @noRd
diaQcStageFigureSizes <- list(pca = c(8, 6), rle = c(10, 6), density = c(8, 6), correlation = c(10, 8))

#' The QC figures action's parameters: the aesthetics the tab draws with, and which figures.
#' @param aesthetics `list(color_var, shape_var)`.
#' @param figures `stage` (the four figures of a stage, with their images), `final_pca` (the correlation tab's PCA,
#'   which the tab draws from the saved product; no image) or `composite` (the normalisation's composite QC figure, its
#'   three files; `saveWorkflowDiaCompositeFigure`).
#' @param composite For `composite` only: `list(product, ruv_mode)`, the entity id of the normalisation run's product
#'   and the run's RUV mode.
#' @return The parameters.
#' @noRd
diaQcFiguresParameters <- function(aesthetics, figures = "stage", composite = NULL) {
    c(list(color_variable = aesthetics$color_var, shape_variable = aesthetics$shape_var, figures = figures), composite)
}

#' Index, in the controller's store handle, the products the head roots, when a QC figures action is delivered. The
#' worker saved them under its own handle and the completion rooted them in the head, but an identity lookup
#' (`resolveSavedWorkflowEntity`) reads the handle's own index, which is otherwise rebuilt from the head's roots only when
#' a store is reopened (`newDIAStateAdapter`). Without this the tab would not find the worker's products, would compute
#' and save them again, and the two saves of one identity would be refused as ambiguous by the next reader.
#' @param manager The session owner.
#' @return NULL, invisibly.
#' @noRd
indexWorkflowHeadProducts <- function(manager) {
    store <- manager$.context$store
    for (reference in manager$.context$snapshot$products) registerSavedWorkflowMetadata(store, reference)
    invisible(NULL)
}

#' The normalisation tab's QC figures on the worker (`normalization-qc-figures`): the products and images of the
#' request's state under the request's aesthetics; for `composite`, the composite QC figure of the normalisation run
#' whose final state the request names.
#' @noRd
computeWorkflowDiaQcFiguresAction <- function(manager, request, variant) {
    source <- request$input_references$state
    parameters <- request$parameters
    if (!isTRUE(parameters$figures %in% c("stage", "final_pca", "composite")))
        cli::cli_abort("the QC figures must be stage, final_pca or composite", class = "workflow_operation_parameters")
    aesthetics <- list(color_var = parameters$color_variable, shape_var = parameters$shape_variable)
    plot <- function(kind) resolveDIAReferenceQcPlot(source, aesthetics, kind, buildProtNormQcPlotProduct,
        plotPca, plotRle, buildProtNormDensityPlot, buildProtNormCorrelationPlot)
    if (identical(parameters$figures, "final_pca")) {
        plot("pca")
    } else if (identical(parameters$figures, "composite")) {
        saveWorkflowDiaCompositeFigure(manager, request, aesthetics)
    } else {
        # KIT-131: in a worker, not under its request-owned temporary folder (on Windows past 259 characters).
        directory <- tempfile("dia-qc-figures-", tmpdir = workflowShortScratchRoot())
        dir.create(directory)
        on.exit(unlink(directory, recursive = TRUE), add = TRUE)
        # KIT-259 (B2): a figure that cannot be drawn is skipped and logged; the publication then draws it, or replaces
        # it by its placeholder.
        for (kind in names(diaQcStageFigureSizes)) {
            size <- diaQcStageFigureSizes[[kind]]
            drawProtNormQcFigure(directory, paste0(kind, ".png"), size[[1L]], size[[2L]], function()
                saveProtNormQcPlotArtifact(directory, paste0(kind, ".png"), plot(kind), width = size[[1L]], height = size[[2L]]),
                placeholder = FALSE, signal = FALSE)
        }
    }
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source)
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-122 (composite): the normalisation's composite QC figure as a worker product. Assembled in the application, its
# rasters are the normalisation step's largest memory event there (as the lipid tab's were,
# `computeWorkflowOmicsCompositeFigureAction`). The publication stays the journaled transaction it is: the figure's
# three files are written inside its `produce`, where the builder wrote them. On a managed session the worker has built
# the figure before the publication half runs (the `composite` figures of `normalization-qc-figures`) and saved the
# files; the publication half finds them by identity and writes their bytes. With no runner, or when a file is not
# found saved, it builds the figure as before and saves nothing.
# ---------------------------------------------------------------------------------------------------------------

#' The composite QC figure's files, by kind, with the media type each is saved under. The store holds each as a
#' rendered asset under an identity of its own, as it holds the interactive volcano pages
#' (`saveDIAInteractiveVolcanoAsset`): `newSavedWorkflowRenderIdentity` names image and page formats only, and the
#' `rds` is the plot object `savePlot` writes beside its renderings.
#' @noRd
diaCompositeFigureMediaTypes <- c(png = "image/png", pdf = "application/pdf", rds = "application/octet-stream")

#' The identities of the composite QC figure's files for the stage images in a QC directory: the normalisation run's
#' product, the RUV mode (the layout), each stage image there by name and content (the builder's inputs,
#' `resolveProtNormCompositeFigureInputs`), and the file's kind. The images' content carries the aesthetics they were
#' drawn with and a missing or older pre-normalisation image, so a figure is found saved only for the images it was
#' built from. The canonical-correlation image is not read here: it is the product's own plot at a fixed size. The
#' worker computes the identities on the directory it builds the figure in, the application on the publication's
#' staged directory.
#' @param manager The session owner.
#' @param product The normalisation run's product reference.
#' @param ruv_mode The run's RUV mode.
#' @param qc_dir The directory that holds the stage images.
#' @return The identities, named by kind (`diaCompositeFigureMediaTypes`).
#' @noRd
diaCompositeFigureIdentities <- function(manager, product, ruv_mode, qc_dir) {
    files <- sort(list.files(qc_dir, pattern = "^(pre_norm|post_norm|ruv_corrected)_(pca|density|rle|correlation)[.]png$"),
        method = "radix")
    images <- stats::setNames(vapply(file.path(qc_dir, files), function(file) digest::digest(file = file, algo = "sha256"), ""),
        files)
    parameters <- digest::digest(list(ruv_mode = ruv_mode, images = images), algo = "sha256")
    lapply(stats::setNames(nm = names(diaCompositeFigureMediaTypes)), function(kind) newSavedWorkflowScientificIdentity(
        product$entity_id, "dia-composite-qc-figure-v1", parameters, kind, manager$.context$software))
}

#' Whether a composite QC figure's file is saved under `identity`, without reading it.
#' @noRd
diaCompositeFigureSaved <- function(manager, identity) {
    identical(resolveSavedWorkflowEntity(manager$.context$store, "rendered_asset", identity, hydrate = FALSE)$status, "hit")
}

#' The composite QC figure on the worker (the `composite` figures of `normalization-qc-figures`). In a temporary
#' directory the worker gathers what the publication half stages before it writes the figure
#' (`publishDIAReferenceNormalization`): the published pre-normalisation images, the images of the normalised state and,
#' unless RUV was skipped, of the final one (the stage figures saved them before this action, so their bytes are
#' written and nothing is drawn again) and the canonical-correlation image. It runs the tab's unchanged builder there
#' and saves the files it wrote, rooted on the run's product. A figure already saved for these images is not built again.
#' @param manager The worker's owner.
#' @param request The request: `parameters$product` names the normalisation run's product in the session,
#'   `parameters$ruv_mode` its RUV mode.
#' @param aesthetics `list(color_var, shape_var)`, the stage figures'.
#' @return NULL, invisibly.
#' @noRd
saveWorkflowDiaCompositeFigure <- function(manager, request, aesthetics) {
    store <- manager$.context$store
    parameters <- request$parameters
    product <- manager$.context$snapshot$products[[parameters$product]]
    if (is.null(product)) cli::cli_abort("the normalisation product is not in the session", class = "workflow_operation_parameters")
    saved <- hydrateSavedWorkflowEntity(store, product)$value
    skipped <- identical(parameters$ruv_mode, "skip")
    directory <- tempfile("dia-qc-composite-", tmpdir = workflowShortScratchRoot())
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    published <- workflowRequestExperimentPaths(request)$protein_qc_dir
    if (!is.null(published)) file.copy(list.files(published, pattern = "^pre_norm_.*[.]png$", full.names = TRUE), directory)
    states <- c(list(post_norm = saved$normalized), if (!skipped) list(ruv_corrected = saved$final))
    # KIT-259 (B2): a figure that cannot be drawn is replaced by the publication's placeholder, so that the composite's
    # identity matches the publication's staged images.
    for (stage in names(states)) for (kind in names(diaQcStageFigureSizes)) {
        size <- diaQcStageFigureSizes[[kind]]
        drawProtNormQcFigure(directory, paste0(stage, "_", kind, ".png"), size[[1L]], size[[2L]], function() {
            plot <- resolveDIAReferenceQcPlot(states[[stage]], aesthetics, kind, buildProtNormQcPlotProduct,
                plotPca, plotRle, buildProtNormDensityPlot, buildProtNormCorrelationPlot)
            saveProtNormQcPlotArtifact(directory, paste0(stage, "_", kind, ".png"), plot, width = size[[1L]], height = size[[2L]])
        }, signal = FALSE)
    }
    identities <- diaCompositeFigureIdentities(manager, product, parameters$ruv_mode, directory)
    absent <- names(identities)[!vapply(identities, function(identity) diaCompositeFigureSaved(manager, identity), NA)]
    if (!length(absent)) return(invisible(NULL))
    if (!skipped && !is.null(saved$ruv_result$best_cancor_plot)) drawProtNormQcFigure(directory, "ruv_corrected_cancor.png",
        8, 6, function() ggplot2::ggsave(file.path(directory, "ruv_corrected_cancor.png"), saved$ruv_result$best_cancor_plot,
            width = 8, height = 6, dpi = 150), signal = FALSE)
    omic_type <- workflowActionSpec(request$operation)$omic_type
    generateProtNormCompositeQcFigure(parameters$ruv_mode, directory, omic_type)
    files <- stats::setNames(file.path(directory, paste0(omic_type, "_composite_QC_figure.", names(identities))), names(identities))
    # The builder reports its own failure and leaves no file; nothing is saved, and the application builds as before.
    if (!all(file.exists(files))) return(invisible(NULL))
    # The figure's size in pixels, from the image's header (the asset envelope needs positive dimensions).
    dimensions <- readBin(readBin(files[["png"]], "raw", n = 24L)[17:24], "integer", n = 2L, size = 4L, endian = "big")
    for (kind in absent) {
        asset <- newSavedWorkflowRenderedAsset(readBin(files[[kind]], "raw", n = file.size(files[[kind]])),
            diaCompositeFigureMediaTypes[[kind]], dimensions)
        reference <- saveWorkflowConsumerEntity(store, asset, basename(tempfile("dia-qc-composite-")), identities[[kind]],
            "rendered_asset")
        rootDIAProduct(manager, reference, list(product))
    }
    invisible(NULL)
}

#' Write the composite QC figure's files from the store into a QC directory, under the builder's names
#' (`generateProtNormCompositeQcFigure`), when all three are saved for the stage images that are there.
#' @param manager The session owner.
#' @param product The normalisation run's product reference.
#' @param ruv_mode The run's RUV mode.
#' @param qc_dir The publication's staged QC directory.
#' @param omic_type The omic, which names the files.
#' @return TRUE when the files were written; FALSE when one is not saved (nothing is written).
#' @noRd
writeDIACompositeFigureFiles <- function(manager, product, ruv_mode, qc_dir, omic_type) {
    if (is.null(qc_dir)) return(FALSE)
    identities <- diaCompositeFigureIdentities(manager, product, ruv_mode, qc_dir)
    if (!all(vapply(identities, function(identity) diaCompositeFigureSaved(manager, identity), NA))) return(FALSE)
    for (kind in names(identities)) {
        asset <- resolveSavedWorkflowEntity(manager$.context$store, "rendered_asset", identities[[kind]])$value
        writeBin(asset$bytes, file.path(qc_dir, paste0(omic_type, "_composite_QC_figure.", kind)))
    }
    TRUE
}
