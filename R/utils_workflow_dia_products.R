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
    breaks <- seq(min(round(data$pearson_correlation - 0.5, 2), na.rm = TRUE),
      1, 0.001)
    return(ggplot2::ggplot(data, ggplot2::aes(.data[["pearson_correlation"]])) +
      ggplot2::geom_histogram(breaks = breaks) +
      ggplot2::scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
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

#' Root a reusable scientific or display product in the durable session head.
#' @param manager Current component owner.
#' @param reference Accepted existing-store product descriptor.
#' @return Descriptor invisibly, without hydrating its source components.
#' @noRd
rootWorkflowProduct <- function(manager, reference, dependencies = list()) {
    snapshot <- manager$.context$snapshot
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
    publishDIAStateHead(manager$.context, snapshot)
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
#' @noRd
diaNormalizationOperation <- function() {
    newWorkflowOperation("normalization-ruv", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        output_class = "ProteinQuantitativeData", product_type = "dia-normalization-ruv",
        parameters = c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
            "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k", "grouping"),
        kernel = function(object, parameters, store) {
            # The kernels record their arguments under the calling name; the legacy path calls them through
            # these aliases, so the aliases keep the recorded args (and every checkpoint hash) identical.
            normaliseBetweenSamplesFn <- normaliseBetweenSamples
            normalized <- normaliseBetweenSamplesFn(object, normalisation_method = parameters$norm_method)
            normalized_reference <- store$checkpoint(normalized, "normalized-science")
            object <- NULL
            if (parameters$ruv_mode == "skip") {
                ruv_result <- buildProtNormSkippedRuvResult()
                final_reference <- normalized_reference
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
                corrected <- ruvIII_C_VaryingFn(normalized,
                    ruv_grouping_variable = parameters$grouping, ruv_number_k = resolved$k, ctrl = resolved$controls)
                corrected <- removeRowsWithMissingValuesPercentFn(theObject = corrected)
                final_reference <- store$checkpoint(corrected, "ruv-science", from = normalized_reference)
                corrected <- NULL
                best_k <- resolved$k
                controls <- resolved$controls
                best_percentage <- resolved$percentage
            }
            normalized <- NULL
            list(normalized = normalized_reference, final = final_reference, ruv_result = ruv_result,
                best_k = best_k, controls = controls, percentage = best_percentage,
                parameters = parameters[c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
                    "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k", "grouping")])
        })
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

runDIAReferenceNormalization <- function(settings, workflow, norm_data, paths, omic_type,
    experiment_label, post_qc, ruv_qc, grouping) {
    prior_plot <- ggplot2::get_last_plot()
    on.exit(ggplot2::set_last_plot(prior_plot), add = TRUE)
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "normalization")
        source <- norm_data$normalization_input_reference %||%
            manager$.context$snapshot$qc_sources$pre_normalization %||%
            getDIAStateReference(manager, "protein_replicate_filtered") %||% getDIAStateReference(manager)
        norm_data$normalization_input_reference <- source
        fields <- c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
            "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k")
        parameters <- stats::setNames(lapply(fields, \(field) settings[[field]]), fields)
        parameters$grouping <- grouping()
        run <- runSavedOperation(manager, diaNormalizationOperation(), source, parameters)
        saved <- run$value
        product <- run$reference
        norm_data$normalized_protein_obj <- saved$normalized
        norm_data$ruv_normalized_obj <- saved$final
        norm_data$best_k <- saved$best_k
        norm_data$control_genes_index <- saved$controls
        norm_data$ruv_optimization_result <- saved$ruv_result
        workflow$ruv_optimization_result <- saved$ruv_result
        # KIT-018: the optimisation grid is a saved display table the module pages on demand.
        norm_data$ruv_optimization_table <- saveDIADisplayTable(manager, saved$ruv_result$optimization_results)
        skipped <- identical(settings$ruv_mode, "skip")
        final <- readDIAStateMetadata(manager, saved$final)$data
        config <- list(norm_method = settings$norm_method, ruv_mode = settings$ruv_mode,
            ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage)
        if (skipped) config$ruv_applied <- FALSE
        state_name <- if (skipped) "normalized" else "ruv_corrected"
        final$products$normalization_science <- product
        previous <- if (manager$current_state == state_name) readDIAStateMetadata(manager) else NULL
        if (is.null(previous) || !identical(previous$data$products$normalization_science, product)) {
            commitDIAReferenceState(manager, state_name, final, config,
                if (skipped) "Post-normalization complete: RUV-III skipped by user" else
                    "Post-normalization complete: RUV-III correction and missing value cleanup completed")
        }
        accepted <- getDIAStateReference(manager)
        display <- publishDIAActionOutputs(manager, accepted, "normalization",
            paths[c("source_dir", "protein_qc_dir", "time_dir")], function(staged) {
                norm_data$publication_paths <- staged
                on.exit({ norm_data$publication_paths <- NULL }, add = TRUE)
                pre_files <- list.files(paths$protein_qc_dir, pattern = "^pre_norm_.*[.]png$", full.names = TRUE)
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
                    if (!is.null(saved$ruv_result$best_cancor_plot)) ggplot2::ggsave(
                        file.path(staged$protein_qc_dir, "ruv_corrected_cancor.png"),
                        saved$ruv_result$best_cancor_plot, width = 8, height = 6, dpi = 150)
                }
                generateProtNormCompositeQcFigure(settings$ruv_mode, staged$protein_qc_dir, omic_type)
                list(qc_paths = norm_data$qc_plot_paths, filtering_plot = norm_data$post_norm_filtering_plot,
                    filtering_text = norm_data$filtering_summary_text, staged_qc_dir = staged$protein_qc_dir)
            }, identity = list(product = product$entity_id, color = settings$color_variable,
                shape = settings$shape_variable, output_schema = 2L))
        norm_data$qc_plot_paths <- lapply(display$qc_paths, function(group) lapply(group, function(path) {
            if (is.character(path)) sub(display$staged_qc_dir, paths$protein_qc_dir, path, fixed = TRUE) else path
        }))
        norm_data$post_norm_filtering_plot <- display$filtering_plot
        norm_data$filtering_summary_text <- display$filtering_text
        norm_data$normalization_complete <- TRUE
        finalizeProtNormWorkflowState(norm_data)
        TRUE
    })
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
#' @noRd
diaRuvCovariatesOperation <- function() {
    newWorkflowOperation("protein_ruv_covariates", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        output_class = "ProteinQuantitativeData", product_type = "dia-ruv-covariates",
        parameters = c("percentage_min", "percentage_max", "num_components_to_impute", "ruv_grouping_variable",
            "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_mode"),
        kernel = function(object, parameters, store) {
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
#' @noRd
resolveDIARuvCovariatesParameters <- function(config_list, grouping = NULL, ruv_mode = "covariates") {
    if (!is.character(ruv_mode) || length(ruv_mode) != 1L || !ruv_mode %in% c("covariates", "correction"))
        cli::cli_abort("ruv_mode must be covariates or correction", class = "workflow_execution_invalid")
    section <- config_list$ruvCovariates %||% list()
    shipped <- readWorkflowShippedDefaults()$ruvCovariates %||% list()
    pick <- function(name, default) {
        value <- section[[name]] %||% shipped[[name]] %||% default
        value
    }
    numeric <- function(name, default) as.numeric(pick(name, default))
    list(percentage_min = numeric("percentage_min", 1), percentage_max = numeric("percentage_max", 30),
        num_components_to_impute = numeric("num_components_to_impute", 5),
        ruv_grouping_variable = grouping %||% as.character(pick("ruv_grouping_variable", "group")),
        ruv_qval_cutoff = numeric("ruv_qval_cutoff", 0.05), ruv_fdr_method = as.character(pick("ruv_fdr_method", "BH")),
        separation_metric = as.character(pick("separation_metric", "max_difference")),
        k_penalty_weight = numeric("k_penalty_weight", 0.5),
        adaptive_k_penalty = tolower(as.character(pick("adaptive_k_penalty", "true"))) %in% c("true", "1", "yes"),
        ruv_mode = ruv_mode)
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
        source <- norm_data$correlation_filtered_obj %||% getDIAStateReference(manager)
        if (!isSavedWorkflowReference(source)) cli::cli_abort("the correlation step must complete before the RUV covariates")
        parameters <- resolveDIARuvCovariatesParameters(workflow$config_list, grouping(), settings$ruv_mode %||% "covariates")
        run <- runSavedOperation(manager, diaRuvCovariatesOperation(), source, parameters)
        saved <- run$value
        product <- run$reference
        norm_data$ruv_optimization_result <- saved$ruv_result
        workflow$ruv_optimization_result <- saved$ruv_result
        norm_data$ruv_optimization_table <- saveDIADisplayTable(manager, saved$ruv_result$optimization_results)
        norm_data$best_k <- saved$best_k
        norm_data$control_genes_index <- saved$controls
        ruv_mode <- saved$ruv_mode %||% "covariates"
        correction <- identical(ruv_mode, "correction")
        state_name <- if (correction) "ruv_corrected" else "ruv_covariates"
        final <- readDIAStateMetadata(manager, saved$final %||% saved$covariates)$data
        final$products$ruv_covariates_science <- product
        config <- list(ruv_mode = ruv_mode, ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage,
            control_count = saved$control_count, separation_score = saved$ruv_result$best_separation_score,
            composite_score = saved$ruv_result$best_composite_score, max_acceptable_k = saved$ruv_result$max_acceptable_k,
            replicate_structure = saved$replicate_structure, block_reason = saved$block_reason, replicate_blocks = saved$blocks,
            ruv_grouping_variable = parameters$ruv_grouping_variable, ruv_applied = TRUE,
            outcomes = list(k = saved$best_k, percentage_as_neg_ctrl = saved$percentage, control_count = saved$control_count,
                separation_score = saved$ruv_result$best_separation_score, composite_score = saved$ruv_result$best_composite_score,
                replicate_structure = saved$replicate_structure, block_reason = saved$block_reason, ruv_mode = ruv_mode))
        previous <- if (identical(manager$current_state, state_name)) readDIAStateMetadata(manager) else NULL
        if (is.null(previous) || !identical(previous$data$products$ruv_covariates_science, product)) {
            commitDIAReferenceState(manager, state_name, final, config, if (correction)
                sprintf("RUV-III correction applied to the quantified matrix and the DPC EList (k = %d; Y - W alpha); no W terms in the design", saved$best_k) else
                sprintf("RUV-III factors appended to the design as covariates (W1..W%d); the quantified matrix unchanged", saved$best_k))
        }
        accepted <- getDIAStateReference(manager)
        norm_data$ruv_normalized_obj <- accepted
        norm_data$correlation_filtered_obj <- accepted
        workflow$ruv_normalised_for_da_analysis_obj <- accepted
        display <- publishDIAActionOutputs(manager, accepted, "ruv-covariates", paths[c("source_dir", "protein_qc_dir")],
            function(staged) {
                norm_data$publication_paths <- staged
                on.exit({ norm_data$publication_paths <- NULL }, add = TRUE)
                writeDIADesignTable(manager, readDIAStateMetadata(manager, saved$covariates)$data$components$design_matrix,
                    file.path(staged$protein_qc_dir, "ruv_covariates_design.tsv"), vroom::vroom_write)
                adjusted_spec <- readDIAStateMetadata(manager, saved$adjusted)$data
                writeDIADesignTable(manager, adjusted_spec$components$protein_quant_table,
                    file.path(staged$protein_qc_dir, "ruv_adjusted_protein_matrix.tsv"), vroom::vroom_write)
                vroom::vroom_write(saved$ruv_result$optimization_results, file.path(staged$protein_qc_dir, "ruv_covariates_optimization.tsv"))
                saveRDS(saved$ruv_result[setdiff(names(saved$ruv_result), "best_cancor_plot")],
                    file.path(staged$source_dir, "ruv_covariates_optimization_results.RDS"))
                if (!is.null(saved$ruv_result$best_cancor_plot)) ggplot2::ggsave(
                    file.path(staged$protein_qc_dir, "ruv_covariates_cancor.png"),
                    saved$ruv_result$best_cancor_plot, width = 8, height = 6, dpi = 150)
                ruv_qc(saved$adjusted)
                list(completed = TRUE, qc_paths = norm_data$qc_plot_paths)
            }, identity = list(product = product$entity_id, color = settings$color_variable,
                shape = settings$shape_variable, output_schema = 1L))
        norm_data$ruv_covariates_complete <- TRUE
        accepted
    })
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
        source <- norm_data$ruv_normalized_obj
        if (!isSavedWorkflowReference(source)) cli::cli_abort("normalization must complete before correlation filtering")
        original <- readDIAStateMetadata(manager, source)$data
        if (skipped) {
            final <- source
        } else {
            pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), source,
                list(group = grouping(), tech_rep_remove_regex = "pool"))
            norm_data$correlation_vector <- pairs$reference
            selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), source,
                list(pairs = pairs$reference$entity_id, threshold = threshold, correlations = pairs$value),
                dependencies = list(pairs$reference))
            final <- selection$value
            pairs <- selection <- NULL
        }
        norm_data$correlation_filtered_obj <- final
        norm_data$correlation_threshold <- if (skipped) 0 else threshold
        spec <- readDIAStateMetadata(manager, final)$data
        config <- if (skipped) list(min_pearson_correlation_threshold = 0, skipped = TRUE) else
            list(min_pearson_correlation_threshold = threshold)
        previous <- if (manager$current_state == "correlation_filtered") readDIAStateMetadata(manager) else NULL
        if (is.null(previous) || !identical(previous$data$components, spec$components) ||
            !identical(previous$config, config)) commitDIAReferenceState(manager,
                "correlation_filtered", spec, config, if (skipped) "Skipped final sample correlation filter" else
                    "Applied final sample correlation filter (chunk 28)")
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
    })
}
