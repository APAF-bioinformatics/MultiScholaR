# KIT-077: the lipidomics and metabolomics figures as store entities. Every figure file is written from a saved, rooted
# rendered asset, on the pattern of the proteomics figure products (utils_workflow_dia_products.R), whose functions are
# bound to the protein class and stay untouched. A hit writes the stored bytes; a miss renders, saves and roots. With no
# omics owner in reach (a script, a test without one) the figure is rendered to its file as before.

#' The omics owner a figure is saved on: the innermost scope's (the kit's own code and, from KIT-078, its workers), else
#' the session's for that omic. One session can run several omics modules, so the session keys its owners by omic.
#' @param omic_type `lipidomics` or `metabolomics`.
#' @return An open owner, or NULL.
#' @noRd
workflowOmicsFigureOwner <- function(omic_type) {
    stack <- .diaActiveScope$stack
    session <- shiny::getDefaultReactiveDomain()
    manager <- if (length(stack)) stack[[length(stack)]] else if (!is.null(session)) session$userData$workflow_omics_owners[[omic_type]]
    if (!inherits(manager, "DIAStateAdapter") || !identical(manager$.context$omic_type, omic_type)) return(NULL)
    if (isTRUE(manager$.context$store$closed)) return(NULL)
    manager
}

#' Identify one figure: the omic, what it draws (the assay table and design for the QC figures, the plotted tables
#' otherwise), the assay, the figure and its plotting parameters, the software.
#' @param omic_type The omic.
#' @param content What the figure draws.
#' @param figure The figure's name (its file name without the extension).
#' @param assay The assay's name, or "" for a figure over several.
#' @param parameters The plotting parameters.
#' @param software Installed software identity.
#' @return The store's scientific identity.
#' @noRd
newWorkflowOmicsFigureIdentity <- function(omic_type, content, figure, assay, parameters, software) {
    newSavedWorkflowScientificIdentity(digest::digest(content, algo = "sha256"), "omics-figure-v1",
        digest::digest(parameters, algo = "sha256"),
        digest::digest(list(omic = omic_type, figure = figure, assay = assay), algo = "sha256"), software)
}

#' What a ggplot draws, for a figure's identity: its data and each layer's own data.
#' @noRd
omicsFigurePlotContent <- function(plot) {
    data_of <- function(data) if (is.data.frame(data)) as.data.frame(data)
    list(data = data_of(plot$data), layers = lapply(plot$layers, function(layer) data_of(layer$data)))
}

#' The plot data a rebuildable figure's product keeps (a render's return value), or TRUE when the plot holds none.
#' @noRd
omicsFigurePlotData <- function(plot) {
    if (is.data.frame(plot$data)) as.data.frame(plot$data) else TRUE
}

#' One figure's plot, built on the first file the store does not hold and only once: a repeat export builds no plot.
#' @param build The plot's generator.
#' @param arguments The generator's arguments, taken now.
#' @return A function returning the plot.
#' @noRd
omicsFigureLazyPlot <- function(build, arguments) {
    force(build)
    force(arguments)
    built <- FALSE
    plot <- NULL
    function() {
        if (!built) {
            plot <<- do.call(build, arguments)
            built <<- TRUE
        }
        plot
    }
}

#' What a DA heatmap draws, known before it is built: the contrast's results and, from the source object, the assay and
#' the design it reads.
#' @param results The contrast's long table for this assay.
#' @param object The source object, or anything else where a test mocks the generator.
#' @param assay_slot The object's slot of assay tables.
#' @param assay The assay's name.
#' @return A hash.
#' @noRd
omicsDaHeatmapContent <- function(results, object, assay_slot, assay) {
    source <- if (isS4(object)) list(assay = methods::slot(object, assay_slot)[[assay]], design = object@design_matrix,
        sample_id = object@sample_id, group_id = object@group_id)
    digest::digest(list(results = results, source = source), algo = "sha256")
}

#' The directory the app writes an omic's filtering figures to: its owner's run directory, where a workbook and the
#' proteomics modules write theirs. The modules pass the omic type only.
#' @param omic_type The omic.
#' @return A directory, or NULL with no owner for this omic in reach.
#' @noRd
workflowOmicsFilteringDir <- function(omic_type) {
    if (!is.character(omic_type) || length(omic_type) != 1L) return(NULL)
    time_dir <- workflowOmicsFigureOwner(omic_type)$.context$experiment_paths$time_dir
    if (is.character(time_dir) && length(time_dir) == 1L && nzchar(time_dir)) time_dir
}

#' The composite QC figure (`savePlot` writes its pdf, png and rds in one call) from the store. It is still assembled from
#' the per-assay files; the `.rds` is a plot object, not a figure, and keeps being written beside it.
#' @param plot,base_path,plot_name,width,height What `savePlot` takes.
#' @param plot_files The per-assay files it was assembled from, whose bytes identify it.
#' @param labels Its row and column labels.
#' @param savePlotFn The renderer.
#' @return NULL, invisibly.
#' @noRd
writeWorkflowOmicsCompositeFigure <- function(plot, base_path, plot_name, width, height, plot_files, labels, omic_type, savePlotFn = savePlot) {
    written <- FALSE
    render <- function(target) {
        if (!written) savePlotFn(plot, base_path, plot_name, width = width, height = height, dpi = 150, limitsize = FALSE)
        written <<- TRUE
        TRUE
    }
    sources <- Filter(function(file) is.character(file) && length(file) == 1L && !is.na(file) && file.exists(file), unlist(plot_files))
    content <- unname(vapply(sources, function(file) digest::digest(file = file, algo = "sha256"), character(1)))
    for (format in c("png", "pdf")) {
        writeWorkflowOmicsFigure(file.path(base_path, paste0(plot_name, ".", format)), render, omic_type, content = content, figure = plot_name,
            type = "composite", parameters = labels, width = width, height = height, dpi = 150)
    }
    if (!written) saveRDS(plot, file.path(base_path, paste0(plot_name, ".rds")))
    invisible(NULL)
}

#' Write one figure file from the store.
#' @param path The file, its extension `png` or `pdf`.
#' @param render `function(path)` that writes the file and returns the plot data (a rebuildable type), TRUE (any other
#'   figure), or NULL when there is nothing to draw and no file was written.
#' @param omic_type,content,figure,assay,parameters The figure's identity (`newWorkflowOmicsFigureIdentity`).
#' @param type `pca`, `rle`, `density` or `correlation` keep their plot data in the product (schema
#'   `omics_qc_plot_product_v1`); any other type keeps none (`omics_figure_v1`): grids and multi-page PDFs are bytes only.
#' @param width,height,dpi The rendered size, in inches and dots per inch.
#' @return `path`, or NULL when there was nothing to draw.
#' @noRd
writeWorkflowOmicsFigure <- function(path, render, omic_type, content, figure, assay = "", type = "figure",
    parameters = list(), width = 7, height = 7, dpi = 150) {
    manager <- workflowOmicsFigureOwner(omic_type)
    if (is.null(manager)) return(if (is.null(render(path))) NULL else path)
    format <- tolower(tools::file_ext(path))
    if (!format %in% c("png", "pdf")) cli::cli_abort("an omics figure is a png or a pdf: {.file {path}}", class = "saved_entity_invalid")
    dimensions <- as.integer(ceiling(c(width, height) * dpi))
    store <- manager$.context$store
    software <- manager$.context$software
    identity <- newWorkflowOmicsFigureIdentity(omic_type, content, figure, assay, parameters, software)
    factors <- list(route = omic_type, figure_type = type)
    found <- resolveSavedWorkflowEntity(store, "scientific_product", identity, factors = factors, hydrate = FALSE)
    product <- found$reference
    render_identity <- function(reference) newSavedWorkflowRenderIdentity(reference, list(dpi = dpi), format, dimensions, software)
    if (identical(found$status, "hit")) {
        asset <- resolveSavedWorkflowAsset(store, product, render_identity(product))
        if (identical(asset$status, "hit")) {
            # KIT-078: an earlier worker's figure whose action never completed has no root yet.
            if (is.list(manager$.context$deferred_roots)) {
                rootWorkflowProduct(manager, product, list())
                rootWorkflowProduct(manager, asset$reference, list(product))
            }
            writeBin(asset$value$bytes, path)
            return(path)
        }
    }
    drawn <- render(path)
    if (is.null(drawn) || !file.exists(path)) return(NULL)
    if (is.null(product)) {
        rebuildable <- type %in% c("pca", "rle", "density", "correlation") && is.data.frame(drawn)
        value <- list(schema = if (rebuildable) "omics_qc_plot_product_v1" else "omics_figure_v1", figure = figure, assay = assay,
            figure_type = type, parameters = parameters, data = if (rebuildable) drawn)
        product <- saveWorkflowConsumerEntity(store, newSavedWorkflowScientificProduct(value, "omics-figure",
            provenance = list(route = omic_type, figure = figure)), basename(tempfile("omics-figure-")), identity,
            "scientific_product", factors = factors)
        rootWorkflowProduct(manager, product, list())
    }
    bytes <- readBin(path, "raw", n = file.size(path))
    reference <- saveWorkflowConsumerEntity(store, newSavedWorkflowRenderedAsset(bytes,
        if (identical(format, "png")) "image/png" else "application/pdf", dimensions),
        basename(tempfile("omics-figure-asset-")), render_identity(product), "rendered_asset")
    rootWorkflowProduct(manager, reference, list(product))
    path
}

#' The normalisation tab's composite QC figure (KIT-105 W2, from `writeLipidNormCompositeFigure` and
#' `runMetabNormCompositeQcFigureStep`, which drew the same grid): per assay and kind (pca, density, rle, correlation), a
#' row of the pre-normalisation, post-normalisation and, unless RUV was skipped, RUV-corrected files, assembled from the
#' files only (no matrix). KIT-078: a worker reads the published files and stages the figure in `output_dir`; in process
#' it is written beside them.
#' @param qc_dir The QC directory the per-assay files are in.
#' @param assay_names The assays, in the order of their rows.
#' @param ruv_mode `skip` draws two columns; any other mode three.
#' @param omic_type The omic, which names the figure (`<omic_type>_composite_QC_figure`).
#' @param output_dir Where the figure is written; `qc_dir` when NULL.
#' @param addLog The tab's log.
#' @param generateCompositeFn,savePlotFn,dirExistsFn,logWarnFn The composite's generator, renderer and checks.
#' @return `list(qcDir, ncolComposite, columnLabels, allPlotFiles, allRowLabels, logEntry, compositeSaved)`, with
#'   `errorMessage` before `compositeSaved` on a failure (metab's fields; lipid read `compositeSaved` and `errorMessage`).
#' @noRd
writeWorkflowOmicsNormComposite <- function(qc_dir, assay_names, ruv_mode, omic_type, output_dir = qc_dir,
    addLog = function(entry) invisible(entry), generateCompositeFn = omicsNormCompositeFromFiles, savePlotFn = savePlot,
    dirExistsFn = dir.exists, logWarnFn = logger::log_warn) {
    log_entry <- "Generating composite QC figure..."
    addLog(log_entry)
    empty <- list(qcDir = qc_dir, ncolComposite = NULL, columnLabels = NULL, allPlotFiles = character(0), allRowLabels = list(),
        logEntry = log_entry)
    tryCatch({
        if (is.null(qc_dir) || !dirExistsFn(qc_dir) || is.null(assay_names)) return(c(empty, list(compositeSaved = FALSE)))
        ncol_composite <- if (identical(ruv_mode, "skip")) 2 else 3
        stages <- c("pre_norm", "post_norm", "ruv_corrected")[seq_len(ncol_composite)]
        column_labels <- c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")[seq_len(ncol_composite)]
        all_plot_files <- c()
        all_row_labels <- list()
        label_counter <- 1L
        for (assay_name in assay_names) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            for (plot_type in c("pca", "density", "rle", "correlation")) {
                all_plot_files <- c(all_plot_files, file.path(qc_dir, sprintf("%s_%s_%s.png", safe_name, stages, plot_type)))
                all_row_labels[[sprintf("%s_%s", safe_name, plot_type)]] <- sprintf("%s)", letters[label_counter + seq_len(ncol_composite) - 1L])
                label_counter <- label_counter + ncol_composite
            }
        }
        composite <- generateCompositeFn(all_plot_files, ncol_composite, all_row_labels, column_labels)
        if (!is.null(composite)) {
            writeWorkflowOmicsCompositeFigure(composite$plot, output_dir %||% qc_dir, paste0(omic_type, "_composite_QC_figure"),
                width = composite$width, height = composite$height, plot_files = all_plot_files,
                labels = list(rows = all_row_labels, columns = column_labels), omic_type = omic_type, savePlotFn = savePlotFn)
            addLog(sprintf("Composite QC figure saved to: %s", file.path(qc_dir, "composite_QC_figure")))
        }
        list(qcDir = qc_dir, ncolComposite = ncol_composite, columnLabels = column_labels, allPlotFiles = all_plot_files,
            allRowLabels = all_row_labels, logEntry = log_entry, compositeSaved = !is.null(composite))
    }, error = function(e) {
        message <- conditionMessage(e)
        addLog(paste("Warning: Could not generate composite QC figure:", message))
        logWarnFn(paste("Composite QC generation failed:", message))
        c(empty, list(errorMessage = message, compositeSaved = FALSE))
    })
}

#' The composite from its per-assay files (KIT-105 W2, from lipid's `buildLipidNormCompositeFromFilesGenerator` and
#' metab's `generateMetabNormCompositeFromFiles`, line for line the same): a row of column titles, then per row its
#' labels above its images; a row none of whose files exists is left out.
#' @param plotFiles The files, row by row.
#' @param ncol The columns.
#' @param rowLabels Each row's labels, named by row; `a)`, `b)`, ... when NULL.
#' @param columnLabels The column titles, drawn when there is one per column.
#' @return `list(plot, width, height)` in inches, or NULL (with a warning) when there is nothing to draw.
#' @noRd
omicsNormCompositeFromFiles <- function(plotFiles, ncol = 3, rowLabels = NULL, columnLabels = NULL) {
    text_plot <- function(title, x, size, hjust, bottom, ...) {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = x, y = 0.5, label = title, size = size, ..., hjust = hjust) +
            ggplot2::xlim(0, 1) +
            ggplot2::theme_void() +
            ggplot2::theme(plot.margin = ggplot2::margin(5, 5, bottom, 5), panel.background = ggplot2::element_blank())
    }
    label_plot <- function(title) text_plot(title, 0, 5, 0, 5)
    title_plot <- function(title) text_plot(title, 0.5, 6, 0.5, 10, fontface = "bold")
    image_plot <- function(path) {
        if (is.na(path) || !file.exists(path)) return(ggplot2::ggplot() + ggplot2::theme_void())
        tryCatch({
            # KIT-079: a native raster, 4 bytes a pixel against 32 for an array of doubles; the rendered figure is byte-identical.
            grob <- grid::rasterGrob(png::readPNG(path, native = TRUE), interpolate = TRUE)
            ggplot2::ggplot() + ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + ggplot2::theme_void()
        }, error = function(e) {
            logger::log_warn(sprintf("[generateCompositeFromFiles] Could not load image: %s", path))
            ggplot2::ggplot() + ggplot2::theme_void()
        })
    }
    logger::log_info(sprintf("[generateCompositeFromFiles] Generating composite from %d files...", length(plotFiles)))
    for (package in c("patchwork", "ggplot2", "png")) {
        if (!requireNamespace(package, quietly = TRUE)) {
            warning(paste(package, "package required for composite generation"))
            return(NULL)
        }
    }
    tryCatch({
        n_files <- length(plotFiles)
        if (is.null(rowLabels)) {
            n_rows <- n_files / ncol
            rowLabels <- split(paste0(letters[seq_len(n_files)], ")"), rep(seq_len(n_rows), each = ncol))
            names(rowLabels) <- paste0("row", seq_len(n_rows))
        }
        plot_sections <- list()
        height_values <- c()
        if (!is.null(columnLabels) && length(columnLabels) == ncol) {
            plot_sections <- append(plot_sections, list(patchwork::wrap_plots(lapply(columnLabels, title_plot), ncol = ncol)))
            height_values <- c(height_values, 0.2)
            logger::log_info("[generateCompositeFromFiles] Added column titles")
        }
        row_names <- names(rowLabels)
        for (i in seq_along(row_names)) {
            row_files <- plotFiles[((i - 1) * ncol + 1):min(i * ncol, n_files)]
            if (any(!is.na(row_files) & vapply(row_files, function(path) !is.na(path) && file.exists(path), logical(1)))) {
                plot_sections <- append(plot_sections, list(
                    patchwork::wrap_plots(lapply(rowLabels[[row_names[i]]], label_plot), ncol = ncol),
                    patchwork::wrap_plots(lapply(row_files, image_plot), ncol = ncol)))
                height_values <- c(height_values, 0.1, 1)
                logger::log_info(sprintf("[generateCompositeFromFiles] Added row: %s", row_names[i]))
            } else {
                logger::log_info(sprintf("[generateCompositeFromFiles] Skipping empty row: %s", row_names[i]))
            }
        }
        if (length(plot_sections) == 0) {
            warning("No valid plot sections to combine")
            return(NULL)
        }
        logger::log_info("[generateCompositeFromFiles] Combining plot sections...")
        combined_plot <- patchwork::wrap_plots(plot_sections, ncol = 1) + patchwork::plot_layout(heights = height_values)
        plot_width <- 4 + (ncol * 3)
        plot_height <- 4 + (length(height_values) * 2)
        rm(plot_sections)
        gc()
        list(plot = combined_plot, width = plot_width, height = plot_height)
    }, error = function(e) {
        logger::log_error(paste("[generateCompositeFromFiles] Error:", e$message))
        NULL
    })
}
