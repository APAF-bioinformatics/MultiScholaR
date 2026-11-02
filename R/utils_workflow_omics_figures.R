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
