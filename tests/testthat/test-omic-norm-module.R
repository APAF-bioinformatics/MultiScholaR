# KIT-105 W1: the normalisation tab through its public entry points, written before the shared
# normalisation core (W2) and the shared frame and panel library (R/mod_omic_norm.R,
# R/mod_omic_norm_panels.R; W4, W5) exist. Lipid and metab are driven through
# `mod_{lipid,metab}_norm_server` in `testServer`, in process (no runner), on the MCI-021/MCI-015
# fixtures, so the blocks hold through phases B to D.
#   (a) Identity: the norm UIs' HTML digests and their namespaced input ids, output ids, test ids and
#       tab titles (lipid, metab, and trx's mount of the lipid UI under its profile).
#   (b) Server captures: the run under RUV skip (ITSD off), manual (a manual ITSD selection) and
#       automatic (ITSD by the state's pattern); correlation apply and skip; reset; the export (and
#       its refusal before the run); the final QC plot; the tab entry's pre-normalisation QC.
#   (c) The twin differences (design D7): skip-correlation when the store commit fails, the final
#       PCA failing, and the composite with a missing figure file, recorded per omic.
#   F105-01 (gated): each QC image slot and assay label shows its own assay, stage and kind.
# (a) to (c) were frozen on the start commit (9a98e765) and pass there
# (docs/qa/workflow-execution/kit/kit-105/w1/test-first.txt). The RUV kernels are the MCI lanes'
# stubs; every other kernel is real. MULTISCHOLAR_KIT105_W1_DUMP=<file> writes the captures as R
# code (how the frozen values were made).

#' Skip a block that fails until a later step lands it, unless MULTISCHOLAR_KIT105_PENDING=1.
omicNormPending <- function(phase) {
    if (!identical(Sys.getenv("MULTISCHOLAR_KIT105_PENDING"), "1")) {
        testthat::skip(sprintf("KIT-105 W1: expected to fail until %s lands it (MULTISCHOLAR_KIT105_PENDING=1 runs it)",
            phase))
    }
}

#' The rendered HTML, with shiny's random tabset ids masked (as the entry suite does, without angle
#' brackets, so that the tags stay readable).
omicNormHtml <- function(x) {
    html <- if (inherits(x, c("shiny.tag", "shiny.tag.list"))) htmltools::renderTags(x)$html else x$html
    html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"ID\"", as.character(html))
    gsub("tab-[0-9]+-([0-9]+)", "tab-ID-\\1", html)
}

#' A value's sha256 (version-3 serialisation).
omicNormSha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

#' Where MULTISCHOLAR_KIT105_W1_DUMP points, or NULL.
omicNormDumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT105_W1_DUMP")
    if (nzchar(path)) path else NULL
}

#' R code for a capture (lists, named vectors and scalars): a list's elements one per line, a long
#' vector's packed into lines, indented.
omicNormCode <- function(x, indent = 0L) {
    pad <- strrep(" ", indent + 4L)
    label <- function(names, i) {
        if (is.null(names) || !nzchar(names[[i]])) return("")
        paste0(if (make.names(names[[i]]) == names[[i]]) names[[i]] else paste0("`", names[[i]], "`"), " = ")
    }
    if (is.list(x)) {
        if (length(x) == 0L) return("list()")
        items <- vapply(seq_along(x), function(i) paste0(pad, label(names(x), i), omicNormCode(x[[i]], indent + 4L)), "")
        return(paste0("list(\n", paste(items, collapse = ",\n"), ")"))
    }
    if (is.null(x)) return("NULL")
    scalar <- function(value) deparse(value, control = c("keepNA", "keepInteger"))
    if (length(x) == 0L) return(if (is.null(names(x))) scalar(x) else sprintf("stats::setNames(%s, character())", scalar(unname(x))))
    if (length(x) == 1L && is.null(names(x))) return(scalar(x))
    items <- paste0(vapply(seq_along(x), function(i) paste0(label(names(x), i), scalar(unname(x[[i]]))), ""), ",")
    items[[length(items)]] <- sub(",$", ")", items[[length(items)]])
    if (sum(nchar(items) + 1L) + indent < 90L) return(paste0("c(", paste(items, collapse = " ")))
    # Items packed into lines of about 110 characters.
    lines <- character()
    for (item in items) {
        last <- length(lines)
        if (last && nchar(lines[[last]]) + nchar(item) < 110L) {
            lines[[last]] <- paste(lines[[last]], item)
        } else {
            lines <- c(lines, paste0(pad, item))
        }
    }
    paste0("c(\n", paste(lines, collapse = "\n"))
}

#' Append `name <- <value>` to the dump file.
omicNormDump <- function(name, value) {
    cat(sprintf("%s <- %s\n\n", name, omicNormCode(value)), file = omicNormDumpPath(), append = TRUE)
}

# ------------------------------------------------------------------------------------------------
# (a) Identity: the norm UIs, frozen at 9a98e765.
# ------------------------------------------------------------------------------------------------

#' One norm UI's HTML digest and its ids: inputs (every namespaced `id` that is not an output, a
#' label or a tab pane), outputs (`shiny-*-output` and widget outputs), `data-testid` values and the
#' tab titles in order.
omicNormUiParts <- function(ui) {
    html <- omicNormHtml(ui)
    tags <- regmatches(html, gregexpr("<[a-zA-Z][^>]*>", html))[[1L]]
    attr_of <- function(name) vapply(tags, function(tag) {
        m <- regmatches(tag, regexec(sprintf("[[:space:]]%s=\"([^\"]*)\"", name), tag))[[1L]]
        if (length(m)) m[[2L]] else NA_character_
    }, "", USE.NAMES = FALSE)
    ids <- attr_of("id")
    classes <- attr_of("class")
    is_output <- !is.na(ids) & grepl("(^| )(shiny-[a-z]+-output|html-widget-output)( |$)", classes)
    inputs <- ids[!is.na(ids) & !is_output & !grepl("-label$", ids) & !startsWith(ids, "tab-ID-")]
    list(html = digest::digest(html, algo = "sha256", serialize = FALSE), inputs = sort(unique(inputs)),
        outputs = sort(unique(ids[is_output])), testids = sort(unique(stats::na.omit(attr_of("data-testid")))),
        tabs = attr_of("data-value")[grepl("data-toggle=\"tab\"", tags)])
}

# The frozen values between the markers are the dump's (MULTISCHOLAR_KIT105_W1_DUMP), pasted.
# OMIC_NORM_FROZEN_UI_START
omicNormFrozenUi <- list(
    lipid = list(
        html = "d27fa2638c9a5943f554fe47eab4e1b8b04c247b0386d1114d78c895c95e4725",
        inputs = c(
            "norm-adaptive_k_penalty", "norm-apply_correlation_filter", "norm-apply_itsd",
            "norm-auto_percentage_max", "norm-auto_percentage_min", "norm-color_variable",
            "norm-export_session", "norm-itsd_aggregation", "norm-k_penalty_weight", "norm-log_offset",
            "norm-min_pearson_correlation_threshold", "norm-norm_method", "norm-norm_qc_tabs",
            "norm-reset_normalization", "norm-run_normalization", "norm-ruv_grouping_variable", "norm-ruv_k",
            "norm-ruv_mode", "norm-ruv_percentage", "norm-separation_metric", "norm-shape_variable",
            "norm-skip_correlation_filter"),
        outputs = c(
            "norm-assay1_label_correlation", "norm-assay1_label_density", "norm-assay1_label_pca",
            "norm-assay1_label_rle", "norm-assay2_label_correlation", "norm-assay2_label_density",
            "norm-assay2_label_pca", "norm-assay2_label_rle", "norm-correlation_filter_summary",
            "norm-correlation_post_filter_assay1", "norm-correlation_post_filter_assay2",
            "norm-correlation_post_norm_assay1", "norm-correlation_post_norm_assay2",
            "norm-correlation_ruv_corrected_assay1", "norm-correlation_ruv_corrected_assay2",
            "norm-density_post_filter_assay1", "norm-density_post_filter_assay2",
            "norm-density_post_norm_assay1", "norm-density_post_norm_assay2",
            "norm-density_ruv_corrected_assay1", "norm-density_ruv_corrected_assay2", "norm-final_qc_plot",
            "norm-itsd_selection_ui", "norm-norm_log", "norm-pca_post_filter_assay1",
            "norm-pca_post_filter_assay2", "norm-pca_post_norm_assay1", "norm-pca_post_norm_assay2",
            "norm-pca_ruv_corrected_assay1", "norm-pca_ruv_corrected_assay2", "norm-rle_post_filter_assay1",
            "norm-rle_post_filter_assay2", "norm-rle_post_norm_assay1", "norm-rle_post_norm_assay2",
            "norm-rle_ruv_corrected_assay1", "norm-rle_ruv_corrected_assay2", "norm-ruv_qc_ui"),
        testids = c(
            "lipid-norm-export-session", "lipid-norm-method", "lipid-norm-run", "lipid-norm-ruv-mode"),
        tabs = c(
            "ITSD Selection", "RUV QC", "PCA", "Density", "RLE", "Correlation", "Correlation Filtering",
            "Log")),
    metab = list(
        html = "51dc7b957cf66b34353b73686b6346e9d87f2392845ef5f15ca115465f95c427",
        inputs = c(
            "norm-adaptive_k_penalty", "norm-apply_correlation_filter", "norm-apply_itsd",
            "norm-auto_percentage_max", "norm-auto_percentage_min", "norm-color_variable",
            "norm-export_session", "norm-itsd_aggregation", "norm-k_penalty_weight", "norm-log_offset",
            "norm-min_pearson_correlation_threshold", "norm-norm_method", "norm-norm_qc_tabs",
            "norm-reset_normalization", "norm-run_normalization", "norm-ruv_grouping_variable", "norm-ruv_k",
            "norm-ruv_mode", "norm-ruv_percentage", "norm-separation_metric", "norm-shape_variable",
            "norm-skip_correlation_filter"),
        outputs = c(
            "norm-assay1_label_correlation", "norm-assay1_label_density", "norm-assay1_label_pca",
            "norm-assay1_label_rle", "norm-assay2_label_correlation", "norm-assay2_label_density",
            "norm-assay2_label_pca", "norm-assay2_label_rle", "norm-correlation_filter_summary",
            "norm-correlation_post_filter_assay1", "norm-correlation_post_filter_assay2",
            "norm-correlation_post_norm_assay1", "norm-correlation_post_norm_assay2",
            "norm-correlation_ruv_corrected_assay1", "norm-correlation_ruv_corrected_assay2",
            "norm-density_post_filter_assay1", "norm-density_post_filter_assay2",
            "norm-density_post_norm_assay1", "norm-density_post_norm_assay2",
            "norm-density_ruv_corrected_assay1", "norm-density_ruv_corrected_assay2", "norm-final_qc_plot",
            "norm-itsd_selection_ui", "norm-norm_log", "norm-pca_post_filter_assay1",
            "norm-pca_post_filter_assay2", "norm-pca_post_norm_assay1", "norm-pca_post_norm_assay2",
            "norm-pca_ruv_corrected_assay1", "norm-pca_ruv_corrected_assay2", "norm-rle_post_filter_assay1",
            "norm-rle_post_filter_assay2", "norm-rle_post_norm_assay1", "norm-rle_post_norm_assay2",
            "norm-rle_ruv_corrected_assay1", "norm-rle_ruv_corrected_assay2", "norm-ruv_qc_ui"),
        testids = c(
            "metab-norm-export-session", "metab-norm-method", "metab-norm-run", "metab-norm-ruv-mode"),
        tabs = c(
            "ITSD Selection", "RUV QC", "PCA", "Density", "RLE", "Correlation", "Correlation Filtering",
            "Log")),
    trx = list(
        html = "2548b720447ed36de5639d3a2dad4b0d5a58744bc78a3700f787737b9718cefa",
        inputs = c(
            "norm-auto_percentage_max", "norm-auto_percentage_min", "norm-auto_percentage_step",
            "norm-k_penalty_weight", "norm-norm_qc_tabs", "norm-reset_normalization",
            "norm-run_normalization", "norm-ruv_k", "norm-ruv_mode", "norm-ruv_percentage"),
        outputs = c(
            "norm-assay1_label_pca", "norm-assay1_label_rle", "norm-norm_log", "norm-pca_post_filter_assay1",
            "norm-pca_post_norm_assay1", "norm-pca_ruv_corrected_assay1", "norm-rle_post_filter_assay1",
            "norm-rle_post_norm_assay1", "norm-rle_ruv_corrected_assay1", "norm-ruv_qc_ui"),
        testids = c("trx-norm-reset", "trx-norm-run", "trx-norm-ruv-mode"),
        tabs = c("RUV QC", "PCA", "RLE", "Log")))
# OMIC_NORM_FROZEN_UI_END

test_that("(a) the lipid, metab and trx-profile norm UIs keep their HTML, ids, test ids and tab titles", {
    actual <- list(lipid = omicNormUiParts(mod_lipid_norm_ui("norm")), metab = omicNormUiParts(mod_metab_norm_ui("norm")),
        trx = omicNormUiParts(mod_lipid_norm_ui("norm", profile = trxOmicProfile())))
    if (!is.null(omicNormDumpPath())) omicNormDump("omicNormFrozenUi", actual)
    expect_named(actual, names(omicNormFrozenUi))
    for (omic in names(omicNormFrozenUi)) {
        for (part in names(omicNormFrozenUi[[omic]])) {
            expect_identical(actual[[omic]][[part]], omicNormFrozenUi[[omic]][[part]], info = paste(omic, part))
        }
    }
})

# ------------------------------------------------------------------------------------------------
# (b) and (c): the server, driven in process.
# ------------------------------------------------------------------------------------------------

#' Each omic's entry points, fixture, QC-directory key and the kernels the lanes stub or record.
omicNormSpec <- function(omic) {
    switch(omic,
        lipid = list(omic_type = "lipidomics", server = mod_lipid_norm_server, ui = mod_lipid_norm_ui,
            object = module_ci_lipid_norm_object, workflow = module_ci_lipid_norm_workflow,
            ruv_results = module_ci_lipid_norm_ruv_results, qc_key = "lipid_qc_dir",
            ruv_kernel = "runLipidPerAssayRuvOptimization", qc_kernel = "generateLipidQcPlots"),
        metab = list(omic_type = "metabolomics", server = mod_metab_norm_server, ui = mod_metab_norm_ui,
            object = module_ci_metab_norm_object, workflow = module_ci_metab_norm_workflow,
            ruv_results = module_ci_metab_norm_ruv_results, qc_key = "metabolite_qc_dir",
            ruv_kernel = "runPerAssayRuvOptimization", qc_kernel = "generateMetabQcPlots"))
}

# The 24 image slots and 8 assay labels of the stage grid, and the file stage each slot shows.
omicNormImageIds <- as.vector(t(outer(c("pca", "density", "rle", "correlation"),
    c("post_filter_assay1", "post_norm_assay1", "ruv_corrected_assay1", "post_filter_assay2", "post_norm_assay2",
        "ruv_corrected_assay2"), paste, sep = "_")))
omicNormLabelIds <- as.vector(outer(1:2, c("pca", "density", "rle", "correlation"),
    function(n, kind) sprintf("assay%d_label_%s", n, kind)))

#' What one image slot shows: the QC files whose bytes its `src` holds (several when figures are
#' byte-identical) with its alt text, or the alt text alone when it has no image.
omicNormImageShown <- function(value, encoded) {
    if (!is.list(value)) return(paste("<value>", paste(value, collapse = " ")))
    src <- value$src %||% ""
    if (!nzchar(src)) return(paste("<none>", value$alt %||% ""))
    shown <- names(encoded)[encoded == sub("^data:[^,]*,", "", src)]
    paste0(if (length(shown)) paste(shown, collapse = "|") else "<unknown>", " alt=", value$alt %||% "")
}

#' The `normData` flags, assay names, ITSD selections and each state slot's committed state.
omicNormFlags <- function(nd, manager) {
    flags <- c("normalization_complete", "ruv_complete", "correlation_filtering_complete", "pre_norm_qc_generated")
    slots <- c("post_filter_obj", "post_itsd_obj", "post_log2_obj", "post_norm_obj", "ruv_corrected_obj",
        "correlation_filtered_obj")
    c(vapply(stats::setNames(flags, flags), function(flag) as.character(isTRUE(nd[[flag]])), ""),
        assay_names = paste(nd$assay_names, collapse = "|"),
        itsd_selections = paste(names(nd$itsd_selections), vapply(nd$itsd_selections, paste, "", collapse = ","),
            sep = "=", collapse = ";"),
        vapply(stats::setNames(slots, slots), function(slot) omicNormSlotState(nd[[slot]], manager), ""))
}

#' The step provenance entries without their build identity and time; a stored product's entity id
#' (the correlation pairs') differs between processes, so it is read as `<ENTITY>`.
omicNormProvenance <- function(manager) {
    lapply(manager$.context$snapshot$step_provenance, function(entry) {
        entry <- entry[intersect(c("operation", "variant", "parameter_names", "fixed", "values", "preceded_by"), names(entry))]
        rapply(entry, function(value) sub("^[0-9a-f]{64}$", "<ENTITY>", value), classes = "character", how = "replace")
    })
}

#' The committed state (or states) a `normData` slot's reference names; an object by its class.
omicNormSlotState <- function(value, manager) {
    if (is.null(value)) return("<NULL>")
    if (!isSavedWorkflowReference(value)) return(paste0("<", class(value)[[1L]], ">"))
    states <- manager$.context$snapshot$states
    named <- names(Filter(function(ref) identical(ref$entity_id, value$entity_id), states))
    if (length(named)) paste(named, collapse = "|") else "<unlisted reference>"
}

#' Run one omic's norm tab through its public server under `steps` (named functions of the context)
#' and capture after each step. The bus is the MCI fixture's, with a real owner; the tab starts on
#' `qc`. Returns `captures` (by step), `extra` (what steps kept) and `choices` (the design-driven
#' select choices the server sent).
omicNormDrive <- function(omic, steps) {
    spec <- omicNormSpec(omic)
    root <- normalizePath(withr::local_tempdir(pattern = paste0("kit105-norm-", omic, "-")))
    source_dir <- file.path(root, "source")
    dir.create(source_dir)
    paths <- stats::setNames(list(file.path(root, "qc"), source_dir, file.path(root, "export")),
        c(spec$qc_key, "source_dir", "export_dir"))
    workflow <- spec$workflow(spec$object(positive_only = TRUE), experiment_paths = list(source_dir = source_dir))
    log <- new.env(parent = emptyenv())
    for (field in c("updates", "notes", "choices")) log[[field]] <- list()
    log$modals <- integer()
    guiEntryInstallBrowserDouble(log)
    withr::local_seed(20261108L)
    acc <- new.env(parent = emptyenv())
    for (field in c("captures", "slots", "extra", "pca", "labels", "previous")) acc[[field]] <- list()
    acc$log <- 0L
    # The MCI lanes' RUV stubs (test-module-ci-{lipid,metab}-norm.R): the optimisation's per-assay
    # results, and a correction that records its arguments.
    ruv_results <- spec$ruv_results
    stubs <- list(function(theObject, ruv_mode, params, experiment_paths) ruv_results(theObject, mode = ruv_mode),
        function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
        theObject@args$ruvIII_C_Varying <- list(ruv_grouping_variable = ruv_grouping_variable,
            ruv_number_k = ruv_number_k, ctrl = ctrl)
        theObject
    })
    names(stubs) <- c(spec$ruv_kernel, "ruvIII_C_Varying")
    # The PCA kernel, called through: the final QC plot's calls are recorded (and refused when
    # `acc$fail` is set).
    plot_pca <- plotPca
    stubs$plotPca <- function(theObject, ...) {
        args <- list(...)
        final <- identical(args$title, "Final QC - PCA")
        if (final && !is.null(acc$fail)) stop(acc$fail)
        plots <- plot_pca(theObject, ...)
        if (final) {
            acc$pca[[length(acc$pca) + 1L]] <- list(class = class(theObject)[[1L]], state = guiEntryDigest(theObject),
                grouping_variable = args$grouping_variable, shape_variable = args$shape_variable, title = args$title,
                data = omicNormSha(lapply(if (inherits(plots, "ggplot")) list(plots) else plots,
                    function(plot) guiEntryCanonical(plot$data))))
        }
        plots
    }
    # The QC figures' kernel, called through; `acc$drop` names a file removed after a stage (the
    # composite with a missing file).
    qc_kernel <- get(spec$qc_kernel, envir = asNamespace("MultiScholaR"))
    stubs[[spec$qc_kernel]] <- function(...) {
        result <- qc_kernel(...)
        if (!is.null(acc$drop) && identical(list(...)$stage, acc$drop$stage)) unlink(file.path(paths[[1L]], acc$drop$file))
        result
    }
    do.call(testthat::local_mocked_bindings, c(stubs, list(.package = "MultiScholaR")))
    # Every annotation label drawn (the placeholders, the composite's labels and titles).
    annotate <- ggplot2::annotate
    testthat::local_mocked_bindings(annotate = function(geom, x = NULL, y = NULL, ..., label = NULL, size = NULL) {
        if (!is.null(label)) acc$labels[[length(acc$labels) + 1L]] <- paste0(paste(label, collapse = " "), " [", size %||% "", "]")
        annotate(geom, x = x, y = y, ..., label = label, size = size)
    }, .package = "ggplot2")
    tab <- shiny::reactiveVal("qc")
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(spec$server, args = list(id = "norm",
        workflow_data = workflow, experiment_paths = paths, omic_type = spec$omic_type, experiment_label = "kit105",
        selected_tab = function() tab()), {
        apply_updates <- function() {
            for (i in seq_len(20L)) {
                if (!length(log$updates)) return(invisible())
                pending <- log$updates
                log$updates <- list()
                values <- stats::setNames(lapply(pending, `[[`, "value"), sub("^norm-", "", vapply(pending, `[[`, "", "id")))
                do.call(session$setInputs, values[!duplicated(names(values), fromLast = TRUE)])
            }
        }
        set <- function(...) {
            session$setInputs(...)
            apply_updates()
        }
        read <- function(name) tryCatch(output[[name]], error = function(e) {
            if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", conditionMessage(e))
        })
        norm_data <- session$returned
        # The UI's own defaults, as the browser sends them.
        defaults <- omicsUiInputDefaults(spec$ui("x"))
        names(defaults) <- sub("^x-", "", names(defaults))
        do.call(set, defaults)
        capture <- function(step) {
            session$flushReact()
            manager <- workflow$state_manager
            nd <- shiny::reactiveValuesToList(norm_data)
            log_lines <- nd$normalization_log
            new_log <- log_lines[seq_along(log_lines) > acc$log]
            acc$log <- length(log_lines)
            notes <- vapply(log$notes, function(note) paste0("[", note$type, "] ", guiEntryMask(note$message, root)), "")
            log$notes <- list()
            qc_files <- sort(list.files(paths[[1L]]))
            pngs <- file.path(paths[[1L]], qc_files[endsWith(qc_files, ".png")])
            encoded <- stats::setNames(vapply(pngs, base64enc::base64encode, ""), basename(pngs))
            out_names <- sub("^norm-", "", sort(names(session$.__enclos_env__$private$outs)))
            others <- setdiff(out_names[!grepl("^(final_qc_plot|cancor_plot_)", out_names)], c(omicNormImageIds, omicNormLabelIds))
            files <- names(guiEntryArtifacts(source_dir, root))
            record <- list(
                history = as.character(manager$getHistory()),
                tab_status = vapply(workflow$tab_status, as.character, ""),
                log = vapply(sub("^\\[[0-9:]+\\] ", "", new_log), guiEntryMask, "", roots = root, USE.NAMES = FALSE),
                notes = notes,
                norm = omicNormFlags(nd, manager),
                operations = vapply(omicNormProvenance(manager), function(entry) {
                    paste0(entry$operation, if (!is.na(entry$variant %||% NA)) paste0(":", entry$variant))
                }, ""),
                states = omicNormSha(guiEntryStateDigests(manager)),
                provenance = omicNormSha(omicNormProvenance(manager)),
                normdata = omicNormSha(guiEntryCanonical(nd[c("ruv_optimization_results", "correlation_results",
                    "itsd_candidate_ids")])),
                images = NULL, labels = NULL,
                outputs = omicNormSha(lapply(stats::setNames(others, others), function(id) guiEntryOutputValue(read(id), root))),
                final_pca = acc$pca,
                annotations = unlist(acc$labels),
                qc_files = qc_files,
                source_files = files[!startsWith(files, "<store")])
            # The grid's slots in full (F105-01 reads them); the capture holds their digests.
            slots <- list(assay_names = nd$assay_names, qc_files = qc_files,
                images = vapply(omicNormImageIds, function(id) omicNormImageShown(read(id), encoded), ""),
                labels = vapply(omicNormLabelIds, function(id) paste(read(id), collapse = " "), ""))
            record$images <- omicNormSha(slots$images)
            record$labels <- omicNormSha(slots$labels)
            # A field the step left as the previous step had it is held as "<unchanged>".
            for (field in c("history", "tab_status", "norm", "operations", "qc_files", "source_files")) {
                value <- record[[field]]
                if (identical(value, acc$previous[[field]])) record[[field]] <- "<unchanged>"
                acc$previous[[field]] <- value
            }
            acc$pca <- list()
            acc$labels <- list()
            acc$captures[[step]] <- record
            acc$slots[[step]] <- slots
        }
        ctx <- list(set = set, read = read, tab = tab, workflow = workflow, paths = paths, acc = acc,
            capture = capture, keep = function(name, value) acc$extra[[name]] <- value,
            flags = function() omicNormFlags(shiny::isolate(shiny::reactiveValuesToList(norm_data)), workflow$state_manager),
            closed = function() session$isClosed(),
            notes = function() {
                notes <- vapply(log$notes, function(note) paste0("[", note$type, "] ", guiEntryMask(note$message, root)), "")
                log$notes <- list()
                notes
            })
        for (step in names(steps)) steps[[step]](ctx)
    }))))
    list(captures = acc$captures, extra = acc$extra, choices = lapply(log$choices, as.character), slots = acc$slots)
}

#' The export's session RDS and summary as the step left them.
omicNormExportRecord <- function(source_dir) {
    out <- list()
    latest <- list.files(source_dir, pattern = "_filtered_session_data_latest\\.rds$", full.names = TRUE)
    if (length(latest)) {
        session <- readRDS(latest[[1L]])
        out$rds_fields <- names(session)
        out$feature_counts <- vapply(names(session$feature_counts), function(assay) {
            counts <- session$feature_counts[[assay]]
            sprintf("%s: %s features, %s samples", assay, counts$features, counts$samples)
        }, "", USE.NAMES = FALSE)
        out$rds_values <- omicNormSha(guiEntryCanonical(session[setdiff(names(session),
            c("export_timestamp", "current_s4_object"))]))
        out$current_s4_object <- guiEntryDigest(session$current_s4_object)
    }
    summary <- list.files(source_dir, pattern = "_filtered_session_summary\\.txt$", full.names = TRUE)
    if (length(summary)) out$summary <- strsplit(guiEntryMask(readLines(summary[[1L]], warn = FALSE), source_dir), "\n")[[1L]]
    out
}

#' The steps of each scenario, by name.
omicNormScenarioSteps <- function(scenario) {
    run <- function(ctx, n, ...) {
        ctx$set(...)
        ctx$set(run_normalization = n)
    }
    switch(scenario,
        # RUV skip, ITSD off, cyclicloess: the tab entry's pre-QC, the export refused before the run,
        # the run, the correlation filter (0.5), the final QC plot, the export and the reset.
        skip = list(
            entry = function(ctx) {
                ctx$tab("norm")
                ctx$capture("entry")
            },
            export_refused = function(ctx) {
                ctx$set(export_session = 1L)
                ctx$capture("export_refused")
            },
            run = function(ctx) {
                run(ctx, 1L, apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "skip")
                ctx$capture("run")
            },
            correlation = function(ctx) {
                ctx$set(min_pearson_correlation_threshold = 0.5)
                ctx$set(apply_correlation_filter = 1L)
                ctx$capture("correlation")
            },
            final_plot = function(ctx) {
                value <- ctx$read("final_qc_plot")
                ctx$keep("final_plot", if (is.list(value) && !is.null(value$src)) "<plot>" else paste(value, collapse = " "))
                ctx$capture("final_plot")
            },
            export = function(ctx) {
                ctx$set(export_session = 2L)
                ctx$keep("export", omicNormExportRecord(ctx$paths$source_dir))
                ctx$capture("export")
            },
            reset = function(ctx) {
                ctx$set(reset_normalization = 1L)
                ctx$capture("reset")
            }),
        # RUV manual (k 2, 5%), a manual ITSD selection (IS_A in the first assay, IS_B in the
        # second), quantile; correlation skipped; reset.
        manual = list(
            run = function(ctx) {
                ctx$set(itsd_table_lcms_pos_rows_selected = 2L, itsd_table_gcms_rows_selected = 1L)
                run(ctx, 1L, apply_itsd = TRUE, norm_method = "quantile", ruv_mode = "manual", ruv_k = 2, ruv_percentage = 5)
                ctx$capture("run")
            },
            skip_correlation = function(ctx) {
                ctx$set(skip_correlation_filter = 1L)
                ctx$capture("skip_correlation")
            },
            reset = function(ctx) {
                ctx$set(reset_normalization = 1L)
                ctx$capture("reset")
            }),
        # RUV automatic, ITSD by the state's pattern (no manual selection), cyclicloess, with the first
        # assay's post-normalisation density figure removed before the composite (c); the final QC
        # plot; (c) the final PCA failing; (c) skip-correlation with the store commit failing.
        automatic = list(
            run = function(ctx) {
                ctx$acc$drop <- list(stage = "post_norm", file = "lcms_pos_post_norm_density.png")
                run(ctx, 1L, apply_itsd = TRUE, norm_method = "cyclicloess", ruv_mode = "automatic")
                ctx$acc$drop <- NULL
                ctx$capture("run")
            },
            final_plot = function(ctx) {
                value <- ctx$read("final_qc_plot")
                ctx$keep("final_plot", if (is.list(value) && !is.null(value$src)) "<plot>" else paste(value, collapse = " "))
                ctx$capture("final_plot")
            },
            final_plot_failing = function(ctx) {
                ctx$acc$fail <- "kit105: the final PCA is refused"
                ctx$set(shape_variable = "batch")
                value <- ctx$read("final_qc_plot")
                ctx$keep("final_plot_failing", if (is.list(value) && !is.null(value$src)) "<plot>" else paste(value, collapse = " "))
                ctx$acc$fail <- NULL
                ctx$capture("final_plot_failing")
            },
            # Last: the observer's error ends the session (shiny closes it), so what the tab held is
            # read when the commit is refused, and the bus and owner after.
            skip_correlation_failing = function(ctx) {
                held <- new.env(parent = emptyenv())
                testthat::local_mocked_bindings(commitWorkflowReferenceState = function(...) {
                    held$at_commit <- ctx$flags()
                    stop("kit105: the store commit is refused")
                }, .package = "MultiScholaR")
                held$errors <- character()
                withCallingHandlers(ctx$set(skip_correlation_filter = 1L), warning = function(w) {
                    if (grepl("kit105", conditionMessage(w), fixed = TRUE)) {
                        held$errors <- c(held$errors, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }
                })
                ctx$keep("skip_correlation_failing", list(at_commit = held$at_commit, errors = held$errors,
                    session_closed = ctx$closed(), history = as.character(ctx$workflow$state_manager$getHistory()),
                    tab_status = vapply(ctx$workflow$tab_status, as.character, ""), notes = ctx$notes()))
            }))
}

# Each scenario runs once per omic in this file; the blocks that read it share the run.
omicNormRuns <- new.env(parent = emptyenv())

omicNormScenario <- function(omic, scenario) {
    key <- paste(omic, scenario)
    if (is.null(omicNormRuns[[key]])) {
        omicNormRuns[[key]] <- omicNormDrive(omic, omicNormScenarioSteps(scenario))
        if (!is.null(omicNormDumpPath())) {
            run <- omicNormRuns[[key]]
            omicNormDump(sprintf("omicNormFrozen[[\"%s\"]][[\"%s\"]]", omic, scenario), run[c("captures", "extra", "choices")])
            omicNormDump(sprintf("omicNormSlotsSeen[[\"%s\"]][[\"%s\"]]", omic, scenario), run$slots)
        }
    }
    omicNormRuns[[key]]
}

# The frozen captures, per omic and scenario (the dump's `omicNormFrozen` entries, pasted; its
# `omicNormSlotsSeen` entries are each step's full image and label slots, for reading only).
# OMIC_NORM_FROZEN_START
omicNormFrozen <- list(lipid = list(), metab = list())
omicNormFrozen[["lipid"]][["skip"]] <- list(
    captures = list(
        entry = list(
            history = c("initial", "lipid_qc_complete"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            log = character(0),
            notes = character(0),
            norm = c(
                normalization_complete = "FALSE", ruv_complete = "FALSE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "<NULL>",
                post_itsd_obj = "<NULL>", post_log2_obj = "<NULL>", post_norm_obj = "<NULL>",
                ruv_corrected_obj = "<NULL>", correlation_filtered_obj = "<NULL>"),
            operations = c(initial = "import", lipid_qc_complete = "legacy"),
            states = "18718ca7e8dc19cb278104f17e926b39a4fed1c3e4da6849abd58f74fb9ba00d",
            provenance = "e5742a3192359a1dcb99c39a6776d7dfc81dcd2ef7ef37d78336723d9722c711",
            normdata = "a1e8493f5e12162b1d94c08c2995c23fcdc9c9deece0d43cc697982c9c3975ec",
            images = "168bd32fe90bada9a2bdc262834b45c767c081475a5683c0aee201d46310b7b5",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "42d09557af7dc50f2f25ce4e50d26dde6f24f737f258e34817e263b99a428111",
            final_pca = list(),
            annotations = NULL,
            qc_files = c("gcms_pre_norm_correlation.png", "lcms_pos_pre_norm_correlation.png"),
            source_files = character(0)),
        export_refused = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = "[warning] Please complete normalization before exporting session data.",
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "18718ca7e8dc19cb278104f17e926b39a4fed1c3e4da6849abd58f74fb9ba00d",
            provenance = "e5742a3192359a1dcb99c39a6776d7dfc81dcd2ef7ef37d78336723d9722c711",
            normdata = "a1e8493f5e12162b1d94c08c2995c23fcdc9c9deece0d43cc697982c9c3975ec",
            images = "168bd32fe90bada9a2bdc262834b45c767c081475a5683c0aee201d46310b7b5",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "42d09557af7dc50f2f25ce4e50d26dde6f24f737f258e34817e263b99a428111",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        run = list(
            history = c("initial", "lipid_qc_complete", "lipid_log2", "lipid_normalized"),
            tab_status = "<unchanged>",
            log = c(
                "Starting normalization pipeline...", "Post-filtering state captured",
                "Pre-normalization QC plots generated", "ITSD normalization skipped",
                "Applying log2 transformation (offset: 1 )", "Log2 transformation complete",
                "Applying between-sample normalization (method: cyclicloess )",
                "Between-sample normalization complete", "Post-normalization QC plots generated",
                "RUV-III skipped", "Generating composite QC figure...",
                "Composite QC figure saved to: <ROOT>/qc/composite_QC_figure",
                "Normalization pipeline complete!"),
            notes = "[message] Normalization pipeline complete!",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "lipid_qc_complete",
                post_itsd_obj = "<NULL>", post_log2_obj = "lipid_log2", post_norm_obj = "lipid_normalized",
                ruv_corrected_obj = "lipid_normalized", correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole"),
            states = "61f40c879fcabd7fd12ac7ab2434bb680513f0d441aa67b0a04f1cbe65f2a040",
            provenance = "5a589b9c8122fd31678616579327660b75dc8cba51ef1938b2fa8a4d414dacb7",
            normdata = "a1e8493f5e12162b1d94c08c2995c23fcdc9c9deece0d43cc697982c9c3975ec",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "7aeee0129a88196d778d76329d167a0f4d89d90f224e097b75a6bfe2752f4de4",
            final_pca = list(
                list(
                    class = "LipidomicsAssayData",
                    state = "814d34425159bfd122c00f94d6e1f83f2af248cb58b7c75a4b9e89c9eb4b08af",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "3d05a34b14714cdffd228befc1557670bdb3526d022459d32d07804a6ba5e851")),
            annotations = c(
                "Pre-Normalisation [6]", "Post-Normalisation [6]", "a) [5]", "b) [5]", "c) [5]", "d) [5]",
                "e) [5]", "f) [5]", "g) [5]", "h) [5]", "i) [5]", "j) [5]", "k) [5]", "l) [5]", "m) [5]",
                "n) [5]", "o) [5]", "p) [5]"),
            qc_files = c(
                "gcms_post_norm_correlation.png", "gcms_post_norm_density.png", "gcms_post_norm_pca.png",
                "gcms_post_norm_rle.png", "gcms_pre_norm_correlation.png",
                "lcms_pos_post_norm_correlation.png", "lcms_pos_post_norm_density.png",
                "lcms_pos_post_norm_pca.png", "lcms_pos_post_norm_rle.png",
                "lcms_pos_pre_norm_correlation.png", "lipidomics_composite_QC_figure.pdf",
                "lipidomics_composite_QC_figure.png", "lipidomics_composite_QC_figure.rds"),
            source_files = "<unchanged>"),
        correlation = list(
            history = c(
                "initial", "lipid_qc_complete", "lipid_log2", "lipid_normalized",
                "lipid_correlation_filtered"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "complete", differential_analysis = "locked"),
            log = c(
                "Applying correlation filter (threshold: 0.5 )", "Correlation filtering complete"),
            notes = c(
                "[default] Applying correlation filter...",
                "[message] Correlation filtering complete! Ready for DE analysis."),
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "TRUE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "lipid_qc_complete",
                post_itsd_obj = "<NULL>", post_log2_obj = "lipid_log2", post_norm_obj = "lipid_normalized",
                ruv_corrected_obj = "lipid_normalized",
                correlation_filtered_obj = "lipid_correlation_filtered"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole",
                lipid_correlation_filtered = "lipid_pearson_selection:whole"),
            states = "5fdeff360d4eb520f3e5526f0bd93ce3356c1b0cde259d712fb4f57bd43329c7",
            provenance = "8310d5ce05971ad3ea733765495f8c50de6da6b45f952e96f3cf4aa36ba4bc73",
            normdata = "c9f2bad45bdd42f860944e6ab72e30ef2424745eb1ed20edb74b2519746cca49",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "ef1c5f8b2a20b5698b201bcf7759ce9a68f30e32bc726ecb6ab2a4898aa84c10",
            final_pca = list(
                list(
                    class = "LipidomicsAssayData",
                    state = "814d34425159bfd122c00f94d6e1f83f2af248cb58b7c75a4b9e89c9eb4b08af",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "3d05a34b14714cdffd228befc1557670bdb3526d022459d32d07804a6ba5e851")),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        final_plot = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = character(0),
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "5fdeff360d4eb520f3e5526f0bd93ce3356c1b0cde259d712fb4f57bd43329c7",
            provenance = "8310d5ce05971ad3ea733765495f8c50de6da6b45f952e96f3cf4aa36ba4bc73",
            normdata = "c9f2bad45bdd42f860944e6ab72e30ef2424745eb1ed20edb74b2519746cca49",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "ef1c5f8b2a20b5698b201bcf7759ce9a68f30e32bc726ecb6ab2a4898aa84c10",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        export = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = "Exported comprehensive session data to: <ROOT>/source/lipid_filtered_session_data_<TIME>.rds",
            notes = "[message] Session data exported successfully!\nSaved as: lipid_filtered_session_data_<TIME>.rds\nSee summary file for details.",
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "5fdeff360d4eb520f3e5526f0bd93ce3356c1b0cde259d712fb4f57bd43329c7",
            provenance = "8310d5ce05971ad3ea733765495f8c50de6da6b45f952e96f3cf4aa36ba4bc73",
            normdata = "c9f2bad45bdd42f860944e6ab72e30ef2424745eb1ed20edb74b2519746cca49",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "6fba3782e9caa36e4caea859ca2e25678d9a52b91f619319d98a7c5f0f390339",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = c(
                "export_complete_dia-export-<HEX>.json", "lipid_filtered_session_data_<TIME>.rds",
                "lipid_filtered_session_data_latest.rds", "lipid_filtered_session_summary.txt",
                "lipid_qc_params.RDS")),
        reset = list(
            history = c(
                "initial", "lipid_qc_complete", "lipid_log2", "lipid_normalized",
                "lipid_correlation_filtered", "lipid_reset"),
            tab_status = "<unchanged>",
            log = "Reset to pre-normalization state",
            notes = "[message] Reset to pre-normalization state",
            norm = c(
                normalization_complete = "FALSE", ruv_complete = "FALSE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "lipid_qc_complete",
                post_itsd_obj = "<NULL>", post_log2_obj = "lipid_log2", post_norm_obj = "<NULL>",
                ruv_corrected_obj = "<NULL>", correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole",
                lipid_correlation_filtered = "lipid_pearson_selection:whole",
                lipid_reset = "lipid_reset:whole"),
            states = "6933f8e4ed8fd73cf488b5f74f4ee27c6dc9ffc8054dfec7e9ac8e3d3a9f1ff3",
            provenance = "32e9ecb92dd79a4f63e102e94d532f30d8481d24293e6cd855e0afc8237a0042",
            normdata = "c9f2bad45bdd42f860944e6ab72e30ef2424745eb1ed20edb74b2519746cca49",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "e3e065ae81eb10a9501d1d228c4a33e0a10fb9437014859b302b82992200ae98",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>")),
    extra = list(
        final_plot = "<plot>",
        export = list(
            rds_fields = c(
                "r6_current_state_name", "current_s4_object", "contrasts_tbl", "design_matrix", "config_list",
                "itsd_selections", "ruv_optimization_results", "correlation_results", "assay_names",
                "export_timestamp", "omic_type", "experiment_label", "normalization_method", "ruv_mode",
                "itsd_applied", "itsd_aggregation", "log_offset", "correlation_threshold",
                "ruv_grouping_variable", "feature_counts", "lipid_counts", "qc_params",
                "normalization_complete", "ruv_complete", "correlation_filtering_complete"),
            feature_counts = c("LCMS_Pos: 14 features, 4 samples", "GCMS: 14 features, 4 samples"),
            rds_values = "c03c89ed6eb4bc9f5daeb547fed1e0a78bbf4f89193013b4ce835d6def74faf3",
            current_s4_object = "814d34425159bfd122c00f94d6e1f83f2af248cb58b7c75a4b9e89c9eb4b08af",
            summary = c(
                "Lipidomics Normalized Session Data Export Summary",
                "===================================================", "", "Export Timestamp: <TIME>",
                "Session File: lipid_filtered_session_data_<TIME>.rds", "", "Data Summary:",
                "  LCMS_Pos: 14 features, 4 samples", "  GCMS: 14 features, 4 samples", "",
                "Normalization Parameters:", "- Method: cyclicloess", "- ITSD applied: No",
                "- ITSD aggregation: N/A", "- Log2 offset: 1", "- RUV mode: skip",
                "- RUV grouping variable: group", "- Correlation threshold: 0.5", "",
                "RUV Optimization Results (per-assay):", "  (RUV skipped or not applied)", "", "Contrasts:",
                "KO_vs_WT", "", "This data is ready for differential expression analysis.",
                "Use 'Load Filtered Session' in the DE tab to import."))),
    choices = list(
        `norm-color_variable` = c("group", "batch"),
        `norm-shape_variable` = c("group", "batch"),
        `norm-ruv_grouping_variable` = c("group", "batch")))
omicNormFrozen[["metab"]][["skip"]] <- list(
    captures = list(
        entry = list(
            history = c("initial", "metab_qc_complete"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            log = character(0),
            notes = character(0),
            norm = c(
                normalization_complete = "FALSE", ruv_complete = "FALSE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "<NULL>",
                post_itsd_obj = "<NULL>", post_log2_obj = "<NULL>", post_norm_obj = "<NULL>",
                ruv_corrected_obj = "<NULL>", correlation_filtered_obj = "<NULL>"),
            operations = c(initial = "import", metab_qc_complete = "legacy"),
            states = "ae77f074a695bb53081b733e93d8f4192000ed74f1ae3eb7ac862b19b1ff3855",
            provenance = "206b5b584d03d53aedd4c17c18c117f3c4987d73e0d10667ac39d7c4860a9e8a",
            normdata = "ca1f0c4aeca22e0dd27fbea0183feb769fa43a36f40a40846f46eecefac62a7d",
            images = "168bd32fe90bada9a2bdc262834b45c767c081475a5683c0aee201d46310b7b5",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "42d09557af7dc50f2f25ce4e50d26dde6f24f737f258e34817e263b99a428111",
            final_pca = list(),
            annotations = NULL,
            qc_files = c("gcms_pre_norm_correlation.png", "lcms_pos_pre_norm_correlation.png"),
            source_files = character(0)),
        export_refused = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = "[warning] Please complete normalization before exporting session data.",
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "ae77f074a695bb53081b733e93d8f4192000ed74f1ae3eb7ac862b19b1ff3855",
            provenance = "206b5b584d03d53aedd4c17c18c117f3c4987d73e0d10667ac39d7c4860a9e8a",
            normdata = "ca1f0c4aeca22e0dd27fbea0183feb769fa43a36f40a40846f46eecefac62a7d",
            images = "168bd32fe90bada9a2bdc262834b45c767c081475a5683c0aee201d46310b7b5",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "42d09557af7dc50f2f25ce4e50d26dde6f24f737f258e34817e263b99a428111",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        run = list(
            history = c("initial", "metab_qc_complete", "metab_log2", "metab_normalized"),
            tab_status = "<unchanged>",
            log = c(
                "Starting normalization pipeline...", "Post-filtering state captured",
                "Pre-normalization QC plots generated", "ITSD normalization skipped",
                "Applying log2 transformation (offset: 1 )", "Log2 transformation complete",
                "Applying between-sample normalization (method: cyclicloess )",
                "Between-sample normalization complete", "Post-normalization QC plots generated",
                "RUV-III skipped", "Generating composite QC figure...",
                "Composite QC figure saved to: <ROOT>/qc/composite_QC_figure",
                "Normalization pipeline complete!"),
            notes = "[message] Normalization pipeline complete!",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "metab_qc_complete",
                post_itsd_obj = "<NULL>", post_log2_obj = "metab_log2", post_norm_obj = "metab_normalized",
                ruv_corrected_obj = "metab_normalized", correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole"),
            states = "7b0ea9f8ad3e6b15ef1dc442f32afc0d7d81c14b13694ef43aa7d9d312ebb86d",
            provenance = "5cbe229880e1540afded09e9ba53bd0e1da48128d14135635be56918660e67ca",
            normdata = "ca1f0c4aeca22e0dd27fbea0183feb769fa43a36f40a40846f46eecefac62a7d",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "7aeee0129a88196d778d76329d167a0f4d89d90f224e097b75a6bfe2752f4de4",
            final_pca = list(
                list(
                    class = "MetaboliteAssayData",
                    state = "60e0303f4e7398a75611b1092b5dda9b4c1da4d40f1253941b372de24d47039f",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "3d05a34b14714cdffd228befc1557670bdb3526d022459d32d07804a6ba5e851")),
            annotations = c(
                "Pre-Normalisation [6]", "Post-Normalisation [6]", "a) [5]", "b) [5]", "c) [5]", "d) [5]",
                "e) [5]", "f) [5]", "g) [5]", "h) [5]", "i) [5]", "j) [5]", "k) [5]", "l) [5]", "m) [5]",
                "n) [5]", "o) [5]", "p) [5]"),
            qc_files = c(
                "gcms_post_norm_correlation.png", "gcms_post_norm_density.png", "gcms_post_norm_pca.png",
                "gcms_post_norm_rle.png", "gcms_pre_norm_correlation.png",
                "lcms_pos_post_norm_correlation.png", "lcms_pos_post_norm_density.png",
                "lcms_pos_post_norm_pca.png", "lcms_pos_post_norm_rle.png",
                "lcms_pos_pre_norm_correlation.png", "metabolomics_composite_QC_figure.pdf",
                "metabolomics_composite_QC_figure.png", "metabolomics_composite_QC_figure.rds"),
            source_files = "<unchanged>"),
        correlation = list(
            history = c(
                "initial", "metab_qc_complete", "metab_log2", "metab_normalized",
                "metab_correlation_filtered"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "complete", differential_analysis = "locked"),
            log = c(
                "Applying correlation filter (threshold: 0.5 )", "Correlation filtering complete"),
            notes = c(
                "[default] Applying correlation filter...",
                "[message] Correlation filtering complete! Ready for DE analysis."),
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "TRUE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "metab_qc_complete",
                post_itsd_obj = "<NULL>", post_log2_obj = "metab_log2", post_norm_obj = "metab_normalized",
                ruv_corrected_obj = "metab_normalized",
                correlation_filtered_obj = "metab_correlation_filtered"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole",
                metab_correlation_filtered = "metab_pearson_selection:whole"),
            states = "9d534995d3baf87f9269f9631119bc1f27e350c87e51c681fde68a3e66819a4a",
            provenance = "85fb639d11822112feef3d817c637fe05689027dc2bb722f0b6e48845a667699",
            normdata = "332c0f56d24a1a46e6f0bf914a98e41cc5efcfd9ada0d027b16b808b0f26948c",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "ef1c5f8b2a20b5698b201bcf7759ce9a68f30e32bc726ecb6ab2a4898aa84c10",
            final_pca = list(
                list(
                    class = "MetaboliteAssayData",
                    state = "60e0303f4e7398a75611b1092b5dda9b4c1da4d40f1253941b372de24d47039f",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "3d05a34b14714cdffd228befc1557670bdb3526d022459d32d07804a6ba5e851")),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        final_plot = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = character(0),
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "9d534995d3baf87f9269f9631119bc1f27e350c87e51c681fde68a3e66819a4a",
            provenance = "85fb639d11822112feef3d817c637fe05689027dc2bb722f0b6e48845a667699",
            normdata = "332c0f56d24a1a46e6f0bf914a98e41cc5efcfd9ada0d027b16b808b0f26948c",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "ef1c5f8b2a20b5698b201bcf7759ce9a68f30e32bc726ecb6ab2a4898aa84c10",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        export = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = "Exported comprehensive session data to: <ROOT>/source/metab_filtered_session_data_<TIME>.rds",
            notes = "[message] Session data exported successfully!\nSaved as: metab_filtered_session_data_<TIME>.rds\nSee summary file for details.",
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "9d534995d3baf87f9269f9631119bc1f27e350c87e51c681fde68a3e66819a4a",
            provenance = "85fb639d11822112feef3d817c637fe05689027dc2bb722f0b6e48845a667699",
            normdata = "332c0f56d24a1a46e6f0bf914a98e41cc5efcfd9ada0d027b16b808b0f26948c",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "b0e80029a5cbceeb6bb770b5edde364ae209ce52dcf0ed7bb58482476ae4d874",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = c(
                "export_complete_dia-export-<HEX>.json", "metab_filtered_session_data_<TIME>.rds",
                "metab_filtered_session_data_latest.rds", "metab_filtered_session_summary.txt",
                "metab_qc_params.RDS")),
        reset = list(
            history = c(
                "initial", "metab_qc_complete", "metab_log2", "metab_normalized",
                "metab_correlation_filtered", "metab_reset"),
            tab_status = "<unchanged>",
            log = "Reset to pre-normalization state",
            notes = "[message] Reset to pre-normalization state",
            norm = c(
                normalization_complete = "FALSE", ruv_complete = "FALSE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "TRUE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "metab_qc_complete",
                post_itsd_obj = "<NULL>", post_log2_obj = "metab_log2", post_norm_obj = "<NULL>",
                ruv_corrected_obj = "<NULL>", correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole",
                metab_correlation_filtered = "metab_pearson_selection:whole",
                metab_reset = "metab_reset:whole"),
            states = "1c36f2cd9a29076c83a5558a69e7ceea5537d6309650a28a2bb65f1ec940101b",
            provenance = "615eedf13522b608459e38030baa010ec309cda6e2f3f79f6d456492d847a83f",
            normdata = "332c0f56d24a1a46e6f0bf914a98e41cc5efcfd9ada0d027b16b808b0f26948c",
            images = "1cfa99b8bcf8de32b177068e75bf057617d7e7d23f564d96ef264f48846f6740",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "b82b6405c8346bca747500a35fb834518f024a69af4e2c9cdaa4886ad0a2231c",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>")),
    extra = list(
        final_plot = "<plot>",
        export = list(
            rds_fields = c(
                "r6_current_state_name", "current_s4_object", "contrasts_tbl", "design_matrix", "config_list",
                "itsd_selections", "ruv_optimization_results", "correlation_results", "assay_names",
                "export_timestamp", "omic_type", "experiment_label", "normalization_method", "ruv_mode",
                "itsd_applied", "itsd_aggregation", "log_offset", "correlation_threshold",
                "ruv_grouping_variable", "feature_counts", "metabolite_counts", "qc_params",
                "normalization_complete", "ruv_complete", "correlation_filtering_complete"),
            feature_counts = c("LCMS_Pos: 14 features, 5 samples", "GCMS: 14 features, 5 samples"),
            rds_values = "8956c7f7d24e3399f9eb87d62a74e55655b09c3a56fb16dbc2b8645e92eb86d1",
            current_s4_object = "60e0303f4e7398a75611b1092b5dda9b4c1da4d40f1253941b372de24d47039f",
            summary = c(
                "Metabolomics Normalized Session Data Export Summary",
                "===================================================", "", "Export Timestamp: <TIME>",
                "Session File: metab_filtered_session_data_<TIME>.rds", "", "Data Summary:",
                "  LCMS_Pos: 14 features, 5 samples", "  GCMS: 14 features, 5 samples", "",
                "Normalization Parameters:", "- Method: cyclicloess", "- ITSD applied: No",
                "- ITSD aggregation: N/A", "- Log2 offset: 1", "- RUV mode: skip",
                "- RUV grouping variable: group", "- Correlation threshold: 0.5", "",
                "RUV Optimization Results (per-assay):", "  (RUV skipped or not applied)", "", "Contrasts:",
                "KO_vs_WT", "", "This data is ready for differential expression analysis.",
                "Use 'Load Filtered Session' in the DE tab to import."))),
    choices = list(
        `norm-color_variable` = c("group", "batch"),
        `norm-shape_variable` = c("group", "batch"),
        `norm-ruv_grouping_variable` = c("group", "batch")))
omicNormFrozen[["lipid"]][["manual"]] <- list(
    captures = list(
        run = list(
            history = c(
                "initial", "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized",
                "lipid_ruv_corrected"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            log = c(
                "Starting normalization pipeline...", "Post-filtering state captured",
                "Pre-normalization QC plots generated", "Applying ITSD normalization (aggregation: median )",
                "Assay LCMS_Pos : 1 ITSD features selected", "Assay GCMS : 1 ITSD features selected",
                "ITSD normalization complete", "Applying log2 transformation (offset: 1 )",
                "Log2 transformation complete", "Applying between-sample normalization (method: quantile )",
                "Between-sample normalization complete", "Post-normalization QC plots generated",
                "Running RUV-III (mode: manual )",
                "RUV optimization complete. Best k per assay: LCMS_Pos = 2, GCMS = 2",
                "RUV-III correction applied", "RUV QC plots generated", "Generating composite QC figure...",
                "Composite QC figure saved to: <ROOT>/qc/composite_QC_figure",
                "Normalization pipeline complete!"),
            notes = "[message] Normalization pipeline complete!",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "LCMS_Pos=2;GCMS=1",
                post_filter_obj = "lipid_qc_complete", post_itsd_obj = "lipid_itsd_norm",
                post_log2_obj = "lipid_log2", post_norm_obj = "lipid_normalized",
                ruv_corrected_obj = "lipid_ruv_corrected", correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy",
                lipid_itsd_norm = "lipid_itsd_normalisation:whole", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole",
                lipid_ruv_corrected = "lipid_ruv:whole"),
            states = "c974c0dfdc29c9567fd76dae509a47f5d552f274a0936dcb64c793ede3636ce3",
            provenance = "a01af382e056b5e26bff87332f81543caae071eea65fecb9e7b4c97d3e897ab2",
            normdata = "d73553f170cf0ca694467a4ff6e7a214cd714c4e50778c6182bb0b4964d8070e",
            images = "c64a1c2a13649637d602c3406a6f3865984a274768deb19baab9285e6bb1b398",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "45543d7e025cd734bf76e1f34192064ed900a2c5eaa0189911b31cc95a3ff0c7",
            final_pca = list(
                list(
                    class = "LipidomicsAssayData",
                    state = "b1ba801f9c3a13a5651eef5f6b3e8fa822a29c0b46cb9e3e14d766cddef2997c",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "afce65984dc7e339bd3b3e914a6f87707c2b2b02675c58c36553de610e5ccf2e")),
            annotations = c(
                "Pre-Normalisation [6]", "Post-Normalisation [6]", "RUV-Corrected [6]", "a) [5]", "b) [5]",
                "c) [5]", "d) [5]", "e) [5]", "f) [5]", "g) [5]", "h) [5]", "i) [5]", "j) [5]", "k) [5]",
                "l) [5]", "m) [5]", "n) [5]", "o) [5]", "p) [5]", "q) [5]", "r) [5]", "s) [5]", "t) [5]",
                "u) [5]", "v) [5]", "w) [5]", "x) [5]"),
            qc_files = c(
                "gcms_post_norm_correlation.png", "gcms_post_norm_density.png", "gcms_post_norm_pca.png",
                "gcms_post_norm_rle.png", "gcms_pre_norm_correlation.png",
                "gcms_ruv_corrected_correlation.png", "gcms_ruv_corrected_density.png",
                "gcms_ruv_corrected_pca.png", "gcms_ruv_corrected_rle.png",
                "lcms_pos_post_norm_correlation.png", "lcms_pos_post_norm_density.png",
                "lcms_pos_post_norm_pca.png", "lcms_pos_post_norm_rle.png",
                "lcms_pos_pre_norm_correlation.png", "lcms_pos_ruv_corrected_correlation.png",
                "lcms_pos_ruv_corrected_density.png", "lcms_pos_ruv_corrected_pca.png",
                "lcms_pos_ruv_corrected_rle.png", "lipidomics_composite_QC_figure.pdf",
                "lipidomics_composite_QC_figure.png", "lipidomics_composite_QC_figure.rds"),
            source_files = character(0)),
        skip_correlation = list(
            history = c(
                "initial", "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized",
                "lipid_ruv_corrected", "lipid_norm_complete"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "complete", differential_analysis = "locked"),
            log = "Correlation filtering skipped - ready for DE analysis",
            notes = "[message] Normalization complete! Proceeding to DE analysis.",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "TRUE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "LCMS_Pos=2;GCMS=1",
                post_filter_obj = "lipid_qc_complete", post_itsd_obj = "lipid_itsd_norm",
                post_log2_obj = "lipid_log2", post_norm_obj = "lipid_normalized",
                ruv_corrected_obj = "lipid_ruv_corrected", correlation_filtered_obj = "lipid_ruv_corrected"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy",
                lipid_itsd_norm = "lipid_itsd_normalisation:whole", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole",
                lipid_ruv_corrected = "lipid_ruv:whole", lipid_norm_complete = "lipid_norm_complete:whole"),
            states = "8d48d9b982b994d92c7880fab1bf91e7ee4efcc7dd0075ee92295181dc19e643",
            provenance = "f668e5302b5b7e4078dacca6434d9e7b1a7f7ca66e19c5184d5e07aa57507d6d",
            normdata = "d73553f170cf0ca694467a4ff6e7a214cd714c4e50778c6182bb0b4964d8070e",
            images = "c64a1c2a13649637d602c3406a6f3865984a274768deb19baab9285e6bb1b398",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "1177c9e8abd01576c6c612424de3869d6a8b4b7b48136fe6c57107219bff690b",
            final_pca = list(
                list(
                    class = "LipidomicsAssayData",
                    state = "b1ba801f9c3a13a5651eef5f6b3e8fa822a29c0b46cb9e3e14d766cddef2997c",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "afce65984dc7e339bd3b3e914a6f87707c2b2b02675c58c36553de610e5ccf2e")),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        reset = list(
            history = c(
                "initial", "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized",
                "lipid_ruv_corrected", "lipid_norm_complete", "lipid_reset"),
            tab_status = "<unchanged>",
            log = "Reset to pre-normalization state",
            notes = "[message] Reset to pre-normalization state",
            norm = c(
                normalization_complete = "FALSE", ruv_complete = "FALSE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "LCMS_Pos=2;GCMS=1",
                post_filter_obj = "lipid_qc_complete", post_itsd_obj = "lipid_itsd_norm",
                post_log2_obj = "lipid_log2", post_norm_obj = "<NULL>", ruv_corrected_obj = "<NULL>",
                correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy",
                lipid_itsd_norm = "lipid_itsd_normalisation:whole", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole",
                lipid_ruv_corrected = "lipid_ruv:whole", lipid_norm_complete = "lipid_norm_complete:whole",
                lipid_reset = "lipid_reset:whole"),
            states = "7d8d63b94ea2195c066e7eb612ede33010b3f932e4e9a861210731f23a921dd2",
            provenance = "8657461e52e8cea6d8d4a5a77edb6b773dbddcdc5c0cfd47aeb8906d7a7537c2",
            normdata = "a1e8493f5e12162b1d94c08c2995c23fcdc9c9deece0d43cc697982c9c3975ec",
            images = "c64a1c2a13649637d602c3406a6f3865984a274768deb19baab9285e6bb1b398",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "85a94a38f2772cc2511058d2c8a727218e5e1a551b49579ca208147a0bf02523",
            final_pca = list(),
            annotations = c("No data [6]", "No data [6]"),
            qc_files = "<unchanged>",
            source_files = "<unchanged>")),
    extra = list(),
    choices = list(
        `norm-color_variable` = c("group", "batch"),
        `norm-shape_variable` = c("group", "batch"),
        `norm-ruv_grouping_variable` = c("group", "batch")))
omicNormFrozen[["metab"]][["manual"]] <- list(
    captures = list(
        run = list(
            history = c(
                "initial", "metab_qc_complete", "metab_itsd_norm", "metab_log2", "metab_normalized",
                "metab_ruv_corrected"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            log = c(
                "Starting normalization pipeline...", "Post-filtering state captured",
                "Pre-normalization QC plots generated", "Applying ITSD normalization (aggregation: median )",
                "Assay LCMS_Pos : 1 ITSD features selected", "Assay GCMS : 1 ITSD features selected",
                "ITSD normalization complete", "Applying log2 transformation (offset: 1 )",
                "Log2 transformation complete", "Applying between-sample normalization (method: quantile )",
                "Between-sample normalization complete", "Post-normalization QC plots generated",
                "Running RUV-III (mode: manual )",
                "RUV optimization complete. Best k per assay: LCMS_Pos = 2, GCMS = 2",
                "RUV-III correction applied", "RUV QC plots generated", "Generating composite QC figure...",
                "Composite QC figure saved to: <ROOT>/qc/composite_QC_figure",
                "Normalization pipeline complete!"),
            notes = "[message] Normalization pipeline complete!",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "LCMS_Pos=2;GCMS=1",
                post_filter_obj = "metab_qc_complete", post_itsd_obj = "metab_itsd_norm",
                post_log2_obj = "metab_log2", post_norm_obj = "metab_normalized",
                ruv_corrected_obj = "metab_ruv_corrected", correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy",
                metab_itsd_norm = "metab_itsd_normalisation:whole", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole",
                metab_ruv_corrected = "metab_ruv:whole"),
            states = "4069f7e526041c28fd55ae17d98f7f9e477992394952feb47ba90d56340241f0",
            provenance = "28f3eaebc6a1eeb064ad638a8c1b30dfa354f6ec17b464cef21a47a0b1e07404",
            normdata = "6b1b9c5dcdc5e0469b79aed84a3e0a6dfd8ba98b9a936a9a33e23be4889180b0",
            images = "c64a1c2a13649637d602c3406a6f3865984a274768deb19baab9285e6bb1b398",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "45543d7e025cd734bf76e1f34192064ed900a2c5eaa0189911b31cc95a3ff0c7",
            final_pca = list(
                list(
                    class = "MetaboliteAssayData",
                    state = "fdb6fffdbfd873d59895f150dd86a3c9dcd8e2a86598f485dc048b16b29ed6b1",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "afce65984dc7e339bd3b3e914a6f87707c2b2b02675c58c36553de610e5ccf2e")),
            annotations = c(
                "Pre-Normalisation [6]", "Post-Normalisation [6]", "RUV-Corrected [6]", "a) [5]", "b) [5]",
                "c) [5]", "d) [5]", "e) [5]", "f) [5]", "g) [5]", "h) [5]", "i) [5]", "j) [5]", "k) [5]",
                "l) [5]", "m) [5]", "n) [5]", "o) [5]", "p) [5]", "q) [5]", "r) [5]", "s) [5]", "t) [5]",
                "u) [5]", "v) [5]", "w) [5]", "x) [5]"),
            qc_files = c(
                "gcms_post_norm_correlation.png", "gcms_post_norm_density.png", "gcms_post_norm_pca.png",
                "gcms_post_norm_rle.png", "gcms_pre_norm_correlation.png",
                "gcms_ruv_corrected_correlation.png", "gcms_ruv_corrected_density.png",
                "gcms_ruv_corrected_pca.png", "gcms_ruv_corrected_rle.png",
                "lcms_pos_post_norm_correlation.png", "lcms_pos_post_norm_density.png",
                "lcms_pos_post_norm_pca.png", "lcms_pos_post_norm_rle.png",
                "lcms_pos_pre_norm_correlation.png", "lcms_pos_ruv_corrected_correlation.png",
                "lcms_pos_ruv_corrected_density.png", "lcms_pos_ruv_corrected_pca.png",
                "lcms_pos_ruv_corrected_rle.png", "metabolomics_composite_QC_figure.pdf",
                "metabolomics_composite_QC_figure.png", "metabolomics_composite_QC_figure.rds"),
            source_files = character(0)),
        skip_correlation = list(
            history = c(
                "initial", "metab_qc_complete", "metab_itsd_norm", "metab_log2", "metab_normalized",
                "metab_ruv_corrected", "metab_norm_complete"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "complete", differential_analysis = "locked"),
            log = "Correlation filtering skipped - ready for DE analysis",
            notes = "[message] Normalization complete! Proceeding to DE analysis.",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "TRUE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "LCMS_Pos=2;GCMS=1",
                post_filter_obj = "metab_qc_complete", post_itsd_obj = "metab_itsd_norm",
                post_log2_obj = "metab_log2", post_norm_obj = "metab_normalized",
                ruv_corrected_obj = "metab_ruv_corrected", correlation_filtered_obj = "metab_ruv_corrected"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy",
                metab_itsd_norm = "metab_itsd_normalisation:whole", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole",
                metab_ruv_corrected = "metab_ruv:whole", metab_norm_complete = "metab_norm_complete:whole"),
            states = "a618753a9fc35cc25d23ba9a8e8df20f0b97505d81aa0712a2f362b89cc72d77",
            provenance = "8efee29e7f4fb67a11cc18b78d0e521ba4a1befd0878d8ecc3a3be67b2f3de04",
            normdata = "6b1b9c5dcdc5e0469b79aed84a3e0a6dfd8ba98b9a936a9a33e23be4889180b0",
            images = "c64a1c2a13649637d602c3406a6f3865984a274768deb19baab9285e6bb1b398",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "1177c9e8abd01576c6c612424de3869d6a8b4b7b48136fe6c57107219bff690b",
            final_pca = list(
                list(
                    class = "MetaboliteAssayData",
                    state = "fdb6fffdbfd873d59895f150dd86a3c9dcd8e2a86598f485dc048b16b29ed6b1",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "afce65984dc7e339bd3b3e914a6f87707c2b2b02675c58c36553de610e5ccf2e")),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        reset = list(
            history = c(
                "initial", "metab_qc_complete", "metab_itsd_norm", "metab_log2", "metab_normalized",
                "metab_ruv_corrected", "metab_norm_complete", "metab_reset"),
            tab_status = "<unchanged>",
            log = "Reset to pre-normalization state",
            notes = "[message] Reset to pre-normalization state",
            norm = c(
                normalization_complete = "FALSE", ruv_complete = "FALSE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "LCMS_Pos=2;GCMS=1",
                post_filter_obj = "metab_qc_complete", post_itsd_obj = "metab_itsd_norm",
                post_log2_obj = "metab_log2", post_norm_obj = "<NULL>", ruv_corrected_obj = "<NULL>",
                correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy",
                metab_itsd_norm = "metab_itsd_normalisation:whole", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole",
                metab_ruv_corrected = "metab_ruv:whole", metab_norm_complete = "metab_norm_complete:whole",
                metab_reset = "metab_reset:whole"),
            states = "7e9b7a51a4c5c7dad726b1daeadeb52fbd10a086ef42f343cb015cb797464b09",
            provenance = "508e46943ae495bc8ed24e0abf3e9f9429e36fca39da53e7401cc09064986314",
            normdata = "ca1f0c4aeca22e0dd27fbea0183feb769fa43a36f40a40846f46eecefac62a7d",
            images = "c64a1c2a13649637d602c3406a6f3865984a274768deb19baab9285e6bb1b398",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "85a94a38f2772cc2511058d2c8a727218e5e1a551b49579ca208147a0bf02523",
            final_pca = list(),
            annotations = c("No data [6]", "No data [6]"),
            qc_files = "<unchanged>",
            source_files = "<unchanged>")),
    extra = list(),
    choices = list(
        `norm-color_variable` = c("group", "batch"),
        `norm-shape_variable` = c("group", "batch"),
        `norm-ruv_grouping_variable` = c("group", "batch")))
omicNormFrozen[["lipid"]][["automatic"]] <- list(
    captures = list(
        run = list(
            history = c(
                "initial", "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized",
                "lipid_ruv_corrected"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            log = c(
                "Starting normalization pipeline...", "Post-filtering state captured",
                "Pre-normalization QC plots generated", "Applying ITSD normalization (aggregation: median )",
                "ITSD normalization complete", "Applying log2 transformation (offset: 1 )",
                "Log2 transformation complete",
                "Applying between-sample normalization (method: cyclicloess )",
                "Between-sample normalization complete", "Post-normalization QC plots generated",
                "Running RUV-III (mode: automatic )",
                "RUV optimization complete. Best k per assay: LCMS_Pos = 1, GCMS = 2",
                "RUV-III correction applied", "RUV QC plots generated", "Generating composite QC figure...",
                "Composite QC figure saved to: <ROOT>/qc/composite_QC_figure",
                "Normalization pipeline complete!"),
            notes = "[message] Normalization pipeline complete!",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "lipid_qc_complete",
                post_itsd_obj = "lipid_itsd_norm", post_log2_obj = "lipid_log2",
                post_norm_obj = "lipid_normalized", ruv_corrected_obj = "lipid_ruv_corrected",
                correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", lipid_qc_complete = "legacy",
                lipid_itsd_norm = "lipid_itsd_normalisation:whole", lipid_log2 = "lipid_log2:partitioned",
                lipid_normalized = "lipid_between_sample_normalisation:whole",
                lipid_ruv_corrected = "lipid_ruv:whole"),
            states = "a11e397e3b5f7858bc3332997d4a69feac0ea516be79f166935a5fa6ddf1f4c1",
            provenance = "6e252e3853d400fa0823e4620e4853880059bb61dab8d9e07e7e6dc3a14c6bc0",
            normdata = "b2502fa3988a015890ef4adc797513da5d000e55d29c5edb15068f6ebdd557ac",
            images = "422c5931e587d69abdb910d0a29b9e5c5e70bf48f9de038b1d8dbea73f03c50b",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "258a56a356180f82d13f062eadee6255a15a19455795ef22ef16b2884fd8c977",
            final_pca = list(
                list(
                    class = "LipidomicsAssayData",
                    state = "52242e72d7d6bdeb5f5697dcddee7f9d0a44f49e1f21636bb3d01f43bed8d63c",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "e563dcfac02f82d815a963f457c3034f797d2eb4764ff77cf47a3a888fe025ac")),
            annotations = c(
                "Pre-Normalisation [6]", "Post-Normalisation [6]", "RUV-Corrected [6]", "a) [5]", "b) [5]",
                "c) [5]", "d) [5]", "e) [5]", "f) [5]", "g) [5]", "h) [5]", "i) [5]", "j) [5]", "k) [5]",
                "l) [5]", "m) [5]", "n) [5]", "o) [5]", "p) [5]", "q) [5]", "r) [5]", "s) [5]", "t) [5]",
                "u) [5]", "v) [5]", "w) [5]", "x) [5]"),
            qc_files = c(
                "gcms_post_norm_correlation.png", "gcms_post_norm_density.png", "gcms_post_norm_pca.png",
                "gcms_post_norm_rle.png", "gcms_pre_norm_correlation.png",
                "gcms_ruv_corrected_correlation.png", "gcms_ruv_corrected_density.png",
                "gcms_ruv_corrected_pca.png", "gcms_ruv_corrected_rle.png",
                "lcms_pos_post_norm_correlation.png", "lcms_pos_post_norm_pca.png",
                "lcms_pos_post_norm_rle.png", "lcms_pos_pre_norm_correlation.png",
                "lcms_pos_ruv_corrected_correlation.png", "lcms_pos_ruv_corrected_density.png",
                "lcms_pos_ruv_corrected_pca.png", "lcms_pos_ruv_corrected_rle.png",
                "lipidomics_composite_QC_figure.pdf", "lipidomics_composite_QC_figure.png",
                "lipidomics_composite_QC_figure.rds"),
            source_files = character(0)),
        final_plot = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = character(0),
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "a11e397e3b5f7858bc3332997d4a69feac0ea516be79f166935a5fa6ddf1f4c1",
            provenance = "6e252e3853d400fa0823e4620e4853880059bb61dab8d9e07e7e6dc3a14c6bc0",
            normdata = "b2502fa3988a015890ef4adc797513da5d000e55d29c5edb15068f6ebdd557ac",
            images = "422c5931e587d69abdb910d0a29b9e5c5e70bf48f9de038b1d8dbea73f03c50b",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "258a56a356180f82d13f062eadee6255a15a19455795ef22ef16b2884fd8c977",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        final_plot_failing = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = character(0),
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "a11e397e3b5f7858bc3332997d4a69feac0ea516be79f166935a5fa6ddf1f4c1",
            provenance = "6e252e3853d400fa0823e4620e4853880059bb61dab8d9e07e7e6dc3a14c6bc0",
            normdata = "b2502fa3988a015890ef4adc797513da5d000e55d29c5edb15068f6ebdd557ac",
            images = "422c5931e587d69abdb910d0a29b9e5c5e70bf48f9de038b1d8dbea73f03c50b",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "258a56a356180f82d13f062eadee6255a15a19455795ef22ef16b2884fd8c977",
            final_pca = list(),
            annotations = "Error: kit105: the final PCA is refused [4]",
            qc_files = "<unchanged>",
            source_files = "<unchanged>")),
    extra = list(
        final_plot = "<plot>",
        final_plot_failing = "<plot>",
        skip_correlation_failing = list(
            at_commit = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "TRUE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "lipid_qc_complete",
                post_itsd_obj = "lipid_itsd_norm", post_log2_obj = "lipid_log2",
                post_norm_obj = "lipid_normalized", ruv_corrected_obj = "lipid_ruv_corrected",
                correlation_filtered_obj = "lipid_ruv_corrected"),
            errors = "Error in commitWorkflowReferenceState: kit105: the store commit is refused",
            session_closed = TRUE,
            history = c(
                "initial", "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized",
                "lipid_ruv_corrected"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            notes = character(0))),
    choices = list(
        `norm-color_variable` = c("group", "batch"),
        `norm-shape_variable` = c("group", "batch"),
        `norm-ruv_grouping_variable` = c("group", "batch")))
omicNormFrozen[["metab"]][["automatic"]] <- list(
    captures = list(
        run = list(
            history = c(
                "initial", "metab_qc_complete", "metab_itsd_norm", "metab_log2", "metab_normalized",
                "metab_ruv_corrected"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            log = c(
                "Starting normalization pipeline...", "Post-filtering state captured",
                "Pre-normalization QC plots generated", "Applying ITSD normalization (aggregation: median )",
                "ITSD normalization complete", "Applying log2 transformation (offset: 1 )",
                "Log2 transformation complete",
                "Applying between-sample normalization (method: cyclicloess )",
                "Between-sample normalization complete", "Post-normalization QC plots generated",
                "Running RUV-III (mode: automatic )",
                "RUV optimization complete. Best k per assay: LCMS_Pos = 1, GCMS = 2",
                "RUV-III correction applied", "RUV QC plots generated", "Generating composite QC figure...",
                "Composite QC figure saved to: <ROOT>/qc/composite_QC_figure",
                "Normalization pipeline complete!"),
            notes = "[message] Normalization pipeline complete!",
            norm = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "metab_qc_complete",
                post_itsd_obj = "metab_itsd_norm", post_log2_obj = "metab_log2",
                post_norm_obj = "metab_normalized", ruv_corrected_obj = "metab_ruv_corrected",
                correlation_filtered_obj = "<NULL>"),
            operations = c(
                initial = "import", metab_qc_complete = "legacy",
                metab_itsd_norm = "metab_itsd_normalisation:whole", metab_log2 = "metab_log2:partitioned",
                metab_normalized = "metab_between_sample_normalisation:whole",
                metab_ruv_corrected = "metab_ruv:whole"),
            states = "c03173d8653f9eae6ee4a36be3be91369e0540cb54b5deb781fbcda99f5112a3",
            provenance = "a4ad3ceca6369b3eaabb5b3b5725db625dd573b444bc631bca3c5aae006ccec7",
            normdata = "260bfecdd01ac4fb72a07ca8fdd80bc9591bf1f371010c75110d292c5e16f754",
            images = "422c5931e587d69abdb910d0a29b9e5c5e70bf48f9de038b1d8dbea73f03c50b",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "258a56a356180f82d13f062eadee6255a15a19455795ef22ef16b2884fd8c977",
            final_pca = list(
                list(
                    class = "MetaboliteAssayData",
                    state = "cb3876ffff3a90915594a59b700a1b22667cd9dda7e3610b9d122df92d7bae02",
                    grouping_variable = "group",
                    shape_variable = "group",
                    title = "Final QC - PCA",
                    data = "e563dcfac02f82d815a963f457c3034f797d2eb4764ff77cf47a3a888fe025ac")),
            annotations = c(
                "Pre-Normalisation [6]", "Post-Normalisation [6]", "RUV-Corrected [6]", "a) [5]", "b) [5]",
                "c) [5]", "d) [5]", "e) [5]", "f) [5]", "g) [5]", "h) [5]", "i) [5]", "j) [5]", "k) [5]",
                "l) [5]", "m) [5]", "n) [5]", "o) [5]", "p) [5]", "q) [5]", "r) [5]", "s) [5]", "t) [5]",
                "u) [5]", "v) [5]", "w) [5]", "x) [5]"),
            qc_files = c(
                "gcms_post_norm_correlation.png", "gcms_post_norm_density.png", "gcms_post_norm_pca.png",
                "gcms_post_norm_rle.png", "gcms_pre_norm_correlation.png",
                "gcms_ruv_corrected_correlation.png", "gcms_ruv_corrected_density.png",
                "gcms_ruv_corrected_pca.png", "gcms_ruv_corrected_rle.png",
                "lcms_pos_post_norm_correlation.png", "lcms_pos_post_norm_pca.png",
                "lcms_pos_post_norm_rle.png", "lcms_pos_pre_norm_correlation.png",
                "lcms_pos_ruv_corrected_correlation.png", "lcms_pos_ruv_corrected_density.png",
                "lcms_pos_ruv_corrected_pca.png", "lcms_pos_ruv_corrected_rle.png",
                "metabolomics_composite_QC_figure.pdf", "metabolomics_composite_QC_figure.png",
                "metabolomics_composite_QC_figure.rds"),
            source_files = character(0)),
        final_plot = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = character(0),
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "c03173d8653f9eae6ee4a36be3be91369e0540cb54b5deb781fbcda99f5112a3",
            provenance = "a4ad3ceca6369b3eaabb5b3b5725db625dd573b444bc631bca3c5aae006ccec7",
            normdata = "260bfecdd01ac4fb72a07ca8fdd80bc9591bf1f371010c75110d292c5e16f754",
            images = "422c5931e587d69abdb910d0a29b9e5c5e70bf48f9de038b1d8dbea73f03c50b",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "258a56a356180f82d13f062eadee6255a15a19455795ef22ef16b2884fd8c977",
            final_pca = list(),
            annotations = NULL,
            qc_files = "<unchanged>",
            source_files = "<unchanged>"),
        final_plot_failing = list(
            history = "<unchanged>",
            tab_status = "<unchanged>",
            log = character(0),
            notes = character(0),
            norm = "<unchanged>",
            operations = "<unchanged>",
            states = "c03173d8653f9eae6ee4a36be3be91369e0540cb54b5deb781fbcda99f5112a3",
            provenance = "a4ad3ceca6369b3eaabb5b3b5725db625dd573b444bc631bca3c5aae006ccec7",
            normdata = "260bfecdd01ac4fb72a07ca8fdd80bc9591bf1f371010c75110d292c5e16f754",
            images = "422c5931e587d69abdb910d0a29b9e5c5e70bf48f9de038b1d8dbea73f03c50b",
            labels = "aaf4247ea00fa60b9ef6d876a18a08971175b031209ba32ca975b9e31310635f",
            outputs = "258a56a356180f82d13f062eadee6255a15a19455795ef22ef16b2884fd8c977",
            final_pca = list(),
            annotations = "Error: kit105: the final PCA is refused [4]",
            qc_files = "<unchanged>",
            source_files = "<unchanged>")),
    extra = list(
        final_plot = "<plot>",
        final_plot_failing = "<plot>",
        skip_correlation_failing = list(
            at_commit = c(
                normalization_complete = "TRUE", ruv_complete = "TRUE",
                correlation_filtering_complete = "FALSE", pre_norm_qc_generated = "FALSE",
                assay_names = "LCMS_Pos|GCMS", itsd_selections = "", post_filter_obj = "metab_qc_complete",
                post_itsd_obj = "metab_itsd_norm", post_log2_obj = "metab_log2",
                post_norm_obj = "metab_normalized", ruv_corrected_obj = "metab_ruv_corrected",
                correlation_filtered_obj = "<NULL>"),
            errors = "Error in commitWorkflowReferenceState: kit105: the store commit is refused",
            session_closed = TRUE,
            history = c(
                "initial", "metab_qc_complete", "metab_itsd_norm", "metab_log2", "metab_normalized",
                "metab_ruv_corrected"),
            tab_status = c(
                setup_import = "complete", design_matrix = "complete", quality_control = "complete",
                normalization = "pending", differential_analysis = "locked"),
            notes = character(0))),
    choices = list(
        `norm-color_variable` = c("group", "batch"),
        `norm-shape_variable` = c("group", "batch"),
        `norm-ruv_grouping_variable` = c("group", "batch")))
# OMIC_NORM_FROZEN_END

#' Compare one scenario's captures with the frozen ones, step by step and field by field.
expectOmicNormScenario <- function(omic, scenario) {
    actual <- omicNormScenario(omic, scenario)
    frozen <- omicNormFrozen[[omic]][[scenario]]
    expect_false(is.null(frozen), info = paste(omic, scenario, "is frozen"))
    expect_identical(names(actual$captures), names(frozen$captures), info = paste(omic, scenario, "steps"))
    for (step in names(frozen$captures)) {
        for (field in union(names(frozen$captures[[step]]), names(actual$captures[[step]]))) {
            expect_identical(actual$captures[[step]][[field]], frozen$captures[[step]][[field]],
                info = paste(omic, scenario, step, field))
        }
    }
    expect_identical(actual$extra, frozen$extra, info = paste(omic, scenario, "extra"))
    expect_identical(actual$choices, frozen$choices, info = paste(omic, scenario, "choices"))
    invisible(actual)
}

test_that("(b) lipid, RUV skip and ITSD off: entry pre-QC, export refused, run, correlation, final QC plot, export, reset", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormScenario("lipid", "skip")
})

test_that("(b) metab, RUV skip and ITSD off: entry pre-QC, export refused, run, correlation, final QC plot, export, reset", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormScenario("metab", "skip")
})

test_that("(b) lipid, RUV manual with a manual ITSD selection: run, skip correlation, reset", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormScenario("lipid", "manual")
})

test_that("(b) metab, RUV manual with a manual ITSD selection: run, skip correlation, reset", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormScenario("metab", "manual")
})

test_that("(b, c) lipid, RUV automatic: run with a figure missing at the composite, final QC plot, the final PCA failing, skip-correlation with the commit failing", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormScenario("lipid", "automatic")
})

test_that("(b, c) metab, RUV automatic: run with a figure missing at the composite, final QC plot, the final PCA failing, skip-correlation with the commit failing", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormScenario("metab", "automatic")
})

# ------------------------------------------------------------------------------------------------
# F105-01: each image slot `<kind>_<stage>_assay<n>` shows `<assay n>_<pre_norm|post_norm|
# ruv_corrected>_<kind>.png` (its alt `<kind> - <assay>`), or says that file is not generated yet;
# each label `assay<n>_label_<kind>` names assay n. Metab's hand-written bindings do; lipid's
# generated ones (registerLipidNormStaticQcImageOutputs, registerLipidNormAssayLabelOutputs) bind
# every slot to the last spec (lazy arguments read after the loop), so every lipid slot shows the
# second assay's RUV-corrected correlation figure and every label the second assay. Trx's mount
# registers the same lipid functions.
# ------------------------------------------------------------------------------------------------

#' The slot and label mismatches of one step: what each shows against what it should.
omicNormSlotMismatches <- function(capture) {
    assays <- capture$assay_names
    stages <- c(post_filter = "pre_norm", post_norm = "post_norm", ruv_corrected = "ruv_corrected")
    bad <- character()
    for (id in omicNormImageIds) {
        parts <- regmatches(id, regexec("^([a-z]+)_(post_filter|post_norm|ruv_corrected)_assay([0-9])$", id))[[1L]]
        assay <- assays[as.integer(parts[[4L]])]
        file <- sprintf("%s_%s_%s.png", gsub("[^A-Za-z0-9]", "_", tolower(assay)), stages[[parts[[3L]]]], parts[[2L]])
        shown <- capture$images[[id]]
        ok <- if (file %in% capture$qc_files) {
            file %in% strsplit(sub(" alt=.*$", "", shown), "|", fixed = TRUE)[[1L]] &&
                identical(sub("^.* alt=", "", shown), paste(parts[[2L]], "-", assay))
        } else identical(shown, paste("<none> Plot not generated yet:", file))
        if (!ok) bad <- c(bad, sprintf("%s shows '%s', not %s", id, shown, file))
    }
    for (id in omicNormLabelIds) {
        assay <- assays[as.integer(sub("^assay([0-9])_.*$", "\\1", id))]
        if (!identical(capture$labels[[id]], paste("Assay:", assay))) {
            bad <- c(bad, sprintf("%s reads '%s', not Assay: %s", id, capture$labels[[id]], assay))
        }
    }
    bad
}

expectOmicNormSlots <- function(omic) {
    for (scenario in c("skip", "manual", "automatic")) {
        slots <- omicNormScenario(omic, scenario)$slots
        for (step in names(slots)) {
            expect_identical(omicNormSlotMismatches(slots[[step]]), character(), info = paste(omic, scenario, step))
        }
    }
}

test_that("F105-01: each metab QC image slot and assay label shows its own assay, stage and kind", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormSlots("metab")
})

test_that("F105-01: each lipid QC image slot and assay label shows its own assay, stage and kind (fails before)", {
    skip_if_not(workflowOwnerCapable())
    expectOmicNormSlots("lipid")
})
