# KIT-110 W1 (b-trx): the transcriptomics normalisation tab through its public entry point,
# `mod_lipid_norm_server("norm", ..., "transcriptomics", ..., profile = trxOmicProfile())`, in
# `testServer`, written before the trx normalisation core (W2) and trx's move onto the shared frame
# (W3). The fixture is MCI-030's (the shipped MCI corpus, `module_ci_trx_project(module_ci_trx_corpus())`),
# so the file runs in check mode as MCI-029 does.
#   In process: an automatic run then a reset; a manual run (k 2, 10 %); a reset before any run then a
#   manual k above the cap (refused).
#   Managed (MCI-030.3's chain mocks): an automatic run; a manual k above the cap (refused).
# Captured after each step: the registered output names, each output's value, the `normData` flags and
# references, the history, the state digests, the record's `trx_ruvs` outcomes, `tab_status`, the
# notifications (type, duration, message), the log and the files under `gene_qc_dir` and
# `publication_graphs_dir`; managed, the worker's operations in order too.
# Frozen on the start commit (d5244ffe; docs/qa/workflow-execution/kit/kit-110/w1/test-first.txt).
# MULTISCHOLAR_KIT110_W1_DUMP=<file> writes the captures as R code (how the frozen values were made).

trxNormW1Sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

#' Where MULTISCHOLAR_KIT110_W1_DUMP points, or NULL.
trxNormW1DumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT110_W1_DUMP")
    if (nzchar(path)) path else NULL
}

#' R code for a capture (lists, named vectors and scalars): a list's elements one per line, a long
#' vector's packed into lines, indented (test-omic-norm-module.R's `omicNormCode`).
trxNormW1Code <- function(x, indent = 0L) {
    pad <- strrep(" ", indent + 4L)
    label <- function(names, i) {
        if (is.null(names) || !nzchar(names[[i]])) return("")
        paste0(if (make.names(names[[i]]) == names[[i]]) names[[i]] else paste0("`", names[[i]], "`"), " = ")
    }
    if (is.list(x)) {
        if (length(x) == 0L) return("list()")
        items <- vapply(seq_along(x), function(i) paste0(pad, label(names(x), i), trxNormW1Code(x[[i]], indent + 4L)), "")
        return(paste0("list(\n", paste(items, collapse = ",\n"), ")"))
    }
    if (is.null(x)) return("NULL")
    scalar <- function(value) deparse(value, control = c("keepNA", "keepInteger"), width.cutoff = 500L)
    if (length(x) == 0L) return(if (is.null(names(x))) scalar(x) else sprintf("stats::setNames(%s, character())", scalar(unname(x))))
    if (length(x) == 1L && is.null(names(x))) return(scalar(x))
    items <- paste0(vapply(seq_along(x), function(i) paste0(label(names(x), i), scalar(unname(x[[i]]))), ""), ",")
    items[[length(items)]] <- sub(",$", ")", items[[length(items)]])
    if (sum(nchar(items) + 1L) + indent < 90L) return(paste0("c(", paste(items, collapse = " ")))
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

#' The notifications, each with its type and duration (a `showNotification` double that records
#' `duration`, shiny's default 5 when the caller gives none), scoped to `env`.
trxNormW1InstallNotifications <- function(log, env = parent.frame()) {
    notify <- function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
        type = c("default", "message", "warning", "error"), session = NULL) {
        log$notes[[length(log$notes) + 1L]] <- list(type = type[[1L]], duration = duration,
            message = paste(as.character(ui), collapse = " "))
        invisible(id %||% "notification")
    }
    bindings <- list(showNotification = notify, removeNotification = function(...) invisible(NULL))
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "shiny", .env = env)))
    imported <- intersect(names(bindings), unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE))
    if (length(imported)) do.call(testthat::local_mocked_bindings, c(bindings[imported], list(.package = "MultiScholaR", .env = env)))
    invisible(log)
}

#' The record's last `trx_ruvs` outcomes (as MCI-030 reads them).
trxNormW1RuvsOutcomes <- function(owner) {
    record <- buildWorkflowAnalysisRecord(owner)
    steps <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)
    if (length(steps)) steps[[length(steps)]]$outcomes else list()
}

#' The committed state (or states) a `normData` slot's reference names.
trxNormW1SlotState <- function(value, manager) {
    if (is.null(value)) return("<NULL>")
    if (!isSavedWorkflowReference(value)) return(paste0("<", class(value)[[1L]], ">"))
    named <- names(Filter(function(ref) identical(ref$entity_id, value$entity_id), manager$.context$snapshot$states))
    if (length(named)) paste(named, collapse = "|") else "<unlisted reference>"
}

#' What one output shows. An image slot of the stage grid: the QC files whose bytes its `src` holds,
#' with its alt text (or the alt text alone when it has no image). Otherwise the entry suite's output
#' value (a plot as its kind and alt text, a UI or table as its digest, text masked).
trxNormW1OutputShown <- function(name, value, encoded, roots) {
    if (grepl("^(pca|density|rle|correlation)_", name) && is.list(value)) {
        if (!nzchar(value$src %||% "")) return(paste("<none>", value$alt %||% ""))
        shown <- names(encoded)[encoded == sub("^data:[^,]*,", "", value$src)]
        return(paste0(if (length(shown)) paste(shown, collapse = "|") else "<unknown>", " alt=", value$alt %||% ""))
    }
    guiEntryOutputValue(value, roots)
}

#' The server the tests mount: the trx norm tab through the lipid wrapper under the transcriptomics
#' profile (as `module_ci_trx_tabs_server`), with the omic's runner when `managed`.
trxNormW1Server <- function(fixture, managed) {
    force(fixture)
    force(managed)
    function(input, output, session) {
        if (managed) {
            runner <- installWorkflowOmicsActions(session, fixture$workflow, "transcriptomics", "kit110", fixture$paths)
            runner$progress <- module_ci_trx_quiet_progress
        }
        session$userData$norm <- mod_lipid_norm_server("norm", fixture$workflow, fixture$paths, "transcriptomics", "kit110",
            profile = trxOmicProfile())
    }
}

#' Drive the tab through `steps` (named functions of the context) and capture after each step.
trxNormW1Drive <- function(steps, managed = FALSE) {
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    paths <- fixture$paths
    roots <- unique(c(normalizePath(dirname(dirname(paths$source_dir)), mustWork = FALSE), dirname(dirname(paths$source_dir))))
    probe <- new.env()
    if (managed) omicsManagedChainMocks(fixture, probe)
    log <- new.env(parent = emptyenv())
    log$notes <- list()
    trxNormW1InstallNotifications(log)
    withr::local_seed(20261108L)
    acc <- new.env(parent = emptyenv())
    acc$captures <- list()
    acc$log <- 0L
    acc$operations <- 0L
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(trxNormW1Server(fixture, managed), {
        capture <- function(step) {
            session$flushReact()
            manager <- fixture$workflow$state_manager
            nd <- shiny::isolate(shiny::reactiveValuesToList(session$userData$norm))
            log_lines <- nd$normalization_log
            new_log <- log_lines[seq_along(log_lines) > acc$log]
            acc$log <- length(log_lines)
            notes <- vapply(log$notes, function(note) {
                sprintf("[%s %s] %s", note$type, if (is.null(note$duration)) "NULL" else format(note$duration),
                    guiEntryMask(note$message, roots))
            }, "")
            log$notes <- list()
            qc_files <- sort(list.files(paths$gene_qc_dir))
            pngs <- file.path(paths$gene_qc_dir, qc_files[endsWith(qc_files, ".png")])
            encoded <- stats::setNames(vapply(pngs, base64enc::base64encode, ""), basename(pngs))
            out_names <- sort(names(session$.__enclos_env__$private$outs))
            out_names <- out_names[startsWith(out_names, "norm-")]
            read <- function(name) tryCatch(output[[name]], error = function(e) {
                if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", conditionMessage(e))
            })
            values <- vapply(out_names, function(name) trxNormW1OutputShown(sub("^norm-", "", name), read(name), encoded, roots), "")
            names(values) <- sub("^norm-", "", names(values))
            ruv <- nd$ruv_optimization_results
            operations <- if (managed) probe$operations[seq_along(probe$operations) > acc$operations] else NULL
            acc$operations <- length(probe$operations)
            acc$captures[[step]] <- list(
                outputs = names(values),
                values = values,
                norm = c(normalization_complete = as.character(isTRUE(nd$normalization_complete)),
                    ruv_complete = as.character(isTRUE(nd$ruv_complete)),
                    plot_refresh_trigger = format(nd$plot_refresh_trigger),
                    assay_names = paste(nd$assay_names, collapse = "|"),
                    post_filter_obj = trxNormW1SlotState(nd$post_filter_obj, manager),
                    ruv_corrected_obj = trxNormW1SlotState(nd$ruv_corrected_obj, manager),
                    ruv_optimization_results = paste(names(ruv), collapse = "|"),
                    cancor_plot_class = paste(class(ruv$counts$cancor_plot), collapse = "|"),
                    ruv_result = trxNormW1Sha(guiEntryCanonical(ruv$counts[setdiff(names(ruv$counts), "cancor_plot")]))),
                history = as.character(manager$getHistory()),
                states = trxNormW1Sha(guiEntryStateDigests(manager)),
                outcomes = trxNormW1Sha(guiEntryCanonical(trxNormW1RuvsOutcomes(manager))),
                tab_status = vapply(fixture$workflow$tab_status, as.character, ""),
                notes = notes,
                log = vapply(sub("^\\[[0-9:]+\\] ", "", new_log), guiEntryMask, "", roots = roots, USE.NAMES = FALSE),
                gene_qc_files = qc_files,
                publication_files = sort(list.files(paths$publication_graphs_dir)),
                operations = operations)
        }
        ctx <- list(set = function(...) session$setInputs(...), capture = capture)
        for (step in names(steps)) steps[[step]](ctx)
    }))))
    acc$captures
}

#' The steps of each case.
trxNormW1Steps <- function(case) {
    run <- function(ctx, n, ...) {
        if (...length()) ctx$set(...)
        ctx$set(`norm-run_normalization` = n)
    }
    reset <- function(ctx, n) ctx$set(`norm-reset_normalization` = n)
    switch(case,
        automatic_reset = list(
            run = function(ctx) {
                run(ctx, 1L, `norm-ruv_mode` = "automatic")
                ctx$capture("run")
            },
            reset = function(ctx) {
                reset(ctx, 1L)
                ctx$capture("reset")
            }),
        manual = list(
            run = function(ctx) {
                run(ctx, 1L, `norm-ruv_mode` = "manual", `norm-ruv_k` = 2L, `norm-ruv_percentage` = 10)
                ctx$capture("run")
            }),
        reset_refused = list(
            reset_before_run = function(ctx) {
                reset(ctx, 1L)
                ctx$capture("reset_before_run")
            },
            refused = function(ctx) {
                run(ctx, 1L, `norm-ruv_mode` = "manual", `norm-ruv_k` = 4L, `norm-ruv_percentage` = 10)
                ctx$capture("refused")
            }),
        managed = list(
            run = function(ctx) {
                run(ctx, 1L, `norm-ruv_mode` = "automatic")
                ctx$capture("run")
            }),
        # Its own session: after a failed action the in-test worker (this process) refuses the next
        # one ("native reference roots contain a live payload"), which a real worker process does not.
        managed_refused = list(
            refused = function(ctx) {
                run(ctx, 1L, `norm-ruv_mode` = "manual", `norm-ruv_k` = 4L, `norm-ruv_percentage` = 10)
                ctx$capture("refused")
            }))
}

trxNormW1Runs <- new.env(parent = emptyenv())

trxNormW1Case <- function(case) {
    if (is.null(trxNormW1Runs[[case]])) {
        testthat::skip_if_not_installed("edgeR")
        trxNormW1Runs[[case]] <- trxNormW1Drive(trxNormW1Steps(case), managed = startsWith(case, "managed"))
        if (!is.null(trxNormW1DumpPath())) {
            cat(sprintf("trxNormW1Frozen[[\"%s\"]] <- %s\n\n", case, trxNormW1Code(trxNormW1Runs[[case]])),
                file = trxNormW1DumpPath(), append = TRUE)
        }
    }
    trxNormW1Runs[[case]]
}

# The frozen captures, per case (the dump's entries, pasted).
# TRX_NORM_W1_FROZEN_START
trxNormW1Frozen <- list()
trxNormW1Frozen[["automatic_reset"]] <- list(
    run = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "cancor_plot_counts", "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui",
            "ruv_summary_counts", "ruv_table_counts"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)", cancor_plot_counts = "<plot> alt=Plot object",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Running TMM, the RUV grid (automatic) and RUVs\n[<TIME>] RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')\n[<TIME>] PCA and RLE (raw) generated\n[<TIME>] PCA and RLE (tmm) generated\n[<TIME>] PCA and RLE (ruvs) generated\n[<TIME>] CCA panels, grid workbook and the normalisation comparison grid saved\n[<TIME>] Normalization complete - ready for differential expression",
            pca_post_filter_assay1 = "trx_raw_pca.png alt=pca - counts",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "trx_tmm_pca.png alt=pca - counts",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "trx_ruvs_pca.png alt=pca - counts",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "trx_raw_rle.png alt=rle - counts",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "trx_tmm_rle.png alt=rle - counts",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "trx_ruvs_rle.png alt=rle - counts",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558",
            ruv_summary_counts = "Best k: 2\nBest %: 10\nSeparation: 0.9892\nControls: 18",
            ruv_table_counts = "<text> 0eda554b452587f838480abb09185d552be1615fa89824563914bbc78b7e4f0c"),
        norm = c(
            normalization_complete = "TRUE", ruv_complete = "TRUE", plot_refresh_trigger = "1",
            assay_names = "counts", post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "trx_ruvs",
            ruv_optimization_results = "counts",
            cancor_plot_class = "ggplot2::ggplot|ggplot|ggplot2::gg|S7_object|gg",
            ruv_result = "7ef0b152b466d471fd5fcfe5b521d5193e8e114aba5e9d89615b644f54d10355"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"),
        states = "bbcc8e91c9921cafbad19de435d2542e4dcbcaccf24addea238253fb56938ddc",
        outcomes = "ed69a6f33d5780c4aa83c1086c9b4f288b616ccc7c6fc8c261ebfea561d99ba1",
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = "[message 5] Normalization complete! Proceeding to DE analysis.",
        log = c(
            "Running TMM, the RUV grid (automatic) and RUVs",
            "RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')", "PCA and RLE (raw) generated",
            "PCA and RLE (tmm) generated", "PCA and RLE (ruvs) generated",
            "CCA panels, grid workbook and the normalisation comparison grid saved",
            "Normalization complete - ready for differential expression"),
        gene_qc_files = c(
            "ruv_normalised_log_values.tsv", "trx_raw_pca.png", "trx_raw_rle.png", "trx_ruvs_pca.png",
            "trx_ruvs_rle.png", "trx_tmm_pca.png", "trx_tmm_rle.png"),
        publication_files = c(
            "canonical_correlation_optimization_summary.xlsx", "canonical_correlation_panels.pdf",
            "canonical_correlation_panels.png", "filtering_qc", "normalization_comparison_grid.pdf",
            "normalization_comparison_grid.png"),
        operations = NULL),
    reset = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "cancor_plot_counts", "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui",
            "ruv_summary_counts", "ruv_table_counts"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)", cancor_plot_counts = "<plot> alt=Plot object",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Running TMM, the RUV grid (automatic) and RUVs\n[<TIME>] RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')\n[<TIME>] PCA and RLE (raw) generated\n[<TIME>] PCA and RLE (tmm) generated\n[<TIME>] PCA and RLE (ruvs) generated\n[<TIME>] CCA panels, grid workbook and the normalisation comparison grid saved\n[<TIME>] Normalization complete - ready for differential expression\n[<TIME>] Reset to pre-normalization state",
            pca_post_filter_assay1 = "trx_raw_pca.png alt=pca - counts",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "trx_tmm_pca.png alt=pca - counts",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "trx_ruvs_pca.png alt=pca - counts",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "trx_raw_rle.png alt=rle - counts",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "trx_tmm_rle.png alt=rle - counts",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "trx_ruvs_rle.png alt=rle - counts",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558",
            ruv_summary_counts = "Not yet computed",
            ruv_table_counts = "{\"x\":null,\"evals\":[],\"jsHooks\":[],\"deps\":[{\"name\":\"jquery\",\"version\":\"3.6.0\",\"src\":{\"href\":\"jquery-3.6.0\"},\"meta\":null,\"script\":\"jquery-3.6.0.min.js\",\"stylesheet\":null,\"head\":null,\"attachment\":null,\"all_files\":true},{\"name\":\"dt-core\",\"version\":\"1.13.6\",\"src\":{\"href\":\"dt-core-1.13.6\"},\"meta\":null,\"script\":\"js/jquery.dataTables.min.js\",\"stylesheet\":[\"css/jquery.dataTables.min.css\",\"css/jquery.dataTables.extra.css\"],\"head\":null,\"attachment\":null,\"package\":null,\"all_files\":false},{\"name\":\"crosstalk\",\"version\":\"1.2.2\",\"src\":{\"href\":\"crosstalk-1.2.2\"},\"meta\":null,\"script\":\"js/crosstalk.min.js\",\"stylesheet\":\"css/crosstalk.min.css\",\"head\":null,\"attachment\":null,\"all_files\":true}]}"),
        norm = c(
            normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "1",
            assay_names = "counts", post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "<NULL>",
            ruv_optimization_results = "", cancor_plot_class = "NULL",
            ruv_result = "1e8d93415411c9caecc3ce4932c05f26923726fcdab89e7c8f63ea3d4f3720fc"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs",
            "trx_reset"),
        states = "a87303d8ac9111423e1c2e38ab88d0a7e4369eb9fef17261f6da47da0a3822ab",
        outcomes = "ed69a6f33d5780c4aa83c1086c9b4f288b616ccc7c6fc8c261ebfea561d99ba1",
        tab_status = c(
            quality_control = "complete", normalization = "pending", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = "[message 5] Reset to pre-normalization state",
        log = "Reset to pre-normalization state",
        gene_qc_files = c(
            "ruv_normalised_log_values.tsv", "trx_raw_pca.png", "trx_raw_rle.png", "trx_ruvs_pca.png",
            "trx_ruvs_rle.png", "trx_tmm_pca.png", "trx_tmm_rle.png"),
        publication_files = c(
            "canonical_correlation_optimization_summary.xlsx", "canonical_correlation_panels.pdf",
            "canonical_correlation_panels.png", "filtering_qc", "normalization_comparison_grid.pdf",
            "normalization_comparison_grid.png"),
        operations = NULL))

trxNormW1Frozen[["manual"]] <- list(
    run = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "cancor_plot_counts", "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui",
            "ruv_summary_counts", "ruv_table_counts"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)", cancor_plot_counts = "<plot> alt=Plot object",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Running TMM, the RUV grid (manual) and RUVs\n[<TIME>] RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')\n[<TIME>] PCA and RLE (raw) generated\n[<TIME>] PCA and RLE (tmm) generated\n[<TIME>] PCA and RLE (ruvs) generated\n[<TIME>] CCA panels, grid workbook and the normalisation comparison grid saved\n[<TIME>] Normalization complete - ready for differential expression",
            pca_post_filter_assay1 = "trx_raw_pca.png alt=pca - counts",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "trx_tmm_pca.png alt=pca - counts",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "trx_ruvs_pca.png alt=pca - counts",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "trx_raw_rle.png alt=rle - counts",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "trx_tmm_rle.png alt=rle - counts",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "trx_ruvs_rle.png alt=rle - counts",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558",
            ruv_summary_counts = "Best k: 2\nBest %: 10\nSeparation: 0.9892\nControls: 18",
            ruv_table_counts = "<text> 0eda554b452587f838480abb09185d552be1615fa89824563914bbc78b7e4f0c"),
        norm = c(
            normalization_complete = "TRUE", ruv_complete = "TRUE", plot_refresh_trigger = "1",
            assay_names = "counts", post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "trx_ruvs",
            ruv_optimization_results = "counts",
            cancor_plot_class = "ggplot2::ggplot|ggplot|ggplot2::gg|S7_object|gg",
            ruv_result = "8f0af161e8e22ab7ea0a8a34e35ad2f487f55c8e71518273551f35927351cd78"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"),
        states = "bbcc8e91c9921cafbad19de435d2542e4dcbcaccf24addea238253fb56938ddc",
        outcomes = "ec85beb9d1dd555a0e788552abd9eb04c487cb5961390f322d856cb612783364",
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = "[message 5] Normalization complete! Proceeding to DE analysis.",
        log = c(
            "Running TMM, the RUV grid (manual) and RUVs",
            "RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')", "PCA and RLE (raw) generated",
            "PCA and RLE (tmm) generated", "PCA and RLE (ruvs) generated",
            "CCA panels, grid workbook and the normalisation comparison grid saved",
            "Normalization complete - ready for differential expression"),
        gene_qc_files = c(
            "ruv_normalised_log_values.tsv", "trx_raw_pca.png", "trx_raw_rle.png", "trx_ruvs_pca.png",
            "trx_ruvs_rle.png", "trx_tmm_pca.png", "trx_tmm_rle.png"),
        publication_files = c(
            "canonical_correlation_optimization_summary.xlsx", "canonical_correlation_panels.pdf",
            "canonical_correlation_panels.png", "filtering_qc", "normalization_comparison_grid.pdf",
            "normalization_comparison_grid.png"),
        operations = NULL))

trxNormW1Frozen[["reset_refused"]] <- list(
    reset_before_run = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Reset to pre-normalization state",
            pca_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_pca.png",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_pca.png",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_pca.png",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_rle.png",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_rle.png",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_rle.png",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558"),
        norm = c(
            normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "0",
            assay_names = "counts", post_filter_obj = "<NULL>", ruv_corrected_obj = "<NULL>",
            ruv_optimization_results = "", cancor_plot_class = "NULL",
            ruv_result = "1e8d93415411c9caecc3ce4932c05f26923726fcdab89e7c8f63ea3d4f3720fc"),
        history = c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete"),
        states = "d6b0c4f539ac28e638197e985adea5895477821122c9dee997f57ebc72f6859c",
        outcomes = "743a1d2c832c72b269ed1dc4181d24657bf9c7296c6f9b1c4934328d5583605e",
        tab_status = c(
            quality_control = "complete", normalization = "pending", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = "[message 5] Reset to pre-normalization state",
        log = "Reset to pre-normalization state",
        gene_qc_files = character(0),
        publication_files = "filtering_qc",
        operations = NULL),
    refused = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Reset to pre-normalization state\n[<TIME>] Running TMM, the RUV grid (manual) and RUVs\n[<TIME>] Normalization error: manual k 4 is above the cap of 3 for 12 samples (trx_ruv_k_above_cap)",
            pca_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_pca.png",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_pca.png",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_pca.png",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_rle.png",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_rle.png",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_rle.png",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558"),
        norm = c(
            normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "0",
            assay_names = "counts", post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "<NULL>",
            ruv_optimization_results = "", cancor_plot_class = "NULL",
            ruv_result = "1e8d93415411c9caecc3ce4932c05f26923726fcdab89e7c8f63ea3d4f3720fc"),
        history = c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete"),
        states = "d6b0c4f539ac28e638197e985adea5895477821122c9dee997f57ebc72f6859c",
        outcomes = "743a1d2c832c72b269ed1dc4181d24657bf9c7296c6f9b1c4934328d5583605e",
        tab_status = c(
            quality_control = "complete", normalization = "pending", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = "[error 15] Normalization error: manual k 4 is above the cap of 3 for 12 samples (trx_ruv_k_above_cap)",
        log = c(
            "Running TMM, the RUV grid (manual) and RUVs",
            "Normalization error: manual k 4 is above the cap of 3 for 12 samples (trx_ruv_k_above_cap)"),
        gene_qc_files = character(0),
        publication_files = "filtering_qc",
        operations = NULL))

trxNormW1Frozen[["managed"]] <- list(
    run = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "cancor_plot_counts", "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui",
            "ruv_summary_counts", "ruv_table_counts"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)", cancor_plot_counts = "<plot> alt=Plot object",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Running TMM, the RUV grid (automatic) and RUVs\n[<TIME>] PCA and RLE (raw) generated\n[<TIME>] TMM factors and logCPM computed\n[<TIME>] RUV grid (automatic): k = 2, 10% inferred negative controls\n[<TIME>] RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')\n[<TIME>] PCA and RLE (tmm) generated\n[<TIME>] PCA and RLE (ruvs) generated\n[<TIME>] CCA panels, grid workbook and the normalisation comparison grid saved\n[<TIME>] Normalization complete - ready for differential expression",
            pca_post_filter_assay1 = "trx_raw_pca.png alt=pca - counts",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "trx_tmm_pca.png alt=pca - counts",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "trx_ruvs_pca.png alt=pca - counts",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "trx_raw_rle.png alt=rle - counts",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "trx_tmm_rle.png alt=rle - counts",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "trx_ruvs_rle.png alt=rle - counts",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558",
            ruv_summary_counts = "Best k: 2\nBest %: 10\nSeparation: 0.9892\nControls: 18",
            ruv_table_counts = "<text> 0eda554b452587f838480abb09185d552be1615fa89824563914bbc78b7e4f0c"),
        norm = c(
            normalization_complete = "TRUE", ruv_complete = "TRUE", plot_refresh_trigger = "1",
            assay_names = "counts", post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "trx_ruvs",
            ruv_optimization_results = "counts",
            cancor_plot_class = "ggplot2::ggplot|ggplot|ggplot2::gg|S7_object|gg",
            ruv_result = "7ef0b152b466d471fd5fcfe5b521d5193e8e114aba5e9d89615b644f54d10355"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"),
        states = "bbcc8e91c9921cafbad19de435d2542e4dcbcaccf24addea238253fb56938ddc",
        outcomes = "ed69a6f33d5780c4aa83c1086c9b4f288b616ccc7c6fc8c261ebfea561d99ba1",
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = c(
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[message 5] Normalization complete! Proceeding to DE analysis."),
        log = c(
            "Running TMM, the RUV grid (automatic) and RUVs", "PCA and RLE (raw) generated",
            "TMM factors and logCPM computed", "RUV grid (automatic): k = 2, 10% inferred negative controls",
            "RUVs applied: W_1..W_2 added to the design (state 'trx_ruvs')", "PCA and RLE (tmm) generated",
            "PCA and RLE (ruvs) generated",
            "CCA panels, grid workbook and the normalisation comparison grid saved",
            "Normalization complete - ready for differential expression"),
        gene_qc_files = c(
            "ruv_normalised_log_values.tsv", "trx_raw_pca.png", "trx_raw_rle.png", "trx_ruvs_pca.png",
            "trx_ruvs_rle.png", "trx_tmm_pca.png", "trx_tmm_rle.png"),
        publication_files = c(
            "canonical_correlation_optimization_summary.xlsx", "canonical_correlation_panels.pdf",
            "canonical_correlation_panels.png", "filtering_qc", "normalization_comparison_grid.pdf",
            "normalization_comparison_grid.png"),
        operations = c(
            "trx_qc_figures", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_qc_figures",
            "trx_qc_figures", "trx_ruv_publication", "trx_composite_figure")))

trxNormW1Frozen[["managed_refused"]] <- list(
    refused = list(
        outputs = c(
            "assay1_label_correlation", "assay1_label_density", "assay1_label_pca", "assay1_label_rle",
            "assay2_label_correlation", "assay2_label_density", "assay2_label_pca", "assay2_label_rle",
            "correlation_post_filter_assay1", "correlation_post_filter_assay2",
            "correlation_post_norm_assay1", "correlation_post_norm_assay2",
            "correlation_ruv_corrected_assay1", "correlation_ruv_corrected_assay2",
            "density_post_filter_assay1", "density_post_filter_assay2", "density_post_norm_assay1",
            "density_post_norm_assay2", "density_ruv_corrected_assay1", "density_ruv_corrected_assay2",
            "norm_log", "pca_post_filter_assay1", "pca_post_filter_assay2", "pca_post_norm_assay1",
            "pca_post_norm_assay2", "pca_ruv_corrected_assay1", "pca_ruv_corrected_assay2",
            "rle_post_filter_assay1", "rle_post_filter_assay2", "rle_post_norm_assay1",
            "rle_post_norm_assay2", "rle_ruv_corrected_assay1", "rle_ruv_corrected_assay2", "ruv_qc_ui"),
        values = c(
            assay1_label_correlation = "Assay: counts", assay1_label_density = "Assay: counts",
            assay1_label_pca = "Assay: counts", assay1_label_rle = "Assay: counts",
            assay2_label_correlation = "Assay 2: (detecting...)",
            assay2_label_density = "Assay 2: (detecting...)", assay2_label_pca = "Assay 2: (detecting...)",
            assay2_label_rle = "Assay 2: (detecting...)",
            correlation_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_correlation.png",
            correlation_post_filter_assay2 = "<none> Assay not detected yet",
            correlation_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_correlation.png",
            correlation_post_norm_assay2 = "<none> Assay not detected yet",
            correlation_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_correlation.png",
            correlation_ruv_corrected_assay2 = "<none> Assay not detected yet",
            density_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_density.png",
            density_post_filter_assay2 = "<none> Assay not detected yet",
            density_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_density.png",
            density_post_norm_assay2 = "<none> Assay not detected yet",
            density_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_density.png",
            density_ruv_corrected_assay2 = "<none> Assay not detected yet",
            norm_log = "[<TIME>] Running TMM, the RUV grid (manual) and RUVs\n[<TIME>] PCA and RLE (raw) generated\n[<TIME>] TMM factors and logCPM computed\n[<TIME>] Normalization error: worker failed: manual k 4 is above the cap of 3 for 12 samples (simpleError)",
            pca_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_pca.png",
            pca_post_filter_assay2 = "<none> Assay not detected yet",
            pca_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_pca.png",
            pca_post_norm_assay2 = "<none> Assay not detected yet",
            pca_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_pca.png",
            pca_ruv_corrected_assay2 = "<none> Assay not detected yet",
            rle_post_filter_assay1 = "<none> Plot not generated yet: trx_raw_rle.png",
            rle_post_filter_assay2 = "<none> Assay not detected yet",
            rle_post_norm_assay1 = "<none> Plot not generated yet: trx_tmm_rle.png",
            rle_post_norm_assay2 = "<none> Assay not detected yet",
            rle_ruv_corrected_assay1 = "<none> Plot not generated yet: trx_ruvs_rle.png",
            rle_ruv_corrected_assay2 = "<none> Assay not detected yet",
            ruv_qc_ui = "<ui> 26fd2510633d7e5170d56cbcca2c6655538f7c1b659d4b7e17b3b86cb7f6c558"),
        norm = c(
            normalization_complete = "FALSE", ruv_complete = "FALSE", plot_refresh_trigger = "0",
            assay_names = "counts", post_filter_obj = "trx_qc_complete", ruv_corrected_obj = "<NULL>",
            ruv_optimization_results = "", cancor_plot_class = "NULL",
            ruv_result = "1e8d93415411c9caecc3ce4932c05f26923726fcdab89e7c8f63ea3d4f3720fc"),
        history = c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete"),
        states = "d6b0c4f539ac28e638197e985adea5895477821122c9dee997f57ebc72f6859c",
        outcomes = "743a1d2c832c72b269ed1dc4181d24657bf9c7296c6f9b1c4934328d5583605e",
        tab_status = c(
            quality_control = "complete", normalization = "pending", differential_analysis = "disabled",
            session_summary = "disabled"),
        notes = c(
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[default NULL] Applying transcriptomics step ...",
            "[error 15] Error applying transcriptomics step: worker failed: manual k 4 is above the cap of 3 for 12 samples",
            "[error 15] Normalization error: worker failed: manual k 4 is above the cap of 3 for 12 samples (simpleError)"),
        log = c(
            "Running TMM, the RUV grid (manual) and RUVs", "PCA and RLE (raw) generated",
            "TMM factors and logCPM computed",
            "Normalization error: worker failed: manual k 4 is above the cap of 3 for 12 samples (simpleError)"),
        gene_qc_files = c("trx_raw_pca.png", "trx_raw_rle.png"),
        publication_files = "filtering_qc",
        operations = c("trx_qc_figures", "trx_tmm_logcpm", "trx_ruv_grid")))
# TRX_NORM_W1_FROZEN_END

expectTrxNormW1Case <- function(case) {
    seen <- trxNormW1Case(case)
    frozen <- trxNormW1Frozen[[case]]
    testthat::expect_identical(names(seen), names(frozen), info = case)
    for (step in names(frozen)) {
        for (field in names(frozen[[step]])) {
            testthat::expect_identical(seen[[step]][[field]], frozen[[step]][[field]], info = paste(case, step, field))
        }
    }
}

test_that("(b-trx) in process: an automatic run, then a reset after it", {
    expectTrxNormW1Case("automatic_reset")
})

test_that("(b-trx) in process: a manual run (k 2, 10 %)", {
    expectTrxNormW1Case("manual")
})

test_that("(b-trx) in process: a reset before any run, then a manual k above the cap (refused)", {
    expectTrxNormW1Case("reset_refused")
})

test_that("(b-trx) managed: an automatic run on the worker, its operations in order", {
    expectTrxNormW1Case("managed")
})

test_that("(b-trx) managed: a manual k above the cap, refused by the worker's grid", {
    expectTrxNormW1Case("managed_refused")
})
