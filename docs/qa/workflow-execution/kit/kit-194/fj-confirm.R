# KIT-194 gate G0 (ticket section 8): F-j executed on the base. For lipidomics and metabolomics, in process, on fixture
# F172: the tab's run (`runOmicNormPipeline`: ITSD by pattern, offset 0, cyclic loess, RUV skip), then one route of 2.2,
# then "Run DA" with no state loaded (`runOmicDaAnalysis` with the real `executeOmicDaRun` and in-process DA;
# `completeOmicDaRun` records the results). Routes: J1 Revert to `<prefix>_qc_complete` through the Intensity Filter
# sub-tab's real observer, then "Apply Correlation Filter" (threshold 0.5); J2 the same with Skip; J3 Apply, then the
# Revert, then the DA; J4 Apply, then "Reset to Pre-Normalization", then the DA. Writes `fj-confirm.tsv` in <out dir>.
# Stop rules (section 8): J1 not reproduced on both omics stops 7.4, T6, T7; neither J3 nor J4 reproduced stops 7.5, T8.
#
# Usage (a run-tests `scripts` entry, from the package root after R CMD INSTALL):
#   Rscript docs/qa/workflow-execution/kit/kit-194/fj-confirm.R <out dir>

suppressPackageStartupMessages(library(testthat))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[[1L]] else "kit-194-g0"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g0 <- new.env(parent = asNamespace("MultiScholaR"))
for (helper in c("helper-workflow-omics-owner.R", "helper-workflow-lipid-declarations.R", "helper-workflow-metab-declarations.R",
    "helper-kit172.R")) {
    sys.source(file.path("tests", "testthat", helper), envir = g0)
}
g0$out_dir <- out_dir

evalq({
    quiet <- function(expr) suppressWarnings(suppressMessages(expr))

    fixture <- function(root, class) {
        out <- kit172Owner(root, kit172Object(kit172Values(172001L)$values, class), class)
        binding <- kit172Binding(class)
        quiet({
            runWorkflowOmicsReferenceIntensityFilter(binding, out$owner, 1, 0.5, out$workflow$config_list)
            commitWorkflowOmicsStoreStep(out$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
                out$workflow$config_list, "QC complete")
        })
        out$class <- class
        out$omic <- if (identical(class, "lipid")) "lipidomics" else "metabolomics"
        out$binding <- binding
        out$prefix <- binding$prefix
        config <- omicGuiConfig(out$omic)
        config$binding$binding$qc_plots <- function(...) invisible(NULL)
        out$config <- config
        out$log <- new.env(parent = emptyenv())
        out$log$lines <- character()
        log <- out$log
        out$addLog <- function(entry) log$lines <- c(log$lines, entry)
        out
    }

    input <- list(apply_itsd = TRUE, itsd_aggregation = "median", log_offset = 0, norm_method = "cyclicloess", ruv_mode = "skip",
        auto_percentage_min = 1, auto_percentage_max = 20, ruv_grouping_variable = "group", separation_metric = "max_difference",
        k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_k = 1L, ruv_percentage = 20, min_pearson_correlation_threshold = 0.5,
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)

    owner <- function(fx) fx$workflow$state_manager

    run <- function(fx, nd) quiet(shiny::isolate(runOmicNormPipeline(input, fx$workflow, list(), nd, fx$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fx$config)))

    revert <- function(fx, until) {
        quiet(utils::capture.output(shiny::testServer(mod_omic_qc_intensity_server, args = list(workflow_data = fx$workflow,
            omic_type = fx$omic, experiment_label = "kit194", config = fx$config), {
            press <- 0L
            while (press < 10L && !identical(fx$workflow$state_manager$current_state, until)) {
                press <- press + 1L
                session$setInputs(revert_filter = press)
            }
        })))
        owner(fx)$current_state
    }

    lcms <- function(fx, reference) {
        object <- readWorkflowOmicsState(fx$binding, owner(fx), reference)
        table <- as.data.frame(methods::slot(object, fx$binding$data_slot)$LCMS_Pos)
        id <- methods::slot(object, fx$binding$id_slot)
        runs <- intersect(kit172Design()$Run, colnames(table))
        matrix(as.matrix(table[, runs, drop = FALSE]), ncol = length(runs), dimnames = list(table[[id]], runs))
    }

    names_of <- function(fx, reference) {
        states <- owner(fx)$.context$snapshot$states
        paste(names(Filter(function(r) identical(r$entity_id, reference$entity_id), states)), collapse = "|")
    }

    # The DA as the tab runs it with no state loaded: the real run, its results recorded.
    da <- function(fx, seen) {
        recorded <- new.env(parent = emptyenv())
        config <- fx$config$da
        in_process <- config$run$in_process
        config$run$in_process <- function(manager, source, settings, contrasts) {
            recorded$source <- source
            in_process(manager, source, settings, contrasts)
        }
        testthat::with_mocked_bindings(
            quiet(shiny::isolate(runOmicDaAnalysis(input, NULL, fx$workflow, new.env(parent = emptyenv()), list(), config, fx$binding))),
            completeOmicDaRun = function(results, ...) recorded$results <- results,
            .package = "MultiScholaR")
        recorded
    }

    route <- function(root, class, name) {
        fx <- fixture(file.path(root, paste(class, name, sep = "-")), class)
        on.exit(try(closeDIAStateAdapter(owner(fx)), silent = TRUE), add = TRUE)
        prefix <- fx$prefix
        seen <- new.env(parent = emptyenv())
        seen$notes <- character()
        row <- list(omic = fx$omic, route = name)
        testthat::with_mocked_bindings({
            nd <- omicNormReactiveState()
            run(fx, nd)
            row$precondition <- identical(owner(fx)$current_state, paste0(prefix, "_normalized")) &&
                isTRUE(shiny::isolate(nd$normalization_complete))
            normalised_max <- max(lcms(fx, paste0(prefix, "_normalized")), na.rm = TRUE)
            qc <- lcms(fx, paste0(prefix, "_qc_complete"))
            apply <- function() quiet(shiny::isolate(applyOmicNormCorrelation(input, fx$workflow, nd, fx$addLog, fx$config)))
            if (name %in% c("J1", "J2")) {
                row$revert_head <- revert(fx, paste0(prefix, "_qc_complete"))
                if (identical(name, "J1")) apply() else quiet(shiny::isolate(skipOmicNormCorrelation(fx$workflow, nd, fx$addLog, fx$config)))
            } else {
                apply()
                if (identical(name, "J3")) row$revert_head <- revert(fx, paste0(prefix, "_qc_complete")) else {
                    quiet(shiny::isolate(resetOmicNorm(fx$workflow, nd, fx$addLog, fx$config)))
                    row$revert_head <- NA_character_
                }
            }
            head <- getWorkflowStateReference(owner(fx))
            spec <- readWorkflowStateMetadata(owner(fx), head)$data
            matrix <- lcms(fx, head)
            row$head_state <- owner(fx)$current_state
            row$head_products <- paste(sort(names(spec$products)), collapse = ",")
            row$head_dependency_labels <- paste(sort(unique(unlist(lapply(spec$dependencies, function(d)
                if (is.list(d)) d$identity$parameters)))), collapse = ",")
            row$head_max <- max(matrix, na.rm = TRUE)
            row$qc_complete_max <- max(qc, na.rm = TRUE)
            row$normalised_max <- normalised_max
            rows <- intersect(rownames(matrix), rownames(qc))
            runs <- intersect(colnames(matrix), colnames(qc))
            row$head_equals_qc_complete <- length(rows) == nrow(matrix) && identical(unname(matrix[rows, runs, drop = FALSE]),
                unname(qc[rows, runs, drop = FALSE]))
            row$normalization_status <- as.character(fx$workflow$tab_status$normalization %||% NA_character_)
            result <- da(fx, seen)
            row$da_source_state <- if (is.null(result$source)) NA_character_ else names_of(fx, result$source)
            long <- result$results[[fx$binding$long_table]]
            row$da_rows <- if (is.null(long)) 0L else nrow(long)
            row$da_input_max <- if (is.null(result$source)) NA_real_ else max(lcms(fx, result$source), na.rm = TRUE)
        }, showNotification = function(ui, ..., type = "default") {
            seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        }, removeNotification = function(...) invisible(NULL), withProgress = function(expr, ...) expr,
        incProgress = function(...) invisible(NULL), updateSelectInput = function(...) invisible(NULL), .package = "shiny")
        row$notes <- paste(gsub("[\t\n]", " ", seen$notes), collapse = " | ")
        row$log_errors <- paste(grep("^ERROR", fx$log$lines, value = TRUE), collapse = " | ")
        row
    }

    root <- tempfile("kit194-g0-")
    dir.create(root)
    rows <- list()
    for (class in c("lipid", "metab")) {
        for (name in c("J1", "J2", "J3", "J4")) {
            row <- tryCatch(route(root, class, name), error = function(e) list(omic = class, route = name,
                error = gsub("[\t\n]", " ", conditionMessage(e))))
            rows[[length(rows) + 1L]] <- row
            cat(sprintf("KIT-194 G0\t%s\t%s\t%s\n", class, name, paste(names(row), unlist(lapply(row, function(v) paste(v, collapse = ","))),
                sep = "=", collapse = "\t")))
        }
    }
    fields <- unique(unlist(lapply(rows, names)))
    table <- do.call(rbind, lapply(rows, function(row) as.data.frame(stats::setNames(lapply(fields, function(f)
        if (is.null(row[[f]])) NA_character_ else paste(row[[f]], collapse = ",")), fields), stringsAsFactors = FALSE)))
    utils::write.table(table, file.path(out_dir, "fj-confirm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("KIT-194 G0: wrote", file.path(out_dir, "fj-confirm.tsv"), "\n")
}, g0)
