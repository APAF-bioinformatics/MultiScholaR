# KIT-259 gate G0-B (ticket r1.1 section 10): does a proteomics normalisation figure that cannot be drawn leave the
# normalisation committed but unfinished, at the default warn, on base? Observation only: no package code changes.
# Binding: v0.9-staging c4505c7f plus this file (the branch point of kit-259).
#
# Part IP (in process; KIT-247's guiEntryRunLane pattern): the label-free lane prot_lfq, its steps up to the
#   normalisation's refusal, then (a) the choices 'Color by' and 'Shape by' offer and the design's columns; (b) a 'Color
#   by' among them with no within-group pair, Run Normalization, then Apply correlation, Skip, DA entry, Reset, Run again,
#   and 'Color by' = group then Run again.
# Part IPD (in process, supplementary): the same (b) sequence on KIT-122's DIA owner with an unset batch, through
#   testServer(mod_prot_norm_server), with Run again before Reset.
# Part M (managed; KIT-122's in-test harness, protNormW3Managed): a lifecycle check first (one managed normalisation
#   reaching normalization_complete with its notice; one worker failure reaching the tab's failed handler with its
#   notice), then the counted (b) rows on the DIA owner with an unset batch.
# Part C: the trigger A2 does not remove: 16 groups x 2 samples with the default 'Shape by' (group), through
#   generateProtNormPostNormalizationQcArtifacts on the S4 object, and through runDIAReferenceNormalization (the tab, in
#   process) after tab entry; records whether the pre-normalisation PCA failed at entry and redraw_pre.
#
# Usage (package root, on a runner, NOT_CRAN=true): Rscript <this file> <out dir>

k259Args <- commandArgs(trailingOnly = TRUE)
k259Out <- if (length(k259Args) >= 1L) k259Args[[1L]] else "kit259-g0b"
dir.create(k259Out, recursive = TRUE, showWarnings = FALSE)
k259Out <- normalizePath(k259Out)

# ---------------------------------------------------------------------------------------------------------------------
# PREDICTIONS (written before any value is read; regular expressions applied to the observed value), by reading at
# c4505c7f (ticket section 2.3). An unset batch in 'Color by' gives zero within-group pairs and the correlation figure
# fails at seq(Inf, ...) ("'from' must be a finite number"); in process the state is committed before the figures, so the
# head is 'normalized', normalization_complete stays FALSE, the publication is pending and every later step is refused
# ("normalization outputs are pending; retry that action"); Reset recovers; another 'Color by' then Run again completes.
# Managed: the action commits on acceptance, the figures action fails, publish never runs (nothing pending), Apply is
# refused by diaCorrelationSource. 16 groups with the default 'Shape by' fail the PCA at ggsave ("Insufficient values in
# manual scale"), at tab entry first, and after the commit at Run Normalization (redraw_pre FALSE on a first run).
# ---------------------------------------------------------------------------------------------------------------------
k259P <- function(key, predicted, pattern) data.frame(key = key, predicted = predicted, pattern = pattern,
    stringsAsFactors = FALSE)
k259Predictions <- rbind(
    # (a)
    k259P("ip.a.color_choices", "the builder's design offers group, factor1, factor2 and batch", "^group,factor1,factor2,batch$"),
    k259P("ip.a.shape_choices", "the same for 'Shape by'", "^group,factor1,factor2,batch$"),
    k259P("ip.a.sample_id_offered", "sample_id is not offered (the builder names it Run)", "^FALSE$"),
    k259P("ip.a.design_columns", "Run ... tech_rep_group, no sample_id column", "^Run,group,.*batch.*tech_rep_group$"),
    k259P("ip.a.zero_pair_variable", "the unset batch gives no within-group pair", "^batch$"),
    k259P("ip.a.warn", "the lane runs at warn = 1 (no warning promoted)", "^1$"),
    # (b) in process, the lane
    k259P("ip.b_color.events", "the pre-normalisation correlation figure fails on regeneration", "plot correlation \\[batch/group\\] = error: .*finite"),
    k259P("ip.b_run.head", "the run commits before the figures", "^normalized$"),
    k259P("ip.b_run.normalization_complete", "the step is not finished", "^FALSE$"),
    k259P("ip.b_run.ruv_normalized_obj", "set before the figures are drawn", "^set$"),
    k259P("ip.b_run.pending", "the publication is left pending", "^normalization:FALSE$"),
    k259P("ip.b_run.notes", "the tab's error notice", "Error in normalization: .*finite"),
    k259P("ip.b_run.events", "the post-normalisation correlation figure is the one that fails", "plot correlation \\[batch/group\\] = error"),
    k259P("ip.b_apply.head", "nothing moves", "^normalized$"),
    k259P("ip.b_apply.notes", "Apply correlation is refused", "normalization outputs are pending; retry that action"),
    k259P("ip.b_skip.notes", "Skip is refused", "normalization outputs are pending; retry that action"),
    k259P("ip.b_da.notes", "DA is refused", "outputs are pending"),
    k259P("ip.b_reset.head", "Reset reverts to the pre-normalisation table", "^protein_replicate_filtered$"),
    k259P("ip.b_reset.notes", "with its notice", "Normalisation reset"),
    k259P("ip.b_rerun.head", "Run again (same 'Color by') commits again", "^normalized$"),
    k259P("ip.b_rerun.normalization_complete", "and fails the same way", "^FALSE$"),
    k259P("ip.b_rerun.notes", "with the same notice", "Error in normalization: .*finite"),
    k259P("ip.b_group_rerun.head", "'Color by' = group, Run again", "^normalized$"),
    k259P("ip.b_group_rerun.normalization_complete", "completes (a new identity)", "^TRUE$"),
    k259P("ip.b_group_rerun.notes", "with the completion notice", "Normalization completed \\(RUV skipped\\)"),
    # (b) in process on the DIA owner (supplementary)
    k259P("ipd.b_run.head", "commits before the figures", "^normalized$"),
    k259P("ipd.b_run.normalization_complete", "not finished", "^FALSE$"),
    k259P("ipd.b_run.pending", "publication pending", "^normalization:FALSE$"),
    k259P("ipd.b_run.notes", "error notice", "Error in normalization: .*finite"),
    k259P("ipd.b_apply.notes", "Apply refused", "normalization outputs are pending; retry that action"),
    k259P("ipd.b_skip.notes", "Skip refused", "normalization outputs are pending; retry that action"),
    k259P("ipd.b_rerun.notes", "Run again from the stranded state fails the same way", "Error in normalization: .*finite"),
    k259P("ipd.b_rerun.normalization_complete", "still unfinished", "^FALSE$"),
    k259P("ipd.b_reset.head", "Reset recovers", "^protein_replicate_filtered$"),
    k259P("ipd.b_group_rerun.normalization_complete", "'Color by' = group then Run completes", "^TRUE$"),
    # Managed: lifecycle, then (b)
    k259P("m.lifecycle.l_complete.normalization_complete", "a managed run completes", "^TRUE$"),
    k259P("m.lifecycle.l_complete.notes", "with its completion notice", "Normalization completed \\(RUV skipped\\)"),
    k259P("m.lifecycle.l_complete.actions", "on the worker", "normalization-ruv"),
    k259P("m.lifecycle.l_failed.notes", "a worker failure reaches the tab's failed handler", "\\[error 10\\] Error in normalization: "),
    k259P("m.lifecycle.ok", "the lifecycle check passes", "^TRUE$"),
    k259P("m.b_run.head", "the action commits on acceptance", "^normalized$"),
    k259P("m.b_run.normalization_complete", "not finished", "^FALSE$"),
    k259P("m.b_run.ruv_normalized_obj", "publish never ran (first run)", "^null$"),
    k259P("m.b_run.pending", "nothing pending", "^none$"),
    k259P("m.b_run.notes", "the figures action's failure reaches the tab", "Error in normalization: .*finite"),
    k259P("m.b_apply.notes", "Apply is refused by diaCorrelationSource", "normalization must complete before correlation filtering"),
    k259P("m.b_rerun.normalization_complete", "Run again fails the same way", "^FALSE$"),
    k259P("m.b_reset.head", "Reset recovers", "^protein_replicate_filtered$"),
    k259P("m.b_group_rerun.normalization_complete", "'Color by' = group then Run completes", "^TRUE$"),
    # (c) 16 groups x 2, default 'Shape by'
    k259P("c.direct.status", "the post-normalisation PCA fails at ggsave", "Insufficient values in manual scale"),
    k259P("c.c_entry.events", "the pre-normalisation PCA fails at tab entry", "save pre_norm_pca.png = error: .*Insufficient values"),
    k259P("c.c_run.events", "the run's request has redraw_pre FALSE", "redraw_pre=FALSE"),
    k259P("c.c_run.head", "committed", "^normalized$"),
    k259P("c.c_run.normalization_complete", "not finished", "^FALSE$"),
    k259P("c.c_run.pending", "publication pending", "^normalization:FALSE$"),
    k259P("c.c_run.notes", "the error notice names the shape limit", "Error in normalization: .*Insufficient values"),
    k259P("c.c_run.warn", "default warn", "^0$"))

# ---------------------------------------------------------------------------------------------------------------------
k259Log <- function(...) cat(sprintf("[KIT-259 G0-B %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)), file = stderr())
k259Text <- function(value, n = 3000L) {
    if (is.null(value)) return("<null>")
    value <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(value)), collapse = " | "))
    if (!nzchar(value)) "<none>" else substr(value, 1L, n)
}
k259Rows <- new.env(parent = emptyenv())
k259Rows$rows <- list()
k259Put <- function(key, value) {
    value <- k259Text(value)
    k259Rows$rows[[length(k259Rows$rows) + 1L]] <- data.frame(key = key, value = value, stringsAsFactors = FALSE)
    k259Log(key, " = ", substr(value, 1L, 400L))
}
k259Write <- function() {
    observed <- do.call(rbind, k259Rows$rows)
    if (is.null(observed)) observed <- data.frame(key = character(), value = character())
    utils::write.table(observed, file.path(k259Out, "g0b-observed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    compare <- k259Predictions
    compare$observed <- vapply(compare$key, function(key) {
        hit <- observed$value[observed$key == key]
        if (length(hit)) hit[[length(hit)]] else "<not observed>"
    }, "")
    compare$held <- ifelse(compare$observed == "<not observed>", NA,
        mapply(function(p, o) grepl(p, o, perl = TRUE), compare$pattern, compare$observed))
    utils::write.table(compare, file.path(k259Out, "g0b-base.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    k259Log("written: ", sum(compare$held %in% TRUE), " held, ", sum(compare$held %in% FALSE), " not held, ",
        sum(is.na(compare$held)), " not observed")
}
k259Chain <- function(cnd) {
    out <- character()
    depth <- 0L
    while (!is.null(cnd) && inherits(cnd, "condition") && depth < 10L) {
        out <- c(out, sprintf("[%d] <%s> %s", depth, paste(class(cnd), collapse = "/"), conditionMessage(cnd)))
        cnd <- cnd$parent
        depth <- depth + 1L
    }
    paste(out, collapse = " ==> ")
}
k259Capture <- function(expr) {
    result <- tryCatch({
        force(expr)
        list(status = "ok")
    }, error = function(e) list(status = "error", error = e))
    list(status = result$status, chain = if (identical(result$status, "error")) k259Chain(result$error) else "<none>")
}

k259Log("loading the package and the helpers")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
k259Put("env.R", R.version.string)
k259Put("env.ggplot2", as.character(utils::packageVersion("ggplot2")))
k259Put("env.warn_outside", getOption("warn"))
k259Put("env.owner_capable", workflowOwnerCapable())

# ---------------------------------------------------------------------------------------------------------------------
# Spies: every figure built from a saved reference, every figure saved, every normalisation request (redraw_pre).
# ---------------------------------------------------------------------------------------------------------------------
k259Spy <- new.env(parent = emptyenv())
k259Spy$events <- character()
k259Event <- function(label, error = NULL) {
    k259Spy$events <- c(k259Spy$events, paste0(label, " = ", if (is.null(error)) "ok" else
        paste("error:", k259Text(conditionMessage(error), 300L))))
}
k259InstallSpies <- function(env = parent.frame()) {
    ns <- asNamespace("MultiScholaR")
    o_resolve <- get("resolveDIAReferenceQcPlot", ns)
    o_save <- get("saveProtNormQcPlotArtifact", ns)
    o_request <- get("diaNormalizationRequest", ns)
    testthat::local_mocked_bindings(
        resolveDIAReferenceQcPlot = function(reference, aesthetics, plot_type, ...) {
            got <- tryCatch(list(value = o_resolve(reference, aesthetics, plot_type, ...)), error = function(e) list(error = e))
            k259Event(sprintf("plot %s [%s/%s]", plot_type, aesthetics$color_var, aesthetics$shape_var), got$error)
            if (!is.null(got$error)) stop(got$error)
            got$value
        },
        saveProtNormQcPlotArtifact = function(qcDir, filename, plotObject, width, height, dpi = 150,
                                              savePlotFn = ggplot2::ggsave) {
            got <- tryCatch(list(value = o_save(qcDir, filename, plotObject, width, height, dpi, savePlotFn)),
                error = function(e) list(error = e))
            k259Event(paste("save", filename), got$error)
            if (!is.null(got$error)) stop(got$error)
            got$value
        },
        diaNormalizationRequest = function(manager, settings, norm_data, grouping) {
            value <- o_request(manager, settings, norm_data, grouping)
            k259Event(paste0("request redraw_pre=", isTRUE(value$redraw_pre)))
            value
        },
        .package = "MultiScholaR", .env = env)
}

#' The tab's flags and the head after a step.
k259Flags <- function(prefix, norm, manager) {
    k259Put(paste0(prefix, ".warn"), getOption("warn"))
    got <- tryCatch(shiny::isolate(list(
        complete = isTRUE(norm$normalization_complete),
        ruv = if (is.null(norm$ruv_normalized_obj)) "null" else "set",
        normalized = if (is.null(norm$normalized_protein_obj)) "null" else "set")), error = function(e) NULL)
    if (!is.null(got)) {
        k259Put(paste0(prefix, ".normalization_complete"), got$complete)
        k259Put(paste0(prefix, ".ruv_normalized_obj"), got$ruv)
        k259Put(paste0(prefix, ".normalized_protein_obj"), got$normalized)
    }
    k259Put(paste0(prefix, ".head"), tryCatch(manager$current_state, error = function(e) paste("<error>", conditionMessage(e))))
    k259Put(paste0(prefix, ".pending"), tryCatch(withDIAStateScope(manager, function() {
        reference <- getDIAStateReference(manager)
        request <- manager$.context$snapshot$action_outputs[[reference$entity_id]]
        if (is.null(request)) "none" else paste0(request$label, ":",
            isTRUE(manager$.context$snapshot$completed_outputs[[reference$entity_id]][[request$label]]))
    }), error = function(e) paste("<error>", conditionMessage(e))))
}

# ---------------------------------------------------------------------------------------------------------------------
# Part IP: the label-free lane.
# ---------------------------------------------------------------------------------------------------------------------
k259PartIp <- function() {
    k259InstallSpies()
    ns <- asNamespace("MultiScholaR")
    o_frame <- get("protNormFrameServer", ns)
    held <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(protNormFrameServer = function(input, output, session, context) {
        held$norm <- context$norm_data
        o_frame(input, output, session, context)
    }, .package = "MultiScholaR")
    step <- function(id, tab, fn) guiEntryStep(id, tab, function(api) {
        k259Spy$events <- character()
        got <- k259Capture({
            fn(api)
            api$session$flushReact()
        })
        if (!identical(got$status, "ok")) k259Put(paste0("ip.", id, ".step_error"), got$chain)
        k259Flags(paste0("ip.", id), held$norm, shiny::isolate(api$workflow_data$state_manager))
        k259Put(paste0("ip.", id, ".events"), k259Spy$events)
    })
    held$variable <- "batch"
    custom <- list(
        step("a", "normalization", function(api) {
            design <- shiny::isolate(api$workflow_data$design_matrix)
            color <- api$choices("normalization-color_variable")
            k259Put("ip.a.color_choices", paste(color, collapse = ","))
            k259Put("ip.a.shape_choices", paste(api$choices("normalization-shape_variable"), collapse = ","))
            k259Put("ip.a.sample_id_offered", "sample_id" %in% color)
            k259Put("ip.a.design_columns", paste(colnames(design), collapse = ","))
            for (column in intersect(c("group", "factor1", "factor2", "factor3", "batch"), colnames(design))) {
                k259Put(paste0("ip.a.design.", column), paste(as.character(design[[column]]), collapse = ","))
            }
            pairs <- vapply(color[color %in% colnames(design)], function(column) {
                values <- as.character(design[[column]])
                values <- values[!is.na(values) & nzchar(values)]
                sum(choose(table(values), 2))
            }, numeric(1))
            k259Put("ip.a.pairs_per_choice", paste(names(pairs), pairs, sep = "=", collapse = ","))
            zero <- names(pairs)[pairs == 0]
            held$variable <- if (length(zero)) zero[[length(zero)]] else "batch"
            k259Put("ip.a.zero_pair_variable", if (length(zero)) paste(zero, collapse = ",") else "<none>")
            k259Put("ip.a.warn", getOption("warn"))
            k259Put("ip.a.color_input", api$input("normalization-color_variable"))
            k259Put("ip.a.shape_input", api$input("normalization-shape_variable"))
        }),
        step("b_color", "normalization", function(api) api$set(stats::setNames(list(held$variable),
            "normalization-color_variable"))),
        step("b_run", "normalization", function(api) {
            api$set(list(`normalization-ruv_mode` = "skip"))
            api$set(list(`normalization-run_normalization` = 1L))
        }),
        step("b_apply", "normalization", function(api) {
            api$set(list(`normalization-norm_qc_tabs` = "Correlation Filtering",
                `normalization-min_pearson_correlation_threshold` = 0.5))
            api$set(list(`normalization-apply_correlation_filter` = 1L))
        }),
        step("b_skip", "normalization", function(api) api$set(list(`normalization-skip_correlation_filter` = 1L))),
        step("b_da", "da", function(api) {
            api$set(list(`differential_expression-load_filtered_session` = 2L))
            api$set(list(`differential_expression-formula_string` = "~ 0 + group",
                `differential_expression-da_q_val_thresh` = 0.05, `differential_expression-treat_lfc_cutoff` = 0))
            api$set(list(`differential_expression-run_da_analysis` = 1L))
            if (length(api$modals())) api$set(list(`differential_expression-acknowledge_qvalue_warning` = 1L))
        }),
        step("b_reset", "normalization", function(api) api$set(list(`normalization-reset_normalization` = 1L))),
        step("b_rerun", "normalization", function(api) api$set(list(`normalization-run_normalization` = 2L))),
        step("b_group_color", "normalization", function(api) api$set(list(`normalization-color_variable` = "group"))),
        step("b_group_rerun", "normalization", function(api) api$set(list(`normalization-run_normalization` = 3L))))
    lane_env <- environment(guiEntryRunLane)
    o_steps <- get("guiEntrySteps", envir = lane_env)
    cut_steps <- function(lane) {
        steps <- o_steps(lane)
        ids <- vapply(steps, `[[`, "", "id")
        if (!"norm_refusal" %in% ids) stop("the lane has no norm_refusal step: ", paste(ids, collapse = ", "))
        c(steps[seq_len(match("norm_refusal", ids))], custom)
    }
    locked <- bindingIsLocked("guiEntrySteps", lane_env)
    if (locked) unlockBinding("guiEntrySteps", lane_env)
    assign("guiEntrySteps", cut_steps, envir = lane_env)
    on.exit({
        assign("guiEntrySteps", o_steps, envir = lane_env)
        if (locked) lockBinding("guiEntrySteps", lane_env)
    }, add = TRUE)
    record <- guiEntryRunLane("prot_lfq")
    k259Put("ip.norm_captured", !is.null(held$norm))
    k259Put("ip.steps", paste(names(record$steps), collapse = ","))
    for (id in names(record$steps)) {
        k259Put(paste0("ip.", id, ".notes"), vapply(record$steps[[id]]$notifications, function(n)
            paste0("[", n$type, "] ", n$message), ""))
        k259Put(paste0("ip.", id, ".lane_head"), record$steps[[id]]$head)
    }
    invisible(record)
}

# ---------------------------------------------------------------------------------------------------------------------
# The DIA owner of KIT-122 (helper-prot-norm-core.R protNormCoreOwner, copied): optionally another sample count, other
# groups and an unset batch.
# ---------------------------------------------------------------------------------------------------------------------
k259Owner <- function(envir = parent.frame(), n_samples = 12L, groups = NULL, unset_batch = FALSE) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_samples <- n_samples
    original <- makeWorkflowCheckpointProtein(recipe)
    design <- original@design_matrix
    if (!is.null(groups)) {
        design$group <- groups
        design$tech_rep_group <- paste(groups, seq_along(groups), sep = "_")
    }
    if (unset_batch) design$batch <- NA_character_
    original@design_matrix <- design
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$design_matrix <- original@design_matrix
    workflow$config_list <- original@args
    initializeDIANativeImport(workflow, paths, session, "diann", "iq")
    manager <- workflow$state_manager
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$saveState("protein_replicate_filtered", original, list(), "source")
    list(root = root, paths = paths, workflow = workflow, manager = manager, original = original)
}

#' Drive testServer(mod_prot_norm_server) through `steps` (named function(session)), recording after each (a reduced
#' copy of test-prot-norm-core.R protNormW3Drive).
k259Drive <- function(prefix, fixture, steps, managed, selected_tab = NULL) {
    k259InstallSpies()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    probe <- new.env()
    if (managed) protNormW3Managed(fixture, probe)
    args <- list(workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core")
    if (!is.null(selected_tab)) args$selected_tab <- selected_tab
    utils::capture.output(shiny::testServer(mod_prot_norm_server, args = args, {
        norm <- session$getReturned()
        for (case in names(steps)) {
            seen <- length(probe$actions)
            log$notes <- character()
            k259Spy$events <- character()
            got <- k259Capture({
                steps[[case]](session)
                session$flushReact()
            })
            if (!identical(got$status, "ok")) k259Put(paste0(prefix, case, ".step_error"), got$chain)
            k259Flags(paste0(prefix, case), norm, fixture$manager)
            k259Put(paste0(prefix, case, ".notes"), log$notes)
            k259Put(paste0(prefix, case, ".events"), k259Spy$events)
            k259Put(paste0(prefix, case, ".actions"), vapply(probe$actions[seq_along(probe$actions) > seen],
                `[[`, "", "operation"))
        }
    }))
    invisible(NULL)
}

k259StepsB <- function() {
    click <- protNormW3Clicker()
    list(
        b_run = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip", color_variable = "batch"))),
        b_apply = click("apply_correlation_filter"),
        b_skip = click("skip_correlation_filter"),
        b_rerun = click("run_normalization"),
        b_reset = click("reset_normalization"),
        b_group_rerun = click("run_normalization", color_variable = "group"))
}

k259PartIpd <- function() {
    fixture <- k259Owner(unset_batch = TRUE)
    k259Drive("ipd.", fixture, k259StepsB(), managed = FALSE)
}

k259PartMLifecycle <- function() {
    fixture <- k259Owner()
    click <- protNormW3Clicker()
    k259Drive("m.lifecycle.", fixture, list(
        l_complete = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip"))),
        l_failed = click("run_normalization", ruv_mode = "manual", ruv_k = 50)), managed = TRUE)
}

k259PartM <- function() {
    if (!workflowOwnerCapable()) {
        k259Put("m.lifecycle.ok", "<owner not capable>")
        return(invisible(NULL))
    }
    k259PartMLifecycle()
    observed <- do.call(rbind, k259Rows$rows)
    value <- function(key) {
        hit <- observed$value[observed$key == key]
        if (length(hit)) hit[[length(hit)]] else ""
    }
    ok <- identical(value("m.lifecycle.l_complete.normalization_complete"), "TRUE") &&
        grepl("Normalization completed", value("m.lifecycle.l_complete.notes")) &&
        grepl("normalization-ruv", value("m.lifecycle.l_complete.actions")) &&
        grepl("[error 10] Error in normalization: ", value("m.lifecycle.l_failed.notes"), fixed = TRUE)
    k259Put("m.lifecycle.ok", ok)
    if (!ok) return(invisible(NULL))
    fixture <- k259Owner(unset_batch = TRUE)
    k259Drive("m.", fixture, k259StepsB(), managed = TRUE)
}

k259Groups16 <- paste0("g", sprintf("%02d", rep(seq_len(16L), each = 2L)))

k259PartCDirect <- function() {
    # Direct: the generator on the S4 object.
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_samples <- 32L
    object <- makeWorkflowCheckpointProtein(recipe)
    design <- object@design_matrix
    design$group <- k259Groups16
    design$tech_rep_group <- paste(k259Groups16, seq_along(k259Groups16), sep = "_")
    object@design_matrix <- design
    directory <- tempfile("k259-c-")
    dir.create(directory)
    k259InstallSpies()
    k259Spy$events <- character()
    got <- k259Capture(generateProtNormPostNormalizationQcArtifacts(object, directory,
        list(color_var = "group", shape_var = "group"), messageFn = function(...) invisible(NULL)))
    k259Put("c.direct.status", if (identical(got$status, "ok")) "ok" else paste("error:", got$chain))
    k259Put("c.direct.events", k259Spy$events)
    files <- list.files(directory)
    k259Put("c.direct.files", if (length(files)) paste0(files, " (", file.size(file.path(directory, files)), " bytes)") else "<none>")
    k259Put("c.direct.warn", getOption("warn"))
}

k259PartC <- function() {
    k259PartCDirect()
    # The tab, in process: tab entry, then Run Normalization with the default 'Shape by'.
    fixture <- k259Owner(n_samples = 32L, groups = k259Groups16)
    tab <- shiny::reactiveVal("setup")
    click <- protNormW3Clicker()
    k259Drive("c.", fixture, list(
        c_entry = function(session) tab("normalization"),
        c_run = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))), managed = FALSE, selected_tab = tab)
    k259Put("c.c_run.qc_files", list.files(fixture$paths$protein_qc_dir))
}

for (part in c("k259PartIp", "k259PartIpd", "k259PartM", "k259PartC")) {
    k259Log("part ", part)
    got <- k259Capture(get(part)())
    if (!identical(got$status, "ok")) k259Put(paste0("part.", part, ".error"), got$chain)
    k259Write()
}
k259Log("done")
