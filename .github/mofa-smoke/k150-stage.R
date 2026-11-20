# KIT-150 (branch run-mofa-measure-k150): one measured stage of the integration module, in a fresh R process.
# Declaration: docs/qa/workflow-execution/kit/kit-150/declaration.md.
#
# Usage: Rscript --vanilla k150-stage.R <stage arguments (JSON file)>
# Started by k150-run-stage.sh through the bundle's launcher (RETICULATE_PYTHON = the bundle's Python), inside its own
# memory cgroup. The arguments name the stage, the build's library, the study and the folders.
#
# Writes <out>.baseline.json as soon as the baseline is taken, then <out>.json when the stage ends (completed or
# error). A stage stopped by the memory or time limit leaves only the baseline file; the driver records the stop.
#
# Figures (declaration, "Terms"):
#   baseline     after MultiScholaR and MOFA2 are loaded and mofapy2 is imported through reticulate (and, for the input
#                write, after the synthetic study is generated): gc(full = TRUE), then gc(reset = TRUE); R heap in use
#                (Ncells + Vcells, MiB) and VmRSS; VmHWM is then reset to the current RSS (/proc/self/clear_refs, 5).
#   heap_peak    gc() "max used" (Ncells + Vcells, MiB) at stage end minus the baseline's heap in use.
#   vmhwm        VmHWM at stage end (MiB): the process peak RSS since the baseline, the embedded Python included.
#   python_peak  tracemalloc's peak, started at the trainer's entry (before MOFA2's create_mofa and the hand-off) and
#                read when training returns.
#   wall         elapsed seconds from the baseline to the stage's end (package loading excluded; load_s apart).

k150_started <- proc.time()[["elapsed"]]
k150_args <- jsonlite::read_json(commandArgs(trailingOnly = TRUE)[[1L]], simplifyVector = TRUE)
k150_rec <- new.env(parent = emptyenv())
k150_rec$marks <- list()

k150_status_mib <- function(field) {
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    if (length(line) != 1L) return(NA_real_)
    round(as.numeric(sub("^[^0-9]*([0-9]+).*$", "\\1", line)) / 1024, 1)
}

k150_write <- function(path, value) {
    tmp <- paste0(path, ".tmp")
    jsonlite::write_json(value, tmp, auto_unbox = TRUE, digits = NA, null = "null", na = "null", pretty = TRUE)
    file.rename(tmp, path)
}

k150_common <- function() {
    list(id = k150_args$id, stage = k150_args$stage, build = k150_args$build, commit = k150_args$commit,
        size = k150_args$size, rep = k150_args$rep, case = k150_args$case)
}

# --- Load: the build's library first; MultiScholaR, MOFA2, mofapy2 ---------------------------------------------------
.libPaths(c(k150_args$lib, .libPaths()))
suppressPackageStartupMessages({
    loadNamespace("MultiScholaR")
    loadNamespace("MOFA2")
    loadNamespace("reticulate")
})
k150_ns <- asNamespace("MultiScholaR")
k150_loaded_from <- normalizePath(find.package("MultiScholaR"), mustWork = TRUE)
if (!identical(dirname(k150_loaded_from), normalizePath(k150_args$lib, mustWork = TRUE))) {
    stop(sprintf("MultiScholaR loaded from %s, not from the build's library %s", k150_loaded_from, k150_args$lib))
}
k150_marker <- file.path(k150_args$lib, "k150-build-commit.txt")
k150_lib_commit <- if (file.exists(k150_marker)) trimws(readLines(k150_marker, n = 1L)) else NA_character_
if (!identical(k150_lib_commit, k150_args$commit)) {
    stop(sprintf("the build's library holds commit %s, the stage names %s", k150_lib_commit, k150_args$commit))
}
invisible(reticulate::import("mofapy2"))
invisible(reticulate::import("numpy"))
k150_tracemalloc <- reticulate::import("tracemalloc")
# KIT-148 item 2's internal switch (builds B3 and C only): off, as INT-22 removes the direct hand-off. Set in every
# build; the builds without the switch do not read it.
options(MultiScholaR.integration_mofa_handoff = "run_mofa")
k150_h <- new.env(parent = k150_ns)

# --- The synthetic study (declaration, "Sizes and the synthetic studies") -------------------------------------------
K150_WORKFLOWS <- c(proteomics = "LFQ", transcriptomics = "TRANSCRIPTOMICS", lipidomics = "LIPIDOMICS",
    metabolomics = "METABOLOMICS")
K150_SCALES <- c(proteomics = "log2", transcriptomics = "log2_cpm", lipidomics = "log2", metabolomics = "log2")
K150_VIEWS <- c(proteomics = "proteomics", transcriptomics = "transcriptomics", lipidomics = "lipidomics_LCMS_Pos",
    metabolomics = "metabolomics_LCMS_Pos")
K150_PREFIX <- c(proteomics = "PROT", transcriptomics = "TRAN", lipidomics = "LIPI", metabolomics = "META")
# Which of the five latent factors each omic carries (factor 1, the group factor, in all).
K150_ACTIVE <- list(proteomics = c(1L, 2L, 4L), transcriptomics = c(1L, 2L), lipidomics = c(1L, 3L, 5L),
    metabolomics = c(1L, 3L))

k150_omics <- function() as.character(k150_args$omics)

k150_study <- function() {
    n <- as.integer(k150_args$n)
    p <- as.integer(k150_args$p)
    omics <- k150_omics()
    companion <- k150_args$companion_p
    set.seed(as.integer(k150_args$seed))
    groups <- rep(c("A", "B"), length.out = n)
    scores <- matrix(stats::rnorm(n * 5L), n, 5L)
    scores[, 1L] <- scores[, 1L] + ifelse(groups == "B", 1, -1)
    out <- list()
    for (i in seq_along(omics)) {
        omic <- omics[[i]]
        rows <- if (!is.null(companion) && !is.na(companion) && i > 1L) as.integer(companion) else p
        set.seed(as.integer(k150_args$seed) + i)
        loadings <- matrix(0, rows, 5L)
        for (k in K150_ACTIVE[[omic]]) {
            hit <- sample.int(rows, ceiling(0.2 * rows))
            loadings[hit, k] <- stats::rnorm(length(hit))
        }
        spread <- stats::runif(rows, 0.3, 1.5)
        centre <- stats::rnorm(rows, 20, 2)
        values <- loadings %*% t(scores) + matrix(stats::rnorm(rows * n), rows, n) * spread + centre
        loadings <- NULL
        values[sample.int(length(values), floor(0.03 * length(values)))] <- NA_real_
        runs <- sprintf("%s_R%04d", K150_PREFIX[[omic]], seq_len(n))
        colnames(values) <- runs
        ids <- sprintf("%s_F%06d", K150_PREFIX[[omic]], seq_len(rows))
        view <- K150_VIEWS[[omic]]
        matrix_table <- data.frame(feature_id = ids, values, check.names = FALSE, stringsAsFactors = FALSE)
        values <- NULL
        samples <- data.frame(Run = runs, group = groups, biological_sample = sprintf("BS%04d", seq_len(n)),
            stringsAsFactors = FALSE)
        features <- data.frame(feature_id = ids, display_name = sprintf("%s feature %d", omic, seq_len(rows)),
            stringsAsFactors = FALSE)
        manifest <- list(omic_type = omic, workflow = K150_WORKFLOWS[[omic]], matrix_role = "normalised",
            scale = K150_SCALES[[omic]], correction = NULL,
            source_state = list(name = "synthetic", synthetic = TRUE, note = "KIT-150 synthetic study"),
            sample_id_column = "Run", group_column = "group", technical_replicate_column = "biological_sample",
            biological_sample_column = "biological_sample")
        out[[omic]] <- list(manifest = manifest, samples = samples, matrices = stats::setNames(list(matrix_table), view),
            features = stats::setNames(list(features), view))
    }
    out
}

# --- The model's parameters, every one passed (declaration, "Model settings") ---------------------------------------
k150_parameters <- function(case) {
    capped <- identical(case, "capped")
    parameters <- list(technical_replicates = "average", internal_standards = "drop", missing_samples = "keep_missing",
        feature_selection = if (capped) "top_variance" else "none",
        features_per_view = if (capped) as.integer(k150_args$cap) else NA_integer_,
        scaling = "none", scale_views = TRUE, center_groups = TRUE, num_factors = 5L, convergence_mode = "slow",
        maxiter = if (capped) 10000L else as.integer(k150_args$floor_rounds), seed = 20261119L, precision = "double")
    # KIT-149's variant (build C only; the other builds have one variant and ignore the name).
    if (nzchar(k150_args$variant %||% "")) parameters$variant <- k150_args$variant
    parameters
}

`%||%` <- function(a, b) if (is.null(a)) b else a

k150_mark <- function(name) {
    k150_rec$marks[[length(k150_rec$marks) + 1L]] <- list(name = name,
        at_s = round(proc.time()[["elapsed"]] - k150_rec$t0, 3), vmhwm_mib = k150_status_mib("VmHWM"),
        vmrss_mib = k150_status_mib("VmRSS"))
}

k150_selected <- function(views) {
    dims <- vapply(views, function(m) c(nrow(m), ncol(m)), numeric(2))
    k150_rec$selected <- data.frame(view = names(views), rows = dims[1L, ], samples = dims[2L, ], stringsAsFactors = FALSE,
        row.names = NULL)
}

# --- Store owners ----------------------------------------------------------------------------------------------------
k150_owner_new <- function(root) newIntegrationStateOwner(root, "integration")
k150_owner_reopen <- function(root) {
    identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
    manager <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id),
        WorkflowState$new(audit_enabled = FALSE), reopen = TRUE)
    manager$.context$reference_states <- TRUE
    manager$.context$omic_type <- "integration"
    manager
}
k150_dirs <- function() file.path(k150_args$inputs, k150_omics())
k150_views <- function() unname(K150_VIEWS[k150_omics()])

# --- The trainers handed to runIntegrationModel (its `trainer` argument; the kernel calls trainer(views, options)) ---
k150_measuring_trainer <- function(views, options) {
    k150_selected(views)
    k150_mark("handoff_start")
    k150_tracemalloc$start()
    fit <- integrationTrainMofa(views, options)
    traced <- k150_tracemalloc$get_traced_memory()
    k150_tracemalloc$stop()
    k150_rec$python_peak_mib <- round(as.numeric(traced[[2L]]) / 2^20, 1)
    k150_mark("training_end")
    fit
}
k150_stop_trainer <- function(views, options) {
    k150_selected(views)
    k150_mark("handoff_reached")
    k150_rec$handoff_reached <- TRUE
    stop(structure(class = c("k150_handoff_reached", "error", "condition"),
        list(message = "k150: the alignment is complete; stopped at the hand-off", call = NULL)))
}

# --- The stages --------------------------------------------------------------------------------------------------------
k150_prepare <- function(stage) {
    if (identical(stage, "write")) return(k150_study())
    NULL
}

k150_run <- function(stage, prepared) {
    switch(stage,
        write = {
            for (omic in names(prepared)) {
                entry <- prepared[[omic]]
                writeIntegrationInput(file.path(k150_args$inputs, omic), entry$manifest, entry$samples, entry$matrices,
                    entry$features)
            }
            "written"
        },
        import = {
            manager <- k150_owner_new(k150_args$store)
            importIntegrationInputs(manager, k150_dirs())
            closeWorkflowStateOwner(manager)
            "imported"
        },
        align = {
            # The model's own path (runIntegrationModel, the build's kernel or variant) up to the trainer's call; the
            # stand-in trainer stops there. Any other error is the stage's.
            manager <- k150_owner_reopen(k150_args$store)
            k150_rec$handoff_reached <- FALSE
            tryCatch(runIntegrationModel(manager, k150_parameters("capped"), trainer = k150_stop_trainer),
                error = function(error) if (!isTRUE(k150_rec$handoff_reached)) stop(error))
            closeWorkflowStateOwner(manager)
            if (!isTRUE(k150_rec$handoff_reached)) stop("the alignment stage did not reach the hand-off")
            "aligned"
        },
        train = ,
        floor = {
            manager <- k150_owner_reopen(k150_args$store)
            run <- runIntegrationModel(manager, k150_parameters(k150_args$case), trainer = k150_measuring_trainer)
            k150_mark("product_written")
            closeWorkflowStateOwner(manager)
            k150_rec$model <- list(status = run$status, iterations = run$value$elbo$iterations,
                elbo_last = run$value$elbo$last,
                factors_sha256 = digest::digest(round(as.matrix(run$value$factors[-1L]), 8), algo = "sha256"),
                features_model = run$value$view_summary$features_model, views = run$value$view_summary$view,
                samples = nrow(run$value$samples))
            run <- NULL
            "trained"
        },
        results = {
            manager <- k150_owner_reopen(k150_args$store)
            export <- manager$.context$snapshot$exports$integration_mofa_model
            if (is.null(export)) stop("no model is accepted on the store's head")
            value <- withWorkflowStateScope(manager, function()
                hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value)
            k150_mark("product_read")
            # As the Results tab calls it (R/mod_integration_results.R of each build).
            results <- if ("manager" %in% names(formals(integrationResults))) {
                integrationResults(value, integrationImportedObject(manager), list(), integrationGuiTexts(), manager = manager)
            } else {
                integrationResults(value, integrationImportedObject(manager), list(), integrationGuiTexts())
            }
            closeWorkflowStateOwner(manager)
            k150_rec$results_parts <- names(results)
            "read"
        },
        app_import = {
            # The Import tab's action (R/mod_integration_import.R): a new store, integrationCommitImport on the request
            # integrationImportRequest gives for every data set selected (the folders in the scan's order).
            manager <- k150_owner_new(k150_args$store)
            request <- list(dirs = normalizePath(k150_dirs(), winslash = "/", mustWork = TRUE), omics = k150_omics(),
                views = k150_views())
            summary <- integrationCommitImport(manager, request)
            closeWorkflowStateOwner(manager)
            k150_rec$app_views <- summary$views
            "imported"
        },
        app_preview = {
            # The Sample Matching tab's preview (R/mod_integration_alignment.R of each build).
            manager <- k150_owner_reopen(k150_args$store)
            preview <- if (exists("INTEGRATION_PREVIEW_FIELDS", envir = k150_ns, inherits = FALSE)) {
                integrationAlignmentPreview(integrationImportedObject(manager, INTEGRATION_PREVIEW_FIELDS), k150_views(), NULL,
                    manager = manager)
            } else {
                integrationAlignmentPreview(integrationImportedObject(manager), k150_views(), NULL)
            }
            closeWorkflowStateOwner(manager)
            k150_rec$preview_samples <- nrow(preview$aligned$samples)
            "previewed"
        },
        fail_memory = {
            # Induced failure (lifecycle): allocate 64 MiB blocks (zero-filled, so resident) until the cgroup's limit
            # stops the process.
            held <- list()
            repeat held[[length(held) + 1L]] <- numeric(8388608)
        },
        fail_time = {
            Sys.sleep(3600)
            "slept"
        },
        fail_error = {
            # Induced failure (lifecycle): the import refuses folders that do not exist.
            manager <- k150_owner_new(k150_args$store)
            importIntegrationInputs(manager, c(file.path(k150_args$store, "absent-a"), file.path(k150_args$store, "absent-b")))
            "imported"
        },
        stop(sprintf("unknown stage %s", stage)))
}

# Post-measures, taken after the stage's figures are recorded (they do not enter the figures).
k150_post <- function(stage) {
    if (identical(stage, "import") && exists("openIntegrationInputView", envir = k150_ns, inherits = FALSE)) {
        # Claim 3's allowance: the feature and sample tables as the import holds them (contract v2 builds).
        tables <- 0
        runs <- integer()
        for (dir in k150_dirs()) {
            input <- readIntegrationInput(dir)
            tables <- tables + as.numeric(utils::object.size(input$samples)) +
                sum(vapply(input$views, function(view) as.numeric(utils::object.size(view$features)), numeric(1)))
            runs <- c(runs, nrow(input$samples))
        }
        k150_rec$tables_mib <- round(tables / 2^20, 1)
        k150_rec$row_group_mib <- round(10000 * max(runs) * 8 / 2^20, 1)
    }
    if (identical(stage, "write")) {
        files <- list.files(k150_args$inputs, recursive = TRUE, full.names = TRUE)
        k150_rec$inputs_mib <- round(sum(file.size(files)) / 2^20, 1)
    }
}

for (name in c("k150_study", "k150_owner_new", "k150_owner_reopen", "k150_measuring_trainer", "k150_stop_trainer",
    "k150_prepare", "k150_run", "k150_post", "k150_dirs", "k150_views", "k150_parameters")) {
    fn <- get(name)
    environment(fn) <- k150_h
    assign(name, fn)
}

# --- Measure ---------------------------------------------------------------------------------------------------------
k150_stage <- k150_args$stage
prepared <- k150_prepare(k150_stage)
invisible(gc(full = TRUE))
k150_g0 <- gc(full = TRUE, reset = TRUE)
writeLines("5", "/proc/self/clear_refs")
k150_rec$t0 <- proc.time()[["elapsed"]]
k150_baseline <- list(heap_mib = round(sum(k150_g0[, 2L]), 1), vmrss_mib = k150_status_mib("VmRSS"),
    vmhwm_after_reset_mib = k150_status_mib("VmHWM"), load_s = round(k150_rec$t0 - k150_started, 2),
    threads = list(OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", NA), OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", NA)),
    loaded_from = k150_loaded_from)
k150_write(paste0(k150_args$out, ".baseline.json"), c(k150_common(), list(status = "baseline", baseline = k150_baseline)))

k150_outcome <- tryCatch(list(status = "completed", value = k150_run(k150_stage, prepared)),
    error = function(error) list(status = "error", class = class(error)[[1L]],
        message = substr(conditionMessage(error), 1L, 600L)))
k150_wall <- round(proc.time()[["elapsed"]] - k150_rec$t0, 2)
k150_g1 <- gc()
k150_figures <- list(heap_peak_mib = round(sum(k150_g1[, 6L]) - k150_baseline$heap_mib, 1),
    vmhwm_mib = k150_status_mib("VmHWM"), vmrss_end_mib = k150_status_mib("VmRSS"), wall_s = k150_wall,
    python_peak_mib = k150_rec$python_peak_mib %||% NA_real_)
prepared <- NULL
if (identical(k150_outcome$status, "completed")) {
    k150_post_error <- tryCatch({
        k150_post(k150_stage)
        NA_character_
    }, error = function(error) substr(conditionMessage(error), 1L, 300L))
} else {
    k150_post_error <- NA_character_
}
k150_record <- c(k150_common(), list(status = k150_outcome$status, error_class = k150_outcome$class %||% NA_character_,
    error_message = k150_outcome$message %||% NA_character_, baseline = k150_baseline, figures = k150_figures,
    marks = k150_rec$marks, selected = k150_rec$selected, model = k150_rec$model,
    tables_mib = k150_rec$tables_mib %||% NA_real_, row_group_mib = k150_rec$row_group_mib %||% NA_real_,
    inputs_mib = k150_rec$inputs_mib %||% NA_real_, app_views = k150_rec$app_views,
    preview_samples = k150_rec$preview_samples %||% NA_integer_, results_parts = k150_rec$results_parts,
    post_error = k150_post_error))
k150_write(paste0(k150_args$out, ".json"), k150_record)
cat(sprintf("k150 %s %s %s rep %s: %s; heap peak %.1f MiB, VmHWM %.1f MiB, %.1f s\n", k150_stage, k150_args$build,
    k150_args$size, k150_args$rep, k150_outcome$status, k150_figures$heap_peak_mib, k150_figures$vmhwm_mib, k150_wall))
if (!identical(k150_outcome$status, "completed")) cat("k150 error:", k150_outcome$message, "\n")
quit(save = "no", status = if (identical(k150_outcome$status, "completed")) 0L else 3L)
