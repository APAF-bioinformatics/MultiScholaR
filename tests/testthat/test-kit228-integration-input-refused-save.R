# KIT-228 (VAL-201): a Summary Save that does not write a new integration input removes the omic's earlier one (at the
# start of the Save, in one place for all four omics) and tells the user; the integration then cannot import an input the
# last Save did not write. Written before the fix: on the tests-only commit each E row finds the earlier folder still
# there (and, as the head's notice is not yet written, no notice), so it fails at the values the ticket states (section
# 5); the G rows and R01's second half pass on both commits. Lines `KIT-228 <row>` are the evidence the ticket's records
# copy (outcome.md: folder presence, file digests with the manifest's `written_at` masked, the notes).

library(testthat)

# ------------------------------------------------------------------------------------------------
# Copied unchanged from tests/testthat/test-kit207-summary-head-guard.R (itself copied from test-omics-norm-state-guards.R
# at 3d4a585d): lines 16-127, 129-175, 181-193, 238-256, 275, 279-286 and 352-371. That file is not edited.
# ------------------------------------------------------------------------------------------------

k194Omic <- function(class) c(lipid = "lipidomics", metab = "metabolomics")[[class]]

#' Copy of `kit172PreNorm` (test-kit172-downstream.R:127-137; A11): an owner at the QC-complete state (intensity filter,
#' QC finalise), as the normalisation tab finds it.
k194PreNorm <- function(root, values, class) {
    fixture <- kit172Owner(root, kit172Object(values, class), class)
    binding <- kit172Binding(class)
    kit172Quiet({
        runWorkflowOmicsReferenceIntensityFilter(binding, fixture$owner, 1, 0.5, fixture$workflow$config_list)
        commitWorkflowOmicsStoreStep(fixture$owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"),
            fixture$workflow$config_list, "QC complete")
    })
    fixture
}

#' Copy of `kit172InProcessRunner` (test-kit172-downstream.R:140-160; A11): the managed actions run in this process, each
#' commit through its declaration's wrapper, its result delivered as the worker's; a product delivers an empty display.
k194InProcessRunner <- function(manager) {
    list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
        run = function(actions, done = function() NULL, failed = function(error) stop(error)) {
            outcome <- tryCatch({
                for (action in actions) {
                    if (is.function(action$before)) action$before()
                    if (is.null(action$context)) {
                        action$deliver(list(display = list()))
                    } else {
                        step <- list(operation = action$operation, parameters = action$parameters, state_name = action$context$state_name,
                            description = action$context$description)
                        result <- runWorkflowOmicsReferenceStep(manager, step, action$context$config)
                        action$deliver(list(stateReference = result$reference, value = result$value))
                    }
                }
                NULL
            }, error = function(e) e)
            if (inherits(outcome, "error")) failed(outcome) else done()
        })
}

#' The recording runner of the managed cases (B04): `step$active()` is NULL while idle and an action while `busy`; a
#' press while `busy` or `hold` is queued, and `start()` runs what was queued; `actions` given as a function are resolved
#' when they start (`R/utils_workflow_omics_actions.R:43`), a refusal then reaching `failed`; the actions run as
#' `k194InProcessRunner` runs them, on the bus's owner of the moment. `operations` records every action run.
k194Runner <- function(workflow, busy = FALSE, hold = FALSE) {
    rec <- new.env(parent = emptyenv())
    rec$operations <- character()
    rec$queued <- list()
    rec$busy <- busy
    rec$hold <- hold
    execute <- function(actions, done, failed) {
        resolved <- tryCatch(if (is.function(actions)) actions() else actions, error = function(error) error)
        if (inherits(resolved, "error")) {
            if (is.function(failed)) failed(resolved) else stop(resolved)
            return(invisible(NULL))
        }
        rec$operations <- c(rec$operations, vapply(resolved, function(action) action$operation, ""))
        k194InProcessRunner(workflow$state_manager)$run(resolved, done = done %||% function() NULL,
            failed = failed %||% function(error) stop(error))
        invisible(NULL)
    }
    rec$start <- function() {
        queued <- rec$queued
        rec$queued <- list()
        rec$hold <- FALSE
        for (waiting in queued) execute(waiting$actions, waiting$done, waiting$failed)
        invisible(NULL)
    }
    runner <- new.env(parent = emptyenv())
    runner$step <- list(active = function() if (isTRUE(rec$busy)) "a running action")
    runner$progress <- function(...) list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL),
        close = function(...) invisible(NULL))
    runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
        if (isTRUE(rec$busy) || isTRUE(rec$hold)) {
            rec$queued[[length(rec$queued) + 1L]] <- list(actions = actions, done = done, failed = failed)
            return(invisible(NULL))
        }
        execute(actions, done, failed)
    }
    rec$runner <- runner
    rec
}

#' A lipid or metab project at QC complete on F172 (seed 172001), with the omic's GUI configuration (its QC figures
#' drawn by nothing) and the tab's log.
k194Fixture <- function(root, class) {
    fixture <- k194PreNorm(root, kit172Values(172001L)$values, class)
    fixture$class <- class
    fixture$omic <- k194Omic(class)
    fixture$binding <- kit172Binding(class)
    fixture$prefix <- fixture$binding$prefix
    config <- omicGuiConfig(fixture$omic)
    config$binding$binding$qc_plots <- function(...) invisible(NULL)
    fixture$config <- config
    log <- new.env(parent = emptyenv())
    log$lines <- character()
    fixture$log <- log
    fixture$addLog <- function(entry) log$lines <- c(log$lines, entry)
    fixture
}

k194Close <- function(fixture) try(closeDIAStateAdapter(fixture$workflow$state_manager), silent = TRUE)

#' Close the bus's owner of the moment (a load replaces it) when the calling test ends; `fixture` is forced here, so a
#' loop's fixtures are each closed.
k194Defer <- function(fixture, envir = parent.frame()) {
    force(fixture)
    withr::defer(k194Close(fixture), envir = envir)
}

k194Owner <- function(fixture) fixture$workflow$state_manager
k194History <- function(fixture) as.character(unlist(k194Owner(fixture)$getHistory()))
k194Head <- function(fixture) k194Owner(fixture)$current_state
k194Ref <- function(fixture, state = NULL) getWorkflowStateReference(k194Owner(fixture), state)

#' Shiny's notifications recorded as "[type] text"; removal, progress and selector updates quiet.
k194Notes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default") {
            seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        },
        removeNotification = function(...) invisible(NULL),
        withProgress = function(expr, ...) expr,
        incProgress = function(...) invisible(NULL),
        updateSelectInput = function(...) invisible(NULL),
        .package = "shiny", .env = env)
    seen
}

k194Noted <- function(seen, text) any(grepl(text, seen$notes, fixed = TRUE))

#' The tab's inputs: ITSD off unless asked (by pattern), offset 0, `method`, RUV skip unless `ruv` (the legacy list).
k194Input <- function(method = "quantile", itsd = FALSE, ruv = NULL, threshold = 0.5) {
    list(apply_itsd = itsd, itsd_aggregation = "median", log_offset = 0, norm_method = method,
        ruv_mode = if (is.null(ruv)) "skip" else ruv$ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = "group", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_k = ruv$manual_k %||% 1L, ruv_percentage = ruv$manual_percentage %||% 20,
        min_pearson_correlation_threshold = threshold)
}

#' "Run Normalization Pipeline": the frame's `runOmicNormPipeline`, in process (`runner` NULL) or managed.
k194Run <- function(fixture, nd, input = k194Input(), runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(runOmicNormPipeline(input, fixture$workflow, list(), nd, fixture$addLog,
        function(...) list(color_var = "group", shape_var = NULL), fixture$config)))
}

#' The run a case builds on completed (an instrument check, VAL-54): head the normalised or RUV state, flag set.
k194ExpectRun <- function(fixture, nd, info = "") {
    expect_true(k194Head(fixture) %in% paste0(fixture$prefix, c("_normalized", "_ruv_corrected")),
        info = paste("run precondition:", info, k194Head(fixture), paste(fixture$log$lines, collapse = " | ")))
    expect_true(isTRUE(shiny::isolate(nd$normalization_complete)), info = paste("run precondition:", info))
}

k194Apply <- function(fixture, nd, runner = NULL, threshold = 0.5) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    kit172Quiet(shiny::isolate(applyOmicNormCorrelation(k194Input(threshold = threshold), fixture$workflow, nd, fixture$addLog,
        fixture$config)))
}

#' Revert: the Intensity Filter sub-tab's real observer under `testServer`, pressed until the head is `until` (or
#' `times` presses).
k194Revert <- function(fixture, until = NULL, times = 10L) {
    kit172Quiet(utils::capture.output(shiny::testServer(mod_omic_qc_intensity_server, args = list(workflow_data = fixture$workflow,
        omic_type = fixture$omic, experiment_label = "kit194", config = fixture$config), {
        press <- 0L
        while (press < times && (is.null(until) || !identical(fixture$workflow$state_manager$current_state, until))) {
            press <- press + 1L
            session$setInputs(revert_filter = press)
        }
    })))
    k194Head(fixture)
}

#' trx's "Run Normalization Pipeline" (`runTrxNormalisation`, as test-trx-norm-core.R drives it) with the notifications
#' recorded in `seen$notes` and the tab's log in `nd`.
k194TrxRun <- function(project, nd, inputs, seen, runner = NULL) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    notify <- function(ui, ..., type = "default", duration = 5) {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    set_status <- function(tab, value) {
        status <- project$workflow$tab_status
        status[[tab]] <- value
        project$workflow$tab_status <- status
    }
    withr::local_seed(20261124L)
    kit172Quiet(utils::capture.output(shiny::isolate(runTrxNormalisation(list(managed = runTrxNormManaged, in_process = runTrxNormInProcess),
        transcriptomicsDeclarationBinding, project$workflow$state_manager, trxNormGridParameters(inputs), nd, project$workflow$config_list,
        project$paths, omicNormAddLog(nd), notify, set_status))))
    invisible(NULL)
}

k207Refused <- function(text) paste("The session could not be saved:", text)

#' A fresh project's directories for the omic (the Summary's `projectDirs`).
k207Dirs <- function(root, omic) {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    stats::setNames(list(paths), omic)
}

k207IntegrationDir <- function(paths) paths$integration_dir %||% file.path(paths$base_dir, "integration")

#' The reference DA on the head (as `kit172Journey`); it commits no state.
k207Da <- function(fx) {
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(fx$binding, k194Owner(fx), k194Ref(fx),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts()))
    invisible(fx)
}

#' A lipid or metab fixture after "Run Normalization Pipeline" (`input`), Apply (unless `apply` is FALSE) and the DA
#' (when `da`); closed when the calling test ends.
k207Normalised <- function(root, class, input = k194Input(), apply = TRUE, da = FALSE, env = parent.frame()) {
    fx <- k194Fixture(root, class)
    k194Defer(fx, envir = env)
    nd <- omicNormReactiveState()
    k194Run(fx, nd, input)
    k194ExpectRun(fx, nd, paste("KIT-207", class, basename(root)))
    if (apply) expect_true(isTRUE(k194Apply(fx, nd)), info = paste("apply precondition", class, basename(root)))
    if (da) k207Da(fx)
    list(fx = fx, nd = nd)
}


# ------------------------------------------------------------------------------------------------
# KIT-228's harness.
# ------------------------------------------------------------------------------------------------

# The notices of section 4 (DRAFTS, VAL-41: the owner's wording pass may change them; Re-check 1 R02 and R03 texts).
k228Texts <- list(
    failed = function(omic) sprintf(paste("This omic's earlier integration input (integration/%s/) was removed because this save",
        "did not write a new one (the reason is in the message above). The integration cannot use %s until this Summary is",
        "saved without errors."), omic, omic),
    failed_unnamed = function(omic) sprintf(paste("This omic's earlier integration input (integration/%s/) was removed because this save",
        "did not write a new one. The integration cannot use %s until this Summary is saved without errors."), omic, omic),
    gate = function(omic) sprintf(paste("This omic's earlier integration input (integration/%s/) was removed because this save",
        "did not write a new one: the RUVIII export is missing. Run the RUVIII export, then save again before integrating."), omic),
    no_table = function(omic) sprintf(paste("This omic's earlier integration input (integration/%s/) was removed, and this save",
        "did not write a new one. Save this Summary again before integrating; if this message returns, finish the analysis first",
        "(the integration input is written from the final normalised table)."), omic),
    retry = function(omic) sprintf(paste("This omic's earlier integration input (integration/%s/) was removed by this save. If this",
        "save only finishes an earlier save's publication, that does not write it again: save once more before integrating."), omic),
    not_removed = function(folder) sprintf(paste("The earlier integration input at %s could not be removed. The integration may",
        "still read it, although this save did not write it: delete that folder, or save again once it can be removed."), folder))

k228Folder <- function(dirs, omic) file.path(k207IntegrationDir(dirs[[omic]]), omic)
k228Manifest <- function(folder) file.path(folder, "integration_input.json")
k228ManifestSha <- function(folder) {
    if (file.exists(k228Manifest(folder))) digest::digest(file = k228Manifest(folder), algo = "sha256") else "<absent>"
}

#' Every file of an input folder by SHA-256, the manifest without its `written_at` line (the only line that differs between
#' two saves of one table): the G1 guard's evidence.
k228Digests <- function(folder) {
    files <- sort(list.files(folder))
    stats::setNames(vapply(files, function(name) {
        path <- file.path(folder, name)
        if (!identical(name, "integration_input.json")) return(digest::digest(file = path, algo = "sha256"))
        lines <- readLines(path, warn = FALSE)
        digest::digest(paste(lines[!grepl("\"written_at\"", lines, fixed = TRUE)], collapse = "\n"), algo = "sha256",
            serialize = FALSE)
    }, ""), files)
}

#' "Save Workflow Arguments" under the omic's group (as `k207Save`), in process (`runner` NULL) or with a runner double;
#' `args` replaces the study-parameter and record writers (the proteomics owner of test-workflow-integration-input.R).
#' Records the notifications ("[type] text"), `done` and `failed`.
k228Save <- function(workflow, omic, dirs, runner = NULL, args = list()) {
    rec <- new.env(parent = emptyenv())
    rec$notes <- character()
    rec$done <- 0L
    rec$failed <- list()
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    call <- c(list(omicGuiConfig(omic)$summary, workflow, dirs, omic, "k228", "",
        notify = function(message, type, duration) rec$notes <- c(rec$notes, paste0("[", type, "] ", message)),
        done = function(saved) rec$done <- rec$done + 1L,
        failed = function(error) rec$failed[[length(rec$failed) + 1L]] <- error,
        saveRDSFn = function(object, path) invisible(NULL), catFn = function(...) NULL), args)
    kit172Quiet(utils::capture.output(do.call(saveWorkflowSummaryParameters, call)))
    rec
}

#' The evidence line of a row (both commits): the folder's presence, its digests and the notes.
k228Line <- function(row, omic, rec, folder) {
    digests <- if (dir.exists(folder)) k228Digests(folder) else character()
    cat(sprintf("KIT-228 %s %s\tfolder=%s\tfiles=%s\tdone=%d\tfailed=%s\tnotes=%s\n", row, omic, dir.exists(folder),
        if (length(digests)) paste0(names(digests), ":", substr(digests, 1L, 12L), collapse = ",") else "<none>", rec$done,
        paste(vapply(rec$failed, function(e) class(e)[[1L]], ""), collapse = ","), paste(rec$notes, collapse = " | ")))
}

#' The precondition of every E row: a successful Save, which leaves the omic's folder with its manifest (an instrument
#' check, VAL-54). Returns the manifest's SHA-256 (M).
k228FirstSave <- function(workflow, omic, dirs, info, runner = NULL, args = list()) {
    folder <- k228Folder(dirs, omic)
    rec <- k228Save(workflow, omic, dirs, runner = runner, args = args)
    k228Line(paste(info, "first"), omic, rec, folder)
    expect_identical(rec$done, 1L, info = paste(info, "first save:", paste(rec$notes, collapse = " | ")))
    expect_true(file.exists(k228Manifest(folder)), info = paste(info, "first save wrote no input:", paste(rec$notes, collapse = " | ")))
    k228ManifestSha(folder)
}

#' What a Save that writes no new input must leave: the earlier folder gone, the notice once.
k228ExpectRemoved <- function(rec, folder, text, info) {
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_false(dir.exists(folder), info = info)
    expect_identical(sum(rec$notes == paste0("[warning] ", text)), 1L, info = info)
}

#' The runner doubles of the managed route: `failed` with a class-less error as `R/utils_workflow_execution.R` hands it on
#' (E5), or neither callback (E6, the application closed mid-write).
k228FailingRunner <- function() {
    list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
        run = function(actions, done = NULL, failed = NULL, dropped = NULL) {
            failed(simpleError("managed action ended: worker_failed"))
            invisible(NULL)
        })
}

k228SilentRunner <- function() {
    list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
        run = function(actions, done = NULL, failed = NULL, dropped = NULL) invisible(NULL))
}

#' The managed route's "publication pending" sequence (I05, C02): the first request's worker writes the input, then the step
#' fails as publication pending (`failed` is called, `done` and `failed` are kept); a second Save finds the runner retrying
#' and resumes the first request: its `done` runs (`resume = "done"`) or, the request superseded, its `failed`
#' (`resume = "dropped"`); nothing is recomputed or rewritten (`R/utils_workflow_omics_actions.R:73-82`).
k228RetryRunner <- function(write_input, resume = c("done", "dropped")) {
    resume <- match.arg(resume)
    state <- new.env(parent = emptyenv())
    state$pending <- FALSE
    list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
        retrying = function() isTRUE(state$pending),
        run = function(actions, done = NULL, failed = NULL, dropped = NULL) {
            if (isTRUE(state$pending)) {
                state$pending <- FALSE
                if (identical(resume, "done")) state$done() else
                    state$failed(simpleError("the earlier publication belongs to a superseded state and was dropped"))
                return(invisible(NULL))
            }
            write_input()
            state$pending <- TRUE
            state$done <- done
            state$failed <- failed
            failed(simpleError("publication pending: press the button again to write the files"))
            invisible(NULL)
        })
}

#' A runner whose omic has ANOTHER step's publication pending (C01: `retrying()` is true for the whole omic): the Save's own
#' export fails as a managed failure does.
k228BusyRunner <- function() {
    list(progress = function(...) list(inc = function(...) invisible(NULL), close = function() invisible(NULL)),
        retrying = function() TRUE,
        run = function(actions, done = NULL, failed = NULL, dropped = NULL) {
            failed(simpleError("managed action ended: worker_failed"))
            invisible(NULL)
        })
}

#' Put an input folder back from a copy taken after a successful Save (a refusal row needs one more earlier folder).
k228Restore <- function(stash, folder) {
    unlink(folder, recursive = TRUE)
    file.copy(file.path(stash, basename(folder)), dirname(folder), recursive = TRUE)
}

#' The managed failure's expectations, base and head alike: the class is lost, so the Summary names no refusal.
k228ExpectClassLost <- function(rec, info) {
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_true(length(rec$failed) == 1L && !inherits(rec$failed[[1L]], "integration_input_invalid"), info = info)
    expect_true("[warning] Study parameters saved with warnings" %in% rec$notes, info = info)
    expect_false(any(grepl("could not be saved", rec$notes, fixed = TRUE)), info = info)
}

#' G1 and the notice guard of a successful Save after a successful one: the folder is there with a manifest, `done` once,
#' and no notice about a removal; the evidence line carries the digests (compared base against head from the logs).
k228ExpectSavedAgain <- function(workflow, omic, dirs, info, runner = NULL, args = list()) {
    folder <- k228Folder(dirs, omic)
    rec <- k228Save(workflow, omic, dirs, runner = runner, args = args)
    k228Line(paste(info, "again"), omic, rec, folder)
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_identical(rec$done, 1L, info = info)
    expect_true(file.exists(k228Manifest(folder)), info = info)
    expect_false(any(grepl("earlier integration input", rec$notes, fixed = TRUE)), info = info)
    expect_error(readIntegrationInput(folder), NA, info = info)
}

#' G4: a refused Save with no earlier folder says nothing about a removal and creates no folder of the omic.
k228ExpectNoNotice <- function(rec, folder, info) {
    info <- paste(info, "| notes:", paste(rec$notes, collapse = " | "))
    expect_false(any(grepl("earlier integration input", rec$notes, fixed = TRUE)), info = info)
    expect_false(dir.exists(folder), info = info)
}

#' The proteomics Summary owner of test-workflow-integration-input.R (a classic route without RUV, copied from
#' `integrationProteinOwner` and the KIT-061 case that uses it), with the writers stubbed as `integrationSummarySave` does.
k228ProtOwner <- function(root) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit228", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$setWorkflowType("LFQ")
    manager$.context$reference_states <- TRUE
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("normalized", protein, list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_applied = FALSE), "skipped")
    manager$saveState("correlation_filtered", protein, list(min_pearson_correlation_threshold = 0.5), "correlation")
    list(manager = manager, workflow = list2env(list(state_manager = manager, contrasts_tbl = NULL, config_list = NULL)))
}

k228ProtWriters <- list(createWorkflowArgsFn = function(...) "study_parameters.txt", writeAnalysisRecordFn = function(...) NULL)

# From test-workflow-integration-input.R:5-24 (R01's fixture).
k228InputFixture <- function() {
    samples <- data.frame(Run = c("s1", "s2", "s3", "s4"), group = c("A", "A", "B", "B"),
        tech_rep_group = c("A_1", "A_2", "B_1", "B_2"), stringsAsFactors = FALSE)
    matrix <- data.frame(feature_id = c("P1", "P2", "P3"), s1 = c(20.125, 21.5, NA), s2 = c(20.25, 21.375, 18.0625),
        s3 = c(22.1, 21.4, 18.3), s4 = c(22.3, 21.6, 1 / 3), stringsAsFactors = FALSE)
    features <- data.frame(feature_id = c("P1", "P2", "P3"), display_name = c("GENE1", "P2", "GENE3"),
        stringsAsFactors = FALSE)
    manifest <- list(omic_type = "proteomics", workflow = "DIA", matrix_role = "normalised", scale = "log2",
        correction = NULL, source_state = list(name = "correlation_filtered", entity_id = strrep("a", 64),
            sha256 = strrep("b", 64)), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group")
    list(manifest = manifest, samples = samples, matrix = matrix, features = features)
}

# ------------------------------------------------------------------------------------------------
# R01 (Re-check 1): the failed publish says a previous input is kept only when there was one.
# ------------------------------------------------------------------------------------------------

test_that("KIT-228 R01: a publish that fails with no previous input does not say a previous input is kept", {
    root <- withr::local_tempdir()
    dir <- file.path(root, "integration", "proteomics")
    fixture <- k228InputFixture()
    write <- function(rename) writeIntegrationInput(dir, fixture$manifest, fixture$samples, list(proteomics = fixture$matrix),
        list(proteomics = fixture$features), renameFn = rename)
    error <- tryCatch(write(function(from, to) FALSE), integration_input_invalid = function(e) e)
    expect_s3_class(error, "integration_input_invalid")
    message <- conditionMessage(error)
    expect_match(message, "the integration input could not be published at", fixed = TRUE, info = message)
    expect_false(grepl("kept", message, fixed = TRUE), info = message)
    expect_false(dir.exists(dir))
    # With a previous input the text is as before (test-workflow-integration-input.R:98 holds it too).
    writeIntegrationInput(dir, fixture$manifest, fixture$samples, list(proteomics = fixture$matrix),
        list(proteomics = fixture$features))
    calls <- 0L
    failPublish <- function(from, to) {
        calls <<- calls + 1L
        if (identical(to, dir) && !grepl("previous", basename(from))) FALSE else file.rename(from, to)
    }
    error <- tryCatch(write(failPublish), integration_input_invalid = function(e) e)
    expect_match(conditionMessage(error), "; the previous input is kept", fixed = TRUE)
    expect_true(file.exists(k228Manifest(dir)))
})

# ------------------------------------------------------------------------------------------------
# E4, G3, G4: proteomics.
# ------------------------------------------------------------------------------------------------

test_that("KIT-228 E4, G1, G3, G4: proteomics, a generic failure removes the earlier input; a removal that fails is named", {
    root <- withr::local_tempdir()
    owner <- k228ProtOwner(root)
    withr::defer(closeDIAStateAdapter(owner$manager))
    dirs <- k207Dirs(file.path(root, "e4"), "proteomics")
    folder <- k228Folder(dirs, "proteomics")
    sha <- k228FirstSave(owner$workflow, "proteomics", dirs, "E4", args = k228ProtWriters)
    # G1: a successful Save after a successful Save.
    k228ExpectSavedAgain(owner$workflow, "proteomics", dirs, "G1", args = k228ProtWriters)
    # E4: the study parameters cannot be written (the generic branch, the fallback file as before).
    failing <- list(createWorkflowArgsFn = function(...) stop("injected: study parameters not written"),
        writeAnalysisRecordFn = function(...) NULL)
    rec <- k228Save(owner$workflow, "proteomics", dirs, args = failing)
    k228Line("E4", "proteomics", rec, folder)
    cat(sprintf("KIT-228 E4 proteomics\tmanifest_sha_before=%s\tafter=%s\n", sha,
        k228ManifestSha(folder)))
    # I04: an in-process error whose reason the Summary never shows: the notice does not point to a message above.
    k228ExpectRemoved(rec, folder, k228Texts$failed_unnamed("proteomics"), "E4")
    expect_true("[warning] Study parameters saved with warnings" %in% rec$notes)
    expect_true(file.exists(file.path(dirs$proteomics$source_dir, "study_parameters.txt")), info = "E4 fallback file")
    expect_true(length(rec$failed) == 1L)
    # G4: the same failure with no earlier folder: no notice, no folder.
    fresh <- k207Dirs(file.path(root, "g4"), "proteomics")
    rec <- k228Save(owner$workflow, "proteomics", fresh, args = failing)
    k228Line("G4", "proteomics", rec, k228Folder(fresh, "proteomics"))
    k228ExpectNoNotice(rec, k228Folder(fresh, "proteomics"), "G4 proteomics")
    # I06 (i): the start removal fails, then the writer's own no-source removal succeeds (the folder is gone at the end):
    # the notice is chosen from the end state, so it says the input was removed (T3), not that it could not be.
    k228FirstSave(owner$workflow, "proteomics", dirs, "I06 prot", args = k228ProtWriters)
    rec <- testthat::with_mocked_bindings(k228Save(owner$workflow, "proteomics", dirs, args = k228ProtWriters),
        removeOmicIntegrationInput = function(integrationDir, omicType) abortIntegrationInput("injected: cannot remove"),
        protIntegrationInputSource = function(manager) NULL, .package = "MultiScholaR")
    k228Line("I06 prot", "proteomics", rec, folder)
    info <- paste("I06 prot | notes:", paste(rec$notes, collapse = " | "))
    expect_false(dir.exists(folder), info = info)
    expect_identical(sum(rec$notes == paste0("[warning] ", k228Texts$no_table("proteomics"))), 1L, info = info)
    expect_false(paste0("[warning] ", k228Texts$not_removed(folder)) %in% rec$notes, info = info)
    # G3: the removal fails (a read-only integration folder, Linux): named, the Save goes on, T4 at its end.
    k228FirstSave(owner$workflow, "proteomics", dirs, "G3", args = k228ProtWriters)
    integration <- k207IntegrationDir(dirs$proteomics)
    Sys.chmod(integration, "555")
    withr::defer(Sys.chmod(integration, "755"))
    skip_if(file.access(integration, 2L) == 0L, "the integration folder is still writable (running as root)")
    rec <- k228Save(owner$workflow, "proteomics", dirs, args = failing)
    k228Line("G3", "proteomics", rec, folder)
    info <- paste("G3 | notes:", paste(rec$notes, collapse = " | "))
    expect_true(any(grepl(folder, rec$notes[startsWith(rec$notes, "[warning] the earlier integration input at")], fixed = TRUE)),
        info = info)
    expect_identical(sum(rec$notes == paste0("[warning] ", k228Texts$not_removed(folder))), 1L, info = info)
    expect_false(any(grepl("was removed", rec$notes, fixed = TRUE)), info = info)
    expect_true("[warning] Study parameters saved with warnings" %in% rec$notes, info = info)
})

# ------------------------------------------------------------------------------------------------
# E1, E2, E5, E6, G1, G4: transcriptomics.
# ------------------------------------------------------------------------------------------------

test_that("KIT-228 E1: transcriptomics, an RUV step run again with no new export: the Save is refused and removes the earlier input", {
    skip_if_not_installed("edgeR")
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("E1 run precondition:", paste(seen$notes, collapse = " | ")))
    kit172Quiet(runTrxReferenceRuviiiExport(owner, "group"))
    dirs <- project$project_dirs
    folder <- k228Folder(dirs, "transcriptomics")
    sha <- k228FirstSave(project$workflow, "transcriptomics", dirs, "E1")
    # G1 on the same project: a successful Save after a successful Save.
    k228ExpectSavedAgain(project$workflow, "transcriptomics", dirs, "G1")
    # A second run (k 1) with no new RUVIII export: the accepted export is the first run's (KIT-207 E9).
    k194TrxRun(project, nd, list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10), seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("E1 second run precondition:", paste(seen$notes, collapse = " | ")))
    expect_false(identical(owner$.context$snapshot$exports$trx_ruviii_export$source$entity_id,
        getWorkflowStateReference(owner, "trx_ruvs")$entity_id), info = "E1 precondition: the accepted export is the first run's")
    rec <- k228Save(project$workflow, "transcriptomics", dirs)
    k228Line("E1", "transcriptomics", rec, folder)
    cat(sprintf("KIT-228 E1 transcriptomics\tmanifest_sha_before=%s\tafter=%s\n", sha,
        k228ManifestSha(folder)))
    k228ExpectRemoved(rec, folder, k228Texts$failed("transcriptomics"), "E1")
    expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "integration_input_invalid"), info = "E1 the refusal's class")
})

test_that("KIT-228 E2, G4: transcriptomics, the RUVIII export missing at the gate removes the earlier input; with none, nothing is said", {
    skip_if_not_installed("edgeR")
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    # A run without RUV (the TMM logCPM, no export needed), saved; then a run with RUV k > 0 and no export (KIT-207 G5).
    kit172Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none")))
    dirs <- project$project_dirs
    folder <- k228Folder(dirs, "transcriptomics")
    sha <- k228FirstSave(project$workflow, "transcriptomics", dirs, "E2")
    stash_e2 <- file.path(withr::local_tempdir(), "stash")
    dir.create(stash_e2)
    file.copy(folder, stash_e2, recursive = TRUE)
    kit172Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("groups_all_features", 2L, "user")))
    expect_null(owner$.context$snapshot$exports$trx_ruviii_export)
    expect_identical(trxSummaryIntegrationRequires(project$workflow), "trx_ruviii_export", info = "E2 gate precondition")
    rec <- k228Save(project$workflow, "transcriptomics", dirs)
    k228Line("E2", "transcriptomics", rec, folder)
    cat(sprintf("KIT-228 E2 transcriptomics\tmanifest_sha_before=%s\tafter=%s\n", sha,
        k228ManifestSha(folder)))
    k228ExpectRemoved(rec, folder, k228Texts$gate("transcriptomics"), "E2")
    missing_message <- omicGuiConfig("transcriptomics")$summary$integration$missing_message
    expect_true(paste0("[warning] ", missing_message) %in% rec$notes, info = "E2 the gate's own warning")
    expect_true("[message] Study parameters saved successfully" %in% rec$notes, info = "E2 study parameters saved")
    expect_identical(rec$done, 1L)
    # C01: the gate while ANOTHER step's publication is pending (the runner retries): the earlier input is removed, T2 shown.
    k228Restore(stash_e2, folder)
    expect_true(file.exists(k228Manifest(folder)), info = "C01 gate precondition")
    rec <- k228Save(project$workflow, "transcriptomics", dirs, runner = k228BusyRunner())
    k228Line("C01 gate", "transcriptomics", rec, folder)
    k228ExpectRemoved(rec, folder, k228Texts$gate("transcriptomics"), "C01 gate")
    # G4: the same gate with no earlier folder (the notes as before the change).
    fresh <- k207Dirs(file.path(withr::local_tempdir(), "g4"), "transcriptomics")
    rec <- k228Save(project$workflow, "transcriptomics", fresh)
    k228Line("G4", "transcriptomics", rec, k228Folder(fresh, "transcriptomics"))
    k228ExpectNoNotice(rec, k228Folder(fresh, "transcriptomics"), "G4 transcriptomics gate")
    expect_true(paste0("[warning] ", missing_message) %in% rec$notes)
})

test_that("KIT-228 E5, E6: transcriptomics on the managed route: a failure that lost its class, and a Save that never ends, remove the earlier input", {
    skip_if_not_installed("edgeR")
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    kit172Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none")))
    dirs <- project$project_dirs
    folder <- k228Folder(dirs, "transcriptomics")
    sha <- k228FirstSave(project$workflow, "transcriptomics", dirs, "E5")
    rec <- k228Save(project$workflow, "transcriptomics", dirs, runner = k228FailingRunner())
    k228Line("E5", "transcriptomics", rec, folder)
    cat(sprintf("KIT-228 E5 transcriptomics\tmanifest_sha_before=%s\tafter=%s\n", sha,
        k228ManifestSha(folder)))
    k228ExpectClassLost(rec, "E5 transcriptomics")
    k228ExpectRemoved(rec, folder, k228Texts$failed("transcriptomics"), "E5 transcriptomics")
    k228FirstSave(project$workflow, "transcriptomics", dirs, "E6")
    rec <- k228Save(project$workflow, "transcriptomics", dirs, runner = k228SilentRunner())
    k228Line("E6", "transcriptomics", rec, folder)
    expect_false(dir.exists(folder), info = paste("E6 transcriptomics | notes:", paste(rec$notes, collapse = " | ")))
    expect_identical(rec$done, 0L)
    expect_identical(length(rec$failed), 0L)
})

# ------------------------------------------------------------------------------------------------
# E3, E5, E6, G1, G2, G4: lipidomics and metabolomics.
# ------------------------------------------------------------------------------------------------

for (class in c("lipid", "metab")) {
    test_that(paste("KIT-228 E3, E5, E6, G1, G2, G4:", class, "a Save that writes no new input removes the earlier one"), {
        skip_if_not(workflowOwnerCapable())
        root <- withr::local_tempdir()
        seen <- k194Notes()
        omic <- k194Omic(class)
        prefix <- kit172Binding(class)$prefix
        q <- k207Normalised(file.path(root, "q"), class, da = TRUE)
        dirs <- k207Dirs(file.path(root, "project"), omic)
        folder <- k228Folder(dirs, omic)
        sha <- k228FirstSave(q$fx$workflow, omic, dirs, "E3")
        # G1: a successful Save after a successful Save.
        k228ExpectSavedAgain(q$fx$workflow, omic, dirs, "G1")
        # E3: the integration source refuses (in process).
        rec <- testthat::with_mocked_bindings(k228Save(q$fx$workflow, omic, dirs),
            omicsIntegrationInputSource = function(binding, manager) abortIntegrationInput("injected"), .package = "MultiScholaR")
        k228Line("E3", omic, rec, folder)
        cat(sprintf("KIT-228 E3 %s\tmanifest_sha_before=%s\tafter=%s\n", omic, sha,
            k228ManifestSha(folder)))
        k228ExpectRemoved(rec, folder, k228Texts$failed(omic), paste("E3", class))
        expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "integration_input_invalid"), info = paste("E3", class))
        # E5: the managed route, the worker's failure with its class lost.
        k228FirstSave(q$fx$workflow, omic, dirs, "E5")
        rec <- k228Save(q$fx$workflow, omic, dirs, runner = k228FailingRunner())
        k228Line("E5", omic, rec, folder)
        k228ExpectClassLost(rec, paste("E5", class))
        k228ExpectRemoved(rec, folder, k228Texts$failed(omic), paste("E5", class))
        # E6: the managed route, a Save that never calls back.
        k228FirstSave(q$fx$workflow, omic, dirs, "E6")
        rec <- k228Save(q$fx$workflow, omic, dirs, runner = k228SilentRunner())
        k228Line("E6", omic, rec, folder)
        expect_false(dir.exists(folder), info = paste("E6", class, "| notes:", paste(rec$notes, collapse = " | ")))
        expect_identical(rec$done, 0L, info = paste("E6", class))
        expect_identical(length(rec$failed), 0L, info = paste("E6", class))
        # I05: the managed route's "publication pending" failure, then Save pressed again: the worker's input must be on disk
        # and what the user is told must be true (the re-press must not delete the input the resumed request published).
        k228FirstSave(q$fx$workflow, omic, dirs, "I05")
        stash <- file.path(root, "stash")
        dir.create(stash)
        file.copy(folder, stash, recursive = TRUE)
        write_input <- function() {
            unlink(folder, recursive = TRUE)
            file.copy(file.path(stash, basename(folder)), dirname(folder), recursive = TRUE)
        }
        retry <- k228RetryRunner(write_input)
        first <- k228Save(q$fx$workflow, omic, dirs, runner = retry)
        k228Line("I05 first press", omic, first, folder)
        expect_true(file.exists(k228Manifest(folder)), info = paste("I05", class, "the worker's input after the first press"))
        again <- k228Save(q$fx$workflow, omic, dirs, runner = retry)
        k228Line("I05 second press", omic, again, folder)
        info <- paste("I05", class, "| first notes:", paste(first$notes, collapse = " | "), "| second notes:",
            paste(again$notes, collapse = " | "))
        # Decision VAL-393 (option (a)): every Save removes the earlier input at its start, the re-press included; the
        # re-press says so once. The resumed request's own "Saved Integration S4 Object" stays (a limit for the wording pass).
        expect_true("[message] Saved Integration S4 Object" %in% first$notes, info = info)
        expect_false(dir.exists(folder), info = info)
        expect_identical(sum(again$notes == paste0("[warning] ", k228Texts$retry(omic))), 1L, info = info)
        # C02: the resumed request is dropped as superseded: no input is left.
        k228FirstSave(q$fx$workflow, omic, dirs, "C02")
        dropped <- k228RetryRunner(write_input, resume = "dropped")
        k228Save(q$fx$workflow, omic, dirs, runner = dropped)
        expect_true(file.exists(k228Manifest(folder)), info = paste("C02", class, "the worker's input after the first press"))
        rec <- k228Save(q$fx$workflow, omic, dirs, runner = dropped)
        k228Line("C02", omic, rec, folder)
        expect_false(dir.exists(folder), info = paste("C02", class, "| notes:", paste(rec$notes, collapse = " | ")))
        # C01: ANOTHER step's publication pending (the runner retries), then a Save that fails: the earlier input is removed
        # and the Save's notice is shown.
        k228FirstSave(q$fx$workflow, omic, dirs, "C01")
        rec <- k228Save(q$fx$workflow, omic, dirs, runner = k228BusyRunner())
        k228Line("C01 failed", omic, rec, folder)
        k228ExpectRemoved(rec, folder, k228Texts$failed(omic), paste("C01 failed", class))
        # I06 (ii): the start removal fails, then the writer with no source removes the folder and reports TRUE (the folder is
        # gone at the end): T3, not silence.
        k228FirstSave(q$fx$workflow, omic, dirs, "I06")
        rec <- testthat::with_mocked_bindings(k228Save(q$fx$workflow, omic, dirs),
            removeOmicIntegrationInput = function(integrationDir, omicType) abortIntegrationInput("injected: cannot remove"),
            omicsIntegrationInputSource = function(binding, manager) NULL, .package = "MultiScholaR")
        k228Line("I06", omic, rec, folder)
        info <- paste("I06", class, "| notes:", paste(rec$notes, collapse = " | "))
        expect_false(dir.exists(folder), info = info)
        expect_identical(sum(rec$notes == paste0("[warning] ", k228Texts$no_table(omic))), 1L, info = info)
        expect_false(paste0("[warning] ", k228Texts$not_removed(folder)) %in% rec$notes, info = info)
        # I07: T1 in the stale-head branch (the commonest refusal) when an earlier folder exists: a Save, a Revert to QC
        # complete, a Save into the same folder.
        k228FirstSave(q$fx$workflow, omic, dirs, "I07")
        stash_i07 <- file.path(root, "stash-i07")
        dir.create(stash_i07)
        file.copy(folder, stash_i07, recursive = TRUE)
        # G2 and G4: KIT-207's refusal (a Revert to QC complete) with no earlier folder: the notes as before, no notice.
        expect_identical(k194Revert(q$fx, paste0(prefix, "_qc_complete")), paste0(prefix, "_qc_complete"), info = "G2 precondition")
        fresh <- k207Dirs(file.path(root, "g2-project"), omic)
        rec <- k228Save(q$fx$workflow, omic, fresh)
        k228Line("G2", omic, rec, k228Folder(fresh, omic))
        k228ExpectNoNotice(rec, k228Folder(fresh, omic), paste("G2 G4", class))
        expect_true(paste0("[error] ", k207Refused(paste("The current table is not the result of the last normalisation: Quality",
            "Control was reverted, or the normalisation was reset, after it ran. Run the normalisation (and the correlation filter",
            "or its Skip) again, or, if you have only just pressed Reset, press Revert once in a Quality Control sub-tab; then",
            "press this button again."))) %in% rec$notes, info = paste("G2", class))
        rec <- k228Save(q$fx$workflow, omic, dirs)
        k228Line("I07", omic, rec, folder)
        k228ExpectRemoved(rec, folder, k228Texts$failed(omic), paste("I07", class))
        expect_true(length(rec$failed) == 1L && inherits(rec$failed[[1L]], "omic_norm_correlation_input_stale"), info = paste("I07", class))
        # C01: the same refusal while ANOTHER step's publication is pending (the runner retries): the earlier input is removed
        # and T1 shown (the refusal comes before the runner is reached).
        k228Restore(stash_i07, folder)
        expect_true(file.exists(k228Manifest(folder)), info = paste("C01 refused precondition", class))
        rec <- k228Save(q$fx$workflow, omic, dirs, runner = k228BusyRunner())
        k228Line("C01 refused", omic, rec, folder)
        k228ExpectRemoved(rec, folder, k228Texts$failed(omic), paste("C01 refused", class))
    })
}
