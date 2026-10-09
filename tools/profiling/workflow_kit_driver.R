# KIT-001 journey driver: the whole DIA journey through the production GUI handlers of a
# managed `gui`-mode application, with per-step phase marks, per-step store snapshots and
# fingerprints. Composed by sourcing the existing browser helpers (workflow_checkpoint_dia_driver.R)
# from the stage's frozen snapshot; no R/ hook and no CORE-pinned file is changed.
#
# The driver runs in the coordinator process, outside the managed slice (as the CORE launcher did);
# the browser it drives is therefore outside the measured slice too. The slice holds the
# application controller, its monitor and any worker.

#' The uptime clock shared with the coordinator's trace.
#' @noRd
kitClock <- function() as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])

#' Open the step-mark file beside the memory trace.
#' @noRd
kitOpenPhases <- function(path) {
    events <- file(path, "wt")
    writeLines("uptime\tstep\tevent", events)
    flush(events)
    events
}

#' Write one uptime-clocked mark.
#' @noRd
kitMark <- function(events, step, event) {
    writeLines(paste(kitClock(), step, event, sep = "\t"), events)
    flush(events)
    invisible(NULL)
}

#' Drive the production design builder for any number of groups (the pinned helper allows two).
#' Factors, metadata assignment and technical replicates follow runCheckpointDesignActions;
#' one contrast is added: `groups[2]` minus `groups[1]`.
#' @noRd
runKitDesignActions <- function(h, driver, design_path) {
    design <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("sample", "group") %in% names(design)), !anyNA(design[c("sample", "group")]),
        !anyDuplicated(design$sample))
    groups <- unique(design$group)
    stopifnot(length(groups) >= 2L)
    prefix <- "proteomics_workflow-design_matrix-builder-"
    set <- function(id, value) h$checkpointSetInput(driver, paste0(prefix, id), value)
    click <- function(id) {
        encoded <- jsonlite::toJSON(paste0(prefix, id), auto_unbox = TRUE)
        h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", encoded),
            timeout = h$checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("document.getElementById(%s).click();true", encoded))
        h$checkpointWaitForIdle(driver, timeout = h$checkpointWaitBudget())
    }
    h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "design")
    set("main_tabs", "Factors")
    for (group in groups) {
        set("new_factor", group)
        click("add_factor")
    }
    set("main_tabs", "Assign Metadata")
    for (group in groups) {
        h$checkpointChooseExisting(driver, paste0(prefix, "selected_runs"),
            design$sample[design$group == group], multiple = TRUE)
        h$checkpointChooseExisting(driver, paste0(prefix, "factor1_select"), group)
        click("assign_metadata")
    }
    if ("replicates" %in% names(design)) {
        stopifnot(!anyNA(design$replicates), all(design$replicates >= 1))
        set("main_tabs", "Technical Replicates")
        set("tech_rep_assignment_mode", "manual")
        for (group in groups) for (replicate in unique(design$replicates[design$group == group])) {
            samples <- design$sample[design$group == group & design$replicates == replicate]
            if (length(samples) > 1L) {
                driver$get_js(h$checkpointDesignAssignmentScript(prefix, samples, replicate, prepare = TRUE))
                expected <- sprintf(paste0("Group: %s, Biological Replicate: %s\n",
                    "  Samples: %s\n  Technical Replicates: %s"), group, replicate,
                    paste(samples, collapse = ", "), paste(seq_along(samples), collapse = ", "))
                h$checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
                    "return !!el && el.textContent.split('\\n\\n').includes(%s) && ",
                    "window.shinytest2 && window.shinytest2.busy!==true;})()"),
                    jsonlite::toJSON(paste0(prefix, "tech_rep_summary"), auto_unbox = TRUE),
                    jsonlite::toJSON(expected, auto_unbox = TRUE)),
                    timeout = h$checkpointWaitBudget(control = TRUE), stable_ms = 200)
            }
        }
    }
    set("main_tabs", "Contrasts")
    h$checkpointChooseExisting(driver, paste0(prefix, "contrast_group1"), groups[[2L]])
    h$checkpointChooseExisting(driver, paste0(prefix, "contrast_group2"), groups[[1L]])
    click("add_contrast")
    click("save_results")
    h$waitCheckpointStepComplete(driver, "Design", timeout = h$checkpointWaitBudget())
    list(groups = groups, contrast = paste(groups[[2L]], "vs", groups[[1L]]))
}

#' Click one QC action and wait for its result text or an error notification (as the pinned helpers).
#' @noRd
runKitQcAction <- function(h, driver, prefix, tab_input, tab, action, result, state, settings = list()) {
    h$checkpointSetInput(driver, paste0(prefix, tab_input), tab)
    for (field in names(settings)) h$checkpointSetInput(driver, paste0(prefix, field), settings[[field]])
    id <- jsonlite::toJSON(paste0(prefix, action), auto_unbox = TRUE)
    h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", id),
        timeout = h$checkpointWaitBudget(control = TRUE))
    driver$get_js(sprintf("document.getElementById(%s).click();true", id))
    output <- jsonlite::toJSON(paste0(prefix, result), auto_unbox = TRUE)
    encoded <- jsonlite::toJSON(state, auto_unbox = TRUE)
    h$checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "return (el && el.textContent.includes(%s)) || ",
        "document.querySelector('.shiny-notification-error') !== null;})()"), output, encoded),
        timeout = h$checkpointWaitBudget())
    error <- driver$get_js(paste0("(function(){const el=document.querySelector(",
        "'.shiny-notification-error');return el?el.textContent:null;})()"))
    if (!is.null(error)) stop(paste("production action failed:", state, error))
    invisible(state)
}

#' Snapshot the adopted store after a step: head hash, current state, entity id and rows.
#' NULL while no store exists (an in-memory owner has no references to snapshot).
#' @noRd
kitSnapshotAfterStep <- function(project, api) {
    store_dir <- kitFindStore(project)
    if (is.null(store_dir)) return(NULL)
    head <- kitReadHead(store_dir, api)
    list(store = store_dir, head_sha256 = head$head_sha256, state = head$current_state,
        entity_id = head$states[[head$current_state]], rows = head$rows[[head$current_state]],
        history = head$history)
}

#' Run the journey. Every step is marked `start` then `complete` or `failed`; after each step the
#' store is snapshotted; the fault (lifecycle failure member) injects the application's own guard
#' stop after `fault$after_step` and expects `fault$failed_step` to fail with the head unchanged.
#' @param h Environment holding the pinned browser helpers.
#' @param driver Managed browser driver attached to the application.
#' @param project_parent Directory the GUI creates the project under.
#' @param corpus `list(report, fasta, design)` paths.
#' @param journey Frozen journey settings from the bindings.
#' @param fault NULL or `list(name, after_step, failed_step)`.
#' @param application_root Managed application root (for the guard injection file).
#' @param events Open phase-mark connection.
#' @param api Accessor into the installed namespace.
#' @return Journey record: populations per step, files, final head, status.
#' @noRd
runKitJourney <- function(h, driver, project_parent, corpus, journey, fault, application_root, events, api) {
    steps <- kitJourneySteps()
    record <- list(status = "running", populations = list(), completed = character(), store = NULL,
        corpus = journey$corpus_key, project = NULL, files = NULL)
    # Everything the GUI writes lands under the member's project parent; that is the project root.
    project <- function() project_parent
    settings <- journey$settings
    peptide_prefix <- "proteomics_workflow-quality_control-peptide_qc-"
    protein_prefix <- "proteomics_workflow-quality_control-protein_qc-"
    norm_prefix <- "proteomics_workflow-normalization-"
    peptide_actions <- stats::setNames(h$checkpointPeptideActions(),
        c("qvalue", "protein_peptide", "rollup", "intensity", "sample", "replicate", "imputation"))
    protein_actions <- list(
        iq = list(tab = "IQ Protein Rollup", action = "rollup-apply_iq_rollup",
            result = "rollup-iq_rollup_results", state = "protein_s4_created"),
        accession = list(tab = "Accession Cleanup", action = "cleanup-apply_accession_cleanup",
            result = "cleanup-accession_cleanup_results", state = "protein_accession_cleaned"),
        protein_intensity = list(tab = "Protein Intensity Filter", action = "intensity_filter-apply_protein_intensity_filter",
            result = "intensity_filter-protein_intensity_filter_results", state = "protein_intensity_filtered"),
        duplicates = list(tab = "Duplicate Removal", action = "duplicate_removal-apply_duplicate_removal",
            result = "duplicate_removal-duplicate_removal_results", state = "duplicates_removed"),
        protein_replicate = list(tab = "Protein Replicate Filter", action = "replicate_filter-apply_protein_replicate_filter",
            result = "replicate_filter-protein_replicate_filter_results", state = "protein_replicate_filtered"))
    run_step <- function(step) {
        id <- step$id
        switch(step$kind,
            store = switch(id,
                import = h$runCheckpointImportActions(driver, project_parent, corpus),
                design = runKitDesignActions(h, driver, corpus$design),
                session_export = h$checkpointClickAndWaitText(driver, paste0(norm_prefix, "export_filtered_session"),
                    ".shiny-notification", "Filtered session data exported successfully!")),
            peptide = {
                if (id == "qvalue") {
                    h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "qc")
                    h$checkpointSetInput(driver, "proteomics_workflow-quality_control-qc_tabs_lfq", "Peptide QC")
                }
                action <- peptide_actions[[id]]
                runKitQcAction(h, driver, peptide_prefix, "peptide_filter_tabs", action$tab, action$action,
                    action$result, action$state, settings[[id]] %||% list())
            },
            protein = {
                if (id == "iq") h$checkpointSetInput(driver, "proteomics_workflow-quality_control-qc_tabs_lfq", "Protein QC")
                action <- protein_actions[[id]]
                runKitQcAction(h, driver, protein_prefix, "protein_filter_tabs", action$tab, action$action,
                    action$result, action$state, settings[[id]] %||% list())
                if (id == "protein_replicate") h$waitCheckpointStepComplete(driver, "QC", timeout = h$checkpointWaitBudget())
            },
            normalisation = switch(id,
                normalisation_ruv = {
                    h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "normalization")
                    h$checkpointSetInput(driver, paste0(norm_prefix, "norm_method"), journey$normalisation$norm_method)
                    h$checkpointSetInput(driver, paste0(norm_prefix, "ruv_mode"), journey$normalisation$ruv_mode)
                    for (name in names(journey$normalisation$ruv)) {
                        h$checkpointSetInput(driver, paste0(norm_prefix, name), journey$normalisation$ruv[[name]])
                    }
                    h$checkpointClickAndWaitText(driver, paste0(norm_prefix, "run_normalization"),
                        ".shiny-notification", "Normalization and RUV correction completed!")
                },
                correlation = {
                    h$checkpointSetInput(driver, paste0(norm_prefix, "norm_qc_tabs"), "Correlation Filtering")
                    h$checkpointSetInput(driver, paste0(norm_prefix, "min_pearson_correlation_threshold"),
                        journey$normalisation$correlation_threshold)
                    h$checkpointClickAndWaitText(driver, paste0(norm_prefix, "apply_correlation_filter"),
                        paste0("#", norm_prefix, "correlation_filter_summary"), "Correlation filtering completed successfully!")
                }),
            da = {
                prefix <- "proteomics_workflow-differential_expression-"
                h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "da")
                h$checkpointClickAndWaitText(driver, paste0(prefix, "load_filtered_session"),
                    ".shiny-notification", "Session loaded successfully!")
                h$checkpointSetInput(driver, paste0(prefix, "formula_string"), journey$da$formula)
                h$checkpointSetInput(driver, paste0(prefix, "da_q_val_thresh"), journey$da$q_value)
                h$checkpointSetInput(driver, paste0(prefix, "treat_lfc_cutoff"), journey$da$treat_lfc)
                h$checkpointClickAndWaitText(driver, paste0(prefix, "run_da_analysis"),
                    paste0("#", prefix, "da_status"), "[OK] Analysis Complete")
            })
        invisible(NULL)
    }
    previous_files <- NULL
    for (step in steps) {
        kitMark(events, step$id, "start")
        outcome <- tryCatch({ run_step(step); NULL }, error = function(error) error)
        if (!is.null(outcome)) {
            kitMark(events, step$id, "failed")
            record$status <- "journey_failed"
            record$failed_step <- step$id
            record$error <- conditionMessage(outcome)
            break
        }
        kitMark(events, step$id, "complete")
        record$completed <- c(record$completed, step$id)
        record$project <- project()
        snapshot <- kitSnapshotAfterStep(record$project, api)
        files <- kitFingerprintFiles(record$project, previous_files)
        changed <- if (is.null(previous_files)) files else
            files[!(paste(files$path, files$sha256) %in% paste(previous_files$path, previous_files$sha256)), , drop = FALSE]
        previous_files <- files
        record$populations[[step$id]] <- c(list(step = step$id, expected_state = step$state,
            clock = kitClock(), files_changed = changed$path), snapshot)
        if (!is.null(snapshot)) record$store <- snapshot$store
        if (!is.null(fault) && identical(step$id, fault$after_step)) {
            # The application's own guard path: the monitor stops the slice on this request.
            saveRDS(list(schema = 1L, reason = "injected_guard_stop"),
                file.path(application_root, "application-guard-inject.rds"), compress = FALSE, version = 3)
            kitMark(events, "guard_injected", "start")
        }
    }
    if (record$status == "running") record$status <- "complete"
    record$files <- previous_files
    fingerprints <- if (is.null(record$store)) NULL else tryCatch(kitFingerprintStore(record$store, api),
        error = function(error) { record$fingerprint_error <<- conditionMessage(error); NULL })
    if (!is.null(fingerprints)) {
        # Fingerprint every state the final head references and attach each to the step that
        # produced it (a state snapshotted only after a later step is marked after_adoption).
        record$final_head_sha256 <- fingerprints$head_sha256
        record$product_fingerprints <- fingerprints$products
        head <- kitReadHead(record$store, api)
        for (id in names(record$populations)) {
            state <- record$populations[[id]]$expected_state
            if (is.null(state) || is.null(fingerprints$states[[state]])) next
            if (is.null(record$populations[[id]]$state)) {
                record$populations[[id]]$state <- state
                record$populations[[id]]$store <- record$store
                record$populations[[id]]$entity_id <- head$states[[state]]
                record$populations[[id]]$rows <- head$rows[[state]]
                record$populations[[id]]$snapshot_time <- "after_adoption"
            } else record$populations[[id]]$snapshot_time <- "after_step"
            record$populations[[id]]$fingerprint <- fingerprints$states[[state]]
        }
        record$populations_complete <- all(vapply(steps, function(step) is.null(step$state) ||
            !is.null(record$populations[[step$id]]$fingerprint), logical(1)))
    } else record$populations_complete <- FALSE
    record
}

#' Attach a browser to the managed application and run the journey, recording the result.
#' @param h Browser helper environment (the pinned driver file sourced from the snapshot).
#' @param member_dir Member directory; `application_root` its managed application root.
#' @noRd
runKitMember <- function(h, member_dir, application_root, library_path, corpus, journey, fault, api) {
    config <- readRDS(file.path(application_root, "application.rds"))
    withr::local_options(multischolar.dia_memory_audit_directory = member_dir,
        multischolar.dia_memory_audit_seconds = 900)
    # The phase file exists before the browser so a driver failure still leaves a readable record.
    events <- kitOpenPhases(file.path(member_dir, "journey-phases.tsv"))
    on.exit(close(events), add = TRUE)
    managed_url <- sprintf("http://127.0.0.1:%d", config$port)
    browser <- h$newCheckpointBrowser(library_path, member_dir, managed_url = managed_url)
    on.exit(try(browser$chrome$close(), silent = TRUE), add = TRUE)
    h$stabilizeCheckpointManagedBrowser(browser$driver, managed_url)
    h$checkpointWaitForJs(browser$driver, "document.querySelector('[data-testid=\"tile-proteomics\"]') !== null",
        timeout = 60000, stable_ms = 200)
    project_parent <- file.path(member_dir, "project")
    dir.create(project_parent)
    record <- tryCatch(runKitJourney(h, browser$driver, project_parent, corpus, journey, fault, application_root,
        events, api), error = function(error) list(status = "driver_failed", error = conditionMessage(error),
            corpus = journey$corpus_key, project = project_parent))
    kitMark(events, "journey", "teardown")
    record$teardown_clock <- kitClock()
    saveRDS(record, file.path(member_dir, "journey-result.rds"), compress = FALSE, version = 3)
    record
}
