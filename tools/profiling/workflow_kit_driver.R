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
runKitDesignActions <- function(h, driver, design_path, prefix = "proteomics_workflow-design_matrix-builder-",
    tab_input = "proteomics_workflow-workflow_tabs", tab_value = "design", save = NULL) {
    # KIT-068: the lipidomics and metabolomics design builders share the proteomics builder's input ids under their own
    # prefix and tab input; `save` clicks the save button and waits on the module's own completion text (the DIA
    # builder's completion is the workflow progress widget, the default below).
    design <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("sample", "group") %in% names(design)), !anyNA(design[c("sample", "group")]),
        !anyDuplicated(design$sample))
    groups <- unique(design$group)
    stopifnot(length(groups) >= 2L)
    set <- function(id, value) h$checkpointSetInput(driver, paste0(prefix, id), value)
    click <- function(id) {
        encoded <- jsonlite::toJSON(paste0(prefix, id), auto_unbox = TRUE)
        h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", encoded),
            timeout = h$checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("document.getElementById(%s).click();true", encoded))
        h$checkpointWaitForIdle(driver, timeout = h$checkpointWaitBudget())
    }
    h$checkpointSetInput(driver, tab_input, tab_value)
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
    if (is.null(save)) {
        click("save_results")
        h$waitCheckpointStepComplete(driver, "Design", timeout = h$checkpointWaitBudget())
    } else save(paste0(prefix, "save_results"))
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
    # `head_entity`: the checkpoint entity whose descriptor is this head, so the step stays reachable
    # after the live head moves on (KIT-005 audit PH-01).
    list(store = store_dir, head_sha256 = head$head_sha256, head_entity = kitFindHeadEntity(store_dir, head$head_sha256),
        state = head$current_state,
        entity_id = head$states[[head$current_state]], rows = head$rows[[head$current_state]],
        history = head$history, product_ids = as.character(unlist(head$products)))
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
    steps <- kitJourneySteps(kitRouteOf(journey), journey$normalisation$ruv_mode %||% "covariates")
    # KIT-018: a single journey member runs import through `journey$stop_after` only.
    if (!is.null(journey$stop_after)) {
        last <- match(journey$stop_after, vapply(steps, `[[`, character(1), "id"))
        stopifnot(!is.na(last))
        steps <- steps[seq_len(last)]
    }
    record <- list(status = "running", populations = list(), completed = character(), store = NULL,
        corpus = journey$corpus_key, project = NULL, files = NULL, stop_after = journey$stop_after)
    # Everything the GUI writes lands under the member's project parent; that is the project root.
    project <- function() project_parent
    settings <- journey$settings
    if (identical(kitRouteOf(journey), "DIA_LIMPA")) {
        # KIT-053: the route's between-sample normalisation step takes the bound tab setting (the ticket's case binding:
        # none on the screen corpus, cyclic loess on the full corpus) through the module's own select before its apply.
        settings$protein_normalisation_limpa <- utils::modifyList(settings$protein_normalisation_limpa %||% list(),
            list("normalisation-norm_method" = journey$normalisation$norm_method %||% "none"))
    }
    peptide_prefix <- "proteomics_workflow-quality_control-peptide_qc-"
    protein_prefix <- "proteomics_workflow-quality_control-protein_qc-"
    norm_prefix <- "proteomics_workflow-normalization-"
    peptide_actions <- stats::setNames(h$checkpointPeptideActions(),
        c("qvalue", "protein_peptide", "rollup", "intensity", "sample", "replicate", "imputation"))
    # KIT-053: the DIA_LIMPA route's own tabs (KIT-047, KIT-041) beside the shared ones.
    peptide_actions$peptide_log2 <- list(tab = "Log2 Matrix", action = "log2_matrix-apply_log2_matrix",
        result = "log2_matrix-log2_results", state = "log2_matrix")
    peptide_actions$peptide_accession <- list(tab = "Peptide Accession Cleanup", action = "accession_cleanup-apply_peptide_accession",
        result = "accession_cleanup-peptide_accession_results", state = "peptide_accession_cleaned")
    protein_actions <- list(
        iq = list(tab = "IQ Protein Rollup", action = "rollup-apply_iq_rollup",
            result = "rollup-iq_rollup_results", state = "protein_s4_created"),
        dpc_quant = list(tab = "DPC Quantification", action = "dpc_quant-apply_dpc_quant",
            result = "dpc_quant-dpc_quant_results", state = "protein_dpc_quant"),
        protein_normalisation_limpa = list(tab = "Between-Sample Normalisation", action = "normalisation-apply_protein_normalisation",
            result = "normalisation-protein_normalisation_results", state = "protein_normalised"),
        peptide_count = list(tab = "Peptide Count Filter", action = "peptide_count-apply_peptide_count",
            result = "peptide_count-peptide_count_results", state = "peptide_count_filtered"),
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
        # KIT-068: the lipidomics and metabolomics routes run their own modules' tab actions.
        if (kitRouteOf(journey) %in% c("lipidomics", "metabolomics"))
            return(runKitOmicsStep(h, driver, project_parent, corpus, journey, step))
        # KIT-098: the transcriptomics route runs its own module ids (kit-098/plan.md section 2).
        if (identical(kitRouteOf(journey), "transcriptomics")) return(runKitTrxStep(h, driver, project_parent, corpus, journey, step))
        switch(step$kind,
            store = switch(id,
                # KIT-053: the limpa route is chosen at the entry point before the import (the pinned import action
                # knows no route; the kit's own import sets the radio and then runs the same actions).
                import = if (identical(kitRouteOf(journey), "DIA_LIMPA")) runKitLimpaImportActions(h, driver, project_parent, corpus) else
                    h$runCheckpointImportActions(driver, project_parent, corpus),
                design = runKitDesignActions(h, driver, corpus$design),
                session_export = kitClickAndWaitText(h, driver, paste0(norm_prefix, "export_filtered_session"),
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
                # The first protein step opens the Protein QC tab: the IQ rollup on the DIA route; on the protein route
                # (KIT-021) there is no peptide stage, so the QC page is opened here too.
                if (id == steps[[match("protein", vapply(steps, `[[`, "", "kind"))]]$id) {
                    h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "qc")
                    # The QC page's tabset is `qc_tabs_lfq` (Peptide QC, Protein QC) on the DIA route and `qc_tabs_tmt`
                    # (Protein QC alone, the design having committed protein_s4_initial) on LFQ/TMT (mod_prot_qc.R).
                    tabset <- if (identical(kitRouteOf(journey), "protein")) "qc_tabs_tmt" else "qc_tabs_lfq"
                    h$checkpointSetInput(driver, paste0("proteomics_workflow-quality_control-", tabset), "Protein QC")
                }
                action <- protein_actions[[id]]
                runKitQcAction(h, driver, protein_prefix, "protein_filter_tabs", action$tab, action$action,
                    action$result, action$state, settings[[id]] %||% list())
                if (id == "protein_replicate") h$waitCheckpointStepComplete(driver, "QC", timeout = h$checkpointWaitBudget())
            },
            normalisation = switch(id,
                ruv_covariates = {
                    # KIT-053: the DIA_LIMPA route's RUV covariates, the tab's own action after the correlation filter.
                    # KIT-056: the bound mode is set on the tab's radio first; the state follows it.
                    limpa_ruv_mode <- journey$normalisation$ruv_mode %||% "covariates"
                    if (!identical(limpa_ruv_mode, "covariates")) h$checkpointSetInput(driver, paste0(norm_prefix, "ruv_mode"), limpa_ruv_mode)
                    h$checkpointSetInput(driver, paste0(norm_prefix, "norm_qc_tabs"), "Correlation Filtering")
                    h$checkpointClickAndWaitText(driver, paste0(norm_prefix, "apply_ruv_covariates"),
                        paste0("#", norm_prefix, "ruv_covariates_summary"),
                        sprintf("State saved as: '%s'", steps[[match(id, vapply(steps, `[[`, "", "id"))]]$state))
                },
                normalisation_ruv = {
                    h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "normalization")
                    h$checkpointSetInput(driver, paste0(norm_prefix, "norm_method"), journey$normalisation$norm_method)
                    h$checkpointSetInput(driver, paste0(norm_prefix, "ruv_mode"), journey$normalisation$ruv_mode)
                    for (name in names(journey$normalisation$ruv)) {
                        h$checkpointSetInput(driver, paste0(norm_prefix, name), journey$normalisation$ruv[[name]])
                    }
                    # The module's completion text names the RUV mode (buildProtNormCompletionNotification): the fixture
                    # corpora skip RUV (KIT-021); the 300-sample corpora and the DIA route run it.
                    kitClickAndWaitText(h, driver, paste0(norm_prefix, "run_normalization"), ".shiny-notification",
                        if (identical(journey$normalisation$ruv_mode, "skip")) "Normalization completed (RUV skipped)!" else
                            "Normalization and RUV correction completed!")
                },
                correlation = {
                    if (identical(kitRouteOf(journey), "DIA_LIMPA")) {
                        # KIT-053: the route's normalisation button is a pass-through (its normalisation is a protein QC
                        # step); it designates the correlation input and completes before the filter runs.
                        h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "normalization")
                        h$checkpointSetInput(driver, paste0(norm_prefix, "ruv_mode"), journey$normalisation$ruv_mode %||% "covariates")
                        kitClickAndWaitText(h, driver, paste0(norm_prefix, "run_normalization"), ".shiny-notification",
                            "Normalisation stage ready (DIA with limpa)")
                    }
                    h$checkpointSetInput(driver, paste0(norm_prefix, "norm_qc_tabs"), "Correlation Filtering")
                    h$checkpointSetInput(driver, paste0(norm_prefix, "min_pearson_correlation_threshold"),
                        journey$normalisation$correlation_threshold)
                    h$checkpointClickAndWaitText(driver, paste0(norm_prefix, "apply_correlation_filter"),
                        paste0("#", norm_prefix, "correlation_filter_summary"), "Correlation filtering completed successfully!")
                }),
            da = {
                prefix <- "proteomics_workflow-differential_expression-"
                h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "da")
                kitClickAndWaitText(h, driver, paste0(prefix, "load_filtered_session"),
                    ".shiny-notification", "Session loaded successfully!")
                h$checkpointSetInput(driver, paste0(prefix, "formula_string"), journey$da$formula)
                h$checkpointSetInput(driver, paste0(prefix, "da_q_val_thresh"), journey$da$q_value)
                h$checkpointSetInput(driver, paste0(prefix, "treat_lfc_cutoff"), journey$da$treat_lfc)
                h$checkpointClickAndWaitText(driver, paste0(prefix, "run_da_analysis"),
                    paste0("#", prefix, "da_status"), "[OK] Analysis Complete")
            },
            report = {
                # KIT-065: the summary tab: the publication copy (which bootstraps the workflow arguments when they were
                # not saved), then the report render with the route's shipped template.
                prefix <- "proteomics_workflow-session_summary-"
                h$checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "session_summary")
                h$checkpointSetInput(driver, paste0(prefix, "experiment_label"), journey$experiment_label)
                h$checkpointClickAndWaitText(driver, paste0(prefix, "copy_to_publication"),
                    paste0("#", prefix, "copy_status"), "Files copied to publication directory successfully")
                kitClickAndWaitText(h, driver, paste0(prefix, "generate_report"), ".shiny-notification",
                    "Report generated successfully!")
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
    # The journey is over: let the application settle, mark teardown (the resident-end window is the
    # 5 s before it) and release the sampler before the fingerprint pass, which reads the immutable
    # store and counts against neither elapsed nor the ceiling.
    Sys.sleep(6)
    kitMark(events, "journey", "teardown")
    record$teardown_clock <- kitClock()
    file.create(file.path(dirname(application_root), "stop-sampling"))
    fingerprints <- if (is.null(record$store)) NULL else tryCatch(kitFingerprintStore(record$store, api),
        error = function(error) { record$fingerprint_error <<- conditionMessage(error); NULL })
    if (!is.null(fingerprints)) {
        # Fingerprint every state the final head references and attach each to the step that
        # produced it (a state snapshotted only after a later step is marked after_adoption).
        record$final_head_sha256 <- fingerprints$head_sha256
        record$product_fingerprints <- fingerprints$products
        head <- kitReadHead(record$store, api)
        for (id in names(record$populations)) {
            # The state the step actually left current (in-run snapshot), else the expected state.
            state <- record$populations[[id]]$state %||% record$populations[[id]]$expected_state
            if (is.null(state) || is.null(fingerprints$states[[state]])) next
            record$populations[[id]]$state_as_expected <- identical(state, record$populations[[id]]$expected_state)
            if (is.null(record$populations[[id]]$state)) {
                record$populations[[id]]$state <- state
                record$populations[[id]]$store <- record$store
                record$populations[[id]]$entity_id <- head$states[[state]]
                record$populations[[id]]$rows <- head$rows[[state]]
                record$populations[[id]]$snapshot_time <- "after_adoption"
            } else record$populations[[id]]$snapshot_time <- "after_step"
            record$populations[[id]]$fingerprint <- fingerprints$states[[state]]
        }
        # Products are attached to the step after which they first appear (canonical fingerprints,
        # ordered by hash so the set compares across runs).
        seen <- character()
        for (id in names(record$populations)) {
            ids <- record$populations[[id]]$product_ids %||% character()
            new_ids <- setdiff(ids, seen)
            seen <- union(seen, ids)
            products <- fingerprints$products[intersect(new_ids, names(fingerprints$products))]
            products <- unname(products[order(vapply(products, `[[`, character(1), "sha256"))])
            names(products) <- vapply(products, function(p) paste0(p$product_type %||% p$kind, ":", substr(p$sha256, 1L, 12L)), character(1))
            if (!is.null(record$populations[[id]]$fingerprint)) record$populations[[id]]$fingerprint$products <- c(
                record$populations[[id]]$fingerprint$products, products)
        }
        record$populations_complete <- all(vapply(steps, function(step) is.null(step$state) ||
            !is.null(record$populations[[step$id]]$fingerprint), logical(1)))
    } else record$populations_complete <- FALSE
    record
}

#' A browser tab rendering a large production table freezes for longer than one evaluation
#' call (the pinned helpers allow 5 s per poll and 30 s per click); at 300 samples the IQ
#' protein table did so for over 30 s (KIT-005 correction 4). A transport timeout while the
#' renderer is busy means "not ready yet", so the evaluation is retried, for up to ten minutes
#' per call, until the tab answers; every other error is raised unchanged.
#' @noRd
kitPatientDriver <- function(driver, patience_seconds = 600) {
    original <- driver$get_js
    driver$get_js <- function(script, timeout = NULL) {
        started <- proc.time()[["elapsed"]]
        repeat {
            outcome <- tryCatch(list(value = original(script, timeout)), error = function(error) error)
            if (!inherits(outcome, "error")) return(outcome$value)
            frozen <- grepl("timed out waiting for response", conditionMessage(outcome), fixed = TRUE)
            if (!frozen || proc.time()[["elapsed"]] - started > patience_seconds) stop(outcome)
            Sys.sleep(1)
        }
    }
    invisible(driver)
}

#' Record every Shiny notification the page ever shows in `window.__kitNotifications` (text and
#' error flag), so a completion message that stays visible for 10 s cannot be missed while the
#' renderer is frozen (KIT-005 correction 5). Idempotent.
#' @noRd
kitInstallNotificationLog <- function(driver) {
    driver$get_js(paste0("(function(){if(window.__kitNotifications)return true;window.__kitNotifications=[];",
        "const push=n=>{if(n&&n.textContent)window.__kitNotifications.push({text:n.textContent,",
        "error:!!(n.classList&&n.classList.contains('shiny-notification-error'))});};",
        "const visit=n=>{if(n.nodeType!==1)return;if(n.classList&&n.classList.contains('shiny-notification'))push(n);",
        "if(n.querySelectorAll)n.querySelectorAll('.shiny-notification').forEach(push);};",
        "new MutationObserver(ms=>{for(const m of ms)for(const n of m.addedNodes)visit(n);})",
        ".observe(document.documentElement,{childList:true,subtree:true});return true;})()"))
    invisible(NULL)
}

#' Click a production action and wait for its completion text in the live DOM or in the
#' notification log; an error notification (live or logged) fails the step.
#' @noRd
kitClickAndWaitText <- function(h, driver, input_id, selector, text) {
    kitInstallNotificationLog(driver)
    id <- jsonlite::toJSON(input_id, auto_unbox = TRUE)
    encoded_selector <- jsonlite::toJSON(selector, auto_unbox = TRUE)
    encoded_text <- jsonlite::toJSON(text, auto_unbox = TRUE)
    h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", id),
        timeout = h$checkpointWaitBudget(control = TRUE))
    driver$get_js(sprintf("window.__kitNotifications=[];document.getElementById(%s).click();true", id))
    h$checkpointWaitForJs(driver, sprintf(paste0("(function(){return Array.from(document.querySelectorAll(%s))",
        ".some(el=>el.textContent.includes(%s)) || (window.__kitNotifications||[]).some(n=>n.text.includes(%s)) || ",
        "document.querySelector('.shiny-notification-error') !== null || ",
        "(window.__kitNotifications||[]).some(n=>n.error);})()"), encoded_selector, encoded_text, encoded_text),
        timeout = h$checkpointWaitBudget())
    error <- driver$get_js(paste0("(function(){const el=document.querySelector('.shiny-notification-error');",
        "if(el)return el.textContent;const logged=(window.__kitNotifications||[]).find(n=>n.error);",
        "return logged?logged.text:null;})()"))
    if (!is.null(error)) stop(paste("production action failed:", input_id, error))
    invisible(TRUE)
}

#' Attach a browser to the managed application and run the journey, recording the result.
#' @param h Browser helper environment (the pinned driver file sourced from the snapshot).
#' @param member_dir Member directory; `application_root` its managed application root.
#' @noRd
runKitMember <- function(h, member_dir, application_root, library_path, corpus, journey, fault, api) {
    config <- readRDS(file.path(application_root, "application.rds"))
    # The pinned helpers' wait budget is the diagnostic ceiling minus this process's elapsed time; the
    # member runs in its own process (workflow_kit_member.R) and the 2 h ceiling covers the journey
    # ceiling (90 min), so no browser wait expires before the member's own deadline (KIT-005 correction 3).
    withr::local_options(multischolar.dia_memory_audit_directory = member_dir,
        multischolar.dia_memory_audit_seconds = 7200, multischolar.dia_extended120 = TRUE)
    # The phase file exists before the browser so a driver failure still leaves a readable record.
    events <- kitOpenPhases(file.path(member_dir, "journey-phases.tsv"))
    on.exit(close(events), add = TRUE)
    managed_url <- sprintf("http://127.0.0.1:%d", config$port)
    browser <- h$newCheckpointBrowser(library_path, member_dir, managed_url = managed_url)
    on.exit(try(browser$chrome$close(), silent = TRUE), add = TRUE)
    kitPatientDriver(browser$driver)
    kitInstallNotificationLog(browser$driver)
    h$stabilizeCheckpointManagedBrowser(browser$driver, managed_url)
    h$checkpointWaitForJs(browser$driver, "document.querySelector('[data-testid=\"tile-proteomics\"]') !== null",
        timeout = 60000, stable_ms = 200)
    project_parent <- file.path(member_dir, "project")
    dir.create(project_parent)
    record <- tryCatch(runKitJourney(h, browser$driver, project_parent, corpus, journey, fault, application_root,
        events, api), error = function(error) list(status = "driver_failed", error = conditionMessage(error),
            corpus = journey$corpus_key, project = project_parent))
    if (is.null(record$teardown_clock)) {
        kitMark(events, "journey", "teardown")
        record$teardown_clock <- kitClock()
    }
    saveRDS(record, file.path(member_dir, "journey-result.rds"), compress = FALSE, version = 3)
    record
}


#' The limpa route's import (KIT-053): the pinned import actions (`runCheckpointImportActions`) with the DIA-NN route
#' radio set to `DIA_LIMPA` before the import runs; the pinned helper is not edited.
#' @noRd
runKitLimpaImportActions <- function(h, driver, run_dir, inputs) {
    h$checkpointClickControl(driver, "tile-proteomics")
    h$checkpointClickControl(driver, "btn-start-analysis")
    driver$set_inputs(experiment_label = "dia_baseline_smoke", project_dir = run_dir)
    h$checkpointClickControl(driver, "btn-confirm-setup")
    h$checkpointWaitForJs(driver, "document.querySelector('[data-testid=\"prot-tab-setup\"]') !== null",
        timeout = h$checkpointWaitBudget())
    roots <- shinyFiles::getVolumes()()
    for (field in c("report", "fasta")) {
        input_name <- if (field == "report") "search_results" else "fasta_file"
        id <- paste0("proteomics_workflow-setup_import-", input_name)
        selection <- h$checkpointFileSelection(inputs[[field]], roots)
        h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null",
            jsonlite::toJSON(id, auto_unbox = TRUE)), timeout = h$checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("Shiny.setInputValue(%s,%s,{priority:'event'});true",
            jsonlite::toJSON(id, auto_unbox = TRUE),
            jsonlite::toJSON(selection, auto_unbox = TRUE)))
        h$checkpointWaitForIdle(driver, timeout = h$checkpointWaitBudget())
    }
    args <- list(inputs$format %||% "diann")
    names(args) <- "proteomics_workflow-setup_import-format_override"
    do.call(driver$set_inputs, args)
    h$checkpointSetInput(driver, "proteomics_workflow-setup_import-dia_route", "DIA_LIMPA")
    h$checkpointClickControl(driver, "prot-import-process")
    h$waitCheckpointStepComplete(driver, "Import", timeout = h$checkpointWaitBudget())
    invisible(NULL)
}

#' KIT-068: the lipidomics or metabolomics import: the omic's tile, the project setup, the import tab's vendor radio,
#' the assay names and files (the module's own shinyFiles buttons, selected as the pinned helpers select the DIA report),
#' sanitised names off (the corpus names are already clean and the design names them as written), then the process
#' button, waiting on the module's own completion notification.
#' @param corpus `list(assays = <named list of paths by assay name>, vendor, design)` from the coordinator.
#' @noRd
runKitOmicsImportActions <- function(h, driver, run_dir, corpus, route) {
    # KIT-098: lipidomics or metabolomics only; any other route is refused, never driven as metabolomics.
    short <- switch(route, lipidomics = "lipid", metabolomics = "metab", stop("the omics import drives lipidomics or metabolomics only, not ", route))
    prefix <- paste0(route, "_workflow-import-")
    stopifnot(length(corpus$assays) >= 1L, length(corpus$assays) <= 2L, !is.null(names(corpus$assays)))
    h$checkpointClickControl(driver, paste0("tile-", route))
    h$checkpointClickControl(driver, "btn-start-analysis")
    driver$set_inputs(experiment_label = "kit_omics_lifecycle", project_dir = run_dir)
    h$checkpointClickControl(driver, "btn-confirm-setup")
    h$checkpointWaitForJs(driver, sprintf("document.querySelector('[data-testid=\"%s-tab-import\"]') !== null", short),
        timeout = h$checkpointWaitBudget())
    h$checkpointSetInput(driver, paste0(prefix, "vendor_format"), corpus$vendor)
    roots <- shinyFiles::getVolumes()()
    for (i in seq_along(corpus$assays)) {
        h$checkpointSetInput(driver, paste0(prefix, "assay", i, "_name"), names(corpus$assays)[[i]])
        id <- paste0(prefix, "assay", i, "_file")
        selection <- h$checkpointFileSelection(corpus$assays[[i]], roots)
        h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null",
            jsonlite::toJSON(id, auto_unbox = TRUE)), timeout = h$checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("Shiny.setInputValue(%s,%s,{priority:'event'});true",
            jsonlite::toJSON(id, auto_unbox = TRUE), jsonlite::toJSON(selection, auto_unbox = TRUE)))
        h$checkpointWaitForIdle(driver, timeout = h$checkpointWaitBudget())
    }
    h$checkpointSetInput(driver, paste0(prefix, "sanitize_names"), FALSE)
    kitClickAndWaitText(h, driver, paste0(prefix, "process_import"), ".shiny-notification", "Data imported successfully!")
    invisible(NULL)
}

#' KIT-068: one step of the lipidomics or metabolomics route through the module's own tab actions (the twin modules
#' share their input ids under `<route>_workflow-<tab>-`). Every action waits on the module's own completion text (a
#' notification or a status output) or an error notification, as `runKitQcAction` does for the proteomics tabs.
#' @noRd
runKitOmicsStep <- function(h, driver, project_parent, corpus, journey, step) {
    route <- kitRouteOf(journey)
    short <- switch(route, lipidomics = "lipid", metabolomics = "metab", stop("the omics steps drive lipidomics or metabolomics only, not ", route))
    tabs <- paste0(route, "_workflow-", route, "_tabs")
    prefix <- function(tab) paste0(route, "_workflow-", tab, "-")
    qc_tabs <- paste0(prefix("qc"), short, "_qc_tabs")
    notify <- function(input_id, text) kitClickAndWaitText(h, driver, input_id, ".shiny-notification", text)
    open_tab <- function(value) h$checkpointSetInput(driver, tabs, value)
    settings <- journey$normalisation
    # KIT-082: a display is observed by what it shows. `settled` ends a wait for something that will not come: the session has
    # been idle for 20 s, which a worker's display never relies on.
    settled <- paste0("(()=>{const busy=document.documentElement.classList.contains('shiny-busy');",
        "if(busy||!window.__kitSettled)window.__kitSettled=Date.now();return Date.now()-window.__kitSettled>20000;})()")
    wait_settled <- function(found) {
        driver$get_js("window.__kitSettled=0;true")
        h$checkpointWaitForJs(driver, sprintf("(%s)||(%s)", found, settled), timeout = h$checkpointWaitBudget())
        isTRUE(driver$get_js(found))
    }
    has_data_row <- function(selector) sprintf("Boolean(document.querySelector(%s))",
        jsonlite::toJSON(paste(selector, "tbody tr td:not(.dataTables_empty)"), auto_unbox = TRUE))
    # The share of an output's image that is not blank, or -1 with no loaded image: an empty plot is still an image.
    ink <- function(output) sprintf(paste0("(function(){const el=document.getElementById(%s);const img=el&&el.querySelector('img');",
        "if(!img||!img.complete||!img.naturalWidth)return -1;const c=document.createElement('canvas');c.width=img.naturalWidth;c.height=img.naturalHeight;",
        "const x=c.getContext('2d');x.drawImage(img,0,0);const d=x.getImageData(0,0,c.width,c.height).data;let n=0;",
        "for(let i=0;i<d.length;i+=4)if(d[i+3]>0&&(d[i]<245||d[i+1]<245||d[i+2]<245))n++;return n/(d.length/4);})()"),
        jsonlite::toJSON(output, auto_unbox = TRUE))
    notifications <- function() driver$get_js(paste0("(function(){const live=Array.from(document.querySelectorAll('.shiny-notification')).map(n=>n.textContent);",
        "const log=(window.__kitNotifications||[]).map(n=>n.text);return live.concat(log).join(' | ');})()"))
    switch(step$id,
        import = runKitOmicsImportActions(h, driver, project_parent, corpus, route),
        design = runKitDesignActions(h, driver, corpus$design, prefix = prefix("design-builder"), tab_input = tabs,
            tab_value = "design", save = function(id) notify(id, "Design matrix and contrasts saved successfully!")),
        intensity = {
            open_tab("qc")
            h$checkpointSetInput(driver, qc_tabs, "Intensity Filter")
            notify(paste0(prefix("qc"), "intensity-apply_filter"), "Intensity filter applied")
        },
        duplicates = {
            h$checkpointSetInput(driver, qc_tabs, "Duplicate Resolution")
            notify(paste0(prefix("qc"), "duplicates-detect_duplicates"), "Detection complete")
            notify(paste0(prefix("qc"), "duplicates-resolve_duplicates"), "Duplicates resolved")
        },
        qc_finalize = {
            # KIT-082: the internal-standard analysis, which commits nothing. Its start notice also says "internal standards",
            # and the summary shows a placeholder from the start, so only the finished count is the observation.
            h$checkpointSetInput(driver, qc_tabs, "Internal Standards")
            h$checkpointSetInput(driver, paste0(prefix("qc"), "itsd-is_pattern"), journey$is_pattern %||% "^IS_")
            notify(paste0(prefix("qc"), "itsd-analyze_is"), "Found ")
            count <- regmatches(notifications(), regexpr("Found [0-9]+ internal standards", notifications()))
            if (!length(count) || identical(count, "Found 0 internal standards")) stop("the internal-standard analysis found none: ", toString(count), call. = FALSE)
            cat("DISPLAY", count, "\n")
            h$checkpointSetInput(driver, qc_tabs, "Finalize QC")
            notify(paste0(prefix("qc"), "s4_finalize-finalize_qc"), "QC complete!")
        },
        normalisation = {
            # The one tab action: ITSD (when bound), log2, between-sample normalisation and RUV-III in the bound mode;
            # the RUV sweep settings are the journey's (the same input ids as the proteomics tab).
            open_tab("norm")
            norm <- prefix("norm")
            h$checkpointSetInput(driver, paste0(norm, "apply_itsd"), isTRUE(settings$apply_itsd))
            h$checkpointSetInput(driver, paste0(norm, "norm_method"), settings$norm_method %||% "cyclicloess")
            h$checkpointSetInput(driver, paste0(norm, "ruv_mode"), settings$ruv_mode %||% "automatic")
            for (name in names(settings$ruv)) h$checkpointSetInput(driver, paste0(norm, name), settings$ruv[[name]])
            # KIT-082: the ITSD candidate tables, built a page at a time: a data row in every assay's table before the run.
            if (isTRUE(settings$apply_itsd)) {
                drawn_tables <- sprintf("Array.from(document.querySelectorAll('[id^=\"%sitsd_table_\"]')).filter(el=>el.querySelector('tbody tr td:not(.dataTables_empty)')).length", norm)
                if (!wait_settled(sprintf("(%s)>=%d", drawn_tables, length(corpus$assays)))) {
                    stop("the ITSD candidate tables did not all draw: ", driver$get_js(drawn_tables), " of ", length(corpus$assays), call. = FALSE)
                }
            }
            notify(paste0(norm, "run_normalization"), "Normalization pipeline complete!")
        },
        correlation = {
            norm <- prefix("norm")
            h$checkpointSetInput(driver, paste0(norm, "norm_qc_tabs"), "Correlation Filtering")
            h$checkpointSetInput(driver, paste0(norm, "min_pearson_correlation_threshold"), settings$correlation_threshold)
            notify(paste0(norm, "apply_correlation_filter"), "Correlation filtering complete!")
        },
        session_export = notify(paste0(prefix("norm"), "export_session"), "Session data exported successfully!"),
        da = {
            open_tab("de")
            de <- prefix("de")
            notify(paste0(de, "load_filtered_session"), "Session loaded successfully!")
            h$checkpointSetInput(driver, paste0(de, "formula_string"), journey$da$formula)
            h$checkpointSetInput(driver, paste0(de, "da_q_val_thresh"), journey$da$q_value)
            h$checkpointSetInput(driver, paste0(de, "treat_lfc_cutoff"), journey$da$treat_lfc)
            h$checkpointClickAndWaitText(driver, paste0(de, "run_da_analysis"), paste0("#", de, "da_status"), "Analysis Complete")
            # KIT-082: the tab's interactive volcano. The assay opens on Combined, whose banner draws no plot, so select the first
            # assay and wait for the plot itself: a worker's page from KIT-082, the widget before it. A worker's display takes as long
            # as it takes, with the session idle, so this wait has its own ends: the plot, or a warning, error or output error.
            volcano <- jsonlite::toJSON(paste0(de, "volcano_glimma"), auto_unbox = TRUE)
            assay_input <- jsonlite::toJSON(paste0(de, "volcano_assay"), auto_unbox = TRUE)
            assay_script <- sprintf("(function(){const el=document.getElementById(%s);if(!el||!el.selectize)return '';return Object.keys(el.selectize.options).filter(k=>k!=='Combined')[0]||'';})()", assay_input)
            h$checkpointWaitForJs(driver, sprintf("%s!==''", assay_script), timeout = h$checkpointWaitBudget())
            assay <- driver$get_js(assay_script)
            h$checkpointSetInput(driver, paste0(de, "volcano_assay"), assay)
            h$checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);return el!==null&&",
                "(el.classList.contains('shiny-output-error')||el.querySelector('iframe,.html-widget,.alert-danger,.alert-warning')!==null);})()"), volcano),
                timeout = h$checkpointWaitBudget())
            failure <- driver$get_js(sprintf("(function(){const el=document.getElementById(%s);if(el.classList.contains('shiny-output-error'))return el.textContent||'output error';const a=el.querySelector('.alert-danger,.alert-warning');return a?a.textContent:'';})()", volcano))
            if (nzchar(failure)) stop("the DA tab's interactive volcano failed: ", failure, call. = FALSE)
            # KIT-082: the tab's other displays. A plot must have ink, a table a data row, the download data rows.
            output_error <- function(output) {
                error <- driver$get_js(sprintf("(function(){const el=document.getElementById(%s);return el&&el.classList.contains('shiny-output-error')?(el.textContent||'output error'):'';})()",
                    jsonlite::toJSON(output, auto_unbox = TRUE)))
                if (nzchar(error)) stop("the DA tab's ", output, " failed: ", error, call. = FALSE)
            }
            plotted <- function(name, least) {
                output <- paste0(de, name)
                wait_settled(sprintf("(%s)>=%s", ink(output), least))
                output_error(output)
                share <- driver$get_js(ink(output))
                cat(sprintf("DISPLAY %s ink %.4f\n", name, share))
                if (share < least) stop("the DA tab's ", name, " drew nothing: ink ", share, call. = FALSE)
            }
            h$checkpointSetInput(driver, paste0(de, "volcano_interactive"), FALSE)
            plotted("volcano_static", 0.005)
            h$checkpointSetInput(driver, paste0(de, "volcano_interactive"), TRUE)
            h$checkpointSetInput(driver, paste0(de, "da_results_tabs"), "Heatmap")
            h$checkpointSetInput(driver, paste0(de, "heatmap_assay"), assay)
            plotted("heatmap_plot", 0.05)
            h$checkpointSetInput(driver, paste0(de, "da_results_tabs"), "DA Results Table")
            table_found <- has_data_row(paste0("#", de, "da_results_table"))
            found <- wait_settled(table_found)
            output_error(paste0(de, "da_results_table"))
            if (!found) stop("the DA tab's results table has no data row", call. = FALSE)
            # The download, through the handler: its text arrives in the page, and only its counts are kept.
            link <- jsonlite::toJSON(paste0(de, "download_da_results"), auto_unbox = TRUE)
            href <- sprintf("(document.getElementById(%s).getAttribute('href')||'').length>1", link)
            if (!wait_settled(href)) stop("the DA tab's download has no link", call. = FALSE)
            driver$get_js(sprintf(paste0("window.__kitDownload=null;fetch(document.getElementById(%s).getAttribute('href')).then(r=>r.ok?r.text():Promise.reject(r.status))",
                ".then(t=>{const lines=t.split('\\n').filter(l=>l.length);window.__kitDownload={rows:lines.length-1,sampleColumns:lines.length?(lines[0].match(/intensity[.]/g)||[]).length:0};})",
                ".catch(e=>{window.__kitDownload={error:String(e)};});true"), link))
            h$checkpointWaitForJs(driver, "window.__kitDownload!==null", timeout = h$checkpointWaitBudget())
            download <- driver$get_js("window.__kitDownload")
            if (!is.null(download$error)) stop("the DA tab's download failed: ", download$error, call. = FALSE)
            cat(sprintf("DISPLAY download rows %d sample columns %d\n", download$rows, download$sampleColumns))
            if (download$rows < 1L || download$sampleColumns < 1L) stop("the DA tab's download has no data rows or no sample columns", call. = FALSE)
        },
        report = {
            # The summary tab: the publication copy (which writes the study parameters when they were not saved), then
            # the report render with the omic's shipped template (KIT-065's step on the proteomics routes).
            open_tab("summary")
            summary <- prefix("summary")
            h$checkpointSetInput(driver, paste0(summary, "experiment_label"), journey$experiment_label)
            h$checkpointSetInput(driver, paste0(summary, "description"), "KIT-068 lifecycle member")
            # KIT-082: the study parameters and the integration object (on a managed session, the worker's).
            notify(paste0(summary, "save_workflow_args"), "Study parameters saved")
            if (grepl("saved with warnings", notifications(), fixed = TRUE)) stop("the study parameters were saved with warnings", call. = FALSE)
            h$checkpointClickAndWaitText(driver, paste0(summary, "copy_to_publication"), paste0("#", summary, "copy_status"),
                "Files copied to publication directory")
            notify(paste0(summary, "generate_report"), "Report generated successfully!")
        },
        stop("unknown omics step ", step$id))
    invisible(NULL)
}

#' KIT-098: the transcriptomics import: the tile, the project setup, the import tab's count-table format, the count file through
#' the module's own shinyFiles button (selected as the pinned helpers select the DIA report), then the process button. The id
#' and name columns are the tab's proposal (`GeneID`, `GeneName`); the profile offers no sanitising control.
#' @param corpus `list(assays = list(counts = <path>), vendor = "trx_counts", design)` from the coordinator.
#' @noRd
runKitTrxImportActions <- function(h, driver, run_dir, corpus) {
    prefix <- "transcriptomics_workflow-import-"
    stopifnot(length(corpus$assays) == 1L)
    h$checkpointClickControl(driver, "tile-transcriptomics")
    h$checkpointClickControl(driver, "btn-start-analysis")
    driver$set_inputs(experiment_label = "kit_omics_lifecycle", project_dir = run_dir)
    h$checkpointClickControl(driver, "btn-confirm-setup")
    h$checkpointWaitForJs(driver, "document.querySelector('[data-testid=\"trx-tab-import\"]') !== null", timeout = h$checkpointWaitBudget())
    h$checkpointSetInput(driver, paste0(prefix, "vendor_format"), corpus$vendor %||% "trx_counts")
    id <- paste0(prefix, "assay1_file")
    selection <- h$checkpointFileSelection(corpus$assays[[1L]], shinyFiles::getVolumes()())
    h$checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", jsonlite::toJSON(id, auto_unbox = TRUE)),
        timeout = h$checkpointWaitBudget(control = TRUE))
    driver$get_js(sprintf("Shiny.setInputValue(%s,%s,{priority:'event'});true", jsonlite::toJSON(id, auto_unbox = TRUE),
        jsonlite::toJSON(selection, auto_unbox = TRUE)))
    h$checkpointWaitForIdle(driver, timeout = h$checkpointWaitBudget())
    kitClickAndWaitText(h, driver, paste0(prefix, "process_import"), ".shiny-notification", "Data imported successfully!")
    invisible(NULL)
}

#' KIT-098: one step of the transcriptomics route through KIT-097's tabs (the lipid modules under `trxOmicProfile()`, ids under
#' `transcriptomics_workflow-<tab>-`). Each action waits on the module's own completion text or an error notification. The
#' settings are the binding's (`journey$normalisation`, `journey$da`, `journey$ruviii_m_source`), never the shared journey's.
#' @noRd
runKitTrxStep <- function(h, driver, project_parent, corpus, journey, step) {
    tabs <- "transcriptomics_workflow-transcriptomics_tabs"
    prefix <- function(tab) paste0("transcriptomics_workflow-", tab, "-")
    notify <- function(input_id, text) kitClickAndWaitText(h, driver, input_id, ".shiny-notification", text)
    open_tab <- function(value) h$checkpointSetInput(driver, tabs, value)
    settings <- journey$normalisation
    switch(step$id,
        import = runKitTrxImportActions(h, driver, project_parent, corpus),
        design = runKitDesignActions(h, driver, corpus$design, prefix = prefix("design-builder"), tab_input = tabs,
            tab_value = "design", save = function(id) notify(id, "Design matrix and contrasts saved successfully!")),
        expression_filter = {
            open_tab("qc")
            notify(paste0(prefix("qc"), "filter-apply_filter"), "Expression filter applied")
        },
        qc_finalize = {
            h$checkpointSetInput(driver, paste0(prefix("qc"), "trx_qc_tabs"), "Finalize QC")
            notify(paste0(prefix("qc"), "s4_finalize-finalize_qc"), "QC complete!")
        },
        normalisation = {
            open_tab("norm")
            norm <- prefix("norm")
            h$checkpointSetInput(driver, paste0(norm, "ruv_mode"), settings$ruv_mode)
            for (name in c("auto_percentage_min", "auto_percentage_max", "auto_percentage_step", "k_penalty_weight"))
                h$checkpointSetInput(driver, paste0(norm, name), settings[[name]])
            notify(paste0(norm, "run_normalization"), "Normalization complete!")
        },
        da = {
            open_tab("de")
            de <- prefix("de")
            h$checkpointSetInput(driver, paste0(de, "da_q_val_thresh"), journey$da$q_value)
            h$checkpointSetInput(driver, paste0(de, "treat_lfc_cutoff"), journey$da$treat_lfc)
            h$checkpointClickAndWaitText(driver, paste0(de, "run_da_analysis"), paste0("#", de, "da_status"), "Analysis Complete")
        },
        ruviii_export = {
            open_tab("summary")
            summary <- prefix("summary")
            h$checkpointSetInput(driver, paste0(summary, "ruviii_m_source"), journey$ruviii_m_source)
            h$checkpointClickAndWaitText(driver, paste0(summary, "run_ruviii_export"), paste0("#", summary, "ruviii_status"),
                paste("M =", journey$ruviii_m_source))
        },
        report = {
            summary <- prefix("summary")
            h$checkpointSetInput(driver, paste0(summary, "experiment_label"), journey$experiment_label)
            h$checkpointSetInput(driver, paste0(summary, "description"), "KIT-098 transcriptomics journey")
            notify(paste0(summary, "save_workflow_args"), "Study parameters saved")
            logged <- driver$get_js("(window.__kitNotifications||[]).map(n=>n.text).join(' | ')")
            if (grepl("saved with warnings", paste(logged, collapse = " "), fixed = TRUE)) stop("the study parameters were saved with warnings", call. = FALSE)
            h$checkpointClickAndWaitText(driver, paste0(summary, "copy_to_publication"), paste0("#", summary, "copy_status"),
                "Files copied to publication directory")
            notify(paste0(summary, "generate_report"), "Report generated successfully!")
        },
        session_export = notify(paste0(prefix("summary"), "export_session_state"), "Session state exported to:"),
        stop("unknown transcriptomics step ", step$id))
    invisible(NULL)
}
