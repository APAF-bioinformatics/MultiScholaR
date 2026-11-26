# KIT-277 (VAL-331, VAL-345, VAL-353): a design saved again with only the model formula, (proteomics import) only the
# organism or taxonomy id, or only the FASTA file changed asks first and, on Save, is carried through to every consumer.
# Each scenario drives one lane's module server from import under shiny::testServer, through its namespaced inputs only
# (the GUI-entry driver, helper-gui-entry.R), as test-kit247-design-resave.R does (its helpers are copied here, file-local,
# never sourced); one scenario's observations serve several rows. Rows are tagged by stage: S1 (formula), S2 (organism),
# S3 (FASTA); MULTISCHOLAR_KIT277_STAGE=1 runs the S1 rows and G2 only (the stage-1 gate). A row FAILS, it does not
# ERROR, on the old code: every step after the Save runs whether or not an earlier expectation held, and a read that the
# old code cannot answer is recorded as an error text, not raised. Texts are drafts (VAL-41).

kit277StageOnly <- function() identical(Sys.getenv("MULTISCHOLAR_KIT277_STAGE"), "1")

#' Diagnostic log (MULTISCHOLAR_KIT277_DIAG set; O4's scenario only; no expectation reads it): the store's state at a point.
kit277Diag <- function(label, api, probe) {
    file <- Sys.getenv("MULTISCHOLAR_KIT277_DIAG_FILE")
    if (!nzchar(Sys.getenv("MULTISCHOLAR_KIT277_DIAG")) || !nzchar(file)) return(invisible(NULL))
    emit <- function(name, value) cat("KIT277DIAG", label, name, "=", paste(utils::capture.output(utils::str(value, give.attr = FALSE)), collapse = " "), "\n", file = file, append = TRUE)
    manager <- probe$manager()
    context <- manager$.context
    snap <- context$snapshot
    emit("manager_current_state", kit277Try(manager$current_state))
    emit("snapshot_current_state", kit277Try(snap$current_state))
    emit("history_entities", kit277Try(vapply(as.character(unlist(snap$state_history)), function(s) snap$states[[s]]$entity_id %||% NA_character_, "")))
    emit("context_head_entity", kit277Try(context$head$entity_id))
    emit("store_identity", kit277Try(context$store$identity))
    emit("da_results_entity", kit277Try(snap$da_results$entity_id %||% "<none>"))
    emit("enrichment_marker", kit277Try(!is.null(snap$enrichment_results)))
    emit("da_manifest", kit277Try(withWorkflowStateScope(manager, function() {
        if (!isSavedWorkflowReference(snap$da_results)) return("<no DA>")
        m <- resolveWorkflowStateValue(context, snap$da_results)
        list(source = m$source$entity_id, owner = m$source$owner, parameters = m$parameters, contrasts = m$contrasts$friendly_names,
            result_names = names(m$results))
    })))
    emit("da_status", kit277Try(workflowDaResultStatus(manager)))
    emit("status_tests", kit277Try(withWorkflowStateScope(manager, function() {
        if (!isSavedWorkflowReference(snap$da_results)) return("<no DA reference>")
        source <- resolveWorkflowStateValue(context, snap$da_results)$source
        history <- as.character(unlist(snap$state_history))
        table <- workflowDaTableState(manager)
        fields <- workflowDaTableFields(table$record$data$class)
        list(owner_same = identical(source$owner, context$store$identity),
            member = any(vapply(history, function(s) identical(snap$states[[s]]$entity_id, source$entity_id), logical(1))),
            differ = workflowDaTablesDiffer(manager, source, table),
            source_components = vapply(readComponentMetadata(manager, source)$data$components[fields], function(c) c$entity_id %||% NA_character_, ""),
            head_components = vapply(table$record$data$components[fields], function(c) c$entity_id %||% NA_character_, ""))
    })))
    emit("daAnalysisUI_head", kit277Try(withWorkflowStateScope(manager, function() readWorkflowStateMetadata(manager)$data$slots$args$daAnalysisUI)))
    emit("candidate_files", kit277Try(sub(probe$paths$base_dir, "", workflowDaCandidateFiles(probe$paths, "proteomics"), fixed = TRUE)))
    emit("made_files", kit277Try(withWorkflowStateScope(manager, function() {
        if (!isSavedWorkflowReference(snap$da_results)) return("<no DA>")
        sub(probe$paths$base_dir, "", workflowDaMadeFiles(manager, probe$paths, "proteomics"), fixed = TRUE)
    })))
    emit("unmade_files", kit277Try(withWorkflowStateScope(manager, function() {
        if (!isSavedWorkflowReference(snap$da_results)) return("<no DA>")
        sub(probe$paths$base_dir, "", workflowDaFilesUnmade(manager, probe$paths, "proteomics"), fixed = TRUE)
    })))
    emit("publication_marker", kit277Try(tryCatch({ assertWorkflowDaPublicationComplete(probe$paths); "complete" }, error = function(e) paste(class(e)[[1L]], conditionMessage(e)))))
    emit("notes", kit277Try(probe$notes()))
    emit("copy_gate", kit277Try(tryCatch({
        shiny::isolate(assertWorkflowSummaryDaCurrent(api$workflow_data, "proteomics", probe$paths, files = TRUE)); "allowed"
    }, error = function(e) paste(paste(class(e), collapse = "/"), conditionMessage(e)))))
    invisible(NULL)
}
kit277Skip <- function(stage) {
    if (kit277StageOnly() && !stage %in% c("S1", "G2")) testthat::skip("stage-1 gate: the S1 rows and G2 only")
    # Diagnostic run (O4's miss): MULTISCHOLAR_KIT277_STAGE=2 runs the S2 rows only; nothing is judged by it.
    if (identical(Sys.getenv("MULTISCHOLAR_KIT277_STAGE"), "2") && !identical(stage, "S2")) testthat::skip("diagnostic: S2 rows only")
}

kit277Cache <- new.env(parent = emptyenv())
kit277Once <- function(key, make) {
    if (!exists(key, envir = kit277Cache, inherits = FALSE)) assign(key, make(), envir = kit277Cache)
    get(key, envir = kit277Cache, inherits = FALSE)
}
kit277Try <- function(expr) tryCatch(expr, error = function(e) paste("<error>", conditionMessage(e)))

kit277Click <- function(api, id) {
    current <- suppressWarnings(as.integer(api$input(id)))
    following <- if (length(current) == 1L && !is.na(current)) current + 1L else 1L
    api$set(stats::setNames(list(following), id))
}

kit277DesignActions <- function(api, design, contrasts, prefix, tech_rep, formula = "~ 0 + group") {
    builder <- function(name) paste0(prefix, name)
    groups <- unique(design$group)
    for (i in seq_along(groups)) {
        api$set(stats::setNames(list(groups[[i]]), builder("new_factor")))
        kit277Click(api, builder("add_factor"))
    }
    for (i in seq_along(groups)) {
        api$set(stats::setNames(list(design$sample[design$group == groups[[i]]], groups[[i]], "", "", 1L),
            builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
        kit277Click(api, builder("assign_metadata"))
    }
    if (tech_rep && !is.null(design$tech_rep_group)) {
        for (set_name in unique(design$tech_rep_group)) {
            samples <- design$sample[design$tech_rep_group == set_name]
            if (length(samples) < 2L) next
            api$set(stats::setNames(list(samples, "manual", as.integer(design$replicates[match(samples[[1L]], design$sample)])),
                builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
            kit277Click(api, builder("assign_tech_reps"))
        }
    }
    for (pair in contrasts) {
        api$set(stats::setNames(list(pair[[1L]], pair[[2L]], formula), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
        kit277Click(api, builder("add_contrast"))
    }
    kit277Click(api, builder("save_results"))
}

kit277Contrast <- function(lane_id) {
    switch(lane_id, trx_screen = c("Sera", "RPMI"), c("KO", "WT"))
}

# Save Design again with the design as it was: kind "same", "formula" (the formula only), "contrast" (a second comparison)
# or "formula_contrast" (both).
kit277Resave <- function(api, lane, kind, formula = "~ 0 + group + replicates") {
    prot <- identical(lane$omic_type, "proteomics")
    prefix <- if (prot) "design_matrix-builder-" else "design-builder-"
    design <- guiEntryReadDesign(api)
    if (prot) design$sample <- e2e_match_available_values(design$sample, api$choices("design_matrix-builder-selected_runs"))
    pair <- kit277Contrast(lane$lane_id)
    contrasts <- if (kind %in% c("contrast", "formula_contrast")) list(pair, rev(pair)) else list(pair)
    chosen <- if (kind %in% c("formula", "formula_contrast")) formula else "~ 0 + group"
    api$set(stats::setNames(list("design"), api$tabs))
    kit277DesignActions(api, design, contrasts, prefix, tech_rep = !prot, formula = chosen)
}

kit277Confirm <- function(api, lane) {
    kit277Click(api, if (identical(lane$omic_type, "proteomics")) "design_matrix-confirm_design_resave" else "design-confirm_design_resave")
}

#' One lane up to the lane step `until`, then `actions(api, probe)`. As test-kit247-design-resave.R's `kit247Lane`, with the
#' text of every modal recorded (`probe$modal_text()`).
kit277Lane <- function(lane_id, until, actions) {
    testthat::skip_if_not(workflowOwnerCapable(), "the session owner is not available (qs2/duckdb)")
    lane <- guiEntryLane(lane_id)
    omic <- lane$omic_type
    spec <- guiEntryServerSpec(omic)
    standard <- guiEntrySteps(lane)
    ids <- vapply(standard, `[[`, "", "id")
    prefix_steps <- standard[seq_len(match(until, ids))]
    root <- normalizePath(withr::local_tempdir(pattern = paste0("k277-", lane_id, "-")))
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
    withr::local_options(warn = 1)
    withr::local_seed(20261107L)
    project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", lane_id), omic_types = omic,
        label = NULL, force = TRUE, reuse_existing = TRUE))
    paths <- project_dirs[[omic]]
    old_dirs <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
    assign("project_dirs", project_dirs, envir = globalenv())
    withr::defer(if (is.null(old_dirs)) rm("project_dirs", envir = globalenv()) else assign("project_dirs", old_dirs, envir = globalenv()))
    if (identical(omic, "proteomics")) {
        guiEntryInstallUniprotCache(paths)
        install_enrichment_doubles()
    }
    log <- new.env(parent = emptyenv())
    log$updates <- list(); log$notes <- list(); log$modals <- integer(); log$choices <- list(); log$modal_text <- character()
    guiEntryInstallBrowserDouble(log)
    recorder <- function(ui, session = shiny::getDefaultReactiveDomain()) {
        log$modals <- c(log$modals, 1L)
        log$modal_text <- c(log$modal_text, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }
    testthat::local_mocked_bindings(showModal = recorder, .package = "shiny")
    if ("showModal" %in% unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE)) {
        testthat::local_mocked_bindings(showModal = recorder, .package = "MultiScholaR")
    }
    inputs <- list(dir = guiEntryFixtureDir(lane), scratch = file.path(root, "inputs"), root = root, lane = lane, paths = paths)
    dir.create(inputs$scratch, showWarnings = FALSE)
    holder <- new.env(parent = emptyenv())
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(spec$server,
        args = list(id = spec$id, project_dirs = project_dirs, omic_type = omic, experiment_label = lane_id), {
        rootSession <- session$rootScope()
        prefix <- paste0(spec$id, "-")
        apply_updates <- function() {
            for (i in seq_len(50L)) {
                if (!length(log$updates)) return(invisible())
                pending <- log$updates; log$updates <- list()
                values <- stats::setNames(lapply(pending, `[[`, "value"), vapply(pending, `[[`, "", "id"))
                values <- values[!duplicated(names(values), fromLast = TRUE)]
                do.call(rootSession$setInputs, values)
            }
        }
        api <- list(
            set = function(values) {
                names(values) <- paste0(prefix, names(values))
                do.call(rootSession$setInputs, values)
                apply_updates()
            },
            defaults = function(ui_fn, id) {
                values <- omicsUiInputDefaults(ui_fn(id))
                if (length(values)) { names(values) <- paste0(prefix, names(values)); do.call(rootSession$setInputs, values); apply_updates() }
            },
            upload = function(id, path) {
                copy <- file.path(inputs$scratch, basename(path))
                file.copy(path, copy, overwrite = TRUE)
                value <- data.frame(name = basename(path), size = file.size(copy), type = "", datapath = copy, stringsAsFactors = FALSE)
                do.call(rootSession$setInputs, stats::setNames(list(value), paste0(prefix, id)))
                apply_updates()
            },
            input = function(id) shiny::isolate(rootSession$input[[paste0(prefix, id)]]),
            choices = function(id) log$choices[[paste0(prefix, id)]] %||% character(),
            modals = function() log$modals,
            workflow_data = workflow_data, session = session, inputs = inputs, tabs = spec$tabs)
        for (step in prefix_steps) {
            if (!is.null(step$tab)) api$set(stats::setNames(list(step$tab), spec$tabs))
            step$run(api)
            session$flushReact()
        }
        log$notes <- list()
        manager <- function() shiny::isolate(workflow_data$state_manager)
        probe <- list(
            lane = lane, paths = paths,
            head = function() manager()$current_state,
            entity = function(name = NULL) getWorkflowStateReference(manager(), name)$entity_id,
            history = function() guiEntryHistory(manager()),
            status = function() lapply(shiny::isolate(workflow_data$tab_status), as.character),
            notes = function() vapply(log$notes, `[[`, "", "message"),
            clear_notes = function() log$notes <- list(),
            modals = function() length(log$modals),
            modal_text = function() log$modal_text,
            bus = function(name) shiny::isolate(workflow_data[[name]]),
            steps = stats::setNames(standard, ids), tabs = spec$tabs,
            manager = manager, flush = function() session$flushReact())
        holder$result <- actions(api, probe)
        shiny::isolate(closeDIAStateAdapter(workflow_data$state_manager))
    }))))
    holder$result
}

# ---------------------------------------------------------------------------------------------------------------------
# What the rows read.
# ---------------------------------------------------------------------------------------------------------------------

kit277StatusKey <- function(lane) if (identical(lane$omic_type, "proteomics")) "differential_expression" else "differential_analysis"
kit277DaTab <- function(lane) if (identical(lane$omic_type, "proteomics")) "da" else "de"
kit277BoxId <- function(lane) if (identical(lane$omic_type, "proteomics")) "differential_expression-formula_string" else "de-formula_string"

kit277DaRef <- function(manager) manager$.context$snapshot$da_results

kit277BusFormula <- function(probe) {
    config <- probe$bus("config_list")
    kit277Try(config[["deAnalysisParameters"]][["formula_string"]])
}

#' The formula the store's DA ran with (its manifest's parameters), or "<no DA>".
kit277DaFormula <- function(manager) {
    kit277Try(withWorkflowStateScope(manager, function() {
        reference <- manager$.context$snapshot$da_results
        if (!isSavedWorkflowReference(reference)) return("<no DA>")
        resolveWorkflowStateValue(manager$.context, reference)$parameters$formula
    }))
}

#' The head's recorded formulas (proteomics): the design settings' and the DA arguments'.
kit277HeadFormulas <- function(manager) {
    kit277Try(withWorkflowStateScope(manager, function() {
        args <- readWorkflowStateMetadata(manager)$data$slots$args
        list(design = args$deAnalysisParameters$formula_string, da_ui = args$daAnalysisUI$formula_string)
    }))
}

#' Whether the Copy and Report gate (`files = TRUE`) or the Summary's (`files = FALSE`) lets the DA through.
kit277Gate <- function(api, probe, files) {
    tryCatch({
        shiny::isolate(assertWorkflowSummaryDaCurrent(api$workflow_data, probe$lane$omic_type, probe$paths, files = files))
        "allowed"
    }, workflow_summary_da_stale = function(e) "refused", error = function(e) paste("<error>", conditionMessage(e)))
}

kit277RecordHasDa <- function(probe) {
    kit277Try("differential_abundance" %in% names(buildWorkflowAnalysisRecord(probe$manager())$analyses))
}

#' Whether a session export written now would list the DA folder (proteomics: `collectDIANativeGuiFiles`).
kit277ExportHasDa <- function(probe) {
    kit277Try(withWorkflowStateScope(probe$manager(), function() {
        specs <- collectDIANativeGuiFiles(probe$manager(), list())
        any(startsWith(vapply(specs, function(spec) spec$destination, ""), "da/"))
    }))
}

kit277RunDa <- function(api, lane) {
    if (identical(lane$omic_type, "proteomics")) {
        api$set(stats::setNames(list("da"), api$tabs))
        api$set(list(`differential_expression-da_q_val_thresh` = 0.05, `differential_expression-treat_lfc_cutoff` = 0))
        modals <- length(api$modals())
        kit277Click(api, "differential_expression-run_da_analysis")
        if (length(api$modals()) > modals) kit277Click(api, "differential_expression-acknowledge_qvalue_warning")
    } else {
        api$set(stats::setNames(list("de"), api$tabs))
        api$set(list(`de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 0))
        kit277Click(api, "de-run_da_analysis")
    }
}

kit277Box <- function(api, lane) {
    api$set(stats::setNames(list(kit277DaTab(lane)), api$tabs))
    api$input(kit277BoxId(lane))
}

#' The Summary's Save Workflow Arguments then Copy, through the tab's own inputs; returns the notifications Copy raised.
kit277CopyClick <- function(api, probe) {
    lane <- probe$lane
    prefix <- if (identical(lane$omic_type, "proteomics")) "session_summary" else "summary"
    api$set(stats::setNames(list(prefix), api$tabs))
    api$set(stats::setNames(list(lane$lane_id, "KIT-277 copy check"), paste0(prefix, c("-experiment_label", "-description"))))
    kit277Click(api, paste0(prefix, "-save_workflow_args"))
    probe$clear_notes()
    kit277Click(api, paste0(prefix, "-copy_to_publication"))
    probe$notes()
}

#' The head rule the Summary's Save (R/utils_workflow_summary.R:297), Report (:524) and Export Session State (:622) all call.
kit277HeadRule <- function(api, probe) {
    tryCatch({
        shiny::isolate(assertWorkflowSummaryHeadFinal(api$workflow_data, probe$lane$omic_type))
        "allowed"
    }, omic_norm_correlation_input_stale = function(e) paste("refused:", conditionMessage(e)),
    error = function(e) paste("<error>", conditionMessage(e)))
}

kit277State <- function(api, probe) {
    lane <- probe$lane
    list(history = probe$history(), head = probe$head(), entity = probe$entity(), status = probe$status(),
        da = !is.null(kit277DaRef(probe$manager())), formula = kit277BusFormula(probe), box = kit277Box(api, lane),
        copy = kit277Gate(api, probe, TRUE), summary = kit277Gate(api, probe, FALSE))
}

#' The lane's exported design folder (the project's source_dir after its design save), copied, with `edit` applied to its
#' design_matrix.tab (as test-kit247-design-resave.R's `kit247ExportedDesign`).
kit277ExportedDesign <- function(probe, edit) {
    source <- probe$paths$source_dir
    target <- file.path(tempfile("k277-export-"))
    dir.create(target)
    files <- list.files(source, pattern = "^(design_matrix|data_cln.*|contrast_strings|assay_manifest|column_mapping|config)\\.(tab|txt|json|ini)$")
    file.copy(file.path(source, files), target)
    design <- utils::read.delim(file.path(target, "design_matrix.tab"), stringsAsFactors = FALSE, check.names = FALSE)
    design <- edit(design)
    utils::write.table(design, file.path(target, "design_matrix.tab"), sep = "\t", row.names = FALSE, quote = FALSE)
    target
}

kit277Import <- function(api, lane, folder, taxon = 9606L, organism = "Homo sapiens") {
    api$set(stats::setNames(list("design"), api$tabs))
    if (identical(lane$omic_type, "proteomics")) {
        api$set(list(`design_matrix-import_dir_text` = folder, `design_matrix-import_taxon_id` = taxon,
            `design_matrix-import_organism_name` = organism))
        kit277Click(api, "design_matrix-confirm_import")
    } else {
        api$set(list(`design-import_dir` = folder))
        kit277Click(api, "design-confirm_import")
    }
}

# ---------------------------------------------------------------------------------------------------------------------
# Stage 1: the formula.
# ---------------------------------------------------------------------------------------------------------------------

kit277NewFormula <- "~ 0 + group + replicates"

kit277FormulaScenario <- function(lane_id) {
    kit277Once(paste0("formula-", lane_id), function() kit277Lane(lane_id, "da", function(api, probe) {
        lane <- probe$lane
        prot <- identical(lane$omic_type, "proteomics")
        out <- list(key = kit277StatusKey(lane))
        out$before <- kit277State(api, probe)
        out$results_before <- if (prot) !is.null(probe$bus("da_analysis_results_list"))
        probe$clear_notes()
        modals <- probe$modals()
        kit277Resave(api, lane, "formula")
        out$asked <- probe$modals() > modals
        out$question <- utils::tail(probe$modal_text(), 1L)
        out$cancelled <- kit277State(api, probe)
        probe$clear_notes()
        kit277Resave(api, lane, "formula")
        if (out$asked) kit277Confirm(api, lane)
        probe$flush()
        out$notes <- probe$notes()
        out$after <- kit277State(api, probe)
        out$head_formulas <- if (prot) kit277HeadFormulas(probe$manager())
        out$results_after <- if (prot) !is.null(probe$bus("da_analysis_results_list"))
        out$design_formula <- kit277Try(omicDaDesignFormula(shiny::isolate(api$workflow_data)))
        out$record_da <- kit277RecordHasDa(probe)
        out$export_da <- if (prot) kit277ExportHasDa(probe)
        out$copy_notes_after <- kit277CopyClick(api, probe)
        probe$clear_notes()
        modals <- probe$modals()
        kit277Resave(api, lane, "formula")
        out$repeat_modals <- probe$modals() - modals
        out$repeat_notes <- probe$notes()
        kit277RunDa(api, lane)
        probe$flush()
        out$rerun_formula <- kit277DaFormula(probe$manager())
        out$rerun_copy <- kit277Gate(api, probe, TRUE)
        out$rerun_copy_notes <- kit277CopyClick(api, probe)
        out$rerun_design_formula <- kit277Try(omicDaDesignFormula(shiny::isolate(api$workflow_data)))
        out
    }))
}

for (lane_id in c("prot_lfq", "lipid_canonical", "metab_lc")) {
    test_that(paste("KIT-277 F1/F2 (S1):", lane_id, "a formula-only re-save asks, and Save puts the formula on the DA tab and sets the DA back"), {
        kit277Skip("S1")
        out <- kit277FormulaScenario(lane_id)
        expect_identical(out$before$status[[out$key]], "complete", info = lane_id)
        expect_identical(out$before$copy, "allowed", info = lane_id)
        expect_true(out$asked, info = lane_id)
        expect_match(paste(out$question, collapse = " "), "~ 0 + group + replicates", fixed = TRUE, info = lane_id)
        expect_match(paste(out$question, collapse = " "), "~ 0 + group", fixed = TRUE, info = lane_id)
        # Cancel (the question left unanswered) changes nothing.
        expect_identical(out$cancelled[c("history", "head", "entity", "status", "da", "formula", "box")],
            out$before[c("history", "head", "entity", "status", "da", "formula", "box")], info = lane_id)
        # Save: the formula on the DA tab, the DA not current for every consumer.
        expect_identical(out$after$box, kit277NewFormula, info = lane_id)
        expect_identical(out$after$status[[out$key]], "pending", info = lane_id)
        expect_false(out$after$da, info = lane_id)
        expect_identical(out$after$copy, "refused", info = lane_id)
        expect_identical(out$after$summary, "allowed", info = lane_id)
        expect_false(out$record_da, info = lane_id)
        # The Copy press itself: refused after the Save, working after the DA ran again.
        expect_true(any(grepl("Copy error", out$copy_notes_after, fixed = TRUE)), info = lane_id)
        if (identical(lane_id, "prot_lfq")) {
            expect_identical(out$head_formulas$design, kit277NewFormula)
            expect_identical(out$head_formulas$da_ui, kit277NewFormula)
            expect_true(out$results_before)
            expect_false(out$results_after)
            expect_false(out$export_da)
        } else {
            expect_identical(out$design_formula, kit277NewFormula, info = lane_id)
        }
        # The next run fits the new formula; Copy works again; the head's marker ends with the run.
        expect_identical(out$rerun_formula, kit277NewFormula, info = lane_id)
        expect_identical(out$rerun_copy, "allowed", info = lane_id)
        expect_false(any(grepl("Copy error", out$rerun_copy_notes, fixed = TRUE)), info = lane_id)
        expect_gt(length(out$rerun_copy_notes), 0L)
        if (!identical(lane_id, "prot_lfq")) expect_null(out$rerun_design_formula, info = lane_id)
    })

    test_that(paste("KIT-277 F5, F6 (S1):", lane_id, "the same formula saved again asks nothing; the notices say what was done"), {
        kit277Skip("S1")
        out <- kit277FormulaScenario(lane_id)
        # F5 holds only after the first Save applied (F1's values), then a repeat is quiet.
        expect_true(out$asked, info = lane_id)
        expect_identical(out$after$box, kit277NewFormula, info = lane_id)
        expect_identical(out$repeat_modals, 0L, info = lane_id)
        expect_true(any(grepl("Design saved; nothing changed.", out$repeat_notes, fixed = TRUE)), info = lane_id)
        # F6: after the applied change the notice is the formula notice, never "nothing changed".
        expect_true(any(grepl("Formula saved. Quality Control and normalisation are kept; run the differential analysis again.",
            out$notes, fixed = TRUE)), info = lane_id)
        expect_false(any(grepl("nothing changed", out$notes, fixed = TRUE)), info = lane_id)
    })
}

test_that("KIT-277 F2b (S1): the formula a lipidomics or metabolomics Save placed is on the head until the next DA run", {
    kit277Skip("S1")
    for (lane_id in c("lipid_canonical", "metab_lc")) {
        out <- kit277FormulaScenario(lane_id)
        expect_identical(out$design_formula, kit277NewFormula, info = lane_id)
        expect_identical(out$rerun_formula, kit277NewFormula, info = lane_id)
        expect_null(out$rerun_design_formula, info = lane_id)
    }
})

kit277TrxScenario <- function() {
    kit277Once("formula-trx", function() kit277Lane("trx_screen", "da", function(api, probe) {
        lane <- probe$lane
        out <- list(before = kit277State(api, probe), before_da = kit277DaFormula(probe$manager()))
        probe$clear_notes()
        modals <- probe$modals()
        kit277Resave(api, lane, "formula")
        probe$flush()
        out$modals <- probe$modals() - modals
        out$notes <- probe$notes()
        out$after <- kit277State(api, probe)
        kit277RunDa(api, lane)
        probe$flush()
        out$rerun_formula <- kit277DaFormula(probe$manager())
        out
    }))
}

test_that("KIT-277 F3 (S1): transcriptomics asks nothing and the notice says the formula does not change its analysis", {
    kit277Skip("S1")
    out <- kit277TrxScenario()
    expect_identical(out$modals, 0L)
    expect_true(any(grepl(paste("Formula saved. The transcriptomics differential analysis models the groups and the RUV",
        "factors from the design, so this formula does not change it."), out$notes, fixed = TRUE)))
    expect_identical(out$after$status, out$before$status)
    expect_identical(out$after$copy, "allowed")
    # The guard: the DA still fits what the design gives, not the builder's formula.
    expect_false(grepl("replicates", out$rerun_formula, fixed = TRUE))
    expect_identical(out$rerun_formula, out$before_da)
})

kit277ImportFormulaScenario <- function(lane_id) {
    kit277Once(paste0("import-formula-", lane_id), function() kit277Lane(lane_id, "da", function(api, probe) {
        lane <- probe$lane
        prot <- identical(lane$omic_type, "proteomics")
        folder <- kit277ExportedDesign(probe, identity)
        # The imported configuration: the folder's config.ini (lipidomics), the project's (proteomics).
        configFile <- if (prot) file.path(probe$paths$source_dir, "config.ini") else file.path(folder, "config.ini")
        config <- ini::read.ini(configFile)
        config$deAnalysisParameters$formula_string <- kit277NewFormula
        ini::write.ini(config, configFile)
        out <- list(key = kit277StatusKey(lane), before = kit277State(api, probe))
        modals <- probe$modals()
        probe$clear_notes()
        kit277Import(api, lane, folder)
        out$asked <- probe$modals() > modals
        out$question <- utils::tail(probe$modal_text(), 1L)
        out$cancelled <- kit277State(api, probe)
        kit277Import(api, lane, folder)
        if (out$asked) kit277Click(api, if (prot) "design_matrix-confirm_design_resave" else "design-confirm_design_resave")
        probe$flush()
        out$after <- kit277State(api, probe)
        out$head_formulas <- if (prot) kit277HeadFormulas(probe$manager())
        out$notes <- probe$notes()
        out
    }))
}

for (lane_id in c("prot_lfq", "lipid_canonical")) {
    test_that(paste("KIT-277 F4 (S1):", lane_id, "Import Existing Design with only the formula changed asks, Cancel leaves the configuration, Save carries it through"), {
        kit277Skip("S1")
        out <- kit277ImportFormulaScenario(lane_id)
        expect_true(out$asked, info = lane_id)
        expect_match(paste(out$question, collapse = " "), "~ 0 + group + replicates", fixed = TRUE, info = lane_id)
        expect_false(grepl("the organism", paste(out$question, collapse = " "), fixed = TRUE), info = lane_id)
        expect_identical(out$cancelled[c("history", "head", "entity", "status", "da", "formula", "box")],
            out$before[c("history", "head", "entity", "status", "da", "formula", "box")], info = lane_id)
        expect_identical(out$after$box, kit277NewFormula, info = lane_id)
        expect_identical(out$after$status[[out$key]], "pending", info = lane_id)
        expect_identical(out$after$copy, "refused", info = lane_id)
        if (identical(lane_id, "prot_lfq")) expect_identical(out$head_formulas$design, kit277NewFormula)
    })
}

test_that("KIT-277 F7 (S1): the formula and the comparisons changed in one save ask, and Save applies both", {
    kit277Skip("S1")
    out <- kit277Lane("prot_lfq", "da", function(api, probe) {
        lane <- probe$lane
        modals <- probe$modals()
        probe$clear_notes()
        kit277Resave(api, lane, "formula_contrast")
        asked <- probe$modals() > modals
        if (asked) kit277Confirm(api, lane)
        probe$flush()
        list(asked = asked, box = kit277Box(api, lane), contrasts = nrow(probe$bus("contrasts_tbl")),
            status = probe$status(), copy = kit277Gate(api, probe, TRUE))
    })
    expect_true(out$asked)
    expect_identical(out$box, kit277NewFormula)
    expect_identical(out$contrasts, 2L)
    expect_identical(out$status$differential_expression, "pending")
    expect_identical(out$copy, "refused")
})

test_that("KIT-277 G2 (guard): a comparisons-only change is KIT-247's: 'Comparisons saved', no question, the DA tab pending", {
    kit277Skip("G2")
    out <- kit277Lane("prot_lfq", "da", function(api, probe) {
        modals <- probe$modals()
        probe$clear_notes()
        kit277Resave(api, probe$lane, "contrast")
        list(modals = probe$modals() - modals, notes = probe$notes(), status = probe$status())
    })
    expect_identical(out$modals, 0L)
    expect_true(any(grepl("Comparisons saved", out$notes, fixed = TRUE)))
    expect_identical(out$status$differential_expression, "pending")
})

# ---------------------------------------------------------------------------------------------------------------------
# Stage 2: the organism.
# ---------------------------------------------------------------------------------------------------------------------

kit277IntegrationDir <- function(paths) paths$integration_dir %||% file.path(paths$base_dir, "integration")

kit277Enrich <- function(api, probe, taxid) {
    lane <- probe$lane
    api$set(stats::setNames(list("enrichment"), api$tabs))
    api$defaults(mod_prot_enrich_ui, "enrichment_analysis")
    api$set(list(`enrichment_analysis-organism_taxid` = taxid, `enrichment_analysis-enrichment_method_tabs` = "gprofiler2"))
    api$set(list(`enrichment_analysis-selected_contrast` = lane$expected_contrasts[[1L]], `enrichment_analysis-up_cutoff` = 0,
        `enrichment_analysis-down_cutoff` = 0, `enrichment_analysis-q_cutoff` = 0.05))
    kit277Click(api, "enrichment_analysis-run_enrichment_analysis")
}

kit277Organisms <- function(table) kit277Try(paste(sort(unique(as.character(table$Organism))), collapse = "|"))

kit277OrganismScenario <- function() {
    kit277Once("organism", function() kit277Lane("prot_dia", "enrichment", function(api, probe) {
        lane <- probe$lane
        paths <- probe$paths
        # The earlier annotation in every place a test can plant it, and an integration input written before the change.
        integration <- kit277IntegrationDir(paths)
        dir.create(file.path(integration, "proteomics"), recursive = TRUE, showWarnings = FALSE)
        writeLines("planted", file.path(integration, "proteomics", "manifest.json"))
        writeLines("planted", file.path(integration, "proteomics_planted_final_s4.RDS"))
        saveRDS(data.frame(Entry = "P00001", Organism = "Planted"), file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
        cache <- file.path(paths$results_dir, "cache", "uniprot_annotations")
        dir.create(cache, recursive = TRUE, showWarnings = FALSE)
        writeLines("planted", file.path(cache, "uniprot_results.tsv"))
        places <- function() list(bus = !is.null(probe$bus("uniprot_dat_cln")),
            project = file.exists(file.path(paths$source_dir, "uniprot_dat_cln.RDS")),
            cache = file.exists(file.path(cache, "uniprot_annotations.RDS")),
            tsv = file.exists(file.path(cache, "uniprot_results.tsv")),
            integration = file.exists(file.path(integration, "proteomics")) ||
                file.exists(file.path(integration, "proteomics_planted_final_s4.RDS")))
        out <- list(places_before = places(), taxon_before = probe$bus("taxon_id"), organism_before = probe$bus("organism_name"),
            before = kit277State(api, probe), enrichment_before = !is.null(probe$manager()$.context$snapshot$enrichment_results))
        folder <- kit277ExportedDesign(probe, identity)
        modals <- probe$modals()
        probe$clear_notes()
        kit277Import(api, lane, folder, 10090L, "Mus musculus")
        out$asked <- probe$modals() > modals
        out$question <- utils::tail(probe$modal_text(), 1L)
        out$cancelled <- list(places = places(), taxon = probe$bus("taxon_id"), organism = probe$bus("organism_name"),
            state = kit277State(api, probe))
        kit277Import(api, lane, folder, 10090L, "Mus musculus")
        if (out$asked) kit277Click(api, "design_matrix-confirm_design_resave")
        probe$flush()
        out$notes <- probe$notes()
        out$places_after <- places()
        diag_file <- tempfile("k277diag-")
        withr::local_envvar(MULTISCHOLAR_KIT277_DIAG_FILE = diag_file)
        out$diag_file <- diag_file
        kit277Diag("A_after_organism_save", api, probe)
        out$taxon_after <- probe$bus("taxon_id")
        out$organism_after <- probe$bus("organism_name")
        out$after <- kit277State(api, probe)
        out$enrichment_after <- !is.null(probe$manager()$.context$snapshot$enrichment_results)
        out$record_da <- kit277RecordHasDa(probe)
        out$export_da <- kit277ExportHasDa(probe)
        out$status_enrichment <- probe$status()$enrichment_analysis
        out$copy_notes_after <- kit277CopyClick(api, probe)
        kit277Diag("B_after_first_save_and_copy", api, probe)
        # The enrichment Run is refused until the DA runs again.
        kit277Enrich(api, probe, "10090")
        out$enrichment_before_da <- !is.null(probe$manager()$.context$snapshot$enrichment_results)
        kit277Diag("C_after_refused_enrichment", api, probe)
        # The DA run fetches the new organism's annotation.
        kit277RunDa(api, lane)
        probe$flush()
        kit277Diag("D_after_da_run", api, probe)
        out$da_organisms <- kit277Organisms(probe$bus("uniprot_dat_cln"))
        # The enrichment run uses it.
        kit277Enrich(api, probe, "10090")
        probe$flush()
        kit277Diag("E_after_enrichment", api, probe)
        out$enrich_organisms <- kit277Organisms(probe$bus("uniprot_dat_cln"))
        out$copy_after <- kit277Gate(api, probe, TRUE)
        out$copy_notes_rerun <- kit277CopyClick(api, probe)
        kit277Diag("F_after_second_save_and_copy", api, probe)
        out
    }))
}

test_that("KIT-277 O1, O3 (S2): an organism-only import asks, Cancel changes nothing, Save applies it and removes the earlier annotation everywhere", {
    kit277Skip("S2")
    out <- kit277OrganismScenario()
    expect_true(all(unlist(out$places_before)))
    expect_true(out$asked)
    expect_match(paste(out$question, collapse = " "), "Homo sapiens", fixed = TRUE)
    expect_match(paste(out$question, collapse = " "), "Mus musculus", fixed = TRUE)
    expect_true(all(unlist(out$cancelled$places)))
    expect_identical(out$cancelled$taxon, out$taxon_before)
    expect_identical(out$cancelled$state[c("history", "head", "entity", "status")], out$before[c("history", "head", "entity", "status")])
    expect_equal(as.numeric(out$taxon_after), 10090)
    expect_identical(out$organism_after, "Mus musculus")
    expect_false(any(unlist(out$places_after)))
    expect_identical(out$after$status$differential_expression, "pending")
    expect_identical(out$status_enrichment, "pending")
    expect_true(any(grepl("Organism saved. Run the differential analysis, then enrichment, again.", out$notes, fixed = TRUE)))
    expect_false(any(grepl("nothing changed", out$notes, fixed = TRUE)))
})

test_that("KIT-277 O4 (S2): after an organism Save the earlier DA and enrichment are not current for any consumer", {
    kit277Skip("S2")
    out <- kit277OrganismScenario()
    if (nzchar(Sys.getenv("MULTISCHOLAR_KIT277_DIAG")) && file.exists(out$diag_file)) cat(readLines(out$diag_file), sep = "\n")
    expect_identical(out$after$copy, "refused")
    expect_identical(out$after$summary, "allowed")
    expect_false(out$after$da)
    expect_false(out$enrichment_after)
    expect_false(out$record_da)
    expect_false(out$export_da)
    expect_false(out$enrichment_before_da)
    expect_identical(out$copy_after, "allowed")
    expect_true(any(grepl("Copy error", out$copy_notes_after, fixed = TRUE)))
    expect_false(any(grepl("Copy error", out$copy_notes_rerun, fixed = TRUE)))
    expect_gt(length(out$copy_notes_rerun), 0L)
})

test_that("KIT-277 O2b, O2 (S2): after an organism Save the DA run and the enrichment run use the new organism's annotation", {
    kit277Skip("S2")
    out <- kit277OrganismScenario()
    expect_identical(out$da_organisms, "Synthetic test organism")
    expect_identical(out$enrich_organisms, "Synthetic test organism")
})

test_that("KIT-277 G4 (S2, guard): an import whose organism name differs only textually, with the same taxonomy id, asks nothing", {
    kit277Skip("S2")
    out <- kit277Lane("prot_lfq", "correlation", function(api, probe) {
        folder <- kit277ExportedDesign(probe, identity)
        modals <- probe$modals()
        probe$clear_notes()
        kit277Import(api, probe$lane, folder, 9606L, "Homo sapiens (human)")
        probe$flush()
        list(modals = probe$modals() - modals, notes = probe$notes(), taxon = probe$bus("taxon_id"),
            organism = probe$bus("organism_name"))
    })
    expect_identical(out$modals, 0L)
    expect_true(any(grepl("Design saved; nothing changed.", out$notes, fixed = TRUE)))
    expect_equal(as.numeric(out$taxon), 9606)
})

# ---------------------------------------------------------------------------------------------------------------------
# Stage 3: the FASTA.
# ---------------------------------------------------------------------------------------------------------------------

kit277FastaFile <- function(api, lane) file.path(api$inputs$dir, lane$fasta_file)

#' The lane's FASTA with the first gene name replaced by `gene` (bytes differ; the table differs).
kit277ModifiedFasta <- function(api, lane, gene, to) {
    lines <- readLines(kit277FastaFile(api, lane), warn = FALSE)
    at <- grep("GN=", lines, fixed = TRUE)[[1L]]
    lines[[at]] <- sub("GN=[A-Za-z0-9_.-]+", paste0("GN=", gene), lines[[at]])
    writeLines(lines, to)
    invisible(to)
}

kit277TableText <- function(table) paste(unlist(table, use.names = FALSE), collapse = " ")

kit277FastaState <- function(probe) {
    context <- probe$manager()$.context
    dependency <- context$snapshot$file_dependencies$fasta
    projectTable <- file.path(probe$paths$source_dir, "aa_seq_tbl_final.RDS")
    list(bus = kit277TableText(probe$bus("aa_seq_tbl_final")),
        project = if (file.exists(projectTable)) kit277TableText(readRDS(projectTable)) else "<none>",
        sha = dependency$sha256 %||% "<none>", name = dependency$name %||% "<none>", history = probe$history(),
        head = probe$head(), modals = probe$modals(), notes = probe$notes())
}

kit277CleanupId <- function(lane_id) {
    ids <- vapply(guiEntrySteps(guiEntryLane(lane_id)), `[[`, "", "id")
    grep("cleanup", ids, fixed = TRUE)[[1L]]
}

kit277PreCleanupId <- function(lane_id) {
    ids <- vapply(guiEntrySteps(guiEntryLane(lane_id)), `[[`, "", "id")
    ids[kit277CleanupId(lane_id) - 1L]
}

kit277CleanupStepId <- function(lane_id) {
    ids <- vapply(guiEntrySteps(guiEntryLane(lane_id)), `[[`, "", "id")
    ids[kit277CleanupId(lane_id)]
}

kit277FastaChain <- function() {
    kit277Once("fasta-chain", function() kit277Lane("prot_dia", kit277PreCleanupId("prot_dia"), function(api, probe) {
        lane <- probe$lane
        out <- list(start = kit277FastaState(probe))
        step <- function(prepare) {
            folder <- kit277ExportedDesign(probe, identity)
            prepare(folder)
            modals <- probe$modals()
            probe$clear_notes()
            kit277Import(api, lane, folder)
            probe$flush()
            state <- kit277FastaState(probe)
            state$asked <- probe$modals() > modals
            state
        }
        new_fasta <- function(folder, gene = "GENEX1") kit277ModifiedFasta(api, lane, gene, file.path(folder, "proteins.fasta"))
        modified_table <- function(gene) {
            table <- probe$bus("aa_seq_tbl_final")
            column <- if ("gene_name" %in% names(table)) "gene_name" else utils::tail(names(table)[vapply(table, is.character, TRUE)], 1L)
            table[[column]][[1L]] <- gene
            table
        }
        expected <- tempfile("k277-expected-", fileext = ".fasta")
        kit277ModifiedFasta(api, lane, "GENEX1", expected)
        out$p1_sha <- digest::digest(file = expected, algo = "sha256")
        out$p1 <- step(function(folder) new_fasta(folder))
        out$p2 <- step(function(folder) new_fasta(folder))
        rds_folder_table <- modified_table("GENEY1")
        out$p3 <- step(function(folder) saveRDS(rds_folder_table, file.path(folder, "aa_seq_tbl_final.RDS")))
        both_table <- modified_table("GENEZ1")
        out$p9 <- step(function(folder) {
            saveRDS(both_table, file.path(folder, "aa_seq_tbl_final.RDS"))
            new_fasta(folder, "GENEQ1")
        })
        out
    }))
}

test_that("KIT-277 P1, P2 (S3): a changed FASTA before Accession Cleanup is used without a question; the same file again changes nothing", {
    kit277Skip("S3")
    out <- kit277FastaChain()
    expect_false(out$p1$asked)
    expect_true(grepl("GENEX1", out$p1$bus, fixed = TRUE))
    expect_true(grepl("GENEX1", out$p1$project, fixed = TRUE))
    expect_identical(out$p1$sha, out$p1_sha)
    expect_true(any(grepl("FASTA file saved", out$p1$notes, fixed = TRUE)))
    expect_true(any(grepl("Accession Cleanup will use it", out$p1$notes, fixed = TRUE)))
    # P2: the same bytes again: nothing asked, nothing changed.
    expect_false(out$p2$asked)
    expect_identical(out$p2$bus, out$p1$bus)
    expect_identical(out$p2$sha, out$p1$sha)
    expect_true(any(grepl("Design saved; nothing changed.", out$p2$notes, fixed = TRUE)))
})

test_that("KIT-277 P3, P9 (S3): a folder table is compared by content; the table applied is the file recorded", {
    kit277Skip("S3")
    out <- kit277FastaChain()
    expect_true(grepl("GENEY1", out$p3$bus, fixed = TRUE))
    expect_identical(out$p3$name, "aa_seq_tbl_final.RDS")
    expect_true(grepl("GENEZ1", out$p9$bus, fixed = TRUE))
    expect_false(grepl("GENEQ1", out$p9$bus, fixed = TRUE))
    expect_identical(out$p9$name, "aa_seq_tbl_final.RDS")
})

#' A lane to its Accession Cleanup, the folder with a changed FASTA, the question and its Save (`kit277FastaCut`).
kit277FastaCut <- function(lane_id, until, rerun) {
    kit277Once(paste0("fasta-cut-", lane_id), function() kit277Lane(lane_id, until, function(api, probe) {
        lane <- probe$lane
        history <- probe$history()
        cleanup <- intersect(c("protein_accession_cleaned", "peptide_accession_cleaned"), history)[[1L]]
        out <- list(history = history, cleanup = cleanup, expected = history[match(cleanup, history) - 1L],
            design = history[2L], cleaned_sha_before = kit277Try(kit277CleanupSha(probe, cleanup)),
            sha_before = probe$manager()$.context$snapshot$file_dependencies$fasta$sha256)
        integration <- kit277IntegrationDir(probe$paths)
        dir.create(file.path(integration, "proteomics"), recursive = TRUE, showWarnings = FALSE)
        writeLines("planted", file.path(integration, "proteomics", "manifest.json"))
        folder <- kit277ExportedDesign(probe, identity)
        kit277ModifiedFasta(api, lane, "GENEX1", file.path(folder, "proteins.fasta"))
        out$new_sha <- digest::digest(file = file.path(folder, "proteins.fasta"), algo = "sha256")
        before <- list(history = probe$history(), head = probe$head(), entity = probe$entity(), status = probe$status(),
            bus = kit277TableText(probe$bus("aa_seq_tbl_final")))
        modals <- probe$modals()
        probe$clear_notes()
        kit277Import(api, lane, folder)
        out$asked <- probe$modals() > modals
        out$question <- utils::tail(probe$modal_text(), 1L)
        out$cancelled <- list(history = probe$history(), head = probe$head(), entity = probe$entity(), status = probe$status(),
            bus = kit277TableText(probe$bus("aa_seq_tbl_final")))
        out$before <- before
        kit277Import(api, lane, folder)
        if (out$asked) kit277Click(api, "design_matrix-confirm_design_resave")
        probe$flush()
        out$notes <- probe$notes()
        out$history_after <- probe$history()
        out$head_after <- probe$head()
        out$status_after <- probe$status()
        out$integration_after <- file.exists(file.path(integration, "proteomics"))
        out$sha_after <- probe$manager()$.context$snapshot$file_dependencies$fasta$sha256
        out$bus_after <- kit277TableText(probe$bus("aa_seq_tbl_final"))
        out$project_after <- kit277Try(kit277TableText(readRDS(file.path(probe$paths$source_dir, "aa_seq_tbl_final.RDS"))))
        out$da_status <- kit277Try(workflowDaResultStatus(probe$manager()))
        out$copy <- kit277Gate(api, probe, TRUE)
        out$head_rule <- kit277HeadRule(api, probe)
        out$record_da <- kit277RecordHasDa(probe)
        if (identical(lane_id, "prot_lfq")) out$copy_notes <- kit277CopyClick(api, probe)
        if (isTRUE(rerun)) {
            id <- kit277CleanupStepId(lane_id)
            api$set(stats::setNames(list(probe$steps[[id]]$tab), probe$tabs))
            # The lane's step presses its button with a fixed count; the press is raised from 0 so that it fires again.
            again <- api
            again$set <- function(values) {
                actions <- grep("-apply_", names(values), value = TRUE)
                if (length(actions)) api$set(stats::setNames(as.list(rep(0L, length(actions))), actions))
                api$set(values)
            }
            kit277Try(probe$steps[[id]]$run(again))
            probe$flush()
            out$rerun_history <- probe$history()
            out$rerun_cleanup_sha <- kit277Try(kit277CleanupSha(probe, cleanup))
            out$rerun_bus <- kit277TableText(probe$bus("aa_seq_tbl_final"))
        }
        out
    }))
}

kit277CleanupSha <- function(probe, state) {
    manager <- probe$manager()
    withWorkflowStateScope(manager, function() {
        readWorkflowStateMetadata(manager, getWorkflowStateReference(manager, state))$config$fasta_sha256
    })
}

kit277FastaRoutes <- list(
    prot_lfq = list(until = "da", rerun = TRUE, design = "protein_s4_initial"),
    prot_tmt = list(until = NULL, rerun = FALSE, design = "protein_s4_initial"),
    prot_dia = list(until = NULL, rerun = FALSE, design = NULL),
    prot_dia_limpa = list(until = NULL, rerun = TRUE, design = NULL))

for (lane_id in names(kit277FastaRoutes)) {
    test_that(paste("KIT-277 P4 (S3):", lane_id, "a changed FASTA after Accession Cleanup asks; Save goes back to the step before the cleanup"), {
        kit277Skip("S3")
        route <- kit277FastaRoutes[[lane_id]]
        until <- route$until %||% kit277CleanupStepId(lane_id)
        out <- kit277FastaCut(lane_id, until, route$rerun)
        expect_true(out$asked, info = lane_id)
        expect_match(paste(out$question, collapse = " "), "proteins.fasta", fixed = TRUE, info = lane_id)
        expect_match(paste(out$question, collapse = " "), "Accession Cleanup", fixed = TRUE, info = lane_id)
        expect_identical(out$cancelled, out$before, info = lane_id)
        expect_identical(out$head_after, out$expected, info = lane_id)
        expect_identical(utils::tail(out$history_after, 1L), out$expected, info = lane_id)
        if (!is.null(route$design)) expect_identical(out$expected, route$design, info = lane_id)
        if (is.null(route$design)) expect_false(identical(out$expected, out$design), info = lane_id)
        expect_identical(out$status_after$quality_control, "pending", info = lane_id)
        expect_false(out$integration_after, info = lane_id)
        expect_true(grepl("GENEX1", out$bus_after, fixed = TRUE), info = lane_id)
        expect_true(grepl("GENEX1", out$project_after, fixed = TRUE), info = lane_id)
        expect_identical(out$sha_after, out$new_sha, info = lane_id)
        expect_true(any(grepl("FASTA file saved", out$notes, fixed = TRUE)), info = lane_id)
        if (identical(lane_id, "prot_lfq")) {
            # The DA's source left the history: Copy refused, the record lists no DA.
            expect_identical(out$da_status, "stale")
            expect_identical(out$copy, "refused")
            expect_false(out$record_da)
            expect_true(any(grepl("Copy error", out$copy_notes, fixed = TRUE)))
        }
        # R15: the head is not one of the DA's four states, so the Summary Save, Report and Export Session State (which
        # call this rule) are refused, with the class and a text.
        expect_match(out$head_rule, "^refused: .+", info = lane_id)
    })
}

for (lane_id in c("prot_lfq", "prot_dia_limpa")) {
    test_that(paste("KIT-277 P5, P6 (S3):", lane_id, "the cleanup run again names and uses the file that was applied"), {
        kit277Skip("S3")
        route <- kit277FastaRoutes[[lane_id]]
        out <- kit277FastaCut(lane_id, route$until %||% kit277CleanupStepId(lane_id), route$rerun)
        # In process the cleanup reads the bus table; on the worker the project copy: both are the new file's table.
        expect_true(grepl("GENEX1", out$rerun_bus, fixed = TRUE), info = lane_id)
        expect_identical(out$project_after, out$bus_after, info = lane_id)
        expect_identical(out$rerun_cleanup_sha, out$new_sha, info = lane_id)
        expect_identical(out$rerun_cleanup_sha, out$sha_after, info = lane_id)
    })
}

test_that("KIT-277 P7 (S3): a full import with a new FASTA file into a project that already holds a table applies and records that file", {
    kit277Skip("S3")
    out <- kit277Lane("prot_lfq", "design", function(api, probe) {
        lane <- probe$lane
        folder <- kit277ExportedDesign(probe, identity)
        kit277ModifiedFasta(api, lane, "GENEX1", file.path(folder, "proteins.fasta"))
        expected <- digest::digest(file = file.path(folder, "proteins.fasta"), algo = "sha256")
        kit277Import(api, lane, folder)
        probe$flush()
        state <- kit277FastaState(probe)
        state$expected <- expected
        state
    })
    expect_true(grepl("GENEX1", out$bus, fixed = TRUE))
    expect_true(grepl("GENEX1", out$project, fixed = TRUE))
    expect_identical(out$sha, out$expected)
})

test_that("KIT-277 P8 (S3): with no application session a limited change after Accession Cleanup is refused and nothing is applied", {
    kit277Skip("S3")
    ran <- character()
    limited <- list(fasta = list(from = "a.fasta", to = "b.fasta", steps = "Accession Cleanup"))
    condition <- tryCatch(runOmicsDesignResaveGate(NULL, "unchanged", TRUE, commit = function(cut) ran <<- c(ran, "commit"),
        unchanged = function() ran <<- c(ran, "unchanged"), limited = limited), error = function(e) e)
    expect_s3_class(condition, "omics_design_resave_unconfirmed")
    expect_identical(ran, character())
})

test_that("KIT-277 R1 (S3, C05): a confirmed import refused at the checkpoint puts back the project's copies and the recorded input", {
    kit277Skip("S3")
    out <- kit277Lane("prot_lfq", "correlation", function(api, probe) {
        lane <- probe$lane
        paths <- probe$paths
        collide <- function(design) {
            row <- which(design$group == design$group[[1L]])[[2L]]
            design$replicates[[row]] <- design$replicates[[1L]]
            if (!"tech_reps" %in% names(design)) design$tech_reps <- NA_integer_
            design$tech_reps[c(1L, row)] <- NA_integer_
            if ("tech_rep_group" %in% names(design)) design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
            design
        }
        folder <- kit277ExportedDesign(probe, collide)
        saveRDS(data.frame(planted = "KIT-277 aa"), file.path(folder, "aa_seq_tbl_final.RDS"))
        saveRDS(data.frame(planted = "KIT-277 meta"), file.path(folder, "fasta_metadata.RDS"))
        saveRDS(data.frame(planted = "KIT-277 uniprot"), file.path(folder, "uniprot_dat_cln.RDS"))
        copies <- function() vapply(c("aa_seq_tbl_final.RDS", "fasta_metadata.RDS", "uniprot_dat_cln.RDS"), function(name) {
            path <- file.path(paths$source_dir, name)
            if (file.exists(path)) digest::digest(file = path, algo = "sha256") else "<none>"
        }, "")
        before <- list(copies = copies(), input = probe$manager()$.context$snapshot$file_dependencies$fasta)
        kit277Import(api, lane, folder)
        kit277Click(api, "design_matrix-confirm_design_resave")
        probe$flush()
        list(before = before, copies = copies(), input = probe$manager()$.context$snapshot$file_dependencies$fasta,
            notes = probe$notes())
    })
    expect_true(any(grepl("Error during import", out$notes, fixed = TRUE)))
    expect_identical(out$copies, out$before$copies)
    expect_identical(out$input, out$before$input)
})

test_that("KIT-277 R1b (S3, E-8): a confirmed FASTA Save refused after the files were written puts them back", {
    kit277Skip("S3")
    lane_id <- "prot_lfq"
    out <- kit277Lane(lane_id, kit277CleanupStepId(lane_id), function(api, probe) {
        lane <- probe$lane
        paths <- probe$paths
        folder <- kit277ExportedDesign(probe, identity)
        kit277ModifiedFasta(api, lane, "GENEX1", file.path(folder, "proteins.fasta"))
        copies <- function() vapply(c("aa_seq_tbl_final.RDS", "fasta_metadata.RDS"), function(name) {
            path <- file.path(paths$source_dir, name)
            if (file.exists(path)) digest::digest(file = path, algo = "sha256") else "<none>"
        }, "")
        before <- list(copies = copies(), input = probe$manager()$.context$snapshot$file_dependencies$fasta,
            bus = kit277TableText(probe$bus("aa_seq_tbl_final")), history = probe$history())
        modals <- probe$modals()
        kit277Import(api, lane, folder)
        asked <- probe$modals() > modals
        if (asked) {
            testthat::with_mocked_bindings(kit277Click(api, "design_matrix-confirm_design_resave"),
                cutOmicsDesignHistory = function(...) stop("KIT-277 planted failure"), .package = "MultiScholaR")
        }
        probe$flush()
        list(asked = asked, before = before, copies = copies(), input = probe$manager()$.context$snapshot$file_dependencies$fasta,
            bus = kit277TableText(probe$bus("aa_seq_tbl_final")), history = probe$history(), notes = probe$notes())
    })
    expect_true(out$asked)
    expect_true(any(grepl("Error during import", out$notes, fixed = TRUE)))
    expect_identical(out$copies, out$before$copies)
    expect_identical(out$input, out$before$input)
    expect_identical(out$bus, out$before$bus)
    expect_identical(out$history, out$before$history)
})

test_that("KIT-277 G3 (guard): a first import into a project with no table processes and applies the FASTA file as before", {
    kit277Skip("S3")
    root <- withr::local_tempdir()
    import_dir <- file.path(root, "import")
    source_dir <- file.path(root, "scripts")
    results_dir <- file.path(root, "results")
    for (path in c(import_dir, source_dir, results_dir)) dir.create(path)
    fasta <- testthat::test_path("..", "testdata", "e2e", "prot_tmt", "proteins.fasta")
    workflow_data <- new.env(parent = emptyenv())
    session <- shiny::MockShinySession$new()
    suppressMessages(shiny::withReactiveDomain(session, hydrateProtDesignImportedFastaSidecar(workflow_data, import_dir,
        source_dir, resultsDir = results_dir, fastaPath = fasta, organismName = "Homo sapiens")))
    expect_true(is.data.frame(workflow_data$aa_seq_tbl_final))
    expect_gt(nrow(workflow_data$aa_seq_tbl_final), 0L)
    expect_true(file.exists(file.path(source_dir, "aa_seq_tbl_final.RDS")))
})
