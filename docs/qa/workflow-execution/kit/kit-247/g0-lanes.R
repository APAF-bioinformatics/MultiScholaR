# KIT-247 gate G0 (ticket r1.2, section 10): observation of the base, no package code. A design saved again after
# Quality Control and normalisation, through the GUI-entry driver (tests/testthat/helper-gui-entry.R: each omic's public
# module server under shiny::testServer, through the namespaced inputs only). Each case is a fresh lane run: the lane's
# own steps up to a named step (the prefix), then this script's presses (the re-save and what follows), each recorded.
#
# Binding: v0.9-staging 16112917 plus this folder only (coordinator's binding, r1.2; KIT-222 changes no lane's design).
#
# PREDICTIONS are written below (k247Predictions) before any value is read. Each has a regular expression applied to the
# observed value; "held" is TRUE/FALSE, or NA for a value recorded with no prediction.
#
# Instruments (observation only; each calls the package function it wraps and returns its value unchanged): spies on
# diaNormalizationRequest (the proteomics run's source and the KIT-182 notice), omicNormApplyRunStart (the lipid, metab
# and trx run's start), revertProtNormStateManagerToPreNormalization, resetOmicNorm, resetTrxNormalisation (the recorded
# start at Reset) and startWorkflowSummaryIntegrationExport (the state the Summary's Save takes). Entities are labelled
# "<step>:<state>" by the step after which they first appear in the owner's states map.
#
# Usage (package root, on a runner): Rscript docs/qa/workflow-execution/kit/kit-247/g0-lanes.R <g01|g02|g03> <out dir>

k247Args <- commandArgs(trailingOnly = TRUE)
k247Group <- if (length(k247Args) >= 1L) k247Args[[1L]] else "g01"
k247Out <- if (length(k247Args) >= 2L) k247Args[[2L]] else "kit247-g0"
dir.create(k247Out, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------------------------------------------------
# Predictions (by reading at 16112917; KIT-247 r1.2 sections 2.2 to 2.8 and Re-check 1 R1-08). Written before any run.
# ---------------------------------------------------------------------------------------------------------------------

k247P <- function(lane, case, key, predicted, pattern) {
    data.frame(lane = lane, case = case, key = key, predicted = predicted, pattern = pattern, stringsAsFactors = FALSE)
}

k247Predictions <- function() {
    rows <- list()
    add <- function(...) rows[[length(rows) + 1L]] <<- k247P(...)
    resave <- function(lane, case, design) {
        add(lane, case, "resave.head", "the head is the re-committed design state", paste0("^", design, "$"))
        add(lane, case, "resave.history_unchanged", "the active history is unchanged: every later state of the previous design stays in place (2.2)", "^TRUE$")
        add(lane, case, "resave.later_entities_unchanged", "every later state keeps its entity", "^TRUE$")
        add(lane, case, "resave.status_unchanged_except_design", "no tab status changes (design_matrix stays complete)", "^TRUE$")
        add(lane, case, "resave.notes", "the success text only", "saved successfully")
    }
    # G0-1 proteomics, label-free (prot_lfq): moved KO_5 from KO to WT.
    l <- "prot_lfq"; d <- "protein_s4_initial"
    resave(l, "main", d)
    add(l, "main", "resave.design_entity_changed", "the design state is a new entity", "^TRUE$")
    add(l, "main", "reset.head", "Reset (KIT-225) returns the previous design's Protein Replicate Filter table P4", "^protein_replicate_filtered$")
    add(l, "main", "reset.head_label", "P4 is the prefix's entity (previous design)", "^qc_protein_[^:]*:protein_replicate_filtered$")
    add(l, "main", "reset.head_group_moved", "the moved sample is in its OLD group in that table", "^KO$")
    add(l, "main", "run.source_label", "the run reads P4", "^qc_protein_[^:]*:protein_replicate_filtered$")
    add(l, "main", "run.source_group_moved", "old group", "^KO$")
    add(l, "main", "run.notice_kit182", "no KIT-182 notice (normalization_complete FALSE after Reset)", "^FALSE$")
    add(l, "main", "corr.head", "the correlation filter applies", "^correlation_filtered$")
    add(l, "main", "da.result_changed", "the DA runs, with no message", "^TRUE$")
    add(l, "main", "da.group_moved", "the DA model's design has the moved sample in its OLD group", "^KO$")
    add(l, "main", "record.design_group_moved", "the record's design rows are the previous design's (INFERRED, P08)", "^KO$")
    resave(l, "b", d)
    add(l, "b", "corr.head", "normalisation and the correlation filter run first", "^correlation_filtered$")
    add(l, "b", "da.result_changed", "the DA refuses (a renamed group is not in the table's design)", "^FALSE$")
    resave(l, "c", d)
    add(l, "c", "da.result_changed", "DA pressed directly after the re-save does not run", "^FALSE$")
    add(l, "c", "da.notes", "refused: no filtered protein data is loaded", "no filtered protein data is loaded")
    resave(l, "d", d)
    add(l, "d", "qcre.head", "Accession Cleanup re-applied: its state, at its old position", "accession")
    add(l, "d", "run.source_label", "the run reads the re-applied Accession Cleanup table", "^qcre:")
    add(l, "d", "record.steps", "the record lists the previous design's later QC steps between the re-applied step and the new normalisation", "accession.* > .*replicate_filtered.* > .*normalized")
    resave(l, "e", d)
    add(l, "e", "run.source_label", "path 1: the run reads the NEW design state, before any Quality Control", "^resave:protein_s4_initial$")
    add(l, "e", "run.source_group_moved", "the new design's group", "^WT$")
    add(l, "e", "run.notice_kit182", "the KIT-182 notice fires ('uses your current quality-control result', false here)", "^TRUE$")
    add(l, "e", "run.prereq_warning", "the tab-entry warning stays silent (QC still complete)", "^FALSE$")
    add(l, "e", "run.status", "Quality Control still shown complete", "quality_control=complete")
    add(l, "e", "run.states_between_design_and_run", "the previous design's QC states sit between the design state and the run's states", "^[1-9]")
    add(l, "e", "record.qc_steps_before_run", "the record lists the previous design's QC steps before the new normalisation", "^[1-9]")
    resave(l, "f0", d)
    add(l, "f0", "summary.refused", "no DA: Save goes ahead (status 'none'; proteomics has no head guard)", "^FALSE$")
    add(l, "f0", "summary.state_used", "on the previous design's correlation_filtered", "^correlation_filtered$")
    add(l, "f0", "summary.object_group_moved", "whose design has the moved sample in its old group", "^KO$")
    add(l, "f1", "before.integration_files", "the Save before the re-save wrote the integration input", "\\S")
    resave(l, "f1", d)
    add(l, "f1", "summary.refused", "a DA before the re-save: Save refused as stale (workflow_summary_da_stale; R1-08)", "^TRUE$")
    add(l, "f1", "summary.integration_files", "and the omic's integration input removed (R1-08)", "^<none>$")
    resave(l, "g", d)
    add(l, "g", "resave.contrast_count", "the bus holds both comparisons", "^2$")
    resave(l, "h", d)
    add(l, "h", "resave.components_all_identical", "an unchanged re-save yields identical table identities (stop rule 5 if not)", "^TRUE$")
    # G0-2 proteomics, DIA with limpa (prot_dia_limpa): moved KO_4 from KO to WT.
    l <- "prot_dia_limpa"; d <- "raw_data_s4"
    resave(l, "main", d)
    add(l, "main", "reset.head", "Reset (KIT-225) returns the previous design's Peptide Count Filter table L7", "^peptide_count_filtered$")
    add(l, "main", "reset.head_group_moved", "old group", "^KO$")
    add(l, "main", "run.source_label", "the pass-through designates L7 (previous design)", ":peptide_count_filtered$")
    add(l, "main", "run.source_group_moved", "old group", "^KO$")
    add(l, "main", "da.result_changed", "numbers: the DA runs (newly reachable through KIT-225, INFERRED)", "^TRUE$")
    add(l, "main", "da.group_moved", "on the previous design's samples and groups", "^KO$")
    resave(l, "e", d)
    add(l, "e", "run.source_label", "path 1: the run starts at the new design's peptide state raw_data_s4", "^resave:raw_data_s4$")
    add(l, "e", "run.notice_kit182", "no KIT-182 notice (a peptide table)", "^FALSE$")
    # G0-3 lipidomics, metabolomics (KO_3 to WT) and transcriptomics (Sera_03's two runs to RPMI).
    for (l in c("lipid_canonical", "metab_lc", "trx_screen")) {
        prefix <- switch(l, lipid_canonical = "lipid", metab_lc = "metab", trx_screen = "trx")
        d <- paste0(prefix, "_raw_data_s4")
        old <- if (identical(l, "trx_screen")) "Sera" else "KO"
        new <- if (identical(l, "trx_screen")) "RPMI" else "WT"
        resave(l, "main", d)
        add(l, "main", "reset.head", "Reset re-commits the recorded start whatever the history holds (2.4)", paste0("^", prefix, "_reset$"))
        add(l, "main", "reset.head_group_moved", "that table is the previous design's", paste0("^", old, "$"))
        add(l, "main", "run.source_group_moved", "the run reads it: old group", paste0("^", old, "$"))
        if (!identical(l, "trx_screen")) add(l, "main", "corr.head", "the correlation filter applies", paste0("^", prefix, "_correlation_filtered$"))
        add(l, "main", "da.result_changed", "the DA runs, with no message", "^TRUE$")
        add(l, "main", "da.group_moved", "on the previous design's groups", paste0("^", old, "$"))
        resave(l, "e", d)
        add(l, "e", "run.source_label", "path 1: the run reads the NEW design state, before any Quality Control", paste0("^resave:", d, "$"))
        add(l, "e", "run.source_group_moved", "the new design's group", paste0("^", new, "$"))
        add(l, "e", "run.status", "Quality Control still shown complete", "quality_control=complete")
        add(l, "e", "run.states_between_design_and_run", "the previous design's QC states sit between the design state and the run's states", "^[1-9]")
        add(l, "e", "record.qc_steps_before_run", "the record lists them before the new normalisation", "^[1-9]")
        resave(l, "f0", d)
        add(l, "f0", "summary.refused", "Save refused at a design head (assertWorkflowSummaryHeadFinal)", "^TRUE$")
        add(l, "f1", "before.integration_files", "the Save before the re-save wrote the integration input", "\\S")
        resave(l, "f1", d)
        add(l, "f1", "summary.refused", "Save refused at a design head", "^TRUE$")
        add(l, "f1", "summary.integration_files", "and the omic's integration input removed", "^<none>$")
        resave(l, "g", d)
        add(l, "g", "resave.contrast_count", "the bus holds both comparisons", "^2$")
        resave(l, "h", d)
        add(l, "h", "resave.components_all_identical", "an unchanged re-save yields identical table identities (stop rule 5 if not)", "^TRUE$")
    }
    do.call(rbind, rows)
}

k247PredictionTable <- k247Predictions()

# ---------------------------------------------------------------------------------------------------------------------
# The package (base) and the helpers.
# ---------------------------------------------------------------------------------------------------------------------

k247Log <- function(...) cat(sprintf("[KIT-247 G0 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)), file = stderr())
k247Log("loading the package")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)

k247Rows <- new.env(parent = emptyenv())
k247Rows$rows <- list()
k247Text <- function(value) {
    if (is.null(value)) return("<null>")
    value <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(value)), collapse = " | "))
    if (!nzchar(value)) "<none>" else substr(value, 1L, 2000L)
}
k247Put <- function(lane, case, key, value) {
    value <- k247Text(value)
    k247Rows$rows[[length(k247Rows$rows) + 1L]] <- data.frame(lane = lane, case = case, key = key, value = value,
        stringsAsFactors = FALSE)
    k247Log(lane, " ", case, " ", key, " = ", substr(value, 1L, 300L))
}
k247Try <- function(expr) tryCatch(expr, error = function(e) paste("<error>", conditionMessage(e)))

#' The lane's constants: the moved runs (fixture names), the old and new group, the comparisons.
k247LaneSpec <- function(lane_id) {
    switch(lane_id,
        prot_lfq = list(moved = "KO_5", old = "KO", new = "WT", contrast = c("KO", "WT")),
        prot_dia_limpa = list(moved = "KO_4", old = "KO", new = "WT", contrast = c("KO", "WT")),
        lipid_canonical = list(moved = "KO_3", old = "KO", new = "WT", contrast = c("KO", "WT")),
        metab_lc = list(moved = "KO_3", old = "KO", new = "WT", contrast = c("KO", "WT")),
        trx_screen = list(moved = c("Sera_03_r1", "Sera_03_r2"), old = "Sera", new = "RPMI", contrast = c("Sera", "RPMI"),
            new_replicate = 4L, new_set = "RPMI_04"))
}

#' The group of `sample` in the design matrix of the table `reference` names (read through the owner).
k247GroupOf <- function(manager, reference, sample) {
    if (is.null(reference)) return("<null>")
    k247Try(withWorkflowStateScope(manager, function() {
        # VAL-54 instrument fix (g03 run 45594081675: rebuilding a trx object for its design alone failed validation,
        # "counts must be numeric"): the design matrix component is read alone, as buildWorkflowAnalysisRecord reads it.
        component <- readComponentMetadata(manager, reference)$data$components$design_matrix
        design <- as.data.frame(materializeWorkflowTable(manager$.context, component), stringsAsFactors = FALSE)
        hit <- which(vapply(seq_len(nrow(design)), function(i) any(as.character(unlist(design[i, ])) %in% sample), NA))
        if (!length(hit)) "<absent>" else paste(unique(as.character(design$group[hit])), collapse = "|")
    }))
}

#' The group of `sample` in an S4 object's design matrix, or in the table a saved reference names.
k247ObjectGroupOf <- function(manager, object, sample) {
    if (isSavedWorkflowReference(object)) return(k247GroupOf(manager, object, sample))
    k247Try({
        design <- as.data.frame(methods::slot(object, "design_matrix"), stringsAsFactors = FALSE)
        hit <- which(vapply(seq_len(nrow(design)), function(i) any(as.character(unlist(design[i, ])) %in% sample), NA))
        if (!length(hit)) "<absent>" else paste(unique(as.character(design$group[hit])), collapse = "|")
    })
}

#' The component descriptors of a state's record (the table identities KIT-212's predicate compares).
k247Components <- function(manager, reference) {
    if (is.null(reference)) return(NULL)
    tryCatch(withWorkflowStateScope(manager, function() readComponentMetadata(manager, reference)$data$components),
        error = function(e) NULL)
}

# ---------------------------------------------------------------------------------------------------------------------
# The case runner (after guiEntryRunLane; the lane's steps up to `until`, then `custom`).
# ---------------------------------------------------------------------------------------------------------------------

k247RunCase <- function(lane_id, case_id, until, custom) {
    started <- Sys.time()
    k247Log("case ", lane_id, " ", case_id, " starts")
    lane <- guiEntryLane(lane_id)
    omic <- lane$omic_type
    spec <- guiEntryServerSpec(omic)
    standard <- guiEntrySteps(lane)
    ids <- vapply(standard, `[[`, "", "id")
    if (!until %in% ids) stop("the lane has no step ", until, ": ", paste(ids, collapse = ", "))
    prefix_steps <- standard[seq_len(match(until, ids))]
    root <- normalizePath(withr::local_tempdir(pattern = paste0("k247-", lane_id, "-")))
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
    withr::local_options(warn = 1)
    withr::local_seed(20261107L)
    project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", lane_id), omic_types = omic,
        label = NULL, force = TRUE, reuse_existing = TRUE))
    paths <- project_dirs[[omic]]
    integration_dir <- paths$integration_dir %||% file.path(paths$base_dir, "integration")
    old_dirs <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
    assign("project_dirs", project_dirs, envir = globalenv())
    withr::defer(if (is.null(old_dirs)) rm("project_dirs", envir = globalenv()) else assign("project_dirs", old_dirs, envir = globalenv()))
    if (identical(omic, "proteomics")) {
        guiEntryInstallUniprotCache(paths)
        install_enrichment_doubles()
    }
    log <- new.env(parent = emptyenv())
    log$updates <- list(); log$notes <- list(); log$modals <- integer(); log$choices <- list()
    guiEntryInstallBrowserDouble(log)

    # The spies (observation only).
    spy <- new.env(parent = emptyenv())
    spy$events <- list()
    push <- function(kind, ...) spy$events[[length(spy$events) + 1L]] <- c(list(kind = kind), list(...))
    ns <- asNamespace("MultiScholaR")
    o_request <- get("diaNormalizationRequest", ns)
    o_start <- get("omicNormApplyRunStart", ns)
    o_revert <- get("revertProtNormStateManagerToPreNormalization", ns)
    o_reset_omic <- get("resetOmicNorm", ns)
    o_reset_trx <- get("resetTrxNormalisation", ns)
    o_export <- get("startWorkflowSummaryIntegrationExport", ns)
    entity <- function(reference) if (isSavedWorkflowReference(reference)) reference$entity_id else NULL
    testthat::local_mocked_bindings(
        diaNormalizationRequest = function(manager, settings, norm_data, grouping) {
            previous <- entity(shiny::isolate(norm_data$normalization_input_reference))
            complete <- isTRUE(shiny::isolate(norm_data$normalization_complete))
            notice <- FALSE
            result <- withCallingHandlers(o_request(manager, settings, norm_data, grouping),
                prot_norm_input_current = function(m) notice <<- TRUE)
            push("prot_run", previous = previous, complete_before = complete, notice = notice,
                source = entity(shiny::isolate(norm_data$normalization_input_reference)))
            result
        },
        omicNormApplyRunStart = function(binding, manager, normData, addLog) {
            before <- manager$current_state
            start <- entity(shiny::isolate(normData$post_filter_obj))
            result <- o_start(binding, manager, normData, addLog)
            push("omic_run", head_before = before, recorded_start = start, head_after = manager$current_state,
                source = entity(getWorkflowStateReference(manager)), refused = inherits(result, "error"))
            result
        },
        revertProtNormStateManagerToPreNormalization = function(workflowData, normData = NULL, messageFn = message) {
            start <- if (is.null(normData)) NULL else entity(shiny::isolate(normData$normalization_input_reference))
            push("prot_reset", recorded_start = start)
            o_revert(workflowData = workflowData, normData = normData, messageFn = messageFn)
        },
        resetOmicNorm = function(workflowData, normData, addLog, config) {
            push("omic_reset", recorded_start = entity(shiny::isolate(normData$post_filter_obj)))
            o_reset_omic(workflowData, normData, addLog, config)
        },
        resetTrxNormalisation = function(manager, normData, config, addLog, notify, setStatus) {
            push("omic_reset", recorded_start = entity(shiny::isolate(normData$post_filter_obj)))
            o_reset_trx(manager, normData, config, addLog, notify, setStatus)
        },
        startWorkflowSummaryIntegrationExport = function(config, workflowData, omicType, path, stateName, object, ...) {
            push("summary_export", state = stateName, object = object, manager = workflowData$state_manager)
            o_export(config, workflowData, omicType, path, stateName, object, ...)
        },
        .package = "MultiScholaR")

    laneSpec <- k247LaneSpec(lane_id)
    route <- if (isTRUE(lane$use_limpa)) "DIA_LIMPA" else lane$workflow_type
    design_state <- if (identical(omic, "proteomics")) {
        if (route %in% c("DIA", "DIA_LIMPA")) "raw_data_s4" else "protein_s4_initial"
    } else omicGuiConfig(omic)$design$state_name
    put <- function(key, value) k247Put(lane_id, case_id, key, value)
    inputs <- list(dir = guiEntryFixtureDir(lane), scratch = file.path(root, "inputs"), root = root, lane = lane, paths = paths)
    dir.create(inputs$scratch, showWarnings = FALSE)
    roots <- c(root, paths$base_dir)
    outcome <- tryCatch(suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(spec$server,
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
        manager <- function() shiny::isolate(workflow_data$state_manager)
        snapshot <- function() manager()$.context$snapshot
        history <- function() guiEntryHistory(manager())
        registry <- new.env(parent = emptyenv())
        register <- function(step_id) {
            states <- snapshot()$states
            for (name in names(states)) {
                id <- states[[name]]$entity_id
                if (!is.null(id) && !exists(id, envir = registry, inherits = FALSE)) assign(id, paste0(step_id, ":", name), envir = registry)
            }
        }
        label <- function(id) {
            if (is.null(id)) return("<null>")
            if (exists(id, envir = registry, inherits = FALSE)) get(id, envir = registry) else paste0("other:", substr(id, 1L, 8L))
        }
        entities <- function() {
            states <- snapshot()$states
            vapply(history(), function(name) states[[name]]$entity_id %||% "<none>", "")
        }
        status_text <- function() {
            status <- shiny::isolate(workflow_data$tab_status)
            paste(paste0(names(status), "=", vapply(status, function(x) paste(as.character(x), collapse = ","), "")), collapse = "; ")
        }
        notes_text <- function() vapply(log$notes, function(note) paste0("[", note$type, "] ", substr(guiEntryMask(note$message, roots), 1L, 400L)), "")
        integration_files <- function() {
            files <- list.files(integration_dir, recursive = TRUE)
            if (!length(files)) "<none>" else paste(files, collapse = ", ")
        }
        da_reference <- function() snapshot()$da_results
        da_source <- function() {
            reference <- da_reference()
            if (!isSavedWorkflowReference(reference)) return(NULL)
            tryCatch(withWorkflowStateScope(manager(), function() resolveWorkflowStateValue(manager()$.context, reference)$source),
                error = function(e) NULL)
        }
        ctx <- new.env(parent = emptyenv())
        ctx$lane <- lane; ctx$omic <- omic; ctx$spec <- laneSpec; ctx$route <- route; ctx$design_state <- design_state
        ctx$moved <- laneSpec$moved

        # The lane's own steps (the prefix).
        for (step in prefix_steps) {
            if (!is.null(step$tab)) api$set(stats::setNames(list(step$tab), spec$tabs))
            step$run(api)
            session$flushReact()
            register(step$id)
            log$notes <- list()
        }
        if (identical(omic, "proteomics")) ctx$moved <- e2e_match_available_values(laneSpec$moved,
            api$choices("design_matrix-builder-selected_runs"))
        # trx sums technical replicates into units named by tech_rep_group at the design (func_trx_design.R:126-149): the
        # moved runs are found in a table's design under their old unit (Sera_03) or their new one (RPMI_04).
        ctx$probe <- if (identical(omic, "transcriptomics")) c(ctx$moved, "Sera_03", laneSpec$new_set) else ctx$moved[[1L]]
        before <- list(head = manager()$current_state, history = history(), entities = entities(),
            status = shiny::isolate(workflow_data$tab_status), design = getWorkflowStateReference(manager(), design_state),
            da = entity(da_reference()))
        ctx$before <- before
        put("before.head", before$head)
        put("before.history", paste(before$history, collapse = " > "))
        put("before.status", status_text())
        put("before.integration_files", integration_files())
        put("before.design_group_moved", k247GroupOf(manager(), before$design, ctx$probe))
        put("moved_runs", ctx$moved)
        put("prefix_seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))

        # This case's presses.
        for (step in custom) {
            spy$events <- list()
            da_before <- entity(da_reference())
            head_before <- manager()$current_state
            head_before_ref <- getWorkflowStateReference(manager())
            if (!is.null(step$tab)) api$set(stats::setNames(list(step$tab), spec$tabs))
            error <- tryCatch({ step$run(api, ctx); NULL }, error = function(e) conditionMessage(e))
            session$flushReact()
            register(step$id)
            k <- function(name) paste0(step$id, ".", name)
            if (!is.null(error)) put(k("error"), error)
            head_ref <- getWorkflowStateReference(manager())
            put(k("head"), manager()$current_state)
            put(k("head_label"), label(entity(head_ref)))
            put(k("head_group_moved"), k247GroupOf(manager(), head_ref, ctx$probe))
            put(k("history"), paste(history(), collapse = " > "))
            put(k("history_labels"), paste(vapply(entities(), label, ""), collapse = " > "))
            put(k("status"), status_text())
            put(k("notes"), notes_text())
            put(k("modals"), length(log$modals))
            if (identical(step$kind, "resave")) {
                now <- entities()
                put(k("head_is_design_state"), identical(manager()$current_state, design_state))
                put(k("history_unchanged"), identical(history(), before$history))
                later <- setdiff(names(before$entities), design_state)
                put(k("later_entities_unchanged"), identical(unname(now[later]), unname(before$entities[later])))
                put(k("design_entity_changed"), !identical(entity(getWorkflowStateReference(manager(), design_state)), entity(before$design)))
                status <- shiny::isolate(workflow_data$tab_status)
                fields <- setdiff(union(names(status), names(before$status)), "design_matrix")
                put(k("status_unchanged_except_design"), identical(status[fields], before$status[fields]))
                old <- k247Components(manager(), before$design)
                new <- k247Components(manager(), getWorkflowStateReference(manager(), design_state))
                fields <- union(names(old), names(new))
                same <- vapply(fields, function(f) identical(old[[f]], new[[f]]), NA)
                put(k("components_identical"), paste0(fields, "=", same))
                put(k("components_all_identical"), length(fields) > 0L && all(same))
                put(k("design_group_moved"), k247GroupOf(manager(), getWorkflowStateReference(manager(), design_state), ctx$probe))
                contrasts <- shiny::isolate(workflow_data$contrasts_tbl)
                put(k("contrast_count"), if (is.data.frame(contrasts)) nrow(contrasts) else 0L)
                put(k("contrasts"), if (is.data.frame(contrasts)) as.character(contrasts[[1L]]) else "<none>")
            }
            if (identical(step$kind, "run")) {
                event <- Filter(function(e) e$kind %in% c("prot_run", "omic_run"), spy$events)
                source <- if (length(event)) event[[length(event)]]$source else entity(head_before_ref)
                put(k("source_from"), if (length(event)) event[[length(event)]]$kind else "head_before_press")
                put(k("source_label"), label(source))
                ref <- NULL
                for (name in names(snapshot()$states)) if (identical(snapshot()$states[[name]]$entity_id, source)) ref <- snapshot()$states[[name]]
                if (is.null(ref) && identical(entity(head_before_ref), source)) ref <- head_before_ref
                put(k("source_group_moved"), k247GroupOf(manager(), ref, ctx$probe))
                if (length(event) && identical(event[[length(event)]]$kind, "prot_run")) {
                    put(k("prev_input_label"), label(event[[length(event)]]$previous))
                    put(k("normalization_complete_before"), event[[length(event)]]$complete_before)
                    put(k("notice_kit182"), event[[length(event)]]$notice)
                } else {
                    put(k("notice_kit182"), FALSE)
                    if (length(event)) put(k("recorded_start_label"), label(event[[length(event)]]$recorded_start))
                }
                put(k("prereq_warning"), any(grepl("Please complete the Quality Control", notes_text(), fixed = TRUE)))
                hist <- history(); labs <- vapply(entities(), label, "")
                first_run <- which(startsWith(labs, paste0(step$id, ":")))
                at_design <- match(design_state, hist)
                put(k("states_between_design_and_run"), if (length(first_run) && !is.na(at_design)) max(0L, first_run[[1L]] - at_design - 1L) else "<no run state>")
            }
            if (identical(step$kind, "reset")) {
                event <- Filter(function(e) e$kind %in% c("prot_reset", "omic_reset"), spy$events)
                if (length(event)) put(k("recorded_start_label"), label(event[[1L]]$recorded_start))
            }
            if (identical(step$kind, "da")) {
                put(k("result_changed"), !is.null(entity(da_reference())) && !identical(entity(da_reference()), da_before))
                source <- da_source()
                put(k("source_label"), label(entity(source)))
                put(k("group_moved"), k247GroupOf(manager(), source, ctx$probe))
            }
            if (identical(step$kind, "summary")) {
                event <- Filter(function(e) identical(e$kind, "summary_export"), spy$events)
                put(k("refused"), any(grepl("could not be saved", notes_text(), fixed = TRUE)))
                if (length(event)) {
                    put(k("state_used"), event[[1L]]$state)
                    put(k("object_group_moved"), k247ObjectGroupOf(event[[1L]]$manager, event[[1L]]$object, ctx$probe))
                } else put(k("state_used"), "<no integration export started>")
                put(k("integration_files"), integration_files())
            }
            if (identical(step$kind, "record")) {
                record <- tryCatch(buildWorkflowAnalysisRecord(manager()), error = function(e) e)
                if (inherits(record, "error")) {
                    put(k("error"), conditionMessage(record))
                } else {
                    states <- vapply(record$steps, function(s) paste0(s$state, "=", s$operation), "")
                    put(k("steps"), paste(states, collapse = " > "))
                    labs <- vapply(record$steps, function(s) label(s$result_entity_id), "")
                    put(k("step_labels"), paste(labs, collapse = " > "))
                    names_only <- vapply(record$steps, `[[`, "", "state")
                    at_design <- match(design_state, names_only)
                    first_run <- which(grepl("^(run|run2):", labs))
                    put(k("qc_steps_before_run"), if (length(first_run) && !is.na(at_design)) max(0L, first_run[[1L]] - at_design - 1L) else "<no run step>")
                    rows <- record$design$rows
                    group <- k247Try({
                        rows <- as.data.frame(rows, stringsAsFactors = FALSE)
                        hit <- which(vapply(seq_len(nrow(rows)), function(i) any(as.character(unlist(rows[i, ])) %in% ctx$probe), NA))
                        if (!length(hit)) "<absent>" else paste(unique(as.character(rows$group[hit])), collapse = "|")
                    })
                    put(k("design_group_moved"), group)
                    put(k("design_builder_label"), paste(labs[grepl("design", vapply(record$steps, `[[`, "", "operation"))], collapse = ","))
                }
            }
            log$notes <- list()
        }
        put("seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
        shiny::isolate(closeDIAStateAdapter(workflow_data$state_manager))
    })))), error = function(e) e)
    if (inherits(outcome, "error")) {
        put("case_error", conditionMessage(outcome))
        k247Log("case ", lane_id, " ", case_id, " FAILED: ", conditionMessage(outcome))
    }
    invisible(NULL)
}

# ---------------------------------------------------------------------------------------------------------------------
# The presses (ids relative to the omic's module; every button is pressed with its last value plus one, P10).
# ---------------------------------------------------------------------------------------------------------------------

k247Click <- function(api, id) {
    current <- suppressWarnings(as.integer(api$input(id)))
    following <- if (length(current) == 1L && !is.na(current)) current + 1L else 1L
    api$set(stats::setNames(list(following), id))
}

#' The design builder through its inputs (guiEntryDesignActions, with the presses raised and several comparisons).
k247DesignActions <- function(api, design, contrasts, prefix, tech_rep) {
    builder <- function(name) paste0(prefix, name)
    groups <- unique(design$group)
    for (i in seq_along(groups)) {
        api$set(stats::setNames(list(groups[[i]]), builder("new_factor")))
        k247Click(api, builder("add_factor"))
    }
    for (i in seq_along(groups)) {
        api$set(stats::setNames(list(design$sample[design$group == groups[[i]]], groups[[i]], "", "", 1L),
            builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
        k247Click(api, builder("assign_metadata"))
    }
    if (tech_rep && !is.null(design$tech_rep_group)) {
        for (set_name in unique(design$tech_rep_group)) {
            samples <- design$sample[design$tech_rep_group == set_name]
            if (length(samples) < 2L) next
            api$set(stats::setNames(list(samples, "manual", as.integer(design$replicates[match(samples[[1L]], design$sample)])),
                builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
            k247Click(api, builder("assign_tech_reps"))
        }
    }
    for (pair in contrasts) {
        api$set(stats::setNames(list(pair[[1L]], pair[[2L]], "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
        k247Click(api, builder("add_contrast"))
    }
    k247Click(api, builder("save_results"))
}

#' Save Design again: `kind` move (one sample, or trx's one technical-replicate set, to the other group), rename (WT to
#' Ctrl, after the builder's own Reset to Initial State), contrast (one comparison added) or same (unchanged).
k247Resave <- function(kind) function(api, ctx) {
    prot <- identical(ctx$omic, "proteomics")
    prefix <- if (prot) "design_matrix-builder-" else "design-builder-"
    design <- guiEntryReadDesign(api)
    if (prot) design$sample <- e2e_match_available_values(design$sample, api$choices("design_matrix-builder-selected_runs"))
    pair <- ctx$spec$contrast
    contrasts <- list(pair)
    if (identical(kind, "move")) {
        rows <- design$sample %in% ctx$moved
        design$group[rows] <- ctx$spec$new
        if (!is.null(ctx$spec$new_set)) {
            design$tech_rep_group[rows] <- ctx$spec$new_set
            design$replicates[rows] <- ctx$spec$new_replicate
        }
    }
    if (identical(kind, "rename")) {
        design$group[design$group == "WT"] <- "Ctrl"
        contrasts <- list(c("KO", "Ctrl"))
        api$set(stats::setNames(list("all"), paste0(prefix, "reset_scope")))
        k247Click(api, paste0(prefix, "reset_changes"))
        k247Click(api, paste0(prefix, "confirm_reset"))
    }
    if (identical(kind, "contrast")) contrasts <- list(pair, rev(pair))
    k247DesignActions(api, design, contrasts, prefix, tech_rep = !prot)
}

k247Ui <- function(omic, name) {
    if (identical(omic, "transcriptomics")) {
        config <- omicGuiConfig("transcriptomics")
        # As the lane's ui(fn, group): ui(mod_omic_norm_ui, "norm"), ui(mod_omic_da_ui, "da"), ui(mod_omic_summary_ui, "summary").
        fn <- switch(name, norm = mod_omic_norm_ui, da = mod_omic_da_ui, summary = mod_omic_summary_ui)
        return(function(id) fn(id, config[[name]]))
    }
    short <- if (identical(omic, "lipidomics")) "lipid" else "metab"
    get(paste0("mod_", short, "_", name, "_ui"), envir = asNamespace("MultiScholaR"))
}

k247Press <- function(what) function(api, ctx) {
    prot <- identical(ctx$omic, "proteomics")
    trx <- identical(ctx$omic, "transcriptomics")
    if (prot) {
        switch(what,
            reset = k247Click(api, "normalization-reset_normalization"),
            run = {
                if (!identical(ctx$route, "DIA_LIMPA")) api$set(list(`normalization-ruv_mode` = "skip"))
                k247Click(api, "normalization-run_normalization")
            },
            corr = {
                api$set(list(`normalization-norm_qc_tabs` = "Correlation Filtering", `normalization-min_pearson_correlation_threshold` = 0.5))
                k247Click(api, "normalization-apply_correlation_filter")
            },
            ruvcov = k247Click(api, "normalization-apply_ruv_covariates"),
            export = k247Click(api, "normalization-export_filtered_session"),
            daload = k247Click(api, "differential_expression-load_filtered_session"),
            da = {
                api$set(list(`differential_expression-formula_string` = "~ 0 + group", `differential_expression-da_q_val_thresh` = 0.05,
                    `differential_expression-treat_lfc_cutoff` = 0))
                modals <- length(api$modals())
                k247Click(api, "differential_expression-run_da_analysis")
                if (length(api$modals()) > modals) k247Click(api, "differential_expression-acknowledge_qvalue_warning")
            },
            summary = {
                api$defaults(mod_prot_summary_ui, "session_summary")
                api$set(list(`session_summary-experiment_label` = ctx$lane$lane_id, `session_summary-description` = "KIT-247 G0"))
                k247Click(api, "session_summary-save_workflow_args")
            },
            stop("no proteomics press ", what))
        return(invisible(NULL))
    }
    switch(what,
        reset = k247Click(api, "norm-reset_normalization"),
        run = {
            if (trx) {
                api$defaults(k247Ui(ctx$omic, "norm"), "norm")
            } else {
                api$set(list(`norm-apply_itsd` = TRUE, `norm-norm_method` = "cyclicloess", `norm-ruv_mode` = "skip"))
            }
            k247Click(api, "norm-run_normalization")
        },
        corr = {
            api$set(list(`norm-min_pearson_correlation_threshold` = 0.5))
            k247Click(api, "norm-apply_correlation_filter")
        },
        export = k247Click(api, "norm-export_session"),
        daload = k247Click(api, "de-load_filtered_session"),
        da = {
            if (trx) api$defaults(k247Ui(ctx$omic, "da"), "de") else
                api$set(list(`de-formula_string` = "~ 0 + group", `de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 0))
            k247Click(api, "de-run_da_analysis")
        },
        summary = {
            api$defaults(k247Ui(ctx$omic, "summary"), "summary")
            api$set(list(`summary-experiment_label` = ctx$lane$lane_id, `summary-description` = "KIT-247 G0"))
            k247Click(api, "summary-save_workflow_args")
        },
        stop("no press ", what))
    invisible(NULL)
}

#' (d): Accession Cleanup re-applied through its sub-tab (the lane's qc_step, the press raised).
k247QcAccession <- function(api, ctx) {
    specs <- lapply(getProtQcProteinModuleSpecs(ctx$route), function(spec) c(spec, list(group = "protein_qc", tabs = "protein_filter_tabs")))
    hit <- Filter(function(spec) grepl("accession|cleanup", paste(spec$id, spec$title), ignore.case = TRUE), specs)
    if (!length(hit)) stop("no Accession Cleanup sub-tab in ", paste(vapply(specs, `[[`, "", "id"), collapse = ", "))
    spec <- hit[[1L]]
    sub_id <- paste0("quality_control-", spec$group, "-", spec$id)
    ui_fn <- get(spec$uiFn, envir = asNamespace("MultiScholaR"))
    api$set(stats::setNames(list("Protein QC"), "quality_control-qc_tabs_tmt"))
    api$set(stats::setNames(list(spec$title), paste0("quality_control-", spec$group, "-", spec$tabs)))
    settings <- guiEntryProtQcSettings()
    rendered <- names(omicsUiInputDefaults(ui_fn(sub_id)))
    overrides <- settings[paste0(sub_id, "-", names(settings)) %in% rendered]
    if (length(overrides)) api$set(stats::setNames(overrides, paste0(sub_id, "-", names(overrides))))
    actions <- grep("-apply_", guiEntryActionIds(ui_fn(sub_id)), value = TRUE)
    if (!length(actions)) stop("no apply action in ", spec$uiFn)
    k247Click(api, actions[[1L]])
}

k247Step <- function(id, tab, kind, run) list(id = id, tab = tab, kind = kind, run = run)
k247Tabs <- function(omic) {
    if (identical(omic, "proteomics")) list(design = "design", qc = "qc", norm = "normalization", da = "da", summary = "session_summary")
    else list(design = "design", qc = "qc", norm = "norm", da = "de", summary = "summary")
}

#' A case's presses after the prefix, by name.
k247Steps <- function(omic, names, resave = "move") {
    t <- k247Tabs(omic)
    lapply(names, function(name) switch(name,
        resave = k247Step("resave", t$design, "resave", k247Resave(resave)),
        reset = k247Step("reset", t$norm, "reset", k247Press("reset")),
        run = k247Step("run", t$norm, "run", k247Press("run")),
        run2 = k247Step("run2", t$norm, "run", k247Press("run")),
        corr = k247Step("corr", t$norm, "corr", k247Press("corr")),
        ruvcov = k247Step("ruvcov", t$norm, "ruvcov", k247Press("ruvcov")),
        export = k247Step("export", t$norm, "export", k247Press("export")),
        daload = k247Step("daload", t$da, "daload", k247Press("daload")),
        da = k247Step("da", t$da, "da", k247Press("da")),
        summary = k247Step("summary", t$summary, "summary", k247Press("summary")),
        qcre = k247Step("qcre", t$qc, "qc", k247QcAccession),
        record = k247Step("record", NULL, "record", function(api, ctx) invisible(NULL)),
        stop("no step ", name)))
}

# ---------------------------------------------------------------------------------------------------------------------
# The cases.
# ---------------------------------------------------------------------------------------------------------------------

k247Cases <- function(group) {
    cases <- list()
    add <- function(lane, case, until, steps) cases[[length(cases) + 1L]] <<- list(lane = lane, case = case, until = until, steps = steps)
    if (identical(group, "g01")) {
        l <- "prot_lfq"; o <- "proteomics"
        add(l, "main", "correlation", k247Steps(o, c("resave", "reset", "run", "corr", "export", "daload", "da", "record")))
        add(l, "b", "correlation", k247Steps(o, c("resave", "reset", "run", "corr", "export", "daload", "da"), resave = "rename"))
        add(l, "c", "da", k247Steps(o, c("resave", "da")))
        add(l, "d", "correlation", k247Steps(o, c("resave", "qcre", "run", "record")))
        add(l, "e", "correlation", k247Steps(o, c("resave", "run", "record")))
        add(l, "f0", "correlation", k247Steps(o, c("resave", "summary")))
        add(l, "f1", "summary_save", k247Steps(o, c("resave", "summary")))
        add(l, "g", "correlation", k247Steps(o, "resave", resave = "contrast"))
        add(l, "h", "correlation", k247Steps(o, "resave", resave = "same"))
    }
    if (identical(group, "g02")) {
        l <- "prot_dia_limpa"; o <- "proteomics"
        add(l, "main", "ruv_covariates", k247Steps(o, c("resave", "reset", "run", "corr", "ruvcov", "export", "daload", "da", "record")))
        add(l, "e", "ruv_covariates", k247Steps(o, c("resave", "run", "record")))
    }
    if (identical(group, "g03t")) {
        # VAL-54 instrument rerun 1 of 1 (coordinator's ruling): trx's main and (e) again, the group lookup fixed.
        add("trx_screen", "main", "norm", k247Steps("transcriptomics", c("resave", "reset", "run", "da", "record")))
        add("trx_screen", "e", "norm", k247Steps("transcriptomics", c("resave", "run", "record")))
    }
    if (identical(group, "g03")) {
        for (l in c("lipid_canonical", "metab_lc", "trx_screen")) {
            o <- guiEntryLane(l)$omic_type
            trx <- identical(o, "transcriptomics")
            until <- if (trx) "norm" else "correlation"
            main <- if (trx) c("resave", "reset", "run", "da", "record") else
                c("resave", "reset", "run", "corr", "export", "daload", "da", "record")
            add(l, "main", until, k247Steps(o, main))
            add(l, "e", until, k247Steps(o, c("resave", "run", "record")))
            add(l, "f0", until, k247Steps(o, c("resave", "summary")))
            add(l, "f1", "summary_save", k247Steps(o, c("resave", "summary")))
            add(l, "g", until, k247Steps(o, "resave", resave = "contrast"))
            add(l, "h", until, k247Steps(o, "resave", resave = "same"))
        }
    }
    cases
}

k247Write <- function() {
    observed <- if (length(k247Rows$rows)) do.call(rbind, k247Rows$rows) else
        data.frame(lane = character(), case = character(), key = character(), value = character())
    utils::write.table(observed, file.path(k247Out, paste0("g0-", k247Group, "-observed.tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)
    predicted <- k247PredictionTable[k247PredictionTable$lane %in% unique(c(observed$lane, k247GroupLanes)), , drop = FALSE]
    if (identical(k247Group, "g03t")) predicted <- predicted[predicted$case %in% c("main", "e"), , drop = FALSE]
    utils::write.table(predicted, file.path(k247Out, paste0("g0-", k247Group, "-predictions.tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)
    compare <- predicted
    compare$observed <- vapply(seq_len(nrow(compare)), function(i) {
        hit <- observed$value[observed$lane == compare$lane[[i]] & observed$case == compare$case[[i]] & observed$key == compare$key[[i]]]
        if (length(hit)) hit[[length(hit)]] else "<not observed>"
    }, "")
    compare$held <- vapply(seq_len(nrow(compare)), function(i) {
        if (identical(compare$observed[[i]], "<not observed>")) return(NA)
        grepl(compare$pattern[[i]], compare$observed[[i]], perl = TRUE)
    }, NA)
    utils::write.table(compare, file.path(k247Out, paste0("g0-", k247Group, "-compare.tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)
    k247Log("written: ", nrow(observed), " observed rows; ", sum(compare$held %in% TRUE), " predictions held, ",
        sum(compare$held %in% FALSE), " not held, ", sum(is.na(compare$held)), " not observed")
}

k247GroupLanes <- switch(k247Group, g01 = "prot_lfq", g02 = "prot_dia_limpa", g03 = c("lipid_canonical", "metab_lc", "trx_screen"),
    g03t = "trx_screen",
    stop("unknown group ", k247Group))
k247Log("group ", k247Group, ": ", paste(k247GroupLanes, collapse = ", "))
if (!workflowOwnerCapable()) stop("the session owner is not available (qs2/duckdb)")
for (case in k247Cases(k247Group)) {
    tryCatch(k247RunCase(case$lane, case$case, case$until, case$steps), error = function(e) {
        k247Put(case$lane, case$case, "case_error", conditionMessage(e))
    })
    k247Write()
}
k247Write()
k247Log("done")
