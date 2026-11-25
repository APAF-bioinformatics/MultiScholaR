# KIT-247 (VAL-239, the owner's OD-1 (D), OD-2 (a), OD-4 (a)): a design saved again after Quality Control or
# normalisation. Each case drives one omic's module server from import under shiny::testServer, through its namespaced
# inputs only (the GUI-entry driver, helper-gui-entry.R: the lane's own steps up to a named step), then saves the design
# again; the runner is G0's (docs/qa/workflow-execution/kit/kit-247/g0-lanes.R). E22 is a function-level row.
# Texts are drafts (VAL-41; kit-247/wording.md).

kit247Click <- function(api, id) {
    current <- suppressWarnings(as.integer(api$input(id)))
    following <- if (length(current) == 1L && !is.na(current)) current + 1L else 1L
    api$set(stats::setNames(list(following), id))
}

#' The design builder through its inputs (guiEntryDesignActions, presses raised, several comparisons).
kit247DesignActions <- function(api, design, contrasts, prefix, tech_rep, formula = "~ 0 + group") {
    builder <- function(name) paste0(prefix, name)
    groups <- unique(design$group)
    for (i in seq_along(groups)) {
        api$set(stats::setNames(list(groups[[i]]), builder("new_factor")))
        kit247Click(api, builder("add_factor"))
    }
    for (i in seq_along(groups)) {
        api$set(stats::setNames(list(design$sample[design$group == groups[[i]]], groups[[i]], "", "", 1L),
            builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
        kit247Click(api, builder("assign_metadata"))
    }
    if (tech_rep && !is.null(design$tech_rep_group)) {
        for (set_name in unique(design$tech_rep_group)) {
            samples <- design$sample[design$tech_rep_group == set_name]
            if (length(samples) < 2L) next
            api$set(stats::setNames(list(samples, "manual", as.integer(design$replicates[match(samples[[1L]], design$sample)])),
                builder(c("tech_rep_samples", "tech_rep_assignment_mode", "manual_replicate_number"))))
            kit247Click(api, builder("assign_tech_reps"))
        }
    }
    for (pair in contrasts) {
        api$set(stats::setNames(list(pair[[1L]], pair[[2L]], formula), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
        kit247Click(api, builder("add_contrast"))
    }
    kit247Click(api, builder("save_results"))
}

kit247LaneSpec <- function(lane_id) {
    switch(lane_id,
        prot_lfq = list(moved = "KO_5", new = "WT", contrast = c("KO", "WT")),
        prot_tmt = list(moved = NULL, new = "WT", contrast = c("KO", "WT")),
        prot_dia = list(moved = "KO_4", new = "WT", contrast = c("KO", "WT")),
        prot_dia_limpa = list(moved = "KO_4", new = "WT", contrast = c("KO", "WT")),
        lipid_canonical = list(moved = "KO_3", new = "WT", contrast = c("KO", "WT")),
        metab_lc = list(moved = "KO_3", new = "WT", contrast = c("KO", "WT")),
        trx_screen = list(moved = c("Sera_03_r1", "Sera_03_r2"), new = "RPMI", contrast = c("Sera", "RPMI"),
            new_replicate = 4L, new_set = "RPMI_04"))
}

#' Save Design again: kind "move" (one sample, or trx's technical-replicate set, to the other group), "contrast" (one
#' comparison added) or "same".
kit247Resave <- function(api, lane, kind) {
    prot <- identical(lane$omic_type, "proteomics")
    prefix <- if (prot) "design_matrix-builder-" else "design-builder-"
    spec <- kit247LaneSpec(lane$lane_id)
    design <- guiEntryReadDesign(api)
    if (prot) design$sample <- e2e_match_available_values(design$sample, api$choices("design_matrix-builder-selected_runs"))
    contrasts <- list(spec$contrast)
    if (identical(kind, "move")) {
        moved <- if (is.null(spec$moved)) design$sample[design$group != spec$new][[1L]] else spec$moved
        if (prot) moved <- e2e_match_available_values(moved, design$sample)
        rows <- design$sample %in% moved
        design$group[rows] <- spec$new
        if (!is.null(spec$new_set)) {
            design$tech_rep_group[rows] <- spec$new_set
            design$replicates[rows] <- spec$new_replicate
        }
    }
    if (identical(kind, "contrast")) contrasts <- list(spec$contrast, rev(spec$contrast))
    formula <- if (identical(kind, "formula")) "~ 0 + group + replicates" else "~ 0 + group"
    api$set(stats::setNames(list("design"), api$tabs))
    kit247DesignActions(api, design, contrasts, prefix, tech_rep = !prot, formula = formula)
}

kit247Confirm <- function(api, lane) {
    kit247Click(api, if (identical(lane$omic_type, "proteomics")) "design_matrix-confirm_design_resave" else "design-confirm_design_resave")
}

kit247NormTab <- function(lane) if (identical(lane$omic_type, "proteomics")) "normalization" else "norm"

kit247Reset <- function(api, lane) {
    api$set(stats::setNames(list(kit247NormTab(lane)), api$tabs))
    kit247Click(api, if (identical(lane$omic_type, "proteomics")) "normalization-reset_normalization" else "norm-reset_normalization")
}

kit247Run <- function(api, lane) {
    api$set(stats::setNames(list(kit247NormTab(lane)), api$tabs))
    if (identical(lane$omic_type, "proteomics")) {
        if (!isTRUE(lane$use_limpa)) api$set(list(`normalization-ruv_mode` = "skip"))
        kit247Click(api, "normalization-run_normalization")
    } else if (identical(lane$omic_type, "transcriptomics")) {
        kit247Click(api, "norm-run_normalization")
    } else {
        api$set(list(`norm-apply_itsd` = TRUE, `norm-norm_method` = "cyclicloess", `norm-ruv_mode` = "skip"))
        kit247Click(api, "norm-run_normalization")
    }
}

#' One lane up to the lane step `until`, then `actions(api, probe)`; `probe` reads the session. Returns what the actions
#' returned. The KIT-182 notice (`prot_norm_input_current`) is counted in `probe$notices()`.
kit247Lane <- function(lane_id, until, actions) {
    testthat::skip_if_not(workflowOwnerCapable(), "the session owner is not available (qs2/duckdb)")
    lane <- guiEntryLane(lane_id)
    omic <- lane$omic_type
    spec <- guiEntryServerSpec(omic)
    standard <- guiEntrySteps(lane)
    ids <- vapply(standard, `[[`, "", "id")
    prefix_steps <- standard[seq_len(match(until, ids))]
    root <- normalizePath(withr::local_tempdir(pattern = paste0("k247-", lane_id, "-")))
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
    log$updates <- list(); log$notes <- list(); log$modals <- integer(); log$choices <- list()
    guiEntryInstallBrowserDouble(log)
    notices <- new.env(parent = emptyenv())
    notices$count <- 0L
    original_request <- get("diaNormalizationRequest", asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(diaNormalizationRequest = function(manager, settings, norm_data, grouping) {
        withCallingHandlers(original_request(manager, settings, norm_data, grouping),
            prot_norm_input_current = function(m) notices$count <- notices$count + 1L)
    }, .package = "MultiScholaR")
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
            note_types = function() vapply(log$notes, `[[`, "", "type"),
            clear_notes = function() log$notes <- list(),
            modals = function() length(log$modals),
            notices = function() notices$count,
            bus = function(name) shiny::isolate(workflow_data[[name]]),
            steps = stats::setNames(standard, ids), tabs = spec$tabs,
            manager = manager, flush = function() session$flushReact())
        holder$result <- actions(api, probe)
        shiny::isolate(closeDIAStateAdapter(workflow_data$state_manager))
    }))))
    holder$result
}

kit247DesignState <- function(lane) {
    switch(lane$omic_type, proteomics = if (isTRUE(lane$use_limpa) || identical(lane$workflow_type, "DIA")) "raw_data_s4" else "protein_s4_initial",
        lipidomics = "lipid_raw_data_s4", metabolomics = "metab_raw_data_s4", transcriptomics = "trx_raw_data_s4")
}

#' E1/E4/E5/E6/E8/E10/E11 with E12, E16, and E2/E9/E11's Reset: a changed re-save confirmed after the normalisation.
kit247ChangedResave <- function(lane_id, until) {
    kit247Lane(lane_id, until, function(api, probe) {
        lane <- probe$lane
        design <- kit247DesignState(lane)
        out <- list(before_history = probe$history(), modals_before = probe$modals())
        kit247Resave(api, lane, "move")
        out$modals_after_save <- probe$modals()
        kit247Confirm(api, lane)
        probe$flush()
        out$design <- design
        out$history <- probe$history()
        out$head <- probe$head()
        out$status <- probe$status()
        out$notes <- probe$notes()
        probe$clear_notes()
        kit247Reset(api, lane)
        out$reset_history <- probe$history()
        out$reset_head <- probe$head()
        out$reset_notes <- probe$notes()
        out
    })
}

kit247ExpectCut <- function(out, info) {
    testthat::expect_gt(out$modals_after_save, out$modals_before)
    testthat::expect_identical(out$head, out$design, info = info)
    testthat::expect_identical(utils::tail(out$history, 1L), out$design, info = info)
    testthat::expect_identical(length(out$history), match(out$design, out$before_history), info = info)
    testthat::expect_identical(out$status$quality_control, "pending", info = info)
    testthat::expect_identical(out$status$normalization, "disabled", info = info)
    testthat::expect_identical(sum(grepl("The design was saved again", out$notes, fixed = TRUE)), 1L, info = info)
    # Reset after the cut brings nothing back.
    testthat::expect_identical(out$reset_head, out$design, info = info)
    testthat::expect_identical(out$reset_history, out$history, info = info)
}

test_that("KIT-247 E1, E2, E12, E16: label-free, a changed design confirmed over later steps sets them aside; Reset brings nothing back", {
    out <- kit247ChangedResave("prot_lfq", "correlation")
    kit247ExpectCut(out, "prot_lfq")
})

test_that("KIT-247 E3, E14: label-free, a run after the cut reads the new design state, with no 'current quality-control result' notice", {
    out <- kit247Lane("prot_lfq", "correlation", function(api, probe) {
        kit247Resave(api, probe$lane, "move")
        kit247Confirm(api, probe$lane)
        design_entity <- probe$entity("protein_s4_initial")
        notices <- probe$notices()
        kit247Run(api, probe$lane)
        list(history = probe$history(), notices = probe$notices() - notices, design_entity = design_entity,
            design_after = probe$entity("protein_s4_initial"))
    })
    testthat::expect_identical(out$history, c("initial", "protein_s4_initial", "normalized"))
    testthat::expect_identical(out$notices, 0L)
    testthat::expect_identical(out$design_after, out$design_entity)
})

test_that("KIT-247 E19: a changed re-save answered Cancel commits nothing and leaves the bus as it was", {
    out <- kit247Lane("prot_lfq", "correlation", function(api, probe) {
        before <- list(history = probe$history(), head = probe$head(), entity = probe$entity(),
            design = probe$bus("design_matrix"), contrasts = probe$bus("contrasts_tbl"))
        kit247Resave(api, probe$lane, "move")
        list(before = before, history = probe$history(), head = probe$head(), entity = probe$entity(),
            design = probe$bus("design_matrix"), contrasts = probe$bus("contrasts_tbl"))
    })
    testthat::expect_identical(out$history, out$before$history)
    testthat::expect_identical(out$head, out$before$head)
    testthat::expect_identical(out$entity, out$before$entity)
    testthat::expect_identical(out$design, out$before$design)
    testthat::expect_identical(out$contrasts, out$before$contrasts)
})

test_that("KIT-247 E17, E18: label-free, a comparison-only or unchanged re-save keeps every later step", {
    for (kind in c("contrast", "same")) {
        out <- kit247Lane("prot_lfq", "correlation", function(api, probe) {
            before <- list(history = probe$history(), head = probe$head(), entity = probe$entity())
            kit247Resave(api, probe$lane, kind)
            list(before = before, history = probe$history(), head = probe$head(), entity = probe$entity(),
                contrasts = probe$bus("contrasts_tbl"), status = probe$status(), notes = probe$notes(),
                pending = probe$manager()$.context$pending_provenance)
        })
        testthat::expect_identical(out$history, out$before$history, info = kind)
        testthat::expect_identical(out$head, out$before$head, info = kind)
        testthat::expect_identical(out$entity, out$before$entity, info = kind)
        testthat::expect_null(out$pending, info = kind)
        if (identical(kind, "contrast")) {
            testthat::expect_identical(nrow(out$contrasts), 2L)
            testthat::expect_identical(out$status$differential_expression, "pending")
            testthat::expect_true(any(grepl("Comparisons saved", out$notes, fixed = TRUE)))
        } else {
            testthat::expect_true(any(grepl("Design saved; nothing changed.", out$notes, fixed = TRUE)))
        }
    }
})

test_that("KIT-247 E21 (OD-4 (a)): the cut removes the omic's earlier integration input, folder and top-level file", {
    out <- kit247Lane("prot_lfq", "summary_save", function(api, probe) {
        integration <- probe$paths$integration_dir %||% file.path(probe$paths$base_dir, "integration")
        listing <- function() list.files(integration, recursive = TRUE)
        before <- listing()
        kit247Resave(api, probe$lane, "move")
        kit247Confirm(api, probe$lane)
        list(before = before, after = listing(), notes = probe$notes())
    })
    testthat::expect_true(length(out$before) > 0L)
    testthat::expect_identical(length(grep("^proteomics", out$after)), 0L)
    testthat::expect_true(any(grepl("integration input of this omic was removed", out$notes, fixed = TRUE)))
})

test_that("KIT-247 E4, E5, E6: TMT, DIA (iq) and DIA with limpa: the cut, and Reset brings nothing back", {
    for (lane_id in c("prot_tmt", "prot_dia", "prot_dia_limpa")) {
        until <- if (identical(lane_id, "prot_dia_limpa")) "ruv_covariates" else "correlation"
        if (identical(lane_id, "prot_dia_limpa")) testthat::skip_if_not_installed("limpa")
        kit247ExpectCut(kit247ChangedResave(lane_id, until), lane_id)
    }
})

test_that("KIT-247 E8, E9, E10, E11, E15: lipidomics, metabolomics, transcriptomics: the cut; Reset brings nothing back (OD-2 (a)); a run reads only the new design state", {
    for (lane_id in c("lipid_canonical", "metab_lc", "trx_screen")) {
        until <- if (identical(lane_id, "trx_screen")) "norm" else "correlation"
        out <- kit247ChangedResave(lane_id, until)
        kit247ExpectCut(out, lane_id)
        testthat::expect_false(any(grepl("_reset$", out$reset_history)), info = lane_id)
        testthat::expect_true(any(grepl("no longer current", out$reset_notes, fixed = TRUE)), info = lane_id)
    }
    out <- kit247Lane("lipid_canonical", "correlation", function(api, probe) {
        kit247Resave(api, probe$lane, "move")
        kit247Confirm(api, probe$lane)
        kit247Run(api, probe$lane)
        probe$history()
    })
    at <- match("lipid_raw_data_s4", out)
    testthat::expect_false(any(grepl("intensity|duplicate|qc|filtered_qc", out[-seq_len(at)])))
})

test_that("KIT-247 E17 (lipidomics): a comparison-only re-save keeps every later step", {
    out <- kit247Lane("lipid_canonical", "correlation", function(api, probe) {
        before <- list(history = probe$history(), head = probe$head(), entity = probe$entity())
        kit247Resave(api, probe$lane, "contrast")
        list(before = before, history = probe$history(), head = probe$head(), entity = probe$entity(),
            contrasts = probe$bus("contrasts_tbl"))
    })
    testthat::expect_identical(out$history, out$before$history)
    testthat::expect_identical(out$head, out$before$head)
    testthat::expect_identical(out$entity, out$before$entity)
    testthat::expect_identical(nrow(out$contrasts), 2L)
})

test_that("KIT-247 G2: a design saved twice before any Quality Control is committed as before, with no question and no cut", {
    out <- kit247Lane("lipid_canonical", "design", function(api, probe) {
        modals <- probe$modals()
        kit247Resave(api, probe$lane, "move")
        list(history = probe$history(), head = probe$head(), modals = probe$modals() - modals, notes = probe$notes())
    })
    testthat::expect_identical(out$head, "lipid_raw_data_s4")
    testthat::expect_identical(out$history, c("initial", "lipid_raw_data_s4"))
    testthat::expect_false(any(grepl("saved again", out$notes, fixed = TRUE)))
})

test_that("KIT-247 E22 (R1-03): a re-save that changes only a batch or factor value is 'changed'", {
    manager <- structure(new.env(parent = emptyenv()), class = "DIAStateAdapter")
    manager$.context <- new.env(parent = emptyenv())
    manager$.context$snapshot <- list(state_history = list("initial", "lipid_raw_data_s4"))
    design <- data.frame(Run = c("A", "B", "C", "D"), group = c("x", "x", "y", "y"), replicates = c(1L, 2L, 1L, 2L),
        factor1 = c("x", "x", "y", "y"), factor2 = NA_character_, batch = c("b1", "b2", "b1", "b2"), stringsAsFactors = FALSE)
    data <- list(assay = data.frame(id = 1:3, A = 1:3, B = 4:6, C = 7:9, D = 10:12))
    batch <- design
    batch$batch[[1L]] <- "b2"
    factor <- design
    factor$factor2[[3L]] <- "z"
    current <- list(design_matrix = design, data = data)
    testthat::expect_identical(compareOmicsDesignResult(manager, list(design_matrix = design, data = data), "lipid_raw_data_s4", current), "unchanged")
    testthat::expect_identical(compareOmicsDesignResult(manager, list(design_matrix = batch, data = data), "lipid_raw_data_s4", current), "changed")
    testthat::expect_identical(compareOmicsDesignResult(manager, list(design_matrix = factor, data = data), "lipid_raw_data_s4", current), "changed")
    testthat::expect_identical(compareOmicsDesignResult(manager, list(design_matrix = design, data = data), "metab_raw_data_s4", current), "first")
})

test_that("KIT-247 E13 (OD-2 (a)): lipidomics, Reset after a Quality Control Revert below the run's start brings nothing back", {
    out <- kit247Lane("lipid_canonical", "correlation", function(api, probe) {
        api$set(stats::setNames(list("qc"), api$tabs))
        # Each Revert selects the state before the head (registerOmicQcRevertObserver): pressed until the head is the
        # design state, below the run's start (R1 run 45604461406: one press only reached lipid_normalized).
        for (i in seq_len(15L)) {
            if (identical(probe$head(), "lipid_raw_data_s4")) break
            kit247Click(api, "qc-intensity-revert_filter")
        }
        reverted <- list(history = probe$history(), head = probe$head())
        probe$clear_notes()
        kit247Reset(api, probe$lane)
        list(reverted = reverted, history = probe$history(), head = probe$head(), notes = probe$notes())
    })
    testthat::expect_identical(out$reverted$head, "lipid_raw_data_s4")
    testthat::expect_false(any(grepl("_reset$", out$history)))
    testthat::expect_identical(out$head, out$reverted$head)
    testthat::expect_true(any(grepl("no longer current", out$notes, fixed = TRUE)))
})

test_that("KIT-247 E23 (VAL-239 (iv), VAL-148's rule completed): a refused Summary Save removes the omic's folder AND its top-level file", {
    for (lane_id in c("lipid_canonical", "metab_lc", "trx_screen")) {
        out <- kit247Lane(lane_id, "summary_save", function(api, probe) {
            integration <- probe$paths$integration_dir %||% file.path(probe$paths$base_dir, "integration")
            omic <- probe$lane$omic_type
            mine <- function() grep(paste0("^", omic), list.files(integration, recursive = TRUE), value = TRUE)
            before <- mine()
            kit247Reset(api, probe$lane)
            api$set(stats::setNames(list("summary"), api$tabs))
            kit247Click(api, "summary-save_workflow_args")
            list(before = before, after = mine(), notes = probe$notes())
        })
        testthat::expect_true(length(out$before) > 1L, info = lane_id)
        testthat::expect_true(any(grepl("could not be saved", out$notes, fixed = TRUE)), info = lane_id)
        testthat::expect_identical(out$after, character(), info = lane_id)
    }
})

# ---------------------------------------------------------------------------------------------------------------------
# The second correction (owner's VAL-277; AUDIT-KIT-247-impl-20261126 I01, I02, I07, G3).
# ---------------------------------------------------------------------------------------------------------------------

#' The lane's exported design folder (the project's source_dir after its design save), copied, with `edit` applied to
#' its design_matrix.tab.
kit247ExportedDesign <- function(probe, edit) {
    source <- probe$paths$source_dir
    target <- file.path(tempfile("k247-export-"))
    dir.create(target)
    files <- list.files(source, pattern = "^(design_matrix|data_cln.*|contrast_strings|assay_manifest|column_mapping|config)\\.(tab|txt|json|ini)$")
    file.copy(file.path(source, files), target)
    design <- utils::read.delim(file.path(target, "design_matrix.tab"), stringsAsFactors = FALSE, check.names = FALSE)
    design <- edit(design)
    utils::write.table(design, file.path(target, "design_matrix.tab"), sep = "\t", row.names = FALSE, quote = FALSE)
    target
}

#' One run moved to the other group (a new replicate number after the group's last).
kit247MoveRun <- function(design, run, to) {
    row <- which(apply(design, 1L, function(r) any(as.character(r) == run)))[[1L]]
    design$group[[row]] <- to
    if ("factor1" %in% names(design)) design$factor1[[row]] <- to
    design$replicates[[row]] <- max(design$replicates[design$group == to], na.rm = TRUE) + 1L
    if ("tech_rep_group" %in% names(design)) design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    design
}

#' A KIT-222 refusal: one run given the replicate number of another of its group, neither marked.
kit247Collide <- function(design, run) {
    row <- which(apply(design, 1L, function(r) any(as.character(r) == run)))[[1L]]
    other <- which(design$group == design$group[[row]] & seq_len(nrow(design)) != row)[[1L]]
    design$replicates[[row]] <- design$replicates[[other]]
    if (!"tech_reps" %in% names(design)) design$tech_reps <- NA_integer_
    design$tech_reps[c(row, other)] <- NA_integer_
    if ("tech_rep_group" %in% names(design)) design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    design
}

kit247Import <- function(api, lane, folder) {
    api$set(stats::setNames(list("design"), api$tabs))
    if (identical(lane$omic_type, "proteomics")) {
        api$set(list(`design_matrix-import_dir_text` = folder, `design_matrix-import_taxon_id` = 9606L,
            `design_matrix-import_organism_name` = "Homo sapiens"))
        kit247Click(api, "design_matrix-confirm_import")
    } else {
        api$set(list(`design-import_dir` = folder))
        kit247Click(api, "design-confirm_import")
    }
}

kit247MovedRun <- function(api, lane) {
    spec <- kit247LaneSpec(lane$lane_id)
    if (identical(lane$omic_type, "proteomics")) {
        e2e_match_available_values(spec$moved, api$choices("design_matrix-builder-selected_runs"))
    } else spec$moved[[1L]]
}

for (lane_id in c("prot_lfq", "lipid_canonical")) {
    test_that(paste("KIT-247 E7 (VAL-277, AUDIT I01):", lane_id, "Import Existing Design of a changed design asks; Save cuts; Cancel changes nothing"), {
        state <- if (identical(lane_id, "prot_lfq")) "protein_s4_initial" else "lipid_raw_data_s4"
        out <- kit247Lane(lane_id, "correlation", function(api, probe) {
            run <- kit247MovedRun(api, probe$lane)
            to <- kit247LaneSpec(lane_id)$new
            folder <- kit247ExportedDesign(probe, function(d) kit247MoveRun(d, run, to))
            before <- list(history = probe$history(), head = probe$head(), design = probe$bus("design_matrix"),
                config = probe$bus("config_list"), modals = probe$modals())
            kit247Import(api, probe$lane, folder)
            cancelled <- list(history = probe$history(), head = probe$head(), design = probe$bus("design_matrix"),
                config = probe$bus("config_list"), modals = probe$modals())
            kit247Import(api, probe$lane, folder)
            if (identical(probe$lane$omic_type, "proteomics")) kit247Click(api, "design_matrix-confirm_design_resave") else
                kit247Click(api, "design-confirm_design_resave")
            probe$flush()
            list(before = before, cancelled = cancelled, history = probe$history(), head = probe$head(),
                status = probe$status(), notes = probe$notes())
        })
        expect_gt(out$cancelled$modals, out$before$modals)
        expect_identical(out$cancelled$history, out$before$history)
        expect_identical(out$cancelled$head, out$before$head)
        expect_identical(out$cancelled$design, out$before$design)
        expect_identical(out$cancelled$config, out$before$config)
        expect_identical(out$head, state)
        expect_identical(utils::tail(out$history, 1L), state)
        expect_identical(out$status$quality_control, "pending")
        expect_true(any(grepl("The design was saved again", out$notes, fixed = TRUE)))
    })

    test_that(paste("KIT-247 E7 error (VAL-277, AUDIT I01):", lane_id, "an import refused after the confirmation is shown, the session lives on, nothing changes"), {
        out <- kit247Lane(lane_id, "correlation", function(api, probe) {
            run <- kit247MovedRun(api, probe$lane)
            folder <- kit247ExportedDesign(probe, function(d) kit247Collide(d, run))
            before <- list(history = probe$history(), head = probe$head(), design = probe$bus("design_matrix"))
            kit247Import(api, probe$lane, folder)
            if (identical(probe$lane$omic_type, "proteomics")) kit247Click(api, "design_matrix-confirm_design_resave") else
                kit247Click(api, "design-confirm_design_resave")
            probe$flush()
            list(before = before, history = probe$history(), head = probe$head(), design = probe$bus("design_matrix"),
                notes = probe$notes(), alive = TRUE)
        })
        expect_true(isTRUE(out$alive))
        expect_true(any(grepl("Error during import", out$notes, fixed = TRUE)))
        expect_identical(out$history, out$before$history)
        expect_identical(out$head, out$before$head)
        expect_identical(out$design, out$before$design)
    })
}

test_that("KIT-247 I02 (VAL-277): a design commit that fails after the bus write, then the same design saved again, is not judged unchanged", {
    out <- kit247Lane("lipid_canonical", "correlation", function(api, probe) {
        # VAL-54 instrument fix (run 45676491133): the planted failure is scoped to these two presses; before, it outlived
        # the row and every later lipidomics lane's design save failed with it.
        testthat::with_mocked_bindings({
            kit247Resave(api, probe$lane, "move")
            kit247Click(api, "design-confirm_design_resave")
        }, writeOmicsDesignArtifacts = function(...) stop("KIT-247 planted failure after the bus write"), .package = "MultiScholaR")
        stopifnot(!any(grepl("planted failure", deparse(get("writeOmicsDesignArtifacts", asNamespace("MultiScholaR"))))))
        failed <- list(history = probe$history(), notes = probe$notes())
        probe$clear_notes()
        kit247Resave(api, probe$lane, "move")
        kit247Click(api, "design-confirm_design_resave")
        probe$flush()
        list(failed = failed, history = probe$history(), head = probe$head(), notes = probe$notes())
    })
    expect_true(any(grepl("Error saving design matrix results", out$failed$notes, fixed = TRUE)))
    expect_false(any(grepl("nothing changed|Comparisons saved", out$notes)))
    expect_identical(out$head, "lipid_raw_data_s4")
    expect_identical(utils::tail(out$history, 1L), "lipid_raw_data_s4")
})

#' The group of `runs` in the design of the DA the store holds (its source table's design component).
kit247DaGroup <- function(manager, runs) {
    withWorkflowStateScope(manager, function() {
        reference <- manager$.context$snapshot$da_results
        if (!isSavedWorkflowReference(reference)) return("<no DA>")
        source <- resolveWorkflowStateValue(manager$.context, reference)$source
        component <- readComponentMetadata(manager, source)$data$components$design_matrix
        design <- as.data.frame(materializeWorkflowTable(manager$.context, component), stringsAsFactors = FALSE)
        hit <- which(apply(design, 1L, function(r) any(as.character(r) %in% runs)))
        paste(unique(as.character(design$group[hit])), collapse = "|")
    })
}

#' I07's sequence after the confirmed re-save: Reset, (Quality Control re-applied when `reapply`), run, correlation,
#' export, load, DA; every notification and the active history's state names are recorded after each press, and whether
#' the normalisation run committed a state (coordinator's ruling: recording only; I07's presses and expectations kept).
kit247I07 <- function(lane_id, reapply = FALSE) {
    kit247Lane(lane_id, "correlation", function(api, probe) {
        lane <- probe$lane
        log <- list()
        note <- function(press) {
            log[[length(log) + 1L]] <<- list(press = press, history = paste(probe$history(), collapse = " > "),
                notes = paste(probe$notes(), collapse = " || "))
            probe$clear_notes()
        }
        run <- kit247MovedRun(api, lane)
        kit247Resave(api, lane, "move")
        kit247Click(api, if (identical(lane$omic_type, "proteomics")) "design_matrix-confirm_design_resave" else "design-confirm_design_resave")
        note("resave + confirm")
        kit247Reset(api, lane)
        note("reset")
        if (reapply) {
            for (id in c("qc_intensity", "qc_duplicates", "qc_finalize")) {
                api$set(stats::setNames(list(probe$steps[[id]]$tab), probe$tabs))
                probe$steps[[id]]$run(api)
                note(id)
            }
        }
        before_run <- probe$history()
        kit247Run(api, lane)
        after_run <- probe$history()
        run_committed <- length(after_run) > length(before_run) || !identical(after_run, before_run)
        note("run")
        prot <- identical(lane$omic_type, "proteomics")
        api$set(stats::setNames(list(if (prot) "normalization" else "norm"), api$tabs))
        if (prot) {
            api$set(list(`normalization-norm_qc_tabs` = "Correlation Filtering", `normalization-min_pearson_correlation_threshold` = 0.5))
            kit247Click(api, "normalization-apply_correlation_filter")
            note("correlation")
            kit247Click(api, "normalization-export_filtered_session")
            note("export")
            api$set(stats::setNames(list("da"), api$tabs))
            kit247Click(api, "differential_expression-load_filtered_session")
            note("load")
            api$set(list(`differential_expression-formula_string` = "~ 0 + group", `differential_expression-da_q_val_thresh` = 0.05,
                `differential_expression-treat_lfc_cutoff` = 0))
            modals <- probe$modals()
            kit247Click(api, "differential_expression-run_da_analysis")
            if (probe$modals() > modals) kit247Click(api, "differential_expression-acknowledge_qvalue_warning")
            note("da")
        } else {
            api$set(list(`norm-min_pearson_correlation_threshold` = 0.5))
            kit247Click(api, "norm-apply_correlation_filter")
            note("correlation")
            kit247Click(api, "norm-export_session")
            note("export")
            api$set(stats::setNames(list("de"), api$tabs))
            kit247Click(api, "de-load_filtered_session")
            note("load")
            api$set(list(`de-formula_string` = "~ 0 + group", `de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 0))
            kit247Click(api, "de-run_da_analysis")
            note("da")
        }
        probe$flush()
        record <- buildWorkflowAnalysisRecord(probe$manager())
        rows <- as.data.frame(record$design$rows, stringsAsFactors = FALSE)
        hit <- which(apply(rows, 1L, function(r) any(as.character(r) == run)))
        builder <- Filter(function(step) identical(step$operation, "design_builder"), record$steps)
        list(da = kit247DaGroup(probe$manager(), run), record_group = unique(as.character(rows$group[hit])),
            builder_entity = if (length(builder)) builder[[1L]]$result_entity_id else NA_character_,
            design_entity = probe$entity(kit247DesignState(lane)), run_committed = run_committed, log = log,
            record_steps = paste(vapply(record$steps, function(step) paste0(step$state, "=", step$operation), ""), collapse = " > "))
    })
}

kit247PrintLog <- function(label, out) {
    cat(sprintf("KIT247-LOG %s run_committed=%s da=%s record_group=%s builder_entity=%s design_entity=%s\n", label,
        out$run_committed, out$da, paste(out$record_group, collapse = "|"), out$builder_entity, out$design_entity))
    cat(sprintf("KIT247-LOG %s record_steps: %s\n", label, out$record_steps))
    for (entry in out$log) cat(sprintf("KIT247-LOG %s [%s] history: %s | notes: %s\n", label, entry$press, entry$history, entry$notes))
}

for (lane_id in c("prot_lfq", "lipid_canonical")) {
    test_that(paste("KIT-247 I07 (VAL-277):", lane_id, "after the confirmed re-save, Reset and a run, the DA is on the NEW groups and the record agrees with itself"), {
        out <- kit247I07(lane_id)
        kit247PrintLog(paste("I07", lane_id), out)
        new <- kit247LaneSpec(lane_id)$new
        expect_identical(out$da, new)
        expect_identical(out$record_group, new)
        expect_identical(out$builder_entity, out$design_entity)
    })
}

test_that("KIT-247 I07b (VAL-277, coordinator's ruling): lipidomics, Quality Control re-applied after the confirmed re-save, then the run: the DA is on the NEW groups and the record agrees with itself", {
    out <- kit247I07("lipid_canonical", reapply = TRUE)
    kit247PrintLog("I07b lipid_canonical", out)
    expect_identical(out$da, "WT")
    expect_identical(out$record_group, "WT")
    expect_false(is.na(out$builder_entity))
    expect_identical(out$builder_entity, out$design_entity)
})

test_that("KIT-247 E15b (coordinator's ruling): after the cut, a lipidomics run with no Quality Control either commits a state or says why it did not; never neither", {
    out <- kit247Lane("lipid_canonical", "correlation", function(api, probe) {
        kit247Resave(api, probe$lane, "move")
        kit247Click(api, "design-confirm_design_resave")
        before <- probe$history()
        probe$clear_notes()
        kit247Run(api, probe$lane)
        list(before = before, after = probe$history(), notes = probe$notes(), types = probe$note_types())
    })
    cat(sprintf("KIT247-LOG E15b before: %s | after: %s | notes: %s\n", paste(out$before, collapse = " > "),
        paste(out$after, collapse = " > "), paste(out$notes, collapse = " || ")))
    committed <- !identical(out$after, out$before)
    explained <- any(out$types %in% c("error", "warning"))
    expect_true(committed || explained)
})

test_that("KIT-247 G3 (VAL-277): after the cut, Quality Control re-applied and a run give the numbers of a project designed that way from the start", {
    qc <- c("qc_intensity", "qc_duplicates", "qc_finalize")
    normalised <- function(probe) guiEntryDigest(withWorkflowStateScope(probe$manager(), function() probe$manager()$getState("lipid_normalized")))
    redone <- kit247Lane("lipid_canonical", "correlation", function(api, probe) {
        kit247Resave(api, probe$lane, "move")
        kit247Click(api, "design-confirm_design_resave")
        for (id in qc) {
            api$set(stats::setNames(list(probe$steps[[id]]$tab), probe$tabs))
            probe$steps[[id]]$run(api)
        }
        kit247Run(api, probe$lane)
        normalised(probe)
    })
    fresh <- kit247Lane("lipid_canonical", "design_refusal", function(api, probe) {
        kit247Resave(api, probe$lane, "move")
        for (id in c(qc, "norm")) {
            api$set(stats::setNames(list(probe$steps[[id]]$tab), probe$tabs))
            probe$steps[[id]]$run(api)
        }
        normalised(probe)
    })
    expect_identical(redone, fresh)
})

# ---------------------------------------------------------------------------------------------------------------------
# VAL-331 (the owner, 2026-11-26): VAL-288 is WITHDRAWN from KIT-247. A formula-only Save and, on the proteomics import,
# an organism, taxon or FASTA-only Save answer "nothing changed": no question, nothing reset, no status moves.
# ---------------------------------------------------------------------------------------------------------------------

for (lane_id in c("prot_lfq", "lipid_canonical")) {
    test_that(paste("KIT-247 W288a (VAL-331):", lane_id, "a formula-only re-save asks nothing, resets nothing and says nothing changed"), {
        key <- if (identical(lane_id, "prot_lfq")) "differential_expression" else "differential_analysis"
        out <- kit247Lane(lane_id, "da", function(api, probe) {
            before <- list(history = probe$history(), head = probe$head(), entity = probe$entity(), status = probe$status(),
                modals = probe$modals())
            probe$clear_notes()
            kit247Resave(api, probe$lane, "formula")
            list(before = before, modals = probe$modals(), history = probe$history(), head = probe$head(),
                entity = probe$entity(), status = probe$status(), notes = probe$notes())
        })
        expect_identical(out$before$status[[key]], "complete")
        expect_identical(out$modals, out$before$modals)
        expect_true(any(grepl("Design saved; nothing changed.", out$notes, fixed = TRUE)))
        expect_identical(out$history, out$before$history)
        expect_identical(out$head, out$before$head)
        expect_identical(out$entity, out$before$entity)
        expect_identical(out$status, out$before$status)
    })
}

test_that("KIT-247 W288b (VAL-331): a proteomics Import Existing Design with only the organism changed asks nothing, applies nothing and says nothing changed", {
    out <- kit247Lane("prot_lfq", "correlation", function(api, probe) {
        folder <- kit247ExportedDesign(probe, identity)
        before <- list(history = probe$history(), head = probe$head(), entity = probe$entity(), modals = probe$modals(),
            status = probe$status(), taxon = probe$bus("taxon_id"), organism = probe$bus("organism_name"))
        probe$clear_notes()
        api$set(stats::setNames(list("design"), api$tabs))
        api$set(list(`design_matrix-import_dir_text` = folder, `design_matrix-import_taxon_id` = 10090L,
            `design_matrix-import_organism_name` = "Mus musculus"))
        kit247Click(api, "design_matrix-confirm_import")
        probe$flush()
        list(before = before, modals = probe$modals(), history = probe$history(), head = probe$head(),
            entity = probe$entity(), status = probe$status(), taxon = probe$bus("taxon_id"),
            organism = probe$bus("organism_name"), notes = probe$notes())
    })
    expect_identical(out$modals, out$before$modals)
    expect_true(any(grepl("Design saved; nothing changed.", out$notes, fixed = TRUE)))
    expect_identical(out$history, out$before$history)
    expect_identical(out$head, out$before$head)
    expect_identical(out$entity, out$before$entity)
    expect_identical(out$status, out$before$status)
    expect_identical(out$taxon, out$before$taxon)
    expect_identical(out$organism, out$before$organism)
})
