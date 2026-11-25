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
kit247DesignActions <- function(api, design, contrasts, prefix, tech_rep) {
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
        api$set(stats::setNames(list(pair[[1L]], pair[[2L]], "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
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
    api$set(stats::setNames(list("design"), api$tabs))
    kit247DesignActions(api, design, contrasts, prefix, tech_rep = !prot)
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
            clear_notes = function() log$notes <- list(),
            modals = function() length(log$modals),
            notices = function() notices$count,
            bus = function(name) shiny::isolate(workflow_data[[name]]),
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
        kit247Click(api, "qc-intensity-revert_filter")
        reverted <- list(history = probe$history(), head = probe$head())
        probe$clear_notes()
        kit247Reset(api, probe$lane)
        list(reverted = reverted, history = probe$history(), head = probe$head(), notes = probe$notes())
    })
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
