# KIT-103 W4: the shared design tab (R/mod_omic_design.R) and builder (R/mod_omic_design_builder.R),
# with lipid and trx on them through the lipid wrappers. Each block pins today's behaviour: the
# frozen values (UI digests, output names, notification texts and their order) were taken from the
# lipid design files at 616e1147, and the wrapper blocks pass there unchanged
# (docs/qa/workflow-execution/kit/kit-103/w4/test-first.txt). KIT-103 W5 adds metab on them (the
# frozen values from the metab design files at de2d891c; w5/test-first.txt) and KIT-101 F-01's two
# tests, retargeted from test-metab-04 to the shared commit with identical expectations.

# The rendered HTML's sha256, with shiny's random tabset ids masked (as the entry suite does).
omicDesignUiDigest <- function(tag) {
    html <- as.character(htmltools::renderTags(tag)$html)
    html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"<ID>\"", html)
    html <- gsub("tab-[0-9]+-([0-9]+)", "tab-<ID>-\\1", html)
    digest::digest(html, algo = "sha256", serialize = FALSE)
}

# The modal, notification and removal calls of a run, in order (shiny's functions mocked in the caller).
omicDesignCaptureEvents <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$events <- character()
    add <- function(event) seen$events <- c(seen$events, event)
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default", duration = 5, id = NULL) {
            add(sprintf("note[%s|%s|%s]: %s", type, if (is.null(duration)) "NULL" else duration, id %||% "",
                paste(as.character(ui), collapse = " ")))
            invisible(NULL)
        },
        removeNotification = function(id, ...) add(paste("remove_note:", id)),
        showModal = function(ui, ...) add("modal"),
        removeModal = function(...) add("remove_modal"),
        .package = "shiny", .env = env)
    seen
}

omicDesignFrozen <- list(
    lipid_frame = "e988e3400e2b9882fb67bd5c69054034e17c7857c9dac965b936de9bc5a7d44a",
    trx_frame = "32e9a8da7eaa7bdf516456f77bce804b1f3b7967e8f2f6e78e62b398c4ef2b07",
    builder = "6a09783df4d081979c18d5156c6775c1ac6240c135209a8cbabfe0dc9df47588",
    modal_test = "b1acb5e622d9c60d2edc7acc66d0d2255129529a717b1273c3d05f2323df7b2f",
    modal_app = "2f465045936f7f8796cf0d16d34b9af52e5a7f6015b9e41375627529f343b3bb",
    reset_modal = "f574bf2d9c160d7d1771a6585a794a64cdf82ab7b3383bc0c6c23e4edb564fc4",
    outputs = c("assays_preview", "builder-available_factors_display", "builder-contrast_factors_info",
        "builder-data_table", "builder-defined_contrasts_display", "builder-range_preview",
        "builder-removed_samples_display", "builder-replicate_inputs", "builder-tech_rep_summary",
        "contrasts_preview", "data_available", "design_matrix_exists", "design_matrix_preview", "import_dir_path"),
    # KIT-103 W5: the metab design files at de2d891c.
    metab_frame = "7f8b554fa4f786d6fd2ee6e10d8c83ad62ff68c1bdc19157fc30b27c20e2516c",
    metab_builder = "fec6c68a518969a96e8037b92c1c32a050582271e8433896c3097d589a17bc65",
    metab_modal_test = "817486e426077ff3c43e02e96e0dd08d7c667a9c6ba57a0226eb1cd62c4978e3",
    metab_modal_app = "40d0cadbecf809cc8c7cd64dce3547d92f2c66c6fb17fa3eca401c9834907b71")

test_that("W4 the design tab's UI through the lipid wrappers is today's bytes (lipid, trx; frame and builder)", {
    expect_identical(omicDesignUiDigest(mod_lipid_design_ui("design")), omicDesignFrozen$lipid_frame)
    expect_identical(omicDesignUiDigest(mod_lipid_design_ui("design", profile = trxOmicProfile())),
        omicDesignFrozen$trx_frame)
    expect_identical(omicDesignUiDigest(mod_lipid_design_builder_ui("builder")), omicDesignFrozen$builder)
})

test_that("W4 the shared frame, builder, import modal and reset modal render today's bytes from the configuration", {
    expect_identical(omicDesignUiDigest(mod_omic_design_ui("design", omicGuiConfig("lipidomics"))),
        omicDesignFrozen$lipid_frame)
    expect_identical(omicDesignUiDigest(mod_omic_design_ui("design", omicGuiConfig("transcriptomics"))),
        omicDesignFrozen$trx_frame)
    expect_identical(omicDesignUiDigest(mod_omic_design_builder_ui("builder", omicGuiConfig("lipidomics"))),
        omicDesignFrozen$builder)
    # Trx renders the lipid builder (its builder texts are lipid's).
    expect_identical(omicDesignUiDigest(mod_omic_design_builder_ui("builder", omicGuiConfig("transcriptomics"))),
        omicDesignFrozen$builder)
    ns <- shiny::NS("design")
    config <- omicGuiConfig("lipidomics")
    expect_identical(omicDesignUiDigest(buildOmicsDesignImportModal(ns, config, testMode = TRUE)),
        omicDesignFrozen$modal_test)
    expect_identical(omicDesignUiDigest(buildOmicsDesignImportModal(ns, config, testMode = FALSE)),
        omicDesignFrozen$modal_app)
    expect_identical(omicDesignUiDigest(omicsDesignResetModal("all", shiny::NS("design-builder"))),
        omicDesignFrozen$reset_modal)
})

test_that("W4 the design tab's server keeps today's outputs; the two status outputs read as today's error", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    for (profile in list(NULL, trxOmicProfile())) {
        payload <- module_ci_lipid_design_payload(layout = "lc_pair")
        workflow <- module_ci_lipid_design_workflow(payload)
        shiny::testServer(mod_lipid_design_server, args = list(workflow_data = workflow,
            experiment_paths = module_ci_lipid_design_paths(), profile = profile), {
            session$flushReact()
            outputs <- sort(sub("^proxy[0-9]+-", "", names(session$.__enclos_env__$private$outs)))
            expect_identical(outputs, omicDesignFrozen$outputs)
            for (name in c("data_available", "design_matrix_exists")) {
                error <- tryCatch(output[[name]], error = identity)
                expect_s3_class(error, "error")
                expect_identical(conditionMessage(error), "unused arguments (self, name)")
            }
            expect_identical(output$assays_preview, "Included assays: LCMS_Pos, LCMS_Neg")
            expect_identical(output$`builder-removed_samples_display`, "No samples have been removed.")
            expect_identical(output$`builder-tech_rep_summary`, "No technical replicates assigned yet.")
        })
    }
})

test_that("W4 the builder keeps its reactive `design_matrix` and no name that masks a test's own variables", {
    payload <- module_ci_lipid_design_payload(layout = "lc_pair")
    # The variables the MCI-027 and W2 builder sessions use inside testServer.
    test_names <- c("design", "levels", "runs", "group", "clicks", "tech", "set", "remove", "contrast",
        "technical", "out", "workflow", "fixture", "i", "result", "samples", "assays", "mapping", "seen")
    # The fresh design's run order: lipid mixedsort, trx the input order.
    cases <- list(list(profile = NULL, runs = c("KO_1", "KO_2", "WT_1", "WT_2")),
        list(profile = trxOmicProfile(), runs = c("WT_1", "WT_2", "KO_1", "KO_2")))
    for (case in cases) {
        shiny::testServer(mod_lipid_design_builder_server, args = list(data_tbl = function() payload$data_cln,
            config_list = shiny::reactive(payload$config_list), column_mapping = shiny::reactive(payload$column_mapping),
            profile = case$profile), {
            module_names <- c(ls(session$env, all.names = TRUE), ls(parent.env(session$env), all.names = TRUE))
            expect_true("design_matrix" %in% ls(session$env))
            expect_s3_class(design_matrix, "reactiveVal")
            expect_identical(intersect(module_names, test_names), character())
            session$flushReact()
            expect_identical(shiny::isolate(design_matrix()$Run), case$runs)
        })
    }
})

test_that("W4 the builder's texts, refusal and save are today's (lipid)", {
    payload <- module_ci_lipid_design_payload(layout = "lc_pair")
    seen <- omicDesignCaptureEvents()
    shiny::testServer(mod_lipid_design_builder_server, args = list(data_tbl = function() payload$data_cln,
        config_list = shiny::reactive(payload$config_list), column_mapping = shiny::reactive(payload$column_mapping)), {
        session$flushReact()
        session$setInputs(save_results = 1L)
        expect_null(session$returned())
        session$setInputs(formula_string = "~ 0 + group", new_factor = "WT")
        session$setInputs(add_factor = 1L)
        session$setInputs(new_factor = "KO")
        session$setInputs(add_factor = 2L)
        session$setInputs(selected_runs = c("WT_1", "WT_2"), factor1_select = "WT", factor2_select = "",
            factor3_select = "", replicate_start = 1L)
        expect_match(output$replicate_inputs$html, "Starting replicate number for 2 selected samples:", fixed = TRUE)
        session$setInputs(assign_metadata = 1L)
        session$setInputs(selected_runs = c("KO_1", "KO_2"), factor1_select = "KO")
        session$setInputs(assign_metadata = 2L)
        session$setInputs(contrast_group1 = "KO", contrast_group2 = "WT")
        session$setInputs(add_contrast = 1L)
        expect_identical(output$contrast_factors_info,
            "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group")
        expect_match(output$defined_contrasts_display$html, "KO_vs_WT=groupKO-groupWT", fixed = TRUE)
        session$setInputs(formula_string = "~ group")
        expect_identical(output$contrast_factors_info, "Note: Contrasts will use group names as-is")
        session$setInputs(formula_string = "~ 0 + group")
        session$setInputs(save_results = 2L)
        result <- shiny::isolate(session$returned())
        expect_identical(result$design_matrix$Run, c("KO_1", "KO_2", "WT_1", "WT_2"))
        expect_identical(result$contrasts_tbl$full_format, "KO_vs_WT=groupKO-groupWT")
    })
    expect_identical(seen$events, c("note[warning|5|]: No samples have been assigned to groups.",
        "note[message|5|]: Design saved successfully!"))
})

test_that("W4 the commit's GUI side is today's under the lipid configuration: success and the two refusals", {
    # KIT-103 W7: retargeted from the lipid wrapper (runLipidDesignBuilderObserverShell, deleted) to the
    # shared commit it called, with the same configuration; expectations unchanged.
    config <- omicGuiConfig("lipidomics")
    # Success: the modal, then removed, then the success text; the QC trigger and the tab complete.
    fixture <- list(payload = module_ci_lipid_design_payload(layout = "lc_pair"), paths = module_ci_lipid_design_paths())
    workflow <- module_ci_lipid_design_workflow(fixture$payload)
    seen <- omicDesignCaptureEvents()
    qc <- NULL
    shiny::withReactiveDomain(shiny::MockShinySession$new(), runOmicsDesignCommit(
        results = fixture$payload, workflowData = workflow, experimentPaths = fixture$paths,
        qcTrigger = function(value) qc <<- value, config = config))
    expect_identical(seen$events, c("modal", "remove_modal",
        "note[message|5|]: Design matrix and contrasts saved successfully!"))
    expect_true(qc)
    expect_identical(workflow$tab_status$design_matrix, "complete")
    expect_identical(setdiff(workflow$state_manager$getHistory(), "initial"), "lipid_raw_data_s4")

    # A missing source directory: the text, then the modal removed; nothing written.
    workflow <- module_ci_lipid_design_workflow(fixture$payload)
    seen <- omicDesignCaptureEvents()
    shiny::withReactiveDomain(shiny::MockShinySession$new(), runOmicsDesignCommit(
        results = fixture$payload, workflowData = workflow,
        experimentPaths = list(source_dir = file.path(tempdir(), "kit103-absent")), qcTrigger = function(value) stop("no"),
        config = config))
    expect_identical(seen$events, c("modal", "note[error|15|]: Could not find source directory to save files.",
        "remove_modal"))
    expect_null(workflow$design_matrix)
    expect_identical(workflow$tab_status$design_matrix, "pending")

    # A failed preflight: "Invalid lipidomics design: ...", before the bus.
    drifted <- fixture$payload
    drifted$data_cln$LCMS_Pos$WT_2 <- NULL
    workflow <- module_ci_lipid_design_workflow(fixture$payload)
    paths <- module_ci_lipid_design_paths()
    seen <- omicDesignCaptureEvents()
    shiny::withReactiveDomain(shiny::MockShinySession$new(), runOmicsDesignCommit(
        results = drifted, workflowData = workflow, experimentPaths = paths, qcTrigger = NULL, config = config))
    expect_length(seen$events, 3L)
    expect_identical(seen$events[c(1L, 3L)], c("modal", "remove_modal"))
    expect_match(seen$events[[2L]], "^note\\[error\\|15\\|\\]: Invalid lipidomics design: Assay 'LCMS_Pos' ")
    expect_null(workflow$design_matrix)
    expect_length(list.files(paths$source_dir), 0L)
    expect_identical(workflow$state_manager$getHistory(), "initial")
})

test_that("W4 the shared commit: a failure after the bus write is today's text, after the modal is removed", {
    payload <- module_ci_lipid_design_payload(layout = "lc_pair")
    workflow <- module_ci_lipid_design_workflow(payload)
    paths <- module_ci_lipid_design_paths()
    config <- omicGuiConfig("lipidomics")
    config$design$build_object <- function(results, columnMapping) stop("object refused")
    seen <- omicDesignCaptureEvents()
    qc <- NULL
    shiny::withReactiveDomain(shiny::MockShinySession$new(),
        runOmicsDesignCommit(payload, workflow, paths, function(value) qc <<- value, config))
    expect_identical(seen$events, c("modal", "remove_modal",
        "note[error|15|]: Error saving design matrix results: object refused"))
    expect_null(qc)
    expect_identical(workflow$tab_status$design_matrix, "pending")
    # The bus was written before the failure, as today.
    expect_identical(workflow$design_matrix, payload$design_matrix)
})

test_that("W4 the import plugin: test mode reads the typed path; the confirmation runs the shared import", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    payload <- module_ci_lipid_design_payload(layout = "lc_pair")
    workflow <- module_ci_lipid_design_workflow(payload)
    calls <- list()
    testthat::local_mocked_bindings(runOmicsDesignImport = function(workflowData, experimentPaths, importPath, qcTrigger,
                                                                    config) {
        calls[[length(calls) + 1L]] <<- list(path = importPath, omic = config$shell$omic_type)
        invisible(NULL)
    })
    shiny::testServer(mod_lipid_design_server, args = list(workflow_data = workflow,
        experiment_paths = module_ci_lipid_design_paths()), {
        expect_error(output$import_dir_path, class = "shiny.silent.error")
        session$setInputs(import_dir = "/some/where")
        expect_identical(output$import_dir_path, "/some/where")
        session$setInputs(confirm_import = 1L)
    })
    expect_identical(calls, list(list(path = "/some/where", omic = "lipidomics")))
})

test_that("W4 the shared import (MCI-019.4's path): the bus, the texts in order, and a missing design file", {
    payload <- module_ci_lipid_design_payload(layout = "combined")
    import_root <- module_ci_lipid_design_write_import_pack(tempfile("omic-design-import-"), payload)
    workflow <- module_ci_lipid_design_workflow(payload)
    config <- omicGuiConfig("lipidomics")
    config$design$design_import$qc_baseline <- function(...) invisible(NULL)
    seen <- omicDesignCaptureEvents()
    qc <- NULL
    shiny::withReactiveDomain(omicsOwnerSession(), runOmicsDesignImport(workflow, list(source_dir = import_root),
        import_root, function(value) qc <<- value, config))
    expect_identical(seen$events, c("remove_modal", "note[default|NULL|importing_design]: Importing design files...",
        "remove_note: importing_design",
        "note[message|5|]: Design imported successfully! Loaded 2 assays with 4 samples."))
    expect_true(qc)
    expect_identical(workflow$tab_status$design_matrix, "complete")
    expect_identical(workflow$design_matrix$Run, payload$design_matrix$Run)
    expect_identical(setdiff(workflow$state_manager$getHistory(), "initial"), "lipid_raw_data_s4")

    empty <- withr::local_tempdir()
    workflow <- module_ci_lipid_design_workflow(payload)
    seen <- omicDesignCaptureEvents()
    shiny::withReactiveDomain(omicsOwnerSession(), runOmicsDesignImport(workflow, list(source_dir = empty), empty,
        NULL, config))
    expect_identical(seen$events, c("remove_modal",
        "note[error|10|]: Import failed: 'design_matrix.tab' not found in the selected directory."))
    expect_null(workflow$design_matrix)
})

test_that("W4 the lipid design wrappers keep their signatures", {
    expect_identical(names(formals(mod_lipid_design_ui)), c("id", "profile"))
    expect_identical(names(formals(mod_lipid_design_server)),
        c("id", "workflow_data", "experiment_paths", "volumes", "qc_trigger", "profile"))
    expect_identical(names(formals(mod_lipid_design_builder_ui)), "id")
    expect_identical(names(formals(mod_lipid_design_builder_server)), c("id", "data_tbl", "config_list",
        "column_mapping", "existing_design_matrix", "existing_contrasts", "state_owner", "profile"))
    expect_identical(names(formals(validateLipidDesignAlignment)), c("designMatrix", "assayList", "columnMapping",
        "sampleIdCol", "groupCol", "requireGroups", "requireReplicates"))
    expect_identical(names(formals(validateLipidDesignContrasts)), c("designMatrix", "contrastsTbl", "formulaString",
        "allowEmpty"))
    expect_identical(names(formals(validateLipidDesignDaPreflight)), c("designMatrix", "assayList", "contrastsTbl",
        "formulaString", "columnMapping", "sampleIdCol", "groupCol", "requireContrasts"))
})

# ------------------------------------------------------------------------------------------------
# KIT-103 W5: metab on the shared design tab and builder.
# ------------------------------------------------------------------------------------------------

test_that("W5 the metab design tab's UI is today's bytes, through its wrappers and the shared modules", {
    expect_identical(omicDesignUiDigest(mod_metab_design_ui("design")), omicDesignFrozen$metab_frame)
    expect_identical(omicDesignUiDigest(mod_metab_design_builder_ui("builder")), omicDesignFrozen$metab_builder)
    config <- omicGuiConfig("metabolomics")
    expect_identical(omicDesignUiDigest(mod_omic_design_ui("design", config)), omicDesignFrozen$metab_frame)
    expect_identical(omicDesignUiDigest(mod_omic_design_builder_ui("builder", config)), omicDesignFrozen$metab_builder)
    ns <- shiny::NS("proxy1")
    expect_identical(omicDesignUiDigest(buildOmicsDesignImportModal(ns, config, testMode = TRUE)),
        omicDesignFrozen$metab_modal_test)
    expect_identical(omicDesignUiDigest(buildOmicsDesignImportModal(ns, config, testMode = FALSE)),
        omicDesignFrozen$metab_modal_app)
})

test_that("W5 the metab design tab's server keeps today's outputs, modal and the two status outputs' error", {
    for (mode in c("true", "false")) {
        withr::local_envvar(MULTISCHOLAR_TEST_MODE = mode)
        testthat::local_mocked_bindings(shinyDirChoose = function(...) invisible(NULL), .package = "shinyFiles")
        seen <- new.env(parent = emptyenv())
        testthat::local_mocked_bindings(showModal = function(ui, ...) seen$modal <- ui, .package = "shiny")
        payload <- module_ci_metab_design_payload(layout = "lc_pair")
        workflow <- module_ci_metab_design_workflow(payload)
        shiny::testServer(mod_metab_design_server, args = list(workflow_data = workflow,
            experiment_paths = module_ci_metab_design_paths()), {
            session$flushReact()
            outputs <- sort(sub("^proxy[0-9]+-", "", names(session$.__enclos_env__$private$outs)))
            expect_identical(outputs, omicDesignFrozen$outputs)
            for (name in c("data_available", "design_matrix_exists")) {
                error <- tryCatch(output[[name]], error = identity)
                expect_s3_class(error, "error")
                expect_identical(conditionMessage(error), "unused arguments (self, name)")
            }
            expect_identical(output$assays_preview, "Included assays: LCMS_Pos, LCMS_Neg")
            expect_identical(output$`builder-removed_samples_display`, "No samples have been removed.")
            expect_identical(output$`builder-tech_rep_summary`, "No technical replicates assigned yet.")
            session$setInputs(show_import_modal = 1L)
        })
        expect_identical(omicDesignUiDigest(seen$modal),
            if (identical(mode, "true")) omicDesignFrozen$metab_modal_test else omicDesignFrozen$metab_modal_app)
    }
})

test_that("W5 the metab builder's texts, refusal, run order and save are today's", {
    payload <- module_ci_metab_design_payload(layout = "lc_pair")
    test_names <- c("design", "levels", "runs", "group", "clicks", "tech", "set", "remove", "contrast",
        "technical", "out", "workflow", "fixture", "i", "result", "samples", "assays", "mapping", "seen")
    seen <- omicDesignCaptureEvents()
    shiny::testServer(mod_metab_design_builder_server, args = list(data_tbl = function() payload$data_cln,
        config_list = shiny::reactive(payload$config_list), column_mapping = shiny::reactive(payload$column_mapping)), {
        module_names <- c(ls(session$env, all.names = TRUE), ls(parent.env(session$env), all.names = TRUE))
        expect_s3_class(design_matrix, "reactiveVal")
        expect_identical(intersect(module_names, test_names), character())
        session$flushReact()
        expect_identical(shiny::isolate(design_matrix()$Run), c("KO_1", "KO_2", "WT_1", "WT_2"))
        session$setInputs(save_results = 1L)
        expect_null(session$returned())
        session$setInputs(formula_string = "~ 0 + group", new_factor = "WT")
        session$setInputs(add_factor = 1L)
        session$setInputs(new_factor = "KO")
        session$setInputs(add_factor = 2L)
        session$setInputs(selected_runs = c("WT_1", "WT_2"), factor1_select = "WT", factor2_select = "",
            factor3_select = "", replicate_start = 1L)
        expect_match(output$replicate_inputs$html, "Starting replicate number for 2 selected samples:", fixed = TRUE)
        expect_no_match(output$replicate_inputs$html, "data-testid", fixed = TRUE)
        session$setInputs(assign_metadata = 1L)
        session$setInputs(selected_runs = c("KO_1", "KO_2"), factor1_select = "KO")
        session$setInputs(assign_metadata = 2L)
        session$setInputs(contrast_group1 = "KO", contrast_group2 = "WT")
        session$setInputs(add_contrast = 1L)
        expect_identical(output$contrast_factors_info,
            "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group")
        expect_match(output$defined_contrasts_display$html, "KO_vs_WT=groupKO-groupWT", fixed = TRUE)
        session$setInputs(formula_string = "~ group")
        expect_identical(output$contrast_factors_info, "Note: Contrasts will use group names as-is")
        session$setInputs(formula_string = "~ 0 + group")
        session$setInputs(save_results = 2L)
        result <- shiny::isolate(session$returned())
        expect_identical(result$design_matrix$Run, c("KO_1", "KO_2", "WT_1", "WT_2"))
        expect_identical(result$contrasts_tbl$full_format, "KO_vs_WT=groupKO-groupWT")
        expect_identical(names(result$data_cln$LCMS_Pos),
            c("Feature.Name", "Annotation", "KO_1", "KO_2", "WT_1", "WT_2"))
        expect_false(any(c("group_levels", "excluded_runs") %in% names(result)))
    })
    expect_identical(seen$events, c("note[warning|5|]: No samples have been assigned to groups.",
        "note[message|5|]: Design saved successfully!"))
})

test_that("W5 the metab import checks its files before closing the modal (P-11), as today", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    empty <- withr::local_tempdir()
    design_only <- withr::local_tempdir()
    writeLines("Run\tgroup", file.path(design_only, "design_matrix.tab"))
    payload <- module_ci_metab_design_payload(layout = "lc_pair")
    cases <- list(
        list(dir = empty, text = "Import failed: 'design_matrix.tab' not found in the selected directory."),
        list(dir = design_only, text = paste("Import failed: 'assay_manifest.txt' not found.",
            "Cannot determine which assay files to load.")))
    for (case in cases) {
        workflow <- module_ci_metab_design_workflow(payload)
        seen <- omicDesignCaptureEvents()
        shiny::testServer(mod_metab_design_server, args = list(workflow_data = workflow,
            experiment_paths = module_ci_metab_design_paths()), {
            session$setInputs(import_dir = case$dir)
            session$setInputs(confirm_import = 1L)
        })
        expect_identical(seen$events, paste0("note[error|10|]: ", case$text))
        expect_null(workflow$design_matrix)
    }
})

test_that("W5 the shared import under the metab configuration: the bus, the state and the texts in order", {
    payload <- module_ci_metab_design_payload(layout = "combined")
    import_root <- module_ci_metab_design_write_import_pack(tempfile("omic-design-metab-import-"), payload)
    workflow <- module_ci_metab_design_workflow(payload)
    config <- omicGuiConfig("metabolomics")
    config$design$design_import$qc_baseline <- function(...) invisible(NULL)
    seen <- omicDesignCaptureEvents()
    qc <- NULL
    shiny::withReactiveDomain(omicsOwnerSession(), runOmicsDesignImport(workflow, list(source_dir = import_root),
        import_root, function(value) qc <<- value, config))
    expect_identical(seen$events, c("remove_modal", "note[default|NULL|importing_design]: Importing design files...",
        "remove_note: importing_design",
        "note[message|5|]: Design imported successfully! Loaded 2 assays with 4 samples."))
    expect_true(qc)
    expect_identical(workflow$tab_status$design_matrix, "complete")
    expect_identical(workflow$design_matrix$Run, payload$design_matrix$Run)
    expect_identical(setdiff(workflow$state_manager$getHistory(), "initial"), "metab_raw_data_s4")
})

# KIT-101 F-01 (test-metab-04-design-module-characterization.R:447, :466 at de2d891c), retargeted by
# KIT-103 W5 from metab's builder-results observer to the shared commit with the metab
# configuration; the expectations are unchanged.
runMetabDesignBuilderResultsDefaultFlow <- function(results, sourceDir) {
  ui_calls <- character()
  notifications <- character()
  testthat::local_mocked_bindings(
    showModal = function(...) ui_calls <<- c(ui_calls, "showModal"),
    removeModal = function(...) ui_calls <<- c(ui_calls, "removeModal"),
    showNotification = function(ui, ...) {
      notifications <<- c(notifications, ui)
      ui_calls <<- c(ui_calls, "showNotification")
    },
    .package = "shiny"
  )
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$column_mapping <- list(sample_columns = c("Sample1", "Sample2"))
  outcome <- tryCatch(
    {
      runOmicsDesignCommit(
        results = results,
        workflowData = workflow_data,
        experimentPaths = list(source_dir = sourceDir),
        qcTrigger = NULL,
        config = omicGuiConfig("metabolomics")
      )
      "returned"
    },
    error = function(e) paste("observer error:", conditionMessage(e))
  )
  list(
    outcome = outcome,
    ui_calls = ui_calls,
    notifications = notifications,
    workflow_data = workflow_data
  )
}

makeMetabDesignBuilderResults <- function(runs) {
  list(
    design_matrix = data.frame(
      Run = runs,
      group = c("Control", "Treatment"),
      replicates = c(1L, 1L),
      stringsAsFactors = FALSE
    ),
    data_cln = list(LCMS_Pos = data.frame(
      feature = c("M1", "M2"),
      Sample1 = c(10, 11),
      Sample2 = c(12, 13),
      stringsAsFactors = FALSE
    )),
    contrasts_tbl = NULL,
    config_list = list(normalization = "none")
  )
}

test_that("KIT-101 F-01: an invalid metab design notifies, closes the modal and leaves the bus unwritten", {
  source_dir <- tempfile("metab-design-f01-")
  dir.create(source_dir)
  on.exit(unlink(source_dir, recursive = TRUE), add = TRUE)

  run <- runMetabDesignBuilderResultsDefaultFlow(
    makeMetabDesignBuilderResults(c("Sample1", "Sample9")),
    source_dir
  )

  expect_identical(run$outcome, "returned")
  expect_identical(run$ui_calls, c("showModal", "showNotification", "removeModal"))
  expect_length(run$notifications, 1)
  expect_match(run$notifications, "^Invalid metabolomics design: ")
  expect_null(run$workflow_data$design_matrix)
  expect_null(run$workflow_data$config_list)
  expect_length(list.files(source_dir), 0)
})

test_that("KIT-101 F-01: a missing metab source directory is reported before the bus is written", {
  run <- runMetabDesignBuilderResultsDefaultFlow(
    makeMetabDesignBuilderResults(c("Sample1", "Sample2")),
    file.path(tempdir(), "metab-design-f01-missing-dir")
  )

  expect_identical(run$outcome, "returned")
  expect_identical(run$ui_calls, c("showModal", "showNotification", "removeModal"))
  expect_identical(run$notifications, "Could not find source directory to save files.")
  expect_null(run$workflow_data$design_matrix)
  expect_null(run$workflow_data$config_list)
})

test_that("W5 the metab design wrappers keep their signatures", {
    expect_identical(names(formals(mod_metab_design_ui)), "id")
    expect_identical(names(formals(mod_metab_design_server)),
        c("id", "workflow_data", "experiment_paths", "volumes", "qc_trigger"))
    expect_identical(names(formals(mod_metab_design_builder_ui)), "id")
    expect_identical(names(formals(mod_metab_design_builder_server)), c("id", "data_tbl", "config_list",
        "column_mapping", "existing_design_matrix", "existing_contrasts", "state_owner"))
    expect_identical(names(formals(registerMetabDesignPreviewOutputs)),
        c("output", "workflowData", "renderDt", "renderText", "req"))
    expect_identical(names(formals(validateMetabDesignAlignment)), c("designMatrix", "assayList", "columnMapping",
        "sampleIdCol", "groupCol", "requireGroups", "requireReplicates"))
    expect_identical(names(formals(validateMetabDesignContrasts)), c("designMatrix", "contrastsTbl", "formulaString",
        "allowEmpty"))
    expect_identical(names(formals(validateMetabDesignDaPreflight)), c("designMatrix", "assayList", "contrastsTbl",
        "formulaString", "columnMapping", "sampleIdCol", "groupCol", "requireContrasts"))
})

# KIT-105 W4: the masking check extended to the shared normalisation frame (R/mod_omic_norm*.R). testServer evaluates a
# test's code over the module's frame, so a frame local named after a variable the norm tests read inside testServer
# (test-omic-norm-module.R's omicNormDrive) would be read in its place.
test_that("KIT-105 W4 the normalisation frame keeps its returned reactive values and no name that masks a test's own variables", {
    skip_if_not(workflowOwnerCapable())
    test_names <- c("spec", "root", "source_dir", "paths", "workflow", "log", "acc", "ruv_results", "stubs", "plot_pca",
        "qc_kernel", "annotate", "tab", "steps", "step", "set", "read", "capture", "ctx", "norm_data", "defaults",
        "apply_updates", "field", "i", "pending", "values", "manager", "nd", "log_lines", "new_log", "notes", "qc_files",
        "pngs", "encoded", "out_names", "others", "files", "record", "slots", "value", "fixture", "managed", "result")
    object <- module_ci_lipid_norm_object(positive_only = TRUE)
    bus <- module_ci_lipid_norm_workflow(object)
    shiny::testServer(mod_lipid_norm_server, args = list(id = "norm", workflow_data = bus, experiment_paths = list(),
        omic_type = "lipidomics", experiment_label = "kit105"), {
        module_names <- c(ls(session$env, all.names = TRUE), ls(parent.env(session$env), all.names = TRUE))
        expect_identical(intersect(module_names, test_names), character())
        session$flushReact()
        expect_true(shiny::is.reactivevalues(session$returned))
        expect_identical(shiny::isolate(session$returned$assay_names), c("LCMS_Pos", "GCMS"))
    })
})
