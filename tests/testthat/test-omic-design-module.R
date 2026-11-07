# KIT-103 W4: the shared design tab (R/mod_omic_design.R) and builder (R/mod_omic_design_builder.R),
# with lipid and trx on them through the lipid wrappers. Each block pins today's behaviour: the
# frozen values (UI digests, output names, notification texts and their order) were taken from the
# lipid design files at 616e1147, and the wrapper blocks pass there unchanged
# (docs/qa/workflow-execution/kit/kit-103/w4/test-first.txt).

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
        "contrasts_preview", "data_available", "design_matrix_exists", "design_matrix_preview", "import_dir_path"))

test_that("W4 the design tab's UI through the lipid wrappers is today's bytes (lipid, trx; frame and builder)", {
    expect_identical(omicDesignUiDigest(mod_lipid_design_ui("design")), omicDesignFrozen$lipid_frame)
    expect_identical(omicDesignUiDigest(mod_lipid_design_ui("design", profile = trxOmicProfile())),
        omicDesignFrozen$trx_frame)
    expect_identical(omicDesignUiDigest(mod_lipid_design_builder_ui("builder")), omicDesignFrozen$builder)
})

test_that("W4 the shared frame, builder, import modal and reset modal render today's bytes from the configuration", {
    expect_identical(omicDesignUiDigest(mod_omic_design_ui("design", lipidDesignGuiConfig())),
        omicDesignFrozen$lipid_frame)
    expect_identical(omicDesignUiDigest(mod_omic_design_ui("design", lipidDesignGuiConfig(trxOmicProfile()))),
        omicDesignFrozen$trx_frame)
    expect_identical(omicDesignUiDigest(mod_omic_design_builder_ui("builder", omicGuiConfig("lipidomics"))),
        omicDesignFrozen$builder)
    # Trx renders the lipid builder (its builder texts are lipid's).
    expect_identical(omicDesignUiDigest(mod_omic_design_builder_ui("builder", omicGuiConfig("transcriptomics"))),
        omicDesignFrozen$builder)
    ns <- shiny::NS("design")
    config <- lipidDesignGuiConfig()
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

test_that("W4 the commit's GUI side is today's through the lipid wrapper: success and the two refusals", {
    # Success: the modal, then removed, then the success text; the QC trigger and the tab complete.
    fixture <- list(payload = module_ci_lipid_design_payload(layout = "lc_pair"), paths = module_ci_lipid_design_paths())
    workflow <- module_ci_lipid_design_workflow(fixture$payload)
    seen <- omicDesignCaptureEvents()
    qc <- NULL
    shiny::withReactiveDomain(shiny::MockShinySession$new(), runLipidDesignBuilderObserverShell(
        results = fixture$payload, workflowData = workflow, experimentPaths = fixture$paths,
        qcTrigger = function(value) qc <<- value))
    expect_identical(seen$events, c("modal", "remove_modal",
        "note[message|5|]: Design matrix and contrasts saved successfully!"))
    expect_true(qc)
    expect_identical(workflow$tab_status$design_matrix, "complete")
    expect_identical(setdiff(workflow$state_manager$getHistory(), "initial"), "lipid_raw_data_s4")

    # A missing source directory: the text, then the modal removed; nothing written.
    workflow <- module_ci_lipid_design_workflow(fixture$payload)
    seen <- omicDesignCaptureEvents()
    shiny::withReactiveDomain(shiny::MockShinySession$new(), runLipidDesignBuilderObserverShell(
        results = fixture$payload, workflowData = workflow,
        experimentPaths = list(source_dir = file.path(tempdir(), "kit103-absent")), qcTrigger = function(value) stop("no")))
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
    shiny::withReactiveDomain(shiny::MockShinySession$new(), runLipidDesignBuilderObserverShell(
        results = drifted, workflowData = workflow, experimentPaths = paths))
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
    config <- lipidDesignGuiConfig()
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
    config <- lipidDesignGuiConfig()
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
