# KIT-138 (the integration shell and its Import tab) and KIT-139 (its Sample Matching tab). The scanned project is the
# KIT-137 tier-1 corpus copied into `<project>/integration/<omic>/` (four omics, six data sets, twenty biological
# samples, KO_7 absent from metabolomics; the lipidomics, metabolomics and transcriptomics runs two per sample), and,
# for the internal standards, the few-samples note and the refusals, two small inputs written by the contract's writer
# in the form of test-integration-model.R. Expected counts are read from the corpus files themselves, not from the
# code under test. The whole surface is behind `MultiScholaR.integration_gui` (off by default).

intGuiCorpus <- function() normalizePath(testthat::test_path("..", "testdata", "integration", "corpus"), mustWork = FALSE)
intGuiOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intGuiPaths <- function(root) {
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = file.path(root, "results", "integration"),
        multiomic_inputs_dir = file.path(root, "results", "integration", "multiomic_inputs"))
    dir.create(paths$integration_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

# A project whose integration folder holds the corpus's inputs for `omics`.
intGuiProject <- function(omics = intGuiOmics, .local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(intGuiCorpus()), "the tier-1 corpus is not present")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    paths <- intGuiPaths(root)
    for (omic in omics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir)
        stopifnot(all(file.copy(list.files(file.path(intGuiCorpus(), omic), full.names = TRUE), dir)))
    }
    list(root = root, paths = paths)
}

# The corpus as read without the code under test: per data set its features and runs, per omic its biological samples.
intGuiCorpusCounts <- function() {
    views <- list()
    samples <- list()
    for (omic in intGuiOmics) {
        manifest <- jsonlite::read_json(file.path(intGuiCorpus(), omic, "integration_input.json"))
        design <- utils::read.delim(file.path(intGuiCorpus(), omic, "samples.tsv"), stringsAsFactors = FALSE)
        samples[[omic]] <- design[[manifest$biological_sample_column]]
        for (view in manifest$views) {
            matrix <- utils::read.delim(file.path(intGuiCorpus(), omic, view$matrix_file), check.names = FALSE)
            views[[view$name]] <- data.frame(view = view$name, omic = omic, features = nrow(matrix),
                runs = ncol(matrix) - 1L, stringsAsFactors = FALSE)
        }
    }
    list(views = do.call(rbind, views), samples = samples)
}

# Two small inputs (test-integration-model.R's form): proteomics one run per biological sample (A_1..A_3, B_1, B_2; P4
# observed once), lipidomics two runs per biological sample (A_3 absent; ITSD1 an internal standard).
intGuiInput <- function(root, omic, views, samples, standards = list()) {
    matrices <- lapply(views, function(m) data.frame(feature_id = rownames(m), m, check.names = FALSE,
        stringsAsFactors = FALSE))
    features <- lapply(stats::setNames(names(views), names(views)), function(view) {
        table <- data.frame(feature_id = rownames(views[[view]]), display_name = tolower(rownames(views[[view]])),
            stringsAsFactors = FALSE)
        if (!is.null(standards[[view]])) table$is_internal_standard <- standards[[view]]
        table
    })
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
        source_state = list(name = "kit138-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, matrices, features)
    dir
}

intGuiSamples <- function(runs, biological) {
    data.frame(Run = runs, group = sub("_.*$", "", biological), replicates = as.integer(sub("^.*_", "", biological)),
        tech_rep_group = biological, stringsAsFactors = FALSE)
}

intGuiSmallProject <- function(lipid_biological = c("A_1", "A_1", "A_2", "A_2", "B_1", "B_1", "B_2", "B_2"),
    shift = 0, .local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    paths <- intGuiPaths(root)
    biological <- c("A_1", "A_2", "A_3", "B_1", "B_2")
    runs <- paste0("p", seq_along(biological))
    values <- rbind(P1 = c(10, 11, 12, 13, 14), P2 = c(5, 7, 6, 9, 8), P3 = c(1, 1.5, 2, 2.5, 3.5),
        P4 = c(NA, NA, 4, NA, NA)) + shift
    colnames(values) <- runs
    intGuiInput(root, "proteomics", list(proteomics = values), intGuiSamples(runs, biological))
    runs <- paste0("l", seq_along(lipid_biological))
    values <- rbind(L1 = c(20, 22, 21, 23, 25, 27, 26, 28), L2 = c(3, NA, 4, 5, 6, 7, NA, NA),
        ITSD1 = c(15, 15, 15, 15, 15, 15, 15, 15))
    colnames(values) <- runs
    intGuiInput(root, "lipidomics", list(lipidomics_LCMS_Pos = values), intGuiSamples(runs, lipid_biological),
        standards = list(lipidomics_LCMS_Pos = c(FALSE, FALSE, TRUE)))
    list(root = root, paths = paths)
}

intGuiBus <- function() {
    do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
}

# An import committed on a store outside the app (the Sample Matching tab's starting point).
intGuiImported <- function(paths, views = NULL, .local_envir = parent.frame()) {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, "store"))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    request <- integrationImportRequest(scan, views %||% scan$views$view)
    bus <- intGuiBus()
    shiny::isolate({
        bus$state_manager <- manager
        bus$integration_import <- integrationCommitImport(manager, request)
        status <- bus$tab_status
        status$integration_import <- "complete"
        status$integration_alignment <- "pending"
        bus$tab_status <- status
    })
    bus
}

intGuiHtml <- function(tag) as.character(htmltools::renderTags(tag)$html)

# ------------------------------------------------------------------------------------------------
# KIT-138: the gate, the shell, the home page.
# ------------------------------------------------------------------------------------------------

test_that("KIT-138: with the gate off the home page and the integration module are as before", {
    withr::local_options(MultiScholaR.integration_gui = NULL)
    expect_false(integrationGuiEnabled())
    page <- intGuiHtml(app_ui(NULL))
    expect_match(page, '<div id="integration_box" class="omic-selection-box disabled" data-testid="tile-integration">',
        fixed = TRUE)
    expect_match(page, "Multi-omics integration (available when 2+ omics selected)", fixed = TRUE)
    expect_match(page, "if (omic === 'integration') return; // Integration is auto-enabled", fixed = TRUE)
    expect_match(page, "if (selectedOmics.length >= 2) {", fixed = TRUE)
    expect_false(grepl("KIT-138", page, fixed = TRUE))
    script <- "if (omic === 'integration') return; // Integration is auto-enabled"
    expect_identical(integrationHomeScript(script), script)
    # The panel app_server shows for an omic without a module (R/app_server.R, "Module not yet implemented").
    expect_identical(mod_integration_ui("integration_workflow"), shiny::div(class = "alert alert-warning",
        shiny::h4("Module not yet implemented"), shiny::p("The integration workflow module is under development.")))
    expect_null(mod_integration_server("integration_workflow", list(integration = list())))
})

test_that("KIT-138: with the gate on the Integration tile is live and the shell has the Import and Sample Matching tabs", {
    withr::local_options(MultiScholaR.integration_gui = TRUE)
    texts <- integrationGuiTexts()
    page <- intGuiHtml(app_ui(NULL))
    expect_match(page, paste0('<div id="integration_box" class="omic-selection-box" ',
        'onclick="toggleOmicSelection(&#39;integration&#39;)" data-testid="tile-integration">'), fixed = TRUE)
    expect_match(page, texts[["home_tile_text"]], fixed = TRUE)
    expect_false(grepl("if (omic === 'integration') return;", page, fixed = TRUE))
    expect_match(page, "if (true) { // KIT-138: the Integration tile stays enabled", fixed = TRUE)
    expect_error(integrationHomeScript("a script without the lines"), class = "integration_gui_invalid")
    config <- integrationGuiConfig()
    expect_identical(vapply(config$shell$tabs, `[[`, "", "status_key"), c("integration_import", "integration_alignment"))
    ui <- intGuiHtml(mod_integration_ui("integration_workflow"))
    for (id in c("integration-workflow-tabs", "integration-tab-import", "integration-tab-alignment",
        "integration-import-commit", "integration-alignment-confirm")) {
        expect_match(ui, sprintf('data-testid="%s"', id), fixed = TRUE)
    }
    expect_match(ui, texts[["tab_alignment"]], fixed = TRUE)
    project <- intGuiProject()
    shiny::testServer(mod_integration_server, args = list(project_dirs = list(integration = project$paths),
        omic_type = "integration", experiment_label = "kit138"), {
        session$flushReact()
        bus <- session$returned
        expect_identical(bus$tab_status, list(integration_import = "pending", integration_alignment = "disabled"))
        progress <- output$workflow_progress$html
        expect_match(progress, texts[["step_import"]], fixed = TRUE)
        expect_match(progress, texts[["step_alignment"]], fixed = TRUE)
        # The Import tab scanned the project when the module started.
        expect_identical(bus$integration_scan$omics$omic, intGuiOmics)
    })
})

# ------------------------------------------------------------------------------------------------
# KIT-138: the Import tab.
# ------------------------------------------------------------------------------------------------

test_that("KIT-138: the scan lists every omic's data sets with their sizes and correction, as saved", {
    project <- intGuiProject()
    counts <- intGuiCorpusCounts()
    scan <- integrationScanProject(project$paths$integration_dir)
    expect_identical(scan$omics$omic, intGuiOmics)
    expect_identical(scan$omics$status, rep("ready", 4L))
    expect_identical(scan$omics$data_type, vapply(intGuiOmics, function(omic) {
        jsonlite::read_json(file.path(intGuiCorpus(), omic, "integration_input.json"))$workflow
    }, character(1), USE.NAMES = FALSE))
    texts <- integrationGuiTexts()
    expect_identical(scan$omics$values, rep(texts[["values_normalised"]], 4L))
    expect_identical(scan$omics$correction, rep(texts[["correction_none"]], 4L))
    expected <- counts$views[match(scan$views$view, counts$views$view), ]
    expect_setequal(scan$views$view, counts$views$view)
    expect_identical(scan$views$omic, expected$omic)
    expect_identical(as.integer(scan$views$features), as.integer(expected$features))
    expect_identical(as.integer(scan$views$runs), as.integer(expected$runs))
    expect_identical(as.integer(scan$views$biological_samples),
        unname(vapply(scan$views$omic, function(omic) length(unique(counts$samples[[omic]])), integer(1))))
    display <- integrationScanDisplay(scan)
    expect_identical(names(display), unname(texts[c("col_omic", "col_data_type", "col_values", "col_correction",
        "col_saved", "col_status")]))
    expect_identical(display[[1L]], c("Lipidomics", "Metabolomics", "Proteomics", "Transcriptomics"))
    expect_identical(display[[5L]], rep("2026-11-17 00:00 UTC", 4L))
})

test_that("KIT-138: an omic not saved, an input that fails the contract and one changed since the import are listed in words", {
    project <- intGuiProject(c("lipidomics", "metabolomics", "transcriptomics"))
    texts <- integrationGuiTexts()
    dir.create(file.path(project$paths$integration_dir, "proteomics"))
    matrix <- file.path(project$paths$integration_dir, "lipidomics", "matrix_lipidomics_LCMS_Pos.tsv")
    cat("\n", file = matrix, append = TRUE)
    metabolomics <- file.path(project$paths$integration_dir, "metabolomics", "integration_input.json")
    scan <- integrationScanProject(project$paths$integration_dir,
        imported = c(metabolomics = strrep("0", 64L), transcriptomics = digest::digest(file = file.path(
            project$paths$integration_dir, "transcriptomics", "integration_input.json"), algo = "sha256")))
    status <- stats::setNames(scan$omics$status, scan$omics$omic)
    expect_identical(status[["proteomics"]], "not_saved")
    expect_identical(scan$omics$reason[scan$omics$omic == "proteomics"], texts[["status_not_saved"]])
    # The reader's own message.
    expect_identical(status[["lipidomics"]], "refused")
    expect_identical(scan$omics$reason[scan$omics$omic == "lipidomics"],
        tryCatch(readIntegrationInput(dirname(matrix)), integration_input_invalid = conditionMessage))
    expect_match(scan$omics$reason[scan$omics$omic == "lipidomics"], "does not match its sha256", fixed = TRUE)
    expect_identical(status[["metabolomics"]], "changed")
    expect_identical(status[["transcriptomics"]], "ready")
    expect_identical(scan$omics$manifest_sha256[scan$omics$omic == "metabolomics"],
        digest::digest(file = metabolomics, algo = "sha256"))
    # A refused omic's data sets are not offered; a changed one's are (importing again is the remedy).
    expect_setequal(integrationScanReadyViews(scan)$view, c("metabolomics_LCMS_Pos", "metabolomics_GCMS", "transcriptomics"))
    expect_error(integrationImportRequest(scan, c("lipidomics_LCMS_Pos", "transcriptomics")),
        class = "integration_input_not_ready")
    expect_error(integrationImportRequest(scan, c("metabolomics_GCMS", "metabolomics_LCMS_Pos")),
        regexp = sprintf(texts[["import_need_two_omics"]], 1L), fixed = TRUE, class = "integration_too_few_omics")
    request <- integrationImportRequest(scan, c("transcriptomics", "metabolomics_GCMS"))
    expect_identical(request$omics, c("metabolomics", "transcriptomics"))
    expect_identical(request$views, c("metabolomics_GCMS", "transcriptomics"))
})

test_that("KIT-138: the correction for unwanted variation each input names is stated in words", {
    texts <- integrationGuiTexts()
    expect_identical(integrationCorrectionText(list(matrix_role = "normalised", correction = NULL)), texts[["correction_none"]])
    declared <- list(matrix_role = "corrected", correction = list(id = "C2", method = "RUVIII", k = 2L,
        ruv_design = "groups_all_features", ruv_k = 2L))
    expect_identical(integrationCorrectionText(declared),
        paste0(sprintf(texts[["correction_design"]], omicGuiRuvRecordTexts()[["design_groups_all_features"]]), "; ",
            sprintf(texts[["correction_k"]], "2")))
    legacy <- list(matrix_role = "corrected", correction = list(id = "C1", method = "RUVIII_C_Varying", k = 3L,
        replicate_structure = NULL, ruv_design = NULL))
    expect_identical(integrationCorrectionText(legacy),
        "applied: the approach of earlier versions (RUV-III-C; replicate sets: the experimental groups); 3 factor(s)")
    per_view <- list(matrix_role = "corrected", correction = list(id = "C4", method = "RUVIII", ruv_design = "technical_replicates",
        views = list(list(assay = "LCMS_Pos", applied = TRUE, k = 2L), list(assay = "LCMS_Neg", applied = FALSE, k = 0L))))
    expect_identical(integrationCorrectionText(per_view),
        paste0(sprintf(texts[["correction_design"]], omicGuiRuvRecordTexts()[["design_technical_replicates"]]),
            "; LCMS_Pos: 2 factor(s); LCMS_Neg: not applied"))
    # Read back from a saved input: the scan shows what the file says.
    project <- intGuiProject(c("lipidomics"))
    input <- readIntegrationInput(file.path(intGuiCorpus(), "proteomics"))
    manifest <- input$manifest[c("omic_type", "workflow", "scale", "source_state", "sample_id_column", "group_column",
        "technical_replicate_column", "biological_sample_column")]
    manifest$matrix_role <- "corrected"
    manifest$correction <- list(state = "ruv_covariates", id = "C2", method = "RUVIII", k = 2L, ruv_design = "groups_all_features",
        ruv_k = 2L, replicate_structure = "group", grouping_variable = "group", ruv_k_estimate = NA)
    writeIntegrationInput(file.path(project$paths$integration_dir, "proteomics"), manifest, input$samples,
        lapply(input$views, `[[`, "matrix"), lapply(input$views, `[[`, "features"))
    scan <- integrationScanProject(project$paths$integration_dir)
    row <- scan$omics[scan$omics$omic == "proteomics", ]
    expect_identical(row$values, texts[["values_corrected"]])
    expect_identical(row$correction, integrationCorrectionText(declared))
})

test_that("KIT-138: the Import tab refuses fewer than two omics and commits the selected omics through the import", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    project <- intGuiProject()
    counts <- intGuiCorpusCounts()
    texts <- integrationGuiTexts()
    bus <- intGuiBus()
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)))
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = project$paths), {
        session$flushReact()
        expect_identical(bus$integration_scan$omics$omic, intGuiOmics)
        expect_match(output$scan_summary$html, sprintf(texts[["import_found"]], 4L, 0L), fixed = TRUE)
        expect_match(output$view_choice$html, "lipidomics_LCMS_Neg (Lipidomics; 150 features, 20 biological samples)",
            fixed = TRUE)
        session$setInputs(views = c("lipidomics_LCMS_Pos", "lipidomics_LCMS_Neg"), import = 1)
        expect_match(output$import_status$html, sprintf(texts[["import_need_two_omics"]], 1L), fixed = TRUE)
        expect_null(bus$integration_import)
        initial <- omicWorkflowInitialBus(integrationGuiConfig()$shell)$tab_status
        expect_identical(bus$tab_status, initial)
        session$setInputs(views = setdiff(counts$views$view, "metabolomics_GCMS"), import = 2)
        manager <- bus$state_manager
        expect_identical(manager$getHistory(), c("initial", "integration_inputs"))
        roles <- manager$.context$snapshot$file_dependencies
        expect_setequal(names(roles), paste0("integration_input_", intGuiOmics))
        summary <- bus$integration_import
        expect_identical(summary$omics, intGuiOmics)
        # The selected data sets in the state's order; the unselected one is imported with its omic but not selected.
        expect_identical(summary$views, c("lipidomics_LCMS_Pos", "lipidomics_LCMS_Neg", "metabolomics_LCMS_Pos",
            "proteomics", "transcriptomics"))
        expect_identical(summary$manifests, vapply(stats::setNames(intGuiOmics, intGuiOmics), function(omic) {
            digest::digest(file = file.path(project$paths$integration_dir, omic, "integration_input.json"), algo = "sha256")
        }, character(1)))
        runs <- sum(vapply(counts$samples, length, integer(1)))
        expect_identical(summary$runs, runs)
        expect_identical(summary$biological_samples, length(unique(unlist(counts$samples))))
        expect_match(output$import_status$html, sprintf(texts[["import_done"]], 4L, 5L, runs,
            length(unique(unlist(counts$samples)))), fixed = TRUE)
        expect_identical(bus$tab_status$integration_import, "complete")
        expect_identical(bus$tab_status, integrationImportStatus(initial))
        # The project searched again after the import: nothing changed.
        session$setInputs(rescan = 1)
        expect_identical(bus$integration_scan$omics$status, rep("ready", 4L))
        # Importing again commits on a new store; the earlier one is closed.
        first <- bus$state_manager
        session$setInputs(import = 3)
        expect_false(identical(bus$state_manager, first))
        expect_identical(bus$state_manager$getHistory(), c("initial", "integration_inputs"))
        expect_identical(bus$integration_import$views, summary$views)
    })
    # After an import the tab after it opens and every later one closes, whatever the tabs (KIT-140 adds one).
    expect_identical(integrationImportStatus(list(integration_import = "complete", second = "complete", third = "pending")),
        list(integration_import = "complete", second = "pending", third = "disabled"))
    expect_identical(integrationImportStatus(list(integration_import = "pending")), list(integration_import = "complete"))
})

# ------------------------------------------------------------------------------------------------
# KIT-139: the Sample Matching tab.
# ------------------------------------------------------------------------------------------------

test_that("KIT-139: the corpus's samples are matched by the biological sample, with missing samples and averaged runs stated", {
    project <- intGuiProject()
    counts <- intGuiCorpusCounts()
    texts <- integrationGuiTexts()
    bus <- intGuiImported(project$paths)
    spec <- integrationMofaParameterSpec()
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = project$paths), {
        preview <- session$returned()
        views <- shiny::isolate(bus$integration_import$views)
        # The parameters as resolved, never set by the tab: the provisional defaults.
        for (name in c("technical_replicates", "internal_standards", "missing_samples")) {
            expect_identical(preview$parameters[[name]], spec[[name]]$default)
        }
        expect_identical(preview$parameters$views, views)
        everyone <- unique(unlist(counts$samples))
        shared <- Reduce(intersect, counts$samples)
        expect_identical(preview$all_samples, length(everyone))
        expect_identical(preview$aligned$shared_samples, length(shared))
        facts <- integrationAlignmentFacts(preview)
        expect_identical(facts[[1L]], sprintf(texts[["fact_samples"]], length(everyone), length(shared)))
        expect_identical(facts[[2L]], sprintf(texts[["fact_missing_keep_missing"]], length(everyone) - length(shared)))
        for (view in views) {
            omic <- counts$views$omic[counts$views$view == view]
            per_sample <- table(counts$samples[[omic]])
            expected <- if (any(per_sample > 1L)) {
                sprintf(texts[["fact_replicates_average"]], view, length(counts$samples[[omic]]), length(per_sample),
                    sum(per_sample > 1L))
            } else {
                sprintf(texts[["fact_replicates_single"]], view, length(counts$samples[[omic]]))
            }
            expect_true(expected %in% facts, info = view)
        }
        expect_false(isTRUE(preview$few_samples))
        expect_match(output$facts$html, texts[["align_provisional"]], fixed = TRUE)
        status <- output$alignment_status$html
        expect_match(status, sprintf(texts[["align_key"]], "Metabolomics", "tech_rep_group"), fixed = TRUE)
        # One row per biological sample; the sample missing from metabolomics is shown as missing there (D4).
        table <- integrationAlignmentSampleTable(preview)
        expect_identical(table[[texts[["col_sample"]]]], sort(everyone, method = "radix"))
        absent <- setdiff(everyone, counts$samples$metabolomics)
        expect_identical(table[table[[1L]] %in% absent, "metabolomics_LCMS_Pos"], rep(texts[["cell_missing"]], length(absent)))
        metabolomics_views <- sum(counts$views$omic[match(views, counts$views$view)] == "metabolomics")
        expect_identical(table[table[[1L]] %in% absent, texts[["col_views_present"]]],
            rep(length(views) - metabolomics_views, length(absent)))
        lipid_runs <- table[table[[1L]] == "WT_7", "lipidomics_LCMS_Pos"]
        expect_identical(length(strsplit(lipid_runs, ", ", fixed = TRUE)[[1L]]), sum(counts$samples$lipidomics == "WT_7"))
        # Not confirmed until the user confirms; then saved as the integration design.
        expect_match(output$confirm_status$html, texts[["not_confirmed"]], fixed = TRUE)
        session$setInputs(confirm = 1)
        confirmed <- bus$integration_alignment
        expect_identical(confirmed$samples, preview$aligned$samples)
        expect_identical(confirmed$parameters, preview$parameters)
        expect_identical(confirmed$manifests, bus$integration_import$manifests)
        expect_identical(confirmed$file, file.path(project$paths$multiomic_inputs_dir, "integration_design.tsv"))
        saved <- utils::read.delim(confirmed$file, stringsAsFactors = FALSE, na.strings = "NA")
        expect_identical(saved$biological_sample, preview$aligned$samples$biological_sample)
        expect_identical(names(saved), names(preview$aligned$samples))
        expect_identical(bus$tab_status$integration_alignment, "complete")
        expect_match(output$confirm_status$html, "Sample table confirmed at", fixed = TRUE)
        # A new import (on its own store, as the Import tab commits it) clears the confirmation: it is for the import it
        # was made on.
        again <- newIntegrationStateOwner(file.path(project$paths$source_dir, "store-again"))
        withr::defer(closeWorkflowStateOwner(again))
        bus$state_manager <- again
        bus$integration_import <- integrationCommitImport(again,
            integrationImportRequest(integrationScanProject(project$paths$integration_dir), views))
        session$flushReact()
        expect_null(bus$integration_alignment)
        expect_match(output$confirm_status$html, texts[["not_confirmed"]], fixed = TRUE)
    })
})

test_that("KIT-139: internal standards, features seen once and few samples are stated with their counts", {
    project <- intGuiSmallProject()
    texts <- integrationGuiTexts()
    bus <- intGuiImported(project$paths)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = project$paths), {
        preview <- session$returned()
        facts <- integrationAlignmentFacts(preview)
        expect_identical(facts, c(sprintf(texts[["fact_samples"]], 5L, 4L), sprintf(texts[["fact_missing_keep_missing"]], 1L),
            sprintf(texts[["fact_replicates_average"]], "lipidomics_LCMS_Pos", 8L, 4L, 4L),
            sprintf(texts[["fact_standards_drop"]], "lipidomics_LCMS_Pos", 1L),
            sprintf(texts[["fact_replicates_single"]], "proteomics", 5L),
            sprintf(texts[["fact_sparse"]], "proteomics", 1L)))
        expect_true(preview$few_samples)
        expect_match(output$alignment_status$html, sprintf(texts[["few_samples"]], 4L, INTEGRATION_FEW_SAMPLES), fixed = TRUE)
        views <- integrationAlignmentViewTable(preview)
        expect_identical(views[[texts[["col_standards_removed"]]]], c(1L, 0L))
        expect_identical(views[[texts[["col_features_used"]]]], c(2L, 3L))
        # A warning (D7), not a refusal: the table can be confirmed.
        session$setInputs(confirm = 1)
        expect_identical(bus$tab_status$integration_alignment, "complete")
    })
})

test_that("KIT-139: no shared sample, fewer than two data sets and a changed input are refused in words", {
    texts <- integrationGuiTexts()
    project <- intGuiSmallProject(lipid_biological = c("C_1", "C_1", "C_2", "C_2", "D_1", "D_1", "D_2", "D_2"))
    bus <- intGuiImported(project$paths)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = project$paths), {
        preview <- session$returned()
        expect_s3_class(preview, "integration_alignment_refused")
        expected <- sprintf(texts[["align_refused"]], sprintf(texts[["refuse_no_shared"]], "lipidomics_LCMS_Pos, proteomics"))
        expect_identical(preview$refusal, expected)
        expect_match(output$alignment_status$html, expected, fixed = TRUE)
        session$setInputs(confirm = 1)
        expect_null(bus$integration_alignment)
        expect_identical(bus$tab_status$integration_alignment, "pending")
        expect_match(output$confirm_status$html, expected, fixed = TRUE)
    })
    project <- intGuiSmallProject()
    bus <- intGuiImported(project$paths)
    object <- integrationImportedObject(shiny::isolate(bus$state_manager))
    refusal <- tryCatch(integrationAlignmentPreview(object, "proteomics"), integration_refused = function(error) error)
    expect_s3_class(refusal, "integration_too_few_views")
    expect_identical(integrationAlignmentRefusal(refusal, "proteomics"),
        sprintf(texts[["align_refused"]], sprintf(texts[["refuse_too_few_views"]], 1L)))
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = project$paths), {
        expect_false(inherits(session$returned(), "integration_alignment_refused"))
        # The proteomics input saved again after the import: the confirmation is refused, nothing is saved.
        root <- project$root
        biological <- c("A_1", "A_2", "A_3", "B_1", "B_2")
        values <- rbind(P1 = c(10, 11, 12, 13, 14), P2 = c(5, 7, 6, 9, 8), P3 = c(1, 1.5, 2, 2.5, 3.5),
            P4 = c(NA, NA, 4, NA, NA)) + 1
        colnames(values) <- paste0("p", seq_along(biological))
        intGuiInput(root, "proteomics", list(proteomics = values), intGuiSamples(colnames(values), biological))
        session$setInputs(confirm = 1)
        expect_null(bus$integration_alignment)
        expect_match(output$confirm_status$html, "A data set has changed since it was imported", fixed = TRUE)
        expect_false(file.exists(file.path(project$paths$multiomic_inputs_dir, "integration_design.tsv")))
        scan <- integrationScanProject(project$paths$integration_dir, bus$integration_import$manifests)
        expect_identical(scan$omics$status[scan$omics$omic == "proteomics"], "changed")
    })
})
