# KIT-140: one real MOFA2 training through the Model tab's server logic, on the KIT-137 tier-1 corpus. The corpus is
# copied into a project's integration folder, imported and matched through the Import and Sample Matching tabs, then
# trained by the Model tab's run handler at the table's defaults (double precision, INT-6). Skipped where MOFA2 or a
# Python with mofapy2 is absent (the package's own CI); run on the MOFA2 lane (branch run-mofa-smoke-k140) with the
# bundle's Python. Declared before the run:
#   - the tab's model equals `runIntegrationModel` called directly on its own store of the same import (identical: one
#     machine, double precision);
#   - its factors, weights and variance explained equal the stored double-precision result
#     (tests/testdata/integration/mofa_reference) within 1e-5 relative, the tolerance of test-integration-mofa-model.R;
#   - the same training as a managed-worker action (the worker's entry point and the controller's acceptance in this
#     process, helper-workflow-omics-worker.R) gives the identical model.

intTabRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
intTabCorpus <- function() file.path(intTabRoot(), "tests", "testdata", "integration", "corpus")
intTabReference <- function() file.path(intTabRoot(), "tests", "testdata", "integration", "mofa_reference")
intTabOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")
intTabFields <- c("factors", "weights", "variance_explained", "variance_total", "elbo", "samples", "view_summary",
    "shared_samples", "notes", "versions")

intTabSkip <- function() {
    testthat::skip_if_not(dir.exists(intTabCorpus()), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intTabProject <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = file.path(root, "results", "integration"),
        multiomic_inputs_dir = file.path(root, "results", "integration", "multiomic_inputs"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intTabOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intTabCorpus(), omic), full.names = TRUE), dir)))
    }
    paths
}

intTabNoProgress <- function(message, value, expr) force(expr)

# The project imported through the Import tab and its sample table confirmed through the Sample Matching tab.
intTabConfirmed <- function(paths, .local_envir = parent.frame()) {
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)),
        envir = .local_envir)
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(views = bus$integration_scan$views$view, import = 1)
    })
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = paths), {
        session$setInputs(confirm = 1)
    })
    stopifnot(!is.null(shiny::isolate(bus$integration_alignment)))
    bus
}

intTabDirect <- function(paths, views, .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0("direct-", basename(tempfile()))))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    integrationCommitImport(manager, integrationImportRequest(integrationScanProject(paths$integration_dir), views))
    manager
}

test_that("KIT-140: MOFA2 trained through the Model tab equals the model step and the stored double-precision result", {
    intTabSkip()
    paths <- intTabProject()
    texts <- integrationGuiTexts()
    bus <- intTabConfirmed(paths)
    started <- Sys.time()
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths,
        withProgressFn = intTabNoProgress), {
        session$flushReact()
        expect_match(output$readiness$html, texts[["readiness_ok"]], fixed = TRUE)
        session$setInputs(run = 1)
        model <- bus$integration_model
        expect_identical(model$status, "computed", info = output$run_status$html)
        expect_identical(model$parameters$precision, "double")
        value <- session$returned()$value
        message(sprintf("KIT-140 tab training: %d iterations, %d factors, %.1f s; model sha256 %s",
            value$elbo$iterations, ncol(value$factors) - 1L, as.numeric(difftime(Sys.time(), started, units = "secs")),
            value$model_file$sha256))
        expect_identical(value$versions$MOFA2, "1.22.0")
        expect_identical(value$versions$mofapy2, "0.7.3")
        expect_match(output$result$html, sprintf(texts[["result_factors"]], 5L, 5L), fixed = TRUE)
        expect_match(output$result$html, sprintf(texts[["result_saved"]], value$model_file$sha256), fixed = TRUE)
        # The model step called directly on its own store of the same import.
        direct <- runIntegrationModel(intTabDirect(paths, bus$integration_import$views), bus$integration_alignment$parameters)
        for (name in intTabFields) expect_identical(value[[name]], direct$value[[name]], info = name)
        # The stored double-precision result (KIT-062, re-recorded at 596e25e2) within 1e-5 relative.
        tables <- list(factors = value$factors, weights = value$weights, variance_explained = value$variance_explained)
        for (name in names(tables)) {
            stored <- utils::read.delim(file.path(intTabReference(), paste0(name, ".tsv")), check.names = FALSE,
                stringsAsFactors = FALSE)
            numeric <- names(stored)[vapply(stored, is.numeric, logical(1))]
            if (identical(dim(stored), dim(tables[[name]]))) {
                a <- as.matrix(tables[[name]][numeric])
                b <- as.matrix(stored[numeric])
                message(sprintf("KIT-140 stored-result comparison, %s: largest relative difference %.3g",
                    name, max(abs(a - b)) / max(abs(b))))
            }
            expect_equal(tables[[name]], stored, tolerance = 1e-5, info = name)
        }
    })
})

test_that("KIT-140: MOFA2 trained as a managed-worker action through the Model tab gives the identical model", {
    intTabSkip()
    testthat::skip_if_not(exists("omicsWorkerAction", mode = "function"), "helper-workflow-omics-worker.R is not loaded")
    paths <- intTabProject()
    bus <- intTabConfirmed(paths)
    manager <- shiny::isolate(bus$state_manager)
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
    fixture <- list(owner = manager, source_dir = paths$source_dir, paths = paths)
    controller <- omicsWorkerController(fixture)
    runner <- list(run = function(actions, done = NULL, failed = NULL, dropped = NULL) {
        for (action in actions) {
            result <- tryCatch(omicsWorkerAction(fixture, controller, action$operation, action$parameters),
                error = function(error) error)
            if (inherits(result, "error")) {
                if (is.function(failed)) failed(result)
                return(invisible(NULL))
            }
            withWorkflowStateScope(manager, function() action$deliver(result))
        }
        if (is.function(done)) done()
    })
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) {
        if (identical(omic_type, "integration")) runner else NULL
    }, .package = "MultiScholaR")
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths), {
        session$setInputs(run = 1)
        model <- bus$integration_model
        expect_identical(model$status, "computed", info = output$run_status$html)
        value <- session$returned()$value
        direct <- runIntegrationModel(intTabDirect(paths, bus$integration_import$views), bus$integration_alignment$parameters)
        for (name in intTabFields) expect_identical(value[[name]], direct$value[[name]], info = name)
        expect_identical(manager$.context$snapshot$exports$integration_mofa_model$reference, model$reference)
    })
})
