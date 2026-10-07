test_that("PERF-107 accepts only its frozen GUI qualification binding", {
    config <- list(mode = "gui", qualification = list(
        programme = "perf107-phase-a",
        variant = "physical",
        entry = "resident",
        member_sha256 = strrep("a", 64L),
        deadline = 900
    ))
    expect_identical(workflowQualificationVariant(config,
        "minimum_peptides_per_sample"), "physical")
    changed <- config
    changed$qualification$entry <- "unknown"
    expect_error(workflowQualificationVariant(changed,
        "minimum_peptides_per_sample"), class = "workflow_execution_invalid")
    expect_error(workflowQualificationVariant(config, "precursor_rollup"),
        class = "workflow_execution_invalid")
})

test_that("PERF-107 acknowledgement handler is registered in the static app head", {
    rendered <- htmltools::renderTags(app_ui(NULL))
    application <- paste(rendered$head, rendered$html)
    module <- htmltools::renderTags(mod_prot_qc_peptide_sample_ui("sample"))$html
    marker <- "multischolar-perf107-ready"
    matches <- gregexpr(marker, application, fixed = TRUE)[[1L]]
    expect_length(matches, 1L)
    expect_gt(matches[[1L]], 0L)
    expect_match(application, "setTimeout(install, 25)", fixed = TRUE)
    expect_match(application, "MultiScholaRPerf107AckState")
    expect_match(paste(deparse(body(mod_prot_qc_peptide_sample_server)),
        collapse = "\n"), "browser acknowledgement timed out")
    expect_false(grepl(marker, module, fixed = TRUE))
})

test_that("PERF-107 stages complete resident history without changing its owner", {
    root <- tempfile("perf107-adoption-")
    dir.create(root)
    paths <- list(source_dir = root)
    native <- WorkflowState$new(audit_enabled = FALSE)
    first <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    second <- first
    second@peptide_data$Precursor.Quantity[[1L]] <- 42
    native$saveState("h0", first, list(version = 0L), "first")
    native$saveState("h1", second, list(version = 1L), "second")
    before <- list(history = native$state_history, current = native$current_state,
        states = digest::digest(native$states, algo = "sha256"),
        audits = digest::digest(native$audit_records, algo = "sha256"))
    manager <- stageDIAGUIReferenceOwner(native, paths)
    withr::defer(closeDIAStateAdapter(manager))
    after <- list(history = native$state_history, current = native$current_state,
        states = digest::digest(native$states, algo = "sha256"),
        audits = digest::digest(native$audit_records, algo = "sha256"))
    expect_identical(after, before)
    expect_identical(as.character(unlist(manager$getHistory())),
        as.character(unlist(native$getHistory())))
    expect_identical(manager$current_state, native$current_state)
    expect_true(usesDIAReferenceStates(manager))
    expect_identical(diaStateOwnership(manager)$full_s4_reads, 0L)
})

test_that("PERF-107 GUI reservation owns the complete action window", {
    root <- tempfile("perf107-reservation-")
    dir.create(root)
    application <- WorkflowExecutionApplication$new()
    application$config <- list(root = root, stem = "multischolarwfxaaaaaaaaaaaaaaaa")
    application$cgroup <- root
    application$pid <- Sys.getpid()
    manager <- WorkflowState$new(audit_enabled = FALSE)
    peak_path <- file.path(root, "peak")
    writeLines("0", peak_path)
    publication <- list(experiment_label = "fixture", experiment_paths = list(
        source_dir = root), omic_type = "proteomics")
    testthat::local_mocked_bindings(
        checkWorkflowApplication = function(...) invisible(TRUE),
        openWorkflowPeakWindow = function(...) file(peak_path, open = "r+"),
        closeWorkflowPeakWindow = function(con) {
            close(con)
            123
        },
        workflowCgroupAccounting = function(...) list(cpu_seconds = 1,
            read_bytes = 2, written_bytes = 3),
        workflowMemoryCounter = function(...) 456,
        .package = "MultiScholaR",
    )
    action <- beginWorkflowGUIAction(application, manager,
        "minimum_peptides_per_sample", publication)
    expect_identical(application$active, action$id)
    expect_error(beginWorkflowGUIAction(application, manager,
        "minimum_peptides_per_sample", publication),
        class = "workflow_execution_busy")
    receipt <- finishWorkflowGUIAction(action, "complete", list(
        request_id = action$id,
        content_id = strrep("b", 64L),
        plot_ready = TRUE,
        stable_seconds = 0.2
    ))
    expect_identical(receipt$status, "complete")
    expect_identical(receipt$peak_bytes, 123)
    expect_true(receipt$browser_acknowledged)
    expect_null(application$active)
    expect_true(file.exists(file.path(root,
        paste0(action$id, "-gui-completion.rds"))))
})
