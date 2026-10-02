test_that("the GUI pair comparer runs in its isolated environment and preserves metadata", {
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    expect_false(exists("object.size", driver, inherits = TRUE))
    root <- withr::local_tempdir()
    control <- file.path(root, "control")
    candidate <- file.path(root, "candidate")
    report <- file.path(root, "report.tsv")
    writeLines("small comparer unit fixture; no scientific execution", report)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    object@protein_quant_table[[2L]][1L] <- NA_real_
    metadata <- data.frame(group = ordered(c("B", "A"), levels = c("B", "A", "unused")))
    attr(metadata, "audit_label") <- "retained"
    object@args$verifier_fixture <- list(first = NULL, metadata = metadata,
        annotated = structure(c("x", "y"), units = "fixture"), last = NULL)
    object@args$removeProteinsWithOnlyOneReplicateFn <- list(
        core_utilisation = newDIAWorkerSpec(4L), grouping_variable = "group")
    records <- list(initial = list(data = NULL, config = NULL),
        protein = list(data = object, config = list(optional = NULL, values = metadata),
            timestamp = as.POSIXct("2026-09-21", tz = "UTC")))
    session <- list(r6_complete_states = records, r6_state_history = as.list(names(records)),
        r6_current_state_name = "protein")
    for (directory in c(control, candidate)) {
        source <- file.path(directory, "dia_baseline_smoke", "scripts", "proteomics")
        dir.create(source, recursive = TRUE)
        saveRDS(session, file.path(source, "filtered_session_data_latest.rds"))
        jsonlite::write_json(list(capacity = list(treatment = if (directory == control) "resident" else "qs2_duckdb"),
            input_sha256 = list(report = "unit-fixture"), inputs = list(report = report),
            installation = list(package_path = find.package("MultiScholaR"))),
            file.path(directory, "smoke-receipt.json"), auto_unbox = TRUE)
        png::writePNG(matrix(c(0, 1, 1, 0), 2), file.path(directory, "dia_baseline_smoke", "plot.png"))
    }
    session$r6_complete_states$protein$data@args$removeProteinsWithOnlyOneReplicateFn$core_utilisation <- NA
    saveRDS(session, file.path(candidate, "dia_baseline_smoke", "scripts", "proteomics",
        "filtered_session_data_latest.rds"))
    before <- serialize(session, NULL)
    result <- driver$compareCheckpointGuiCapacityMembers(control, candidate, session)
    expect_identical(result$status, "passed_small_gui_science")
    expect_true(result$png_contents_identical)
    expect_gt(result$sizes$protein$s4_object_size_bytes, 0)
    expect_gt(result$sizes$protein$table_object_size_bytes, 0)
    expect_identical(serialize(session, NULL), before)

    changed <- session
    changed$r6_complete_states$protein$data@args$removeProteinsWithOnlyOneReplicateFn$grouping_variable <- "condition"
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, changed), "state differs")
    changed <- session
    changed$r6_complete_states$protein$data@protein_quant_table[[2L]][2L] <- 1000
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, changed), "state differs")
    changed <- session
    attr(changed$r6_complete_states$protein$data@args$verifier_fixture$metadata, "audit_label") <- "changed"
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, changed), "state differs")
    changed <- session
    attr(changed$r6_complete_states$protein$data@args$verifier_fixture$annotated, "units") <- "changed"
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, changed), "state differs")
    changed <- session
    changed$r6_complete_states$protein$data@args$verifier_fixture$first <- 1L
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, changed), "state differs")
    changed <- session
    changed$r6_complete_states$protein$data@protein_quant_table[[2L]][1L] <- NaN
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, changed), "differs")
    png::writePNG(matrix(1, 2, 2), file.path(candidate, "dia_baseline_smoke", "plot.png"))
    expect_error(driver$compareCheckpointGuiCapacityMembers(control, candidate, session), "PNG outputs differ")
})

test_that("only the approved DA UI clock value is normalized in comparison copies", {
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    clock <- as.POSIXct("2026-09-23 07:00:00", tz = "UTC")
    object@args$daAnalysisUI <- list(q_value_threshold = 0.05,
        log_fold_change_cutoff = 0, formula_string = "~ 0 + group", timestamp = clock)
    object@args$unrelated <- list(timestamp = clock)
    original <- list(data = object, config = list(timestamp = clock), description = "retained")
    later <- original
    later$data@args$daAnalysisUI$timestamp <- clock + 3600
    before <- serialize(list(original, later), NULL)
    expected <- driver$comparableCheckpointGuiRecord(original, "fixture")
    expect_identical(driver$comparableCheckpointGuiRecord(later, "fixture"), expected)
    expect_identical(serialize(list(original, later), NULL), before)
    expect_identical(attributes(expected$data@args$daAnalysisUI$timestamp), attributes(clock))
    for (field in c("q_value_threshold", "log_fold_change_cutoff", "formula_string")) {
        changed <- later
        changed$data@args$daAnalysisUI[[field]] <- switch(field,
            q_value_threshold = 0.1, log_fold_change_cutoff = 1, formula_string = "~ group")
        expect_false(isTRUE(all.equal(driver$comparableCheckpointGuiRecord(changed, "fixture"),
            expected, tolerance = 1e-12)), info = field)
    }
    for (clock_value in list(NULL, NA_real_, as.numeric(clock), as.character(clock))) {
        changed <- later
        changed$data@args$daAnalysisUI$timestamp <- clock_value
        expect_false(identical(driver$comparableCheckpointGuiRecord(changed, "fixture"), expected))
    }
    changed <- later
    changed$data@args$unrelated$timestamp <- clock + 3600
    expect_false(identical(driver$comparableCheckpointGuiRecord(changed, "fixture"), expected))
    changed <- later
    changed$config$timestamp <- clock + 3600
    expect_false(identical(driver$comparableCheckpointGuiRecord(changed, "fixture"), expected))
    changed <- later
    changed$data@protein_quant_table[[2L]][1L] <- 1000
    expect_false(identical(driver$comparableCheckpointGuiRecord(changed, "fixture"), expected))
})
