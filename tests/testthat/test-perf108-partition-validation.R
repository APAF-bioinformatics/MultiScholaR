test_that("partitioned admission preserves the conservative estimate and declared identities", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_STRUCTURAL") != "V")
    root <- tempfile("partition-contract-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    store <- newWorkflowEntityStore(root, "partition-contract", "DIA")
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    identity <- newSavedWorkflowScientificIdentity(input = "partition-contract", design = "small",
        parameters = "unchanged", selectors = "all", software = "native-api-v1")
    input <- saveWorkflowPrecursorRollupInput(store, object, "input", identity)
    scope <- newDIARelationalScope(store, input, profile = "P14", native = TRUE)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    validateDIARelationalInput(scope, input, NULL)
    original <- withr::with_options(list(multischolar.perf108.partition_validation = FALSE),
        estimateDIARollupWorkingMemory(scope, input, 16 * 1024^3))
    expect_identical(estimateDIARollupWorkingMemory(scope, input, 16 * 1024^3), original)
    con <- scope$connection
    workflowQueryExecute(con, paste0("CREATE TEMP TABLE generic_measurements AS SELECT ",
        "range % 2 AS specimen, range AS feature FROM range(20005)"))
    expect_false(hasWorkflowPartitionDuplicates(scope, "generic_measurements", "specimen",
        c("specimen", "feature")))
    workflowQueryExecute(con, "INSERT INTO generic_measurements VALUES (0, 0)")
    expect_true(hasWorkflowPartitionDuplicates(scope, "generic_measurements", "specimen",
        c("specimen", "feature")))
    expect_error(hasWorkflowPartitionDuplicates(scope, "generic_measurements", "specimen",
        "feature"), class = "workflow_execution_unsupported")
    workflowQueryExecute(con, "DELETE FROM generic_measurements")
    expect_false(hasWorkflowPartitionDuplicates(scope, "generic_measurements", "specimen",
        c("specimen", "feature")))
    predicates <- list()
    visitWorkflowCompletePartitions(scope, "generic_measurements", "specimen", function(sql, params) {
        predicates[[length(predicates) + 1L]] <<- list(sql, params)
    })
    expect_identical(predicates, list(list("FALSE", list())))
})
