# Extracted from test-workflow-native-session.R:1349

# prequel ----------------------------------------------------------------------
makeDIANativeSessionFixture <- function(root) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "native-test", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$setWorkflowType("DIA")
    raw <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw_data_s4", raw, list(threshold = 0.01), "raw")
    final <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("inactive", final, list(), "retained branch")
    manager$revertToState("raw_data_s4")
    manager$saveState("correlation_filtered", final, list(), "final")
    list(manager = manager, raw = raw, final = final,
        metadata = list(design_matrix = final@design_matrix,
            contrasts_tbl = data.frame(contrasts = "B-A", friendly_names = "test"),
            config_list = final@args, final_protein_count = nrow(final@protein_quant_table),
            final_sample_count = nrow(final@design_matrix), workflow_type = "DIA",
            normalization_method = "cyclicloess", ruv_mode = "skip", ruv_applied = FALSE))
}
makeDIANativeGuiPaths <- function(root) {
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    for (path in c(unlist(paths), file.path(dirname(paths$protein_qc_dir), "publication_graphs"))) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    paths
}
makeDIANativeFixtureSession <- function() {
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
makeKit013CheckpointWorkflow <- function(native, object) {
    shiny::reactiveValues(state_manager = native, data_tbl = object@peptide_data,
        data_cln = object@peptide_data, data_type = "peptide", data_format = "diann",
        config_list = list(globalParameters = list(workflow_type = "DIA")),
        design_matrix = object@design_matrix,
        column_mapping = list(run_col = "Run", protein_col = "Protein.Group",
            peptide_col = "Stripped.Sequence", raw_quantity_col = "Precursor.Quantity",
            norm_quantity_col = "Precursor.Normalised"))
}
makeKit013CheckpointPaths <- function(root, name) {
    paths <- list(source_dir = file.path(root, paste0(name, "-source")),
        peptide_qc_dir = file.path(root, paste0(name, "-peptide")))
    lapply(paths, dir.create)
    paths
}

# test -------------------------------------------------------------------------
skip_if_not(workflowOwnerCapable())
root <- withr::local_tempdir()
source_dir <- file.path(root, "source")
dir.create(source_dir)
session <- new.env(parent = emptyenv())
session$userData <- new.env(parent = emptyenv())
session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
workflow <- list2env(list(state_update_trigger = NULL))
manager <- createDIASessionOwner(workflow, list(source_dir = source_dir), session)
withr::defer(closeDIAStateAdapter(manager))
expect_s3_class(manager, "DIAStateAdapter")
expect_true(isDIASessionOwnerPending(manager))
expect_identical(manager$.context$store$identity$project_id, "PROTEOMICS-pending")
expect_identical(manager$current_state, "initial")
expect_true(usesDIAReferenceStates(manager))
expect_identical(manager$.context$experiment_paths$source_dir, source_dir)
expect_identical(session$userData$dia_capacity_manager, manager)
workflow$state_manager <- manager
bound <- initializeDIANativeImport(workflow, list(source_dir = source_dir), session, "maxquant", data_type = "lfq")
expect_identical(bound, manager)
expect_false(isDIASessionOwnerPending(manager))
expect_identical(manager$.context$store$identity$workflow_id, "LFQ")
expect_identical(manager$.context$store$identity$project_id, workflowNativeStoreLabel("maxquant", "LFQ"))
expect_identical(jsonlite::read_json(file.path(manager$.context$store$directory, "store.json"))$workflow_id, "LFQ")
expect_identical(manager$.context$snapshot$identity, manager$.context$store$identity)
expect_identical(manager$workflow_type, "LFQ")
expect_identical(manager$getHistory(), list("initial"))
expect_error(initializeDIANativeImport(workflow, list(source_dir = source_dir), session, "diann", method = "iq"), class = "saved_entity_scope")
expect_identical(initializeDIANativeImport(workflow, list(source_dir = source_dir), session, "maxquant", data_type = "lfq"), manager)
