# Extracted from test-workflow-native-session.R:426

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

# test -------------------------------------------------------------------------
withr::local_options(multischolar.dia_native_sessions = TRUE)
paths <- makeDIANativeGuiPaths(withr::local_tempdir())
fixture <- makeDIANativeSessionFixture(withr::local_tempdir())
withr::defer(closeDIAStateAdapter(fixture$manager))
legacy <- c(fixture$metadata, list(r6_complete_states = fixture$manager$states,
        r6_state_history = fixture$manager$getHistory(), r6_current_state_name = fixture$manager$current_state,
        current_s4_object = fixture$final, correlation_filtered_s4 = fixture$final))
workflow <- new.env(parent = emptyenv())
workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
previous <- workflow$state_manager
bad <- legacy
bad$r6_complete_states$raw_data_s4$data <- fixture$manager$.context$snapshot$states$raw_data_s4
session <- makeDIANativeFixtureSession()
expect_error(ingestDIANativeLegacySession(bad, workflow, session, paths$source_dir, paths),
        class = "saved_entity_corrupt")
expect_identical(workflow$state_manager, previous)
expect_true(ingestDIANativeLegacySession(legacy, workflow, session, paths$source_dir, paths))
