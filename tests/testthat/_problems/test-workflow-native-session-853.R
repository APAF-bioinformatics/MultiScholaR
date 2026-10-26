# Extracted from test-workflow-native-session.R:853

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
root <- withr::local_tempdir()
paths <- makeDIANativeGuiPaths(root)
output <- withr::local_tempdir()
fixture <- makeDIANativeSessionFixture(paths$source_dir)
manager <- fixture$manager
withr::defer(closeDIAStateAdapter(manager))
manager$.context$session_format <- "native_v1"
peptide <- fixture$raw
keys <- c(peptide@protein_id_column, peptide@sample_id, peptide@peptide_sequence_column)
peptide@peptide_data <- peptide@peptide_data[!duplicated(peptide@peptide_data[keys]), ]
peptide@peptide_data$Peptide.Imputed <- 1000 + seq_len(nrow(peptide@peptide_data))
manager$saveState("imputed", peptide, list(), "synthetic post-imputation IQ input")
workflow <- list2env(list(state_manager = manager), parent = emptyenv())
apply_iq <- function(workflow, paths) withDIAStateScope(workflow$state_manager,
        function() runProteinIqRollupApplyStep(workflow, paths,
            captureCheckpointFn = function(...) NULL, showNotificationFn = function(...) NULL))
