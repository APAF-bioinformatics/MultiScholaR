# Extracted from test-workflow-native-session.R:466

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
paths <- makeDIANativeGuiPaths(withr::local_tempdir())
session <- makeDIANativeFixtureSession()
fixture <- testthat::test_path("..", "testdata", "e2e", "prot_lfq")
memory <- new.env(parent = emptyenv())
memory$state_manager <- WorkflowState$new(audit_enabled = FALSE)
memory$config_list <- list(globalParameters = list(workflow_type = "LFQ"))
imported <- importMaxQuantData(file.path(fixture, "proteinGroups.txt"))
applyProtImportResultToWorkflow(memory, imported, "maxquant", file.path(fixture, "proteins.fasta"),
        logInfo = function(...) NULL, showNotification = function(...) NULL)
memory$design_matrix <- data.frame(Run = c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"),
        group = rep(c("WT", "KO"), each = 3L), replicates = rep(1:3, 2L), stringsAsFactors = FALSE)
expect_identical(buildProtDesignStateCheckpoint(memory, "LFQ", "Import"), "protein_s4_initial")
