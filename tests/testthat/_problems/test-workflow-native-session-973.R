# Extracted from test-workflow-native-session.R:973

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
skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
root <- withr::local_tempdir()
output <- withr::local_tempdir()
restored_root <- withr::local_tempdir()
fixture <- makeDIANativeSessionFixture(root)
manager <- fixture$manager
withr::defer(closeDIAStateAdapter(manager))
report <- file.path(root, "report.tsv")
fasta <- file.path(root, "proteins.fasta")
design <- file.path(root, "design_matrix.tab")
writeLines(c("Run\tProtein.Group", "s1\tP1"), report)
writeLines(c(">sp|P1|ONE", "MPEPTIDE"), fasta)
writeLines(c("Run\tgroup", "s1\tA"), design)
legacy <- file.path(root, "iq_output_file.txt")
writeLines("legacy", legacy)
manager$.context$reference_states <- TRUE
registerDIAInputFiles(manager, c(diann_report = report, fasta = fasta))
completeDIAInputFiles(manager)
initial <- buildDIAReferenceDesignState(manager, saveDIAOwnedValue(data.frame(Run = "s1", Quantity = 1), "input", manager),
        list(), list(run_col = "Run", quantity_col = "Quantity"))
