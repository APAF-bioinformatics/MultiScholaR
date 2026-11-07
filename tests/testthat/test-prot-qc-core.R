# KIT-104 W3: proteomics' QC back-references #6 to #8. `.formatPeptideRollupSummary` and
# `.peptideImputationAuditSummary` moved, unchanged and with their names kept, beside
# buildDIAReferencePeptideHandoff (their caller); `.proteinCleanupActiveKey` beside
# computeDIAReferenceAccessionCleanup (its caller). The roxygen tags of #7 and #8 moved with them, so they stay
# exported. These blocks show where they live and that both callers' output is unchanged through them.

protQcCoreDefinitions <- c(
    .formatPeptideRollupSummary = "utils_workflow_dia_qc_rollup.R",
    .peptideImputationAuditSummary = "utils_workflow_dia_qc_rollup.R",
    .proteinCleanupActiveKey = "utils_workflow_prot_declarations.R"
)

protQcCoreSourceFile <- function(name) {
    ns <- asNamespace("MultiScholaR")
    if (!exists(name, envir = ns, inherits = FALSE)) return(NA_character_)
    file <- utils::getSrcFilename(get(name, envir = ns, inherits = FALSE))
    if (length(file)) file else NA_character_
}

newProtQcCoreManager <- function(label, envir = parent.frame()) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = envir),
        label, "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager
}

# A component checkpoint of the precursor fixture, as the reference path builds it
# (test-workflow-checkpoint-dia-state.R), carrying the step output's arguments and row count as the step's
# checkpoint does when it reaches the handoff.
protQcCoreCheckpoint <- function(manager, precursor, output) {
    reference <- saveDIATableComponent(manager$.context, precursor@peptide_data)
    checkpoint <- newDIAPeptideCheckpoint(manager, reference, precursor@design_matrix, precursor@args)
    checkpoint$slots$args <- output@args
    checkpoint$dimensions$rows <- nrow(output@peptide_data)
    checkpoint
}

# The handoff without its wall-clock stamp.
protQcCoreHandoffDigest <- function(handoff) {
    handoff$qc$timestamp <- NULL
    digest::digest(handoff, algo = "sha256")
}

# Counts the calls made by name to a namespace function while `code` runs.
protQcCoreCountCalls <- function(name, code) {
    original <- get(name, envir = asNamespace("MultiScholaR"))
    calls <- 0L
    counting <- function(...) {
        calls <<- calls + 1L
        original(...)
    }
    value <- rlang::inject(testthat::with_mocked_bindings(code, !!name := counting, .package = "MultiScholaR"))
    list(value = value, calls = calls)
}

# Recorded at b423b764 / 5f8cb531 (before W3) through buildDIAReferencePeptideHandoff and
# computeDIAReferenceAccessionCleanup on these inputs.
protQcCoreRollupDigest <- "7557865443651b81c5be8209e386123847b2cfe0318818e027205b5b6112a317"
protQcCoreImputationDigest <- "1c4118d8e289822a76cf451873b1b843619a689d03e8f24b68b4bfd1696b8f81"
protQcCoreCleanupDigest <- "a94568f6e3b648fa4b472282b80a30008e81bae6e21d2a82bc4ec03c337c7595"

test_that("back-references #6 to #8 are defined in the proteomics core files, with their names kept", {
    if (is.null(utils::getSrcref(buildDIAReferencePeptideHandoff))) {
        testthat::skip("source references are kept under load_all only")
    }
    found <- vapply(names(protQcCoreDefinitions), protQcCoreSourceFile, character(1))
    expect_identical(found, protQcCoreDefinitions)
    # The GUI side keeps its reader in the module file.
    expect_identical(protQcCoreSourceFile(".peptideRollupSummaryFromObject"), "mod_prot_qc_peptide_rollup.R")
})

test_that("#7 and #8 stay exported", {
    exports <- getNamespaceExports("MultiScholaR")
    expect_true(all(c(".peptideImputationAuditSummary", ".proteinCleanupActiveKey") %in% exports))
})

test_that("#6: the precursor rollup handoff is unchanged through the moved summary formatter", {
    manager <- newProtQcCoreManager("qc-core-rollup")
    withr::defer(closeDIAStateAdapter(manager))
    precursor <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    rolled <- rollUpPrecursorToPeptide(precursor, core_utilisation = NA)
    checkpoint <- protQcCoreCheckpoint(manager, precursor, rolled)
    runs <- unique(as.character(rolled@design_matrix$Run))
    counted <- protQcCoreCountCalls(".formatPeptideRollupSummary",
        buildDIAReferencePeptideHandoff(checkpoint, "rollUpPrecursorToPeptide", list(), 3L, runs, character(),
            "Rolled up precursors"))
    expect_identical(counted$calls, 1L)
    handoff <- counted$value
    expect_identical(handoff$qc_name, "precursor_rollup")
    expect_match(handoff$text, "Aggregation: ", fixed = TRUE)
    expect_identical(protQcCoreHandoffDigest(handoff), protQcCoreRollupDigest)
    # A legacy result without the audit keeps its wording.
    expect_identical(.formatPeptideRollupSummary(NULL), "Rollup audit: unavailable for this legacy result")
})

test_that("#7: the imputation handoff is unchanged through the moved audit summary", {
    manager <- newProtQcCoreManager("qc-core-imputation")
    withr::defer(closeDIAStateAdapter(manager))
    precursor <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    rolled <- filterMinNumPeptidesPerSample(rollUpPrecursorToPeptide(precursor, core_utilisation = NA),
        peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    rlang::local_bindings(config_list = list(peptideMissingValueImputation = list(proportion_missing_values = 0.5)),
        .env = .GlobalEnv)
    imputed <- peptideMissingValueImputation(rolled, proportion_missing_values = 0.5, core_utilisation = NA)
    checkpoint <- protQcCoreCheckpoint(manager, precursor, imputed)
    runs <- unique(as.character(imputed@design_matrix$Run))
    parameters <- list(proportion_missing_values = 0.5, exclusion_column = NULL)
    counted <- protQcCoreCountCalls(".peptideImputationAuditSummary",
        buildDIAReferencePeptideHandoff(checkpoint, "peptideMissingValueImputation", parameters, 3L, runs,
            character(), "Applied missing value imputation using technical replicates"))
    expect_identical(counted$calls, 1L)
    handoff <- counted$value
    expect_identical(handoff$qc_name, "imputation")
    expect_identical(handoff$config$imputation_summary,
        .peptideImputationAuditSummary(checkpoint, 0.5, NULL))
    expect_identical(protQcCoreHandoffDigest(handoff), protQcCoreImputationDigest)
})

test_that("#8: the accession cleanup is unchanged through the moved active-key resolver", {
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- newProtQcCoreManager("qc-core-cleanup")
    withr::defer(closeDIAStateAdapter(manager))
    manager$saveState("protein_s4_created", protein, list(), "fixture")
    now <- function() as.POSIXct("2026-11-08 00:00:00", tz = "UTC")
    counted <- protQcCoreCountCalls(".proteinCleanupActiveKey",
        computeDIAReferenceAccessionCleanup(manager, getDIAStateReference(manager),
            list(delimiter = ";", aggregation_method = "mean"), now_fn = now))
    expect_identical(counted$calls, 2L)
    value <- counted$value
    expect_identical(value$config$active_protein_key, protein@protein_id_column)
    expect_identical(digest::digest(value[c("state_name", "description", "config", "section", "protein_count",
        "execution", "text", "qc_name", "qc", "accession_cleanup_results", "cleanupApplied")], algo = "sha256"),
        protQcCoreCleanupDigest)

    # The resolver: the declared key when present, else the first known column, else a refusal.
    expect_identical(.proteinCleanupActiveKey(protein), protein@protein_id_column)
    inferred <- protein
    inferred@protein_id_column <- "absent"
    expect_identical(.proteinCleanupActiveKey(inferred),
        c("Protein.Group", "Protein.Ids", "Protein.IDs", "protein_id")[
            c("Protein.Group", "Protein.Ids", "Protein.IDs", "protein_id") %in% names(protein@protein_quant_table)][1])
    keyless <- inferred
    keyless@protein_quant_table <- keyless@protein_quant_table[, setdiff(names(keyless@protein_quant_table),
        c(protein@protein_id_column, "Protein.Group", "Protein.Ids", "Protein.IDs", "protein_id")), drop = FALSE]
    expect_error(.proteinCleanupActiveKey(keyless), "No active protein identity column is available.", fixed = TRUE)
})
