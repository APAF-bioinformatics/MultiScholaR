# KIT-103 W3: the proteomics design core. Back-reference #3 (buildDIAReferenceDesignState) and its callee
# buildProtDesignInitialState moved first, then the save payload, the rename edits and the state checkpoint,
# all with their names kept, beside transformDIADesignTable; the design step's UniProt retrieval and the
# imported-FASTA processing became core functions in the annotation files. These blocks show where they
# live, that the core caller of #3 works through it, and that the annotation step's output is unchanged.

protDesignCoreDefinitions <- c(
    buildDIAReferenceDesignState = "utils_workflow_dia_state.R",
    buildProtDesignInitialState = "utils_workflow_dia_state.R",
    buildProtDesignSaveResultsContrastsTable = "utils_workflow_dia_state.R",
    buildProtDesignSaveResultsPayload = "utils_workflow_dia_state.R",
    transformProtDesignSampleNames = "utils_workflow_dia_state.R",
    applyProtDesignBulkRenameUpdates = "utils_workflow_dia_state.R",
    applyProtDesignSingleRenameUpdate = "utils_workflow_dia_state.R",
    buildProtDesignStateCheckpoint = "utils_workflow_dia_state.R",
    retrieveProtDesignUniprotAnnotations = "func_prot_annotation_uniprot_full_helpers.R",
    processProtDesignImportedFasta = "func_prot_annotation_fasta_processing_helpers.R"
)

protDesignCoreSourceFile <- function(name) {
    ns <- asNamespace("MultiScholaR")
    if (!exists(name, envir = ns, inherits = FALSE)) return(NA_character_)
    file <- utils::getSrcFilename(get(name, envir = ns, inherits = FALSE))
    if (length(file)) file else NA_character_
}

newProtDesignCoreManager <- function(label, envir = parent.frame()) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = envir),
        label, "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager
}

# The table of test-workflow-checkpoint-dia-state.R's design-metadata block, with a protein column.
protDesignCoreTable <- function() {
    prepareDIAParserTable(readr::read_tsv(I(paste0("Run\tProtein.Group\tPrecursor.Normalised\n",
        "A\tP00001;P00002\t1\nB\tP00003\t2\nA\tP00001;P00002\t3\nB\tP00004\t4\n")), show_col_types = FALSE))
}

# Table content without parser pointers or row-name storage, stable across processes.
protDesignCoreDigest <- function(x) {
    digest::digest(list(names = names(x), class = class(x), cols = lapply(as.list(as.data.frame(x)), identity)),
        algo = "sha256")
}

# Recorded at 3c657c14 / 70af7bb4 (before W3) through completeProtDesignPostCheckpoint and
# hydrateProtDesignImportedFastaSidecar on these inputs.
protDesignCoreUniprotDigest <- "3c5b81689db4ad40de89571cfc27503e2422bdeb562c42dab62b230cdd584037"
protDesignCoreAaSeqDigest <- "1b29406fe149a49a084190dec6cf69db443e720205e0a09653489f4e4cfcdc7a"
protDesignCoreFastaMetadataDigest <- "51f5cd89ef9eaf3b515fd3b1c2fc72abfcf25b75e41c87e6a6e8a62f16551ab3"

protDesignCoreMetadataDigest <- function(metadata) {
    digest::digest(metadata[setdiff(names(metadata), "processing_timestamp")], algo = "sha256")
}

test_that("the moved design functions are defined in the proteomics core files, with their names kept", {
    if (is.null(utils::getSrcref(transformDIADesignTable))) {
        testthat::skip("source references are kept under load_all only")
    }
    found <- vapply(names(protDesignCoreDefinitions), protDesignCoreSourceFile, character(1))
    expect_identical(found, protDesignCoreDefinitions)
    # The GUI sides keep their names in the modules: the diagnostic traces them.
    expect_identical(protDesignCoreSourceFile("completeProtDesignPostCheckpoint"), "mod_prot_design_state_helpers.R")
    expect_identical(protDesignCoreSourceFile("hydrateProtDesignImportedFastaSidecar"),
        "mod_prot_design_import_helpers.R")
})

test_that("back-reference #3: prepareDIAGUIReferences builds the design state through the moved function", {
    manager <- newProtDesignCoreManager("design-references")
    withr::defer(closeDIAStateAdapter(manager))
    table <- protDesignCoreTable()
    mapping <- list(run_col = "Run", protein_col = "Protein.Group", quantity_col = "Precursor.Normalised")
    workflow_data <- list(data_tbl = table, data_cln = table, data_type = "peptide",
        column_mapping = mapping, config_list = list())

    calls <- list()
    original <- buildDIAReferenceDesignState
    references <- testthat::with_mocked_bindings(
        prepareDIAGUIReferences(manager, workflow_data),
        buildDIAReferenceDesignState = function(...) {
            calls[[length(calls) + 1L]] <<- list(...)
            original(...)
        })

    expect_length(calls, 1L)
    expect_identical(calls[[1L]][[2L]], references$data_tbl)
    initial <- manager$.context$snapshot$design_initial
    expect_identical(initial$source_reference, references$data_tbl)
    expect_identical(initial$column_mapping, mapping)
    expect_identical(manager$.context$snapshot$file_dependencies$design,
        list(role = "design_builder", entity_id = initial$data_reference$entity_id))

    # The settled design state: the samples from the import metadata, in mixedsort order, and a data
    # reference in place of an assay.
    state <- buildDIAReferenceDesignState(manager, references$data_tbl, list(), mapping)
    expect_null(state$data_cln)
    expect_true("data_cln" %in% names(state))
    expect_identical(state$data_reference, initial$data_reference)
    expect_identical(state$sample_metadata, getDIAImportMetadata(manager)$samples)
    expect_identical(state$design_matrix$Run, c("A", "B"))
    expect_identical(names(state$design_matrix),
        c("Run", "group", "factor1", "factor2", "factor3", "batch", "replicates", "tech_reps"))
})

test_that("run order: a fresh proteomics design is mixedsort on both initial-state paths", {
    runs <- c("S10", "S2", "S1", "S2")
    plain <- buildProtDesignInitialState(data.frame(Run = runs, Quantity = c("1", "2", "3", "4")),
        list(deAnalysisParameters = list(formula_string = "~ 0 + group")), list(quantity_col = "Quantity"))
    expect_identical(plain$design_matrix$Run, c("S1", "S2", "S10"))
    expect_identical(plain$data_cln$Quantity, c(1, 2, 3, 4))
    expect_identical(plain$formula, "~ 0 + group")

    manager <- newProtDesignCoreManager("design-run-order")
    withr::defer(closeDIAStateAdapter(manager))
    table <- prepareDIAParserTable(readr::read_tsv(I("Run\tPrecursor.Normalised\nS10\t1\nS2\t2\nS1\t3\n"),
        show_col_types = FALSE))
    mapping <- list(run_col = "Run", quantity_col = "Precursor.Normalised")
    references <- prepareDIAGUIReferences(manager, list(data_tbl = table, data_cln = table,
        data_type = "peptide", column_mapping = mapping, config_list = list()))
    state <- buildDIAReferenceDesignState(manager, references$data_tbl, list(), mapping)
    expect_identical(state$design_matrix$Run, c("S1", "S2", "S10"))
})

test_that("the design step's UniProt annotation is unchanged through the offline UniProt double", {
    manager <- newProtDesignCoreManager("design-uniprot")
    withr::defer(closeDIAStateAdapter(manager))
    source <- saveDIATableComponent(manager$.context, protDesignCoreTable())
    root <- withr::local_tempdir()
    paths <- list(results_dir = file.path(root, "results"), source_dir = file.path(root, "scripts"))
    dir.create(paths$results_dir)
    dir.create(paths$source_dir)
    guiEntryInstallUniprotCache(paths)

    # The GUI side (modal, progress, notifications) over the core computation.
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- manager
    workflow_data$data_cln <- source
    workflow_data$column_mapping <- list(protein_col = "Protein.Group")
    workflow_data$taxon_id <- 9606
    workflow_data$tab_status <- list(design_matrix = "pending")
    session <- shiny::MockShinySession$new()
    utils::capture.output(suppressMessages(shiny::withReactiveDomain(session,
        completeProtDesignPostCheckpoint(workflow_data, paths, session, successMessage = "saved"))))
    expect_identical(protDesignCoreDigest(workflow_data$uniprot_dat_cln), protDesignCoreUniprotDigest)
    expect_identical(protDesignCoreDigest(readRDS(file.path(paths$source_dir, "uniprot_dat_cln.RDS"))),
        protDesignCoreUniprotDigest)
    expect_identical(workflow_data$tab_status$design_matrix, "complete")

    # The core function alone: the same table, the same cache directory and scripts copy.
    unlink(file.path(paths$source_dir, "uniprot_dat_cln.RDS"))
    progress <- list()
    annotations <- NULL
    utils::capture.output(annotations <- suppressMessages(retrieveProtDesignUniprotAnnotations(
        manager, source, "Protein.Group", 9606, resultsDir = paths$results_dir, sourceDir = paths$source_dir,
        progressCallback = function(current, total) progress[[length(progress) + 1L]] <<- c(current, total))))
    expect_identical(protDesignCoreDigest(annotations), protDesignCoreUniprotDigest)
    expect_identical(protDesignCoreDigest(readRDS(file.path(paths$source_dir, "uniprot_dat_cln.RDS"))),
        protDesignCoreUniprotDigest)
    expect_true(file.exists(file.path(paths$results_dir, "cache", "uniprot_annotations", "uniprot_annotations.RDS")))
    expect_error(retrieveProtDesignUniprotAnnotations(manager, source, NULL, 9606),
        "Protein column not found in column_mapping. Cannot get annotations.", fixed = TRUE)
})

test_that("the imported-FASTA annotation step is unchanged", {
    root <- withr::local_tempdir()
    import_dir <- file.path(root, "import")
    source_dir <- file.path(root, "scripts")
    results_dir <- file.path(root, "results")
    for (path in c(import_dir, source_dir, results_dir)) dir.create(path)
    fasta <- testthat::test_path("..", "testdata", "e2e", "prot_tmt", "proteins.fasta")

    # The GUI side over the core computation: no sidecar, so the FASTA is processed.
    workflow_data <- new.env(parent = emptyenv())
    session <- shiny::MockShinySession$new()
    suppressMessages(shiny::withReactiveDomain(session,
        hydrateProtDesignImportedFastaSidecar(workflow_data, import_dir, source_dir, resultsDir = results_dir,
            fastaPath = fasta, organismName = "Homo sapiens")))
    expect_identical(protDesignCoreDigest(workflow_data$aa_seq_tbl_final), protDesignCoreAaSeqDigest)
    expect_identical(protDesignCoreMetadataDigest(workflow_data$fasta_metadata), protDesignCoreFastaMetadataDigest)
    expect_identical(protDesignCoreDigest(readRDS(file.path(source_dir, "aa_seq_tbl_final.RDS"))),
        protDesignCoreAaSeqDigest)
    expect_identical(protDesignCoreDigest(readRDS(file.path(results_dir, "cache", "aa_seq_tbl.RDS"))),
        protDesignCoreAaSeqDigest)

    # The core function alone: the cache under results_dir, the same tables.
    unlink(file.path(results_dir, "cache"), recursive = TRUE)
    result <- processProtDesignImportedFasta(fasta, resultsDir = results_dir, organismName = "Homo sapiens")
    expect_identical(protDesignCoreDigest(result$aa_seq_tbl_final), protDesignCoreAaSeqDigest)
    expect_identical(protDesignCoreMetadataDigest(result$fasta_metadata), protDesignCoreFastaMetadataDigest)
    expect_identical(protDesignCoreDigest(readRDS(file.path(results_dir, "cache", "aa_seq_tbl.RDS"))),
        protDesignCoreAaSeqDigest)
})
