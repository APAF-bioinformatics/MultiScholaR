# KIT-102 W3: the proteomics import core (scoping P1 to P6) and the Spectronaut and PD-TMT reader parity.
# The moved functions keep their names; these blocks show where they now live, that the core callers of
# back-references #1 and #2 work through them, and that the readers convert text-numeric quantities.

protImportCoreDefinitions <- c(
    readProtImportParquetSchema = "func_prot_import_detection.R",
    readProtImportHeaders = "func_prot_import_detection.R",
    importProteomicsData = "func_prot_import_readers.R",
    sanitizeProtImportTable = "utils_workflow_dia_sessions.R",
    describeProtImportData = "utils_workflow_dia_sessions.R",
    recordDIAImportMetadata = "utils_workflow_dia_sessions.R",
    getDIAImportMetadata = "utils_workflow_dia_sessions.R",
    processProtImportFastaData = "func_prot_annotation_fasta_processing_helpers.R",
    prepareProtImportOrganismSelectionChoices = "func_prot_import_organisms.R",
    filterProtImportDataToOrganism = "func_prot_import_organisms.R",
    readProtGuiConfig = "func_prot_import_detection.R",
    loadProtImportConfiguration = "func_prot_import_detection.R",
    readProtImportOptionalMapping = "func_prot_import_detection.R",
    loadProtImportOptionalMappings = "func_prot_import_detection.R"
)

protImportCoreSourceFile <- function(name) {
    ns <- asNamespace("MultiScholaR")
    if (!exists(name, envir = ns, inherits = FALSE)) return(NA_character_)
    file <- utils::getSrcFilename(get(name, envir = ns, inherits = FALSE))
    if (length(file)) file else NA_character_
}

newProtImportCoreManager <- function(label) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = parent.frame()),
        label, "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager
}

protImportCoreTable <- function() {
    prepareDIAParserTable(readr::read_tsv(I("Run\tPrecursor.Normalised\nA\t1\nB\t2\nA\t3\n"),
        show_col_types = FALSE))
}

# Reader output without parser pointers, for a digest that is stable across processes.
protImportReaderDigest <- function(result) {
    digest::digest(list(data = lapply(as.list(result$data), identity), names = names(result$data),
        class = class(result$data), data_type = result$data_type,
        column_mapping = result$column_mapping), algo = "sha256")
}

test_that("P1 to P6 are defined in the proteomics core files, with their names kept", {
    if (is.null(utils::getSrcref(detectProteomicsFormat))) {
        testthat::skip("source references are kept under load_all only")
    }
    found <- vapply(names(protImportCoreDefinitions), protImportCoreSourceFile, character(1))
    expect_identical(found, protImportCoreDefinitions)
    # The GUI-side reader name stays in the module: the diagnostic traces and MCI-004 use it.
    expect_identical(protImportCoreSourceFile("importProtImportDataByFormat"),
        "mod_prot_import_processing_helpers.R")
})

test_that("back-reference #1: prepareDIAGUIReferences records the import metadata through the moved function", {
    manager <- newProtImportCoreManager("gui-references")
    withr::defer(closeDIAStateAdapter(manager))
    table <- protImportCoreTable()
    mapping <- list(run_col = "Run", quantity_col = "Precursor.Normalised")
    workflow_data <- list(data_tbl = table, data_cln = table, data_type = "peptide",
        column_mapping = mapping, config_list = list())

    references <- prepareDIAGUIReferences(manager, workflow_data)

    metadata <- getDIAImportMetadata(manager)
    expect_identical(metadata$reference_id, references$data_tbl$entity_id)
    expect_identical(metadata$column_mapping, mapping)
    expect_identical(metadata$data_type, "peptide")
    expect_identical(as.character(metadata$samples$Run), c("A", "B"))
    expect_identical(metadata$import_summary,
        list(n_rows = 3L, n_runs = 2L, n_proteins = NA, n_peptides = NA))
    expect_identical(manager$.context$design_adoption_pending, references$data_tbl$entity_id)
})

test_that("back-reference #2: transformDIADesignTable reads the samples through the moved metadata reader", {
    manager <- newProtImportCoreManager("design-samples")
    withr::defer(closeDIAStateAdapter(manager))
    table <- protImportCoreTable()
    source <- saveDIATableComponent(manager$.context, table)
    recordDIAImportMetadata(manager, table, list(data = table, data_type = "peptide",
        column_mapping = list(run_col = "Run", quantity_col = "Precursor.Normalised")), source)
    manager$.context$import_metadata_cache <- NULL

    binding <- testthat::with_mocked_bindings(
        transformDIADesignTable(manager, source, samples = c("B", "A"), strip_parser = TRUE),
        streamSavedWorkflowTable = function(...) stop("the sample catalogue was scanned from the assay"))

    expect_true(isWorkflowTableComponent(binding))
    expect_identical(manager$.context$import_metadata_cache$value$reference_id, source$entity_id)
    expect_identical(readSavedWorkflowTable(manager$.context$store, binding, "Run", limit = 3L)$data$Run,
        table$Run)
})

test_that("importProteomicsData reads by format with the options as a list; the GUI name delegates to it", {
    calls <- list()
    reader <- function(format) function(filepath, ...) {
        calls[[length(calls) + 1L]] <<- c(list(format = format, filepath = filepath), list(...))
        list(data = data.frame(x = 1), data_type = format, column_mapping = list())
    }
    readers <- stats::setNames(lapply(c("diann", "spectronaut", "fragpipe", "maxquant", "pd_tmt"), reader),
        c("diann", "spectronaut", "fragpipe", "maxquant", "pd_tmt"))

    for (format in names(readers)) importProteomicsData(format, "report.tsv", list(), readers = readers)
    expect_identical(calls, list(
        list(format = "diann", filepath = "report.tsv", use_precursor_norm = TRUE),
        list(format = "spectronaut", filepath = "report.tsv", quantity_type = "pg"),
        list(format = "fragpipe", filepath = "report.tsv", use_maxlfq = TRUE),
        list(format = "maxquant", filepath = "report.tsv", use_lfq = TRUE, filter_contaminants = TRUE),
        list(format = "pd_tmt", filepath = "report.tsv")))

    calls <- list()
    options <- list(use_precursor_norm = FALSE, quantity_type = "pep", use_maxlfq = FALSE,
        use_lfq = FALSE, filter_contaminants = FALSE)
    input <- list(diann_use_precursor_norm = FALSE, spectronaut_quantity = "pep", fragpipe_use_maxlfq = FALSE,
        maxquant_use_lfq = FALSE, maxquant_filter_contaminants = FALSE)
    for (format in names(readers)) {
        core <- importProteomicsData(format, "report.tsv", options, readers = readers)
        gui <- importProtImportDataByFormat(format, "report.tsv", input,
            importDiann = readers$diann, importSpectronaut = readers$spectronaut,
            importFragpipe = readers$fragpipe, importMaxquant = readers$maxquant,
            importPdTmt = readers$pd_tmt)
        expect_identical(gui, core)
    }
    expect_identical(calls[c(1L, 3L, 5L, 7L, 9L)], calls[c(2L, 4L, 6L, 8L, 10L)])
    expect_identical(calls[[7L]], list(format = "maxquant", filepath = "report.tsv",
        use_lfq = FALSE, filter_contaminants = FALSE))

    logged <- character()
    expect_error(importProteomicsData("unknown", "report.tsv", list(), readers = readers,
        logError = function(message) logged <<- c(logged, message)),
        "Failed to import data: Unsupported format: unknown", fixed = TRUE)
    expect_identical(logged, "Error in data import function: Unsupported format: unknown")
})

test_that("PD-TMT parity: text-numeric Excel abundances import as doubles and a text 0 is NA (F-26)", {
    text_path <- withr::local_tempfile(fileext = ".xlsx")
    writexl::write_xlsx(data.frame(Accession = c("P1", "P2"),
        `Abundance: F1: 126, Sample_A` = c("10", "0"),
        `Abundance: F1: 127N, Sample_B` = c(" 5.5", NA),
        check.names = FALSE), text_path)
    text <- suppressMessages(importProteomeDiscovererTMTData(text_path))
    expect_type(text$data$Abundance, "double")
    expect_identical(text$data$Run, c("126_Sample_A", "127N_Sample_B", "126_Sample_A", "127N_Sample_B"))
    expect_identical(text$data$Abundance, c(10, 5.5, NA, NA))

    # A mix of numeric and text cells across the abundance columns pivots once both are doubles.
    mixed_path <- withr::local_tempfile(fileext = ".xlsx")
    writexl::write_xlsx(data.frame(Accession = c("P1", "P2"),
        `Abundance: F1: 126, Sample_A` = c(10, 0),
        `Abundance: F1: 127N, Sample_B` = c("7", "0"),
        check.names = FALSE), mixed_path)
    mixed <- suppressMessages(importProteomeDiscovererTMTData(mixed_path))
    expect_identical(mixed$data$Abundance, c(10, 7, NA, NA))
})

test_that("PD-TMT leaves a non-numeric abundance column unconverted (P-05), so the pivot refuses the mix", {
    path <- withr::local_tempfile(fileext = ".xlsx")
    writexl::write_xlsx(data.frame(Accession = c("P1", "P2"),
        `Abundance: F1: 126, Sample_A` = c("10", "n.q."),
        `Abundance: F1: 127N, Sample_B` = c("7", "0"),
        check.names = FALSE), path)
    expect_error(suppressMessages(importProteomeDiscovererTMTData(path)), "Failed to pivot data")
})

test_that("Spectronaut parity: a character, numeric-like quantity column is converted to doubles", {
    report <- tibble::tibble(R.FileName = c("a", "b", "c"), PG.ProteinGroups = c("P1", "P1", "P2"),
        PG.Quantity = c("12.5", " 3", ""), PG.Qvalue = c(0.001, 0.002, 0.003),
        EG.ModifiedSequence = c("PEP", "PEP", "TIDE"), PEP.Quantity = c("4", NA, "0.5"),
        EG.Qvalue = c(0.01, 0.01, 0.01))
    read <- function(quantity_type) testthat::with_mocked_bindings(
        importSpectronautData("report.tsv", quantity_type = quantity_type),
        vroom = function(...) report, .package = "vroom")

    protein <- read("pg")
    expect_identical(protein$data$PG.Quantity, c(12.5, 3, NA))
    expect_identical(protein$data$PEP.Quantity, report$PEP.Quantity)
    expect_identical(protein$column_mapping$quantity_col, "PG.Quantity")

    peptide <- read("pep")
    expect_identical(peptide$data$PEP.Quantity, c(4, NA, 0.5))
    expect_identical(peptide$data$PG.Quantity, report$PG.Quantity)

    # Non-numeric text stays as the reader leaves it today (P-05).
    report$PG.Quantity <- c("12.5", "Filtered", "3")
    expect_identical(read("pg")$data$PG.Quantity, report$PG.Quantity)
})

test_that("DIA-NN, FragPipe and MaxQuant give identical output on their test inputs", {
    diann <- suppressMessages(importDIANNData(test_path("..", "testdata", "mock_import", "report.tsv")))
    fragpipe <- suppressMessages(importFragPipeData(
        test_path("..", "testdata", "mock_import", "combined_protein.tsv")))
    maxquant_path <- withr::local_tempfile(fileext = ".txt")
    utils::write.table(data.frame(Protein.IDs = c("P1", "P2", "P3"), Gene.names = c("G1", "G2", "G3"),
        LFQ.intensity.WT_1 = c(10, 0, 30), LFQ.intensity.KO_1 = c(40, 50, 60),
        Potential.contaminant = c("", "", "+"), Reverse = c("", "", ""), check.names = FALSE),
        maxquant_path, sep = "\t", row.names = FALSE, quote = FALSE)
    maxquant <- suppressMessages(suppressWarnings(importMaxQuantData(maxquant_path)))

    # Digests recorded at 60334e4b, before the KIT-102 W3 change.
    expect_identical(protImportReaderDigest(diann),
        "8db7d1a2173f4a39cdba330e860c85c6c304a61b408451bb98280127146c71c7")
    expect_identical(protImportReaderDigest(fragpipe),
        "adb18ab700016299e15334b5f01765359a9436ce112eca794a2b460dc463157c")
    expect_identical(protImportReaderDigest(maxquant),
        "289b06bb4ecefdd3af48bd3b293fe50abb60b1aad960bd68e0596498e217c966")
})
