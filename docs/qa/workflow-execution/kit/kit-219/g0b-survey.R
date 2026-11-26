# The survey of g0b.R, as a test file so that the helpers and the test context of the W1 fixtures are available.
out_dir <- Sys.getenv("KIT219_G0B_OUT")

# The W1 setup functions of test-omic-summary-module.R, evaluated here without the rest of that file.
w1 <- parse(file.path("tests", "testthat", "test-omic-summary-module.R"))
wanted <- c("omicSummaryW1Roots", "omicSummaryW1FeatureSetup", "omicSummaryW1ProtSetup")
for (expr in w1) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) && as.character(expr[[2L]]) %in% wanted) eval(expr)
}

rows <- list()
survey <- function(case, omic_type, bus, paths) {
    manager <- bus$state_manager
    files <- if (dir.exists(paths$pathway_dir)) sort(list.files(paths$pathway_dir, recursive = TRUE)) else character()
    current <- workflowEnrichmentCurrentFiles(bus, omic_type)
    unmade <- tryCatch(workflowDaFilesUnmade(manager, paths, omic_type), error = function(e) paste("error:", conditionMessage(e)))
    omitted <- if (is.null(current)) character() else setdiff(files, current)
    runs <- length(unmade) == 0L
    rows[[length(rows) + 1L]] <<- data.frame(case = case, reference_state_owner = usesWorkflowReferenceStates(manager),
        checked = !is.null(current), pathway_files = length(files), kit212_guard_lets_copy_run = runs,
        omitted_if_copy_runs = if (runs) length(omitted) else 0L, omitted_names = paste(omitted, collapse = ","),
        listed_leaf = runs && length(omitted) > 0L, stringsAsFactors = FALSE)
}

test_that("G0 (iii): the W1 fixtures", {
    for (omic in c("lipid", "metab")) {
        for (inputs in c("complete", "optional_missing", "none")) {
            setup <- omicSummaryW1FeatureSetup(omic, inputs)
            survey(paste(omic, inputs), c(lipid = "lipidomics", metab = "metabolomics")[[omic]], setup$bus(), setup$paths)
        }
    }
    for (branch in c("DIA", "DIA_limpa", "TMT", "MaxQuant_LFQ")) {
        setup <- omicSummaryW1ProtSetup(branch)
        survey(paste("prot in-memory", branch), "proteomics", setup$bus(), setup$paths)
    }
    # The proteomics owner (the tab's server builds it; the same lines here).
    config <- module_ci_prot_summary_branch_config("DIA")
    bus <- module_ci_prot_summary_workflow(config)
    project_dirs <- module_ci_prot_summary_paths()
    paths <- project_dirs$proteomics
    withr::defer(unlink(paths$base_dir, recursive = TRUE))
    module_ci_prot_summary_write_publication_inputs(project_dirs)
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    object <- module_ci_prot_da_object(workflow_type = config$workflow_type, use_limpa = config$use_limpa,
        report_template = config$report_template, design = design, protein_quant_table = quant,
        args = module_ci_prot_summary_args(config))
    store_root <- file.path(paths$base_dir, "store")
    dir.create(store_root, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, "kit219-prot", config$workflow_type),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager$saveState("correlation_filtered", object, list(), "fixture state")
    bus$state_manager <- manager
    survey("prot owner DIA", "proteomics", bus, paths)
    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out_dir, "g0b-iii.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(table[, c("case", "reference_state_owner", "kit212_guard_lets_copy_run", "pathway_files", "omitted_if_copy_runs", "listed_leaf")])
    expect_true(nrow(table) == 11L)
})

test_that("G0 (iv): the assertions of the R2 files on the enrichment outputs", {
    pattern <- "Enrichment_Plots|Pathway_enrichment_results|Enrich_Sheet|Enrichment_Index|_enrichment_results"
    files <- sort(list.files(file.path("tests", "testthat"), pattern = "^test-.*[.]R$"))
    files <- setdiff(files, "test-kit219-enrichment-files-current.R")
    found <- list()
    for (file in files) {
        lines <- readLines(file.path("tests", "testthat", file), warn = FALSE)
        hit <- grep(pattern, lines)
        if (length(hit)) found[[file]] <- data.frame(file = file, line = hit, text = trimws(substr(lines[hit], 1L, 200L)),
            stringsAsFactors = FALSE)
    }
    table <- do.call(rbind, found)
    utils::write.table(table, file.path(out_dir, "g0b-iv.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("G0 (iv): ", nrow(table), " lines in ", length(found), " files\n")
    expect_true(nrow(table) > 0L)
})
