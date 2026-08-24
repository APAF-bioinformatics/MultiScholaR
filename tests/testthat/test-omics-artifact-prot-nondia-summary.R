nondiaSummary027SkipDependencies <- function() {
    for (package in c(
        "arrow", "DBI", "duckdb", "filelock", "openxlsx", "rmarkdown"
    )) {
        testthat::skip_if_not_installed(package)
    }
    testthat::skip_if(Sys.which("pandoc") == "")
}

nondiaSummary027RepoPath <- function(...) {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    file.path(root, ...)
}

nondiaSummary027Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaSummary027RepoPath(
            "tests", "testdata", "omics-parity", "proteomics",
            "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaSummary027Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaSummary027WorkflowType <- function(format) {
    if (identical(format, "pd_tmt")) "TMT" else "LFQ"
}

nondiaSummary027Template <- function(format) {
    if (identical(format, "pd_tmt")) "TMT_report.rmd" else "LFQ_report.rmd"
}

nondiaSummary027Paths <- function(root, format) {
    paths <- list(
        base_dir = root,
        project_id = paste0("nondia-027-", format),
        omic_type = "proteomics",
        omic_label = paste0("summary-", format),
        results_dir = file.path(root, "results"),
        results_summary_dir = file.path(root, "results_summary"),
        publication_graphs_dir = file.path(root, "publication_graphs"),
        time_dir = file.path(root, "time"),
        qc_dir = file.path(root, "qc"),
        da_output_dir = file.path(root, "da_output"),
        pathway_dir = file.path(root, "pathway_enrichment"),
        source_dir = file.path(root, "source"),
        feature_qc_dir = file.path(root, "feature_qc"),
        integration_dir = file.path(root, "integration")
    )
    paths$qc_figures_dir <- file.path(paths$results_summary_dir, "QC_figures")
    paths$publication_figures_dir <- file.path(
        paths$results_summary_dir,
        "Publication_figures"
    )
    paths$publication_tables_dir <- file.path(
        paths$results_summary_dir,
        "Publication_tables"
    )
    directories <- unique(unlist(paths[grepl("_dir$", names(paths))]))
    invisible(lapply(
        directories,
        dir.create,
        recursive = TRUE,
        showWarnings = FALSE
    ))
    paths
}

nondiaSummary027Annotations <- function(protein_ids) {
    data.frame(
        Entry = protein_ids,
        gene_names = paste0("GENE", seq_along(protein_ids)),
        Organism = "Homo sapiens",
        stringsAsFactors = FALSE
    )
}

nondiaSummary027DaIndex <- function(generation_id, descriptor) {
    structure(
        list(
            schema = .PROT_DIA_DA_INDEX_SCHEMA,
            schema_version = .PROT_DIA_DA_INDEX_VERSION,
            backend = "artifact",
            run_id = paste0("da-", descriptor$identity$input_format),
            source_generation_id = generation_id,
            manifest_relative_path = "state/da/current.json",
            manifest_digest = strrep("a", 64L),
            contrasts = list(list(
                contrast_id = "contrast-ko-wt",
                contrast_name = "groupKO-groupWT",
                full_format = "KO_vs_WT=groupKO-groupWT",
                friendly_name = "KO_vs_WT",
                manifest_relative_path = "state/da/contrast.json",
                manifest_digest = strrep("b", 64L),
                long_table = NULL,
                query_specification = NULL,
                summary = list(
                    rows = 5L,
                    significant = 4L,
                    up = 2L,
                    down = 2L
                )
            ))
        ),
        class = c("MultiScholaRProtDiaDaIndex", "list")
    )
}

nondiaSummary027EnrichmentIndex <- function(descriptor) {
    structure(
        list(
            schema = .PROT_DIA_ENRICH_INDEX_SCHEMA,
            schema_version = .PROT_DIA_ENRICH_INDEX_VERSION,
            backend = "artifact",
            run_id = paste0("enrichment-", descriptor$identity$input_format),
            source = list(
                da_run_id = paste0("da-", descriptor$identity$input_format),
                contrast_id = "contrast-ko-wt",
                contrast_name = "groupKO-groupWT"
            ),
            parameters = list(
                backend = "gprofiler2",
                organism_taxid = "9606",
                q_cutoff = 0.05
            ),
            software = list(
                multischolar = "0.5.0",
                backend = "gprofiler2",
                backend_package = "gprofiler2",
                backend_version = "deterministic",
                r = as.character(getRversion())
            ),
            manifest_relative_path = "state/enrichment/current.json",
            manifest_digest = strrep("c", 64L),
            result_table = NULL,
            query_specification = NULL,
            products = list(list(
                name = "groupKO-groupWT_up_enrichment_results.tsv",
                relative_path = "results/enrichment.tsv",
                byte_digest = strrep("d", 64L),
                bytes = 128
            )),
            requests = list(list(
                request_id = "request-ko-wt-up",
                request_digest = strrep("e", 64L),
                backend = "gprofiler2",
                contrast = "groupKO-groupWT",
                direction = "up",
                status = "succeeded",
                execution_state = "replay",
                attempts = 0L,
                response = list(
                    request_id = "request-ko-wt-up",
                    response_digest = strrep("f", 64L),
                    rows = 2L,
                    manifest_relative_path = "state/enrichment/response.json",
                    manifest_digest = strrep("1", 64L)
                )
            )),
            interactive_plots = list(runtime_plot = new.env())
        ),
        class = c("MultiScholaRProtDiaEnrichIndex", "list")
    )
}

nondiaSummary027WritePng <- function(path, label) {
    grDevices::png(path, width = 320, height = 240)
    graphics::plot.new()
    graphics::text(0.5, 0.5, label)
    grDevices::dev.off()
}

nondiaSummary027WriteProducts <- function(built) {
    paths <- built$paths
    object <- built$finalObject
    nondiaSummary027WritePng(
        file.path(paths$time_dir, "12_correlation_filtered_combined_plots.png"),
        "correlation"
    )
    nondiaSummary027WritePng(
        file.path(paths$feature_qc_dir, "composite_QC_figure.png"),
        "QC"
    )
    writeLines(
        "%PDF-1.4 deterministic QC",
        file.path(paths$feature_qc_dir, "composite_QC_figure.pdf")
    )
    for (directory in c(
        "Interactive_Volcano_Plots", "NumSigDaMolecules", "Volcano_Plots",
        "Heatmap"
    )) {
        path <- file.path(paths$publication_graphs_dir, directory)
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        writeLines("deterministic figure", file.path(path, "figure.txt"))
    }
    writeLines(
        "contrast\tcount\nKO_vs_WT\t4",
        file.path(
            paths$publication_graphs_dir,
            "NumSigDaMolecules",
            "KO_vs_WT_num_sig_da_molecules.tab"
        )
    )
    normalized_prefix <- if (built$ruvRun) {
        "ruv_normalised_results_cln_with_replicates"
    } else {
        "normalised_results_cln_with_replicates"
    }
    readr::write_tsv(
        object@protein_quant_table,
        file.path(paths$feature_qc_dir, paste0(normalized_prefix, ".tsv"))
    )
    saveRDS(
        object,
        file.path(paths$feature_qc_dir, paste0(normalized_prefix, ".RDS"))
    )
    openxlsx::write.xlsx(
        data.frame(
            Protein.Ids = object@protein_quant_table[[object@protein_id_column]],
            log2FC = c(2, -2, 1, -1, 0),
            fdr_qvalue = c(0.01, 0.02, 0.03, 0.04, 0.5)
        ),
        file.path(paths$da_output_dir, "da_KO_vs_WT_long_annot.xlsx")
    )
    if (built$enrichmentRun) {
        readr::write_tsv(
            data.frame(
                term_id = c("GO:0006915", "GO:0007049"),
                p_value = c(0.001, 0.01),
                directionality = "positive"
            ),
            file.path(
                paths$pathway_dir,
                "groupKO-groupWT_up_enrichment_results.tsv"
            )
        )
    }
    invisible(paths)
}

nondiaSummary027WriteTemplate <- function(path) {
    writeLines(c(
        "---",
        "title: 'Artifact summary'",
        "output: html_document",
        "params:",
        "  omic_type: proteomics",
        "  experiment_label: default",
        "  workflow_name: default",
        "  timestamp: default",
        "  project_paths: !r NULL",
        "  report_dependencies: !r NULL",
        "---",
        "Descriptor: `r params$report_dependencies$descriptor_contract$descriptor_id`",
        "State: `r params$report_dependencies$final_state$state_name`",
        "Template: `r basename(knitr::current_input())`"
    ), path)
}

nondiaSummary027Build <- function(root, scenario) {
    format <- scenario$format
    workflow_type <- nondiaSummary027WorkflowType(format)
    descriptor <- protNonDiaReadthroughDescriptor(scenario$capability_id)
    paths <- nondiaSummary027Paths(root, format)
    project_dirs <- list(proteomics = paths)
    fixture <- nondiaSummary027RepoPath(scenario$fixture_path)
    source <- file.path(paths$source_dir, basename(fixture))
    stopifnot(file.copy(fixture, source))
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        paths$omic_label,
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "complete",
        differential_expression = "complete",
        differential_abundance = "complete",
        enrichment_analysis = "pending"
    )
    imported <- suppressMessages(nondiaSummary027Importer(format)(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType(workflow_type)
    stopifnot(isTRUE(persistProtNonDiaImportArtifacts(
        workflow,
        imported,
        source
    )$ok))
    runs <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$run_col]]
    ))
    groups <- ifelse(grepl("KO", runs, fixed = TRUE), "KO", "WT")
    workflow$design_matrix <- data.frame(
        Run = runs,
        group = groups,
        replicates = seq_along(runs),
        stringsAsFactors = FALSE
    )
    template <- nondiaSummary027Template(format)
    workflow$config_list <- list(globalParameters = list(
        workflow_type = workflow_type,
        report_template = if (identical(format, "fragpipe")) {
            "TMT_report.rmd"
        } else {
            template
        },
        use_limpa = FALSE
    ))
    workflow$contrasts_tbl <- data.frame(
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT",
        stringsAsFactors = FALSE
    )
    protein_ids <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$protein_col]]
    ))
    workflow$uniprot_dat_cln <- nondiaSummary027Annotations(protein_ids)
    workflow$aa_seq_tbl_final <- data.frame(
        accession = protein_ids,
        sequence = "PEPTIDE",
        stringsAsFactors = FALSE
    )
    suppressMessages(buildProtDesignStateCheckpoint(
        workflow,
        workflow_type,
        "non-DIA summary design fixture",
        validateColumnMapping = TRUE
    ))
    stopifnot(isTRUE(persistProtNonDiaDesignArtifacts(workflow)$ok))
    manager <- newProtNonDiaResumeStateManager(
        workflow$workflow_context,
        descriptor
    )
    workflow$state_manager <- manager
    object <- manager$getState()
    args <- object@args
    args$globalParameters <- workflow$config_list$globalParameters
    args$summary_provenance <- list(
        descriptor_id = descriptor$descriptor_id,
        source_format = format,
        final_stage = "correlation_filtered"
    )
    object@args <- args
    ruv_run <- !identical(format, "fragpipe")
    if (ruv_run) {
        manager$saveState(
            "ruv_corrected",
            object,
            args,
            "RUV summary fixture",
            audit_metadata = list(
                record_id = paste0("summary:ruv:", format),
                stage_id = "ruv"
            )
        )
    }
    manager$saveState(
        "correlation_filtered",
        object,
        args,
        "Final non-DIA summary handoff",
        audit_metadata = list(
            record_id = paste0("summary:correlation:", format),
            stage_id = "correlation_filter"
        )
    )
    workflow$config_list <- args
    workflow$organism_name <- "Homo sapiens"
    workflow$taxon_id <- 9606L
    workflow$da_ui_params <- list(selected_contrast = "KO_vs_WT")
    workflow$enrichment_ui_params <- list(selected_contrast = "KO_vs_WT")
    workflow$da_analysis_results_list <- nondiaSummary027DaIndex(
        manager$getCurrentGenerationId(),
        descriptor
    )
    enrichment_run <- !identical(format, "fragpipe")
    if (enrichment_run) {
        workflow$enrichment_artifact_index <-
            nondiaSummary027EnrichmentIndex(descriptor)
        workflow$tab_status$enrichment_analysis <- "complete"
    }
    readr::write_tsv(
        workflow$design_matrix,
        file.path(paths$source_dir, "design_matrix.tab")
    )
    readr::write_tsv(
        workflow$contrasts_tbl,
        file.path(paths$source_dir, "contrasts_tbl.tab")
    )
    writeLines(c(
        "Study Parameters",
        "Workflow Name: Non-DIA summary",
        "Timestamp: 2026-08-24 04:00:00"
    ), file.path(paths$source_dir, "study_parameters.txt"))
    template_dir <- file.path(paths$base_dir, "scripts", "proteomics")
    dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
    template_path <- file.path(template_dir, template)
    nondiaSummary027WriteTemplate(template_path)
    built <- list(
        projectDirs = project_dirs,
        paths = paths,
        workflowData = workflow,
        finalObject = object,
        descriptor = descriptor,
        scenario = scenario,
        template = template,
        templatePath = template_path,
        ruvRun = ruv_run,
        enrichmentRun = enrichment_run
    )
    nondiaSummary027WriteProducts(built)
    workflowStateReleaseHydrationCache(manager)
    built
}

nondiaSummary027PreserveGlobals <- function() {
    names <- c(
        "project_dirs", "config_list", "contrasts_tbl", "design_matrix",
        "organism_name", "taxon_id"
    )
    existing <- names[vapply(
        names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )]
    values <- mget(existing, envir = .GlobalEnv, inherits = FALSE)
    if (length(existing) > 0L) rm(list = existing, envir = .GlobalEnv)
    withr::defer({
        present <- intersect(names, ls(envir = .GlobalEnv))
        if (length(present) > 0L) rm(list = present, envir = .GlobalEnv)
        list2env(values, envir = .GlobalEnv)
    })
    names
}

test_that("all non-DIA summaries render exact reports and final S4 exports", {
    nondiaSummary027SkipDependencies()
    global_names <- nondiaSummary027PreserveGlobals()
    for (scenario in nondiaSummary027Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-027-summary-", scenario$format, "-")
        )
        built <- nondiaSummary027Build(root, scenario)
        workflow <- built$workflowData
        manager <- workflow$state_manager
        withr::defer(manager$close())
        workflow$config_list$globalParameters$workflow_type <-
            if (identical(built$template, "TMT_report.rmd")) "LFQ" else "TMT"
        workflow$config_list$globalParameters$report_template <-
            if (identical(built$template, "TMT_report.rmd")) {
                "LFQ_report.rmd"
            } else {
                "TMT_report.rmd"
            }
        template <- resolveProtSummaryReportTemplate(
            workflow,
            catFn = \(...) invisible(NULL)
        )
        expect_identical(template$templateFilename, built$template)
        expect_identical(template$dataStateUsed, "correlation_filtered")
        expect_identical(manager$getCacheInfo()$entries, 0L)

        report_dependencies <- prepareProtSummaryDependencies(
            workflow,
            built$projectDirs,
            kind = "report_reads",
            hydrate_final = FALSE
        )
        expect_s3_class(
            report_dependencies,
            "MultiScholaRProtDiaSummaryDependencies"
        )
        expect_null(report_dependencies$finalS4Object)
        expect_identical(manager$getCacheInfo()$entries, 0L)
        expect_identical(
            report_dependencies$manifest$descriptor_contract$descriptor_id,
            built$descriptor$descriptor_id
        )
        expect_identical(
            is.null(report_dependencies$manifest$analysis$enrichment),
            !built$enrichmentRun
        )
        report_dependencies <- prepareProtReportDependencies(
            report_dependencies,
            built$templatePath
        )
        expect_true(
            report_dependencies$manifest$dependencies$study_parameters$available
        )
        expect_false(
            report_dependencies$manifest$dependencies$qc_figures$required
        )
        report_path <- RenderReport(
            omic_type = "proteomics",
            experiment_label = paste0("summary-", scenario$format),
            rmd_filename = built$template,
            output_format = "html_document",
            project_dirs = built$projectDirs,
            report_dependencies = report_dependencies$manifest
        )
        expect_true(file.exists(report_path))
        expect_identical(
            basename(report_path),
            paste0(
                tools::file_path_sans_ext(built$template),
                "_proteomics_summary-",
                scenario$format,
                ".html"
            )
        )
        report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
        expect_true(grepl(
            built$descriptor$descriptor_id,
            report_text,
            fixed = TRUE
        ))
        recorded <- recordProtSummaryReportProduct(
            report_dependencies,
            report_path
        )
        expect_true(file.exists(recorded$metadataPath))
        expect_identical(recorded$manifest$report_product$byte_digest,
            artifactByteDigest(report_path)
        )

        final_dependencies <- prepareProtSummaryDependencies(
            workflow,
            built$projectDirs,
            kind = "final_export",
            hydrate_final = TRUE
        )
        expect_identical(final_dependencies$finalS4Object, built$finalObject)
        expect_identical(
            final_dependencies$designMatrix,
            built$finalObject@design_matrix
        )
        expect_identical(
            final_dependencies$contrastsTbl,
            workflow$contrasts_tbl
        )
        expect_identical(
            final_dependencies$manifest$scientific$design_digest,
            artifactSemanticDigest(workflow$design_matrix)
        )
        expect_identical(
            final_dependencies$manifest$scientific$protein_quant_digest,
            artifactSemanticDigest(built$finalObject@protein_quant_table)
        )
        expect_identical(
            final_dependencies$manifest$scientific$protein_id_digest,
            artifactSemanticDigest(built$finalObject@protein_id_table)
        )
        expect_identical(
            final_dependencies$manifest$scientific$args_digest,
            artifactSemanticDigest(built$finalObject@args)
        )
        final_path <- file.path(
            built$paths$integration_dir,
            paste0("proteomics_", scenario$format, "_final_s4.RDS")
        )
        final <- writeProtSummaryFinalExport(
            built$finalObject,
            final_path,
            final_dependencies
        )
        exported <- readRDS(final_path)
        expect_identical(exported, built$finalObject)
        expect_identical(
            exported@protein_quant_table,
            built$finalObject@protein_quant_table
        )
        expect_identical(
            exported@protein_id_table,
            built$finalObject@protein_id_table
        )
        expect_identical(
            exported@args$summary_provenance$descriptor_id,
            built$descriptor$descriptor_id
        )
        metadata <- workflowSessionDecodeMetadata(
            paste(readLines(final$metadataPath, warn = FALSE), collapse = "\n")
        )
        expect_identical(
            metadata$descriptor_contract$descriptor_digest,
            built$descriptor$descriptor_digest
        )
        expect_true(releaseProtSummaryDependencies(final_dependencies))
        expect_identical(manager$getCacheInfo()$entries, 0L)

        output <- new.env(parent = emptyenv())
        values <- new.env(parent = emptyenv())
        expect_true(completeProtSummaryWorkflowArgsSave(
            output = output,
            values = values,
            projectDirs = built$projectDirs,
            workflowData = workflow,
            experimentLabel = scenario$format,
            description = "descriptor-dispatched final export",
            renderTextFn = \(value) value,
            showNotificationFn = \(...) invisible(NULL),
            catFn = \(...) invisible(NULL)
        ))
        integrated_path <- file.path(
            built$paths$integration_dir,
            paste0("proteomics_", scenario$format, "_final_s4.RDS")
        )
        expect_identical(readRDS(integrated_path), built$finalObject)
        expect_true(file.exists(paste0(integrated_path, ".artifact.json")))
        expect_identical(manager$getCacheInfo()$entries, 0L)
        expect_false(any(vapply(
            global_names,
            exists,
            logical(1),
            envir = .GlobalEnv,
            inherits = FALSE
        )))
    }
})

test_that("non-DIA publication workbooks and session products preserve lanes", {
    nondiaSummary027SkipDependencies()
    global_names <- nondiaSummary027PreserveGlobals()
    for (scenario in nondiaSummary027Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-027-products-", scenario$format, "-")
        )
        built <- nondiaSummary027Build(root, scenario)
        workflow <- built$workflowData
        withr::defer(workflow$state_manager$close())
        copy_inputs <- prepareProtSummaryCopyInputs(
            workflow,
            built$projectDirs
        )
        output <- new.env(parent = emptyenv())
        values <- new.env(parent = emptyenv())
        values$files_copied <- FALSE
        runProtSummaryPublicationCopy(
            output = output,
            values = values,
            projectDirs = built$projectDirs,
            experimentLabel = scenario$format,
            description = "deterministic publication",
            contrastsTbl = copy_inputs$contrastsTbl,
            designMatrix = copy_inputs$designMatrix,
            artifactDependencies = copy_inputs$artifactDependencies,
            renderTextFn = \(value) value,
            showNotificationFn = \(...) invisible(NULL),
            catFn = \(...) invisible(NULL)
        )
        expect_true(values$files_copied)
        tables <- built$paths$publication_tables_dir
        da_workbook <- file.path(tables, "DA_results_proteomics.xlsx")
        enrich_workbook <- file.path(
            tables,
            "Pathway_enrichment_results_proteomics.xlsx"
        )
        expect_true(file.exists(da_workbook))
        expect_true(file.exists(enrich_workbook))
        expect_setequal(
            openxlsx::getSheetNames(da_workbook),
            c("DA_Results_Index", "DA_Sheet1")
        )
        enrichment_sheets <- openxlsx::getSheetNames(enrich_workbook)
        expect_true("Enrichment_Index" %in% enrichment_sheets)
        expect_identical("Enrich_Sheet1" %in% enrichment_sheets,
            built$enrichmentRun
        )
        da_index <- openxlsx::read.xlsx(da_workbook, "DA_Results_Index")
        expect_identical(da_index$Description, "KO_vs_WT")
        expect_true(file.exists(file.path(
            tables,
            if (built$ruvRun) {
                "RUV_normalised_results.tsv"
            } else {
                "normalised_results.tsv"
            }
        )))

        prepared <- prepareProtSummarySessionStateExport(
            projectDirs = built$projectDirs,
            experimentLabel = scenario$format,
            description = "deterministic session",
            workflowArgsSaved = TRUE,
            filesCopied = TRUE,
            reportGenerated = FALSE,
            reportPath = NULL,
            workflowData = workflow,
            exportDate = as.Date("2026-08-24"),
            timestamp = as.POSIXct("2026-08-24 05:00:00", tz = "UTC")
        )
        expect_identical(
            prepared$sessionState$report_template,
            built$template
        )
        audit <- prepared$sessionState$parameter_payload$artifact_summary_audit
        expect_identical(
            audit$descriptor_contract$descriptor_id,
            built$descriptor$descriptor_id
        )
        expect_identical(
            is.null(audit$analysis$enrichment),
            !built$enrichmentRun
        )
        expect_true(completeProtSummarySessionStateExport(
            projectDirs = built$projectDirs,
            experimentLabel = scenario$format,
            description = "deterministic session",
            workflowArgsSaved = TRUE,
            filesCopied = TRUE,
            reportGenerated = FALSE,
            reportPath = NULL,
            workflowData = workflow,
            prepareExportFn = \(...) prepared,
            showNotificationFn = \(...) invisible(NULL),
            logInfoFn = \(...) invisible(NULL)
        ))
        expect_identical(
            readRDS(prepared$sessionExportPath),
            prepared$sessionState
        )
        expect_false(any(vapply(
            global_names,
            exists,
            logical(1),
            envir = .GlobalEnv,
            inherits = FALSE
        )))
    }
})

test_that("summary switches and templates retain exact legacy rollback", {
    nondiaSummary027SkipDependencies()
    for (name in c("LFQ_report.rmd", "TMT_report.rmd")) {
        lines <- readLines(nondiaSummary027RepoPath(
            "inst", "reports", "proteomics", name
        ), warn = FALSE, n = 80L)
        expect_true(any(grepl("^[[:space:]]+project_paths:", lines)))
        expect_true(any(grepl("^[[:space:]]+report_dependencies:", lines)))
    }
    supported <- names(protNonDiaReadthroughDescriptors())
    expect_false(any(grepl("spectronaut", supported, ignore.case = TRUE)))

    scenario <- nondiaSummary027Scenarios()[[1L]]
    root <- withr::local_tempdir(pattern = "nondia-027-switch-")
    built <- nondiaSummary027Build(root, scenario)
    withr::defer(built$workflowData$state_manager$close())
    option <- paste0(
        "multischolar.",
        built$descriptor$identity$workflow_slug,
        ".summary_artifact_reads"
    )
    withr::with_options(stats::setNames(list("disabled"), option), {
        expect_null(prepareProtSummaryDependencies(
            built$workflowData,
            built$projectDirs,
            kind = "report_reads"
        ))
        legacy <- resolveProtSummaryFinalS4State(
            built$workflowData,
            catFn = \(...) invisible(NULL)
        )
        expect_identical(legacy$finalS4Object, built$finalObject)
    })
    export_option <- sub("artifact_reads$", "final_export", option)
    withr::with_options(stats::setNames(list("disabled"), export_option), {
        expect_null(prepareProtSummaryDependencies(
            built$workflowData,
            built$projectDirs,
            kind = "final_export"
        ))
    })
})

test_that("failed final export cannot publish product or sidecar", {
    nondiaSummary027SkipDependencies()
    scenario <- nondiaSummary027Scenarios()[[1L]]
    root <- withr::local_tempdir(pattern = "nondia-027-failure-")
    built <- nondiaSummary027Build(root, scenario)
    withr::defer(built$workflowData$state_manager$close())
    dependencies <- prepareProtSummaryDependencies(
        built$workflowData,
        built$projectDirs,
        kind = "final_export"
    )
    withr::defer(releaseProtSummaryDependencies(dependencies))
    path <- file.path(built$paths$integration_dir, "failed-final.RDS")
    expect_error(
        writeProtSummaryFinalExport(
            built$finalObject,
            path,
            dependencies,
            save_rds_fn = \(object, target) {
                saveRDS(object, target)
                stop("deterministic final export failure")
            }
        ),
        "deterministic final export failure",
        fixed = TRUE
    )
    expect_false(file.exists(path))
    expect_false(file.exists(paste0(path, ".artifact.json")))
})
