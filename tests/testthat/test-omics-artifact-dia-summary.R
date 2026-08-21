omics018SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics018Paths <- function(root) {
    paths <- list(
        base_dir = root,
        project_id = "omics-018-project",
        omic_type = "proteomics",
        omic_label = "dia-summary-study",
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
    invisible(lapply(
        unique(unlist(paths[grepl("_dir$", names(paths))])),
        dir.create,
        recursive = TRUE,
        showWarnings = FALSE
    ))
    list(proteomics = paths)
}

omics018Workflow <- function(project_dirs) {
    paths <- project_dirs$proteomics
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
        normalization = "disabled",
        differential_expression = "pending",
        differential_abundance = "pending",
        enrichment_analysis = "pending"
    )
    workflow
}

omics018Indexes <- function(generation_id) {
    da_index <- structure(
        list(
            schema = .PROT_DIA_DA_INDEX_SCHEMA,
            schema_version = .PROT_DIA_DA_INDEX_VERSION,
            backend = "artifact",
            run_id = "da-run-018",
            source_generation_id = generation_id,
            manifest_relative_path = "data/da/current.json",
            manifest_digest = paste(rep("a", 64L), collapse = ""),
            contrasts = list(list(
                contrast_id = "contrast-018",
                contrast_name = "groupKO-groupWT",
                full_format = "KO_vs_WT=groupKO-groupWT",
                friendly_name = "KO_vs_WT",
                manifest_relative_path = "data/da/contrast.json",
                manifest_digest = paste(rep("b", 64L), collapse = ""),
                long_table = NULL,
                query_specification = NULL,
                summary = list(
                    rows = 3L,
                    significant = 2L,
                    up = 1L,
                    down = 1L
                )
            ))
        ),
        class = c("MultiScholaRProtDiaDaIndex", "list")
    )
    enrichment_index <- structure(
        list(
            schema = .PROT_DIA_ENRICH_INDEX_SCHEMA,
            schema_version = .PROT_DIA_ENRICH_INDEX_VERSION,
            backend = "artifact",
            run_id = "enrich-run-018",
            source = list(
                contrast_name = "groupKO-groupWT",
                contrast_id = "contrast-018",
                da_run_id = "da-run-018"
            ),
            parameters = list(
                backend = "gprofiler2",
                organism_taxid = "9606",
                q_cutoff = 0.05
            ),
            manifest_relative_path = "data/enrichment/current.json",
            manifest_digest = paste(rep("c", 64L), collapse = ""),
            result_table = NULL,
            query_specification = NULL,
            products = list(list(
                name = "KO_vs_WT_up_enrichment_results.tsv",
                byte_digest = paste(rep("d", 64L), collapse = ""),
                bytes = 100
            )),
            requests = list(list(
                request_id = "request-018",
                backend = "gprofiler2",
                direction = "up",
                status = "succeeded",
                attempt = 1L,
                response = list(rows = 2L, digest = "response-018"),
                service = list(name = "g:Profiler", source = "live"),
                software = list(package = "gprofiler2", version = "0.2")
            )),
            interactive_plots = list(large_runtime_plot = new.env())
        ),
        class = c("MultiScholaRProtDiaEnrichIndex", "list")
    )
    list(da = da_index, enrichment = enrichment_index)
}

omics018Build <- function(root) {
    project_dirs <- omics018Paths(root)
    paths <- project_dirs$proteomics
    source <- file.path(paths$source_dir, "report.tsv")
    fixture <- testthat::test_path(
        "..", "testdata", "e2e", "prot_dia", "report.tsv"
    )
    stopifnot(file.copy(fixture, source))
    workflow <- omics018Workflow(project_dirs)
    imported <- suppressMessages(importDIANNData(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- "diann"
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType("DIA")
    stopifnot(isTRUE(persistProtDiaImportArtifacts(
        workflow,
        imported,
        source
    )$ok))

    runs <- unique(workflow$data_cln$Run)
    design <- data.frame(
        Run = runs,
        group = sub("_.*", "", runs),
        replicates = seq_along(runs),
        tech_rep_group = runs,
        stringsAsFactors = FALSE
    )
    contrasts <- data.frame(
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT",
        stringsAsFactors = FALSE
    )
    workflow$design_matrix <- design
    workflow$contrasts_tbl <- contrasts
    workflow$config_list <- list(globalParameters = list(
        workflow_type = "DIA",
        report_template = "DIANN_report.rmd",
        use_limpa = FALSE
    ))
    proteins <- unique(workflow$data_cln$Protein.Group)
    workflow$uniprot_dat_cln <- data.frame(
        Protein.Group = proteins,
        Gene = paste0("GENE", seq_along(proteins)),
        stringsAsFactors = FALSE
    )
    workflow$aa_seq_tbl_final <- data.frame(
        accession = proteins,
        sequence = rep("PEPTIDE", length(proteins)),
        stringsAsFactors = FALSE
    )
    peptide <- PeptideQuantitativeDataDiann(
        workflow$data_cln,
        design,
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$saveState(
        "raw_data_s4",
        peptide,
        workflow$config_list,
        "DIA-NN design memory checkpoint."
    )
    stopifnot(isTRUE(persistProtDiaDesignArtifacts(workflow)$ok))

    manager <- newProtDiaArtifactStateManager(workflow$workflow_context)
    workflow$state_manager <- manager
    protein_ids <- paste0("P", seq_len(4L))
    quantities <- matrix(
        seq_len(length(protein_ids) * length(runs)),
        nrow = length(protein_ids),
        ncol = length(runs),
        dimnames = list(protein_ids, runs)
    )
    protein_table <- data.frame(
        Protein.Group = protein_ids,
        as.data.frame(quantities, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    args <- list(
        globalParameters = workflow$config_list$globalParameters,
        deAnalysisParameters = list(contrasts = "groupKO-groupWT"),
        peptide_qc_audit = list(
            records = list(list(record_id = "peptide-audit-018")),
            current_record_id = "peptide-audit-018"
        ),
        protein_qc_audit = list(
            records = list(list(record_id = "protein-audit-018")),
            current_record_id = "protein-audit-018"
        ),
        normalization_artifact_audit = list(
            records = list(list(record_id = "normalization-audit-018")),
            current_record_id = "normalization-audit-018"
        )
    )
    protein <- ProteinQuantitativeData(
        protein_quant_table = protein_table,
        protein_id_column = "Protein.Group",
        protein_id_table = data.frame(
            Protein.Group = protein_ids,
            Gene = paste0("GENE", seq_along(protein_ids)),
            stringsAsFactors = FALSE
        ),
        design_matrix = design,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = args
    )
    manager$saveState(
        "correlation_filtered",
        protein,
        args,
        "Final DIA-NN summary handoff.",
        audit_metadata = list(
            record_id = "summary:correlation_filtered",
            stage_id = "correlation_filter"
        )
    )
    workflow$config_list <- args
    workflow$organism_name <- "Homo sapiens"
    workflow$taxon_id <- 9606L
    workflow$da_ui_params <- list(selected_contrast = "KO_vs_WT")
    workflow$enrichment_ui_params <- list(
        selected_contrast = "KO_vs_WT",
        method = "gprofiler2"
    )
    indexes <- omics018Indexes(manager$getCurrentGenerationId())
    workflow$da_analysis_results_list <- indexes$da
    workflow$enrichment_artifact_index <- indexes$enrichment
    workflow$tab_status$differential_expression <- "complete"
    workflow$tab_status$differential_abundance <- "complete"
    workflow$tab_status$enrichment_analysis <- "complete"

    readr::write_tsv(design, file.path(paths$source_dir, "design_matrix.tab"))
    readr::write_tsv(contrasts, file.path(paths$source_dir, "contrasts_tbl.tab"))
    writeLines("Study Parameters\n================", file.path(
        paths$source_dir,
        "study_parameters.txt"
    ))
    template_dir <- file.path(paths$base_dir, "scripts", "proteomics")
    dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(
        testthat::test_path(
            "..", "..", "inst", "reports", "proteomics", "DIANN_report.rmd"
        ),
        file.path(template_dir, "DIANN_report.rmd")
    )
    list(
        projectDirs = project_dirs,
        workflowData = workflow,
        finalObject = protein,
        templatePath = file.path(template_dir, "DIANN_report.rmd")
    )
}

test_that("OMICS-ART-018 reconstructs and exports the exact final DIA S4", {
    omics018SkipDependencies()
    fixture <- omics018Build(tempfile("omics-018-final-"))
    workflow <- fixture$workflowData
    manager <- workflow$state_manager
    on.exit(manager$close(), add = TRUE)

    global_names <- c(
        "project_dirs", "config_list", "contrasts_tbl", "organism_name",
        "taxon_id"
    )
    previous <- mget(
        global_names,
        envir = .GlobalEnv,
        ifnotfound = rep(list(NULL), length(global_names)),
        inherits = FALSE
    )
    existed <- vapply(global_names, exists, logical(1), envir = .GlobalEnv,
        inherits = FALSE)
    on.exit({
        present <- intersect(global_names, ls(envir = .GlobalEnv))
        if (length(present) > 0L) rm(list = present, envir = .GlobalEnv)
        for (name in global_names[existed]) {
            assign(name, previous[[name]], envir = .GlobalEnv)
        }
    }, add = TRUE)
    rm(list = global_names[existed], envir = .GlobalEnv)

    dependencies <- prepareProtDiaSummaryDependencies(
        workflow,
        fixture$projectDirs,
        kind = "final_export"
    )
    expect_s3_class(dependencies, "MultiScholaRProtDiaSummaryDependencies")
    expect_identical(dependencies$finalS4Object, fixture$finalObject)
    expect_identical(dependencies$designMatrix, fixture$finalObject@design_matrix)
    expect_identical(dependencies$contrastsTbl, workflow$contrasts_tbl)
    expect_identical(
        dependencies$manifest$analysis$da$contrasts[[1L]]$friendly_name,
        "KO_vs_WT"
    )
    expect_identical(
        dependencies$manifest$analysis$enrichment$requests[[1L]]$request_id,
        "request-018"
    )
    expect_identical(
        dependencies$manifest$scientific$peptide_audit_record_id,
        "peptide-audit-018"
    )
    expect_identical(
        dependencies$manifest$scientific$protein_audit_record_id,
        "protein-audit-018"
    )
    expect_false(is.environment(
        dependencies$manifest$analysis$enrichment$interactive_plots
    ))

    export_path <- file.path(
        fixture$projectDirs$proteomics$integration_dir,
        "proteomics_summary_final_s4.RDS"
    )
    result <- writeProtDiaSummaryFinalExport(
        fixture$finalObject,
        export_path,
        dependencies
    )
    exported <- readRDS(export_path)
    expect_identical(exported, fixture$finalObject)
    expect_identical(exported@protein_quant_table,
        fixture$finalObject@protein_quant_table)
    expect_identical(exported@protein_id_table,
        fixture$finalObject@protein_id_table)
    expect_identical(exported@args$peptide_qc_audit$current_record_id,
        "peptide-audit-018")
    expect_true(file.exists(result$metadataPath))
    metadata <- workflowSessionDecodeMetadata(
        paste(readLines(result$metadataPath), collapse = "\n")
    )
    expect_identical(metadata$analysis$da$run_id, "da-run-018")
    expect_identical(metadata$analysis$enrichment$run_id, "enrich-run-018")
    expect_identical(metadata$product$byte_digest, artifactByteDigest(export_path))

    expect_true(manager$getCacheInfo()$entries <= 1L)
    expect_true(releaseProtDiaSummaryDependencies(dependencies))
    expect_identical(manager$getCacheInfo()$entries, 0L)

    output <- new.env(parent = emptyenv())
    values <- new.env(parent = emptyenv())
    expect_true(completeProtSummaryWorkflowArgsSave(
        output = output,
        values = values,
        projectDirs = fixture$projectDirs,
        workflowData = workflow,
        experimentLabel = "summary-final",
        description = "artifact final export",
        renderTextFn = function(expr) expr,
        showNotificationFn = function(...) invisible(NULL),
        catFn = function(...) invisible(NULL)
    ))
    integrated_path <- file.path(
        fixture$projectDirs$proteomics$integration_dir,
        "proteomics_summary-final_final_s4.RDS"
    )
    expect_identical(readRDS(integrated_path), fixture$finalObject)
    expect_true(file.exists(paste0(integrated_path, ".artifact.json")))
    expect_identical(manager$getCacheInfo()$entries, 0L)
    expect_false(any(vapply(global_names, exists, logical(1), envir = .GlobalEnv,
        inherits = FALSE)))
})

test_that("OMICS-ART-018 declares report dependencies without global state", {
    omics018SkipDependencies()
    fixture <- omics018Build(tempfile("omics-018-report-"))
    workflow <- fixture$workflowData
    on.exit(workflow$state_manager$close(), add = TRUE)
    if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
        previous <- get("project_dirs", envir = .GlobalEnv)
        on.exit(assign("project_dirs", previous, envir = .GlobalEnv), add = TRUE)
        rm("project_dirs", envir = .GlobalEnv)
    }

    copy_inputs <- prepareProtSummaryCopyInputs(
        workflow,
        fixture$projectDirs
    )
    expect_identical(copy_inputs$designMatrix, workflow$design_matrix)
    expect_identical(copy_inputs$contrastsTbl, workflow$contrasts_tbl)
    expect_s3_class(
        copy_inputs$artifactDependencies,
        "MultiScholaRProtDiaSummaryDependencies"
    )

    recorded <- new.env(parent = emptyenv())
    copied <- runProtSummaryPublicationCopy(
        output = new.env(parent = emptyenv()),
        values = new.env(parent = emptyenv()),
        projectDirs = fixture$projectDirs,
        experimentLabel = "summary",
        description = "artifact report",
        contrastsTbl = copy_inputs$contrastsTbl,
        designMatrix = copy_inputs$designMatrix,
        artifactDependencies = copy_inputs$artifactDependencies,
        copyFn = function(omic_type, experiment_label, contrasts_tbl,
                          design_matrix, force, project_dirs) {
            recorded$args <- as.list(environment())
            invisible(list())
        },
        renderTextFn = function(expr) expr,
        showNotificationFn = function(...) invisible(NULL),
        catFn = function(...) invisible(NULL)
    )
    expect_identical(copied$project_dirs, fixture$projectDirs)
    expect_identical(recorded$args$design_matrix, workflow$design_matrix)
    expect_false(exists("project_dirs", envir = .GlobalEnv, inherits = FALSE))

    report_dependencies <- prepareProtDiaReportDependencies(
        copy_inputs$artifactDependencies,
        fixture$templatePath
    )
    expect_identical(
        report_dependencies$manifest$dependencies$study_parameters$available,
        TRUE
    )
    expect_identical(
        report_dependencies$manifest$dependencies$publication_tables$required,
        FALSE
    )
    report_path <- file.path(
        fixture$projectDirs$proteomics$results_summary_dir,
        "DIANN_report_proteomics_summary.html"
    )
    writeLines(
        "<html><body><h1>DIA summary</h1><p>KO_vs_WT</p></body></html>",
        report_path
    )
    report_result <- recordProtDiaSummaryReportProduct(
        report_dependencies,
        report_path
    )
    expect_true(file.exists(report_result$metadataPath))
    expect_identical(
        report_result$manifest$analysis$da$contrasts[[1L]]$friendly_name,
        "KO_vs_WT"
    )

    generated_path <- file.path(
        fixture$projectDirs$proteomics$results_summary_dir,
        "DIANN_report_proteomics_generated.html"
    )
    expect_true(runProtSummaryReportGeneration(
        output = new.env(parent = emptyenv()),
        values = new.env(parent = emptyenv()),
        experimentLabel = "generated",
        description = "artifact report",
        templateFilename = "DIANN_report.rmd",
        projectDirs = fixture$projectDirs,
        reportDependencies = report_dependencies,
        renderReportAvailableFn = function() TRUE,
        renderReportFn = function(omic_type, experiment_label, rmd_filename,
                                  project_dirs, report_dependencies) {
            recorded$render <- as.list(environment())
            writeLines(
                "<html><body><h1>DIA summary</h1></body></html>",
                generated_path
            )
            generated_path
        },
        activateReportFn = function(...) TRUE,
        showNotificationFn = function(...) invisible(NULL),
        logInfoFn = function(...) invisible(NULL)
    ))
    expect_identical(recorded$render$project_dirs, fixture$projectDirs)
    expect_identical(
        recorded$render$report_dependencies$manifest_digest,
        report_dependencies$manifest$manifest_digest
    )
    expect_true(file.exists(paste0(generated_path, ".artifact.json")))
    expect_false(exists("project_dirs", envir = .GlobalEnv, inherits = FALSE))

    unlink(
        fixture$projectDirs$proteomics$publication_figures_dir,
        recursive = TRUE
    )
    optional_missing <- prepareProtDiaReportDependencies(
        copy_inputs$artifactDependencies,
        fixture$templatePath
    )
    expect_false(
        optional_missing$manifest$dependencies$publication_figures$available
    )

    design_path <- file.path(
        fixture$projectDirs$proteomics$source_dir,
        "design_matrix.tab"
    )
    unlink(design_path)
    expect_error(
        prepareProtDiaReportDependencies(
            copy_inputs$artifactDependencies,
            fixture$templatePath
        ),
        class = "multischolar_missing_prot_dia_summary_dependency"
    )
})

test_that("OMICS-ART-018 switches and session audit remain independent", {
    omics018SkipDependencies()
    fixture <- omics018Build(tempfile("omics-018-switches-"))
    workflow <- fixture$workflowData
    on.exit(workflow$state_manager$close(), add = TRUE)

    withr::local_options(list(
        multischolar.prot_dia.summary_artifact_reads = "disabled",
        multischolar.prot_dia.summary_final_export = "enabled"
    ))
    expect_null(prepareProtDiaSummaryDependencies(
        workflow,
        fixture$projectDirs,
        kind = "report_reads"
    ))
    expect_null(protDiaSummarySessionAudit(workflow))
    final_dependencies <- prepareProtDiaSummaryDependencies(
        workflow,
        fixture$projectDirs,
        kind = "final_export"
    )
    expect_s3_class(
        final_dependencies,
        "MultiScholaRProtDiaSummaryDependencies"
    )
    releaseProtDiaSummaryDependencies(final_dependencies)

    withr::local_options(list(
        multischolar.prot_dia.summary_artifact_reads = "enabled",
        multischolar.prot_dia.summary_final_export = "disabled"
    ))
    expect_null(prepareProtDiaSummaryDependencies(
        workflow,
        fixture$projectDirs,
        kind = "final_export"
    ))
    report_dependencies <- prepareProtDiaSummaryDependencies(
        workflow,
        fixture$projectDirs,
        kind = "report_reads"
    )
    expect_s3_class(
        report_dependencies,
        "MultiScholaRProtDiaSummaryDependencies"
    )

    metadata <- summariseProtSummaryWorkflowMetadata(workflow)
    expect_identical(
        metadata$parameters$artifact_summary_audit$analysis$da$run_id,
        "da-run-018"
    )
    expect_identical(
        metadata$parameters$artifact_summary_audit$analysis$enrichment$run_id,
        "enrich-run-018"
    )
    expect_null(
        metadata$parameters$artifact_summary_audit$analysis$enrichment$
            interactive_plots
    )
})

test_that("OMICS-ART-018 explicit publication paths preserve product parity", {
    testthat::skip_if_not_installed("openxlsx")
    explicit_dirs <- module_ci_prot_summary_paths()
    legacy_dirs <- module_ci_prot_summary_paths()
    module_ci_prot_summary_write_publication_inputs(explicit_dirs)
    module_ci_prot_summary_write_publication_inputs(legacy_dirs)
    html_relative <- file.path(
        "Interactive_Volcano_Plots",
        "KO_vs_WT.html"
    )
    writeLines(
        "<html><body>KO_vs_WT</body></html>",
        file.path(
            explicit_dirs$proteomics$publication_graphs_dir,
            html_relative
        )
    )
    writeLines(
        "<html><body>KO_vs_WT</body></html>",
        file.path(
            legacy_dirs$proteomics$publication_graphs_dir,
            html_relative
        )
    )
    if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
        previous <- get("project_dirs", envir = .GlobalEnv)
        on.exit(assign("project_dirs", previous, envir = .GlobalEnv), add = TRUE)
        rm("project_dirs", envir = .GlobalEnv)
    }

    copyToResultsSummary(
        omic_type = "proteomics",
        experiment_label = "parity",
        contrasts_tbl = module_ci_prot_enrich_contrasts(),
        design_matrix = module_ci_prot_da_design(),
        force = TRUE,
        project_dirs = explicit_dirs
    )
    expect_false(exists("project_dirs", envir = .GlobalEnv, inherits = FALSE))
    assign("project_dirs", legacy_dirs, envir = .GlobalEnv)
    on.exit(rm("project_dirs", envir = .GlobalEnv), add = TRUE)
    copyToResultsSummary(
        omic_type = "proteomics",
        experiment_label = "parity",
        contrasts_tbl = module_ci_prot_enrich_contrasts(),
        design_matrix = module_ci_prot_da_design(),
        force = TRUE
    )

    explicit_root <- explicit_dirs$proteomics$results_summary_dir
    legacy_root <- legacy_dirs$proteomics$results_summary_dir
    explicit_files <- sort(list.files(explicit_root, recursive = TRUE))
    legacy_files <- sort(list.files(legacy_root, recursive = TRUE))
    expect_identical(explicit_files, legacy_files)

    workbook_names <- c(
        "DA_results_proteomics.xlsx",
        "Pathway_enrichment_results_proteomics.xlsx"
    )
    for (workbook in workbook_names) {
        explicit_path <- file.path(
            explicit_root,
            "Publication_tables",
            workbook
        )
        legacy_path <- file.path(
            legacy_root,
            "Publication_tables",
            workbook
        )
        sheets <- openxlsx::getSheetNames(explicit_path)
        expect_identical(sheets, openxlsx::getSheetNames(legacy_path))
        for (sheet in sheets) {
            expect_equal(
                openxlsx::read.xlsx(explicit_path, sheet = sheet),
                openxlsx::read.xlsx(legacy_path, sheet = sheet)
            )
        }
    }

    exact_products <- c(
        file.path("QC_figures", "proteomics_composite_QC_figure.pdf"),
        file.path("QC_figures", "proteomics_composite_QC_figure.png"),
        file.path("Publication_figures", html_relative),
        file.path("Publication_tables", "ruv_normalised_results.RDS"),
        file.path("Publication_tables", "RUV_normalised_results.tsv")
    )
    for (relative_path in exact_products) {
        expect_identical(
            readBin(file.path(explicit_root, relative_path), "raw", n = 1e7),
            readBin(file.path(legacy_root, relative_path), "raw", n = 1e7)
        )
    }
})

test_that("OMICS-ART-018 renders DIA templates from explicit parameters", {
    testthat::skip_if_not_installed("rmarkdown")
    project_dirs <- module_ci_prot_summary_paths()
    paths <- project_dirs$proteomics
    template_dir <- file.path(paths$base_dir, "scripts", "proteomics")
    dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(
        c(
            "Study Parameters",
            "Workflow Name: explicit-summary",
            "Timestamp: 2026-08-21 00:00:00"
        ),
        file.path(paths$source_dir, "study_parameters.txt")
    )
    template <- c(
        "---",
        "title: DIA explicit summary",
        "output: html_document",
        "params:",
        "  omic_type: proteomics",
        "  experiment_label: explicit",
        "  workflow_name: explicit-summary",
        "  timestamp: 2026-08-21 00:00:00",
        "  project_paths: !r NULL",
        "  report_dependencies: !r NULL",
        "---",
        "",
        "# `r params$workflow_name`",
        "",
        "Contrast: `r params$report_dependencies$analysis$contrast`"
    )
    for (name in c("DIANN_report.rmd", "DIANN_limpa_report.rmd")) {
        writeLines(template, file.path(template_dir, name))
    }
    if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
        previous <- get("project_dirs", envir = .GlobalEnv)
        on.exit(assign("project_dirs", previous, envir = .GlobalEnv), add = TRUE)
        rm("project_dirs", envir = .GlobalEnv)
    }

    for (name in c("DIANN_report.rmd", "DIANN_limpa_report.rmd")) {
        rendered <- suppressMessages(RenderReport(
            omic_type = "proteomics",
            experiment_label = "explicit",
            rmd_filename = name,
            output_format = "html_document",
            project_dirs = project_dirs,
            report_dependencies = list(
                analysis = list(contrast = "KO_vs_WT")
            )
        ))
        expect_true(file.exists(rendered))
        html <- paste(readLines(rendered, warn = FALSE), collapse = "\n")
        expect_match(html, "explicit-summary", fixed = TRUE)
        expect_match(html, "KO_vs_WT", fixed = TRUE)
        expect_false(exists("project_dirs", envir = .GlobalEnv,
            inherits = FALSE))
    }
})
