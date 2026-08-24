nondiaEnrich056SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaEnrich056RepoPath <- function(...) {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    file.path(root, ...)
}

nondiaEnrich056Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaEnrich056RepoPath(
            "tests", "testdata", "omics-parity", "proteomics",
            "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaEnrich056Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaEnrich056CapabilityId <- function(format) {
    lane <- switch(
        format,
        maxquant = "maxquant",
        fragpipe = "fragpipe",
        pd_tmt = "pd_tmt"
    )
    level <- if (identical(format, "pd_tmt")) "tmt" else "lfq"
    paste0("proteomics.", lane, ".protein.", level, ".v1")
}

nondiaEnrich056Paths <- function(root, format) {
    paths <- list(
        base_dir = root,
        project_id = paste0("nondia-056-", format),
        omic_type = "proteomics",
        omic_label = "enrichment-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

nondiaEnrich056Annotations <- function(protein_ids) {
    data.frame(
        Entry = protein_ids,
        gene_names = paste0("GENE", seq_along(protein_ids)),
        Organism = "Homo sapiens",
        go_id_go_biological_process = "GO:0006915; GO:0007049",
        go_term_go_biological_process = "apoptotic process; cell cycle",
        go_id_go_molecular_function = "GO:0005515",
        go_term_go_molecular_function = "protein binding",
        go_id_go_cellular_compartment = "GO:0005739",
        go_term_go_cellular_compartment = "mitochondrion",
        stringsAsFactors = FALSE
    )
}

nondiaEnrich056Build <- function(root, scenario) {
    format <- scenario$format
    descriptor <- protNonDiaReadthroughDescriptor(
        nondiaEnrich056CapabilityId(format)
    )
    paths <- nondiaEnrich056Paths(root, format)
    fixture <- nondiaEnrich056RepoPath(scenario$fixture_path)
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
        normalization = "complete",
        differential_expression = "complete",
        differential_abundance = "complete",
        enrichment_analysis = "pending"
    )
    imported <- suppressMessages(nondiaEnrich056Importer(format)(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow_type <- descriptor$identity$workflow_type
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
    workflow$config_list <- list(globalParameters = list(
        workflow_type = workflow_type
    ))
    workflow$contrasts_tbl <- data.frame(
        contrasts = c("groupKO-groupWT", "groupWT-groupKO"),
        friendly_names = c("KO_vs_WT", "WT_vs_KO"),
        full_format = c(
            "KO_vs_WT=groupKO-groupWT",
            "WT_vs_KO=groupWT-groupKO"
        ),
        stringsAsFactors = FALSE
    )
    imported_ids <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$protein_col]]
    ))
    workflow$uniprot_dat_cln <- nondiaEnrich056Annotations(imported_ids)
    workflow$aa_seq_tbl_final <- data.frame(
        accession = imported_ids,
        sequence = "PEPTIDE",
        stringsAsFactors = FALSE
    )
    suppressMessages(buildProtDesignStateCheckpoint(
        workflow,
        workflow_type,
        "non-DIA enrichment design fixture",
        validateColumnMapping = TRUE
    ))
    stopifnot(isTRUE(persistProtNonDiaDesignArtifacts(workflow)$ok))

    manager <- newProtNonDiaResumeStateManager(
        workflow$workflow_context,
        descriptor
    )
    workflow$state_manager <- manager
    protein_ids <- paste0("P", seq_len(12L))
    values <- matrix(
        seq_len(length(protein_ids) * length(runs)),
        nrow = length(protein_ids),
        dimnames = list(protein_ids, runs)
    )
    table <- data.frame(
        Protein.Ids = protein_ids,
        as.data.frame(values, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    args <- list(
        globalParameters = list(workflow_type = workflow_type),
        deAnalysisParameters = list(formula_string = "~ 0 + group")
    )
    protein <- ProteinQuantitativeData(
        protein_quant_table = table,
        protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(
            Protein.Ids = protein_ids,
            gene_name = paste0("GENE", seq_along(protein_ids)),
            stringsAsFactors = FALSE
        ),
        design_matrix = workflow$design_matrix,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "replicates",
        args = args
    )
    manager$saveState(
        "correlation_filtered",
        protein,
        args,
        "Final non-DIA enrichment handoff",
        audit_metadata = list(
            record_id = "enrichment:correlation_filtered",
            stage_id = "correlation_filter"
        )
    )
    workflow$config_list <- args
    list(
        workflow = workflow,
        paths = paths,
        protein = protein,
        descriptor = descriptor,
        scenario = scenario
    )
}

nondiaEnrich056DaResult <- function(
    object,
    comparison,
    full_format,
    reverse = FALSE
) {
    ids <- object@protein_quant_table[[object@protein_id_column]]
    effect <- rep(c(2, -2, 1.5, -1.5, 0.2, -0.2),
        length.out = length(ids)
    )
    if (isTRUE(reverse)) effect <- -effect
    q_values <- rep(c(0.01, 0.02, 0.03, 0.04, 0.2, 0.3),
        length.out = length(ids)
    )
    table <- data.frame(
        Protein.Ids = ids,
        gene_name = paste0("GENE", seq_along(ids)),
        comparison = comparison,
        numerator = sub("_vs_.*", "", comparison),
        denominator = sub(".*_vs_", "", comparison),
        log2FC = effect,
        raw_pvalue = q_values / 2,
        fdr_qvalue = q_values,
        fdr_value_bh_adjustment = q_values,
        stringsAsFactors = FALSE
    )
    model <- list(
        results = stats::setNames(list(table), full_format),
        fit.eb = list(coefficients = matrix(effect, ncol = 1L)),
        qvalue_warnings = list()
    )
    list(
        theObject = object,
        contrasts_results = model,
        contrasts_results_table = table,
        significant_rows = table[table$fdr_qvalue < 0.05, , drop = FALSE],
        norm_counts = object@protein_quant_table,
        da_proteins_wide = table,
        da_proteins_long = table,
        qvalue_warnings = list(),
        num_sig_da_molecules_first_go = list(table = data.frame(
            comparison = comparison,
            expression = full_format,
            status = c("Up", "Down"),
            counts = c(4L, 4L),
            analysis_type = "RUV skipped",
            stringsAsFactors = FALSE
        ))
    )
}

nondiaEnrich056PublishDa <- function(built) {
    contrasts <- built$workflow$contrasts_tbl
    results <- list(
        nondiaEnrich056DaResult(
            built$protein,
            contrasts$friendly_names[[1L]],
            contrasts$full_format[[1L]]
        ),
        nondiaEnrich056DaResult(
            built$protein,
            contrasts$friendly_names[[2L]],
            contrasts$full_format[[2L]],
            reverse = TRUE
        )
    )
    names(results) <- contrasts$contrasts
    prepared <- prepareProtDaArtifactRun(
        built$workflow,
        contrasts,
        results,
        protDaParameters("~ 0 + group", 0.05, 0, TRUE, TRUE),
        now = as.POSIXct("2026-08-24 02:00:00", tz = "UTC")
    )
    index <- publishProtDaArtifactRun(built$workflow, prepared)$index
    built$workflow$da_analysis_results_list <- index
    index
}

nondiaEnrich056GprofilerRequest <- function(direction) {
    identifiers <- if (identical(direction, "up")) {
        c("GENE1", "GENE1", NA, "", "GENE3", "GENE9")
    } else {
        c("GENE2", "GENE2", NA, "", "GENE4", "GENE10")
    }
    protDiaEnrichCanonicalRequest(
        backend = "gprofiler2",
        contrast = "groupKO-groupWT",
        direction = direction,
        identifiers = identifiers,
        background = c("GENE1", "GENE1", NA, "GENE2", "GENE12"),
        parameters = list(
            organism = "hsapiens",
            ordered_query = FALSE,
            sources = as.list(c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")),
            user_threshold = 0.05,
            correction_method = "gSCS",
            exclude_iea = FALSE,
            evcodes = TRUE,
            domain_scope = "custom",
            significant = TRUE,
            highlight = TRUE
        ),
        mapping = list(
            protein_id_column = "gene_name",
            transform = "remove_colon_suffix"
        )
    )
}

nondiaEnrich056ClusterRequest <- function(direction) {
    term2gene <- data.frame(
        TERM = c("GO:0006915", "GO:0007049", "GO:0005739"),
        Entry = c("P1", "P2", "P3"),
        stringsAsFactors = FALSE
    )
    term2name <- data.frame(
        TERM = term2gene$TERM,
        NAME = c("apoptotic process", "cell cycle", "mitochondrion"),
        stringsAsFactors = FALSE
    )
    genes <- if (identical(direction, "up")) {
        c("P1", "P1", NA, "", "P9")
    } else {
        c("P2", "P2", NA, "", "P10")
    }
    sent <- protDiaEnrichRequestIdentifiers(genes)
    mapped <- unique(as.character(term2gene$Entry))
    protDiaEnrichCanonicalRequest(
        backend = "clusterprofiler",
        contrast = "groupKO-groupWT",
        direction = direction,
        identifiers = genes,
        background = c("P1", "P1", NA, "P2", "P12"),
        parameters = list(
            organism_taxid = "12345",
            sources = as.list(c("GO:BP", "GO:CC", "GO:MF")),
            user_threshold = 0.05,
            correction_method = "BH",
            method = "enricher",
            exclude_iea = FALSE
        ),
        mapping = list(
            protein_id_column = "Protein.Ids",
            term2gene_digest = artifactSemanticDigest(term2gene),
            term2name_digest = artifactSemanticDigest(term2name),
            term2gene_rows = as.integer(nrow(term2gene)),
            term_count = as.integer(length(unique(term2gene$TERM))),
            mapped_identifier_count = as.integer(length(intersect(sent, mapped))),
            missing_mapping_count = as.integer(length(setdiff(sent, mapped)))
        )
    )
}

nondiaEnrich056RunRequests <- function(
    built,
    backend,
    mode = "live",
    execution = NULL
) {
    if (is.null(execution)) {
        execution <- newProtEnrichExecutionContext(
            built$workflow,
            mode = mode,
            sleep_fn = \(...) NULL
        )
    }
    for (direction in c("up", "down")) {
        request <- if (identical(backend, "gprofiler2")) {
            nondiaEnrich056GprofilerRequest(direction)
        } else {
            nondiaEnrich056ClusterRequest(direction)
        }
        response <- if (identical(backend, "gprofiler2")) {
            value <- buildProtEnrichDeterministicGprofilerResult(
                if (identical(direction, "up")) "positive" else "negative"
            )
            value$meta$version <- "deterministic-2026-08-24"
            value
        } else {
            buildProtEnrichDeterministicClusterProfilerResult(direction)
        }
        protDiaEnrichExecuteRequest(
            execution,
            request,
            call = \() response,
            max_retries = 1L,
            wait_time = 0
        )
    }
    execution
}

nondiaEnrich056RunParameters <- function(backend) {
    supported <- identical(backend, "gprofiler2")
    method_info <- list(
        method = backend,
        supported = supported,
        species_name = if (supported) "Homo sapiens" else "Taxon ID 12345"
    )
    input <- list(
        organism_taxid = if (supported) "9606" else "12345",
        up_cutoff = 1,
        down_cutoff = 1,
        q_cutoff = 0.05,
        correction_method = "gSCS",
        enable_organism_filter = FALSE
    )
    protDiaEnrichParameters(input, method_info)
}

nondiaEnrich056WritePlotProducts <- function(pathway, contrast) {
    for (direction in c("up", "down")) {
        prefix <- file.path(
            pathway,
            paste0(contrast, "_", direction, "_enrichment_plot")
        )
        writeLines("<html>deterministic enrichment plot</html>", paste0(prefix, ".html"))
        grDevices::png(paste0(prefix, ".png"), width = 320, height = 240)
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste("deterministic", direction))
        grDevices::dev.off()
    }
}

nondiaEnrich056PrepareRun <- function(built, backend, execution) {
    setup <- protEnrichExplicitSetup(
        built$workflow,
        "KO_vs_WT",
        built$protein
    )
    pathway <- file.path(
        built$paths$results_dir,
        paste0("pathway_", backend)
    )
    result <- runProtEnrichDeterministicProcessEnrichments(
        setup$da_results,
        taxon_id = if (identical(backend, "gprofiler2")) 9606 else 12345,
        up_cutoff = 1,
        down_cutoff = 1,
        q_cutoff = 0.05,
        pathway_dir = pathway,
        go_annotations = setup$annotations,
        protein_id_column = if (identical(backend, "gprofiler2")) {
            "gene_name"
        } else {
            "Protein.Ids"
        },
        contrast_names = setup$entry$contrast_name
    )
    nondiaEnrich056WritePlotProducts(
        pathway,
        setup$entry$contrast_name
    )
    all_results <- buildProtEnrichAllContrastResults(
        result,
        list(method = backend)
    )
    selected <- all_results[[setup$entry$contrast_name]]
    result_table <- if (identical(backend, "gprofiler2")) {
        selected$gprofiler_results
    } else {
        selected$clusterprofiler_results
    }
    prepared <- protDiaEnrichPersistRun(
        built$workflow,
        list(index = setup$index, entry = setup$entry),
        nondiaEnrich056RunParameters(backend),
        execution$getRecords(),
        result_table,
        pathway,
        now = as.POSIXct("2026-08-24 03:00:00", tz = "UTC")
    )
    list(
        setup = setup,
        pathway = pathway,
        result = result,
        result_table = result_table,
        prepared = prepared
    )
}

nondiaEnrich056PublishRun <- function(built, backend, execution) {
    run <- nondiaEnrich056PrepareRun(built, backend, execution)
    published <- publishProtDiaEnrichRun(built$workflow, run$prepared)
    index <- protEnrichRunIndex(
        published$manifest,
        built$workflow$workflow_context
    )
    index <- protDiaEnrichIndexWithPlots(index, run$result)
    built$workflow$enrichment_artifact_index <- index
    built$workflow$enrichment_analysis_results <- protEnrichReactivePayload(
        built$workflow,
        index
    )
    c(run, list(published = published, index = index))
}

nondiaEnrich056ResponseManifest <- function(built, record) {
    path <- artifactResolveContainedPath(
        built$paths$base_dir,
        record$response$manifest_relative_path,
        must_exist = TRUE
    )
    protDiaEnrichReadJson(
        path,
        \(value) protDiaEnrichValidateResponseManifest(
            value,
            built$workflow$workflow_context
        )
    )
}

nondiaEnrich056MovedWorkflow <- function(built, moved) {
    entries <- list.files(
        built$paths$base_dir,
        all.files = TRUE,
        no.. = TRUE,
        full.names = TRUE
    )
    stopifnot(all(file.copy(entries, moved, recursive = TRUE)))
    paths <- nondiaEnrich056Paths(moved, built$scenario$format)
    prepared <- createProtNonDiaResumeContext(
        paths,
        "enrichment-study",
        built$descriptor$descriptor_id
    )
    stopifnot(isTRUE(prepared$enabled))
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- prepared$context
    workflow$state_manager <- newProtNonDiaResumeStateManager(
        prepared$context,
        built$descriptor
    )
    workflow$artifact_stage_results <- NULL
    workflow$design_matrix <- built$workflow$design_matrix
    workflow$contrasts_tbl <- built$workflow$contrasts_tbl
    workflow$uniprot_dat_cln <- built$workflow$uniprot_dat_cln
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "disabled",
        differential_expression = "pending",
        differential_abundance = "pending",
        enrichment_analysis = "pending"
    )
    workflow$da_analysis_results_list <- restoreProtDaArtifactIndex(workflow)
    list(workflow = workflow, paths = paths)
}

test_that("all non-DIA lanes persist exact enrichment provenance and outputs", {
    nondiaEnrich056SkipDependencies()
    previous_index <- NULL
    gprofiler_result_digests <- character()
    cluster_result_digests <- character()
    gprofiler_product_digests <- list()
    cluster_product_digests <- list()
    for (scenario in nondiaEnrich056Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-056-enrichment-", scenario$format, "-")
        )
        built <- nondiaEnrich056Build(root, scenario)
        withr::defer(built$workflow$state_manager$close())
        da_index <- nondiaEnrich056PublishDa(built)
        setup <- protEnrichExplicitSetup(
            built$workflow,
            "KO_vs_WT",
            built$protein
        )
        expect_identical(setup$entry$contrast_name, "groupKO-groupWT")
        expect_identical(names(setup$da_results@da_data), "groupKO-groupWT")
        expect_identical(setup$index$run_id, da_index$run_id)
        da_manifest_path <- artifactResolveContainedPath(
            root,
            da_index$manifest_relative_path,
            must_exist = TRUE
        )
        da_manifest <- protDiaDaReadJson(
            da_manifest_path,
            \(value) protDiaDaValidateRunManifest(
                value,
                built$workflow$workflow_context
            )
        )
        expect_identical(
            da_manifest$descriptor_contract$descriptor_id,
            built$descriptor$descriptor_id
        )
        expect_identical(
            da_manifest$descriptor_contract$descriptor_digest,
            built$descriptor$descriptor_digest
        )
        if (!is.null(previous_index)) {
            built$workflow$da_analysis_results_list <- previous_index
            expect_error(
                protEnrichExplicitSetup(
                    built$workflow,
                    "KO_vs_WT",
                    built$protein
                ),
                class = "multischolar_stale_prot_dia_enrichment_da_index"
            )
            built$workflow$da_analysis_results_list <- da_index
        }

        gprofiler_execution <- nondiaEnrich056RunRequests(
            built,
            "gprofiler2"
        )
        gprofiler <- nondiaEnrich056PublishRun(
            built,
            "gprofiler2",
            gprofiler_execution
        )
        expect_identical(
            gprofiler$published$manifest$source$da_run_id,
            da_index$run_id
        )
        expect_identical(
            gprofiler$published$manifest$source$contrast_id,
            setup$entry$contrast_id
        )
        expect_true(grepl(
            built$descriptor$identity$workflow_slug,
            gprofiler$index$query_specification$query_id,
            fixed = TRUE
        ))
        expect_identical(gprofiler$published$manifest$parameters$backend,
            "gprofiler2"
        )
        expect_identical(gprofiler$published$manifest$parameters$organism_taxid,
            "9606"
        )
        expect_identical(length(gprofiler$index$products), 6L)
        gprofiler_result_digests <- c(
            gprofiler_result_digests,
            gprofiler$index$result_table$semantic_digest
        )
        gprofiler_product_digests[[scenario$format]] <- sort(vapply(
            gprofiler$index$products,
            `[[`,
            character(1),
            "byte_digest"
        ))
        expect_identical(
            protEnrichCompleteTable(built$workflow, gprofiler$index),
            gprofiler$result_table
        )
        page <- queryProtEnrichPage(
            built$workflow,
            gprofiler$index,
            direction = "positive",
            limit = 2L
        )
        expect_identical(nrow(page$data), 2L)
        expect_true(page$has_more)
        expect_false("data" %in% names(gprofiler$index$result_table))
        expect_lte(nrow(built$workflow$enrichment_analysis_results$gprofiler_results),
            .PROT_DIA_ENRICH_MAX_RESULT_ROWS
        )
        for (record in gprofiler_execution$getRecords()) {
            response <- nondiaEnrich056ResponseManifest(built, record)
            expect_identical(response$service_version,
                "deterministic-2026-08-24"
            )
            expect_identical(
                response$request$parameters$sources,
                as.list(c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"))
            )
            expect_identical(response$request$mapping$identifier_stats$missing_count,
                1L
            )
            expect_identical(response$request$mapping$identifier_stats$blank_count,
                1L
            )
            expect_identical(response$request$mapping$identifier_stats$duplicate_count,
                1L
            )
        }

        cluster_execution <- nondiaEnrich056RunRequests(
            built,
            "clusterprofiler"
        )
        cluster <- nondiaEnrich056PublishRun(
            built,
            "clusterprofiler",
            cluster_execution
        )
        expect_identical(cluster$published$manifest$parameters$backend,
            "clusterprofiler"
        )
        expect_identical(cluster$published$manifest$parameters$organism_taxid,
            "12345"
        )
        expect_identical(length(cluster$index$products), 6L)
        cluster_result_digests <- c(
            cluster_result_digests,
            cluster$index$result_table$semantic_digest
        )
        cluster_product_digests[[scenario$format]] <- sort(vapply(
            cluster$index$products,
            `[[`,
            character(1),
            "byte_digest"
        ))
        expect_identical(
            protEnrichCompleteTable(built$workflow, cluster$index),
            cluster$result_table
        )
        expect_true(all(vapply(
            cluster_execution$getRecords(),
            \(record) record$execution_state == "local",
            logical(1)
        )))
        cluster_manifest <- nondiaEnrich056ResponseManifest(
            built,
            cluster_execution$getRecords()[[1L]]
        )
        expect_identical(
            cluster_manifest$request$mapping$mapped_identifier_count,
            1L
        )
        expect_identical(
            cluster_manifest$request$mapping$missing_mapping_count,
            1L
        )
        expect_identical(
            cluster_manifest$request$parameters$sources,
            as.list(c("GO:BP", "GO:CC", "GO:MF"))
        )
        expect_identical(
            cluster_manifest$software$backend_package,
            "clusterProfiler"
        )

        current_text <- paste(
            readLines(cluster$published$current_path, warn = FALSE),
            collapse = "\n"
        )
        expect_false(grepl(root, current_text, fixed = TRUE))
        expect_false(grepl("credential|password|secret|token|api_key",
            current_text,
            ignore.case = TRUE
        ))
        tuple_option <- paste0(
            "multischolar.",
            built$descriptor$identity$workflow_slug,
            ".enrichment_persistence"
        )
        withr::with_options(stats::setNames(list("disabled"), tuple_option), {
            expect_false(protEnrichArtifactEligible(built$workflow))
            expect_true(protDaArtifactEligible(built$workflow))
        })
        query_option <- sub("persistence$", "queries", tuple_option)
        withr::with_options(stats::setNames(list("disabled"), query_option), {
            expect_false(protEnrichArtifactEligible(built$workflow, "queries"))
            expect_error(
                queryProtEnrichPage(
                    built$workflow,
                    cluster$index,
                    limit = 1L
                ),
                class = "multischolar_prot_dia_enrichment_queries_disabled"
            )
        })
        previous_index <- da_index
    }
    expect_identical(length(unique(gprofiler_result_digests)), 1L)
    expect_identical(length(unique(cluster_result_digests)), 1L)
    expect_true(all(vapply(
        gprofiler_product_digests[-1L],
        identical,
        logical(1),
        gprofiler_product_digests[[1L]]
    )))
    expect_true(all(vapply(
        cluster_product_digests[-1L],
        identical,
        logical(1),
        cluster_product_digests[[1L]]
    )))
})

test_that("non-DIA gprofiler cache replay and moved restore are portable", {
    nondiaEnrich056SkipDependencies()
    for (scenario in nondiaEnrich056Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-056-replay-", scenario$format, "-")
        )
        built <- nondiaEnrich056Build(root, scenario)
        withr::defer(built$workflow$state_manager$close())
        nondiaEnrich056PublishDa(built)
        live <- nondiaEnrich056RunRequests(built, "gprofiler2")
        published <- nondiaEnrich056PublishRun(built, "gprofiler2", live)

        replay <- nondiaEnrich056RunRequests(
            built,
            "gprofiler2",
            mode = "replay"
        )
        expect_true(all(vapply(
            replay$getRecords(),
            \(record) record$execution_state == "replay",
            logical(1)
        )))
        nondiaEnrich056RunRequests(
            built,
            "gprofiler2",
            mode = "replay",
            execution = replay
        )
        states <- vapply(
            replay$getRecords(),
            `[[`,
            character(1),
            "execution_state"
        )
        expect_identical(tail(states, 2L), c("cache", "cache"))

        moved <- withr::local_tempdir(
            pattern = paste0("nondia-056-moved-", scenario$format, "-")
        )
        resumed <- nondiaEnrich056MovedWorkflow(built, moved)
        withr::defer(resumed$workflow$state_manager$close())
        restored <- restoreProtEnrichArtifactIndex(resumed$workflow)
        expect_true(isProtEnrichArtifactIndex(restored))
        expect_identical(restored$run_id, published$index$run_id)
        expect_identical(
            protEnrichCompleteTable(resumed$workflow, restored),
            published$result_table
        )
        expect_identical(
            resumed$workflow$workflow_context$getProjectRoot(),
            normalizePath(moved, winslash = "/")
        )
    }
})

test_that("non-DIA enrichment empty and failure paths cannot corrupt current", {
    nondiaEnrich056SkipDependencies()
    scenario <- nondiaEnrich056Scenarios()[[1L]]
    root <- withr::local_tempdir(pattern = "nondia-056-failure-")
    built <- nondiaEnrich056Build(root, scenario)
    withr::defer(built$workflow$state_manager$close())
    nondiaEnrich056PublishDa(built)
    live <- nondiaEnrich056RunRequests(built, "gprofiler2")
    baseline <- nondiaEnrich056PublishRun(built, "gprofiler2", live)
    current_before <- readLines(baseline$published$current_path, warn = FALSE)

    failed <- newProtEnrichExecutionContext(
        built$workflow,
        mode = "live",
        sleep_fn = \(...) NULL
    )
    expect_error(
        protDiaEnrichExecuteRequest(
            failed,
            nondiaEnrich056GprofilerRequest("up"),
            call = \() stop("deterministic service failure"),
            max_retries = 1L,
            wait_time = 0
        ),
        "deterministic service failure",
        fixed = TRUE
    )
    expect_identical(failed$getRecords()[[1L]]$status, "failed")
    expect_identical(
        readLines(baseline$published$current_path, warn = FALSE),
        current_before
    )

    cluster_execution <- nondiaEnrich056RunRequests(
        built,
        "clusterprofiler"
    )
    cluster <- nondiaEnrich056PrepareRun(
        built,
        "clusterprofiler",
        cluster_execution
    )
    expect_error(
        publishProtDiaEnrichRun(
            built$workflow,
            cluster$prepared,
            failure_injector = \(stage, ...) {
                if (identical(stage, "before_enrichment_current_publication")) {
                    stop("deterministic publication failure")
                }
            }
        ),
        "deterministic publication failure",
        fixed = TRUE
    )
    expect_identical(
        readLines(baseline$published$current_path, warn = FALSE),
        current_before
    )

    empty_execution <- newProtEnrichExecutionContext(
        built$workflow,
        mode = "live"
    )
    for (direction in c("up", "down")) {
        request <- nondiaEnrich056GprofilerRequest(direction)
        request$identifiers <- character()
        request$mapping$identifier_stats <- protDiaEnrichIdentifierStats(
            character(),
            character()
        )
        empty_execution$recordSkip(request)
    }
    setup <- protEnrichExplicitSetup(
        built$workflow,
        "KO_vs_WT",
        built$protein
    )
    empty <- protDiaEnrichPersistRun(
        built$workflow,
        list(index = setup$index, entry = setup$entry),
        nondiaEnrich056RunParameters("gprofiler2"),
        empty_execution$getRecords(),
        data.frame(
            directionality = character(),
            analysis_method = character(),
            stringsAsFactors = FALSE
        ),
        withr::local_tempdir(pattern = "nondia-056-empty-products-")
    )
    empty <- publishProtDiaEnrichRun(built$workflow, empty)
    expect_identical(length(empty$manifest$products), 0L)
    expect_true(all(vapply(
        empty$manifest$requests,
        \(record) record$status == "skipped_empty_input",
        logical(1)
    )))
})

test_that("non-DIA enrichment rejects unsafe metadata and renders legacy state", {
    expect_error(
        protDiaEnrichCanonicalRequest(
            "gprofiler2",
            "contrast",
            "up",
            "P1",
            "P1",
            list(),
            mapping = list(api_token = "secret")
        ),
        class = "multischolar_unsafe_artifact_state_metadata"
    )
    expect_error(
        protDiaEnrichCanonicalRequest(
            "gprofiler2",
            "contrast",
            "up",
            "/private/source/P1",
            "P1",
            list(),
            mapping = list()
        ),
        class = "multischolar_absolute_path_in_artifact_state"
    )

    contrasts <- tibble::tibble(
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT"
    )
    legacy <- createEnrichmentResults(contrasts)
    legacy@enrichment_plotly <- list(
        "groupKO-groupWT" = list(
            up = buildProtEnrichDeterministicPlot("gprofiler2", "up"),
            down = buildProtEnrichDeterministicPlot("gprofiler2", "down")
        )
    )
    expect_identical(
        protEnrichInteractivePlots(legacy, "groupKO-groupWT"),
        legacy@enrichment_plotly[["groupKO-groupWT"]]
    )
})
