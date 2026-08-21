omics017SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics017TempRoot <- function() {
    root <- tempfile("omics-017-")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    root
}

omics017Build <- function(root) {
    paths <- list(
        base_dir = root,
        project_id = "omics-017-project",
        omic_type = "proteomics",
        omic_label = "dia-enrichment-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    context <- createWorkflowContext(
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
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$tab_status <- list(enrichment_analysis = "pending")
    source <- file.path(paths$source_dir, "report.tsv")
    fixture <- testthat::test_path(
        "..", "testdata", "e2e", "prot_dia", "report.tsv"
    )
    stopifnot(file.copy(fixture, source))
    imported <- suppressMessages(importDIANNData(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- "diann"
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType("DIA")
    stopifnot(isTRUE(persistProtDiaImportArtifacts(
        workflow, imported, source
    )$ok))

    runs <- unique(workflow$data_cln$Run)
    design <- data.frame(
        Run = runs,
        group = sub("_.*", "", runs),
        replicates = seq_along(runs),
        tech_rep_group = runs,
        stringsAsFactors = FALSE
    )
    workflow$design_matrix <- design
    workflow$config_list <- list(globalParameters = list(workflow_type = "DIA"))
    workflow$contrasts_tbl <- data.frame(
        contrasts = c("groupKO-groupWT", "groupWT-groupKO"),
        friendly_names = c("KO_vs_WT", "WT_vs_KO"),
        full_format = c(
            "KO_vs_WT=groupKO-groupWT",
            "WT_vs_KO=groupWT-groupKO"
        ),
        stringsAsFactors = FALSE
    )
    protein_ids <- paste0("P", seq_len(12L))
    workflow$uniprot_dat_cln <- data.frame(
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
    workflow$aa_seq_tbl_final <- data.frame(
        accession = protein_ids,
        sequence = "PEPTIDE",
        stringsAsFactors = FALSE
    )
    peptide <- PeptideQuantitativeDataDiann(
        workflow$data_cln,
        design,
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$saveState(
        "raw_data_s4", peptide, workflow$config_list, "DIA design checkpoint."
    )
    stopifnot(isTRUE(persistProtDiaDesignArtifacts(workflow)$ok))

    manager <- newProtDiaArtifactStateManager(context)
    workflow$state_manager <- manager
    values <- matrix(
        seq_len(length(protein_ids) * length(runs)),
        nrow = length(protein_ids),
        dimnames = list(protein_ids, runs)
    )
    table <- data.frame(
        Protein.Group = protein_ids,
        as.data.frame(values, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    args <- list(
        globalParameters = list(workflow_type = "DIA"),
        deAnalysisParameters = list(formula_string = "~ 0 + group")
    )
    protein <- ProteinQuantitativeData(
        protein_quant_table = table,
        protein_id_column = "Protein.Group",
        protein_id_table = data.frame(
            Protein.Group = protein_ids,
            gene_names = paste0("GENE", seq_along(protein_ids)),
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
        "DIA enrichment handoff.",
        audit_metadata = list(
            record_id = "enrichment:correlation_filtered",
            stage_id = "correlation_filter"
        )
    )
    list(workflow = workflow, paths = paths, protein = protein)
}

omics017DaResult <- function(object, comparison, full_format, reverse = FALSE) {
    ids <- object@protein_quant_table[[object@protein_id_column]]
    effect <- rep(c(2, -2, 1.5, -1.5, 0.2, -0.2), length.out = length(ids))
    if (isTRUE(reverse)) effect <- -effect
    q_values <- rep(c(0.01, 0.02, 0.03, 0.04, 0.2, 0.3),
                    length.out = length(ids))
    table <- data.frame(
        Protein.Group = ids,
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
            analysis_type = "RUV applied",
            stringsAsFactors = FALSE
        ))
    )
}

omics017PublishDa <- function(built) {
    contrasts <- built$workflow$contrasts_tbl
    results <- list(
        omics017DaResult(
            built$protein,
            contrasts$friendly_names[[1L]],
            contrasts$full_format[[1L]]
        ),
        omics017DaResult(
            built$protein,
            contrasts$friendly_names[[2L]],
            contrasts$full_format[[2L]],
            reverse = TRUE
        )
    )
    names(results) <- contrasts$contrasts
    prepared <- prepareProtDiaDaArtifactRun(
        built$workflow,
        contrasts,
        results,
        protDiaDaParameters("~ 0 + group", 0.05, 0, TRUE, TRUE),
        now = as.POSIXct("2026-08-21 07:00:00", tz = "UTC")
    )
    publishProtDiaDaArtifactRun(built$workflow, prepared)$index
}

omics017Gost <- function(...) {
    result <- buildProtEnrichDeterministicGprofilerResult("positive")
    result$meta$version <- "deterministic-2026-08-21"
    result
}

omics017Process <- function(built, setup, execution, q_cutoff = 0.05) {
    pathway <- file.path(built$paths$results_dir, "pathway_enrichment")
    dir.create(pathway, recursive = TRUE, showWarnings = FALSE)
    result <- suppressWarnings(suppressMessages(processEnrichments(
        setup$da_results,
        taxon_id = 9606,
        up_cutoff = 1,
        down_cutoff = 1,
        q_cutoff = q_cutoff,
        pathway_dir = pathway,
        go_annotations = setup$annotations,
        protein_id_column = "gene_name",
        contrast_names = names(setup$da_results@da_data),
        correction_method = "gSCS",
        execution_context = execution,
        gost_fn = omics017Gost,
        enricher_fn = runProtEnrichDeterministicEnricher
    )))
    list(result = result, pathway = pathway)
}

omics017ClusterProcess <- function(built, setup, execution) {
    pathway <- file.path(
        built$paths$results_dir,
        "pathway_enrichment_cluster"
    )
    dir.create(pathway, recursive = TRUE, showWarnings = FALSE)
    result <- suppressWarnings(suppressMessages(processEnrichments(
        setup$da_results,
        taxon_id = 12345,
        up_cutoff = 1,
        down_cutoff = 1,
        q_cutoff = 0.05,
        pathway_dir = pathway,
        go_annotations = setup$annotations,
        protein_id_column = "Protein.Group",
        contrast_names = names(setup$da_results@da_data),
        correction_method = "gSCS",
        execution_context = execution,
        gost_fn = omics017Gost,
        enricher_fn = runProtEnrichDeterministicEnricher
    )))
    list(result = result, pathway = pathway)
}

omics017Parameters <- function() {
    list(
        backend = "gprofiler2",
        organism_taxid = "9606",
        organism_supported = TRUE,
        organism_name = "Homo sapiens",
        up_cutoff = 1,
        down_cutoff = 1,
        q_cutoff = 0.05,
        correction_method = "gSCS",
        exclude_iea = FALSE,
        organism_filter_enabled = FALSE,
        organism_filter_applied = FALSE,
        organism_filter_stats = list(
            proteins_before = 0L,
            proteins_after = 0L,
            proteins_removed = 0L
        )
    )
}

omics017ResponseManifest <- function(built, record) {
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

omics017ResponseTables <- function(built, manifest) {
    store <- protDiaEnrichStore(built$workflow$workflow_context)
    list(
        identifiers = protDiaEnrichReadTable(store, manifest$identifiers),
        response = protDiaEnrichReadTable(store, manifest$response)
    )
}

omics017TemporaryQueryEntries <- function(built) {
    store <- protDiaEnrichStore(built$workflow$workflow_context)
    path <- artifactStoreResolveFile(store, store$relative_paths$duckdb_tmp)
    if (!dir.exists(path)) return(character())
    setdiff(list.files(path, all.files = TRUE), c(".", ".."))
}

test_that("DIA enrichment binds one exact DA contrast and ignores stale globals", {
    omics017SkipDependencies()
    built <- omics017Build(omics017TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    index <- omics017PublishDa(built)
    built$workflow$da_analysis_results_list <- index
    assign("contrasts_tbl", data.frame(contrasts = "WRONG"), .GlobalEnv)
    assign("design_matrix", data.frame(Run = "WRONG"), .GlobalEnv)
    assign("uniprot_dat_cln", data.frame(Entry = "WRONG"), .GlobalEnv)
    on.exit(rm(
        list = c("contrasts_tbl", "design_matrix", "uniprot_dat_cln"),
        envir = .GlobalEnv
    ), add = TRUE)

    setup <- protDiaEnrichExplicitSetup(built$workflow, "KO_vs_WT", built$protein)
    expect_identical(setup$entry$contrast_name, "groupKO-groupWT")
    expect_identical(names(setup$da_results@da_data), "groupKO-groupWT")
    expect_identical(setup$da_results@design_matrix, built$workflow$design_matrix)
    expect_false(any(setup$da_results@da_data[[1L]]$Protein.Group == "WRONG"))
    expect_identical(setup$annotations, built$workflow$uniprot_dat_cln)
})

test_that("gprofiler responses persist exact requests and support cache and replay", {
    omics017SkipDependencies()
    built <- omics017Build(omics017TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    index <- omics017PublishDa(built)
    built$workflow$da_analysis_results_list <- index
    setup <- protDiaEnrichExplicitSetup(built$workflow, "KO_vs_WT", built$protein)
    execution <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "auto",
        sleep_fn = function(...) NULL
    )
    processed <- omics017Process(built, setup, execution)
    records <- protDiaEnrichValidateExecutionRecords(
        execution, "gprofiler2", "groupKO-groupWT"
    )
    expect_setequal(vapply(records, `[[`, character(1), "execution_state"), "live")
    expect_setequal(vapply(records, `[[`, character(1), "direction"), c("up", "down"))

    response_manifests <- lapply(
        records,
        \(record) omics017ResponseManifest(built, record)
    )
    names(response_manifests) <- vapply(
        records,
        `[[`,
        character(1),
        "direction"
    )
    response_tables <- lapply(
        response_manifests,
        \(manifest) omics017ResponseTables(built, manifest)
    )
    expect_identical(
        response_tables$up$identifiers$identifier[
            response_tables$up$identifiers$role == "query"
        ],
        c("GENE1", "GENE3", "GENE7", "GENE9")
    )
    expect_identical(
        response_tables$down$identifiers$identifier[
            response_tables$down$identifiers$role == "query"
        ],
        c("GENE2", "GENE4", "GENE8", "GENE10")
    )
    expect_identical(
        response_tables$up$identifiers$identifier[
            response_tables$up$identifiers$role == "background"
        ],
        paste0("GENE", seq_len(12L))
    )
    expect_true(all(vapply(response_manifests, \(manifest) {
        identical(manifest$request$parameters$user_threshold, 0.05) &&
            identical(manifest$request$parameters$correction_method, "gSCS") &&
            identical(manifest$request$mapping$protein_id_column, "gene_name") &&
            identical(manifest$request$mapping$identifier_stats$sent_count, 4L) &&
            identical(manifest$request$mapping$background_stats$sent_count, 12L) &&
            identical(manifest$software$backend_package, "gprofiler2") &&
            identical(manifest$service_version, "deterministic-2026-08-21")
    }, logical(1))))

    all_results <- buildProtEnrichAllContrastResults(
        processed$result,
        list(method = "gprofiler2")
    )[["groupKO-groupWT"]]$gprofiler_results
    parameters <- omics017Parameters()
    prepared <- protDiaEnrichPersistRun(
        built$workflow,
        list(index = setup$index, entry = setup$entry),
        parameters,
        records,
        all_results,
        processed$pathway,
        now = as.POSIXct("2026-08-21 08:00:00", tz = "UTC")
    )
    published <- publishProtDiaEnrichRun(built$workflow, prepared)
    restored <- restoreProtDiaEnrichArtifactIndex(built$workflow)
    expect_true(published$published)
    expect_true(isProtDiaEnrichArtifactIndex(restored))
    expect_true(all(vapply(restored$interactive_plots, inherits,
                           logical(1), "plotly")))
    expect_setequal(
        vapply(published$manifest$products, `[[`, character(1), "name"),
        c(
            "groupKO-groupWT_up_enrichment_results.tsv",
            "groupKO-groupWT_up_enrichment_plot.html",
            "groupKO-groupWT_up_enrichment_plot.png",
            "groupKO-groupWT_down_enrichment_results.tsv",
            "groupKO-groupWT_down_enrichment_plot.html",
            "groupKO-groupWT_down_enrichment_plot.png"
        )
    )
    expect_true(all(vapply(published$manifest$products, \(product) {
        product$bytes > 0 && !startsWith(product$relative_path, "/")
    }, logical(1))))
    expect_identical(
        protDiaEnrichCompleteTable(built$workflow, restored),
        all_results
    )
    temporary_entries <- omics017TemporaryQueryEntries(built)
    page <- queryProtDiaEnrichPage(
        built$workflow,
        restored,
        projections = c("term_id", "directionality"),
        direction = "positive",
        limit = 2L
    )
    expect_identical(nrow(page$data), 2L)
    expect_true(page$has_more)
    expect_true(all(page$data$directionality == "positive"))
    expect_identical(omics017TemporaryQueryEntries(built), temporary_entries)
    expect_error(
        queryProtDiaEnrichPage(
            built$workflow,
            restored,
            limit = 2L,
            resource_policy = list(max_result_bytes = 1L)
        ),
        class = "multischolar_artifact_query_byte_limit_exceeded"
    )

    replay <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "replay",
        sleep_fn = function(...) NULL
    )
    replayed <- omics017Process(built, setup, replay)
    expect_s4_class(replayed$result, "EnrichmentResults")
    expect_true(all(vapply(
        replay$getRecords(),
        function(record) identical(record$execution_state, "replay"),
        logical(1)
    )))
    replay_results <- buildProtEnrichAllContrastResults(
        replayed$result,
        list(method = "gprofiler2")
    )[["groupKO-groupWT"]]$gprofiler_results
    replay_prepared <- protDiaEnrichPersistRun(
        built$workflow,
        list(index = setup$index, entry = setup$entry),
        parameters,
        replay$getRecords(),
        replay_results,
        replayed$pathway
    )
    expect_false(identical(
        replay_prepared$manifest$run_id,
        published$manifest$run_id
    ))
    expect_true(all(vapply(
        replay_prepared$manifest$requests,
        \(record) identical(record$execution_state, "replay"),
        logical(1)
    )))

    cached <- omics017Process(built, setup, replay)
    expect_s4_class(cached$result, "EnrichmentResults")
    states <- vapply(replay$getRecords(), `[[`, character(1), "execution_state")
    expect_true(all(tail(states, 2L) == "cache"))

    raw_manifest <- paste(readLines(published$current_path), collapse = "\n")
    expect_false(grepl(built$paths$base_dir, raw_manifest, fixed = TRUE))
    expect_false(grepl("password|token|authorization", raw_manifest, ignore.case = TRUE))
})

test_that("service and publication failures preserve current enrichment state", {
    omics017SkipDependencies()
    built <- omics017Build(omics017TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    index <- omics017PublishDa(built)
    built$workflow$da_analysis_results_list <- index
    setup <- protDiaEnrichExplicitSetup(built$workflow, "KO_vs_WT", built$protein)
    execution <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "auto",
        sleep_fn = function(...) NULL
    )
    processed <- omics017Process(built, setup, execution)
    records <- execution$getRecords()
    results <- buildProtEnrichAllContrastResults(
        processed$result,
        list(method = "gprofiler2")
    )[[1L]]$gprofiler_results
    parameters <- omics017Parameters()
    expect_error(
        protDiaEnrichPersistRun(
            built$workflow,
            list(index = setup$index, entry = setup$entry),
            parameters,
            records,
            results,
            processed$pathway,
            failure_injector = function(stage, value) {
                if (identical(
                    stage,
                    "before_enrichment_product_publication"
                )) {
                    stop("injected product failure")
                }
            }
        ),
        "injected product failure"
    )
    expect_false(file.exists(artifactResolveContainedPath(
        built$paths$base_dir,
        protDiaEnrichPaths(built$workflow$workflow_context)$current
    )))
    first <- protDiaEnrichPersistRun(
        built$workflow,
        list(index = setup$index, entry = setup$entry),
        parameters,
        records,
        results,
        processed$pathway
    )
    published <- publishProtDiaEnrichRun(built$workflow, first)
    before <- readLines(published$current_path, warn = FALSE)

    failed_context <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "live",
        sleep_fn = function(...) NULL
    )
    request <- protDiaEnrichCanonicalRequest(
        "gprofiler2",
        "groupKO-groupWT",
        "up",
        "GENE1",
        "GENE1",
        list(
            organism = "hsapiens",
            ordered_query = FALSE,
            sources = list("GO:BP"),
            user_threshold = 0.05,
            correction_method = "gSCS",
            exclude_iea = FALSE,
            evcodes = TRUE,
            domain_scope = "custom",
            significant = TRUE,
            highlight = TRUE
        )
    )
    expect_error(
        failed_context$execute(
            request,
            function() stop("deterministic outage"),
            max_retries = 1L
        ),
        "deterministic outage"
    )
    expect_error(
        failed_context$assertComplete(),
        class = "multischolar_prot_dia_enrichment_service_failed"
    )
    expect_identical(readLines(published$current_path, warn = FALSE), before)

    failed_persistence <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "live",
        sleep_fn = function(...) NULL,
        failure_injector = function(stage, value) {
            if (identical(stage, "before_enrichment_response_publication")) {
                stop("injected response failure")
            }
        }
    )
    persistence_request <- protDiaEnrichCanonicalRequest(
        "gprofiler2",
        "groupKO-groupWT",
        "up",
        "UNIQUE_FAILURE_ID",
        "UNIQUE_FAILURE_ID",
        request$parameters,
        mapping = protDiaEnrichGprofilerContext(
            "groupKO-groupWT",
            "up",
            "gene_name"
        )$mapping
    )
    expect_error(
        failed_persistence$execute(
            persistence_request,
            \() omics017Gost()
        ),
        "injected response failure"
    )
    expect_identical(
        failed_persistence$getRecords()[[1L]]$status,
        "failed"
    )
    expect_identical(readLines(published$current_path, warn = FALSE), before)

    forged <- first$manifest
    forged$source$contrast_manifest_digest <- strrep("0", 64L)
    forged$manifest_digest <- protDiaEnrichJsonDigest(forged)
    expect_error(
        protDiaEnrichValidateRunManifest(
            forged,
            built$workflow$workflow_context
        ),
        class = "multischolar_prot_dia_enrichment_source_mismatch"
    )
    forged <- first$manifest
    forged$products[[1L]]$bytes <- forged$products[[1L]]$bytes + 1
    forged$manifest_digest <- protDiaEnrichJsonDigest(forged)
    expect_error(
        protDiaEnrichValidateRunManifest(
            forged,
            built$workflow$workflow_context
        ),
        class = "multischolar_prot_dia_enrichment_product_mismatch"
    )
    response_manifest <- omics017ResponseManifest(built, records[[1L]])
    response_manifest$response$semantic_digest <- strrep("0", 64L)
    response_manifest$manifest_digest <- protDiaEnrichJsonDigest(
        response_manifest
    )
    expect_error(
        protDiaEnrichValidateResponseManifest(
            response_manifest,
            built$workflow$workflow_context
        ),
        class = "multischolar_prot_dia_enrichment_table_mismatch"
    )

    expect_error(
        publishProtDiaEnrichRun(
            built$workflow,
            first,
            failure_injector = function(stage, value) {
                if (identical(stage, "before_enrichment_current_publication")) {
                    stop("injected publication failure")
                }
            }
        ),
        "injected publication failure"
    )
    expect_identical(readLines(published$current_path, warn = FALSE), before)
    built$workflow$state_manager$saveState(
        "correlation_filtered",
        built$protein,
        built$protein@args,
        "DIA enrichment stale-generation test.",
        audit_metadata = list(
            record_id = "enrichment:stale-generation",
            stage_id = "correlation_filter"
        )
    )
    expect_error(
        publishProtDiaEnrichRun(built$workflow, first),
        class = "multischolar_stale_prot_dia_enrichment_da_index"
    )
    expect_error(
        protDiaEnrichCanonicalRequest(
            "gprofiler2", "contrast", "up", "/private/identifier", "P1", list()
        ),
        class = "multischolar_absolute_path_in_artifact_state"
    )
    expect_error(
        protDiaEnrichCanonicalRequest(
            "gprofiler2", "contrast", "up", "P1", "P1",
            list(token = "secret")
        ),
        class = "multischolar_unsafe_artifact_state_metadata"
    )
})

test_that("clusterProfiler provenance records mappings, duplicates, and empties", {
    omics017SkipDependencies()
    built <- omics017Build(omics017TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    index <- omics017PublishDa(built)
    built$workflow$da_analysis_results_list <- index
    setup <- protDiaEnrichExplicitSetup(
        built$workflow,
        "KO_vs_WT",
        built$protein
    )
    execution <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "auto",
        sleep_fn = function(...) NULL
    )
    term2gene <- data.frame(
        TERM = c("GO:1", "GO:1", "GO:2"),
        Entry = c("P1", "P2", "P3")
    )
    term2name <- data.frame(
        TERM = c("GO:1", "GO:2"),
        NAME = c("one", "two")
    )
    result <- protDiaEnrichExecuteClusterRequest(
        execution,
        genes = c("P1", "P1", "P9", NA_character_),
        background = c("P1", "P2", "P2"),
        contrast = "custom_contrast",
        direction = "up",
        taxon_id = "12345",
        q_cutoff = 0.05,
        exclude_iea = TRUE,
        protein_id_column = "Protein.Group",
        term2gene = term2gene,
        term2name = term2name,
        enricher_fn = runProtEnrichDeterministicEnricher
    )
    expect_s4_class(result, "ProtEnrichDeterministicResult")
    protDiaEnrichExecuteClusterRequest(
        execution,
        genes = character(),
        background = c("P1", "P2"),
        contrast = "custom_contrast",
        direction = "down",
        taxon_id = "12345",
        q_cutoff = 0.05,
        exclude_iea = TRUE,
        protein_id_column = "Protein.Group",
        term2gene = term2gene,
        term2name = term2name,
        enricher_fn = function(...) stop("must not run")
    )
    records <- protDiaEnrichValidateExecutionRecords(
        execution, "clusterprofiler", "custom_contrast"
    )
    expect_identical(records[[1L]]$execution_state, "local")
    expect_identical(records[[2L]]$status, "skipped_empty_input")
    expect_identical(records[[1L]]$response$rows, 3L)
    manifest <- omics017ResponseManifest(built, records[[1L]])
    tables <- omics017ResponseTables(built, manifest)
    query_ids <- tables$identifiers$identifier[
        tables$identifiers$role == "query"
    ]
    expect_identical(query_ids, c("P1", "P9"))
    expect_identical(manifest$request$mapping$identifier_stats$input_count, 4L)
    expect_identical(manifest$request$mapping$identifier_stats$missing_count, 1L)
    expect_identical(manifest$request$mapping$identifier_stats$duplicate_count, 1L)
    expect_identical(manifest$request$mapping$mapped_identifier_count, 1L)
    expect_identical(manifest$request$mapping$missing_mapping_count, 1L)
    expect_identical(manifest$request$mapping$background_stats$duplicate_count, 1L)
    expect_identical(manifest$request$parameters$correction_method, "BH")
    expect_identical(manifest$response_meta$local_execution$organism_taxid, "12345")
    expect_identical(manifest$software$backend_package, "clusterProfiler")

    zero_execution <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "auto",
        sleep_fn = function(...) NULL
    )
    zero_result <- buildProtEnrichDeterministicClusterProfilerResult("up")
    zero_result@result <- zero_result@result[0, , drop = FALSE]
    protDiaEnrichExecuteClusterRequest(
        zero_execution,
        genes = "P1",
        background = c("P1", "P2"),
        contrast = "zero_contrast",
        direction = "up",
        taxon_id = "12345",
        q_cutoff = 0.05,
        exclude_iea = TRUE,
        protein_id_column = "Protein.Group",
        term2gene = term2gene,
        term2name = term2name,
        enricher_fn = \(...) zero_result
    )
    protDiaEnrichExecuteClusterRequest(
        zero_execution,
        genes = character(),
        background = c("P1", "P2"),
        contrast = "zero_contrast",
        direction = "down",
        taxon_id = "12345",
        q_cutoff = 0.05,
        exclude_iea = TRUE,
        protein_id_column = "Protein.Group",
        term2gene = term2gene,
        term2name = term2name,
        enricher_fn = \(...) stop("must not run")
    )
    zero_records <- protDiaEnrichValidateExecutionRecords(
        zero_execution,
        "clusterprofiler",
        "zero_contrast"
    )
    expect_identical(zero_records[[1L]]$status, "succeeded")
    expect_identical(zero_records[[1L]]$response$rows, 0L)
    expect_identical(zero_records[[2L]]$status, "skipped_empty_input")
    expect_true(resolveEnrichmentOrganism(9606)$supported)
    expect_false(resolveEnrichmentOrganism(12345)$supported)

    route_execution <- newProtDiaEnrichExecutionContext(
        built$workflow,
        mode = "auto",
        sleep_fn = function(...) NULL
    )
    routed <- omics017ClusterProcess(built, setup, route_execution)
    route_records <- protDiaEnrichValidateExecutionRecords(
        route_execution,
        "clusterprofiler",
        "groupKO-groupWT"
    )
    route_results <- buildProtEnrichAllContrastResults(
        routed$result,
        list(method = "clusterprofiler")
    )[["groupKO-groupWT"]]$clusterprofiler_results
    route_parameters <- omics017Parameters()
    route_parameters$backend <- "clusterprofiler"
    route_parameters$organism_taxid <- "12345"
    route_parameters$organism_supported <- FALSE
    route_parameters$organism_name <- "Taxon ID 12345"
    prepared <- protDiaEnrichPersistRun(
        built$workflow,
        list(index = setup$index, entry = setup$entry),
        route_parameters,
        route_records,
        route_results,
        routed$pathway
    )
    published <- publishProtDiaEnrichRun(built$workflow, prepared)
    restored <- restoreProtDiaEnrichArtifactIndex(built$workflow)
    expect_true(published$published)
    expect_identical(
        protDiaEnrichCompleteTable(built$workflow, restored),
        route_results
    )
    expect_true(all(vapply(restored$interactive_plots, inherits,
                           logical(1), "plotly")))
    expect_identical(length(restored$products), 6L)
})

test_that("enrichment artifact flag is independent of DA state", {
    old <- options(multischolar.prot_dia.enrichment_persistence = "disabled")
    on.exit(options(old), add = TRUE)
    expect_identical(protDiaEnrichArtifactMode("persistence"), "disabled")
    expect_identical(protDiaDaArtifactMode("persistence"), "enabled")
})
