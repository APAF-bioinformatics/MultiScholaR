nondiaArtifact026SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaArtifact026RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

nondiaArtifact026Manifest <- function() {
    jsonlite::read_json(
        nondiaArtifact026RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
}

nondiaArtifact026Scenarios <- function() {
    nondiaArtifact026Manifest()$fixture_scenarios
}

nondiaArtifact026Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaArtifact026WorkflowType <- function(format) {
    if (identical(format, "pd_tmt")) "TMT" else "LFQ"
}

nondiaArtifact026CapabilityId <- function(format) {
    paste0(
        "proteomics.",
        switch(format, maxquant = "maxquant", fragpipe = "fragpipe", pd_tmt = "pd_tmt"),
        ".protein.",
        if (identical(format, "pd_tmt")) "tmt" else "lfq",
        ".v1"
    )
}

nondiaArtifact026VendorObject <- function(scenario) {
    imported <- suppressMessages(
        nondiaArtifact026Importer(scenario$format)(
            nondiaArtifact026RepoPath(scenario$fixture_path)
        )
    )
    data <- as.data.frame(imported$data)
    mapping <- imported$column_mapping
    runs <- unique(as.character(data[[mapping$run_col]]))
    groups <- ifelse(grepl("KO", runs, fixed = TRUE), "KO", "WT")
    replicates <- ave(seq_along(runs), groups, FUN = seq_along)
    workflow_type <- nondiaArtifact026WorkflowType(scenario$format)
    manager <- new.env(parent = emptyenv())
    manager$saveState <- \(state_name, s4_data_object, ...) {
        manager$object <- s4_data_object
        invisible(state_name)
    }
    workflow <- list2env(list(
        data_cln = data,
        column_mapping = mapping,
        design_matrix = data.frame(
            Run = runs,
            group = groups,
            replicates = paste0("R", replicates),
            stringsAsFactors = FALSE
        ),
        config_list = list(
            globalParameters = list(workflow_type = workflow_type)
        ),
        state_manager = manager
    ), parent = emptyenv())
    suppressMessages(buildProtDesignStateCheckpoint(
        workflow,
        workflow_type,
        "non-DIA DA scientific fixture",
        validateColumnMapping = TRUE
    ))
    manager$object
}

nondiaArtifact026RunDa <- function(object, treat_lfc_cutoff = 0) {
    matrix <- as.matrix(object@protein_quant_table[
        ,
        setdiff(names(object@protein_quant_table), object@protein_id_column),
        drop = FALSE
    ])
    rownames(matrix) <- object@protein_quant_table[[object@protein_id_column]]
    design <- object@design_matrix
    rownames(design) <- design$Run
    result <- suppressMessages(runTestsContrasts(
        data = matrix,
        contrast_strings = "KO_vs_WT=groupKO-groupWT",
        design_matrix = design,
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = treat_lfc_cutoff,
        eBayes_trend = TRUE,
        eBayes_robust = TRUE
    ))$results[[1L]]
    result$feature_id <- rownames(result)
    result[order(result$feature_id, method = "radix"), , drop = FALSE]
}

nondiaArtifact026Paths <- function(root, project_id) {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "proteomics",
        omic_label = "da-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

nondiaArtifact026Build <- function(root, scenario, protein_count = 14L) {
    format <- scenario$format
    capability_id <- nondiaArtifact026CapabilityId(format)
    descriptor <- protNonDiaReadthroughDescriptor(capability_id)
    paths <- nondiaArtifact026Paths(root, paste0("nondia-026-", format))
    fixture <- nondiaArtifact026RepoPath(scenario$fixture_path)
    source <- file.path(paths$source_dir, basename(fixture))
    stopifnot(file.copy(fixture, source))
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "da-study",
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
        differential_expression = "pending",
        differential_abundance = "pending",
        enrichment_analysis = "locked"
    )
    imported <- suppressMessages(nondiaArtifact026Importer(format)(source))
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
    proteins <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$protein_col]]
    ))
    workflow$uniprot_dat_cln <- data.frame(
        Protein.Ids = proteins,
        gene_name = paste0("GENE 'quoted' ", seq_along(proteins)),
        stringsAsFactors = FALSE
    )
    workflow$aa_seq_tbl_final <- data.frame(
        accession = proteins,
        sequence = "PEPTIDE",
        stringsAsFactors = FALSE
    )
    suppressMessages(buildProtDesignStateCheckpoint(
        workflow,
        workflow_type,
        "non-DIA DA design fixture",
        validateColumnMapping = TRUE
    ))
    stopifnot(isTRUE(persistProtNonDiaDesignArtifacts(workflow)$ok))
    manager <- newProtNonDiaResumeStateManager(
        workflow$workflow_context,
        descriptor
    )
    workflow$state_manager <- manager
    protein_ids <- c(
        paste0("P", seq_len(protein_count - 1L)),
        "P%_';--"
    )
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
            gene_name = paste0("GENE 'quoted' ", seq_along(protein_ids)),
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
        "Final non-DIA DA handoff",
        audit_metadata = list(
            record_id = "da:correlation_filtered",
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

nondiaArtifact026Result <- function(
    object,
    comparison,
    full_format,
    reverse = FALSE
) {
    ids <- object@protein_quant_table[[object@protein_id_column]]
    count <- length(ids)
    effect <- rep(c(2, -2, 0.5, -0.5), length.out = count)
    if (isTRUE(reverse)) effect <- -effect
    q_values <- rep(c(0.01, 0.01, 0.04, 0.2), length.out = count)
    q_values[[5L]] <- NA_real_
    table <- data.frame(
        Protein.Ids = ids,
        gene_name = paste0("GENE 'quoted' ", seq_along(ids)),
        comparison = comparison,
        numerator = sub("_vs_.*", "", comparison),
        denominator = sub(".*_vs_", "", comparison),
        log2FC = effect,
        raw_pvalue = pmin(1, q_values / 2),
        fdr_qvalue = q_values,
        fdr_value_bh_adjustment = q_values,
        stringsAsFactors = FALSE
    )
    model <- list(
        results = stats::setNames(list(table), full_format),
        fit.eb = list(
            coefficients = matrix(effect, ncol = 1L),
            trend = TRUE,
            robust = TRUE,
            treat = FALSE
        ),
        qvalue_warnings = list(
            fallback = "BH used when qvalue is unavailable"
        )
    )
    list(
        theObject = object,
        contrasts_results = model,
        contrasts_results_table = table,
        significant_rows = table[0L, , drop = FALSE],
        norm_counts = object@protein_quant_table,
        da_proteins_wide = table,
        da_proteins_long = table,
        qvalue_warnings = model$qvalue_warnings,
        num_sig_da_molecules_first_go = list(table = data.frame(
            comparison = comparison,
            expression = full_format,
            status = c("Up", "Down"),
            counts = c(2L, 2L),
            analysis_type = "RUV skipped",
            stringsAsFactors = FALSE
        ))
    )
}

nondiaArtifact026Results <- function(built) {
    contrasts <- built$workflow$contrasts_tbl
    results <- list(
        nondiaArtifact026Result(
            built$protein,
            contrasts$friendly_names[[1L]],
            contrasts$full_format[[1L]]
        ),
        nondiaArtifact026Result(
            built$protein,
            contrasts$friendly_names[[2L]],
            contrasts$full_format[[2L]],
            reverse = TRUE
        )
    )
    names(results) <- contrasts$contrasts
    results
}

test_that("actual non-DIA DA statistics remain at each frozen lane oracle", {
    manifest <- nondiaArtifact026Manifest()
    oracle <- jsonlite::read_json(
        nondiaArtifact026RepoPath(
            "tests", "testdata", "omics-parity", "proteomics",
            "memory-oracle.json"
        ),
        simplifyVector = FALSE
    )
    expected <- oracle$scientific_oracles[[1L]]$da_rows
    expected <- stats::setNames(
        expected,
        vapply(expected, `[[`, character(1), "feature_id")
    )
    tolerance <- manifest$comparison_contract$da_numeric_tolerance$absolute
    for (scenario in nondiaArtifact026Scenarios()) {
        object <- nondiaArtifact026VendorObject(scenario)
        normalized <- suppressMessages(normaliseBetweenSamples(
            object,
            normalisation_method = "cyclicloess"
        ))
        actual <- nondiaArtifact026RunDa(normalized)
        treated <- nondiaArtifact026RunDa(
            normalized,
            treat_lfc_cutoff = 1
        )
        expect_identical(treated$feature_id, actual$feature_id)
        expect_equal(treated$logFC, actual$logFC, tolerance = tolerance)
        expect_true(all(is.finite(treated$raw_pvalue)))
        for (feature_id in names(expected)) {
            row <- actual[actual$feature_id == feature_id, , drop = FALSE]
            expected_row <- expected[[feature_id]]
            for (field in c(
                "logFC", "AveExpr", "t", "raw_pvalue", "adj.P.Val",
                "fdr_value_bh_adjustment"
            )) {
                expect_equal(
                    row[[field]],
                    expected_row[[field]],
                    tolerance = tolerance,
                    info = paste(scenario$format, feature_id, field)
                )
            }
            expect_true(is.na(row$fdr_qvalue))
        }
    }
})

test_that("non-DIA DA persistence and bounded queries are tuple exact", {
    nondiaArtifact026SkipDependencies()
    for (scenario in nondiaArtifact026Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-026-da-", scenario$format, "-")
        )
        built <- nondiaArtifact026Build(root, scenario)
        withr::defer(built$workflow$state_manager$close())
        results <- nondiaArtifact026Results(built)
        parameters <- protDaParameters("~ 0 + group", 0.05, 0, TRUE, TRUE)
        prepared <- prepareProtDaArtifactRun(
            built$workflow,
            built$workflow$contrasts_tbl,
            results,
            parameters,
            now = as.POSIXct("2026-08-24 00:00:00", tz = "UTC")
        )
        expect_true(prepared$enabled)
        expect_null(restoreProtDaArtifactIndex(built$workflow))
        published <- publishProtDaArtifactRun(built$workflow, prepared)
        expect_true(published$published)
        index <- restoreProtDaArtifactIndex(built$workflow)
        expect_true(isProtDiaDaArtifactIndex(index))
        expect_identical(
            published$manifest$descriptor_contract$descriptor_id,
            built$descriptor$descriptor_id
        )
        expect_identical(
            published$manifest$source$generation_id,
            built$workflow$state_manager$getCurrentGenerationId()
        )
        expect_setequal(
            names(published$manifest$source$stage_refs),
            c("import", "design")
        )
        expect_identical(published$manifest$parameters, parameters)
        expect_false(grepl(
            root,
            paste(readLines(published$current_path, warn = FALSE), collapse = "\n"),
            fixed = TRUE
        ))
        entry <- protDiaDaIndexEntry(index, "KO_vs_WT")
        expect_true(grepl(
            built$descriptor$identity$workflow_slug,
            entry$query_specification$query_id,
            fixed = TRUE
        ))
        expect_identical(
            entry$query_specification$max_rows,
            .PROT_DIA_DA_MAX_QUERY_ROWS
        )
        expect_identical(
            entry$query_specification$max_bytes,
            .PROT_DIA_DA_MAX_QUERY_BYTES
        )
        contrast_manifest <- protDiaDaReadJson(
            artifactResolveContainedPath(
                root,
                entry$manifest_relative_path,
                must_exist = TRUE
            ),
            \(value) protDiaDaValidateContrastManifest(
                value,
                protDiaDaStore(built$workflow$workflow_context)
            )
        )
        expect_true(file.exists(artifactResolveContainedPath(
            root,
            contrast_manifest$model$relative_path,
            must_exist = TRUE
        )))
        first <- queryProtDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            projections = c(
                "Protein.Ids", "gene_name", "log2FC", "fdr_qvalue"
            ),
            sort_id = "absolute_effect",
            direction = "desc",
            limit = 3L
        )
        expect_identical(nrow(first$data), 3L)
        expect_true(first$has_more)
        second <- queryProtDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            projections = c(
                "Protein.Ids", "gene_name", "log2FC", "fdr_qvalue"
            ),
            sort_id = "absolute_effect",
            direction = "desc",
            cursor = first$next_cursor,
            limit = 3L
        )
        expect_false(any(first$data$Protein.Ids %in% second$data$Protein.Ids))
        quoted <- queryProtDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            filters = list(protein_search = list(
                operator = "contains",
                value = "%_';--"
            )),
            limit = 10L
        )
        expect_identical(quoted$data$Protein.Ids, "P%_';--")
        empty <- queryProtDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            filters = list(protein_search = list(
                operator = "contains",
                value = "not-present"
            )),
            limit = 10L
        )
        expect_identical(nrow(empty$data), 0L)
        heatmap <- protDaSelectedResults(
            built$workflow,
            list2env(list(
                da_results_list = index,
                current_s4_object = built$protein
            ), parent = emptyenv()),
            "KO_vs_WT",
            view = "heatmap",
            q_value_threshold = 0.05,
            top_n = 4L
        )
        expect_lte(nrow(heatmap$da_proteins_long), 4L)
        withr::with_options(list(
            multischolar.prot_nondia.da_queries = "disabled"
        ), {
            expect_error(
                queryProtDaPage(
                    built$workflow,
                    index,
                    "KO_vs_WT",
                    limit = 1L
                ),
                class = "multischolar_prot_nondia_da_queries_disabled"
            )
        })
    }
})

test_that("non-DIA DA failures and globals cannot advance current", {
    nondiaArtifact026SkipDependencies()
    scenario <- nondiaArtifact026Scenarios()[[1L]]
    root <- withr::local_tempdir(pattern = "nondia-026-failure-")
    built <- nondiaArtifact026Build(root, scenario)
    withr::defer(built$workflow$state_manager$close())
    results <- nondiaArtifact026Results(built)
    parameters <- protDaParameters("~ 0 + group", 0.05, 0, TRUE, TRUE)
    old_global <- get0("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
    had_global <- exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
    withr::defer({
        if (had_global) {
            assign("contrasts_tbl", old_global, envir = .GlobalEnv)
        } else if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
            rm("contrasts_tbl", envir = .GlobalEnv)
        }
    })
    assign(
        "contrasts_tbl",
        data.frame(contrasts = "hostile-hostile"),
        envir = .GlobalEnv
    )
    expect_identical(
        protDaResolveContrasts(built$workflow),
        normaliseProtDaContrastsTable(built$workflow$contrasts_tbl)
    )
    expect_identical(
        protDaResolveAnnotations(built$workflow),
        built$workflow$uniprot_dat_cln
    )
    expect_error(
        prepareProtDaArtifactRun(
            built$workflow,
            built$workflow$contrasts_tbl,
            results,
            parameters,
            failure_injector = \(stage, ...) {
                if (identical(stage, "after_contrast_commit")) {
                    stop("injected partial DA failure")
                }
            }
        ),
        "injected partial DA failure",
        fixed = TRUE
    )
    current_path <- artifactResolveContainedPath(
        root,
        protDiaDaPaths(built$workflow$workflow_context)$current
    )
    expect_false(file.exists(current_path))
    prepared <- prepareProtDaArtifactRun(
        built$workflow,
        built$workflow$contrasts_tbl,
        results,
        parameters
    )
    child <- built$protein
    child@args$stale_source <- TRUE
    built$workflow$state_manager$saveState(
        "da_stale_child",
        child,
        child@args,
        "stale source child"
    )
    expect_error(
        publishProtDaArtifactRun(built$workflow, prepared),
        class = "multischolar_prot_dia_da_source_mismatch"
    )
    expect_false(file.exists(current_path))
    withr::with_options(list(
        multischolar.prot_nondia.da_persistence = "disabled"
    ), {
        disabled <- prepareProtDaArtifactRun(
            built$workflow,
            built$workflow$contrasts_tbl,
            results,
            parameters
        )
        expect_false(disabled$enabled)
    })
})
