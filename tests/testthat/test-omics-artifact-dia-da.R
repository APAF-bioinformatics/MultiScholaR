omics016SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics016TempRoot <- function() {
    root <- tempfile("omics-016-")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    root
}

omics016Build <- function(root, protein_count = 14L) {
    paths <- list(
        base_dir = root,
        project_id = "omics-016-project",
        omic_type = "proteomics",
        omic_label = "dia-da-study",
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
    workflow$tab_status <- list(
        normalization = "complete",
        differential_expression = "pending",
        differential_abundance = "pending",
        enrichment_analysis = "locked"
    )
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
    workflow$uniprot_dat_cln <- data.frame(
        Protein.Group = unique(workflow$data_cln$Protein.Group),
        Gene = paste0("GENE", seq_along(unique(
            workflow$data_cln$Protein.Group
        ))),
        stringsAsFactors = FALSE
    )
    workflow$aa_seq_tbl_final <- data.frame(
        accession = workflow$uniprot_dat_cln$Protein.Group,
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
        "raw_data_s4",
        peptide,
        workflow$config_list,
        "DIA-NN design memory checkpoint."
    )
    stopifnot(isTRUE(persistProtDiaDesignArtifacts(workflow)$ok))

    manager <- newProtDiaArtifactStateManager(context)
    workflow$state_manager <- manager
    stopifnot(protein_count >= 2L)
    protein_ids <- c(
        paste0("P", seq_len(protein_count - 1L)),
        "P%_';--"
    )
    values <- matrix(
        seq_len(length(protein_ids) * length(runs)),
        nrow = length(protein_ids),
        dimnames = list(protein_ids, runs)
    )
    protein_table <- data.frame(
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
        protein_quant_table = protein_table,
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
        "Final DIA-NN DA handoff.",
        audit_metadata = list(
            record_id = "da:correlation_filtered",
            stage_id = "correlation_filter"
        )
    )
    workflow$config_list <- args
    list(workflow = workflow, paths = paths, protein = protein)
}

omics016Result <- function(object, comparison, full_format, reverse = FALSE) {
    ids <- object@protein_quant_table[[object@protein_id_column]]
    count <- length(ids)
    effect <- rep(c(2, -2, 0.5, -0.5), length.out = count)
    if (isTRUE(reverse)) effect <- -effect
    q_values <- rep(c(0.01, 0.01, 0.04, 0.2), length.out = count)
    q_values[[5L]] <- NA_real_
    table <- data.frame(
        Protein.Group = ids,
        gene_names = paste0("GENE", seq_along(ids)),
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
            analysis_type = "RUV applied",
            stringsAsFactors = FALSE
        ))
    )
}

omics016Results <- function(built) {
    contrasts <- built$workflow$contrasts_tbl
    results <- list(
        omics016Result(
            built$protein,
            contrasts$friendly_names[[1L]],
            contrasts$full_format[[1L]]
        ),
        omics016Result(
            built$protein,
            contrasts$friendly_names[[2L]],
            contrasts$full_format[[2L]],
            reverse = TRUE
        )
    )
    names(results) <- contrasts$contrasts
    results
}

omics016Parameters <- function(treat = 0) {
    protDiaDaParameters("~ 0 + group", 0.05, treat, TRUE, TRUE)
}

omics016ContainsHeavyPayload <- function(value) {
    if (is.data.frame(value) || is.matrix(value) || isS4(value) ||
        is.environment(value) || inherits(value, "DBIConnection")) {
        return(TRUE)
    }
    if (!is.list(value)) return(FALSE)
    any(vapply(value, omics016ContainsHeavyPayload, logical(1)))
}

omics016TemporaryQueryEntries <- function(store) {
    path <- artifactStoreResolveFile(store, store$relative_paths$duckdb_tmp)
    if (!dir.exists(path)) return(character())
    setdiff(list.files(path, all.files = TRUE), c(".", ".."))
}

test_that("DIA-NN DA run persists exact contrast tables and model provenance", {
    omics016SkipDependencies()
    built <- omics016Build(omics016TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    results <- omics016Results(built)

    prepared <- prepareProtDiaDaArtifactRun(
        built$workflow,
        built$workflow$contrasts_tbl,
        results,
        omics016Parameters(),
        now = as.POSIXct("2026-08-21 06:00:00", tz = "UTC")
    )
    expect_true(prepared$enabled)
    expect_null(restoreProtDiaDaArtifactIndex(built$workflow))
    published <- publishProtDiaDaArtifactRun(built$workflow, prepared)
    expect_true(published$published)
    index <- restoreProtDiaDaArtifactIndex(built$workflow)
    expect_true(isProtDiaDaArtifactIndex(index))
    expect_length(index$contrasts, 2L)

    raw_manifest <- paste(
        readLines(published$current_path, warn = FALSE),
        collapse = "\n"
    )
    expect_false(grepl(built$paths$base_dir, raw_manifest, fixed = TRUE))
    expect_false(grepl("protein_quant_table", raw_manifest, fixed = TRUE))
    expect_false(grepl("coefficients", raw_manifest, fixed = TRUE))
    expect_false(omics016ContainsHeavyPayload(index))

    store <- protDiaDaStore(built$workflow$workflow_context)
    for (contrast_index in seq_along(index$contrasts)) {
        entry <- index$contrasts[[contrast_index]]
        expected <- results[[contrast_index]]$da_proteins_long
        ref <- artifactStoreNormalizeRef(entry$long_table$ref)
        sidecar <- artifactStoreReadSidecar(
            store,
            artifactStoreManagedPaths(
                store, ref$logical_key, ref$artifact_id
            )$sidecar,
            validate_payload = TRUE
        )
        actual <- decodeArtifactRectangular(
            arrow::read_parquet(
                artifactStoreResolveFile(
                    store, ref$relative_path, must_exist = TRUE
                ),
                as_data_frame = FALSE
            ),
            sidecar$codec_metadata
        )
        expect_identical(actual, expected)
        contrast_entry <- prepared$manifest$contrasts[[contrast_index]]
        contrast_path <- artifactResolveContainedPath(
            built$paths$base_dir,
            contrast_entry$manifest_relative_path,
            must_exist = TRUE
        )
        contrast_manifest <- protDiaDaReadJson(
            contrast_path,
            function(value) protDiaDaValidateContrastManifest(value, store)
        )
        model_path <- artifactResolveContainedPath(
            built$paths$base_dir,
            contrast_manifest$model$relative_path,
            must_exist = TRUE
        )
        expect_identical(
            readRDS(model_path),
            results[[contrast_index]]$contrasts_results
        )
    }
    expect_identical(
        prepared$manifest$source$generation_id,
        built$workflow$state_manager$getCurrentGenerationId()
    )
    expect_setequal(
        names(prepared$manifest$source$stage_refs$design$refs),
        .PROT_DIA_SESSION_DESIGN_ROLES
    )
    forged <- prepared$manifest
    forged$run_id <- paste0("darun_", strrep("0", 64L))
    forged$manifest_digest <- protDiaDaJsonDigest(forged)
    expect_error(
        protDiaDaValidateRunManifest(
            forged, built$workflow$workflow_context
        ),
        class = "multischolar_prot_dia_da_manifest_digest_mismatch"
    )
    forged <- prepared$manifest
    forged$descriptor_contract$descriptor_id <- "forged.descriptor"
    forged$manifest_digest <- protDiaDaJsonDigest(forged)
    expect_error(
        protDiaDaValidateRunManifest(
            forged, built$workflow$workflow_context
        ),
        class = "multischolar_prot_dia_da_identity_mismatch"
    )
    forged <- prepared$manifest
    forged$contrasts[[1L]]$manifest_relative_path <-
        forged$contrasts[[2L]]$manifest_relative_path
    forged$contrasts[[1L]]$manifest_digest <-
        forged$contrasts[[2L]]$manifest_digest
    forged$manifest_digest <- protDiaDaJsonDigest(forged)
    expect_error(
        protDiaDaValidateRunManifest(
            forged, built$workflow$workflow_context
        ),
        class = "multischolar_prot_dia_da_contrast_mismatch"
    )
})

test_that("bounded pages preserve ties and treat SQL-like text as data", {
    omics016SkipDependencies()
    built <- omics016Build(omics016TempRoot(), protein_count = 513L)
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    results <- omics016Results(built)
    prepared <- prepareProtDiaDaArtifactRun(
        built$workflow,
        built$workflow$contrasts_tbl,
        results,
        omics016Parameters()
    )
    index <- publishProtDiaDaArtifactRun(built$workflow, prepared)$index
    store <- protDiaDaStore(built$workflow$workflow_context)
    temporary_entries <- omics016TemporaryQueryEntries(store)

    cursor <- NULL
    first_cursor <- NULL
    pages <- list()
    repeat {
        page <- queryProtDiaDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            projections = c("Protein.Group", "fdr_qvalue"),
            sort_id = "q_value",
            direction = "asc",
            cursor = cursor,
            limit = 37L
        )
        pages[[length(pages) + 1L]] <- page$data
        if (length(pages) == 1L) first_cursor <- page$next_cursor
        if (!page$has_more) break
        cursor <- page$next_cursor
    }
    actual <- dplyr::bind_rows(pages)
    expected <- results[[1L]]$da_proteins_long
    expected <- expected[order(
        is.na(expected$fdr_qvalue),
        expected$fdr_qvalue,
        seq_len(nrow(expected))
    ), c("Protein.Group", "fdr_qvalue")]
    rownames(expected) <- NULL
    expect_identical(actual, expected)
    expect_identical(anyDuplicated(actual$Protein.Group), 0L)
    forged_cursor <- first_cursor
    forged_cursor$artifact_id <- paste0("artifact_", strrep("0", 64L))
    expect_error(
        queryProtDiaDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            projections = c("Protein.Group", "fdr_qvalue"),
            sort_id = "q_value",
            direction = "asc",
            cursor = forged_cursor,
            limit = 37L
        ),
        class = "multischolar_invalid_artifact_query_cursor"
    )

    special <- queryProtDiaDaPage(
        built$workflow,
        index,
        "KO_vs_WT",
        projections = "Protein.Group",
        filters = list(protein_search = list(
            operator = "contains",
            value = "%_';--"
        )),
        limit = 10L
    )
    expect_identical(special$data$Protein.Group, "P%_';--")
    empty <- queryProtDiaDaPage(
        built$workflow,
        index,
        "KO_vs_WT",
        projections = "Protein.Group",
        filters = list(protein_search = list(
            operator = "contains",
            value = "definitely-absent-protein"
        )),
        limit = 10L
    )
    expect_identical(nrow(empty$data), 0L)
    expect_false(empty$has_more)
    expect_null(empty$next_cursor)
    expect_error(
        queryProtDiaDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            limit = 11L,
            resource_policy = list(max_result_rows = 10L)
        ),
        class = "multischolar_artifact_query_row_limit_exceeded"
    )
    expect_error(
        queryProtDiaDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            limit = 10L,
            resource_policy = list(max_result_bytes = 1L)
        ),
        class = "multischolar_artifact_query_byte_limit_exceeded"
    )
    expect_error(
        queryProtDiaDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            projections = "not_a_column"
        ),
        class = "multischolar_invalid_artifact_query_projection"
    )
    expect_error(
        queryProtDiaDaPage(
            built$workflow,
            index,
            "KO_vs_WT",
            sort_id = "raw_sql"
        ),
        class = "multischolar_invalid_artifact_query_sort"
    )
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- built$protein
    da_data$da_results_list <- index
    volcano <- protDiaDaSelectedResults(
        built$workflow,
        da_data,
        "KO_vs_WT",
        view = "volcano"
    )
    expect_identical(volcano$da_proteins_long, results[[1L]]$da_proteins_long)
    heatmap <- protDiaDaSelectedResults(
        built$workflow,
        da_data,
        "KO_vs_WT",
        view = "heatmap",
        q_value_threshold = 0.05,
        top_n = 25L
    )
    expect_identical(nrow(heatmap$da_proteins_long), 25L)
    expect_true(all(heatmap$da_proteins_long$fdr_qvalue < 0.05))
    expect_identical(
        protDiaDaArtifactSummary(
            built$workflow,
            index,
            "KO_vs_WT",
            q_value = 0.025,
            lfc = 1
        ),
        protDiaDaSummary(
            results[[1L]]$da_proteins_long,
            list(da_q_val_thresh = 0.025, treat_lfc_cutoff = 1)
        )
    )
    expect_identical(
        omics016TemporaryQueryEntries(store),
        temporary_entries
    )
})

test_that("interrupted contrast commit does not advance current and resumes", {
    omics016SkipDependencies()
    built <- omics016Build(omics016TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    results <- omics016Results(built)
    first <- prepareProtDiaDaArtifactRun(
        built$workflow,
        built$workflow$contrasts_tbl,
        results,
        omics016Parameters()
    )
    first <- publishProtDiaDaArtifactRun(built$workflow, first)
    prior_digest <- artifactByteDigest(first$current_path)

    interrupted <- function(stage, context) {
        if (identical(stage, "after_contrast_commit")) {
            stop("injected contrast interruption", call. = FALSE)
        }
        invisible(context)
    }
    expect_error(
        prepareProtDiaDaArtifactRun(
            built$workflow,
            built$workflow$contrasts_tbl,
            results,
            omics016Parameters(treat = 0.5),
            failure_injector = interrupted
        ),
        "injected contrast interruption"
    )
    expect_identical(artifactByteDigest(first$current_path), prior_digest)
    expect_identical(
        restoreProtDiaDaArtifactIndex(built$workflow)$run_id,
        first$index$run_id
    )

    resumed <- prepareProtDiaDaArtifactRun(
        built$workflow,
        built$workflow$contrasts_tbl,
        results,
        omics016Parameters(treat = 0.5)
    )
    expect_true(resumed$enabled)
    publication_interruption <- function(stage, context) {
        if (identical(stage, "before_da_current_publication")) {
            stop("injected publication interruption", call. = FALSE)
        }
        invisible(context)
    }
    expect_error(
        publishProtDiaDaArtifactRun(
            built$workflow,
            resumed,
            failure_injector = publication_interruption
        ),
        "injected publication interruption"
    )
    expect_identical(artifactByteDigest(first$current_path), prior_digest)
    expect_identical(
        restoreProtDiaDaArtifactIndex(built$workflow)$run_id,
        first$index$run_id
    )
    resumed <- publishProtDiaDaArtifactRun(built$workflow, resumed)
    expect_false(identical(resumed$index$run_id, first$index$run_id))
    expect_identical(
        restoreProtDiaDaArtifactIndex(built$workflow)$run_id,
        resumed$index$run_id
    )
})

test_that("artifact science ignores stale globals and flags roll back independently", {
    omics016SkipDependencies()
    built <- omics016Build(omics016TempRoot())
    on.exit(built$workflow$state_manager$close(), add = TRUE)
    stale <- data.frame(
        contrasts = "groupSTALE-groupWRONG",
        stringsAsFactors = FALSE
    )
    old <- if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
        get("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
    } else {
        NULL
    }
    assign("contrasts_tbl", stale, envir = .GlobalEnv)
    on.exit({
        rm("contrasts_tbl", envir = .GlobalEnv)
        if (!is.null(old)) assign("contrasts_tbl", old, envir = .GlobalEnv)
    }, add = TRUE)
    expect_identical(
        protDiaDaResolveContrasts(built$workflow),
        normaliseProtDaContrastsTable(built$workflow$contrasts_tbl)
    )

    withr::local_options(list(
        multischolar.prot_dia.da_persistence = "disabled",
        multischolar.prot_dia.da_queries = "enabled"
    ))
    expect_false(protDiaDaArtifactEligible(built$workflow, "persistence"))
    expect_true(protDiaDaArtifactEligible(built$workflow, "queries"))
    expect_identical(
        protDiaDaResolveContrasts(built$workflow),
        normaliseProtDaContrastsTable(built$workflow$contrasts_tbl)
    )
    withr::local_options(list(
        multischolar.prot_dia.da_queries = "disabled"
    ))
    expect_false(protDiaDaArtifactEligible(built$workflow, "queries"))
})
