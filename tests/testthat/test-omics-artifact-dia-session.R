omics015SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

omics015Paths <- function(root, project_id = "omics-015-project") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "proteomics",
        omic_label = "dia-session-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

omics015Workflow <- function(paths) {
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
        differential_expression = "disabled"
    )
    workflow
}

omics015Build <- function(root) {
    paths <- omics015Paths(root)
    source <- file.path(paths$source_dir, "report.tsv")
    fixture <- testthat::test_path(
        "..", "testdata", "e2e", "prot_dia", "report.tsv"
    )
    stopifnot(file.copy(fixture, source))
    workflow <- omics015Workflow(paths)
    imported <- suppressMessages(importDIANNData(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- "diann"
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType("DIA")
    import_result <- persistProtDiaImportArtifacts(
        workflow,
        imported,
        source
    )
    stopifnot(isTRUE(import_result$ok))

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
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT",
        stringsAsFactors = FALSE
    )
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
    design_result <- persistProtDiaDesignArtifacts(workflow)
    stopifnot(isTRUE(design_result$ok))

    manager <- newProtDiaArtifactStateManager(workflow$workflow_context)
    workflow$state_manager <- manager
    protein_ids <- paste0("P", seq_len(14L))
    values <- matrix(
        seq_len(length(protein_ids) * length(runs)),
        nrow = length(protein_ids),
        ncol = length(runs),
        dimnames = list(protein_ids, runs)
    )
    protein_table <- data.frame(
        Protein.Group = rownames(values),
        as.data.frame(values, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(protein_table) <- NULL
    args <- list(
        globalParameters = list(
            workflow_type = "DIA",
            report_template = "DIANN_report.rmd",
            use_limpa = FALSE
        ),
        deAnalysisParameters = list(formula_string = "~ 0 + group")
    )
    protein <- ProteinQuantitativeData(
        protein_quant_table = protein_table,
        protein_id_column = "Protein.Group",
        protein_id_table = data.frame(
            Protein.Group = protein_table$Protein.Group,
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
            record_id = "session:correlation_filtered",
            stage_id = "correlation_filter"
        )
    )
    workflow$config_list <- args
    workflow$fasta_metadata <- list(
        fasta_format = "UniProt",
        num_sequences = length(proteins)
    )
    workflow$accession_cleanup_results <- list(
        cleanup_applied = TRUE,
        aggregation_method = "sum"
    )
    workflow$ruv_optimization_result <- list(
        best_k = 2L,
        best_percentage = 10,
        optimization_results = data.frame(
            percentage = c(5, 10),
            score = c(0.2, 0.5)
        )
    )
    workflow$qc_params <- list(protein_qc = list(minimum = 1))
    workflow$protein_counts <- list(final = nrow(protein_table))
    workflow$mixed_species_analysis <- list(enabled = FALSE)
    norm_data <- new.env(parent = emptyenv())
    norm_data$correlation_filtered_obj <- protein
    norm_data$correlation_filtering_complete <- TRUE
    norm_data$best_k <- 2L
    norm_data$correlation_threshold <- 0.75
    norm_data$qc_plot_paths <- list()
    list(
        paths = paths,
        workflow = workflow,
        norm_data = norm_data,
        protein = protein,
        annotations = workflow$uniprot_dat_cln
    )
}

omics015Export <- function(built) {
    runProtNormExportSessionWorkflow(
        workflowData = built$workflow,
        normData = built$norm_data,
        input = list(norm_method = "none", ruv_mode = "manual"),
        sourceDir = built$paths$source_dir,
        withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) invisible(NULL),
        collectSessionDataFn = function(...) {
            collectProtNormExportSessionData(
                workflowData = built$workflow,
                normData = built$norm_data,
                input = list(norm_method = "none", ruv_mode = "manual"),
                timeFn = function() {
                    as.POSIXct("2026-08-21 03:04:05", tz = "UTC")
                },
                messageFn = function(...) invisible(NULL)
            )
        },
        messageFn = function(...) invisible(NULL)
    )
}

omics015ExpectS4Exact <- function(expected, actual) {
    expect_identical(class(actual), class(expected))
    for (slot_name in methods::slotNames(expected)) {
        expect_identical(
            methods::slot(actual, slot_name),
            methods::slot(expected, slot_name),
            info = slot_name
        )
    }
    expect_identical(methods::validObject(actual, test = TRUE), TRUE)
}

test_that("portable metadata and paths reject runtime-only values", {
    expect_error(
        workflowSessionEncodeMetadata(list(path = "/private/data.tsv")),
        class = "multischolar_absolute_path_in_artifact_state"
    )
    expect_error(
        workflowSessionEncodeMetadata(list(cache = new.env())),
        class = "multischolar_unsafe_artifact_state_metadata"
    )
    root <- withr::local_tempdir(pattern = "omics-015-path-")
    expect_error(
        workflowSessionProjectRelativePath(tempdir(), root),
        class = "multischolar_workflow_session_path_escape"
    )
})

test_that("DIA-NN export publishes a bounded generation-pinned session", {
    omics015SkipDependencies()
    root <- withr::local_tempdir(pattern = "omics-015-export-")
    built <- omics015Build(root)
    withr::defer(built$workflow$state_manager$close())
    exported <- omics015Export(built)
    artifact <- exported$exportArtifacts$artifactSession

    expect_true(isTRUE(artifact$ok))
    expect_true(file.exists(exported$exportArtifacts$latestFilepath))
    expect_true(file.exists(artifact$sessionFilepath))
    expect_true(file.exists(artifact$latestFilepath))
    manifest <- readWorkflowSessionManifest(artifact$latestFilepath)
    expect_identical(
        manifest$workflow_state$current_generation_id,
        built$workflow$state_manager$getCurrentGenerationId()
    )
    expect_identical(
        manifest$compatibility$byte_digest,
        artifactByteDigest(exported$exportArtifacts$sessionFilepath)
    )
    expect_setequal(names(manifest$stage_refs), c("import", "design"))
    expect_identical(
        names(manifest$stage_refs$design$refs),
        .PROT_DIA_SESSION_DESIGN_ROLES
    )
    raw_json <- paste(readLines(artifact$latestFilepath), collapse = "\n")
    expect_false(grepl(root, raw_json, fixed = TRUE))
    expect_false(grepl("r6_complete_states", raw_json, fixed = TRUE))
    expect_false(grepl("current_s4_object", raw_json, fixed = TRUE))
    expect_false(grepl("connection", raw_json, ignore.case = TRUE))
    legacy <- readRDS(exported$exportArtifacts$latestFilepath)
    expect_false("workflow_state_manifest" %in% names(legacy))
    expect_true(isS4(
        legacy$r6_complete_states[[legacy$r6_current_state_name]]
    ))
})

test_that("moved artifact session restores exact DA input without transforms", {
    omics015SkipDependencies()
    parent <- withr::local_tempdir(pattern = "omics-015-move-parent-")
    original <- file.path(parent, "original")
    moved <- file.path(parent, "moved")
    dir.create(original)
    built <- omics015Build(original)
    exported <- omics015Export(built)
    oracle <- built$protein
    annotations <- built$annotations
    legacy <- readRDS(exported$exportArtifacts$latestFilepath)
    built$workflow$state_manager$close()
    expect_true(file.rename(original, moved))

    moved_paths <- omics015Paths(moved, built$paths$project_id)
    manifest_path <- file.path(
        moved_paths$source_dir,
        "filtered_session_artifact_latest.json"
    )
    child_output <- file.path(
        moved_paths$source_dir,
        "fresh_process_restore.rds"
    )
    child_script <- file.path(parent, "restore_session.R")
    repository <- normalizePath(
        testthat::test_path("..", ".."),
        winslash = "/",
        mustWork = TRUE
    )
    writeLines(c(
        sprintf("pkgload::load_all(%s, quiet = TRUE)", dQuote(repository)),
        sprintf("manifest_path <- %s", dQuote(manifest_path)),
        sprintf("project_root <- %s", dQuote(moved)),
        sprintf("source_dir <- %s", dQuote(moved_paths$source_dir)),
        sprintf("results_dir <- %s", dQuote(moved_paths$results_dir)),
        sprintf("project_id <- %s", dQuote(moved_paths$project_id)),
        "paths <- list(base_dir = project_root, project_id = project_id,",
        "    omic_type = 'proteomics', omic_label = 'dia-session-study',",
        "    source_dir = source_dir, results_dir = results_dir)",
        paste0(
            "bundle <- MultiScholaR:::restoreProtDiaSessionManifest(",
            "manifest_path, paths)"
        ),
        "object <- bundle$session_data$current_s4_object",
        "result <- list(object = object, design = bundle$session_data$design_matrix,",
        "    contrasts = bundle$session_data$contrasts_tbl,",
        "    config = bundle$session_data$config_list,",
        "    annotations = bundle$uniprot_dat_cln,",
        "    valid = identical(methods::validObject(object, test = TRUE), TRUE))",
        sprintf("saveRDS(result, %s)", dQuote(child_output)),
        "bundle$state_manager$close()"
    ), child_script)
    child_status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", child_script),
        stdout = TRUE,
        stderr = TRUE
    )
    expect_null(attr(child_status, "status"), info = paste(child_status, collapse = "\n"))
    child <- readRDS(child_output)
    expect_true(child$valid)
    omics015ExpectS4Exact(oracle, child$object)
    expect_identical(child$design, oracle@design_matrix)
    expect_identical(child$contrasts, legacy$contrasts_tbl)
    expect_identical(child$config, oracle@args)
    expect_identical(child$annotations, annotations)

    restored <- restoreProtDiaSessionManifest(manifest_path, moved_paths)
    withr::defer(restored$state_manager$close())
    omics015ExpectS4Exact(oracle, restored$session_data$current_s4_object)
    expect_identical(restored$session_data$design_matrix, oracle@design_matrix)
    expect_identical(restored$session_data$contrasts_tbl, legacy$contrasts_tbl)
    expect_identical(restored$session_data$config_list, oracle@args)
    expect_identical(restored$uniprot_dat_cln, annotations)
    da_validation <- validateProtDaModelAndContrasts(
        restored$session_data$current_s4_object,
        "~ 0 + group",
        restored$session_data$contrasts_tbl
    )
    expect_identical(da_validation$contrast_strings, "KO_vs_WT=groupKO-groupWT")

    compatibility_path <- file.path(
        moved_paths$source_dir,
        "reconstructed_session.rds"
    )
    fingerprint <- writeProtDiaCompatibilitySession(
        restored,
        compatibility_path
    )
    expect_identical(
        fingerprint$generation_id,
        restored$manifest$workflow_state$current_generation_id
    )
    reconstructed <- readRDS(compatibility_path)
    memory <- WorkflowState$new()
    restored_snapshot <- restoreWorkflowStateFromSession(memory, reconstructed)
    expect_identical(
        restored_snapshot$r6_current_state_name,
        "correlation_filtered"
    )
    omics015ExpectS4Exact(oracle, memory$getState())
})

test_that("artifact failure falls back to usable legacy without mutation", {
    session_data <- module_ci_prot_da_session_payload()
    object <- session_data$current_s4_object
    paths <- module_ci_prot_da_paths()
    withr::defer(unlink(dirname(paths$source_dir), recursive = TRUE))
    module_ci_prot_da_write_session(paths$source_dir, session_data)
    writeLines(
        "{not valid json",
        file.path(paths$source_dir, "filtered_session_artifact_latest.json")
    )
    bundle <- resolveProtDaSessionBundle(paths)
    expect_identical(bundle$format, "legacy")
    expect_s3_class(bundle$artifact_error, "error")
    expect_identical(bundle$session_data$current_s4_object, object)

    withr::local_options(list(
        multischolar.prot_dia.artifact_session_restore = "disabled"
    ))
    disabled <- resolveProtDaSessionBundle(paths)
    expect_identical(disabled$format, "legacy")
    expect_null(disabled$artifact_error)
})

test_that("failed session application restores manager, reactives, and globals", {
    session_data <- module_ci_prot_da_session_payload()
    manager <- new.env(parent = emptyenv())
    manager$states <- list(original = list(marker = TRUE))
    manager$state_history <- "original"
    manager$current_state <- "original"
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$tab_status <- list(
        normalization = "pending",
        differential_expression = "locked"
    )
    workflow$state_update_trigger <- NULL
    workflow$contrasts_tbl <- NULL
    workflow$design_matrix <- NULL
    workflow$config_list <- list(original = TRUE)
    workflow$uniprot_dat_cln <- NULL
    workflow$aa_seq_tbl_final <- NULL
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- NULL
    da_data$contrasts_available <- NULL
    da_data$formula_from_s4 <- NULL
    bundle <- list(
        format = "legacy",
        session_data = session_data,
        state_manager = NULL,
        workflow_context = NULL,
        uniprot_dat_cln = NULL,
        aa_seq_tbl_final = NULL
    )
    old_config <- if (exists("config_list", envir = .GlobalEnv)) {
        get("config_list", envir = .GlobalEnv)
    } else {
        NULL
    }
    had_config <- exists("config_list", envir = .GlobalEnv)
    withr::defer({
        if (exists("config_list", envir = .GlobalEnv)) {
            rm("config_list", envir = .GlobalEnv)
        }
        if (had_config) assign("config_list", old_config, envir = .GlobalEnv)
    })
    expect_error(
        applyProtDaSessionBundle(
            workflow,
            da_data,
            bundle,
            failure_injector = function(stage, context) {
                if (identical(stage, "before_session_apply_commit")) {
                    stop("injected apply failure")
                }
                invisible(context)
            }
        ),
        "injected apply failure"
    )
    expect_identical(manager$current_state, "original")
    expect_identical(manager$state_history, "original")
    expect_identical(manager$states, list(original = list(marker = TRUE)))
    expect_null(da_data$current_s4_object)
    expect_identical(workflow$config_list, list(original = TRUE))
    expect_identical(workflow$tab_status$normalization, "pending")
    expect_identical(workflow$tab_status$differential_expression, "locked")
})

test_that("failed latest publication preserves the previous portable session", {
    omics015SkipDependencies()
    root <- withr::local_tempdir(pattern = "omics-015-rollback-")
    built <- omics015Build(root)
    withr::defer(built$workflow$state_manager$close())
    exported <- omics015Export(built)
    latest <- exported$exportArtifacts$artifactSession$latestFilepath
    prior_digest <- artifactByteDigest(latest)
    prior_manifest <- readWorkflowSessionManifest(latest)
    second_rds <- file.path(root, "source", "compatibility_second.rds")
    saveRDS(exported$sessionData, second_rds)
    second_manifest <- buildProtDiaSessionManifest(
        built$workflow,
        built$norm_data,
        exported$sessionData,
        second_rds
    )
    expect_error(
        workflowSessionReplaceLatest(
            second_manifest,
            latest,
            failure_injector = function(stage, context) {
                if (identical(stage, "before_latest_publication")) {
                    stop("injected publication failure")
                }
                invisible(context)
            }
        ),
        "injected publication failure"
    )
    expect_identical(artifactByteDigest(latest), prior_digest)
    expect_identical(readWorkflowSessionManifest(latest), prior_manifest)
})
