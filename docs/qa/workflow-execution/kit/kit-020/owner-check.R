# KIT-020 owner check (in process): the three protein-level fixtures (MaxQuant and FragPipe from
# tests/testdata/e2e/prot_lfq, Proteome Discoverer TMT from tests/testdata/e2e/prot_tmt) imported on the on-disk
# owner as the GUI does (initializeDIANativeImport, registerDIAInputFiles, the reader, applyProtImportResultToWorkflow,
# completeDIAInputFiles), the design committed through buildProtDesignStateCheckpoint, then the shared protein QC
# declarations (intensity, duplicates, replicate), normalisation (RUV skipped), correlation and DA driven as the
# replay drives them; the record validated and replayed exactly; a native export attempted and reopened.
#   identity: protein_s4_initial hydrated equals ProteinQuantitativeData built in memory from the same table.
# Usage: R_LIBS=<proof library> Rscript docs/qa/workflow-execution/kit/kit-020/owner-check.R [maxquant fragpipe pd_tmt]
suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
repo <- normalizePath(".")
h <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
only <- commandArgs(TRUE)
fixtures <- file.path(repo, "tests/testdata/e2e")
cases <- list(
    maxquant = list(file = file.path(fixtures, "prot_lfq", "proteinGroups.txt"), fasta = file.path(fixtures, "prot_lfq", "proteins.fasta"),
        design = file.path(fixtures, "prot_lfq", "design.tsv"), read = function(p) importMaxQuantData(p), type = "LFQ"),
    fragpipe = list(file = file.path(fixtures, "prot_lfq", "seed_combined_protein.tsv"), fasta = file.path(fixtures, "prot_lfq", "proteins.fasta"),
        design = file.path(fixtures, "prot_lfq", "design.tsv"), read = function(p) importFragPipeData(p), type = "LFQ"),
    pd_tmt = list(file = file.path(fixtures, "prot_tmt", "pd_tmt_peptides.tsv"), fasta = file.path(fixtures, "prot_tmt", "proteins.fasta"),
        design = file.path(fixtures, "prot_tmt", "design.tsv"), read = function(p) importProteomeDiscovererTMTData(p), type = "TMT"))
if (length(only)) cases <- cases[only]
fake_session <- function() {
    session <- new.env(parent = emptyenv()); session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit020-owner-check")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
config <- readConfigFile(system.file("config", "config.ini", package = "MultiScholaR"), create_workers = FALSE)
results <- list()
run_case <- function(name, case) {
    out <- list(format = name, type = case$type)
    base <- file.path(work, name, "project"); dir.create(base, recursive = TRUE)
    project_dirs <- setupDirectories(base_dir = base, omic_types = "proteomics", label = name, force = TRUE, reuse_existing = TRUE)
    assign("project_dirs", project_dirs, envir = globalenv())
    paths <- project_dirs[[1L]]
    paths$time_dir <- paths$time_dir %||% file.path(paths$publication_graphs_dir, "filtering_qc")
    dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$config_list <- config
    workflow$config_list$globalParameters$workflow_type <- case$type
    workflow$qc_params <- list(protein_qc = list())
    workflow$filtering_progress <- NULL
    initializeDIANativeImport(workflow, paths, fake_session(), name, "iq")
    manager <- workflow$state_manager
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    stopifnot(isDIANativeSession(manager))
    out$store_identity <- manager$.context$store$identity[c("project_id", "workflow_id")]
    design_in <- utils::read.delim(case$design, stringsAsFactors = FALSE)
    design <- NULL
    in_memory <- NULL
    withDIAStateScope(manager, function() {
        inputs <- c(case$file, case$fasta); names(inputs) <- c(paste0(name, "_report"), "fasta")
        registerDIAInputFiles(manager, inputs)
        imported <- case$read(case$file)
        type <- applyProtImportResultToWorkflow(workflow, imported, name, case$fasta, logInfo = function(...) NULL,
            showNotification = function(...) NULL)
        out$workflow_type <<- type
        metadata <- getDIAImportMetadata(manager)
        out$import_metadata <<- metadata[c("data_type", "column_mapping")]
        completeDIAInputFiles(manager)
        # The design is built on the imported runs, as the GUI's design builder does (the TMT reader names runs
        # "126_WT_1"; the fixture design names samples "WT_1", a suffix of the run).
        runs <- as.character(metadata$samples$Run)
        group <- vapply(runs, function(run) design_in$group[[which(endsWith(run, design_in$sample))[[1L]]]], "")
        design <<- data.frame(Run = runs, group = unname(group), stringsAsFactors = FALSE)
        design$replicates <<- stats::ave(seq_len(nrow(design)), design$group, FUN = seq_along)
        workflow$design_matrix <- design
        # The in-memory construction on the same table, for the identity check.
        table <- getDIAOwnedPayload(workflow, "data_cln")
        mapping <- workflow$column_mapping
        wide <- tidyr::pivot_wider(table, id_cols = !!rlang::sym(mapping$protein_col), names_from = !!rlang::sym(mapping$run_col),
            values_from = !!rlang::sym(mapping$quantity_col))
        wide <- dplyr::mutate(wide, dplyr::across(-!!rlang::sym(mapping$protein_col), ~ log2(.x + 1)))
        design_tr <- design; design_tr$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
        in_memory <<- ProteinQuantitativeData(protein_quant_table = wide, protein_id_column = mapping$protein_col,
            protein_id_table = dplyr::distinct(wide, !!rlang::sym(mapping$protein_col)), design_matrix = design_tr,
            sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group", args = workflow$config_list)
        state <- buildProtDesignStateCheckpoint(workflow, case$type, "Import")
        out$design_state <<- state
    })
    restored <- manager$getState("protein_s4_initial")
    strip <- function(x) { x@args <- list(); x }
    out$design_identity <- list(
        tables_identical = identical(restored@protein_quant_table, in_memory@protein_quant_table),
        design_identical = identical(as.data.frame(restored@design_matrix), as.data.frame(in_memory@design_matrix)),
        slots_identical = identical(strip(restored), strip(in_memory)),
        args_identical = identical(restored@args, in_memory@args))
    out$design_provenance <- manager$.context$snapshot$step_provenance$protein_s4_initial[c("operation", "fixed")]
    out$inputs <- lapply(Filter(function(d) !is.null(d$role), manager$.context$snapshot$file_dependencies), function(d) d[c("role", "sha256")])
    # Protein QC through the registry handlers, as the replay drives them.
    steps <- list(
        list(operation = "protein_intensity_filter", parameters = list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L,
            proteins_intensity_cutoff_percentile = 1L)),
        list(operation = "protein_duplicate_removal", parameters = list(aggregation_method = "mean")),
        list(operation = "protein_replicate_filter", parameters = list(grouping_variable = "group", parallel_cores = 1L)))
    out$steps <- list()
    for (step in steps) {
        entry <- workflowOperationEntry(step$operation)
        variant <- entry$default_variant
        withDIAStateScope(manager, function() {
            directory <- file.path(work, name, "requests", step$operation); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
            request <- list(operation = step$operation, input_references = list(state = getDIAStateReference(manager)),
                parameters = step$parameters, worker_cap = 12 * 1024^3, directory = directory)
            value <- get(entry$handler, mode = "function")(manager, request, variant)
            value$operation <- step$operation
            setDIAStepProvenance(manager, step$operation, variant, entry$parameters, entry$declaration$fixed, values = step$parameters)
            completeDIAProteinAction(workflow, value)
        })
        record <- readDIAStateMetadata(manager)
        out$steps[[step$operation]] <- list(state = manager$current_state, rows = record$dimensions$rows %||% NA, columns = record$dimensions$columns %||% NA)
    }
    # Normalisation (RUV skipped on a five-protein fixture), correlation and DA, as the replay commits them.
    norm_parameters <- list(norm_method = "cyclicloess", ruv_mode = "skip", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 5, ruv_k = NA, grouping = "group")
    normalization_final <- NULL
    withDIAStateScope(manager, function() {
        run <- runSavedOperation(manager, diaNormalizationOperation(), getDIAStateReference(manager), norm_parameters)
        saved <- run$value
        final <- readDIAStateMetadata(manager, saved$final)$data
        config_n <- list(norm_method = "cyclicloess", ruv_mode = "skip", ruv_k = saved$best_k, percentage_as_neg_ctrl = saved$percentage, ruv_applied = FALSE)
        final$products$normalization_science <- run$reference
        commitDIAReferenceState(manager, "normalized", final, config_n, "Post-normalization complete: RUV-III skipped by user")
        normalization_final <<- saved$final
        pairs_parameters <- list(group = "group", tech_rep_remove_regex = "pool")
        pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), normalization_final, pairs_parameters)
        selection_parameters <- list(pairs = pairs$reference$entity_id, threshold = 0.5, correlations = pairs$value)
        selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), normalization_final, selection_parameters,
            dependencies = list(pairs$reference))
        spec <- readDIAStateMetadata(manager, selection$value)$data
        commitDIAReferenceState(manager, "correlation_filtered", spec, list(min_pearson_correlation_threshold = 0.5),
            "Applied final sample correlation filter (chunk 28)")
    })
    out$normalisation_state <- "normalized"; out$correlation_state <- manager$current_state
    contrasts <- data.frame(contrasts = "groupKO-groupWT", friendly_names = "KO_vs_WT", full_format = "KO_vs_WT=groupKO-groupWT", stringsAsFactors = FALSE)
    da_data <- new.env(parent = emptyenv()); da_data$current_s4_object <- getDIAStateReference(manager); da_data$contrasts_available <- contrasts$contrasts
    da_workflow <- new.env(parent = emptyenv()); da_workflow$state_manager <- manager; da_workflow$contrasts_tbl <- contrasts
    workflow$contrasts_tbl <- contrasts
    da <- tryCatch({
        runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da_data, da_workflow, paths,
            contrasts = contrasts, publish = FALSE)
        manifest <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
        list(ok = TRUE, contrasts = names(manifest$results), rows = vapply(manifest$results, function(r) r$rows %||% NA_integer_, 1))
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    out$da <- da
    # A native export of the owner carries the head's type in its manifest and reopens (KIT-020 item 4).
    export <- tryCatch({
        norm_data <- new.env(parent = emptyenv()); norm_data$best_k <- NA; norm_data$correlation_threshold <- 0.5
        metadata <- collectDIANativeSessionMetadata(workflow, norm_data, list(norm_method = "cyclicloess", ruv_mode = "skip"))
        output <- file.path(work, name, "export"); dir.create(output)
        exported <- runDIANativeSessionExport(manager, metadata, output)
        archive <- exported$publication$archive_path
        reopen_root <- file.path(work, name, "reopen"); dir.create(reopen_root)
        reopened <- openDIANativeSession(archive, reopen_root)
        on.exit(closeDIAStateAdapter(reopened$manager), add = TRUE)
        manifest <- jsonlite::read_json(file.path(dirname(archive), "filtered_session_latest.json"))
        list(ok = TRUE, archive = basename(archive), manifest_workflow = manifest$workflow %||% NA, manifest_route = manifest$route %||% NA,
            reopened_type = reopened$manager$workflow_type, reopened_state = reopened$manager$current_state)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    out$export <- export
    # The record, validated; its replay on a pruned copy of the store.
    json <- file.path(paths$source_dir, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    check <- validateWorkflowAnalysisRecord(record)
    out$record <- list(valid = check$valid, errors = check$errors, legacy = any(grepl("legacy step", check$incomplete)),
        workflow = record$workflow, states = vapply(record$steps, function(s) s$state, ""), operations = vapply(record$steps, function(s) s$operation, ""))
    store <- manager$.context$store
    out$store_fields <- names(store)
    origin_root <- store$root %||% dirname(store$directory)
    closeDIAStateAdapter(manager); on.exit()
    replay <- tryCatch({
        copy <- file.path(work, name, "replay"); dir.create(copy)
        stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(origin_root, "."), copy)) == 0L)
        states <- out$record$states; operations <- out$record$operations
        preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
        root <- h$kitPrepareReplayProject(copy, preconditions, api)
        scratch <- file.path(work, name, "replay-scratch"); dir.create(scratch)
        inputs <- c(case$file, case$fasta); names(inputs) <- c(paste0(name, "_report"), "fasta")
        r <- h$kitReplayJourney(record, root, origin_root, inputs, paths, scratch, api)
        list(ok = TRUE, states = length(r$states), mismatches = r$mismatches)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    out$replay <- replay
    out
}
for (name in names(cases)) {
    results[[name]] <- tryCatch(run_case(name, cases[[name]]), error = function(e) list(format = name, error = conditionMessage(e)))
    r <- results[[name]]
    if (!is.null(r$error)) { cat(name, "ERROR:", r$error, "\n"); next }
    cat(name, "| type", r$workflow_type, "| store", r$store_identity$project_id, r$store_identity$workflow_id, "| design", r$design_state,
        "identity tables", r$design_identity$tables_identical, "slots", r$design_identity$slots_identical, "args", r$design_identity$args_identical,
        "| provenance", r$design_provenance$operation, "| steps", paste(vapply(r$steps, function(s) s$state, ""), collapse = ","),
        "| corr", r$correlation_state, "| DA", if (isTRUE(r$da$ok)) paste(r$da$contrasts, r$da$rows) else paste("ERR", r$da$error),
        "| record valid", r$record$valid, "legacy", r$record$legacy, "workflow", r$record$workflow,
        "| replay", if (isTRUE(r$replay$ok)) paste(r$replay$states, "states", length(r$replay$mismatches), "mismatches") else paste("ERR", r$replay$error),
        "| export", if (isTRUE(r$export$ok)) paste(r$export$manifest_workflow, r$export$manifest_route, "reopened", r$export$reopened_type, r$export$reopened_state) else paste("ERR", r$export$error), "\n")
}
results$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
results$head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
jsonlite::write_json(results, file.path(repo, "docs/qa/workflow-execution/kit/kit-020", "owner-check.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("done", round(results$seconds), "s\n")
