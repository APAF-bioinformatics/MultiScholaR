# Protein QC declarations (KIT-015). Each DIA protein QC step is a compute-only entry point that hydrates
# the saved protein state (2,500 proteins wide at 300 samples), runs the unchanged kernel once and returns
# the compute-only protein contract `finishDIAReferenceProteinStep` and the controller's completion consume:
# `checkpoint`, `reference`, `state_name`, `description`, `config`, `section`, `protein_count`, `execution`,
# the display `text`, the step's QC parameter record (`qc_name`, `qc`) and the fields its output product
# carries. The GUI runs the same entry point in process; the managed worker runs it from a frozen request.

#' Hydrate a saved protein state for a whole-placement kernel.
#' @noRd
readDIAReferenceProteinObject <- function(manager, reference) {
    object <- readDIAExecutionComponents(manager, reference)
    if (!methods::is(object, "ProteinQuantitativeData"))
        cli::cli_abort("a protein state is required", class = "workflow_execution_unsupported")
    object
}

#' Checkpoint a kernel's protein object as the step's output, inheriting the saved input as its source.
#' @noRd
checkpointDIAReferenceProteinObject <- function(manager, reference, object) {
    object <- inheritDIAExecutionSource(object, manager, reference)
    checkpoint <- checkpointDIAExecutionResult(manager, object)
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, checkpoint$dependencies,
        checkpoint$components, checkpoint$support, checkpoint$products))
    checkpoint
}

#' The configuration environment a kernel's parameter updates write to: the session's bus in the GUI, the
#' global copy in a script that has one, else a private copy seeded from the saved state's configuration
#' (KIT-009: a worker process creates no process-wide binding).
#' @noRd
diaReferenceConfigEnvironment <- function(before) {
    env <- workflowConfigEnvironment()
    if (exists("config_list", envir = env, inherits = FALSE)) return(env)
    env <- new.env(parent = emptyenv())
    env$config_list <- before$config %||% list()
    env
}

#' IQ protein rollup on a saved imputed peptide state (KIT-015): the peptide table streamed to IQ's input
#' file with the aliased sample ids, `iq::process_long_format` (MaxLFQ, `normalization = "none"`), the
#' protein object built as the module builds it. Placement `whole`: IQ reads its own file and estimates
#' each protein from its own rows; the per-group check on the fixture is recorded in the ticket's evidence.
#' @param manager,reference Owner and the saved imputed peptide state.
#' @param parameters None declared; `peptide_qc_dir` and `protein_qc_dir` name the working directories
#'   the IQ input and output files are written to (the GUI's science request directories, the worker's
#'   request directory).
#' @param worker_cap Admitted worker memory in bytes (recorded; IQ bounds its own memory).
#' @param variant `whole`.
#' @param process_long_format,write_fn,read_fn Existing seams of the module step.
#' @return Compute-only protein contract.
#' @noRd
computeDIAReferenceIQRollup <- function(manager, reference, parameters = list(), worker_cap = 12 * 1024^3,
    variant = "whole", process_long_format = iq::process_long_format, write_fn = vroom::vroom_write,
    read_fn = vroom::vroom) {
    if (!identical(variant, "whole")) cli::cli_abort("unsupported IQ rollup variant", class = "workflow_execution_unsupported")
    paths <- parameters[c("peptide_qc_dir", "protein_qc_dir")]
    if (!is.character(paths$peptide_qc_dir) || !is.character(paths$protein_qc_dir))
        cli::cli_abort("IQ working directories required", class = "workflow_execution_invalid")
    for (path in unlist(paths)) dir.create(path, showWarnings = FALSE, recursive = TRUE)
    before <- readDIAStateMetadata(manager, reference)
    if (!identical(before$data$class, "PeptideQuantitativeData"))
        cli::cli_abort("IQ rollup needs a peptide state", class = "workflow_execution_unsupported")
    peptideS4 <- readDIAExecutionComponents(manager, reference, fields = "design_matrix")
    proteinColumn <- peptideS4@protein_id_column
    peptideColumn <- peptideS4@peptide_sequence_column
    sampleColumn <- peptideS4@sample_id
    intensityColumn <- "Peptide.Imputed"
    missing <- setdiff(c(proteinColumn, peptideColumn, sampleColumn, intensityColumn), names(peptideS4@peptide_data))
    if (length(missing)) cli::cli_abort("IQ rollup input is missing required column(s): {missing}", class = "workflow_execution_invalid")
    if (!sampleColumn %in% names(peptideS4@design_matrix))
        cli::cli_abort("the declared sample column is absent from the design matrix", class = "workflow_execution_invalid")
    peptideValuesImputedFile <- file.path(paths$peptide_qc_dir, "peptide_values_imputed.tsv")
    originalSamples <- unique(as.character(peptideS4@design_matrix[[sampleColumn]]))
    sampleMapping <- tibble::tibble(Original = originalSamples, Alias = paste0("S_", sprintf("%03d", seq_along(originalSamples))))
    peptideReference <- before$data$components$peptide_data
    writeDIAIQInput(manager, peptideReference, peptideValuesImputedFile, sampleMapping, sampleColumn, intensityColumn, write_fn)
    lineage <- readDIAProteinLineage(manager, peptideReference, proteinColumn)
    iqOutputFile <- file.path(paths$protein_qc_dir, "iq_output_file.txt")
    process_long_format(peptideValuesImputedFile, output_filename = iqOutputFile, sample_id = sampleColumn,
        primary_id = proteinColumn, secondary_id = peptideColumn, intensity_col = intensityColumn,
        filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"), normalization = "none")
    if (!file.exists(iqOutputFile)) cli::cli_abort("IQ output file not created", class = "workflow_execution_invalid")
    proteinLog2Quant <- read_fn(iqOutputFile, .name_repair = "minimal")
    restored <- vapply(colnames(proteinLog2Quant), function(col) {
        if (col == proteinColumn) return(col)
        match <- sampleMapping$Original[sampleMapping$Alias == col]
        if (length(match)) match[[1L]] else col
    }, character(1), USE.NAMES = FALSE)
    colnames(proteinLog2Quant) <- restored
    if (!proteinColumn %in% names(proteinLog2Quant))
        cli::cli_abort("IQ output is missing the active protein key {.val {proteinColumn}}", class = "workflow_execution_invalid")
    survivingSamples <- intersect(originalSamples, setdiff(colnames(proteinLog2Quant), proteinColumn))
    finalDesignMatrix <- dplyr::filter(peptideS4@design_matrix, !!rlang::sym(sampleColumn) %in% survivingSamples)
    droppedSamples <- setdiff(originalSamples, survivingSamples)
    lineage_view <- peptideS4
    lineage_view@args$protein_accession_provenance <- lineage
    proteinObj <- ProteinQuantitativeData(protein_quant_table = proteinLog2Quant, protein_id_column = proteinColumn,
        protein_id_table = .proteinIdTableFromPeptideLineage(lineage_view, proteinLog2Quant[[proteinColumn]]),
        design_matrix = finalDesignMatrix, sample_id = sampleColumn, group_id = peptideS4@group_id,
        technical_replicate_id = peptideS4@technical_replicate_id, args = peptideS4@args)
    attr(proteinObj, "dia_component_source") <- attr(peptideS4, "dia_component_source")
    attr(proteinObj, "dia_execution_fields") <- c("protein_quant_table", "protein_id_table", "design_matrix")
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, proteinObj)
    protein_count <- length(unique(proteinLog2Quant[[proteinColumn]]))
    samples <- ncol(proteinLog2Quant) - 1L
    design_columns <- colnames(finalDesignMatrix)
    proteinObj <- peptideS4 <- proteinLog2Quant <- NULL
    text <- paste(
        "IQ Protein Rollup & S4 Object Creation Completed Successfully\n",
        "============================================================\n",
        sprintf("Proteins quantified: %d\n", protein_count),
        sprintf("Samples: %d\n", samples),
        sprintf("Active protein key: %s\n", proteinColumn),
        "Algorithm: MaxLFQ (via IQ tool)\n",
        "S4 Class: ProteinQuantitativeData\n",
        sprintf("Design matrix: %s\n", paste(design_columns, collapse = ", ")),
        sprintf("Output file: %s\n", basename(iqOutputFile)),
        "State saved as: 'protein_s4_created'\n",
        "\nReady for protein accession cleanup.")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "protein_s4_created", description = "IQ protein rollup completed and ProteinQuantitativeData S4 object created",
        config = list(iq_output_file = iqOutputFile, peptide_input_file = peptideValuesImputedFile,
            s4_class = "ProteinQuantitativeData", protein_id_column = proteinColumn),
        section = "protein_rollup", protein_count = protein_count,
        execution = list(variant = variant, worker_cap = worker_cap, rows = protein_count, samples = samples,
            dropped_samples = droppedSamples),
        text = text, field = "proteinObj", qc_name = NULL, qc = NULL, dropped_samples = droppedSamples)
}

#' Accession cleanup on a saved protein state (KIT-015): the unchanged `chooseBestProteinAccession` kernel
#' (the DIA branch annotates `protein_id_table`; the quant table is unchanged) with the FASTA table, through
#' the module's injectable alias so the recorded arguments keep their name.
#' @param manager,reference Owner and the saved protein state.
#' @param parameters `delimiter`, `aggregation_method` and optionally `fasta_table` (path to the project's
#'   `aa_seq_tbl_final.RDS`, read when no table is handed in).
#' @param seqinr The FASTA table (`aa_seq_tbl_final`) when the caller holds it (the GUI's bus); NULL reads
#'   `parameters$fasta_table`, else the owner's source directory, else the cleanup is skipped as the module
#'   skips it without FASTA data.
#' @param fasta_metadata The import's FASTA metadata when the caller holds it (recorded in the result only).
#' @return Compute-only protein contract.
#' @noRd
computeDIAReferenceAccessionCleanup <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "whole", seqinr = NULL, fasta_metadata = NULL, now_fn = Sys.time) {
    if (!identical(variant, "whole")) cli::cli_abort("unsupported accession cleanup variant", class = "workflow_execution_unsupported")
    delimiter <- parameters$delimiter
    aggregation <- parameters$aggregation_method
    if (is.null(seqinr)) {
        path <- parameters$fasta_table %||% file.path(manager$.context$experiment_paths$source_dir %||% "", "aa_seq_tbl_final.RDS")
        if (is.character(path) && length(path) == 1L && nzchar(path) && file.exists(path)) seqinr <- readRDS(path)
        metadata_path <- file.path(dirname(path), "fasta_metadata.RDS")
        if (is.null(fasta_metadata) && file.exists(metadata_path)) fasta_metadata <- readRDS(metadata_path)
    }
    object <- readDIAReferenceProteinObject(manager, reference)
    key_before <- .proteinCleanupActiveKey(object)
    proteins_before <- nrow(dplyr::distinct(object@protein_quant_table, !!rlang::sym(key_before)))
    chooseBestProteinAccessionFn <- chooseBestProteinAccession
    if (!is.null(seqinr)) {
        seqinr <- dplyr::rename(seqinr, uniprot_acc = "database_id")
        cleaned <- chooseBestProteinAccessionFn(theObject = object, delim = delimiter, seqinr_obj = seqinr,
            seqinr_accession_column = "uniprot_acc", replace_zero_with_na = TRUE, aggregation_method = aggregation)
        applied <- TRUE
    } else {
        cleaned <- object
        applied <- FALSE
    }
    key_after <- .proteinCleanupActiveKey(cleaned)
    protein_count <- nrow(dplyr::distinct(cleaned@protein_quant_table, !!rlang::sym(key_after)))
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, cleaned)
    object <- cleaned <- NULL
    results <- list(cleanup_applied = applied, delimiter_used = delimiter, aggregation_method = aggregation,
        active_protein_key_before = key_before, active_protein_key_after = key_after,
        proteins_before = proteins_before, proteins_after = protein_count,
        had_full_metadata = if (!is.null(fasta_metadata)) isTRUE(fasta_metadata$has_protein_evidence) &&
            isTRUE(fasta_metadata$has_gene_names) else FALSE,
        timestamp = now_fn())
    text <- if (applied) paste(
        "Protein Accession Cleanup Applied Successfully\n",
        "==============================================\n",
        sprintf("Proteins remaining: %d\n", protein_count),
        sprintf("Active protein key: %s\n", key_after),
        sprintf("Delimiter: %s\n", delimiter),
        sprintf("Aggregation method: %s\n", aggregation),
        "State saved as: 'protein_accession_cleaned'\n") else paste(
        "Protein Accession Cleanup Skipped\n",
        "=================================\n",
        sprintf("Proteins remaining: %d\n", protein_count),
        sprintf("Active protein key: %s\n", key_after),
        "Reason: No FASTA data available for cleanup\n",
        "State saved as: 'protein_accession_cleaned'\n")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "protein_accession_cleaned",
        description = if (applied) "Applied protein accession cleanup" else "Skipped accession cleanup (no FASTA data)",
        config = list(delimiter = delimiter, aggregation_method = aggregation, cleanup_applied = applied,
            active_protein_key = key_after),
        section = "chooseBestProteinAccession", protein_count = protein_count,
        execution = list(variant = variant, worker_cap = worker_cap, rows = protein_count, cleanup_applied = applied),
        text = text, field = "cleanedS4", qc_name = "accession_cleanup", qc = results,
        accession_cleanup_results = results, cleanupApplied = applied)
}

#' Protein intensity (missing-value) filter on a saved protein state (KIT-015, whole matrix, audit A2): the
#' cutoffs from `updateMissingValueParameters` (flexible) or zero (strict), the percentile, then the
#' unchanged `removeRowsWithMissingValuesPercent` kernel (dataset-wide quantile) through its alias.
#' @param parameters `strict_mode`, `min_reps_per_group`, `min_groups`, `proteins_intensity_cutoff_percentile`.
#' @return Compute-only protein contract.
#' @noRd
computeDIAReferenceProteinIntensityFilter <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "whole", now_fn = Sys.time) {
    if (!identical(variant, "whole")) cli::cli_abort("unsupported protein intensity filter variant", class = "workflow_execution_unsupported")
    before <- readDIAStateMetadata(manager, reference)
    strict <- isTRUE(parameters$strict_mode)
    percentile <- parameters$proteins_intensity_cutoff_percentile
    object <- readDIAReferenceProteinObject(manager, reference)
    env <- diaReferenceConfigEnvironment(before)
    section <- "removeRowsWithMissingValuesPercent"
    if (strict) {
        object <- updateConfigParameter(object, function_name = section, parameter_name = "groupwise_percentage_cutoff", new_value = 0, env = env)
        object <- updateConfigParameter(object, function_name = section, parameter_name = "max_groups_percentage_cutoff", new_value = 0, env = env)
    } else {
        object <- updateMissingValueParameters(object, min_reps_per_group = parameters$min_reps_per_group,
            min_groups = parameters$min_groups, env = env)
    }
    object <- updateConfigParameter(object, function_name = section, parameter_name = "proteins_intensity_cutoff_percentile",
        new_value = percentile, env = env)
    groupwise <- object@args[[section]]$groupwise_percentage_cutoff
    max_groups <- object@args[[section]]$max_groups_percentage_cutoff
    removeRowsWithMissingValuesPercentFn <- removeRowsWithMissingValuesPercent
    filtered <- removeRowsWithMissingValuesPercentFn(object)
    protein_count <- nrow(dplyr::distinct(filtered@protein_quant_table, .data[[filtered@protein_id_column]]))
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, filtered)
    object <- filtered <- NULL
    config <- list(strict_mode = strict, min_reps_per_group = if (!strict) parameters$min_reps_per_group else NA,
        min_groups = if (!strict) parameters$min_groups else NA, groupwise_percentage_cutoff = groupwise,
        max_groups_percentage_cutoff = max_groups, proteins_intensity_cutoff_percentile = percentile)
    qc <- c(config, list(timestamp = now_fn()))
    text <- if (strict) paste(
        "Protein Intensity Filter Applied Successfully\n",
        "============================================\n",
        "Mode: STRICT (No Missing Values)\n",
        sprintf("Proteins remaining: %d\n", protein_count),
        "Groupwise % cutoff: 0.000% (strict - no missing allowed)\n",
        "Max groups % cutoff: 0.000% (strict - all groups must pass)\n",
        sprintf("Intensity cutoff percentile: %.1f%%\n", percentile),
        "State saved as: 'protein_intensity_filtered'\n") else paste(
        "Protein Intensity Filter Applied Successfully\n",
        "============================================\n",
        "Mode: FLEXIBLE (Adaptive Thresholds)\n",
        sprintf("Proteins remaining: %d\n", protein_count),
        sprintf("Min replicates per group: %d\n", as.integer(parameters$min_reps_per_group)),
        sprintf("Min groups required: %d\n", as.integer(parameters$min_groups)),
        sprintf("Groupwise %% cutoff: %.3f%% (calculated)\n", groupwise),
        sprintf("Max groups %% cutoff: %.3f%% (calculated)\n", max_groups),
        sprintf("Intensity cutoff percentile: %.1f%%\n", percentile),
        "State saved as: 'protein_intensity_filtered'\n")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "protein_intensity_filtered",
        description = if (strict) "Applied STRICT protein intensity filter (no missing values)" else
            "Applied FLEXIBLE protein intensity filter (adaptive thresholds)",
        config = config, section = section, protein_count = protein_count,
        execution = list(variant = variant, worker_cap = worker_cap, rows = protein_count,
            groupwise_percentage_cutoff = groupwise, max_groups_percentage_cutoff = max_groups),
        text = text, field = "filteredS4", qc_name = "intensity_filter", qc = qc)
}

#' The duplicate-protein resolution kernel, extracted verbatim from `runProteinDuplicateRemovalStep`
#' (KIT-015, audit A5): the duplicated protein ids, then one row per protein id with every numeric column
#' aggregated by `aggregation_fn(x, na.rm = TRUE)`.
#' @param table Wide protein table.
#' @param protein_id_column Grouping column (the active protein key).
#' @param aggregation_fn Aggregation function (`mean`, `median`, ...).
#' @return `list(table, duplicates)`.
#' @noRd
resolveProteinDuplicateRows <- function(table, protein_id_column, aggregation_fn) {
    duplicates <- table |>
        dplyr::group_by(.data[[protein_id_column]]) |>
        dplyr::filter(dplyr::n() > 1) |>
        dplyr::select(dplyr::all_of(protein_id_column)) |>
        dplyr::distinct() |>
        dplyr::pull(dplyr::all_of(protein_id_column))
    quantColumns <- setdiff(colnames(table), protein_id_column)
    numericQuantColumns <- quantColumns[vapply(table[quantColumns], is.numeric, logical(1))]
    if (length(numericQuantColumns) == 0L) {
        stop("No numeric protein intensity columns available for duplicate removal.")
    }
    table <- table |>
        dplyr::group_by(!!rlang::sym(protein_id_column)) |>
        dplyr::summarise(
            dplyr::across(dplyr::all_of(numericQuantColumns), ~ aggregation_fn(.x, na.rm = TRUE)),
            .groups = "drop"
        ) |>
        dplyr::ungroup()
    list(table = table, duplicates = duplicates)
}

#' Duplicate-protein removal on a saved protein state (KIT-015): the extracted kernel once on the whole
#' table (unit feature group by the active protein key; 2,500 rows at 300 samples, placed whole).
#' @param parameters `aggregation_method` (the name of the aggregation function).
#' @return Compute-only protein contract.
#' @noRd
computeDIAReferenceProteinDuplicateRemoval <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "whole", aggregation_resolver = base::get) {
    if (!identical(variant, "whole")) cli::cli_abort("unsupported duplicate removal variant", class = "workflow_execution_unsupported")
    method <- parameters$aggregation_method
    object <- readDIAReferenceProteinObject(manager, reference)
    key <- if ("protein_id_column" %in% methods::slotNames(object) && !is.null(object@protein_id_column) &&
        nzchar(object@protein_id_column)) object@protein_id_column else "Protein.Ids"
    resolved <- resolveProteinDuplicateRows(object@protein_quant_table, key, aggregation_resolver(method))
    object@protein_quant_table <- resolved$table
    protein_count <- nrow(dplyr::distinct(object@protein_quant_table, .data[[object@protein_id_column]]))
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, object)
    object <- NULL
    duplicates <- resolved$duplicates
    text <- paste(
        "Duplicate Protein Removal Completed Successfully\n",
        "===============================================\n",
        sprintf("Proteins remaining: %d\n", protein_count),
        sprintf("Duplicates found: %d\n", length(duplicates)),
        sprintf("Aggregation method: %s\n", method),
        "State saved as: 'duplicates_removed'\n")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "duplicates_removed", description = "Removed duplicate proteins by aggregation",
        config = list(aggregation_method = method, duplicates_found = duplicates, num_duplicates = length(duplicates)),
        section = "resolveProteinDuplicateRows", protein_count = protein_count,
        execution = list(variant = variant, worker_cap = worker_cap, rows = protein_count, duplicates = length(duplicates)),
        text = text, field = "deduplicatedS4", qc_name = NULL, qc = NULL, duplicates = duplicates)
}

#' Commit a computed protein result and record its bounded display product (KIT-015): the GUI in process
#' (`publish = TRUE`, the files and figures published at once) and the controller's completion of a worker
#' result (`publish = FALSE`, published afterwards by `publishWorkflowProteinResult`) share this step.
#' @param workflow Session or completion owner (`state_manager`, `qc_params`, optional `filtering_progress`,
#'   `accession_cleanup_results`).
#' @param value Compute-only protein contract.
#' @return The display handoff with `stateReference`.
#' @noRd
finishDIAReferenceProteinStep <- function(workflow, value, publish = TRUE) {
    manager <- workflow$state_manager
    commitDIAReferenceState(manager, value$state_name, value$checkpoint, value$config, value$description)
    object <- readDIAExecutionComponents(manager)
    if (!is.null(value$qc_name)) {
        qc_params <- shiny::isolate(workflow$qc_params) %||% list()
        if (is.null(qc_params$protein_qc)) qc_params$protein_qc <- list()
        qc_params$protein_qc[[value$qc_name]] <- value$qc
        workflow$qc_params <- qc_params
    }
    if (!is.null(value$accession_cleanup_results)) workflow$accession_cleanup_results <- value$accession_cleanup_results
    result <- list(resultText = value$text)
    result[[value$field]] <- object
    for (name in intersect(c("outputFile", "proteinCount", "cleanupApplied", "duplicates"), names(value))) result[[name]] <- value[[name]]
    if (is.null(result$proteinCount)) result$proteinCount <- value$protein_count
    finishDIAProteinAction(workflow, result, value$field, publish = publish)
}

#' Recover the bounded display product of one accepted managed protein action (KIT-015): the output
#' product the finishing step recorded on the head, with the step's QC entries merged into the bus.
#' @param manager Reference-backed owner.
#' @param result Accepted worker/controller result.
#' @param workflow_data Existing GUI workflow bus.
#' @return The display handoff (`stateReference`, text, summary, `qcParams`, `accessionCleanupResults`).
#' @noRd
readManagedProteinResult <- function(manager, result, workflow_data) {
    reference <- result$state_reference
    product <- manager$.context$snapshot$output_products[[reference$entity_id]]
    if (is.null(product)) cli::cli_abort("protein output product is missing", class = "workflow_execution_invalid")
    envelope <- hydrateSavedWorkflowEntity(manager$.context$store, product)
    checkSavedWorkflowConsumerPayload("scientific_product", envelope)
    value <- c(list(stateReference = reference), envelope$value[setdiff(names(envelope$value), "stateReference")])
    qc_params <- shiny::isolate(workflow_data$qc_params) %||% list()
    if (is.null(qc_params$protein_qc)) qc_params$protein_qc <- list()
    for (name in names(value$qcParams$protein_qc)) qc_params$protein_qc[[name]] <- value$qcParams$protein_qc[[name]]
    workflow_data$qc_params <- qc_params
    if (!is.null(value$accessionCleanupResults)) workflow_data$accession_cleanup_results <- value$accessionCleanupResults
    value
}

# Declarations ---------------------------------------------------------------------------------------------

#' The IQ protein rollup as a declaration (KIT-015): unit feature group (MaxLFQ estimates each protein from
#' its own rows; the surviving sample set is the global statistic), placed whole: IQ reads the streamed
#' input file itself. The kernel runs the compute-only entry point.
#' @noRd
diaProteinIQRollupOperation <- function() {
    newWorkflowOperation(name = "protein_iq_rollup", input_class = "PeptideQuantitativeData",
        fields = "design_matrix", unit = "feature_group",
        kernel = function(object, parameters, store) computeDIAReferenceIQRollup(store$manager, store$source,
            parameters, worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData", parameters = character(),
        product_type = "dia.protein_iq_rollup_result",
        global_statistics = c(surviving_sample_set =
            "the samples with a value in IQ's output over the whole run; the design is filtered once to them"),
        merge = c(design_matrix = "design", protein_id_table = "concatenate", protein_quant_table = "concatenate"))
}

#' The accession cleanup as a declaration (KIT-015): unit row of `protein_id_table` (the DIA branch ranks
#' accessions per protein group and annotates the id table; the quant table is unchanged), placed whole.
#' @noRd
diaProteinAccessionCleanupOperation <- function() {
    newWorkflowOperation(name = "protein_accession_cleanup", input_class = "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "row",
        kernel = function(object, parameters, store) computeDIAReferenceAccessionCleanup(store$manager, store$source,
            parameters[c("delimiter", "aggregation_method", "fasta_table")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole",
            seqinr = parameters$seqinr, fasta_metadata = parameters$fasta_metadata),
        output_class = "ProteinQuantitativeData", parameters = c("delimiter", "aggregation_method"),
        product_type = "dia.protein_accession_cleanup_result",
        merge = c(protein_id_table = "concatenate", protein_quant_table = "concatenate"))
}

#' The protein intensity filter as a declaration (KIT-015): whole matrix (dataset-wide quantile and
#' design-derived cutoffs, audit A2), the unchanged kernel once.
#' @noRd
diaProteinIntensityFilterOperation <- function() {
    newWorkflowOperation(name = "protein_intensity_filter", input_class = "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        kernel = function(object, parameters, store) computeDIAReferenceProteinIntensityFilter(store$manager, store$source,
            parameters[c("strict_mode", "min_reps_per_group", "min_groups", "proteins_intensity_cutoff_percentile")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData",
        parameters = c("strict_mode", "min_reps_per_group", "min_groups", "proteins_intensity_cutoff_percentile"),
        product_type = "dia.protein_intensity_filter_result",
        global_statistics = c(min_protein_intensity_threshold = "the percentile of every finite value in the whole table",
            group_sizes = "samples per group and the number of groups from the design"))
}

#' The duplicate-protein removal as a declaration (KIT-015): unit feature group by the active protein key
#' (the extracted kernel aggregates within a protein id), placed whole.
#' @noRd
diaProteinDuplicateRemovalOperation <- function() {
    newWorkflowOperation(name = "protein_duplicate_removal", input_class = "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "feature_group",
        kernel = function(object, parameters, store) computeDIAReferenceProteinDuplicateRemoval(store$manager, store$source,
            parameters["aggregation_method"], worker_cap = parameters$worker_cap %||% 12 * 1024^3,
            variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData", parameters = "aggregation_method",
        product_type = "dia.protein_duplicate_removal_result",
        merge = c(protein_quant_table = "concatenate", counts = "sum"))
}
