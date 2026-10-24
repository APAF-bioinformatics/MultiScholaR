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
    # The declaration passes no QC plot file (audit A6, KIT-015 pre-execution audit PE-01): IQ's default
    # `qc-plots.pdf` would land in the process working directory (the worker's home); the plot is no
    # component or product of the step.
    fixed <- diaProteinIQRollupFixedParameters()
    process_long_format(peptideValuesImputedFile, output_filename = iqOutputFile, sample_id = sampleColumn,
        primary_id = proteinColumn, secondary_id = peptideColumn, intensity_col = intensityColumn,
        filter_double_less = fixed$filter_double_less, normalization = fixed$normalization, pdf_out = NULL)
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
    # KIT-026: the FASTA the import hashed (the head's `fasta` input) unless the caller declares one.
    fasta_sha256 <- parameters$fasta_sha256 %||% manager$.context$snapshot$file_dependencies$fasta$sha256 %||% NA_character_
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
            active_protein_key = key_after, fasta_sha256 = fasta_sha256),
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
    setDIAStepProvenanceForSection(manager, value$section) # KIT-025: the in-process path names its operation
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
#' The IQ rollup's fixed arguments (KIT-025): the values the kernel always passes to `iq::process_long_format`,
#' declared once so the head's provenance records them (KIT-027 labels them `fixed`); not request parameters.
#' @noRd
diaProteinIQRollupFixedParameters <- function() list(filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
    normalization = "none")

diaProteinIQRollupOperation <- function() {
    newWorkflowOperation(name = "protein_iq_rollup", input_class = "PeptideQuantitativeData",
        fields = "design_matrix", unit = "feature_group",
        kernel = function(object, parameters, store) computeDIAReferenceIQRollup(store$manager, store$source,
            parameters, worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData", parameters = character(), fixed = diaProteinIQRollupFixedParameters(),
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
            parameters[c("delimiter", "aggregation_method", "fasta_table", "fasta_sha256")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole",
            seqinr = parameters$seqinr, fasta_metadata = parameters$fasta_metadata),
        output_class = "ProteinQuantitativeData", parameters = c("delimiter", "aggregation_method", "fasta_sha256"),
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

# ---------------------------------------------------------------------------------------------------------------
# The DIA_LIMPA peptide stage (KIT-047): the limpa script's chunks 13, 16 and 17 as declarations on the owner. Each
# compute-only entry point hydrates the saved peptide state, runs the unchanged kernel once and returns the contract
# the shared finishing steps consume (`finishDIAReferencePeptideStep` for the two peptide outputs,
# `finishDIAReferenceProteinStep` for the DPC quantification); the GUI runs the same entry point in process, the
# managed worker from a frozen request.

#' The log2 peptide matrix as a declaration (chunk 13): unit row (`calcPeptideMatrix` pivots each feature's own rows;
#' `log2TransformPeptideMatrix` maps each cell, values <= 0 to NA), placed whole. No parameters; the kernel's constant
#' is fixed.
#' @noRd
diaPeptideLog2MatrixOperation <- function() {
    newWorkflowOperation(name = "peptide_log2_matrix", input_class = "PeptideQuantitativeData",
        fields = c("peptide_data", "design_matrix"), unit = "row",
        kernel = function(object, parameters, store) computeDIAReferencePeptideLog2Matrix(store$manager, store$source,
            parameters, worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "streamed"),
        output_class = "PeptideQuantitativeData", parameters = character(), fixed = list(nonpositive_to_na = TRUE),
        product_type = "dia.peptide_log2_matrix_result",
        merge = c(design_matrix = "design", peptide_data = "concatenate"))
}

#' The peptide-level accession cleanup as a declaration (chunk 16): unit feature group by the protein group (the
#' kernel resolves one accession per group and aggregates duplicate peptides within it), placed whole.
#' @noRd
diaPeptideAccessionCleanupOperation <- function() {
    newWorkflowOperation(name = "peptide_accession_cleanup", input_class = "PeptideQuantitativeData",
        fields = c("peptide_data", "design_matrix"), unit = "feature_group",
        kernel = function(object, parameters, store) computeDIAReferencePeptideAccessionCleanup(store$manager, store$source,
            parameters[c("delimiter", "aggregation_method", "seqinr_accession_column", "fasta_table", "fasta_sha256")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "partitioned",
            seqinr = parameters$seqinr, fasta_metadata = parameters$fasta_metadata),
        output_class = "PeptideQuantitativeData",
        parameters = c("delimiter", "aggregation_method", "seqinr_accession_column", "fasta_sha256"),
        product_type = "dia.peptide_accession_cleanup_result",
        merge = c(design_matrix = "design", peptide_data = "concatenate"))
}

#' The DPC imputation and quantification as a declaration (chunk 17): whole matrix (the detection probability curve
#' is fitted over every peptide and sample; the hyperparameters over every protein), the unchanged peptide method of
#' `proteinMissingValueImputationLimpa` once with the declared engine.
#' @noRd
diaProteinDpcQuantOperation <- function() {
    newWorkflowOperation(name = "protein_dpc_quant", input_class = "PeptideQuantitativeData",
        fields = c("peptide_data", "design_matrix", "peptide_matrix"), unit = "whole_matrix",
        kernel = function(object, parameters, store) computeDIAReferenceProteinDpcQuant(store$manager, store$source,
            parameters[c("engine", "dpc_slope", "chunk")], worker_cap = parameters$worker_cap %||% 12 * 1024^3,
            variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData", parameters = c("engine", "dpc_slope", "chunk"),
        product_type = "dia.protein_dpc_quant_result")
}

#' Checkpoint a kernel's peptide object as a peptide step's output, inheriting the saved input as its source.
#' @param loaded The slots the kernel rewrote (the writer keeps the stored component of any other slot).
#' @noRd
checkpointDIAReferencePeptideObject <- function(manager, reference, object, loaded) {
    attr(object, "dia_component_source") <- reference
    attr(object, "dia_execution_fields") <- loaded
    checkpoint <- checkpointDIAExecutionResult(manager, object)
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, checkpoint$dependencies,
        checkpoint$components, checkpoint$support, checkpoint$products))
    checkpoint
}

#' @noRd
readDIAReferencePeptideObject <- function(manager, reference, fields = c("peptide_data", "design_matrix")) {
    before <- readDIAStateMetadata(manager, reference)
    if (!identical(before$data$class, "PeptideQuantitativeData"))
        cli::cli_abort("a peptide state is required", class = "workflow_execution_unsupported")
    readDIAExecutionComponents(manager, reference, fields = fields)
}

#' Log2 peptide matrix on a saved peptide state (chunk 13): `calcPeptideMatrix` then `log2TransformPeptideMatrix`.
#' @return Compute-only peptide contract (`checkpoint`, `state_name`, `description`, `config`, `section`, `counts`).
#' @noRd
computeDIAReferencePeptideLog2Matrix <- function(manager, reference, parameters = list(), worker_cap = 12 * 1024^3,
    variant = "streamed") {
    if (!is.character(variant) || length(variant) != 1L || !variant %in% c("whole", "streamed"))
        cli::cli_abort("unsupported log2 matrix variant", class = "workflow_execution_unsupported")
    if (identical(variant, "streamed")) return(computeDIAReferencePeptideLog2MatrixStreamed(manager, reference, worker_cap))
    object <- readDIAReferencePeptideObject(manager, reference)
    if (isTRUE(object@is_logged_data))
        cli::cli_abort("the peptide state is already log2 transformed", class = "workflow_execution_invalid")
    calcPeptideMatrixFn <- calcPeptideMatrix
    log2TransformPeptideMatrixFn <- log2TransformPeptideMatrix
    object <- calcPeptideMatrixFn(object)
    nonpositive <- sum(object@peptide_matrix <= 0, na.rm = TRUE)
    missing_before <- sum(is.na(object@peptide_matrix))
    transformed <- log2TransformPeptideMatrixFn(object)
    object <- NULL
    peptides <- nrow(transformed@peptide_matrix)
    samples <- ncol(transformed@peptide_matrix)
    missing_after <- sum(is.na(transformed@peptide_matrix))
    quantity <- transformed@norm_quantity_column
    checkpoint <- checkpointDIAReferencePeptideObject(manager, reference, transformed,
        c("peptide_data", "design_matrix", "peptide_matrix"))
    transformed <- NULL
    config <- list(peptides = peptides, samples = samples, quantity_column = quantity, nonpositive_set_to_na = nonpositive,
        missing_before = missing_before, missing_after = missing_after,
        outcomes = list(peptides = peptides, samples = samples, nonpositive_set_to_na = nonpositive,
            missing_cells = missing_after))
    # The saved peptide reference beside the checkpoint: the controller accepts a worker's result only with it (the
    # DIA peptide steps' convention; KIT-053's lifecycle run found the worker path refusing the result without it).
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data, state_name = "log2_matrix",
        description = "Peptide matrix calculated and log2 transformed (values <= 0 set to NA)",
        config = config, section = "log2TransformPeptideMatrix", counts = NULL,
        execution = list(variant = variant, worker_cap = worker_cap, rows = peptides, samples = samples))
}

#' The log2 matrix on the kit's feature-group partitioner (KIT-057, audit F-01): the row unit realised. The whole
#' variant hydrates the long peptide table (5.5 million rows on the 300-sample corpus) and pivots it in one
#' `pivot_wider` (worker peak 4.09 GiB for a 142 MiB matrix). Here the table stays on disk: one pass over the
#' identity columns gives the feature-key map and the sample order exactly as the kernel builds them on the whole
#' table (the map from the distinct protein/peptide pairs in first-appearance order; `pivot_wider`'s columns in the
#' samples' first-appearance order); then batches of whole protein groups run the unchanged kernels
#' (`calcPeptideMatrix`, `log2TransformPeptideMatrix`) and their rows are placed into the matrix by name while the
#' kernel's updated quantity column streams through the bounded table writer as a change map over the source rows.
#' The design is cleaned once at the end as the whole kernel does. Admission as the partitioner records it.
#' @return The whole variant's compute-only contract (the checkpoint's `peptide_data` is the streamed table).
#' @noRd
computeDIAReferencePeptideLog2MatrixStreamed <- function(manager, reference, worker_cap = 12 * 1024^3) {
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) || worker_cap <= 0)
        cli::cli_abort("invalid admitted worker capacity", class = "workflow_execution_unsupported")
    context <- manager$.context
    before <- readDIAStateMetadata(manager, reference)
    if (!identical(before$data$class, "PeptideQuantitativeData"))
        cli::cli_abort("a peptide state is required", class = "workflow_execution_unsupported")
    object <- readDIAReferencePeptideObject(manager, reference, fields = "design_matrix")
    if (isTRUE(object@is_logged_data))
        cli::cli_abort("the peptide state is already log2 transformed", class = "workflow_execution_invalid")
    source <- before$data$components$peptide_data
    scope <- newWorkflowTableReadScope(context$store, source)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    con <- scope$connection
    engine_bytes <- min(1e9, floor(worker_cap / 8))
    DBI::dbExecute(con, paste0("SET memory_limit = '", floor(engine_bytes / 1048576), "MB'"))
    DBI::dbExecute(con, "SET threads = 2")
    if (any(c(".dia_source_order", ".dia_design_order", ".dia_feature_order") %in% scope$metadata$columns))
        cli::cli_abort("reserved DIA grouping columns are present", class = "saved_entity_invalid")
    protein <- object@protein_id_column
    peptide <- object@peptide_sequence_column
    sample <- object@sample_id
    quantity <- object@norm_quantity_column
    identity_columns <- intersect(unique(c(protein, peptide, sample, "Protein.Ids")), scope$metadata$columns)
    if (!all(c(protein, peptide, sample, quantity) %in% scope$metadata$columns))
        cli::cli_abort("the peptide table lacks the log2 step's columns", class = "saved_entity_invalid")
    # Prelude (one pass over the identity columns): the distinct identities in first-appearance order, as the
    # kernel's feature-key map and pivot see them on the whole table.
    identities <- NULL
    samples <- character()
    streamWorkflowTableRows(scope, identity_columns, function(page, offset, row_names) {
        if (!nrow(page)) return(invisible(NULL))
        rows <- page[identity_columns]
        for (name in identity_columns) rows[[name]] <- as.character(rows[[name]])
        rows <- rows[!duplicated(rows[c(protein, peptide)]), , drop = FALSE]
        identities <<- if (is.null(identities)) rows else {
            joined <- rbind(identities, rows)
            joined[!duplicated(joined[c(protein, peptide)]), , drop = FALSE]
        }
        samples <<- unique(c(samples, as.character(page[[sample]])))
        invisible(NULL)
    })
    if (is.null(identities) || !nrow(identities)) cli::cli_abort("the peptide table is empty", class = "saved_entity_invalid")
    rownames(identities) <- NULL
    # The kernel's own map and provenance on the distinct identities (the same functions the whole kernel calls on
    # the whole table; they read the two identity columns and `Protein.Ids` only).
    proxy <- object
    proxy@peptide_data <- identities
    proxy <- .ensurePeptideFeatureKeyMap(proxy)
    feature_map <- proxy@args$peptide_feature_key_map
    object@args <- proxy@args
    matrix_rows <- feature_map$feature_key
    assembled <- matrix(NA_real_, nrow = length(matrix_rows), ncol = length(samples),
        dimnames = list(matrix_rows, samples))
    offset <- 0L
    nonpositive <- 0L
    na_rows <- 0L
    kernel <- function(part, index) {
        positions <- attr(part, "workflow_source_order") %||% integer()
        attr(part, "workflow_source_order") <- NULL
        batch <- proxy
        batch@peptide_data <- part
        batch <- calcPeptideMatrix(batch)
        raw <- batch@peptide_matrix
        batch <- log2TransformPeptideMatrix(batch)
        list(raw = raw, logged = batch@peptide_matrix, data = batch@peptide_data, positions = positions, part = part)
    }
    consume <- function(result, index) {
        raw <- result$raw
        nonpositive <<- nonpositive + sum(raw <= 0, na.rm = TRUE)
        # The whole raw matrix's NA cells: every absent feature/sample combination plus the observed rows whose
        # quantity is NA; a batch's own matrix sees only its samples, so the count is closed at the end.
        na_rows <<- na_rows + sum(is.na(suppressWarnings(as.numeric(result$part[[quantity]]))))
        logged <- result$logged
        assembled[rownames(logged), colnames(logged)] <<- logged
        # The kernel's updated quantity column for the batch's source rows: the change map over the source.
        part <- result$part
        updated <- result$data
        key <- paste(as.character(part[[protein]]), as.character(part[[peptide]]), as.character(part[[sample]]), sep = "\r")
        updated_key <- paste(as.character(updated[[protein]]), as.character(updated[[peptide]]), as.character(updated[[sample]]), sep = "\r")
        values <- updated[[quantity]][match(key, updated_key)]
        changes <- data.frame(.dia_source_order = as.double(result$positions), stringsAsFactors = FALSE)
        changes[[quantity]] <- values
        changes <- changes[order(changes$.dia_source_order), , drop = FALSE]
        for (start in seq.int(1L, max(1L, nrow(changes)), by = 10000L)) {
            at <- if (nrow(changes)) seq.int(start, min(nrow(changes), start + 9999L)) else integer()
            offset <<- appendWorkflowTableBatch(scope, "dia_log2_data", changes[at, , drop = FALSE], offset,
                key_column = ".dia_source_order")
        }
        result <- raw <- logged <- part <- updated <- changes <- NULL
        invisible(NULL)
    }
    partition <- partitionByFeatureGroup(scope, c(protein, peptide, sample, quantity), protein, kernel, consume, worker_cap,
        working_bytes_per_row = diaPeptideGroupWorkingBytes,
        overhead_bytes = engine_bytes + as.numeric(object.size(assembled)) + 268435456 + 2 * 67108864,
        policy = "log2-feature-group-batches-v1")
    if (offset != scope$metadata$rows) cli::cli_abort("the log2 change map does not cover the source", class = "saved_entity_corrupt")
    identity <- newSavedWorkflowScientificIdentity(source$entity_id, before$data$components$design_matrix$entity_id,
        digest::digest(list("log2TransformPeptideMatrix", quantity), algo = "sha256"), "complete-declared-groups", context$software)
    map_path <- finishWorkflowTableBatch(scope, "dia_log2_data", identity)
    metadata <- scope$metadata
    streamed <- saveWorkflowTable(scope$store, NULL, basename(tempfile("dia-data-")), identity, metadata$key_column,
        metadata = metadata, write_fn = function(destination, restoration) writeWorkflowTableChanges(scope,
            map_path, NULL, restoration, destination, source_column = ".dia_source_order"))
    # The object as the whole kernel leaves it: the log2 matrix, the flag, the design cleaned to the table's samples.
    object@peptide_matrix <- assembled
    object@is_logged_data <- TRUE
    object@peptide_data <- object@peptide_data[0, , drop = FALSE]
    # `cleanDesignMatrixPeptide` reads the table's distinct samples in first-appearance order: the prelude's list.
    sample_proxy <- object
    sample_proxy@peptide_data <- stats::setNames(data.frame(samples, stringsAsFactors = FALSE), sample)
    object@design_matrix <- cleanDesignMatrixPeptide(sample_proxy)@design_matrix
    sample_proxy <- NULL
    peptides <- nrow(assembled)
    missing_before <- length(assembled) - metadata$rows + na_rows
    missing_after <- sum(is.na(assembled))
    checkpoint <- checkpointDIAReferencePeptideObject(manager, reference, object, c("design_matrix", "peptide_matrix"))
    checkpoint$components$peptide_data <- streamed
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, checkpoint$dependencies, checkpoint$components,
        checkpoint$support, checkpoint$products, partition$dependencies))
    checkpoint$dimensions <- list(rows = metadata$rows, columns = metadata$columns, samples = nrow(object@design_matrix))
    assembled <- NULL
    config <- list(peptides = peptides, samples = length(samples), quantity_column = quantity, nonpositive_set_to_na = nonpositive,
        missing_before = missing_before, missing_after = missing_after,
        outcomes = list(peptides = peptides, samples = length(samples), nonpositive_set_to_na = nonpositive,
            missing_cells = missing_after))
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data, state_name = "log2_matrix",
        description = "Peptide matrix calculated and log2 transformed (values <= 0 set to NA)",
        config = config, section = "log2TransformPeptideMatrix", counts = NULL,
        execution = c(list(variant = "streamed", worker_cap = worker_cap, rows = peptides, samples = length(samples)),
            partition[c("batches", "batch_sizes", "groups", "admission")]))
}

#' The FASTA table a peptide accession cleanup joins (as the protein cleanup resolves it: the bound path or the
#' owner's source directory), with the accession column the step names (the GUI's table carries `database_id`).
#' @noRd
resolveDIAPeptideAccessionFasta <- function(manager, parameters, seqinr, fasta_metadata) {
    if (is.null(seqinr)) {
        path <- parameters$fasta_table %||% file.path(manager$.context$experiment_paths$source_dir %||% "", "aa_seq_tbl_final.RDS")
        if (is.character(path) && length(path) == 1L && nzchar(path) && file.exists(path)) seqinr <- readRDS(path)
        metadata_path <- file.path(dirname(path), "fasta_metadata.RDS")
        if (is.null(fasta_metadata) && file.exists(metadata_path)) fasta_metadata <- readRDS(metadata_path)
    }
    column <- parameters$seqinr_accession_column %||% "uniprot_acc"
    if (!is.null(seqinr) && !column %in% names(seqinr)) {
        if (!"database_id" %in% names(seqinr))
            cli::cli_abort("the FASTA table carries neither {.val {column}} nor {.val database_id}", class = "workflow_execution_invalid")
        seqinr <- dplyr::rename(seqinr, !!column := "database_id")
    }
    list(seqinr = seqinr, fasta_metadata = fasta_metadata, column = column)
}

#' Peptide-level accession cleanup on a saved peptide state (chunk 16): the peptide `chooseBestProteinAccession`
#' (one accession per protein group from the FASTA metadata, duplicate peptides aggregated, the matrix recalculated).
#' Without a FASTA table the step records a skip, as the protein cleanup does.
#' @noRd
computeDIAReferencePeptideAccessionCleanup <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "partitioned", seqinr = NULL, fasta_metadata = NULL, now_fn = Sys.time) {
    if (!is.character(variant) || length(variant) != 1L || !variant %in% c("whole", "partitioned"))
        cli::cli_abort("unsupported peptide accession cleanup variant", class = "workflow_execution_unsupported")
    if (identical(variant, "partitioned")) return(computeDIAReferencePeptideAccessionCleanupPartitioned(manager, reference,
        parameters, worker_cap, seqinr, fasta_metadata, now_fn))
    delimiter <- parameters$delimiter
    aggregation <- parameters$aggregation_method
    fasta_sha256 <- parameters$fasta_sha256 %||% manager$.context$snapshot$file_dependencies$fasta$sha256 %||% NA_character_
    fasta <- resolveDIAPeptideAccessionFasta(manager, parameters, seqinr, fasta_metadata)
    object <- readDIAReferencePeptideObject(manager, reference)
    key <- object@protein_id_column
    groups_before <- length(unique(object@peptide_data[[key]]))
    rows_before <- nrow(object@peptide_data)
    chooseBestProteinAccessionFn <- chooseBestProteinAccession
    if (!is.null(fasta$seqinr)) {
        cleaned <- chooseBestProteinAccessionFn(theObject = object, delim = delimiter, seqinr_obj = fasta$seqinr,
            seqinr_accession_column = fasta$column, aggregation_method = aggregation)
        applied <- TRUE
        loaded <- c("peptide_data", "design_matrix", "peptide_matrix")
    } else {
        cleaned <- object
        applied <- FALSE
        loaded <- c("peptide_data", "design_matrix")
    }
    object <- NULL
    groups_after <- length(unique(cleaned@peptide_data[[key]]))
    rows_after <- nrow(cleaned@peptide_data)
    checkpoint <- checkpointDIAReferencePeptideObject(manager, reference, cleaned, loaded)
    cleaned <- NULL
    config <- list(delimiter = delimiter, aggregation_method = aggregation, seqinr_accession_column = fasta$column,
        cleanup_applied = applied, active_protein_key = key, fasta_sha256 = fasta_sha256,
        protein_groups_before = groups_before, protein_groups_after = groups_after,
        peptide_rows_before = rows_before, peptide_rows_after = rows_after,
        had_full_metadata = if (!is.null(fasta$fasta_metadata)) isTRUE(fasta$fasta_metadata$has_protein_evidence) &&
            isTRUE(fasta$fasta_metadata$has_gene_names) else FALSE,
        outcomes = list(cleanup_applied = applied, protein_groups_before = groups_before, protein_groups_after = groups_after,
            peptide_rows_aggregated = rows_before - rows_after))
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data, state_name = "peptide_accession_cleaned",
        description = if (applied) "Applied peptide-level protein accession cleanup" else
            "Skipped peptide-level accession cleanup (no FASTA data)",
        config = config, section = "chooseBestProteinAccession", counts = NULL,
        execution = list(variant = variant, worker_cap = worker_cap, rows = rows_after, cleanup_applied = applied),
        accession_cleanup_results = c(config[c("cleanup_applied", "aggregation_method", "active_protein_key",
            "had_full_metadata")], list(delimiter_used = delimiter, timestamp = now_fn())))
}

#' The peptide accession cleanup on the kit's feature-group partitioner (KIT-057, audit F-02): the feature-group unit
#' realised. The method keeps the DIA-NN group key and replaces only the `Protein.Ids` annotation, and its
#' duplicate aggregation groups by (sample, group, sequence), so no row ever leaves its protein group: batches of
#' whole groups run the unchanged method (`chooseBestProteinAccession`: the helper on the batch's groups gives the
#' same mapping as on the whole table, which depends on each group's candidate list and the FASTA scores only;
#' the aggregation per group; the matrix per batch) while the aggregated rows stream through the bounded table
#' writer and are saved in the whole method's row order (sorted by sample, group and sequence, as `summarise`
#' orders its keys), and the matrix rows are placed by name in the whole pivot's order (features in their first
#' appearance in that sorted table; samples sorted). The feature-key map and the accession provenance are the
#' kernel's own functions on the distinct identities in that order. Admission as the partitioner records it.
#' @return The whole variant's compute-only contract.
#' @noRd
computeDIAReferencePeptideAccessionCleanupPartitioned <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    seqinr = NULL, fasta_metadata = NULL, now_fn = Sys.time) {
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) || worker_cap <= 0)
        cli::cli_abort("invalid admitted worker capacity", class = "workflow_execution_unsupported")
    delimiter <- parameters$delimiter
    aggregation <- parameters$aggregation_method
    fasta_sha256 <- parameters$fasta_sha256 %||% manager$.context$snapshot$file_dependencies$fasta$sha256 %||% NA_character_
    fasta <- resolveDIAPeptideAccessionFasta(manager, parameters, seqinr, fasta_metadata)
    if (is.null(fasta$seqinr)) {
        # No FASTA table: the whole variant records a skip on the unchanged object; the same here.
        return(computeDIAReferencePeptideAccessionCleanup(manager, reference, parameters, worker_cap, variant = "whole",
            seqinr = seqinr, fasta_metadata = fasta_metadata, now_fn = now_fn))
    }
    context <- manager$.context
    before <- readDIAStateMetadata(manager, reference)
    if (!identical(before$data$class, "PeptideQuantitativeData"))
        cli::cli_abort("a peptide state is required", class = "workflow_execution_unsupported")
    object <- readDIAReferencePeptideObject(manager, reference, fields = "design_matrix")
    source <- before$data$components$peptide_data
    scope <- newWorkflowTableReadScope(context$store, source)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    con <- scope$connection
    engine_bytes <- min(1e9, floor(worker_cap / 8))
    DBI::dbExecute(con, paste0("SET memory_limit = '", floor(engine_bytes / 1048576), "MB'"))
    DBI::dbExecute(con, "SET threads = 2")
    if (any(c(".dia_source_order", ".dia_design_order", ".dia_feature_order") %in% scope$metadata$columns))
        cli::cli_abort("reserved DIA grouping columns are present", class = "saved_entity_invalid")
    key <- object@protein_id_column
    peptide <- object@peptide_sequence_column
    sample <- object@sample_id
    accession <- "Protein.Ids"
    columns <- scope$metadata$columns
    if (!all(c(key, peptide, sample, accession, object@raw_quantity_column, object@norm_quantity_column) %in% columns))
        cli::cli_abort("the peptide table lacks the accession cleanup's columns", class = "saved_entity_invalid")
    # Prelude (one pass over the identity columns): the distinct features with the first sample each appears in, the
    # distinct group/candidate pairs (the helper's input) and the samples.
    # Per feature (group, sequence) only its smallest sample in the C-locale order is kept: each page is reduced to
    # one row per feature before it meets the running table (one row per feature, never the table's rows).
    features <- NULL
    candidates <- NULL
    samples <- character()
    reduce <- function(f) { f <- f[order(f$sample, method = "radix"), , drop = FALSE]; f[!duplicated(f[c("group", "sequence")]), , drop = FALSE] }
    streamWorkflowTableRows(scope, c(key, peptide, sample, accession), function(page, offset, row_names) {
        if (!nrow(page)) return(invisible(NULL))
        f <- reduce(data.frame(group = as.character(page[[key]]), sequence = as.character(page[[peptide]]),
            sample = as.character(page[[sample]]), stringsAsFactors = FALSE))
        features <<- if (is.null(features)) f else reduce(rbind(features, f))
        c2 <- data.frame(group = as.character(page[[key]]), candidates = as.character(page[[accession]]), stringsAsFactors = FALSE)
        c2 <- c2[!duplicated(c2), , drop = FALSE]
        candidates <<- if (is.null(candidates)) c2 else { joined <- rbind(candidates, c2); joined[!duplicated(joined), , drop = FALSE] }
        samples <<- unique(c(samples, as.character(page[[sample]])))
        invisible(NULL)
    })
    if (is.null(features) || !nrow(features)) cli::cli_abort("the peptide table is empty", class = "saved_entity_invalid")
    rows_before <- scope$metadata$rows
    groups_before <- length(unique(features$group))
    # The whole method's output order: `summarise` sorts its keys (sample, group, sequence) in the C locale.
    sample_order <- stats::setNames(seq_along(samples), samples[order(samples, method = "radix")])
    first_sample <- data.frame(group = features$group, sequence = features$sequence, sample_rank = unname(sample_order[features$sample]),
        stringsAsFactors = FALSE)
    # The matrix's rows: each feature's first appearance in the sorted table (sample block, then group, sequence).
    first_sample <- first_sample[order(first_sample$sample_rank, first_sample$group, first_sample$sequence, method = "radix"), , drop = FALSE]
    # The table's rows: the plain key order (sample, group, sequence); the feature rank is the (group, sequence) rank.
    key_rank <- order(first_sample$group, first_sample$sequence, method = "radix")
    feature_order <- stats::setNames(integer(nrow(first_sample)), paste(first_sample$group, first_sample$sequence, sep = "\r"))
    feature_order[key_rank] <- seq_along(key_rank)
    features <- NULL
    # The mapping once, on the distinct group/candidate pairs, as the method computes it row by row.
    mapping_input <- stats::setNames(candidates, c(key, accession))
    mapping <- chooseBestProteinAccessionHelper(input_tbl = mapping_input, acc_detail_tab = fasta$seqinr,
        accessions_column = !!rlang::sym(accession), row_id_column = fasta$column, group_id = !!rlang::sym(key), delim = delimiter)
    mapped <- stats::setNames(as.character(mapping[[fasta$column]]), as.character(mapping[[key]]))
    # The kernel's feature-key map and provenance on the distinct identities in the output's order: the identities
    # carry the new annotation, as the aggregated table does.
    # A tibble, as the kernel's table is: its provenance (the distinct group/accession pairs) then carries the whole
    # variant's class and no row names (the screen oracle comparison of 2026-10-24 found both differing).
    identities <- tibble::tibble(first_sample$group, first_sample$sequence, .name_repair = "minimal")
    names(identities) <- c(key, peptide)
    identities[[accession]] <- unname(mapped[identities[[key]]])
    proxy <- object
    proxy@peptide_data <- identities
    proxy <- .ensurePeptideFeatureKeyMap(proxy)
    feature_map <- proxy@args$peptide_feature_key_map
    object@args <- proxy@args
    sorted_samples <- names(sample_order)
    assembled <- matrix(NA_real_, nrow = nrow(feature_map), ncol = length(sorted_samples),
        dimnames = list(feature_map$feature_key, sorted_samples))
    offset <- 0L
    restoration <- NULL
    method_args <- NULL
    kernel <- function(part, index) {
        attr(part, "workflow_source_order") <- NULL
        batch <- object
        batch@peptide_data <- part
        # Called through the same alias as the whole variant: the method records its parameters under the caller's name.
        chooseBestProteinAccessionFn <- chooseBestProteinAccession
        cleaned <- chooseBestProteinAccessionFn(theObject = batch, delim = delimiter, seqinr_obj = fasta$seqinr,
            seqinr_accession_column = fasta$column, aggregation_method = aggregation)
        list(data = cleaned@peptide_data, matrix = cleaned@peptide_matrix, args = cleaned@args)
    }
    consume <- function(result, index) {
        m <- result$matrix
        assembled[rownames(m), colnames(m)] <<- m
        # The method's recorded parameters (the same every batch): carried onto the object as the whole variant's are.
        method_args <<- result$args[["chooseBestProteinAccessionFn"]]
        table <- result$data # the method's tibble, kept as it is: the saved table's class is part of its metadata
        table$.dia_design_order <- unname(sample_order[as.character(table[[sample]])])
        table$.dia_feature_order <- unname(feature_order[paste(as.character(table[[key]]), as.character(table[[peptide]]), sep = "\r")])
        if (anyNA(table$.dia_design_order) || anyNA(table$.dia_feature_order))
            cli::cli_abort("an aggregated row has no place in the output order", class = "saved_entity_corrupt")
        if (is.null(restoration)) restoration <<- workflowSavedTableMetadata(table, ".rownames", 3L)
        for (start in seq.int(1L, max(1L, nrow(table)), by = 10000L)) {
            at <- if (nrow(table)) seq.int(start, min(nrow(table), start + 9999L)) else integer()
            offset <<- appendWorkflowTableBatch(scope, "dia_accession_data", table[at, , drop = FALSE], offset, key_column = ".rownames")
        }
        result <- m <- table <- NULL
        invisible(NULL)
    }
    partition <- partitionByFeatureGroup(scope, columns, key, kernel, consume, worker_cap,
        working_bytes_per_row = diaPeptideGroupWorkingBytes,
        overhead_bytes = engine_bytes + as.numeric(object.size(assembled)) + 268435456 + 2 * 67108864,
        policy = "accession-feature-group-batches-v1")
    rows_after <- offset
    identity <- newSavedWorkflowScientificIdentity(source$entity_id, before$data$components$design_matrix$entity_id,
        digest::digest(list("chooseBestProteinAccession", delimiter, aggregation, fasta$column, fasta_sha256), algo = "sha256"),
        "complete-declared-groups", context$software)
    path <- finishWorkflowTableBatch(scope, "dia_accession_data", identity)
    internal <- c(".dia_design_order", ".dia_feature_order")
    metadata <- restoration
    metadata$rows <- rows_after
    output_columns <- setdiff(metadata$columns, internal)
    metadata$columns <- output_columns
    metadata$column_metadata <- metadata$column_metadata[output_columns]
    metadata$name_columns <- metadata$name_columns[intersect(names(metadata$name_columns), output_columns)]
    metadata$schema <- 3L
    metadata$codec <- "parquet_duckdb_v3"
    order_columns <- c(internal, ".multischolar_order")
    source_scope <- newWorkflowTableReadScope(context$store, path, reader = "arrow", metadata_only = TRUE)
    ordering <- saveWorkflowOrderIndex(source_scope, order_columns)
    releaseWorkflowTableReadScope(source_scope)
    streamed <- saveWorkflowTable(context$store, NULL, basename(tempfile("dia-data-")), identity, metadata$key_column,
        metadata = metadata, write_fn = function(destination, restoration) {
            cursor <- openWorkflowOrderedRows(context$store, path, output_columns, order_columns)
            on.exit(cursor$close(), add = TRUE)
            writer <- openWorkflowTableWriter(destination, restoration, context$store$progress)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            repeat {
                data <- cursor$next_batch()
                if (is.null(data)) break
                appendWorkflowTableWriter(writer, completeWorkflowTableBatch(data, restoration))
            }
            if (is.null(writer$output)) appendWorkflowTableWriter(writer,
                completeWorkflowTableBatch(readSavedWorkflowTablePrototype(scope), restoration))
            finishWorkflowTableWriter(writer)
        })
    object@peptide_matrix <- assembled
    object@peptide_data <- object@peptide_data[0, , drop = FALSE]
    if (!is.null(method_args)) object@args[["chooseBestProteinAccessionFn"]] <- method_args
    checkpoint <- checkpointDIAReferencePeptideObject(manager, reference, object, c("design_matrix", "peptide_matrix"))
    checkpoint$components$peptide_data <- streamed
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, checkpoint$dependencies, checkpoint$components,
        checkpoint$support, checkpoint$products, partition$dependencies, ordering$reference))
    checkpoint$dimensions <- list(rows = rows_after, columns = output_columns, samples = nrow(object@design_matrix))
    assembled <- NULL
    groups_after <- length(unique(first_sample$group))
    config <- list(delimiter = delimiter, aggregation_method = aggregation, seqinr_accession_column = fasta$column,
        cleanup_applied = TRUE, active_protein_key = key, fasta_sha256 = fasta_sha256,
        protein_groups_before = groups_before, protein_groups_after = groups_after,
        peptide_rows_before = rows_before, peptide_rows_after = rows_after,
        had_full_metadata = if (!is.null(fasta$fasta_metadata)) isTRUE(fasta$fasta_metadata$has_protein_evidence) &&
            isTRUE(fasta$fasta_metadata$has_gene_names) else FALSE,
        outcomes = list(cleanup_applied = TRUE, protein_groups_before = groups_before, protein_groups_after = groups_after,
            peptide_rows_aggregated = rows_before - rows_after))
    list(checkpoint = checkpoint, reference = checkpoint$components$peptide_data, state_name = "peptide_accession_cleaned",
        description = "Applied peptide-level protein accession cleanup",
        config = config, section = "chooseBestProteinAccession", counts = NULL,
        execution = c(list(variant = "partitioned", worker_cap = worker_cap, rows = rows_after, cleanup_applied = TRUE),
            partition[c("batches", "batch_sizes", "groups", "admission")]),
        accession_cleanup_results = c(config[c("cleanup_applied", "aggregation_method", "active_protein_key",
            "had_full_metadata")], list(delimiter_used = delimiter, timestamp = now_fn())))
}

#' DPC imputation and protein quantification on a saved log2 peptide state (chunk 17): the peptide method of
#' `proteinMissingValueImputationLimpa` with the declared engine; the results list it attaches to the protein
#' object's `args` is stored as the class's owned entry (`args_as_owned`).
#' @param parameters `engine` (`fastlimpadpc` or `limpa`), `dpc_slope`, `chunk`.
#' @return Compute-only protein contract.
#' @noRd
#' KIT-057: one streamed pass over the peptide table's identity columns for the bounded DPC worker: the distinct
#' (protein, sequence) identities in first-appearance order with `Protein.Ids` (the kernel's feature-key map and
#' accession provenance are rebuilt from them exactly as from the whole table) and the kernel's per-protein counts
#' (`peptide_count = n_distinct(sequence)`, `peptidoform_count` = the rows per protein), computed by the kernel's own
#' grouping over one row per distinct pair carrying its row count.
#' @return List `identities` (data frame), `summary` (the kernel's `peptide_summary`), `rows` (table rows).
#' @noRd
streamDIAPeptideIdentitySummary <- function(manager, reference, object, worker_cap = 12 * 1024^3) {
    before <- readDIAStateMetadata(manager, reference)
    source <- before$data$components$peptide_data
    scope <- newWorkflowTableReadScope(manager$.context$store, source)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    con <- scope$connection
    engine_bytes <- min(1e9, floor(worker_cap / 8))
    DBI::dbExecute(con, paste0("SET memory_limit = '", floor(engine_bytes / 1048576), "MB'"))
    DBI::dbExecute(con, "SET threads = 2")
    protein <- object@protein_id_column
    peptide <- object@peptide_sequence_column
    if (!all(c(protein, peptide) %in% scope$metadata$columns))
        cli::cli_abort("the peptide table lacks the DPC step's identity columns", class = "saved_entity_invalid")
    identity_columns <- intersect(unique(c(protein, peptide, "Protein.Ids")), scope$metadata$columns)
    identities <- NULL
    proteins <- character()
    sequences <- character()
    counts <- integer()
    rows <- 0L
    streamWorkflowTableRows(scope, identity_columns, function(page, offset, row_names) {
        if (!nrow(page)) return(invisible(NULL))
        rows <<- rows + nrow(page)
        page_protein <- as.character(page[[protein]])
        page_sequence <- as.character(page[[peptide]])
        key <- paste(page_protein, page_sequence, sep = "\r")
        first <- !duplicated(key)
        n <- tabulate(match(key, key[first]), sum(first))
        index <- match(key[first], paste(proteins, sequences, sep = "\r"))
        new <- is.na(index)
        counts[index[!new]] <<- counts[index[!new]] + n[!new]
        if (any(new)) {
            proteins <<- c(proteins, page_protein[first][new])
            sequences <<- c(sequences, page_sequence[first][new])
            counts <<- c(counts, n[new])
        }
        rows_page <- page[identity_columns]
        for (name in identity_columns) rows_page[[name]] <- as.character(rows_page[[name]])
        rows_page <- rows_page[!duplicated(rows_page[c(protein, peptide)]), , drop = FALSE]
        identities <<- if (is.null(identities)) rows_page else {
            joined <- rbind(identities, rows_page)
            joined[!duplicated(joined[c(protein, peptide)]), , drop = FALSE]
        }
        invisible(NULL)
    })
    if (is.null(identities) || !nrow(identities)) cli::cli_abort("the peptide table is empty", class = "saved_entity_invalid")
    rownames(identities) <- NULL
    pairs <- stats::setNames(list(proteins, sequences, counts), c(protein, peptide, "n_rows"))
    pairs <- tibble::as_tibble(pairs)
    summary <- pairs |>
        dplyr::group_by(!!rlang::sym(protein)) |>
        dplyr::summarise(peptide_count = dplyr::n_distinct(!!rlang::sym(peptide)), peptidoform_count = sum(n_rows), .groups = "drop")
    list(identities = identities, summary = summary, rows = rows)
}

computeDIAReferenceProteinDpcQuant <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "bounded", now_fn = Sys.time) {
    if (!variant %in% c("whole", "bounded")) cli::cli_abort("unsupported DPC quantification variant", class = "workflow_execution_unsupported")
    engine <- resolveLimpaDpcEngine(parameters$engine)
    dpc_slope <- parameters$dpc_slope %||% 0.8
    chunk <- parameters$chunk %||% 1000
    # KIT-057 `bounded` (user request 2026-10-25): the worker hydrates the matrix and the design only. The kernel reads
    # the long table for two things, the feature-key map and provenance (rebuilt by `.ensurePeptideFeatureKeyMap` from
    # the distinct identities, as the streamed log2 and partitioned accession preludes do) and the per-protein counts
    # (`n_distinct(sequence)` and the rows per protein), which a streamed pass over the identity columns computes and
    # hands to the kernel as `peptide_summary`. `whole` hydrates the table and lets the kernel compute both.
    fields <- if (identical(variant, "bounded")) c("design_matrix", "peptide_matrix") else c("peptide_data", "design_matrix", "peptide_matrix")
    object <- readDIAReferencePeptideObject(manager, reference, fields = fields)
    if (!length(object@peptide_matrix))
        cli::cli_abort("DPC quantification needs the log2 peptide matrix state", class = "workflow_execution_invalid")
    if (!isTRUE(object@is_logged_data))
        cli::cli_abort("DPC quantification needs a log2 transformed peptide state", class = "workflow_execution_invalid")
    peptides <- nrow(object@peptide_matrix)
    peptide_summary <- NULL
    execution <- list(variant = variant, worker_cap = worker_cap)
    if (identical(variant, "bounded")) {
        prelude <- streamDIAPeptideIdentitySummary(manager, reference, object, worker_cap)
        object@peptide_data <- prelude$identities
        peptide_summary <- prelude$summary
        execution <- c(execution, list(table_rows = prelude$rows, identities = nrow(prelude$identities),
            proteins_counted = nrow(prelude$summary)))
        prelude <- NULL
    }
    proteinMissingValueImputationLimpaFn <- proteinMissingValueImputationLimpa
    quantified <- proteinMissingValueImputationLimpaFn(object, dpc_slope = dpc_slope, verbose = FALSE, chunk = chunk,
        engine = engine, peptide_summary = peptide_summary)
    object <- NULL
    peptide_summary <- NULL
    results <- quantified@args$limpa_dpc_quant_results
    engine_version <- as.character(utils::packageVersion(engine))
    quantified@args$limpaDpcQuant <- list(engine = engine, engine_version = engine_version, dpc_slope = dpc_slope,
        chunk = chunk, limpa_version = as.character(utils::packageVersion("limpa")))
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, quantified)
    protein_count <- length(unique(quantified@protein_quant_table[[quantified@protein_id_column]]))
    samples <- nrow(quantified@design_matrix)
    key <- quantified@protein_id_column
    quantified <- NULL
    dpc <- unname(as.numeric(results$dpc_parameters_used))
    config <- list(engine = engine, engine_version = engine_version, dpc_slope = dpc_slope, chunk = chunk,
        dpc_intercept = dpc[[1L]], dpc_slope_realised = dpc[[2L]], slope_interpretation = results$slope_interpretation,
        proteins_quantified = results$total_proteins_quantified, peptides_used = results$total_peptides_used,
        missing_percentage_peptides = results$missing_percentage_peptides,
        missing_percentage_proteins = results$missing_percentage_proteins, protein_id_column = key,
        outcomes = list(engine = engine, engine_version = engine_version, limpa_version = as.character(utils::packageVersion("limpa")), dpc_intercept = dpc[[1L]],
            dpc_slope = dpc[[2L]], slope_interpretation = results$slope_interpretation,
            proteins_quantified = results$total_proteins_quantified, peptides_used = results$total_peptides_used,
            missing_percentage_peptides = results$missing_percentage_peptides))
    results <- NULL
    text <- paste(
        "Limpa DPC Quantification Completed Successfully\n",
        "===============================================\n",
        sprintf("Proteins quantified: %d\n", protein_count),
        sprintf("Samples: %d\n", samples),
        sprintf("Peptides used: %d\n", peptides),
        sprintf("Engine: %s %s\n", engine, engine_version),
        sprintf("DPC intercept / slope: %.4f / %.4f (%s)\n", config$dpc_intercept, config$dpc_slope_realised, config$slope_interpretation),
        sprintf("Missing values in peptides: %.1f%%; in proteins: %.1f%%\n", config$missing_percentage_peptides, config$missing_percentage_proteins),
        sprintf("Active protein key: %s\n", key),
        "S4 Class: ProteinQuantitativeData\n",
        "State saved as: 'protein_dpc_quant'\n")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "protein_dpc_quant",
        description = "Limpa DPC imputation and protein quantification completed (ProteinQuantitativeData created)",
        config = config, section = "limpaDpcQuant", protein_count = protein_count,
        execution = c(execution, list(rows = protein_count, samples = samples, peptides = peptides)),
        text = text, field = "quantifiedS4", qc_name = "dpc_quant",
        qc = c(config[setdiff(names(config), "outcomes")], list(timestamp = now_fn())))
}

#' The in-process entry points of the three steps (the GUI without a managed application, the kit's checks): the
#' registered declaration, the provenance it names, the shared finishing step.
#' @noRd
runDIAReferencePeptideLog2Matrix <- function(workflow) {
    manager <- workflow$state_manager
    before <- readDIAStateMetadata(manager)
    pending <- if (identical(manager$current_state, "log2_matrix")) resumeDIAQvalueOutput(manager, before) else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    entry <- workflowActionSpec("peptide_log2_matrix")
    value <- computeDIAReferencePeptideLog2Matrix(manager, getDIAStateReference(manager), list(), variant = entry$default_variant)
    setDIAStepProvenance(manager, entry$name, entry$default_variant, entry$parameters, entry$declaration$fixed)
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
        value$section, value$counts)
}

#' @noRd
runDIAReferencePeptideAccessionCleanup <- function(workflow, parameters) {
    manager <- workflow$state_manager
    before <- readDIAStateMetadata(manager)
    pending <- if (identical(manager$current_state, "peptide_accession_cleaned")) resumeDIAQvalueOutput(manager, before) else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    entry <- workflowActionSpec("peptide_accession_cleanup")
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    value <- computeDIAReferencePeptideAccessionCleanup(manager, getDIAStateReference(manager), parameters,
        variant = entry$default_variant, seqinr = shiny::isolate(workflow$aa_seq_tbl_final),
        fasta_metadata = shiny::isolate(workflow$fasta_metadata))
    setDIAStepProvenance(manager, entry$name, entry$default_variant, entry$parameters, entry$declaration$fixed,
        values = c(parameters, list(fasta_sha256 = value$config$fasta_sha256, seqinr_accession_column = value$config$seqinr_accession_column)))
    workflow$accession_cleanup_results <- value$accession_cleanup_results
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
        value$section, value$counts)
}

#' @noRd
runDIAReferenceProteinDpcQuant <- function(workflow, parameters) {
    manager <- workflow$state_manager
    pending <- if (identical(manager$current_state, "protein_dpc_quant")) resumeDIAProteinAction(workflow, "protein_dpc_quant") else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    entry <- workflowActionSpec("protein_dpc_quant")
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    value <- computeDIAReferenceProteinDpcQuant(manager, getDIAStateReference(manager), parameters, variant = entry$default_variant)
    setDIAStepProvenance(manager, entry$name, entry$default_variant, entry$parameters, entry$declaration$fixed,
        values = value$config[c("engine", "dpc_slope", "chunk")])
    finishDIAReferenceProteinStep(workflow, value)
}

# ---------------------------------------------------------------------------------------------------------------
# The DIA_LIMPA protein stage (KIT-041): the limpa script's chunks 17b and 18 as declarations on the owner (the
# correlation steps are the shared `pearson-pairs`/`pearson-selection`; the RUV-III factors as covariates are the
# in-process `protein_ruv_covariates` declaration in utils_workflow_dia_products.R).

#' The between-sample normalisation of the quantified protein matrix as a declaration (chunk 17b): whole matrix
#' (every normaliser other than `none` uses all columns), the unchanged protein `normaliseBetweenSamples` once;
#' `none` (the script's choice and the shipped default) commits the kernel's unchanged output with `skipped = TRUE`.
#' @noRd
diaProteinBetweenSampleNormalisationOperation <- function() {
    newWorkflowOperation(name = "protein_between_sample_normalisation", input_class = "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        kernel = function(object, parameters, store) computeDIAReferenceProteinNormalisation(store$manager, store$source,
            parameters["norm_method"], worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData", parameters = "norm_method",
        product_type = "dia.protein_between_sample_normalisation_result")
}

#' The minimum-peptide-evidence filter on the quantified protein object as a declaration (chunk 18): unit row
#' (each protein kept by its own peptide and peptidoform counts), the unchanged protein
#' `filterMinNumPeptidesPerProtein` once, placed whole.
#' @noRd
diaProteinPeptideCountFilterOperation <- function() {
    newWorkflowOperation(name = "protein_peptide_count_filter", input_class = "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "row",
        kernel = function(object, parameters, store) computeDIAReferenceProteinPeptideCountFilter(store$manager, store$source,
            parameters[c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, variant = parameters$variant %||% "whole"),
        output_class = "ProteinQuantitativeData", parameters = c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"),
        product_type = "dia.protein_peptide_count_filter_result",
        merge = c(protein_quant_table = "concatenate", design_matrix = "design", protein_id_table = "concatenate"))
}

#' @noRd
computeDIAReferenceProteinNormalisation <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "whole", now_fn = Sys.time) {
    if (!identical(variant, "whole")) cli::cli_abort("unsupported protein normalisation variant", class = "workflow_execution_unsupported")
    method <- parameters$norm_method %||% "none"
    object <- readDIAReferenceProteinObject(manager, reference)
    proteins_before <- nrow(object@protein_quant_table)
    # The kernel records its arguments under the calling name (the DIA normalisation declaration calls it through
    # the same alias).
    normaliseBetweenSamplesFn <- normaliseBetweenSamples
    normalised <- normaliseBetweenSamplesFn(object, normalisation_method = method)
    object <- NULL
    protein_count <- nrow(normalised@protein_quant_table)
    samples <- nrow(normalised@design_matrix)
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, normalised)
    normalised <- NULL
    skipped <- identical(method, "none")
    config <- list(norm_method = method, normalisation_method = method, skipped = skipped,
        proteins = protein_count, samples = samples,
        outcomes = list(norm_method = method, proteins = protein_count, samples = samples))
    text <- paste(
        if (skipped) "Between-Sample Normalisation Skipped (none)\n" else "Between-Sample Normalisation Applied Successfully\n",
        "================================================\n",
        sprintf("Proteins remaining: %d\n", protein_count), sprintf("Samples: %d\n", samples),
        sprintf("Method: %s\n", method),
        if (skipped) "The quantified log2 matrix is carried unchanged (the limpa script's choice).\n" else "",
        "State saved as: 'protein_normalised'\n")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "protein_normalised",
        description = if (skipped) "Between-sample normalisation skipped (none): the quantified matrix carried unchanged" else
            sprintf("Between-sample normalisation applied (%s)", method),
        config = config, section = "limpaProteinNormalisation", protein_count = protein_count,
        execution = list(variant = variant, worker_cap = worker_cap, rows = protein_count, samples = samples, skipped = skipped),
        text = text, field = "normalisedS4", qc_name = "protein_normalisation",
        qc = c(config[setdiff(names(config), "outcomes")], list(timestamp = now_fn())))
}

#' @noRd
computeDIAReferenceProteinPeptideCountFilter <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    variant = "whole", now_fn = Sys.time) {
    if (!identical(variant, "whole")) cli::cli_abort("unsupported peptide count filter variant", class = "workflow_execution_unsupported")
    peptides <- parameters$num_peptides_per_protein_thresh
    forms <- parameters$num_peptidoforms_per_protein_thresh
    object <- readDIAReferenceProteinObject(manager, reference)
    if (is.null(object@args$limpa_dpc_quant_results$peptide_counts_per_protein))
        cli::cli_abort("the peptide count filter needs the DPC quantification's peptide counts", class = "workflow_execution_invalid")
    proteins_before <- nrow(object@protein_quant_table)
    filterMinNumPeptidesPerProteinFn <- filterMinNumPeptidesPerProtein
    filtered <- suppressMessages(filterMinNumPeptidesPerProteinFn(object,
        num_peptides_per_protein_thresh = peptides, num_peptidoforms_per_protein_thresh = forms, verbose = FALSE))
    object <- NULL
    protein_count <- nrow(filtered@protein_quant_table)
    elist_rows <- nrow(filtered@args$limpa_dpc_quant_results$quantified_elist)
    checkpoint <- checkpointDIAReferenceProteinObject(manager, reference, filtered)
    filtered <- NULL
    config <- list(num_peptides_per_protein_thresh = peptides, num_peptidoforms_per_protein_thresh = forms,
        peptides_per_protein_cutoff = peptides, peptidoforms_per_protein_cutoff = forms,
        proteins_before = proteins_before, proteins_after = protein_count, quantified_elist_rows = elist_rows,
        outcomes = list(proteins_before = proteins_before, proteins_after = protein_count,
            proteins_removed = proteins_before - protein_count))
    text <- paste(
        "Peptide Count Filter Applied Successfully\n",
        "=========================================\n",
        sprintf("Proteins remaining: %d (removed %d)\n", protein_count, proteins_before - protein_count),
        sprintf("Min unique peptides per protein: %d\n", as.integer(peptides)),
        sprintf("Min peptidoforms per protein: %d\n", as.integer(forms)),
        sprintf("Quantified EList rows: %s\n", if (is.null(elist_rows)) "n/a" else elist_rows),
        "State saved as: 'peptide_count_filtered'\n")
    list(checkpoint = checkpoint, reference = checkpoint$components$protein_quant_table,
        state_name = "peptide_count_filtered",
        description = "Applied minimum peptide evidence filter on the quantified protein object",
        config = config, section = "filterMinNumPeptidesPerProtein", protein_count = protein_count,
        execution = list(variant = variant, worker_cap = worker_cap, rows = protein_count, removed = proteins_before - protein_count),
        text = text, field = "filteredS4", qc_name = "peptide_count_filter",
        qc = c(config[setdiff(names(config), "outcomes")], list(timestamp = now_fn())))
}

#' @noRd
runDIAReferenceProteinNormalisation <- function(workflow, parameters) {
    manager <- workflow$state_manager
    pending <- if (identical(manager$current_state, "protein_normalised")) resumeDIAProteinAction(workflow, "protein_normalised") else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    entry <- workflowActionSpec("protein_between_sample_normalisation")
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    value <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), parameters, variant = entry$default_variant)
    setDIAStepProvenance(manager, entry$name, entry$default_variant, entry$parameters, entry$declaration$fixed,
        values = value$config["norm_method"])
    finishDIAReferenceProteinStep(workflow, value)
}

#' @noRd
runDIAReferenceProteinPeptideCountFilter <- function(workflow, parameters) {
    manager <- workflow$state_manager
    pending <- if (identical(manager$current_state, "peptide_count_filtered")) resumeDIAProteinAction(workflow, "peptide_count_filtered") else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    entry <- workflowActionSpec("protein_peptide_count_filter")
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    value <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager), parameters, variant = entry$default_variant)
    setDIAStepProvenance(manager, entry$name, entry$default_variant, entry$parameters, entry$declaration$fixed,
        values = value$config[c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh")])
    finishDIAReferenceProteinStep(workflow, value)
}
