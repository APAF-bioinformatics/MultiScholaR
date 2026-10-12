#' Save the scientific columns required by the bounded DIA rollup prototype.
#' @param store Existing request-owned entity store.
#' @param object Current peptide S4; never retained by the returned specification.
#' @param request_id,identity Existing durable request and scientific identity.
#' @return A table reference and scalar column roles, not an S4 or assay.
#' @noRd
saveWorkflowPrecursorRollupInput <- function(store, object, request_id, identity) {
    if (!methods::is(object, "PeptideQuantitativeData") ||
        isTRUE(methods::slot(object, "is_logged_data"))) {
        cli::cli_abort("saved rollup requires linear peptide data", class = "saved_entity_unsupported")
    }
    roles <- list(sample = methods::slot(object, "sample_id"),
        protein = methods::slot(object, "protein_id_column"),
        peptide = methods::slot(object, "peptide_sequence_column"),
        modified = "Modified.Sequence", precursor = "Precursor.Id",
        raw = methods::slot(object, "raw_quantity_column"),
        normalised = methods::slot(object, "norm_quantity_column"))
    table <- methods::slot(object, "peptide_data")
    if (any(!unlist(roles) %in% names(table)) || ".qc_source_row" %in% names(table)) {
        cli::cli_abort("saved rollup requires declared precursor identities and column roles",
            class = "saved_entity_unsupported")
    }
    evidence <- intersect(c("identification_peptide_count", "identification_peptidoform_count",
        "peptides_for_protein_count", "peptidoforms_for_protein_count"), names(table))
    annotation <- if (roles$protein != "Protein.Ids" && "Protein.Ids" %in% names(table))
        "Protein.Ids" else character()
    columns <- unique(c(unlist(roles), evidence, annotation))
    table <- table[, columns, drop = FALSE]
    table$.qc_source_row <- seq_len(nrow(table))
    reference <- saveWorkflowTable(store, table, request_id, identity,
        key_column = ".qc_source_row", metadata_version = 2L)
    list(reference = reference, roles = roles, evidence = evidence, annotation = annotation)
}

#' Fetch one complete scientific result in bounded batches, not a display lookup.
#' @param connection Active action-owned DuckDB connection.
#' @param sql Fixed internal query with quoted identifiers.
#' @param max_rows Input row bound; aggregation cannot expand this population.
#' @return Materialized scientific data frame; total residency is measured by the service.
#' @noRd
collectWorkflowRollupQuery <- function(connection, sql, max_rows) {
    # The action owns this connection. Materialize the query once to know its
    # exact result size, rather than retaining all batches and rbind's copy.
    DBI::dbExecute(connection, paste0("CREATE TEMP TABLE qc_collection AS ", sql))
    on.exit(DBI::dbExecute(connection, "DROP TABLE IF EXISTS qc_collection"), add = TRUE)
    size <- as.numeric(DBI::dbGetQuery(connection,
        "SELECT count(*) AS n FROM qc_collection")$n)
    if (length(size) != 1L || !is.finite(size) || size < 0 || size > max_rows) {
        cli::cli_abort("rollup expanded the source population", class = "saved_entity_limit")
    }
    result <- DBI::dbSendQuery(connection, "SELECT * FROM qc_collection")
    on.exit(DBI::dbClearResult(result), add = TRUE, after = FALSE)
    prototype <- DBI::dbFetch(result, n = 0L)
    if (ncol(prototype) > 64L) {
        cli::cli_abort("rollup fetch exceeds collection limits", class = "saved_entity_limit")
    }
    columns <- lapply(prototype, \(column) {
        value <- vector(typeof(column), size)
        attributes(value) <- attributes(column)
        value
    })
    rows <- 0
    repeat {
        part <- DBI::dbFetch(result, n = 10000L)
        if (ncol(part) > 64L || as.numeric(object.size(part)) > 67108864) {
            cli::cli_abort("rollup fetch exceeds collection limits", class = "saved_entity_limit")
        }
        if (!nrow(part)) break
        next_rows <- rows + nrow(part)
        if (next_rows > size || !identical(names(part), names(prototype))) {
            cli::cli_abort("rollup expanded the source population", class = "saved_entity_limit")
        }
        at <- seq.int(rows + 1, next_rows)
        for (column in seq_along(columns)) columns[[column]][at] <- part[[column]]
        rows <- next_rows
    }
    if (rows != size) cli::cli_abort("rollup collection is incomplete", class = "saved_entity_corrupt")
    attributes(columns) <- attributes(prototype)
    attr(columns, "row.names") <- .set_row_names(size)
    columns
}

#' Calculate one full-population precursor rollup before checkpoint hydration.
#' @param store Existing owner; metadata and integrity are validated once per action.
#' @param input Reference/roles from saveWorkflowPrecursorRollupInput.
#' @param non_exact Explicitly select the diagnostic DuckDB sum() arm, whose
#'   summation order is unspecified and therefore not bit-identical to the kernel.
#' @return Native-shaped data/summary after all database resources have closed.
#' @noRd
computeSavedWorkflowPrecursorRollup <- function(store, input, persist = NULL, non_exact = FALSE) {
    profile <- getOption("multischolar.dia_rollup_profile")
    if (!is.null(persist) && !is.null(profile)) return(computeDIARelationalRollup(store, input, persist,
        arm = "F", profile = profile, phase = getOption("multischolar.dia_rollup_phase", function(label, fn) fn())))
    if (!is.null(persist)) return(persistIndexedWorkflowPrecursorRollup(store, input, persist))
    if (!isTRUE(non_exact)) cli::cli_abort(paste("the DuckDB sum() rollup arm is not order-exact;",
        "select it only with non_exact = TRUE"), class = "dia_rollup_non_exact")
    scope <- newWorkflowTableReadScope(store, input$reference)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    assertWorkflowPhysicalTableScope(scope)
    connection <- scope$connection
    metadata <- scope$metadata
    roles <- input$roles
    fields <- unique(c(unlist(roles), input$evidence, input$annotation))
    if (any(!fields %in% metadata$columns) ||
        any(vapply(metadata$column_metadata[c(roles$raw, roles$normalised)],
            function(column) !identical(column$typeof, "double"), logical(1)))) {
        cli::cli_abort("saved rollup requires the declared double-precision quantities",
            class = "saved_entity_unsupported")
    }
    DBI::dbExecute(connection, "SET threads = 2")
    DBI::dbExecute(connection, "SET memory_limit = ?", params = list("1GB"))
    quote_name <- function(value) as.character(DBI::dbQuoteIdentifier(connection, value))
    DBI::dbExecute(connection, paste0("CREATE TEMP VIEW qc_input AS SELECT ",
        paste(quote_name(fields), collapse = ", "), " FROM read_parquet(",
        DBI::dbQuoteString(connection, scope$path), ")"))
    keys <- c(roles$sample, roles$protein, roles$peptide, roles$modified, roles$precursor)
    bad_keys <- vapply(keys, function(name) {
        quoted <- quote_name(name)
        blank <- if (identical(metadata$column_metadata[[name]]$typeof, "character"))
            paste0(" OR regexp_full_match(", quoted, ", '[ \\t\\r\\n]*')") else ""
        paste0("(", quoted, " IS NULL", blank, ")")
    }, character(1))
    if (nrow(DBI::dbGetQuery(connection, paste0("SELECT 1 FROM qc_input WHERE ",
        paste(bad_keys, collapse = " OR "), " LIMIT 1")))) {
        cli::cli_abort("saved rollup contains missing or blank identity keys")
    }
    quantities <- quote_name(c(roles$raw, roles$normalised))
    bad_values <- paste0("(", quantities, " IS NOT NULL AND (NOT isfinite(",
        quantities, ") OR ", quantities, " < 0))")
    if (nrow(DBI::dbGetQuery(connection, paste0("SELECT 1 FROM qc_input WHERE ",
        paste(bad_values, collapse = " OR "), " LIMIT 1")))) {
        cli::cli_abort("saved rollup contains negative or non-finite abundance")
    }
    identity_keys <- quote_name(c(roles$sample, roles$protein, roles$precursor))
    if (nrow(DBI::dbGetQuery(connection, paste0("SELECT 1 FROM qc_input GROUP BY ",
        paste(identity_keys, collapse = ", "), " HAVING count(*) > 1 LIMIT 1")))) {
        cli::cli_abort("duplicate precursor identities would inflate abundance")
    }
    group <- c(roles$sample, roles$protein, roles$peptide)
    group_sql <- paste(quote_name(group), collapse = ", ")
    modified <- quote_name(roles$modified)
    aggregate_sql <- paste0("SELECT ", group_sql,
        ", sum(", quantities[[1L]], ") AS \"Peptide.RawQuantity\"",
        ", sum(", quantities[[2L]], ") AS \"Peptide.Normalised\"",
        ", CAST(count(DISTINCT CASE WHEN regexp_full_match(CAST(", modified,
        " AS VARCHAR), '[ \\t\\r\\n]*') THEN NULL ELSE CAST(", modified,
        " AS VARCHAR) END) AS INTEGER) AS peptidoform_count",
        ", CAST(count(", quantities[[1L]], ") AS INTEGER) AS raw_observed",
        ", CAST(count(", quantities[[2L]], ") AS INTEGER) AS normalised_observed",
        ", CAST(count(*) AS INTEGER) AS source_count",
        " FROM qc_input GROUP BY ", group_sql)
    aggregated <- collectWorkflowRollupQuery(connection, aggregate_sql, metadata$rows)
    missingness <- list(
        raw_partial_missing_outputs = sum(aggregated$raw_observed > 0L &
            aggregated$raw_observed < aggregated$source_count),
        raw_all_missing_outputs = sum(aggregated$raw_observed == 0L),
        normalised_partial_missing_outputs = sum(aggregated$normalised_observed > 0L &
            aggregated$normalised_observed < aggregated$source_count),
        normalised_all_missing_outputs = sum(aggregated$normalised_observed == 0L))
    # Grouping creates new rows: restore column types, not source element/row names.
    derived <- metadata
    derived$name_columns <- character()
    for (field in names(derived$column_metadata)) {
        if (!is.null(derived$column_metadata[[field]]$attribute_order))
            derived$column_metadata[[field]]$attribute_order <- setdiff(derived$column_metadata[[field]]$attribute_order, "names")
    }
    derived$attributes <- list()
    derived$class <- c("tbl_df", "tbl", "data.frame")
    if (!is.null(derived$attribute_order)) derived$attribute_order <- c("names", "row.names", "class")
    evidence <- input$evidence
    counts <- annotations <- NULL
    if (length(evidence)) {
        evidence_sql <- paste(quote_name(c(roles$protein, evidence)), collapse = ", ")
        conflict_sql <- paste0("SELECT 1 FROM (SELECT DISTINCT ", evidence_sql,
            " FROM qc_input) AS evidence GROUP BY ", quote_name(roles$protein),
            " HAVING count(*) <> 1 LIMIT 1")
        if (nrow(DBI::dbGetQuery(connection, conflict_sql))) {
            cli::cli_abort("conflicting frozen identification evidence for protein groups")
        }
        counts <- collectWorkflowRollupQuery(connection, paste0("SELECT DISTINCT ",
            evidence_sql, " FROM qc_input"), metadata$rows)
    }
    if (length(input$annotation)) {
        annotations <- collectWorkflowRollupQuery(connection, paste0("SELECT DISTINCT ",
            paste(quote_name(c(roles$protein, "Protein.Ids")), collapse = ", "),
            " FROM qc_input"), metadata$rows)
    }
    # No query needs the input now. Release the engine before R-side joins/sort.
    releaseWorkflowTableReadScope(scope)
    connection <- NULL
    output <- restoreSavedWorkflowTable(aggregated, derived, group)
    for (name in c("Peptide.RawQuantity", "Peptide.Normalised", "peptidoform_count")) {
        output[[name]] <- aggregated[[name]]
    }
    aggregated <- NULL
    if (length(evidence)) {
        counts <- restoreSavedWorkflowTable(counts, derived, c(roles$protein, evidence))
        output <- dplyr::left_join(output, counts, by = roles$protein)
        counts <- NULL
    }
    if (length(input$annotation)) {
        annotations <- restoreSavedWorkflowTable(annotations, derived, c(roles$protein, "Protein.Ids"))
        collapse_ids <- function(values) {
            tokens <- unlist(strsplit(as.character(values), ";", fixed = TRUE))
            tokens <- trimws(tokens)
            tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
            paste(sort(unique(tokens)), collapse = ";")
        }
        annotations <- annotations |>
            dplyr::group_by(.data[[roles$protein]]) |>
            dplyr::summarise(Protein.Ids = collapse_ids(.data$Protein.Ids), .groups = "drop")
        output <- dplyr::left_join(output, annotations, by = roles$protein)
        annotations <- NULL
    }
    output <- output |>
        dplyr::arrange(.data[[roles$sample]], .data[[roles$protein]], .data[[roles$peptide]])
    summary <- list(aggregation_rule = "linear_sum_modified_unmodified_and_charge_states",
        active_protein_key = roles$protein, precursor_identity_mode = "declared_precursor_id",
        input_precursor_rows = as.integer(metadata$rows),
        unique_precursor_identities = as.integer(metadata$rows),
        output_stripped_peptide_identities = nrow(output),
        raw_partial_missing_outputs = missingness$raw_partial_missing_outputs,
        raw_all_missing_outputs = missingness$raw_all_missing_outputs,
        normalised_partial_missing_outputs = missingness$normalised_partial_missing_outputs,
        normalised_all_missing_outputs = missingness$normalised_all_missing_outputs,
        frozen_identification_evidence = if (length(evidence)) "preserved" else "absent",
        protein_ids_provenance = if (length(input$annotation)) "preserved" else
            if (roles$protein == "Protein.Ids") "active_key" else "absent")
    list(data = output, summary = summary,
        execution = list(aggregation_frontend = "duckdb_sum_non_exact"))
}

#' Roll up the current saved precursor table and commit references only.
#' @param workflow Session workflow owner.
#' @return Reference-backed rollup progress, using genuine existing sum semantics.
#' @noRd
runDIAReferencePrecursorRollup <- function(workflow) {
    manager <- workflow$state_manager
    before <- readDIAStateMetadata(manager)
    pending <- if (manager$current_state == "precursor_rollup") resumeDIAQvalueOutput(manager, before) else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    # KIT-008: the in-process path runs the declaration's own backend (P14, the exact ordered-sum aggregate
    # the worker runs), no longer the legacy row-indexed path; the two were checked page-identical.
    value <- computeDIAReferencePrecursorRollup(manager, getDIAStateReference(manager), "P14")
    manager$.context$last_rollup_execution <- value$execution
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name,
        value$description, value$config, value$section, value$counts)
}

#' Apply the sample filter in process through the same declaration the managed worker runs (KIT-008):
#' the KIT-006 physical mode (adaptive complete-sample batches), then the shared finishing step.
#' @param workflow Session workflow owner.
#' @param parameters `peptides_per_sample_cutoff` and optional `inclusion_list`.
#' @return References, bounded progress and text; no scientific table escapes.
#' @noRd
runDIAReferenceSampleFilter <- function(workflow, parameters) {
    manager <- workflow$state_manager
    before <- readDIAStateMetadata(manager)
    pending <- if (identical(manager$current_state, "sample_filtered")) resumeDIAQvalueOutput(manager, before) else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    entry <- workflowActionSpec("minimum_peptides_per_sample")
    reference <- getDIAStateReference(manager)
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    parameters <- resolveWorkflowSampleFilterParameters(manager, reference, parameters)
    value <- computeDIAReferenceSampleFilter(manager, reference, parameters, output_mode = entry$default_variant)
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description,
        value$config, value$section, value$counts,
        observed_inclusion_absence = value$summary$inclusion_entries_absent_from_data)
}

#' Apply a verified full-population rollup result with the legacy S4 bookkeeping.
#' @param theObject Complete hydrated peptide state required by the unchanged audit.
#' @param rollup_result Native-shaped data/summary from the saved-table operation.
#' @param core_utilisation Existing recorded configuration; not an SQL engine setting.
#' @return Complete S4 state. Called as rollupFn to preserve legacy parameter provenance.
#' @noRd
applySavedWorkflowPrecursorRollup <- function(theObject, rollup_result, core_utilisation = NULL) {
    core_utilisation <- checkParamsObjectFunctionSimplify(theObject, "core_utilisation", NA)
    theObject <- updateParamInObject(theObject, "core_utilisation")
    methods::slot(theObject, "peptide_data") <- rollup_result$data
    methods::slot(theObject, "raw_quantity_column") <- "Peptide.RawQuantity"
    methods::slot(theObject, "norm_quantity_column") <- "Peptide.Normalised"
    args <- methods::slot(theObject, "args")
    args$rollUpPrecursorToPeptide$rollup_summary <- rollup_result$summary
    methods::slot(theObject, "args") <- args
    cleanDesignMatrixPeptide(theObject)
}

#' Apply the existing DIA peptide kernel to complete scientific groups.
#' @param workflow Session workflow owner.
#' @param section Existing peptide-method configuration section.
#' @param parameters Explicit GUI scientific parameters.
#' @param state_name,description Existing checkpoint identity and description.
#' @return References, bounded progress and text; no scientific table escapes.
#' @noRd
runDIAReferencePeptideGroups <- function(workflow, section, parameters, state_name, description) {
    manager <- workflow$state_manager
    before <- readDIAStateMetadata(manager)
    pending <- if (identical(manager$current_state, state_name)) resumeDIAQvalueOutput(manager, before) else NULL
    if (!is.null(pending)) return(pending)
    assertDIAExportReady(manager)
    context <- manager$.context
    checkpoint <- before$data
    checkpoint$dependencies <- collectDIANativeReferences(list(getDIAStateReference(manager),
        checkpoint$dependencies, checkpoint$components, checkpoint$support, checkpoint$products))
    source <- checkpoint$components$peptide_data
    object <- readDIAExecutionComponents(manager, fields = "design_matrix")
    args <- object@args[[section]] %||% list()
    if (section == "filterMinNumPeptidesPerProtein") {
        args$peptides_per_protein_cutoff <- parameters$num_peptides_per_protein_thresh
        args$peptidoforms_per_protein_cutoff <- parameters$num_peptidoforms_per_protein_thresh
    }
    updates <- parameters
    if (section == "peptideMissingValueImputation") {
        exclusion <- parameters$exclusion_column
        exclusion <- if (!length(exclusion) || is.na(exclusion[[1L]]) ||
            !nzchar(trimws(as.character(exclusion[[1L]])))) NULL else trimws(as.character(exclusion[[1L]]))
        parameters["exclusion_column"] <- list(exclusion)
        updates <- parameters
        if (is.null(exclusion)) updates$exclusion_column <- NULL
    }
    args[names(updates)] <- updates
    object@args[[section]] <- args
    ui_updates <- switch(section,
        filterMinNumPeptidesPerProtein = list(peptides_per_protein_cutoff = parameters$num_peptides_per_protein_thresh,
            peptidoforms_per_protein_cutoff = parameters$num_peptidoforms_per_protein_thresh),
        peptideIntensityFiltering = parameters[c("strict_mode", "min_reps_per_group", "min_groups",
            "peptides_intensity_cutoff_percentile")],
        filterMinNumPeptidesPerSample = parameters["peptides_per_sample_cutoff"],
        peptideMissingValueImputation = updates, list())
    for (parameter in names(ui_updates)) object <- updateConfigParameter(object,
        function_name = section, parameter_name = parameter, new_value = ui_updates[[parameter]])
    args <- object@args[[section]]
    design <- object@design_matrix
    scope <- newWorkflowTableReadScope(context$store, source)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    con <- scope$connection
    DBI::dbExecute(con, "SET memory_limit = '1GB'")
    DBI::dbExecute(con, "SET threads = 2")
    q <- \(name) as.character(DBI::dbQuoteIdentifier(con, name))
    if (any(c(".dia_source_order", ".dia_design_order", ".dia_feature_order") %in% scope$metadata$columns)) {
        cli::cli_abort("reserved DIA grouping columns are present", class = "saved_entity_invalid")
    }
    group <- if (section == "filterMinNumPeptidesPerSample") object@sample_id else object@protein_id_column
    scientific_columns <- intersect(scope$metadata$columns, unique(c(
        unlist(checkpoint$slots[c("sample_id", "protein_id_column", "peptide_sequence_column",
            "raw_quantity_column", "norm_quantity_column")]),
        "Modified.Sequence", "identification_peptide_count", "identification_peptidoform_count",
        args$peptide_quantity_column)))
    threshold <- NULL
    if (section == "peptideIntensityFiltering") {
        quantity <- args$peptide_quantity_column %||% object@norm_quantity_column
        values <- materializeDIATable(context, source, quantity, scope = scope)[[quantity]]
        values <- values[is.finite(values)]
        percentile <- parameters$peptides_intensity_cutoff_percentile
        if (length(percentile) != 1L || !is.finite(percentile) || percentile < 0 || percentile > 100) {
            cli::cli_abort("invalid intensity percentile", class = "saved_entity_invalid")
        }
        threshold <- if (length(values)) unname(stats::quantile(values,
            probs = percentile / 100, na.rm = TRUE)[[1L]]) else 0
        values <- NULL
    }
    offsets <- list(data = 0L)
    restoration <- list()
    summary <- NULL
    sum_fields <- c("input_rows", "output_rows", "observed_rows", "generated_rows", "imputed_rows",
        "missing_not_imputed_rows", "eligible_group_features", "input_feature_count",
        "retained_feature_count", "removed_feature_count", "non_finite_measurement_count",
        "input_sample_count", "retained_sample_count", "removed_sample_count")
    result_args <- args
    cursor <- openWorkflowTableGroups(scope, scientific_columns, group)
    checkpoint$dependencies <- collectDIANativeReferences(list(checkpoint$dependencies, cursor$dependencies))
    on.exit(cursor$close(), add = TRUE)
    generated_rows <- 0L
    repeat {
        part <- cursor$next_group()
        if (is.null(part)) break
        positions <- attr(part, "workflow_source_order") %||% integer()
        attr(part, "workflow_source_order") <- NULL
        part$.dia_source_order <- positions
        object@peptide_data <- part
        object@design_matrix <- design
        if (section == "peptideIntensityFiltering") {
            result <- peptideIntensityFilteringHelper(part, design,
                min_peptide_intensity_threshold = threshold, sample_id_column = object@sample_id,
                grouping_variable = parameters$grouping_variable %||% "group",
                protein_id_column = object@protein_id_column,
                peptide_sequence_column = object@peptide_sequence_column,
                peptide_quantity_column = args$peptide_quantity_column %||% object@norm_quantity_column,
                core_utilisation = args$core_utilisation %||% NA,
                min_reps_per_group = parameters$min_reps_per_group,
                min_groups = parameters$min_groups, strict_mode = parameters$strict_mode,
                groupwise_percentage_cutoff = args$groupwise_percentage_cutoff %||% 1,
                max_groups_percentage_cutoff = args$max_groups_percentage_cutoff %||% 50,
                return_filter_result = TRUE, require_complete_design_runs = FALSE)
            result_args <- utils::modifyList(args, list(
                grouping_variable = parameters$grouping_variable %||% "group",
                groupwise_percentage_cutoff = args$groupwise_percentage_cutoff %||% 1,
                max_groups_percentage_cutoff = args$max_groups_percentage_cutoff %||% 50,
                core_utilisation = args$core_utilisation %||% NA,
                peptide_quantity_column = args$peptide_quantity_column %||% object@norm_quantity_column,
                min_peptide_intensity_threshold = threshold))
            result_args$support_table <- result$support_table
            result_args$removal_ledger <- result$removal_ledger
            result_args$filter_summary <- result$summary
        } else {
            object <- do.call(get(section, mode = "function"), c(list(theObject = object), parameters))
            result_args <- object@args[[section]]
            result <- list(data = object@peptide_data, support_table = result_args$support_table,
                removal_ledger = result_args$removal_ledger,
                summary = result_args$filter_summary %||% result_args$imputation_summary)
        }
        if (!is.null(result$summary)) {
            if (is.null(summary)) summary <- result$summary else {
                for (name in intersect(sum_fields, names(result$summary))) {
                    summary[[name]] <- (summary[[name]] %||% 0L) + result$summary[[name]]
                }
            }
        }
        for (field in c("data", "support_table", "removal_ledger")) {
            table <- result[[field]]
            if (!is.data.frame(table) || !ncol(table)) next
            if (field == "data") {
                table$.dia_design_order <- match(as.character(table[[object@sample_id]]),
                    as.character(design[[object@sample_id]]))
                table$.dia_feature_order <- rep(0L, nrow(table))
                technical <- object@technical_replicate_id
                if (section == "peptideMissingValueImputation" && length(technical) == 1L &&
                    technical %in% names(design) && nrow(part)) {
                    observed <- data.frame(.group = as.character(design[[technical]][match(
                        as.character(part[[object@sample_id]]), as.character(design[[object@sample_id]]))]),
                        .protein = as.character(part[[object@protein_id_column]]),
                        .peptide = as.character(part[[object@peptide_sequence_column]]),
                        .position = part$.dia_source_order)
                    observed <- observed |>
                        dplyr::group_by(.data$.group, .data$.protein, .data$.peptide) |>
                        dplyr::summarise(.position = min(.data$.position), .groups = "drop")
                    selectors <- data.frame(.group = as.character(design[[technical]][table$.dia_design_order]),
                        .protein = as.character(table[[object@protein_id_column]]),
                        .peptide = as.character(table[[object@peptide_sequence_column]]))
                    table$.dia_feature_order <- dplyr::left_join(selectors, observed,
                        by = c(".group", ".protein", ".peptide"))$.position
                }
            }
            if (is.null(restoration[[field]])) {
                restoration[[field]] <- workflowSavedTableMetadata(table, ".rownames", 3L)
                offsets[[field]] <- 0L
            }
            if (field == "data") {
                generated <- is.na(table$.dia_source_order)
                if (any(generated)) {
                    additions <- table[generated, , drop = FALSE]
                    for (start in seq.int(1L, nrow(additions), by = 10000L)) {
                        at <- seq.int(start, min(nrow(additions), start + 9999L))
                        generated_rows <- appendWorkflowTableBatch(scope, "dia_group_generated",
                            additions[at, , drop = FALSE], generated_rows)
                    }
                }
                extra <- unique(c(setdiff(names(table), c(scope$metadata$columns,
                    ".dia_design_order", ".dia_feature_order")),
                    if (section == "peptideMissingValueImputation") intersect(names(table),
                        c(result_args$imputed_value_column, "is_imputed"))))
                table <- table[!generated, extra, drop = FALSE]
            }
            for (start in seq.int(1L, max(1L, nrow(table)), by = 10000L)) {
                at <- if (nrow(table)) seq.int(start, min(nrow(table), start + 9999L)) else integer()
                offsets[[field]] <- appendWorkflowTableBatch(scope, paste0("dia_group_", field),
                    table[at, , drop = FALSE], offsets[[field]],
                    key_column = if (field == "data") ".dia_source_order" else ".rownames")
            }
        }
        object@peptide_data <- object@peptide_data[0, , drop = FALSE]
        result <- part <- table <- NULL
    }
    cursor$close()
    identity <- newSavedWorkflowScientificIdentity(source$entity_id,
        checkpoint$components$design_matrix$entity_id, digest::digest(list(section, parameters), algo = "sha256"),
        "complete-declared-groups", context$software)
    references <- list()
    counts <- NULL
    for (field in names(restoration)) {
        path <- finishWorkflowTableBatch(scope, paste0("dia_group_", field), identity)
        if (field == "data") {
            metadata <- restoration[[field]]
            original <- scope$metadata
            internal <- c(".dia_source_order", ".dia_design_order", ".dia_feature_order")
            metadata$columns <- unique(c(original$columns, setdiff(metadata$columns, internal)))
            specs <- c(original$column_metadata[setdiff(original$columns, names(metadata$column_metadata))],
                metadata$column_metadata)
            metadata$column_metadata <- specs[metadata$columns]
            named <- union(setdiff(names(original$name_columns), scientific_columns),
                intersect(names(metadata$name_columns), metadata$columns))
            metadata$name_columns <- stats::setNames(paste0(".multischolar_names_", seq_along(named)), named)
            if (!length(named)) metadata$name_columns <- stats::setNames(character(), character())
            metadata$rows <- offsets[[field]] + generated_rows
            restoration[[field]] <- metadata
            generated_path <- if (generated_rows) finishWorkflowTableBatch(scope, "dia_group_generated", identity) else NULL
            if (!generated_rows && workflowTableChangesAreIdentity(scope, path, metadata,
                ".dia_source_order")) {
                references[[field]] <- source
                next
            }
            accumulator <- newDIAFilteringAccumulator(readSavedWorkflowTablePrototype(scope),
                checkpoint$slots$protein_id_column)
            references[[field]] <- saveWorkflowTable(scope$store, NULL, basename(tempfile("dia-data-")),
                identity, metadata$key_column, metadata = metadata,
                write_fn = function(destination, restoration) writeWorkflowTableChanges(scope,
                    path, generated_path, restoration, destination, source_column = ".dia_source_order",
                    generated_order = c(".dia_design_order", ".dia_feature_order"),
                    consume = accumulator$push))
            counts <- accumulator$finish()
            next
        }
        metadata <- restoration[[field]]
        metadata$rows <- offsets[[field]]
        order <- if (field == "data")
            '".dia_source_order" NULLS LAST, ".dia_design_order", ".dia_feature_order", ".multischolar_order"' else {
                order_fields <- switch(section,
                    filterMinNumPeptidesPerSample = object@sample_id,
                    peptideIntensityFiltering = c(object@protein_id_column, object@peptide_sequence_column,
                        parameters$grouping_variable %||% "group"),
                    removePeptidesWithOnlyOneReplicate = c(object@protein_id_column,
                        object@peptide_sequence_column, parameters$replicate_group_column),
                    peptideMissingValueImputation = c("technical_replicate_group", "protein_identity", "peptide_identity"))
                paste(q(c(intersect(order_fields, metadata$columns), ".multischolar_order")), collapse = ", ")
            }
        columns <- setdiff(metadata$columns, c(".dia_source_order", ".dia_design_order", ".dia_feature_order"))
        metadata$columns <- columns
        metadata$column_metadata <- metadata$column_metadata[columns]
        metadata$name_columns <- metadata$name_columns[intersect(names(metadata$name_columns), columns)]
        metadata$schema <- 3L
        metadata$codec <- "parquet_duckdb_v3"
        source_scope <- newWorkflowTableReadScope(scope$store, path, reader = "arrow", metadata_only = TRUE)
        order_columns <- c(intersect(order_fields, source_scope$metadata$columns), ".multischolar_order")
        ordering <- saveWorkflowOrderIndex(source_scope, order_columns)
        if (is.null(ordering$value$positions)) {
            references[[field]] <- saveWorkflowTableComponent(source_scope, metadata,
                basename(tempfile(paste0("dia-", field, "-"))), identity)
            releaseWorkflowTableReadScope(source_scope)
            next
        }
        releaseWorkflowTableReadScope(source_scope)
        references[[field]] <- saveWorkflowTable(scope$store, NULL,
            basename(tempfile(paste0("dia-", field, "-"))), identity, metadata$key_column,
            metadata = metadata, write_fn = function(destination, restoration) {
                cursor <- openWorkflowOrderedRows(scope$store, path, columns, order_columns)
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
    }
    for (field in intersect(names(result_args), setdiff(names(references), "data")))
        result_args[field] <- list(NULL)
    if (!is.null(summary)) result_args[[if (section == "peptideMissingValueImputation")
        "imputation_summary" else "filter_summary"]] <- summary
    checkpoint$slots$args[[section]] <- result_args
    checkpoint$support[[section]] <- references[setdiff(names(references), "data")]
    for (field in names(checkpoint$support[[section]])) checkpoint$slots$args[[section]][field] <- list(NULL)
    checkpoint$components$peptide_data <- references$data
    checkpoint$dimensions <- list(rows = restoration$data$rows,
        columns = restoration$data$columns)
    releaseWorkflowTableReadScope(scope)
    finishDIAReferencePeptideStep(workflow, checkpoint, state_name, description,
        c(parameters, list(filter_summary = summary)), section, counts)
}

#' Commit a peptide table and its bounded progress without reading the assay.
#' @param workflow Session owner.
#' @param checkpoint Accepted component result.
#' @param state_name,description,config Existing scientific state bookkeeping.
#' @param section Existing configuration section.
#' @param observed_inclusion_absence Optional globally validated absence report from the worker.
#' @return Display result containing references and summaries only.
#' @noRd
finishDIAReferencePeptideStep <- function(workflow, checkpoint, state_name, description, config, section,
    counts = NULL, observed_inclusion_absence = NULL) {
    manager <- workflow$state_manager
    context <- manager$.context
    scope <- newWorkflowTableReadScope(context$store, checkpoint$components$peptide_data, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    prototype <- readDIAExecutionComponents(manager, fields = "design_matrix")
    counted <- summarizeSavedDIAFiltering(scope, prototype, checkpoint$slots$protein_id_column, counts)
    checkpoint$dependencies <- collectDIANativeReferences(list(checkpoint$dependencies, counted$dependencies))
    summary <- counted$summary
    runs <- counted$runs
    sample <- checkpoint$slots$sample_id
    checkpoint$dependencies <- collectDIANativeReferences(list(checkpoint$dependencies, counted$dependencies))
    design <- materializeDIATable(context, checkpoint$components$design_matrix)
    removed <- setdiff(as.character(design[[sample]]), as.character(runs))
    prototype <- readDIAExecutionComponents(manager, fields = "design_matrix")
    probe <- prototype@peptide_data[rep(NA_integer_, length(runs)), , drop = FALSE]
    probe[[sample]] <- runs
    prototype@peptide_data <- probe
    prototype@design_matrix <- design
    design <- cleanDesignMatrixPeptide(prototype)@design_matrix
    prototype <- probe <- NULL
    checkpoint$components$design_matrix <- saveDIATableComponent(context, design)
    if (section == "filterMinNumPeptidesPerSample") {
        checkpoint$slots$args[[section]]$samples_removed <- removed
        checkpoint$slots$args[[section]]$samples_removed_count <- length(removed)
        checkpoint$slots$args[[section]]$samples_before_count <- length(runs) + length(removed)
        checkpoint$slots$args[[section]]$samples_after_count <- length(runs)
        inclusion <- checkpoint$slots$args[[section]]$inclusion_list %||% character()
        checkpoint$slots$args[[section]]$filter_summary$inclusion_entries_absent_from_data <-
            observed_inclusion_absence %||% setdiff(inclusion, c(as.character(runs), removed))
    }
    protein_count <- counted$proteins
    handoff <- buildDIAReferencePeptideHandoff(checkpoint, section, config, protein_count,
        runs, removed, description)
    config <- handoff$config
    description <- handoff$description
    result_text <- handoff$text
    workflow$qc_params$peptide_qc[[handoff$qc_name]] <- handoff$qc
    prior <- shiny::isolate(workflow$filtering_progress)
    prior_slots <- if (methods::is(prior, "FilteringProgress")) stats::setNames(
        lapply(methods::slotNames(prior), \(name) methods::slot(prior, name)), methods::slotNames(prior)) else NULL
    product <- list(filteringSummary = summary, resultText = result_text, priorProgress = prior_slots,
        qcParams = workflow$qc_params)
    identity <- newSavedWorkflowScientificIdentity(checkpoint$components$peptide_data$entity_id,
        checkpoint$components$design_matrix$entity_id, digest::digest(config, algo = "sha256"),
        paste0(state_name, "-progress"), context$software)
    checkpoint$products$qc_progress <- saveWorkflowConsumerEntity(context$store,
        newSavedWorkflowScientificProduct(product, paste0("dia-", state_name, "-progress")),
        basename(tempfile("dia-qc-progress-")), identity, "scientific_product")
    releaseWorkflowTableReadScope(scope)
    reference <- commitDIAReferenceState(manager, state_name, checkpoint, config, description)
    c(list(stateReference = reference), product)
}

#' Preserve existing GUI configuration, summaries and wording for saved peptide steps.
#' @param checkpoint,section,parameters Accepted metadata and existing method inputs.
#' @param proteins,runs,removed Exact accepted population summaries.
#' @param description Original checkpoint description.
#' @return Existing config/QC records and scalar display text.
#' @noRd
buildDIAReferencePeptideHandoff <- function(checkpoint, section, parameters, proteins, runs, removed, description) {
    args <- checkpoint$slots$args[[section]]
    summary <- args$filter_summary
    stamp <- Sys.time()
    if (section == "filterMinNumPeptidesPerProtein") {
        peptides <- parameters$num_peptides_per_protein_thresh
        forms <- parameters$num_peptidoforms_per_protein_thresh
        aliases <- c(peptides = peptides, peptidoforms = forms)
        config <- list(min_peptides_per_protein = peptides, min_peptidoforms_per_protein = forms,
            peptides_per_protein_cutoff = peptides, peptidoforms_per_protein_cutoff = forms,
            num_peptides_per_protein_thresh = peptides, num_peptidoforms_per_protein_thresh = forms,
            parameter_alias_resolution = list(displayed = aliases, configured = aliases,
                resolved = aliases, executed = aliases))
        qc <- list(min_peptides_per_protein = peptides, min_peptidoforms_per_protein = forms, timestamp = stamp)
        text <- paste("Protein Identification Evidence Filter Applied Successfully\n",
            "========================================================\n", sprintf("Proteins remaining: %d\n", proteins),
            sprintf("Min distinct peptides per protein: %d\n", peptides),
            sprintf("Min distinct peptidoforms per protein: %d\n", forms), "State saved as: 'protein_peptide_filtered'\n")
        name <- "protein_peptide_filter"
    } else if (section == "peptideIntensityFiltering") {
        strict <- parameters$strict_mode
        config <- list(strict_mode = strict, min_reps_per_group = parameters$min_reps_per_group,
            min_groups = parameters$min_groups, intensity_cutoff_percentile = parameters$peptides_intensity_cutoff_percentile,
            filter_summary = summary)
        qc <- config
        if (strict) qc$min_reps_per_group <- qc$min_groups <- NA
        qc$timestamp <- stamp
        description <- if (strict) "Applied STRICT peptide intensity filter" else "Applied FLEXIBLE peptide intensity filter"
        text <- paste("Intensity Filter Applied Successfully\n", "====================================\n",
            sprintf("Mode: %s\n", if (strict) "STRICT" else "FLEXIBLE"), sprintf("Proteins remaining: %d\n", proteins),
            sprintf("Intensity cutoff percentile: %.1f%%\n", parameters$peptides_intensity_cutoff_percentile),
            if (strict) "Rule: every design run must be finite and at/above threshold\n" else
                sprintf("Group rule: >= %d distinct observed runs\n", parameters$min_reps_per_group),
            if (strict) "" else sprintf("Peptide rule: >= %d passing groups\n", parameters$min_groups),
            sprintf("Executed threshold: %g from %s\n", summary$intensity_threshold, summary$intensity_quantity_column),
            sprintf("Peptide identities retained/removed: %d / %d\n", summary$retained_feature_count, summary$removed_feature_count),
            "State saved as: 'intensity_filtered'\n")
        name <- "intensity_filter"
    } else if (section == "filterMinNumPeptidesPerSample") {
        config <- list(min_peptides_per_sample = parameters$peptides_per_sample_cutoff,
            samples_removed = removed, samples_removed_count = length(removed), filter_summary = summary)
        qc <- list(min_peptides_per_sample = parameters$peptides_per_sample_cutoff, samples_removed = removed,
            samples_removed_count = length(removed), samples_before_count = length(runs) + length(removed),
            samples_after_count = length(runs), filter_summary = summary, timestamp = stamp)
        text <- paste("Sample Quality Filter Applied Successfully\n", "=========================================\n",
            sprintf("Proteins remaining: %d\n", proteins), sprintf("Samples remaining: %d\n", length(runs)),
            sprintf("Samples removed: %d\n", length(removed)),
            sprintf("Min distinct protein-group/stripped-peptide identities per sample: %d\n", parameters$peptides_per_sample_cutoff),
            "State saved as: 'sample_filtered'\n")
        if (length(removed)) text <- paste0(text, "\nRemoved samples:\n", paste(removed, collapse = ", "))
        name <- "sample_filter"
    } else if (section == "removePeptidesWithOnlyOneReplicate") {
        config <- list(replicate_group_column = parameters$replicate_group_column, minimum_distinct_runs = 2L,
            retention_policy = "supported_in_any_group", filter_summary = summary)
        qc <- c(config, list(timestamp = stamp))
        text <- paste("Replicate Filter Applied Successfully\n", "====================================\n",
            sprintf("Proteins remaining: %d\n", proteins), sprintf("Replicate group column: %s\n", parameters$replicate_group_column),
            "Rule: >= 2 distinct runs in any replicate group\n",
            "Retention: globally keep supported features, including singleton observations in other groups\n",
            sprintf("Peptide identities retained/removed: %d / %d\n", summary$retained_feature_count, summary$removed_feature_count),
            "State saved as: 'replicate_filtered'\n")
        name <- "replicate_filter"
    } else if (section == "peptideMissingValueImputation") {
        audit <- .peptideImputationAuditSummary(checkpoint, parameters$proportion_missing_values, parameters$exclusion_column)
        config <- list(proportion_missing_values = parameters$proportion_missing_values,
            exclusion_column = parameters$exclusion_column, imputation_summary = audit)
        qc <- c(config, list(timestamp = stamp))
        skipped <- identical(audit$status, "skipped")
        if (skipped) description <- paste("Skipped technical-replicate imputation:", audit$skip_reason)
        lines <- c(sprintf("Eligibility: missing fraction %s %.1f\n", audit$eligibility_operator, parameters$proportion_missing_values),
            sprintf("Denominator: %s\n", audit$eligibility_denominator), sprintf("Quantity column: %s\n", audit$quantity_column),
            sprintf("Exclusion source: %s\n", audit$exclusion_source))
        if (!is.na(audit$imputed_rows)) lines <- c(lines, sprintf("Rows imputed: %d\n", audit$imputed_rows))
        if (skipped) lines <- c(lines, sprintf("Skip reason: %s\n", audit$skip_reason))
        text <- paste(if (skipped) "Missing Value Imputation Skipped\n" else "Missing Value Imputation Applied Successfully\n",
            "============================================\n", sprintf("Proteins remaining: %d\n", proteins),
            sprintf("Max proportion missing: %.1f\n", parameters$proportion_missing_values), paste(lines, collapse = ""),
            "State saved as: 'imputed'\n", "\nReady for peptide-to-protein rollup step.")
        name <- "imputation"
    } else {
        summary <- args$rollup_summary
        config <- list(rollup_summary = summary)
        qc <- list(applied = TRUE, timestamp = stamp, summary = summary)
        text <- paste("Precursor Rollup Applied Successfully\n", "====================================\n",
            sprintf("Proteins remaining: %d\n", proteins), paste0(.formatPeptideRollupSummary(summary), "\n"),
            "State saved as: 'precursor_rollup'\n")
        name <- "precursor_rollup"
    }
    list(config = config, qc = qc, qc_name = name, text = text, description = description)
}

#' Stream the exact IQ input population, aliases and compatibility confidence values.
#' @param manager,reference Current owner and peptide table descriptor.
#' @param path Existing required IQ TSV destination.
#' @param mapping Ordered original/alias sample table.
#' @param sample_column,intensity_column Existing IQ roles.
#' @param write_fn Existing TSV writer.
#' @return Destination path after the completed file is durably published.
#' @noRd
writeDIAIQInput <- function(manager, reference, path, mapping, sample_column,
    intensity_column, write_fn) {
    scope <- newWorkflowTableReadScope(manager$.context$store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    temporary <- tempfile("iq-input-", tmpdir = dirname(path), fileext = ".tsv")
    on.exit(unlink(temporary), add = TRUE)
    streamWorkflowTableRows(scope, scope$metadata$columns, function(data, offset, row_names) {
        data$Q.Value <- rep(0.0009, nrow(data))
        data$PG.Q.Value <- rep(0.009, nrow(data))
        data[[intensity_column]][is.na(data[[intensity_column]])] <- 0
        selected <- match(as.character(data[[sample_column]]), mapping$Original)
        if (anyNA(selected)) cli::cli_abort("peptide data contains runs absent from the design matrix")
        data[[sample_column]] <- mapping$Alias[selected]
        write_fn(data, temporary, append = offset > 0L, col_names = offset == 0L)
    })
    releaseWorkflowTableReadScope(scope)
    syncWorkflowExportPaths(c(temporary, dirname(temporary)), dirname(path))
    if (!file.rename(temporary, path)) cli::cli_abort("cannot publish IQ input", class = "saved_entity_write")
    syncWorkflowExportPaths(c(path, dirname(path)), dirname(path))
    invisible(path)
}

#' Read only distinct protein-lineage annotations for IQ output construction.
#' @param manager,reference Current owner and peptide assay descriptor.
#' @param protein_column Active stable protein key.
#' @return Ordered annotation table, assembled from bounded reads.
#' @noRd
readDIAProteinLineage <- function(manager, reference, protein_column) {
    spec <- readDIAStateMetadata(manager)$data
    if (identical(spec$components$peptide_data, reference)) {
        lineage <- spec$support$.arguments$protein_accession_provenance
        if (isSavedWorkflowReference(lineage)) return(materializeDIATable(manager$.context, lineage))
        lineage <- spec$slots$args$protein_accession_provenance
        if (is.data.frame(lineage) && all(c(protein_column, "Protein.Ids") %in% names(lineage))) return(lineage)
    }
    scope <- newWorkflowTableReadScope(manager$.context$store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    fields <- intersect(unique(c(protein_column, "Protein.Ids")), scope$metadata$columns)
    cursor <- openWorkflowTableGroups(scope, fields, protein_column)
    on.exit(cursor$close(), add = TRUE)
    parts <- list()
    repeat {
        part <- cursor$next_group()
        if (is.null(part)) break
        positions <- attr(part, "workflow_source_order")
        attr(part, "workflow_source_order") <- NULL
        selected <- which(!duplicated(part[fields]))
        piece <- as.data.frame(part[selected, fields, drop = FALSE])
        for (field in fields) piece[[field]] <- as.character(piece[[field]])
        piece$.source_order <- positions[selected]
        parts[[length(parts) + 1L]] <- piece
        if (as.numeric(object.size(parts)) > 67108864)
            cli::cli_abort("protein lineage exceeds lookup limits", class = "saved_entity_limit")
    }
    result <- do.call(rbind, parts)
    result <- result[order(result$.source_order), fields, drop = FALSE]
    rownames(result) <- NULL
    result
}

#' Release a protein scientific result before delivery to its GUI observer.
#' @param workflow Session workflow owner.
#' @param result Existing result list containing one already-committed protein object.
#' @param field Name of that temporary S4 field.
#' @return Bounded display result and state reference; legacy callers are unchanged.
#' @noRd
finishDIAProteinAction <- function(workflow, result, field) {
    manager <- workflow$state_manager
    if (!usesDIAReferenceStates(manager)) return(result)
    object <- result[[field]]
    summary <- structure(list(is_protein_quant = TRUE,
        proteins = countUniqueProteins(object), proteins_per_run = countProteinsPerRun(object),
        total_peptides = NA_integer_, peptides_per_protein = data.frame(),
        peptides_per_run = data.frame()), class = "DIAFilteringSummary")
    result[[field]] <- NULL
    result$stateReference <- getDIAStateReference(manager)
    summary <- readDIAProteinProgress(manager, result$stateReference, summary)
    result$filteringSummary <- summary
    result$qcParams <- workflow$qc_params
    result$accessionCleanupResults <- workflow$accession_cleanup_results
    prior <- shiny::isolate(workflow$filtering_progress)
    result$priorProgress <- if (methods::is(prior, "FilteringProgress")) stats::setNames(
        lapply(methods::slotNames(prior), \(name) methods::slot(prior, name)), methods::slotNames(prior)) else NULL
    identity <- newSavedWorkflowScientificIdentity(result$stateReference$entity_id,
        "accepted-protein-state", "required-qc-outputs", manager$current_state, manager$.context$software)
    product <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowScientificProduct(result, "dia-protein-qc-output"),
        basename(tempfile("dia-protein-output-")), identity, "scientific_product")
    snapshot <- manager$.context$snapshot
    snapshot$output_products[[result$stateReference$entity_id]] <- product
    publishDIAStateHead(manager$.context, snapshot)
    publishDIAProteinTables(manager, result)
    result
}

#' Complete protein file publication using accepted results, without rerunning science.
#' @param manager Current owner.
#' @param result Accepted output product and state reference.
#' @return NULL invisibly; previous completed files remain recoverable on failure.
#' @noRd
publishDIAProteinTables <- function(manager, result) {
    state <- manager$current_state
    if (!state %in% c("protein_s4_created", "protein_accession_cleaned", "protein_replicate_filtered"))
        return(invisible(NULL))
    record <- readDIAStateMetadata(manager, result$stateReference)
    paths <- manager$.context$experiment_paths
    destinations <- if (state == "protein_s4_created") paths[c("protein_qc_dir", "peptide_qc_dir")] else
        if (state == "protein_accession_cleaned") paths["source_dir"] else paths[c("protein_qc_dir", "source_dir")]
    publishDIAActionOutputs(manager, result$stateReference, "protein-qc", destinations, function(staged) {
        if (state == "protein_s4_created") {
            files <- list(iq_output_file = file.path(staged$protein_qc_dir, "iq_output_file.txt"),
                peptide_input_file = file.path(staged$peptide_qc_dir, "peptide_values_imputed.tsv"))
            for (field in names(files)) if (!file.copy(record$config[[field]], files[[field]], overwrite = TRUE))
                cli::cli_abort("cannot stage IQ output publication", class = "saved_entity_outputs_pending")
        } else if (state == "protein_accession_cleaned") {
            saveRDS(result$accessionCleanupResults, file.path(staged$source_dir, "accession_cleanup_results.RDS"))
        } else {
            writeDIADesignTable(manager, record$data$components$protein_quant_table,
                file.path(staged$protein_qc_dir, "remove_proteins_with_only_one_rep.tsv"), vroom::vroom_write)
            saveRDS(result$qcParams, file.path(staged$source_dir, "qc_params.RDS"))
        }
        list(completed = TRUE)
    })
    invisible(NULL)
}

#' Retry a committed protein action's required files and figures only.
#' @param workflow Current workflow owner.
#' @param state_name Expected accepted protein action.
#' @return Pending reference product or NULL for a new scientific action.
#' @noRd
resumeDIAProteinAction <- function(workflow, state_name) {
    manager <- workflow$state_manager
    if (!usesDIAReferenceStates(manager) || manager$current_state != state_name) return(NULL)
    result <- resumeDIAQvalueOutput(manager, readDIAStateMetadata(manager))
    if (is.null(result)) return(NULL)
    workflow$qc_params <- result$qcParams
    if (!is.null(result$accessionCleanupResults)) workflow$accession_cleanup_results <- result$accessionCleanupResults
    publishDIAProteinTables(manager, result)
    result
}

#' Publish a protein QC display from its accepted bounded result.
#' @param result Accepted protein references and progress.
#' @param step_name,omic_type,experiment_label Existing output identities.
#' @param render_fn Existing filtering-progress renderer.
#' @param legacy_object Existing native object when reference mode is inactive.
#' @return Display grid, after required files have been acknowledged.
#' @noRd
publishDIAProteinProgress <- function(result, step_name, omic_type, experiment_label,
    render_fn, legacy_object = NULL) {
    arguments <- list(data = result$filteringSummary %||% legacy_object, step_name = step_name,
        omic_type = omic_type, experiment_label = experiment_label, return_grid = TRUE, overwrite = TRUE)
    if (is.null(result$stateReference)) return(do.call(render_fn, arguments))
    manager <- getDIAActiveManager()
    directories <- if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE))
        get("project_dirs", envir = .GlobalEnv) else list()
    key <- paste(omic_type, experiment_label, sep = "_")
    if (is.null(directories[[key]])) directories[[key]] <- manager$.context$experiment_paths
    arguments$project_dirs <- directories
    publishDIAQvalueOutput(manager, result, arguments, render_fn)
}

#' Roll up complete samples with native kernels and no full-assay SQL aggregation table.
#' @param store Existing immutable store.
#' @param input Validated source/column roles and provenance fields.
#' @param persist Native output prototype and scientific identity.
#' @return Accepted Parquet, native summary and reusable index dependencies.
#' @noRd
persistIndexedWorkflowPrecursorRollup <- function(store, input, persist) {
    scope <- newWorkflowTableReadScope(store, input$reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    roles <- input$roles
    fields <- unique(c(unlist(roles), input$evidence, input$annotation))
    if (!all(fields %in% scope$metadata$columns) ||
        any(vapply(scope$metadata$column_metadata[c(roles$raw, roles$normalised)],
            function(column) !identical(column$typeof, "double"), logical(1))))
        cli::cli_abort("saved rollup requires declared double quantities", class = "saved_entity_unsupported")
    dependencies <- list(input$reference)
    annotations <- NULL
    if (length(input$evidence) || length(input$annotation)) {
        protein_fields <- unique(c(roles$protein, input$evidence, input$annotation))
        protein_cursor <- openWorkflowTableGroups(scope, protein_fields, roles$protein)
        on.exit(protein_cursor$close(), add = TRUE)
        dependencies <- c(dependencies, protein_cursor$dependencies)
        annotation_parts <- list()
        repeat {
            part <- protein_cursor$next_group()
            if (is.null(part)) break
            attr(part, "workflow_source_order") <- NULL
            if (length(input$evidence)) {
                evidence <- dplyr::distinct(dplyr::select(part, dplyr::all_of(c(roles$protein, input$evidence))))
                if (anyDuplicated(evidence[[roles$protein]]))
                    cli::cli_abort("conflicting frozen identification evidence for protein groups")
            }
            if (length(input$annotation)) {
                annotation <- part |>
                    dplyr::group_by(.data[[roles$protein]]) |>
                    dplyr::summarise(Protein.Ids = collapseRollupProteinIds(.data$Protein.Ids), .groups = "drop")
                if (is.factor(annotation[[roles$protein]])) annotation[[roles$protein]] <- as.character(annotation[[roles$protein]])
                annotation_parts[[length(annotation_parts) + 1L]] <- annotation
                if (as.numeric(object.size(annotation_parts)) > 67108864)
                    cli::cli_abort("rollup annotation lookup exceeds limits", class = "saved_entity_limit")
            }
        }
        protein_cursor$close()
        if (length(input$annotation)) annotations <- dplyr::bind_rows(annotation_parts)
    }
    sample_fields <- unique(c(unlist(roles), input$evidence))
    cursor <- openWorkflowTableGroups(scope, sample_fields, roles$sample, group_order = "key")
    on.exit(cursor$close(), add = TRUE)
    dependencies <- c(dependencies, cursor$dependencies)
    metadata <- workflowSavedTableMetadata(persist$prototype, ".rownames", 3L)
    metadata$rows <- NULL
    count_fields <- unique(c("Run", "Protein.Ids", "Stripped.Sequence", roles$protein))
    counts <- if (all(count_fields %in% names(persist$prototype)))
        newDIAFilteringAccumulator(persist$prototype, roles$protein) else NULL
    summary <- NULL
    sums <- c("input_precursor_rows", "unique_precursor_identities", "output_stripped_peptide_identities",
        "raw_partial_missing_outputs", "raw_all_missing_outputs", "normalised_partial_missing_outputs",
        "normalised_all_missing_outputs")
    reference <- saveWorkflowTable(store, NULL, basename(tempfile("dia-rollup-")), persist$identity,
        metadata$key_column, metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, store$progress)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            repeat {
                part <- cursor$next_group()
                if (is.null(part)) break
                attr(part, "workflow_source_order") <- NULL
                result <- rollUpPrecursorToPeptideHelper(part,
                    sample_id_column = roles$sample, protein_id_column = roles$protein,
                    peptide_sequence_column = roles$peptide, modified_peptide_sequence_column = roles$modified,
                    precursor_quantity_column = roles$raw, precursor_normalised_column = roles$normalised,
                    precursor_id_column = roles$precursor, precursor_charge_column = "Precursor.Charge",
                    core_utilisation = NA,
                    is_logged_data = FALSE, return_rollup_result = TRUE)
                if (!is.null(annotations)) result$data$Protein.Ids <-
                    annotations$Protein.Ids[match(result$data[[roles$protein]], annotations[[roles$protein]])]
                if (is.null(summary)) summary <<- result$summary else
                    for (field in sums) summary[[field]] <<- summary[[field]] + result$summary[[field]]
                for (start in seq.int(1L, max(1L, nrow(result$data)), by = 10000L)) {
                    at <- if (nrow(result$data)) seq.int(start, min(nrow(result$data), start + 9999L)) else integer()
                    if (!is.null(counts)) counts$push(result$data[at, , drop = FALSE])
                    appendWorkflowTableWriter(writer, result$data, at)
                }
                part <- result <- NULL
            }
            finishWorkflowTableWriter(writer)
        })
    cursor$close()
    summary$protein_ids_provenance <- if (length(input$annotation)) "preserved" else
        if (roles$protein == "Protein.Ids") "active_key" else "absent"
    list(reference = reference, summary = summary, columns = metadata$columns,
        dependencies = dependencies, counts = if (!is.null(counts)) counts$finish() else NULL)
}

#' Compute a precursor checkpoint without mutating workflow history.
#' @param manager Existing read-only component owner.
#' @param reference Frozen input checkpoint.
#' @param backend Finite legacy or qualified P14 implementation.
#' @param worker_cap Admitted worker memory in bytes, separately from the engine limit.
#' @return Bounded checkpoint/result descriptors for controller acceptance.
#' @noRd
computeDIAReferencePrecursorRollup <- function(manager, reference, backend = "P14", worker_cap = NULL) {
    if (!backend %in% c("legacy", "P14")) cli::cli_abort("unsupported rollup backend")
    # Managed actions pass the dispatched cap; only legacy unmanaged callers use the default.
    worker_cap_source <- if (is.null(worker_cap)) "default" else "dispatched"
    if (is.null(worker_cap)) worker_cap <- 12 * 1024^3
    before <- readDIAStateMetadata(manager, reference)
    checkpoint <- before$data
    checkpoint$dependencies <- collectDIANativeReferences(list(reference,
        checkpoint$dependencies, checkpoint$components, checkpoint$support, checkpoint$products))
    object <- readDIAExecutionComponents(manager, reference, fields = "design_matrix")
    if (isTRUE(object@is_logged_data)) cli::cli_abort("saved rollup requires linear peptide data")
    roles <- list(sample = object@sample_id, protein = object@protein_id_column,
        peptide = object@peptide_sequence_column, modified = "Modified.Sequence", precursor = "Precursor.Id",
        raw = object@raw_quantity_column, normalised = object@norm_quantity_column)
    columns <- names(object@peptide_data)
    input <- list(reference = checkpoint$components$peptide_data, roles = roles,
        evidence = intersect(c("identification_peptide_count", "identification_peptidoform_count",
            "peptides_for_protein_count", "peptidoforms_for_protein_count"), columns),
        annotation = if (roles$protein != "Protein.Ids" && "Protein.Ids" %in% columns) "Protein.Ids" else character())
    rollupFn <- rollUpPrecursorToPeptide
    template <- rollupFn(object)
    identity <- newSavedWorkflowScientificIdentity(input$reference$entity_id,
        checkpoint$components$design_matrix$entity_id, digest::digest(roles, algo = "sha256"),
        "precursor-linear-sum-v1", manager$.context$software)
    persist <- list(prototype = template@peptide_data, identity = identity)
    result <- if (identical(backend, "P14"))
        computeDIARelationalRollup(manager$.context$store, input, persist, arm = "F", profile = "P14",
            worker_cap = worker_cap,
            phase = getOption("multischolar.dia_rollup_phase", function(label, fn) fn())) else
        computeSavedWorkflowPrecursorRollup(manager$.context$store, input, persist)
    if (!is.null(result$execution$admission)) result$execution$admission$worker_cap_source <- worker_cap_source
    result$counts <- saveDIALargeFilteringCountsReference(manager$.context$store,
        result, object@protein_id_column)
    checkpoint$dependencies <- collectDIANativeReferences(list(checkpoint$dependencies, result$dependencies))
    checkpoint$slots$args <- template@args
    checkpoint$slots$args$rollUpPrecursorToPeptide$rollup_summary <- result$summary
    checkpoint$slots$raw_quantity_column <- "Peptide.RawQuantity"
    checkpoint$slots$norm_quantity_column <- "Peptide.Normalised"
    checkpoint$components$peptide_data <- result$reference
    checkpoint$dimensions <- list(rows = result$summary$output_stripped_peptide_identities, columns = result$columns)
    list(checkpoint = checkpoint, reference = result$reference, summary = result$summary,
        dependencies = checkpoint$dependencies, counts = result$counts, execution = result$execution,
        state_name = "precursor_rollup", description = "Applied precursor to peptide rollup",
        config = list(rollup_summary = result$summary), section = "rollUpPrecursorToPeptide")
}

#' Persist an explicit sample-filter candidate with existing saved-table machinery.
#' @param scope Configured physical input lease.
#' @param metadata Native result restoration metadata.
#' @param passing Retained sample labels, including explicit overrides.
#' @param roles Operation-declared sample and protein columns.
#' @param identity Frozen input/design/parameter/software identity.
#' @param output_mode Explicit physical or selection candidate.
#' @param batch_rows Pre-admitted row count bounded by expanded transfer bytes.
#' @return Accepted reference, producer counts and complete dependencies.
#' @noRd
persistDIAFilteredSamples <- function(scope, metadata, passing, roles,
    identity, output_mode, batch_rows) {
    sample <- roles$sample
    protein <- roles$protein
    prototype <- readSavedWorkflowTablePrototype(scope)
    prototype[[sample]] <- as.character(prototype[[sample]])
    counts <- newDIAFilteringAccumulator(prototype, protein)
    reader <- newWorkflowTableReadScope(scope$store, scope$reference, reader = "fastduckdb")
    reader$read_threads <- 8L
    on.exit(releaseWorkflowTableReadScope(reader), add = TRUE)
    fields <- if (output_mode == "physical") metadata$columns else
        unique(c(sample, "Run", "Protein.Ids", "Stripped.Sequence", protein))
    # Ordinal ranges preserve native source order. No wide engine result or sort
    # exists: the only result resident in R is this admitted read/filter batch.
    stream <- function(consume) {
        retained <- 0
        if (reader$metadata$rows > 0L) for (offset in seq.int(0L,
            reader$metadata$rows - 1L, by = batch_rows)) {
            page <- readSavedWorkflowTable(scope$store, scope$reference, fields,
                offset = offset, limit = min(batch_rows, reader$metadata$rows - offset),
                scope = reader)
            keep <- as.character(page$data[[sample]]) %in% passing
            data <- page$data[keep, , drop = FALSE]
            data[[sample]] <- as.character(data[[sample]])
            rownames(data) <- NULL
            counts$push(data)
            consume(data, as.double(page$source_order[keep]), retained)
            retained <- retained + nrow(data)
        }
        if (retained != metadata$rows)
            cli::cli_abort("sample filter output population differs", class = "saved_entity_corrupt")
    }
    if (output_mode == "physical") {
        output <- saveWorkflowTable(scope$store, NULL, basename(tempfile("sample-physical-")),
            identity, metadata$key_column, metadata = metadata,
            write_fn = function(path, restoration) {
                writer <- openWorkflowTableWriter(path, restoration, backend = "fast_r")
                on.exit(abortWorkflowTableWriter(writer), add = TRUE)
                appendWorkflowTableWriter(writer, prototype)
                stream(function(data, ordinals, retained) appendWorkflowTableWriter(writer, data))
                finishWorkflowTableWriter(writer)
            })
        dependencies <- list(scope$reference, output)
    } else {
        map_prototype <- data.frame(logical_ordinal = double(), source_ordinal = double())
        map_metadata <- workflowSavedTableMetadata(map_prototype, "logical_ordinal", 3L)
        map_metadata$rows <- metadata$rows
        map_identity <- identity
        map_identity$selectors <- paste0(identity$selectors, "-ordinals")
        map <- saveWorkflowTable(scope$store, NULL, basename(tempfile("sample-map-")),
            map_identity, "logical_ordinal", metadata = map_metadata,
            write_fn = function(path, restoration) {
                writer <- openWorkflowTableWriter(path, restoration, backend = "fast_r")
                on.exit(abortWorkflowTableWriter(writer), add = TRUE)
                appendWorkflowTableWriter(writer, map_prototype)
                stream(function(data, ordinals, retained) appendWorkflowTableWriter(writer,
                    data.frame(logical_ordinal = retained + as.double(seq_along(ordinals)),
                        source_ordinal = ordinals)))
                finishWorkflowTableWriter(writer)
            })
        columns <- data.frame(output = metadata$columns, source = metadata$columns,
            cast = rep("identity", length(metadata$columns)))
        if (!is.null(scope$metadata$column_metadata[[sample]]$levels))
            columns$cast[columns$source == sample] <- "factor_to_character"
        output <- saveWorkflowTableSelection(scope, map, columns, metadata,
            basename(tempfile("sample-selection-")), identity)
        dependencies <- list(scope$reference, map, output)
    }
    list(reference = output, counts = counts$finish(), dependencies = dependencies)
}

#' Admit and count complete samples without retaining full-population key state.
#' @param scope Existing configured physical read lease.
#' @param roles Declared sample/protein/peptide keys.
#' @param design Native design population and order.
#' @param output_mode Explicit physical or selection candidate.
#' @param worker_cap Actual admitted worker bytes.
#' @return Bounded counts and the recorded working-set admission.
#' @noRd
prepareDIASampleFilterCounts <- function(scope, roles, design, output_mode, worker_cap) {
    con <- scope$connection
    q <- \(x) as.character(DBI::dbQuoteIdentifier(con, x))
    table <- paste0("read_parquet(", DBI::dbQuoteString(con, scope$path), ")")
    keys <- unlist(roles, use.names = FALSE)
    fields <- if (output_mode == "physical") scope$metadata$columns else
        unique(c(roles$sample, "Run", "Protein.Ids", "Stripped.Sequence", roles$protein))
    storage <- workflowTableStorageColumns(scope$metadata, fields)
    widths <- vapply(storage, function(field) {
        spec <- scope$metadata$column_metadata[[field]]
        if (is.null(spec) || !is.null(spec$levels) || identical(spec$typeof, "character"))
            return(paste0("80 + coalesce(octet_length(encode(CAST(", q(field),
                " AS VARCHAR))),0)"))
        "16"
    }, character(1))
    key_width <- paste(vapply(keys, \(field) paste0(
        "80 + coalesce(octet_length(encode(CAST(", q(field), " AS VARCHAR))),0)"),
        character(1)), collapse = " + ")
    # Only per-sample scalar reductions precede admission. No wide sort/result
    # and no full-population distinct-key aggregation is materialized.
    counts <- collectWorkflowBulkMetadata(scope, paste0("SELECT CAST(", q(roles$sample),
        " AS VARCHAR) AS sample, count(*) AS rows, max(", key_width,
        ") AS key_width, max(96 + ", paste(widths, collapse = " + "),
        ") AS transfer_width FROM ", table, " GROUP BY ", q(roles$sample)))
    .validatePeptideRunDesign(stats::setNames(data.frame(counts$sample), roles$sample),
        design, roles$sample, roles$sample, caller = "computeDIAReferenceSampleFilter")
    order <- as.character(design[[roles$sample]])
    present <- order[order %in% counts$sample]
    width <- max(c(0, as.numeric(counts$transfer_width)))
    batch_rows <- max(1L, as.integer(min(10000, floor(33554432 / max(1, width)))))
    if (width > 67108864) cli::cli_abort("one sample-filter row exceeds its admitted transfer",
        class = "workflow_execution_unsupported")
    # Complete-sample batches (KIT-006): the fewest, at least three, whose per-batch identity state
    # (bounded by the largest samples any one batch can hold) fits the admitted cap; the fixed three
    # of the earlier policy are the floor, so the count at the normal cap stays three.
    state <- as.numeric(counts$rows) * (as.numeric(counts$key_width) + 128)
    estimate_for <- function(b) 768 * 1024^2 + 2 * completeGroupBatchBound(state, b) +
        2 * 67108864 + 268435456 + 2 * 8388608
    admitted <- tryCatch(chooseCompleteGroupBatches(seq.int(3L, max(3L, length(present))), estimate_for,
        worker_cap, "sample-adaptive-complete-batches-and-expanded-transfer-v4"),
        workflow_execution_unsupported = function(error) cli::cli_abort(
            "sample filter working-memory estimate exceeds admitted worker capacity",
            class = "workflow_execution_unsupported"))
    counts$identities <- numeric(nrow(counts))
    partition <- partitionBySamples(present, admitted$batches, stage = function(selected, index) {
        if (!length(selected)) return(NULL)
        part <- collectWorkflowBulkMetadata(scope, paste0("SELECT CAST(", q(roles$sample),
            " AS VARCHAR) AS sample, count(DISTINCT (", q(roles$protein), ",", q(roles$peptide),
            ")) AS identities FROM ", table, " WHERE ", q(roles$sample), " IN (",
            paste(DBI::dbQuoteString(con, selected), collapse = ","), ") GROUP BY ", q(roles$sample)))
        at <- match(part$sample, counts$sample)
        if (anyNA(at) || nrow(part) != length(selected) ||
            any(as.numeric(part$identities) > .Machine$integer.max))
            cli::cli_abort("sample identity count population differs", class = "saved_entity_corrupt")
        counts$identities[at] <<- as.numeric(part$identities)
        NULL
    })
    list(counts = counts, admission = list(estimated_bytes = admitted$estimated_bytes,
        worker_cap = worker_cap, max_row_bytes = width, batch_rows = batch_rows,
        batches = admitted$batches, count_batch_sizes = partition$batch_sizes,
        estimates = admitted$estimates, identity_groups = sum(counts$identities),
        policy = admitted$policy, limitation = admitted$limitation))
}

#' The minimum-peptides-per-sample filter as a declaration on the complete-sample partitioner
#' (KIT-006): per-sample identity counts in complete-sample batches after the global design
#' validation, one global decision, the unchanged native method resolving parameters and output
#' types. The kernel runs the reference entry point in its batched (`physical`) output mode.
#' @return Declaration from `newWorkflowOperation`.
#' @noRd
diaSampleFilterOperation <- function() {
    newWorkflowOperation(name = "minimum_peptides_per_sample", input_class = "PeptideQuantitativeData",
        fields = "design_matrix", unit = "complete_sample",
        kernel = function(object, parameters, store) computeDIAReferenceSampleFilter(store$manager, store$source,
            parameters[c("peptides_per_sample_cutoff", "inclusion_list")],
            worker_cap = parameters$worker_cap %||% 12 * 1024^3, output_mode = "physical"),
        output_class = "PeptideQuantitativeData", parameters = c("peptides_per_sample_cutoff", "inclusion_list"),
        product_type = "dia.sample_filter_result",
        global_statistics = c(design_validation_and_inclusion_absence =
            "unique design runs, input runs a subset of the design, and the inclusion-absence set, computed once over the whole input"),
        merge = c(design_matrix = "design", peptide_data = "concatenate", support_table = "concatenate",
            removal_ledger = "concatenate", counts = "sum"))
}

#' Execute the unchanged native method with the shared bounded persistence layer.
#' @param manager,reference,parameters Existing operation inputs.
#' @param worker_cap Actual admitted worker capacity.
#' @return Existing compute-only checkpoint/result contract.
#' @noRd
computeDIANativeSampleFilter <- function(manager, reference, parameters, worker_cap) {
    checkpoint <- readDIAStateMetadata(manager, reference)$data
    scope <- newWorkflowTableReadScope(manager$.context$store, checkpoint$components$peptide_data)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    settings <- configureWorkflowBulkQueryScope(scope)
    con <- scope$connection
    q <- \(name) as.character(DBI::dbQuoteIdentifier(con, name))
    table <- paste0("read_parquet(", DBI::dbQuoteString(con, scope$path), ")")
    strings <- scope$metadata$columns[vapply(scope$metadata$column_metadata,
        \(spec) identical(spec$typeof, "character"), logical(1))]
    dictionary <- sum(vapply(strings, function(field) {
        as.numeric(collectWorkflowBulkMetadata(scope, paste0(
            "SELECT coalesce(sum(80 + octet_length(encode(v))),0) AS n FROM (",
            "SELECT DISTINCT CAST(", q(field), " AS VARCHAR) AS v FROM ", table, ")"),
            max_rows = 1L, max_bytes = 65536)$n)
    }, numeric(1)))
    resident <- scope$metadata$rows * workflowProjectedRowBytes(scope$metadata,
        scope$metadata$columns) + dictionary + as.numeric(object.size(scope$metadata))
    estimate <- 768 * 1024^2 + 8 * resident + 2 * 67108864
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) ||
        worker_cap <= 0 || estimate > worker_cap)
        cli::cli_abort("native sample filter exceeds admitted worker capacity",
            class = "workflow_execution_unsupported")
    releaseWorkflowTableReadScope(scope)
    object <- readDIAExecutionComponents(manager, reference,
        fields = c("peptide_data", "design_matrix"))
    object <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    args <- object@args$filterMinNumPeptidesPerSample
    summary <- args$filter_summary
    checkpoint$support$filterMinNumPeptidesPerSample <- list(
        support_table = saveDIATableComponent(manager$.context, args$support_table),
        removal_ledger = saveDIATableComponent(manager$.context, args$removal_ledger))
    args$support_table <- args$removal_ledger <- NULL
    checkpoint$slots$args$filterMinNumPeptidesPerSample <- args
    counts <- newDIAFilteringAccumulator(object@peptide_data[0, , drop = FALSE],
        object@protein_id_column)
    output <- saveDIATableComponent(manager$.context, object@peptide_data,
        backend = "fast_r", consume = counts$push)
    produced <- list(reference = output, counts = counts$finish())
    checkpoint$components$peptide_data <- output
    checkpoint$dimensions <- list(rows = nrow(object@peptide_data), columns = names(object@peptide_data))
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, checkpoint$dependencies,
        checkpoint$components, checkpoint$support, checkpoint$products))
    list(checkpoint = checkpoint, reference = output, summary = summary,
        dependencies = checkpoint$dependencies,
        counts = saveDIALargeFilteringCountsReference(manager$.context$store, produced,
            object@protein_id_column),
        execution = list(strategy = "unchanged_native_shared_persistence", writer_backend = "fast_r",
            bulk_query = list(settings = settings),
            admission = list(estimated_bytes = estimate, worker_cap = worker_cap,
                resident_bytes = resident, policy = "native-resident-copy-slack-v1")),
        state_name = "sample_filtered", description = "Applied sample quality filter",
        config = parameters, section = "filterMinNumPeptidesPerSample")
}

#' Compute sample identity filtering without advancing a workflow head.
#' @param manager Read-only component owner.
#' @param reference Frozen input checkpoint.
#' @param parameters Explicit cutoff and inclusion list.
#' @param worker_cap Admitted worker memory in bytes.
#' @param output_mode Explicit experimental physical/selection output, or legacy default.
#' @return Bounded checkpoint, summary and ledger references for the controller.
#' @noRd
computeDIAReferenceSampleFilter <- function(manager, reference, parameters, worker_cap = 12 * 1024^3,
    output_mode = "legacy") {
    if (!is.character(output_mode) || length(output_mode) != 1L || is.na(output_mode) ||
        !output_mode %in% c("legacy", "native", "physical", "selection"))
        cli::cli_abort("unsupported sample filter output mode", class = "workflow_execution_unsupported")
    if (output_mode == "native")
        return(computeDIANativeSampleFilter(manager, reference, parameters, worker_cap))
    checkpoint <- readDIAStateMetadata(manager, reference)$data
    object <- readDIAExecutionComponents(manager, reference, fields = "design_matrix")
    sample <- object@sample_id
    protein <- object@protein_id_column
    peptide <- object@peptide_sequence_column
    # Resolve parameters and output types through the unchanged native method on zero rows.
    template <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    args <- template@args$filterMinNumPeptidesPerSample
    threshold <- args$filter_summary$peptides_per_sample_cutoff
    inclusion <- args$filter_summary$inclusion_list
    source <- checkpoint$components$peptide_data
    scope <- newWorkflowTableReadScope(manager$.context$store, source)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    assertWorkflowPhysicalTableScope(scope)
    con <- scope$connection
    bulk_settings <- configureWorkflowBulkQueryScope(scope)
    q <- function(name) as.character(DBI::dbQuoteIdentifier(con, name))
    table <- paste0("read_parquet(", DBI::dbQuoteString(con, scope$path), ")")
    keys <- c(sample, protein, peptide)
    if (!all(keys %in% scope$metadata$columns)) cli::cli_abort("missing sample filter keys")
    if (!all(vapply(scope$metadata$column_metadata[keys], function(spec)
        identical(spec$typeof, "character") || !is.null(spec$levels), logical(1))))
        cli::cli_abort("sample filtering requires character or factor identity keys",
            class = "workflow_execution_unsupported")
    bad <- paste0("(", q(keys), " IS NULL OR regexp_full_match(CAST(", q(keys),
        " AS VARCHAR), '[ \\t\\r\\n]*'))")
    if (nrow(DBI::dbGetQuery(con, paste("SELECT 1 FROM", table, "WHERE",
        paste(bad, collapse = " OR "), "LIMIT 1"))))
        cli::cli_abort("sample/protein/peptide identity is missing")
    if (output_mode != "legacy") {
        prepared <- prepareDIASampleFilterCounts(scope,
            list(sample = sample, protein = protein, peptide = peptide),
            object@design_matrix, output_mode, worker_cap)
        counts <- prepared$counts
        admission <- prepared$admission
    } else {
    key_width <- paste(vapply(keys, \(field) paste0(
        "coalesce(octet_length(encode(CAST(", q(field), " AS VARCHAR))),0) + 80"),
        character(1)), collapse = " + ")
    projected_width <- if (output_mode == "legacy") "" else
        paste0(", max(", key_width, ") AS key_width")
    counts <- collectWorkflowBulkMetadata(scope, paste0("SELECT CAST(", q(sample), " AS VARCHAR) AS sample, ",
        "count(DISTINCT (", q(protein), ",", q(peptide), ")) AS identities, count(*) AS rows ",
        projected_width, " FROM ", table, " GROUP BY ", q(sample)), max_rows = 10000L,
        max_bytes = 8388608)
    width <- if (output_mode == "legacy") {
        widths <- vapply(scope$metadata$columns, function(field) {
        spec <- scope$metadata$column_metadata[[field]]
        if (!is.null(spec$levels) || identical(spec$typeof, "character")) return(paste0(
            "coalesce(octet_length(encode(CAST(", q(field), " AS VARCHAR))),0) + 16"))
        "8"
    }, character(1))
        as.numeric(collectWorkflowBulkMetadata(scope, paste0("SELECT coalesce(max(",
        paste(widths, collapse = " + "), "),0) AS n FROM ", table),
        max_rows = 1L, max_bytes = 65536)$n)
    } else max(c(0, as.numeric(counts$key_width)))
    if (width > 67108864) cli::cli_abort("one filter row exceeds the bounded codec transfer",
        class = "workflow_execution_unsupported")
    transfer <- if (output_mode == "legacy")
        2 * min(scope$metadata$rows, 10000) * width else 2 * 67108864
    estimate <- 768 * 1024^2 + transfer +
        2 * sum(as.numeric(counts$identities)) * (width + 128) + 2 * 8388608
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) ||
        worker_cap <= 0 || estimate > worker_cap) cli::cli_abort(
        "sample filter working-memory estimate exceeds admitted worker capacity",
        class = "workflow_execution_unsupported")
    admission <- list(estimated_bytes = estimate, worker_cap = worker_cap,
        max_row_bytes = width, identity_groups = sum(as.numeric(counts$identities)),
        policy = if (output_mode == "legacy") "sample-identity-state-and-bounded-transfer-v1" else
            "sample-projected-identity-and-codec-transfer-v2",
        limitation = "admission estimate; kernel cgroup guard remains authoritative")
    }
    observed <- stats::setNames(data.frame(counts$sample), sample)
    .validatePeptideRunDesign(observed, object@design_matrix, sample, sample,
        caller = "computeDIAReferenceSampleFilter")
    support <- args$support_table[rep(NA_integer_, nrow(counts)), , drop = FALSE]
    rownames(support) <- NULL
    names(support)[names(support) == sample] <- "sample"
    support$sample <- counts$sample
    support$distinct_peptide_identity_count <- as.integer(counts$identities)
    support$inclusion_override <- support$sample %in% inclusion
    support$sample_passes <- support$distinct_peptide_identity_count >= threshold | support$inclusion_override
    support$decision_reason <- rep("below_distinct_identity_threshold", nrow(support))
    support$decision_reason[support$inclusion_override] <- "explicit_inclusion_override"
    support$decision_reason[support$distinct_peptide_identity_count >= threshold] <-
        "meets_distinct_identity_threshold"
    support <- dplyr::arrange(support, .data$sample)
    names(support)[1L] <- sample
    removal <- support[!support$sample_passes, c(sample, "distinct_peptide_identity_count")]
    removal$required_distinct_peptide_identities <- rep(threshold, nrow(removal))
    removal$failure_reason <- support$decision_reason[!support$sample_passes]
    rownames(removal) <- NULL
    passing <- support[[sample]][support$sample_passes]
    summary <- args$filter_summary
    summary$input_sample_count <- nrow(counts)
    summary$retained_sample_count <- sum(support$sample_passes)
    summary$removed_sample_count <- sum(!support$sample_passes)
    summary$inclusion_entries_absent_from_data <- setdiff(inclusion, counts$sample)
    identity <- newSavedWorkflowScientificIdentity(source$entity_id,
        checkpoint$components$design_matrix$entity_id, digest::digest(parameters, algo = "sha256"),
        paste0("distinct-pair-sample-filter-", output_mode, "-v1"), manager$.context$software)
    if (output_mode == "legacy") identity$selectors <- "distinct-pair-sample-filter-v1"
    metadata <- workflowSavedTableMetadata(template@peptide_data, ".rownames", 3L)
    metadata$rows <- as.integer(sum(counts$rows[counts$sample %in% passing]))
    produced <- if (output_mode != "legacy") persistDIAFilteredSamples(scope, metadata,
        passing, list(sample = sample, protein = protein), identity, output_mode,
        admission$batch_rows) else NULL
    output <- if (!is.null(produced)) produced$reference else
        saveWorkflowTable(scope$store, NULL, basename(tempfile("wfx-sample-filter-")), identity,
        ".rownames", metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            streamSavedWorkflowTable(scope, scope$metadata$columns, function(data, offset, row_names) {
                keep <- as.character(data[[sample]]) %in% passing
                data <- data[keep, , drop = FALSE]
                data[[sample]] <- as.character(data[[sample]])
                rownames(data) <- NULL
                appendWorkflowTableWriter(writer, data)
            })
            if (is.null(writer$output)) appendWorkflowTableWriter(writer, template@peptide_data)
            finishWorkflowTableWriter(writer)
        })
    checkpoint$dependencies <- collectDIANativeReferences(list(reference, produced$dependencies,
        checkpoint$dependencies, checkpoint$components, checkpoint$support, checkpoint$products))
    checkpoint$components$peptide_data <- output
    checkpoint$dimensions <- list(rows = metadata$rows, columns = metadata$columns)
    args$support_table <- args$removal_ledger <- NULL
    args$filter_summary <- summary
    checkpoint$slots$args$filterMinNumPeptidesPerSample <- args
    checkpoint$support$filterMinNumPeptidesPerSample <- list(
        support_table = saveDIATableComponent(manager$.context, support),
        removal_ledger = saveDIATableComponent(manager$.context, removal))
    list(checkpoint = checkpoint, reference = output, summary = summary,
        dependencies = checkpoint$dependencies, counts = if (is.null(produced)) NULL else
            saveDIALargeFilteringCountsReference(scope$store, produced, protein),
        execution = list(admission = admission, output_mode = output_mode,
            writer_backend = if (output_mode == "legacy") "arrow" else "fast_r",
            output_transport = "bounded_R",
            bulk_query = list(settings = bulk_settings,
            global_statistics = if (output_mode == "legacy") "single_complete_population_pass" else
                "global_validation_before_complete_sample_counts",
            independence_unit = if (output_mode == "legacy") "sample_decision_no_batching" else
                "complete_sample_after_global_validation",
            declared_unit = if (output_mode == "legacy") NULL else diaSampleFilterOperation()$unit)),
        state_name = "sample_filtered", description = "Applied sample quality filter",
        config = parameters, section = "filterMinNumPeptidesPerSample")
}
