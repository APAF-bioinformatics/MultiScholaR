#' Read DA and Summary prerequisites without constructing a saved S4 state.
#' @param manager Current state owner.
#' @param state_name Optional existing state.
#' @return Bounded slot metadata and design, plus a reference on the component path.
#' @noRd
readDIAStateGuiMetadata <- function(manager, state_name = NULL) {
    if (!usesDIAReferenceStates(manager)) {
        object <- manager$getState(state_name)
        if (is.null(object)) return(NULL)
        return(c(stats::setNames(lapply(methods::slotNames(object),
            \(name) methods::slot(object, name)), methods::slotNames(object)), list(legacy_object = object)))
    }
    reference <- getDIAStateReference(manager, state_name)
    if (is.null(reference)) return(NULL)
    spec <- readDIAStateMetadata(manager, reference)$data
    if (!is.null(spec$slots$args$daAnalysisUI$formula_string)) {
        spec$slots$args$deAnalysisParameters$formula_string <- spec$slots$args$daAnalysisUI$formula_string
    }
    c(spec$slots, list(reference = reference, dimensions = spec$dimensions,
        design_matrix = materializeDIATable(manager$.context, spec$components$design_matrix)))
}

#' Obtain a protein key from a GUI reference without loading quantities.
#' @param value Existing S4 or component reference.
#' @return Declared protein identifier column.
#' @noRd
diaGuiProteinKey <- function(value) {
    if (isSavedWorkflowReference(value)) readDIAStateMetadata(getDIAActiveManager(), value)$data$slots$protein_id_column else
        methods::slot(value, "protein_id_column")
}

#' Invalidate GUI result references when the selected scientific population changes.
#' @param da_data,workflow Existing module owners.
#' @return NULL invisibly; immutable cached products remain reusable for their own inputs.
#' @noRd
invalidateDIAVisibleResults <- function(da_data, workflow) {
    manager <- workflow$state_manager
    if (!usesDIAReferenceStates(manager)) return(invisible(NULL))
    # KIT-053: the DIA_LIMPA route's final state is ruv_covariates; without it here the state-update observer that
    # fires after the DA cleared the accepted results and the status text the browser acknowledges never held.
    valid_state <- manager$current_state %in% c("ruv_covariates", "normalized", "ruv_corrected", "correlation_filtered")
    result <- da_data$da_results_list
    changed <- !valid_state
    if (!changed && isDIAResultManifest(result)) {
        changed <- !identical(result$source$owner, manager$.context$store$identity)
        if (!changed) {
            fields <- c("protein_quant_table", "design_matrix", "protein_id_table")
            changed <- !identical(readDIAStateMetadata(manager)$data$components[fields],
                readDIAStateMetadata(manager, result$source)$data$components[fields])
        }
    }
    if (changed) {
        da_data$da_results_list <- NULL
        da_data$analysis_complete <- FALSE
        da_data$current_heatmap_plot <- NULL
        da_data$current_row_clusters <- da_data$current_col_clusters <- NULL
        workflow$da_analysis_results_list <- NULL
        if (!valid_state) da_data$current_s4_object <- NULL
    }
    invisible(NULL)
}

#' Reuse the exact uncontrasted native limma fit inside an explicit DA action.
#' @param data,design_m,design,samples Existing runTestsContrasts execution inputs.
#' @param formula,weights Existing model identity, including the scalar NA policy.
#' @return Existing base-fit artifact, never a replacement numerical implementation.
#' @noRd
fitDIAExecutionBaseModel <- function(data, design_m, design, samples, formula, weights,
    group_id = "group", technical_replicate_id = "replicates") {
    manager <- getDIAActiveManager()
    if (!usesDIAReferenceStates(manager) || !isTRUE(manager$.context$da_reference_execution)) {
        return(fitProtDaBaseModel(data, design_m, design, samples, group_id = group_id,
            technical_replicate_id = technical_replicate_id))
    }
    selected_design <- design[samples, , drop = FALSE]
    # KIT-035: the block from the declaration (the same string as the group_replicates block it replaces, so the
    # base-fit identity of an existing store is unchanged).
    replicate <- resolveProtTechnicalReplicateBlock(selected_design, group_id, technical_replicate_id)
    block <- replicate$block
    identity <- newWorkflowDiaBaseFitIdentity(data, selected_design, formula, block, weights,
        "accepted-normalized-component", "accepted-ruv-component", manager$.context$software)
    found <- findSavedWorkflowEntities(manager$.context$store, list(kind = "da_fit", identity = identity), max_matches = 2L)
    if (identical(found$status, "multiple")) cli::cli_abort("ambiguous DA base-fit identity")
    hit <- if (identical(found$status, "hit")) resolveSavedWorkflowEntity(manager$.context$store,
        "da_fit", identity, factors = found$matches[[1L]]$metadata$factors) else list(status = "miss")
    if (identical(hit$status, "hit")) {
        value <- hit$value
        manager$.context$da_base_fit_reference <- hit$reference
        return(list(schema = "prot_da_base_fit_v1", fit = value$fit,
            feature_ids = value$feature_ids, sample_ids = value$sample_ids,
            design_columns = value$design_columns, block = block,
            group_id = group_id, technical_replicate_id = technical_replicate_id, block_reason = replicate$reason,
            uses_duplicate_correlation = isTRUE(replicate$technical),
            consensus_correlation = hit$reference$metadata$factors$consensus_correlation))
    }
    fit <- fitProtDaBaseModel(data, design_m, design, samples, group_id = group_id,
        technical_replicate_id = technical_replicate_id)
    reference <- saveWorkflowDiaBaseFit(manager$.context$store, fit,
        basename(tempfile("dia-base-fit-")), identity)
    rootDIAProduct(manager, reference, list(manager$.context$da_source_reference))
    manager$.context$da_base_fit_reference <- reference
    fit
}

#' Identify the bounded reference manifest delivered to DA consumers.
#' @param value GUI-held result manifest.
#' @return Logical scalar, without resolving a payload.
#' @noRd
isDIAResultManifest <- function(value) {
    is.list(value) && identical(value[["schema", exact = TRUE]], "dia.da_results.v1")
}

#' Select an existing contrast reference by its original or friendly label.
#' @param manifest Bounded result manifest.
#' @param contrast Selected UI contrast.
#' @return One immutable contrast manifest.
#' @noRd
selectDIAContrastReference <- function(manifest, contrast) {
    labels <- manifest$contrasts
    selected <- which(labels$contrasts == contrast | labels$friendly_names == contrast |
        labels$full_format == contrast | labels$friendly_names == sub("=.*$", "", contrast))
    if (length(selected) != 1L) cli::cli_abort("contrast is absent or ambiguous", class = "saved_entity_invalid")
    manifest$results[[labels$contrasts[[selected]]]]
}

#' Query table pages and counts from already adjusted full-population results.
#' @param manifest,contrast Accepted results and selected contrast.
#' @param significance,qvalue,lfc Existing display filters.
#' @param offset,limit Bounded page coordinates.
#' @return Display data and full-population counts; statistics are never recomputed.
#' @noRd
readDIAResultPage <- function(manifest, contrast, significance = "all", qvalue = 0.05,
    lfc = 0, offset = 0L, limit = 100L, manager = getDIAActiveManager()) {
    # KIT-017 screen run: the widget's page requests carry no reactive domain; the render passes the owner.
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("no state owner for the result page", class = "dia_owner_missing")
    result <- selectDIAContrastReference(manifest, contrast)
    reference <- result$tables$da_proteins_long
    if (length(limit) != 1L || !is.finite(limit) || limit < 1 || limit > 10000 || limit != floor(limit) ||
        length(offset) != 1L || !is.finite(offset) || offset < 0 || offset != floor(offset) ||
        !significance %in% c("all", "significant", "up", "down")) cli::cli_abort("invalid result page")
    view <- resolveDIAResultSelection(manager, reference, qvalue, lfc)
    columns <- c(manifest$protein_id_column, "log2FC", "raw_pvalue", "fdr_qvalue")
    if (significance == "all") {
        page <- readSavedWorkflowTable(manager$.context$store, reference, columns,
            offset = offset, limit = limit)$data
    } else {
        reader <- if (requireNamespace("fastduckdb", quietly = TRUE)) "fastduckdb" else "duckdb"
        scope <- newWorkflowTableReadScope(manager$.context$store, view$rows, reader = reader)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        selected <- integer()
        skipped <- 0L
        at <- 0L
        fields <- if (significance == "significant") "source_order" else c("source_order", significance)
        while (at < scope$metadata$rows && length(selected) < limit) {
            count <- min(10000L, scope$metadata$rows - at)
            part <- readSavedWorkflowTable(scope$store, view$rows, fields,
                offset = at, limit = count, scope = scope, reader = reader)$data
            matches <- if (significance == "significant") part$source_order else
                part$source_order[part[[significance]]]
            remaining <- max(0, offset - skipped)
            skipped <- skipped + length(matches)
            if (remaining < length(matches)) selected <- c(selected,
                utils::head(matches[seq.int(remaining + 1L, length(matches))], limit - length(selected)))
            at <- at + count
        }
        releaseWorkflowTableReadScope(scope)
        page <- readWorkflowSelectedOrdinals(manager$.context$store, reference, columns, selected)
    }
    list(data = page, counts = view$counts, full_population_rows = reference$metadata$shape[[1L]])
}

#' Provide one selected DA display with an explicit expression boundary.
#' @param manifest,contrast Accepted reference collection and selected contrast.
#' @param expression TRUE only for Glimma's linked expression display.
#' @return Projected results and, when required, a transient expression view.
#' @noRd
readDIAResultDisplay <- function(manifest, contrast, expression = FALSE) {
    manager <- getDIAActiveManager()
    selected <- selectDIAContrastReference(manifest, contrast)
    columns <- intersect(unique(c(manifest$protein_id_column, "comparison", "log2FC",
        "raw_pvalue", "fdr_qvalue", "gene_name", "gene_names", "Protein.Ids", "Entry", "description")),
        selected$columns)
    result <- readDIAContrastOutput(selected, "da_proteins_long", columns)
    if (expression) result$theObject <- readDIAExecutionComponents(manager, manifest$source,
        fields = c("protein_quant_table", "design_matrix"))
    result
}

#' Read the selected heatmap matrix and reuse clustering across render changes.
#' @param manifest,contrast Accepted DA results and current contrast.
#' @param parameters Existing selection, scaling and clustering settings.
#' @param color_scheme,show_gene_names Render-only choices.
#' @return Original heatmap handoff from a bounded saved scientific product.
#' @noRd
generateDIAReferenceHeatmap <- function(manifest, contrast, parameters, color_scheme, show_gene_names) {
    manager <- getDIAActiveManager()
    selected <- selectDIAContrastReference(manifest, contrast)
    identity <- diaComponentScientificIdentity(manager, manifest$source, "DA-heatmap-clustering",
        c(parameters, list(result = selected$tables$da_proteins_long$entity_id)))
    hit <- resolveSavedWorkflowEntity(manager$.context$store, "scientific_product", identity)
    if (identical(hit$status, "hit")) return(renderProtDAHeatmapProduct(hit$value$value,
        color_scheme, show_gene_names))
    result_scope <- newWorkflowTableReadScope(manager$.context$store, selected$tables$da_proteins_long)
    on.exit(releaseWorkflowTableReadScope(result_scope), add = TRUE)
    assertWorkflowPhysicalTableScope(result_scope)
    con <- result_scope$connection
    q <- \(value) as.character(DBI::dbQuoteIdentifier(con, value))
    columns <- intersect(unique(c(manifest$protein_id_column, "comparison", parameters$qvalue_column,
        parameters$log2fc_column, "gene_names")), selected$columns)
    if (!is.numeric(parameters$top_n_genes) || length(parameters$top_n_genes) != 1L ||
        !is.finite(parameters$top_n_genes) || parameters$top_n_genes < 1L || parameters$top_n_genes > 10000L) {
        cli::cli_abort("heatmap selection exceeds read limits", class = "saved_entity_limit")
    }
    fc <- q(parameters$log2fc_column)
    top <- DBI::dbGetQuery(con, paste0("SELECT ", paste(q(workflowTableStorageColumns(result_scope$metadata, columns)), collapse = ", "),
        " FROM read_parquet(?) WHERE ", q(parameters$qvalue_column), " < ? ORDER BY ",
        "CASE WHEN isnan(", fc, ") THEN NULL ELSE abs(", fc, ") END DESC NULLS LAST, ",
        '".multischolar_order" LIMIT ?'), params = list(result_scope$path,
            parameters$da_q_val_thresh, as.integer(parameters$top_n_genes)))
    if (as.numeric(object.size(top)) > 67108864) cli::cli_abort("heatmap selection exceeds byte limits")
    top <- restoreSavedWorkflowTable(top, result_scope$metadata, columns)
    releaseWorkflowTableReadScope(result_scope)
    if (!nrow(top)) return(NULL)
    spec <- readDIAStateMetadata(manager, manifest$source)$data
    expression_scope <- newWorkflowTableReadScope(manager$.context$store, spec$components$protein_quant_table)
    on.exit(releaseWorkflowTableReadScope(expression_scope), add = TRUE)
    assertWorkflowPhysicalTableScope(expression_scope)
    DBI::dbExecute(expression_scope$connection, paste0("CREATE TEMP VIEW dia_heatmap_expression AS SELECT * FROM read_parquet(",
        DBI::dbQuoteString(expression_scope$connection, expression_scope$path), ")"))
    quant <- materializeDIATable(manager$.context, spec$components$protein_quant_table,
        scope = expression_scope, relation = "dia_heatmap_expression", filter_column = manifest$protein_id_column,
        filter_value = unique(as.character(top[[manifest$protein_id_column]])))
    attr(quant, "dia_source_order") <- NULL
    releaseWorkflowTableReadScope(expression_scope)
    object <- readDIAExecutionComponents(manager, manifest$source, fields = "design_matrix")
    object@protein_quant_table <- quant
    product <- do.call(computeProtDAHeatmapProduct, c(list(da_results_list = list(da_proteins_long = top,
        theObject = object), selected_contrast = contrast), parameters))
    object <- quant <- NULL
    if (is.null(product)) return(NULL)
    reference <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowScientificProduct(product, "dia-heatmap-clustering"),
        basename(tempfile("dia-heatmap-")), identity, "scientific_product")
    rootDIAProduct(manager, reference, list(manifest$source, selected$tables$da_proteins_long))
    renderProtDAHeatmapProduct(product, color_scheme, show_gene_names)
}

#' Resolve one contrast's saved tables or display products for an output consumer.
#' @param result Bounded contrast manifest or native result list.
#' @param fields Required result fields.
#' @param columns Optional projected long-result columns.
#' @return Only the requested values; the caller owns their temporary lifetime.
#' @noRd
readDIAContrastOutput <- function(result, fields = c("da_proteins_long", "num_sig_da_molecules_first_go"),
    columns = NULL) {
    if (!is.list(result) || !identical(result$schema, "dia.contrast_reference.v1")) return(result)
    manager <- getDIAActiveManager()
    values <- list()
    for (field in fields) {
        reference <- result$tables[[field]] %||% result$products[[field]]
        if (is.null(reference)) next
        values[[field]] <- if (isWorkflowTableReference(reference))
            materializeDIATable(manager$.context, reference,
                if (field == "da_proteins_long") columns else NULL) else
            hydrateSavedWorkflowEntity(manager$.context$store, reference)$value
    }
    values
}

#' The per-contrast DA operation (KIT-004): the unchanged `differentialAbundanceAnalysis` kernel on one
#' projected object; tables become components and display values derived products, assembled into
#' the existing `dia.contrast_reference.v1` record saved under the existing identity and factors.
#' @noRd
#' The route factor of the DA identities (KIT-053): `limpa` on the `DIA_LIMPA` type, `iq` otherwise (the constant
#' every DIA, LFQ and TMT identity carried before).
#' @noRd
diaActiveRoute <- function(manager = getDIAActiveManager()) {
    type <- tryCatch(manager$workflow_type, error = function(e) NULL)
    if (identical(type, "DIA_LIMPA")) "limpa" else "iq"
}

diaContrastOperation <- function(route = "iq") {
    newWorkflowOperation("DA-contrast-output", "ProteinQuantitativeData",
        fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
        product_type = "dia-contrast-references", factors = list(route = route, layer = "contrast_references"),
        parameters = c("formula", "qvalue", "lfc", "weights", "contrast", "annotations"),
        kernel = function(object, parameters, store) {
            single <- data.frame(contrasts = parameters$contrast, stringsAsFactors = FALSE)
            arguments <- list(contrasts_tbl = single, formula_string = parameters$formula,
                da_q_val_thresh = parameters$qvalue, treat_lfc_cutoff = parameters$lfc,
                eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
                args_row_id = object@protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
            for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) {
                if (is.null(object@args[[section]])) object@args[[section]] <- arguments
            }
            result <- differentialAbundanceAnalysis(object, contrasts_tbl = single,
                formula_string = parameters$formula, da_q_val_thresh = parameters$qvalue,
                treat_lfc_cutoff = parameters$lfc, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
            object <- NULL
            accepted <- list(schema = "dia.contrast_reference.v1", source = store$source,
                contrast = parameters$name, comparison = parameters$comparison,
                tables = list(), products = list(), qvalue_warnings = result$qvalue_warnings)
            for (field in c("da_proteins_long", "da_proteins_wide", "contrasts_results_table",
                "significant_rows", "norm_counts")) {
                if (is.data.frame(result[[field]])) accepted$tables[[field]] <- store$table(result[[field]])
            }
            for (field in c("num_sig_da_molecules_first_go", "num_sig_da_molecules", "qvalue_warnings")) {
                if (is.null(result[[field]])) next
                value <- result[[field]]
                if (field == "num_sig_da_molecules_first_go" ||
                    field == "num_sig_da_molecules" && is.list(value) && !is.data.frame(value))
                    value <- list(table = value$table)
                if (is.list(value)) for (entry in names(value)) {
                    if (inherits(value[[entry]], c("gg", "ggplot", "patchwork")))
                        value[[entry]] <- stripProtNormQcPlotOwners(value[[entry]])
                }
                accepted$products[[field]] <- store$product(value, field, paste0("dia-", field))
            }
            accepted$columns <- names(result$da_proteins_long)
            accepted$rows <- nrow(result$da_proteins_long)
            accepted$model <- result$da_model # KIT-053: limma or dpcDE, the formula used, the W terms, the block
            accepted
        })
}

#' Run native DA one contrast at a time through the DA declaration and retain only accepted references.
#' @param settings Plain list: `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @param da_data,workflow,paths Existing module owners and publication paths.
#' @param contrasts Contrasts table; NULL uses `workflow$contrasts_tbl`, then the available contrasts.
#' @param annotations Bounded annotation table for the published outputs; NULL uses
#'   `workflow$uniprot_dat_cln`, then an empty table.
#' @return Reference manifest after required outputs have completed; the module edge updates the
#'   selectors and notifies.
#' @noRd
runDIAReferenceDaAnalysis <- function(settings, da_data, workflow, paths, contrasts = NULL, annotations = NULL,
    publish = TRUE) {
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "da")
        source <- da_data$current_s4_object
        if (!isSavedWorkflowReference(source)) cli::cli_abort("load a reference-backed protein state before DA")
        if (is.null(contrasts)) contrasts <- workflow$contrasts_tbl
        if (is.null(contrasts)) contrasts <- data.frame(contrasts = da_data$contrasts_available)
        contrasts <- normaliseProtDaContrastsTable(contrasts)
        if (is.null(annotations)) annotations <- workflow$uniprot_dat_cln
        if (is.null(annotations)) annotations <- data.frame()
        spec <- readDIAStateMetadata(manager, source)$data
        prototype <- readDIAExecutionComponents(manager, source, fields = "design_matrix")
        validateProtDaModelAndContrasts(prototype, settings$formula_string, contrasts)
        prototype <- NULL
        manifest <- list(schema = "dia.da_results.v1", source = source,
            protein_id_column = spec$slots$protein_id_column, sample_id = spec$slots$sample_id,
            contrasts = contrasts, results = list(), parameters = list(formula = settings$formula_string,
                qvalue = settings$da_q_val_thresh, lfc = settings$treat_lfc_cutoff, weights = NA))
        # KIT-022: the DA UI arguments are committed onto the current state before the contrast loop, so the
        # manifest's source is the head state the module and the worker see (one reference per scientific state).
        # The contrast identity strips these arguments, so the contrast products are unchanged by the commit.
        workflow$da_ui_params <- list(q_value_threshold = settings$da_q_val_thresh,
            log_fold_change_cutoff = settings$treat_lfc_cutoff, treat_enabled = settings$treat_lfc_cutoff > 0,
            formula_string = settings$formula_string, timestamp = Sys.time())
        state_record <- readDIAStateMetadata(manager)
        checkpoint <- state_record$data
        da_parameters <- workflow$da_ui_params
        recorded <- checkpoint$slots$args$daAnalysisUI
        if (!identical(recorded[setdiff(names(recorded), "timestamp")],
            da_parameters[setdiff(names(da_parameters), c("timestamp", "treat_enabled"))])) {
            first <- data.frame(contrasts = contrasts$full_format[[1L]], stringsAsFactors = FALSE)
            defaults <- list(contrasts_tbl = first, formula_string = settings$formula_string,
                da_q_val_thresh = settings$da_q_val_thresh, treat_lfc_cutoff = settings$treat_lfc_cutoff,
                eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
                args_row_id = manifest$protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
            if (is.null(checkpoint$slots$args$differentialAbundanceAnalysis))
                checkpoint$slots$args$differentialAbundanceAnalysis <- defaults
            checkpoint$slots$args$daAnalysisUI <- da_parameters[setdiff(names(da_parameters), "treat_enabled")]
            if (nrow(contrasts) > 1L && is.null(checkpoint$slots$args$differentialAbundanceAnalysisHelper))
                checkpoint$slots$args$differentialAbundanceAnalysisHelper <- defaults
            commitDIAReferenceState(manager, manager$current_state, checkpoint,
                state_record$config, state_record$description)
            source <- getDIAStateReference(manager)
            manifest$source <- source
            # The module's holder follows the head (a reference, KIT-012), so a rerun starts from the head too.
            da_data$current_s4_object <- source
        }
        context <- manager$.context
        context$da_reference_execution <- TRUE
        context$da_source_reference <- source
        on.exit({ context$da_reference_execution <- FALSE; context$da_source_reference <- NULL }, add = TRUE)
        operation <- diaContrastOperation(diaActiveRoute(manager))
        for (i in seq_len(nrow(contrasts))) {
            name <- contrasts$contrasts[[i]]
            run <- runSavedOperation(manager, operation, source, c(manifest$parameters,
                list(contrast = contrasts$full_format[[i]], annotations = spec$components$protein_id_table$entity_id,
                    name = name, comparison = contrasts$friendly_names[[i]])))
            if (identical(run$status, "computed")) {
                snapshot <- context$snapshot
                snapshot$da_dependencies[[run$reference$entity_id]] <- collectDIANativeReferences(run$value)
                publishDIAStateHead(context, snapshot)
            }
            manifest$results[[name]] <- run$value
        }
        clearDIAStepProvenance(manager) # KIT-025 audit PH-01: the contrasts commit no state
        snapshot <- context$snapshot
        # KIT-035: the model that ran, from the base fit's recorded factors (provenance; not part of any identity).
        fit_factors <- manager$.context$da_base_fit_reference$metadata$factors
        # KIT-053: the kernel's own account of the model (limma or dpcDE with the W terms) first; the limma base fit's
        # factors beside it (NULL on the dpcDE branch, which fits no saved base fit).
        kernel_model <- manifest$results[[1L]]$model %||% list(model = "limma")
        manifest$model <- c(kernel_model, list(route = diaActiveRoute(manager),
            group_id = fit_factors$group_id, technical_replicate_id = fit_factors$technical_replicate_id,
            block_reason = kernel_model$block_reason %||% fit_factors$block_reason,
            uses_duplicate_correlation = kernel_model$uses_duplicate_correlation %||% fit_factors$uses_duplicate_correlation,
            consensus_correlation = fit_factors$consensus_correlation))
        manifest_reference <- saveDIAStateValue(context, manifest, "da-result-references", codec = "rds")
        snapshot$da_results <- manifest_reference
        snapshot$da_dependencies[[manifest_reference$entity_id]] <- collectDIANativeReferences(
            list(manifest$source, manifest$results))
        publishDIAStateHead(context, snapshot)
        # KIT-019: under a managed application the DA module requests the publication from the worker
        # (`da_publication`); in process the publication runs here as before.
        if (publish) publishDIAActionOutputs(manager, getDIAStateReference(manager), "da",
            paths[c("da_output_dir", "publication_graphs_dir")], function(staged)
                produceDIADaPublication(manager, manifest, annotations, settings$da_q_val_thresh, staged),
            identity = diaDaPublicationIdentity(manifest))
        if (publish) completeDIADaDisplay(da_data, workflow, manifest)
        manifest
    })
}

#' The publication identity of one accepted DA manifest (the request key of its action outputs).
#' @noRd
diaDaPublicationIdentity <- function(manifest) digest::digest(list(manifest, output_schema = 2L), algo = "sha256")

#' The DA tab's status text once a manifest is accepted (the render and the managed step's expected text).
#' @noRd
diaDaStatusText <- function(contrasts) paste("[OK] Analysis Complete\n",
    sprintf("Contrasts analyzed: %d\n", contrasts), "Results available in all tabs")

#' Hand an accepted DA manifest to the module's holders and the tab status (after its publication).
#' @noRd
completeDIADaDisplay <- function(da_data, workflow, manifest) {
    da_data$da_results_list <- manifest
    da_data$analysis_complete <- TRUE
    workflow$da_analysis_results_list <- manifest
    status <- workflow$tab_status
    status$differential_expression <- status$differential_abundance <- "complete"
    status$enrichment_analysis <- "pending"
    workflow$tab_status <- status
    invisible(manifest)
}

#' Produce the DA publication into staged directories (KIT-019: the one body the in-process path and the
#' worker share): the annotated exports, the static volcano and NumSig figures through the unchanged
#' `outputDaResultsAllContrasts` kernel, and the interactive Glimma pages; every required file must exist.
#' @param manager Owner (the active manager of the session scope).
#' @param manifest Accepted `dia.da_results.v1` manifest.
#' @param annotations UniProt annotation table or NULL.
#' @param qvalue Display significance threshold.
#' @param staged `list(da_output_dir, publication_graphs_dir)` staging directories.
#' @param save_assets Save each page as a rendered asset here (the in-process owner); the worker leaves that
#'   to the controller's publication, which holds the head.
#' @return `list(completed = TRUE)`.
#' @noRd
produceDIADaPublication <- function(manager, manifest, annotations, qvalue, staged, save_assets = TRUE) {
    metadata_object <- readDIAExecutionComponents(manager, manifest$source, fields = "design_matrix")
    metadata_object@args$outputDaResultsAllContrasts <- list(uniprot_tbl = annotations,
        da_output_dir = staged$da_output_dir, publication_graphs_dir = staged$publication_graphs_dir,
        file_prefix = "da_proteins", args_row_id = manifest$protein_id_column,
        gene_names_column = "gene_names", uniprot_id_column = "Entry",
        da_q_val_thresh = qvalue, fdr_column = "fdr_qvalue", log2fc_column = "log2FC")
    if (!isTRUE(outputDaResultsAllContrasts(metadata_object, manifest$results,
        uniprot_tbl = annotations, da_output_dir = staged$da_output_dir,
        publication_graphs_dir = staged$publication_graphs_dir, file_prefix = "da_proteins",
        args_row_id = manifest$protein_id_column, gene_names_column = "gene_names", uniprot_id_column = "Entry"))) {
        cli::cli_abort("required DA outputs did not complete", class = "saved_entity_outputs_pending")
    }
    interactive_files <- writeDIAInteractiveVolcanoOutputs(manifest, staged$publication_graphs_dir, annotations, qvalue,
        manager = if (save_assets) manager else NULL)
    assertDIADaPublicationFiles(manifest, staged, interactive_files)
    list(completed = TRUE)
}

#' The files a DA publication must have staged, by contrast (the exports, the figures, the pages).
#' @noRd
assertDIADaPublicationFiles <- function(manifest, staged, interactive_files) {
    safe <- gsub("[^A-Za-z0-9_-]", "_", names(manifest$results))
    required <- c(file.path(staged$da_output_dir, as.vector(outer(
        paste0("da_proteins_", safe, "_long_annot"), c("tsv", "xlsx"), paste, sep = "."))),
        file.path(staged$publication_graphs_dir, "Volcano_Plots", as.vector(outer(safe, c("png", "pdf"), paste, sep = "."))),
        file.path(staged$publication_graphs_dir, "Volcano_Plots", "all_volcano_plots_combined.pdf"),
        file.path(staged$publication_graphs_dir, "NumSigDaMolecules", paste0("da_proteins_num_sig_da_molecules.", c("tab", "xlsx"))),
        interactive_files)
    missing <- required[!file.exists(required) | file.size(required) <= 0]
    if (length(missing)) cli::cli_abort(c("required DA files are incomplete", i = paste(basename(missing), collapse = ", ")),
        class = "saved_entity_outputs_pending")
    invisible(required)
}


#' The service identity an enrichment run reports (KIT-016): the g:Profiler database version carried by each
#' gost result's metadata, the deterministic backend's name under the test double, or the clusterProfiler
#' package version; `unavailable` when a run reports none.
#' @param results `EnrichmentResults`.
#' @return Sorted unique character vector.
#' @noRd
diaEnrichmentServiceIdentity <- function(results) {
    found <- character()
    for (contrast in names(results@enrichment_data)) {
        for (direction in names(results@enrichment_data[[contrast]])) {
            entry <- results@enrichment_data[[contrast]][[direction]]
            if (is.null(entry)) next
            if (is.list(entry) && !isS4(entry)) {
                # A live gost result reports the database version in its metadata; the deterministic double marks
                # every result row with its backend name.
                version <- entry$meta$version
                if (is.null(version) && is.data.frame(entry$result) && "deterministic_backend" %in% names(entry$result))
                    version <- unique(as.character(entry$result$deterministic_backend))
                version <- version[!is.na(version) & nzchar(version)]
                if (is.character(version) && length(version)) found <- c(found, paste0("gprofiler2:", version))
            } else if (isS4(entry) && methods::is(entry, "enrichResult")) {
                found <- c(found, paste0("clusterProfiler:", as.character(utils::packageVersion("clusterProfiler"))))
            }
        }
    }
    if (!length(found)) "unavailable" else sort(unique(found))
}

#' The enrichment analysis as a declaration (KIT-016): whole matrix over the unchanged `processEnrichments`
#' kernel (or the deterministic test double the module selects in test mode), the DA results, annotations
#' and pathway directory as extra inputs beside the declared parameters; the output class is the KIT-003
#' `EnrichmentResults`. The declaration is not cacheable (`factors$cacheable`): the g:Profiler database the
#' call reaches is known only after the call, so its identity is saved with the product afterwards.
#' @noRd
diaEnrichmentOperation <- function() {
    newWorkflowOperation("enrichment", "ProteinQuantitativeData", fields = "design_matrix", unit = "whole_matrix",
        output_class = "EnrichmentResults", product_type = "dia-enrichment-results",
        factors = list(cacheable = FALSE, service = "enrichment"),
        parameters = c("taxon_id", "up_cutoff", "down_cutoff", "q_cutoff", "exclude_iea", "correction_method",
            "protein_id_column", "contrast_names", "results"),
        global_statistics = c(per_contrast_da_table =
            "the background universe and the significant set are the complete DA table of each contrast"),
        kernel = function(object, parameters, store) {
            process <- parameters$process_fn %||% processEnrichments
            arguments <- list(da_results = parameters$da_results, taxon_id = parameters$taxon_id,
                up_cutoff = parameters$up_cutoff, down_cutoff = parameters$down_cutoff, q_cutoff = parameters$q_cutoff,
                pathway_dir = parameters$pathway_dir, go_annotations = parameters$go_annotations,
                exclude_iea = parameters$exclude_iea, protein_id_column = parameters$protein_id_column,
                contrast_names = parameters$contrast_names, correction_method = parameters$correction_method)
            results <- do.call(process, arguments)
            if (!methods::is(results, "EnrichmentResults"))
                cli::cli_abort("enrichment kernel returned no EnrichmentResults", class = "workflow_operation_output")
            list(results = results, service = diaEnrichmentServiceIdentity(results))
        })
}

#' Run the enrichment declaration on a DA manifest's source state and record its product (KIT-016).
#' @param workflow Session owner (`state_manager`).
#' @param manifest Accepted `dia.da_results.v1` manifest.
#' @param arguments The module's `processEnrichments` arguments (`da_results`, `taxon_id`, cutoffs,
#'   `pathway_dir`, `go_annotations`, `exclude_iea`, `protein_id_column`, `contrast_names`,
#'   `correction_method`).
#' @param process_fn The selected kernel: `processEnrichments`, or the deterministic double in test mode.
#' @return `list(results, service, reference)`: the `EnrichmentResults`, the reported service identity and
#'   the saved product reference (identity: source state, the contrast result tables, the declared
#'   parameters and the service identity; rooted on the source and the result tables).
#' @noRd
runDIAReferenceEnrichment <- function(workflow, manifest, arguments, process_fn = NULL) {
    manager <- workflow$state_manager
    if (!isDIAResultManifest(manifest)) cli::cli_abort("enrichment needs a DA result manifest", class = "saved_entity_invalid")
    declaration <- diaEnrichmentOperation()
    results <- vapply(names(manifest$results), function(name)
        manifest$results[[name]]$tables$da_proteins_long$entity_id %||% NA_character_, character(1))
    parameters <- c(arguments[setdiff(declaration$parameters, "results")], list(results = results,
        da_results = arguments$da_results, pathway_dir = arguments$pathway_dir,
        go_annotations = arguments$go_annotations, process_fn = process_fn))
    value <- runProjected(manager, declaration, manifest$source, parameters)
    identity <- diaComponentScientificIdentity(manager, manifest$source, "enrichment",
        c(parameters[declaration$parameters], list(service = value$service)))
    product <- list(schema = "dia_enrichment_product_v1", service = value$service,
        parameters = parameters[setdiff(declaration$parameters, "results")], results = results,
        contrasts = names(value$results@enrichment_data), summaries = value$results@enrichment_summaries)
    reference <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowScientificProduct(product, declaration$product_type,
            provenance = list(service = value$service, cacheable = FALSE)),
        basename(tempfile("dia-enrichment-")), identity, "scientific_product", factors = declaration$factors)
    tables <- lapply(manifest$results, function(result) result$tables$da_proteins_long)
    rootDIAProduct(manager, reference, c(list(manifest$source), tables))
    # KIT-026: the service each reported identity names, in the head beside the product.
    for (service in setdiff(value$service, "unavailable")) {
        parts <- strsplit(service, ":", fixed = TRUE)[[1L]]
        recordWorkflowExternalResource(manager, parts[[1L]], paste(parts[-1L], collapse = ":"),
            if (identical(parts[[1L]], "clusterProfiler")) "package_version" else "database_version",
            query_summary = sprintf("enrichment of %d contrasts, taxon %s", length(results), arguments$taxon_id %||% NA))
    }
    list(results = value$results, service = value$service, reference = reference)
}

#' The static volcano display as a saved product (KIT-016): the selected display columns of one contrast's
#' result table, saved once per result table and rebuilt into the plot on demand by the unchanged
#' `generateProtDAVolcanoStatic` kernel (thresholds and labels are render choices, not identity).
#' @param manifest,contrast Accepted results and selected contrast.
#' @return The product (`schema`, `data`, `columns`, `rows`, `comparison`, `protein_id_column`).
#' @noRd
generateDIAReferenceVolcanoProduct <- function(manifest, contrast) {
    manager <- getDIAActiveManager()
    selected <- selectDIAContrastReference(manifest, contrast)
    reference <- selected$tables$da_proteins_long
    identity <- diaComponentScientificIdentity(manager, manifest$source, "DA-volcano-display",
        list(result = reference$entity_id))
    hit <- resolveSavedWorkflowEntity(manager$.context$store, "scientific_product", identity)
    if (identical(hit$status, "hit")) return(hit$value$value)
    columns <- intersect(unique(c(manifest$protein_id_column, "comparison", "log2FC", "raw_pvalue", "fdr_qvalue",
        "gene_name", "gene_names", "Protein.Ids", "Entry", "description")), selected$columns)
    data <- readDIAContrastOutput(selected, "da_proteins_long", columns)$da_proteins_long
    product <- list(schema = "dia_volcano_display_v1", contrast = contrast, comparison = selected$comparison,
        protein_id_column = manifest$protein_id_column, columns = columns, rows = nrow(data), data = data)
    if (as.numeric(object.size(product)) > 67108864)
        cli::cli_abort("volcano display exceeds product limits", class = "saved_entity_limit")
    saved <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowScientificProduct(product, "dia-volcano-display"),
        basename(tempfile("dia-volcano-")), identity, "scientific_product")
    rootDIAProduct(manager, saved, list(manifest$source, reference))
    product
}

#' Identity of one contrast's interactive volcano asset (KIT-016): the result table and the display
#' threshold it was rendered at; the annotation table joins at publication and is not part of the identity.
#' @noRd
diaInteractiveVolcanoIdentity <- function(manager, manifest, contrast, qvalue) {
    selected <- selectDIAContrastReference(manifest, contrast)
    diaComponentScientificIdentity(manager, manifest$source, "DA-volcano-interactive",
        list(result = selected$tables$da_proteins_long$entity_id, qvalue = as.numeric(qvalue)))
}

#' Save one contrast's self-contained interactive volcano page as a rendered asset (KIT-016), rooted on the
#' source state and the contrast's result table.
#' @noRd
saveDIAInteractiveVolcanoAsset <- function(manager, manifest, contrast, qvalue, bytes) {
    identity <- diaInteractiveVolcanoIdentity(manager, manifest, contrast, qvalue)
    # The page's display box in the module (the asset envelope needs positive dimensions).
    asset <- newSavedWorkflowRenderedAsset(bytes, "text/html", c(1200L, 820L))
    reference <- saveWorkflowConsumerEntity(manager$.context$store, asset,
        basename(tempfile("dia-volcano-html-")), identity, "rendered_asset")
    rootDIAProduct(manager, reference, list(manifest$source, selectDIAContrastReference(manifest, contrast)$tables$da_proteins_long))
    reference
}

#' Resolve one contrast's saved interactive volcano page (KIT-016).
#' @return The HTML as a string, or NULL when no asset exists for the result and threshold.
#' @noRd
resolveDIAInteractiveVolcanoAsset <- function(manifest, contrast, qvalue) {
    manager <- getDIAActiveManager()
    identity <- diaInteractiveVolcanoIdentity(manager, manifest, contrast, qvalue)
    hit <- resolveSavedWorkflowEntity(manager$.context$store, "rendered_asset", identity)
    if (!identical(hit$status, "hit") || !identical(hit$value$status, "available")) return(NULL)
    rawToChar(hit$value$bytes)
}

#' Export linked-expression volcano widgets from accepted result references.
#' @param manifest Accepted DIA result manifest.
#' @param directory Owned publication staging directory.
#' @param annotations Existing bounded annotation table.
#' @param qvalue Display significance threshold.
#' @return Paths to HTML files and their portable dependency files.
#' @noRd
writeDIAInteractiveVolcanoOutputs <- function(manifest, directory, annotations, qvalue, manager = NULL) {
    target <- file.path(directory, "Interactive_Volcano_Plots")
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("dia-volcano-assets-")
    dir.create(temporary)
    on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
    for (contrast in manifest$contrasts$friendly_names) {
        widget <- generateProtDAVolcanoPlotGlimma(manifest, contrast,
            uniprot_tbl = annotations, args_row_id = manifest$protein_id_column,
            da_q_val_thresh = qvalue)
        if (is.null(widget)) cli::cli_abort("interactive volcano output is missing",
            class = "saved_entity_outputs_pending")
        filename <- paste0(gsub("[^A-Za-z0-9_-]", "_", contrast), ".html")
        # One self-contained file, written in R (KIT-036): saveWidget's pandoc pass over this page is the DA
        # step's memory peak (KIT-017 audit PH-01). Only the HTML belongs to the durable publication inventory.
        html <- file.path(temporary, filename)
        saveSelfContainedWidget(widget, html)
        if (!file.copy(html, file.path(target, filename), overwrite = TRUE))
            cli::cli_abort("cannot stage interactive volcano HTML",
                class = "saved_entity_outputs_pending")
        # KIT-016: the self-contained page is also a saved rendered asset, so the module renders it on demand
        # without rebuilding the widget (and without loading Glimma) in the application process.
        if (!is.null(manager)) saveDIAInteractiveVolcanoAsset(manager, manifest, contrast, qvalue,
            readBin(html, "raw", n = file.size(html)))
        widget <- NULL
    }
    list.files(target, recursive = TRUE, full.names = TRUE)
}

#' Build the immutable identity for an uncontrasted DIA limma fit.
#' @param data Ordered feature-by-sample quantities.
#' @param design Ordered sample metadata.
#' @param formula_string Existing model formula string.
#' @param block Effective duplicate-correlation block or NULL.
#' @param weights Existing weights argument, including scalar NA.
#' @param normalization,ruv Effective upstream scientific identities.
#' @param software Installed software/build identity.
#' @return Shared-store scientific identity.
#' @noRd
newWorkflowDiaBaseFitIdentity <- function(data, design, formula_string, block,
    weights, normalization, ruv, software) {
    if (!is.matrix(data) || !is.numeric(data) || is.null(rownames(data)) ||
        is.null(colnames(data)) || !is.data.frame(design) ||
        !identical(colnames(data), rownames(design)) ||
        !(is.null(block) || is.character(block) && length(block) == ncol(data))) {
        cli::cli_abort("invalid DIA base-fit identity inputs", class = "saved_entity_invalid")
    }
    for (value in list(formula_string, normalization, ruv, software)) {
        checkSavedEntityString(value, "DIA fit identity")
    }
    input <- digest::digest(list(values = data, feature_order = rownames(data),
        sample_order = colnames(data)), algo = "sha256")
    design_id <- digest::digest(list(formula = formula_string, design = design),
        algo = "sha256")
    parameters <- digest::digest(list(block = block, weights = weights,
        normalization = normalization, ruv = ruv), algo = "sha256")
    selectors <- digest::digest(list(features = rownames(data),
        samples = colnames(data)), algo = "sha256")
    newSavedWorkflowScientificIdentity(input, design_id, parameters, selectors, software)
}

#' Build a downstream result identity from one immutable base-fit descriptor.
#' @param base_fit_reference Accepted da_fit descriptor.
#' @param contrast_strings Ordered contrast definitions.
#' @param treat_lfc_cutoff,eBayes_trend,eBayes_robust Existing statistics settings.
#' @param adjustment Existing qvalue/BH contract identity.
#' @param software Installed software/build identity.
#' @return Shared-store scientific identity.
#' @noRd
newWorkflowDiaResultIdentity <- function(base_fit_reference, contrast_strings,
    treat_lfc_cutoff, eBayes_trend, eBayes_robust, adjustment, software) {
    if (!is.list(base_fit_reference) || !identical(base_fit_reference$kind, "da_fit") ||
        !grepl("^[a-f0-9]{64}$", base_fit_reference$entity_id) ||
        !is.character(contrast_strings) || !length(contrast_strings) ||
        anyNA(contrast_strings) || anyDuplicated(contrast_strings) ||
        !is.logical(eBayes_trend) || length(eBayes_trend) != 1L ||
        !is.logical(eBayes_robust) || length(eBayes_robust) != 1L) {
        cli::cli_abort("invalid DIA result identity inputs", class = "saved_entity_invalid")
    }
    checkSavedEntityString(adjustment, "adjustment")
    checkSavedEntityString(software, "software")
    newSavedWorkflowScientificIdentity(base_fit_reference$entity_id,
        "dia-contrast-result-v1",
        digest::digest(list(treat_lfc_cutoff = treat_lfc_cutoff,
            eBayes_trend = eBayes_trend, eBayes_robust = eBayes_robust,
            adjustment = adjustment), algo = "sha256"),
        digest::digest(contrast_strings, algo = "sha256"), software)
}

#' Bind Limpa's actual quantified EList, including uncertainty support.
#' @param quantified_elist Actual EList passed to limpa::dpcDE.
#' @param design Ordered design metadata.
#' @param formula_string Existing model formula.
#' @param software Limpa/MultiScholaR build identity.
#' @return Identity that cannot be authorized by displayed quantities alone.
#' @noRd
newWorkflowLimpaFitIdentity <- function(quantified_elist, design,
    formula_string, software) {
    if (!methods::is(quantified_elist, "EList") ||
        !is.matrix(quantified_elist$E) || !is.data.frame(quantified_elist$genes) ||
        !is.list(quantified_elist$other) || !is.data.frame(design) ||
        is.null(rownames(quantified_elist$E)) || is.null(colnames(quantified_elist$E))) {
        cli::cli_abort("invalid Limpa fit identity inputs", class = "saved_entity_invalid")
    }
    checkSavedEntityString(formula_string, "formula_string")
    checkSavedEntityString(software, "software")
    newSavedWorkflowScientificIdentity(
        digest::digest(list(E = quantified_elist$E, genes = quantified_elist$genes,
            other = quantified_elist$other), algo = "sha256"),
        digest::digest(list(formula = formula_string, design = design),
            algo = "sha256"),
        digest::digest(list(uncertainty_members = names(quantified_elist$other)),
            algo = "sha256"),
        digest::digest(list(features = rownames(quantified_elist$E),
            samples = colnames(quantified_elist$E)), algo = "sha256"),
        software)
}

#' Save a real uncontrasted fit through the shared da_fit contract.
#' @param store Shared store.
#' @param base_fit_artifact Artifact returned by fitProtDaBaseModel.
#' @param request_id Idempotent request ID.
#' @param identity Complete base-fit identity.
#' @return Accepted da_fit descriptor.
#' @noRd
saveWorkflowDiaBaseFit <- function(store, base_fit_artifact, request_id, identity) {
    if (!is.list(base_fit_artifact) ||
        !identical(base_fit_artifact$schema, "prot_da_base_fit_v1")) {
        cli::cli_abort("invalid DIA base-fit artifact", class = "saved_entity_invalid")
    }
    envelope <- newSavedWorkflowDaFit(base_fit_artifact$fit,
        base_fit_artifact$feature_ids, base_fit_artifact$sample_ids,
        base_fit_artifact$design_columns)
    saveWorkflowConsumerEntity(store, envelope, request_id, identity, "da_fit",
        factors = list(route = diaActiveRoute(), layer = "base_fit",
            uses_duplicate_correlation = base_fit_artifact$uses_duplicate_correlation,
            consensus_correlation = base_fit_artifact$consensus_correlation,
            # KIT-035: which model ran and why (provenance for the analysis record).
            group_id = base_fit_artifact$group_id, technical_replicate_id = base_fit_artifact$technical_replicate_id,
            block_reason = base_fit_artifact$block_reason))
}

#' Resolve one exact DIA base fit without fitting.
#' @param store Shared store.
#' @param identity Complete base-fit identity.
#' @param uses_duplicate_correlation,consensus_correlation Effective provenance.
#' @return Explicit fit hit or miss.
#' @noRd
resolveWorkflowDiaBaseFit <- function(store, identity,
    uses_duplicate_correlation, consensus_correlation, group_id = "group", technical_replicate_id = "replicates",
    block_reason = if (isTRUE(uses_duplicate_correlation)) "technical_replicates" else "no_duplicates") {
    # KIT-035: the saved factors carry the model's provenance, so an exact resolve names it too.
    resolveSavedWorkflowEntity(store, "da_fit", identity,
        factors = list(route = diaActiveRoute(), layer = "base_fit",
            uses_duplicate_correlation = uses_duplicate_correlation,
            consensus_correlation = consensus_correlation,
            group_id = group_id, technical_replicate_id = technical_replicate_id, block_reason = block_reason))
}

#' Save the complete contrasted result without recomputing statistics on hits.
#' @param store Shared store.
#' @param result Complete result returned by buildProtDaResultTables.
#' @param request_id Idempotent request ID.
#' @param identity Complete result identity.
#' @return Accepted scientific-product descriptor.
#' @noRd
saveWorkflowDiaResult <- function(store, result, request_id, identity) {
    if (!is.list(result) || !identical(result$schema, "prot_da_result_tables_v1")) {
        cli::cli_abort("invalid DIA result artifact", class = "saved_entity_invalid")
    }
    envelope <- newSavedWorkflowScientificProduct(result, "dia_contrast_results",
        provenance = list(route = diaActiveRoute(), full_population = TRUE,
            contrasts = names(result$results)))
    saveWorkflowConsumerEntity(store, envelope, request_id, identity,
        "scientific_product", factors = list(route = diaActiveRoute(), layer = "full_results"))
}

#' Resolve one exact complete DA result without statistical computation.
#' @param store Shared store.
#' @param identity Complete result identity.
#' @return Explicit result hit or miss.
#' @noRd
resolveWorkflowDiaResult <- function(store, identity) {
    resolveSavedWorkflowEntity(store, "scientific_product", identity,
        factors = list(route = diaActiveRoute(), layer = "full_results"))
}

#' Return one deterministic page from already adjusted full-population statistics.
#' @param result_table Result-table envelope.
#' @param page,page_size Positive integers.
#' @param order_by Existing result column.
#' @param decreasing Sort direction.
#' @return Page plus full population metadata; never reruns tests/FDR.
#' @noRd
pageWorkflowDiaResult <- function(result_table, page, page_size,
    order_by = "fdr_qvalue", decreasing = FALSE) {
    checkSavedWorkflowConsumerPayload("result_table", result_table)
    if (!is.numeric(page) || length(page) != 1L || page < 1 || page != floor(page) ||
        !is.numeric(page_size) || length(page_size) != 1L || page_size < 1 ||
        page_size != floor(page_size) || !order_by %in% names(result_table$table) ||
        !is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
        cli::cli_abort("invalid DIA result page request", class = "saved_entity_invalid")
    }
    values <- result_table$table[[order_by]]
    order <- order(values, seq_along(values), decreasing = decreasing, na.last = TRUE)
    start <- (page - 1L) * page_size + 1L
    indices <- if (start > length(order)) integer() else
        order[seq.int(start, min(length(order), start + page_size - 1L))]
    list(data = result_table$table[indices, , drop = FALSE], page = as.integer(page),
        page_size = as.integer(page_size), full_population_rows = nrow(result_table$table),
        statistics_recomputed = FALSE)
}

#' Filter a volcano view from already adjusted full-population results.
#' @param result_table Result-table envelope.
#' @param qvalue_threshold,absolute_lfc_threshold Display thresholds.
#' @return View rows and immutable full-population count.
#' @noRd
filterWorkflowDiaVolcano <- function(result_table, qvalue_threshold,
    absolute_lfc_threshold) {
    checkSavedWorkflowConsumerPayload("result_table", result_table)
    table <- result_table$table
    if (!all(c("fdr_qvalue", "logFC") %in% names(table)) ||
        !is.numeric(qvalue_threshold) || length(qvalue_threshold) != 1L ||
        !is.numeric(absolute_lfc_threshold) || length(absolute_lfc_threshold) != 1L) {
        cli::cli_abort("invalid DIA volcano view request", class = "saved_entity_invalid")
    }
    selected <- !is.na(table$fdr_qvalue) & is.finite(table$fdr_qvalue) &
        table$fdr_qvalue <= qvalue_threshold & !is.na(table$logFC) &
        is.finite(table$logFC) & abs(table$logFC) >= absolute_lfc_threshold
    list(data = table[selected, , drop = FALSE], selected = selected,
        full_population_rows = nrow(table), statistics_recomputed = FALSE)
}

#' Select heatmap rows from full DA statistics and the complete tested matrix.
#' @param result_table Result-table envelope.
#' @param matrix Complete tested feature-by-sample matrix.
#' @param top_n Positive number of rows.
#' @param scaling,clustering Existing display settings.
#' @param cluster_method,distance_method,tree_cut_method,n_clusters,color_scheme,show_labels
#'   Existing clustering and display inputs retained for the renderer.
#' @return Ordered view matrix and exact selection inputs.
#' @noRd
selectWorkflowDiaHeatmap <- function(result_table, matrix, top_n,
    scaling = "row", clustering = "both", cluster_method = "ward.D2",
    distance_method = "euclidean", tree_cut_method = "k_clusters",
    n_clusters = 4L, color_scheme = "RdBu", show_labels = TRUE) {
    checkSavedWorkflowConsumerPayload("result_table", result_table)
    table <- result_table$table
    ids <- if (identical(result_table$feature_id, ".rownames")) rownames(table) else
        as.character(table[[result_table$feature_id]])
    if (!is.matrix(matrix) || is.null(rownames(matrix)) || is.null(colnames(matrix)) ||
        !all(ids %in% rownames(matrix)) || !is.numeric(top_n) || length(top_n) != 1L ||
        top_n < 1 || top_n != floor(top_n) || !scaling %in% c("none", "row") ||
        !clustering %in% c("none", "row", "column", "both") ||
        !is.character(cluster_method) || length(cluster_method) != 1L ||
        !is.character(distance_method) || length(distance_method) != 1L ||
        !is.character(tree_cut_method) || length(tree_cut_method) != 1L ||
        !is.integer(n_clusters) || length(n_clusters) != 1L || n_clusters < 1L ||
        !is.character(color_scheme) || length(color_scheme) != 1L ||
        !is.logical(show_labels) || length(show_labels) != 1L || is.na(show_labels)) {
        cli::cli_abort("invalid DIA heatmap view request", class = "saved_entity_invalid")
    }
    ordering <- order(table$fdr_qvalue, -abs(table$logFC), seq_len(nrow(table)),
        na.last = TRUE)
    selected <- ordering[seq_len(min(length(ordering), top_n))]
    selected_ids <- ids[selected]
    view <- matrix[selected_ids, , drop = FALSE]
    if (identical(scaling, "row")) view <- t(scale(t(view)))
    list(matrix = view, feature_ids = selected_ids, sample_ids = colnames(matrix),
        source_rows = selected, full_population_rows = nrow(table),
        full_tested_features = nrow(matrix), scaling = scaling,
        clustering = clustering, cluster_method = cluster_method,
        distance_method = distance_method, tree_cut_method = tree_cut_method,
        n_clusters = n_clusters, color_scheme = color_scheme,
        show_labels = show_labels, statistics_recomputed = FALSE)
}

#' Resolve an interactive feature selection against the complete tested matrix.
#' @param matrix Complete tested feature-by-sample matrix.
#' @param feature_id Selected feature identity.
#' @return Ordered sample/value projection with no statistical recomputation.
#' @noRd
selectWorkflowDiaExpression <- function(matrix, feature_id) {
    if (!is.matrix(matrix) || is.null(rownames(matrix)) || is.null(colnames(matrix)) ||
        !is.character(feature_id) || length(feature_id) != 1L || is.na(feature_id) ||
        !feature_id %in% rownames(matrix)) {
        cli::cli_abort("invalid DIA expression selection", class = "saved_entity_invalid")
    }
    data.frame(sample_id = colnames(matrix), value = unname(matrix[feature_id, ]),
        feature_id = rep(feature_id, ncol(matrix)), stringsAsFactors = FALSE)
}

#' Reuse a narrow DA display-selection index without repeating fit or FDR work.
#' @param manager Current saved-result owner.
#' @param reference One immutable complete contrast-result table.
#' @param qvalue,lfc Existing display cutoffs, included in selection identity.
#' @return Selected ordinals/flags and full-population counts, never result collections.
#' @noRd
resolveDIAResultSelection <- function(manager, reference, qvalue, lfc) {
    if (length(qvalue) != 1L || !is.finite(qvalue) || qvalue < 0 || qvalue > 1 ||
        length(lfc) != 1L || !is.finite(lfc) || lfc < 0) cli::cli_abort("invalid result cutoffs")
    identity <- newSavedWorkflowScientificIdentity(reference$entity_id, "full-testing-population",
        digest::digest(list(qvalue, lfc), algo = "sha256"), "DA-display-selection-v1", manager$.context$software)
    store <- manager$.context$store
    hit <- resolveSavedWorkflowEntity(store, "scientific_product", identity)
    if (identical(hit$status, "hit")) return(hit$value$value)
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    metadata <- workflowSavedTableMetadata(data.frame(source_order = integer(), up = logical(), down = logical()),
        "source_order", 3L)
    metadata$rows <- NULL
    counts <- data.frame(total = as.double(scope$metadata$rows), significant = 0, up = 0, down = 0)
    rows <- saveWorkflowTable(store, NULL, basename(tempfile("da-selection-")), identity, "source_order",
        metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, store$progress)
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            streamSavedWorkflowTable(scope, c("fdr_qvalue", "log2FC"), function(data, offset, row_names) {
                significant <- !is.na(data$fdr_qvalue) & data$fdr_qvalue < qvalue
                # Preserve the existing DuckDB display comparison for NaN effects.
                up <- significant & ((!is.na(data$log2FC) & data$log2FC > lfc) | is.nan(data$log2FC))
                down <- significant & !is.na(data$log2FC) & data$log2FC < -lfc
                counts$significant <<- counts$significant + sum(significant)
                counts$up <<- counts$up + sum(up)
                counts$down <<- counts$down + sum(down)
                selected <- which(significant)
                appendWorkflowTableWriter(writer, data.frame(source_order = as.integer(offset + selected),
                    up = up[selected], down = down[selected]))
            })
            finishWorkflowTableWriter(writer)
        })
    value <- list(rows = rows, counts = counts)
    product <- saveWorkflowConsumerEntity(store, newSavedWorkflowScientificProduct(value, "DA-display-selection"),
        basename(tempfile("da-selection-product-")), identity, "scientific_product")
    rootDIAProduct(manager, product, list(reference, rows))
    value
}
