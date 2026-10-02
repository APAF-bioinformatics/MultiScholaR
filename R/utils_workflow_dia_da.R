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
    valid_state <- manager$current_state %in% c("normalized", "ruv_corrected", "correlation_filtered")
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
fitDIAExecutionBaseModel <- function(data, design_m, design, samples, formula, weights) {
    manager <- getDIAActiveManager()
    if (!usesDIAReferenceStates(manager) || !isTRUE(manager$.context$da_reference_execution)) {
        return(fitProtDaBaseModel(data, design_m, design, samples))
    }
    selected_design <- design[samples, , drop = FALSE]
    block <- if ("replicates" %in% names(design))
        paste(selected_design$group, selected_design$replicates, sep = "_") else NULL
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
            uses_duplicate_correlation = !is.null(block) && any(duplicated(block)),
            consensus_correlation = hit$reference$metadata$factors$consensus_correlation))
    }
    fit <- fitProtDaBaseModel(data, design_m, design, samples)
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
    lfc = 0, offset = 0L, limit = 100L) {
    manager <- getDIAActiveManager()
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

#' Run native DA one contrast at a time and retain only accepted references.
#' @param input,output,session,da_data,workflow,paths Existing module owners.
#' @return Reference manifest after required outputs have completed.
#' @noRd
runDIAReferenceDaAnalysis <- function(input, output, session, da_data, workflow, paths) {
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "da")
        source <- da_data$current_s4_object
        if (!isSavedWorkflowReference(source)) cli::cli_abort("load a reference-backed protein state before DA")
        contrasts <- workflow$contrasts_tbl
        if (is.null(contrasts) && exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE))
            contrasts <- get("contrasts_tbl", envir = .GlobalEnv)
        if (is.null(contrasts)) contrasts <- data.frame(contrasts = da_data$contrasts_available)
        contrasts <- normaliseProtDaContrastsTable(contrasts)
        spec <- readDIAStateMetadata(manager, source)$data
        prototype <- readDIAExecutionComponents(manager, source, fields = "design_matrix")
        validateProtDaModelAndContrasts(prototype, input$formula_string, contrasts)
        manifest <- list(schema = "dia.da_results.v1", source = source,
            protein_id_column = spec$slots$protein_id_column, sample_id = spec$slots$sample_id,
            contrasts = contrasts, results = list(), parameters = list(formula = input$formula_string,
                qvalue = input$da_q_val_thresh, lfc = input$treat_lfc_cutoff, weights = NA))
        context <- manager$.context
        context$da_reference_execution <- TRUE
        context$da_source_reference <- source
        on.exit({ context$da_reference_execution <- FALSE; context$da_source_reference <- NULL }, add = TRUE)
        object <- NULL
        for (i in seq_len(nrow(contrasts))) {
            name <- contrasts$contrasts[[i]]
            identity <- diaComponentScientificIdentity(manager, source, "DA-contrast-output",
                c(manifest$parameters, list(contrast = contrasts$full_format[[i]],
                    annotations = spec$components$protein_id_table$entity_id)))
            hit <- resolveSavedWorkflowEntity(context$store, "scientific_product", identity,
                factors = list(route = "iq", layer = "contrast_references"))
            if (identical(hit$status, "hit")) {
                manifest$results[[name]] <- hit$value$value
                next
            }
            if (is.null(object)) object <- readDIAExecutionComponents(manager, source)
            single <- data.frame(contrasts = contrasts$full_format[[i]], stringsAsFactors = FALSE)
            parameters <- list(contrasts_tbl = single, formula_string = input$formula_string,
                da_q_val_thresh = input$da_q_val_thresh, treat_lfc_cutoff = input$treat_lfc_cutoff,
                eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
                args_row_id = object@protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
            for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) {
                if (is.null(object@args[[section]])) object@args[[section]] <- parameters
            }
            result <- differentialAbundanceAnalysis(object, contrasts_tbl = single,
                formula_string = input$formula_string, da_q_val_thresh = input$da_q_val_thresh,
                treat_lfc_cutoff = input$treat_lfc_cutoff, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
            accepted <- list(schema = "dia.contrast_reference.v1", source = source,
                contrast = name, comparison = contrasts$friendly_names[[i]],
                tables = list(), products = list(), qvalue_warnings = result$qvalue_warnings)
            for (field in c("da_proteins_long", "da_proteins_wide", "contrasts_results_table",
                "significant_rows", "norm_counts")) {
                if (is.data.frame(result[[field]])) accepted$tables[[field]] <- saveDIATableComponent(context, result[[field]])
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
                product_identity <- identity
                product_identity$selectors <- digest::digest(list(identity$selectors, field), algo = "sha256")
                accepted$products[[field]] <- saveWorkflowConsumerEntity(context$store,
                    newSavedWorkflowScientificProduct(value, paste0("dia-", field)),
                    basename(tempfile("dia-da-display-")), product_identity, "scientific_product")
            }
            accepted$columns <- names(result$da_proteins_long)
            accepted$rows <- nrow(result$da_proteins_long)
            result <- NULL
            reference <- saveWorkflowConsumerEntity(context$store,
                newSavedWorkflowScientificProduct(accepted, "dia-contrast-references"),
                basename(tempfile("dia-contrast-")), identity, "scientific_product",
                factors = list(route = "iq", layer = "contrast_references"))
            snapshot <- context$snapshot
            snapshot$da_dependencies[[reference$entity_id]] <- collectDIANativeReferences(accepted)
            publishDIAStateHead(context, snapshot)
            rootDIAProduct(manager, reference, list(source, accepted$tables, accepted$products))
            manifest$results[[name]] <- accepted
        }
        object <- prototype <- NULL
        workflow$da_ui_params <- list(q_value_threshold = input$da_q_val_thresh,
            log_fold_change_cutoff = input$treat_lfc_cutoff, treat_enabled = input$treat_lfc_cutoff > 0,
            formula_string = input$formula_string, timestamp = Sys.time())
        state_record <- readDIAStateMetadata(manager)
        checkpoint <- state_record$data
        da_parameters <- workflow$da_ui_params
        recorded <- checkpoint$slots$args$daAnalysisUI
        if (!identical(recorded[setdiff(names(recorded), "timestamp")],
            da_parameters[setdiff(names(da_parameters), c("timestamp", "treat_enabled"))])) {
            first <- data.frame(contrasts = contrasts$full_format[[1L]], stringsAsFactors = FALSE)
            defaults <- list(contrasts_tbl = first, formula_string = input$formula_string,
                da_q_val_thresh = input$da_q_val_thresh, treat_lfc_cutoff = input$treat_lfc_cutoff,
                eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
                args_row_id = manifest$protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
            if (is.null(checkpoint$slots$args$differentialAbundanceAnalysis))
                checkpoint$slots$args$differentialAbundanceAnalysis <- defaults
            checkpoint$slots$args$daAnalysisUI <- da_parameters[setdiff(names(da_parameters), "treat_enabled")]
            if (nrow(contrasts) > 1L && is.null(checkpoint$slots$args$differentialAbundanceAnalysisHelper))
                checkpoint$slots$args$differentialAbundanceAnalysisHelper <- defaults
            commitDIAReferenceState(manager, manager$current_state, checkpoint,
                state_record$config, state_record$description)
        }
        snapshot <- context$snapshot
        manifest_reference <- saveDIAStateValue(context, manifest, "da-result-references", codec = "rds")
        snapshot$da_results <- manifest_reference
        snapshot$da_dependencies[[manifest_reference$entity_id]] <- collectDIANativeReferences(
            list(manifest$source, manifest$results))
        publishDIAStateHead(context, snapshot)
        publishDIAActionOutputs(manager, getDIAStateReference(manager), "da",
            paths[c("da_output_dir", "publication_graphs_dir")], function(staged) {
                metadata_object <- readDIAExecutionComponents(manager, source, fields = "design_matrix")
                annotations <- workflow$uniprot_dat_cln
                if (is.null(annotations) && exists("uniprot_dat_cln", envir = .GlobalEnv, inherits = FALSE))
                    annotations <- get("uniprot_dat_cln", envir = .GlobalEnv)
                if (is.null(annotations)) annotations <- data.frame()
                metadata_object@args$outputDaResultsAllContrasts <- list(uniprot_tbl = annotations,
                    da_output_dir = staged$da_output_dir, publication_graphs_dir = staged$publication_graphs_dir,
                    file_prefix = "da_proteins", args_row_id = manifest$protein_id_column,
                    gene_names_column = "gene_names", uniprot_id_column = "Entry",
                    da_q_val_thresh = input$da_q_val_thresh, fdr_column = "fdr_qvalue", log2fc_column = "log2FC")
                if (!isTRUE(outputDaResultsAllContrasts(metadata_object, manifest$results,
                    uniprot_tbl = annotations, da_output_dir = staged$da_output_dir,
                    publication_graphs_dir = staged$publication_graphs_dir, file_prefix = "da_proteins",
                    args_row_id = manifest$protein_id_column, gene_names_column = "gene_names", uniprot_id_column = "Entry"))) {
                    cli::cli_abort("required DA outputs did not complete", class = "saved_entity_outputs_pending")
                }
                interactive_files <- writeDIAInteractiveVolcanoOutputs(manifest,
                    staged$publication_graphs_dir, annotations, input$da_q_val_thresh)
                safe <- gsub("[^A-Za-z0-9_-]", "_", names(manifest$results))
                required <- c(file.path(staged$da_output_dir, as.vector(outer(
                    paste0("da_proteins_", safe, "_long_annot"), c("tsv", "xlsx"), paste, sep = "."))),
                    file.path(staged$publication_graphs_dir, "Volcano_Plots", as.vector(outer(safe,
                        c("png", "pdf"), paste, sep = "."))),
                    file.path(staged$publication_graphs_dir, "Volcano_Plots", "all_volcano_plots_combined.pdf"),
                    file.path(staged$publication_graphs_dir, "NumSigDaMolecules",
                        paste0("da_proteins_num_sig_da_molecules.", c("tab", "xlsx"))))
                required <- c(required, interactive_files)
                if (any(!file.exists(required)) || any(file.size(required) <= 0)) cli::cli_abort(
                    "required DA files are incomplete", class = "saved_entity_outputs_pending")
                list(completed = TRUE)
            }, identity = digest::digest(list(manifest, output_schema = 2L), algo = "sha256"))
        da_data$da_results_list <- manifest
        da_data$analysis_complete <- TRUE
        workflow$da_analysis_results_list <- manifest
        choices <- stats::setNames(contrasts$friendly_names, contrasts$friendly_names)
        for (id in c("volcano_contrast", "heatmap_contrast", "table_contrast"))
            shiny::updateSelectInput(session, id, choices = choices)
        status <- workflow$tab_status
        status$differential_expression <- status$differential_abundance <- "complete"
        status$enrichment_analysis <- "pending"
        workflow$tab_status <- status
        warnings <- unlist(lapply(manifest$results, \(result) result$qvalue_warnings), recursive = FALSE)
        if (length(warnings)) shiny::showNotification(
            "Some contrasts used the existing BH fallback after q-value estimation failed.", type = "warning", duration = NULL)
        shiny::showNotification("Differential abundance analysis completed successfully!", type = "message", duration = 5)
        manifest
    })
}

#' Export linked-expression volcano widgets from accepted result references.
#' @param manifest Accepted DIA result manifest.
#' @param directory Owned publication staging directory.
#' @param annotations Existing bounded annotation table.
#' @param qvalue Display significance threshold.
#' @return Paths to HTML files and their portable dependency files.
#' @noRd
writeDIAInteractiveVolcanoOutputs <- function(manifest, directory, annotations, qvalue) {
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
        # saveWidget retains dependency files even with selfcontained = TRUE.
        # Only the embedded HTML belongs to the durable publication inventory.
        html <- file.path(temporary, filename)
        htmlwidgets::saveWidget(widget, html, selfcontained = TRUE)
        if (!file.copy(html, file.path(target, filename), overwrite = TRUE))
            cli::cli_abort("cannot stage interactive volcano HTML",
                class = "saved_entity_outputs_pending")
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
        factors = list(route = "iq", layer = "base_fit",
            uses_duplicate_correlation = base_fit_artifact$uses_duplicate_correlation,
            consensus_correlation = base_fit_artifact$consensus_correlation))
}

#' Resolve one exact DIA base fit without fitting.
#' @param store Shared store.
#' @param identity Complete base-fit identity.
#' @param uses_duplicate_correlation,consensus_correlation Effective provenance.
#' @return Explicit fit hit or miss.
#' @noRd
resolveWorkflowDiaBaseFit <- function(store, identity,
    uses_duplicate_correlation, consensus_correlation) {
    resolveSavedWorkflowEntity(store, "da_fit", identity,
        factors = list(route = "iq", layer = "base_fit",
            uses_duplicate_correlation = uses_duplicate_correlation,
            consensus_correlation = consensus_correlation))
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
        provenance = list(route = "iq", full_population = TRUE,
            contrasts = names(result$results)))
    saveWorkflowConsumerEntity(store, envelope, request_id, identity,
        "scientific_product", factors = list(route = "iq", layer = "full_results"))
}

#' Resolve one exact complete DA result without statistical computation.
#' @param store Shared store.
#' @param identity Complete result identity.
#' @return Explicit result hit or miss.
#' @noRd
resolveWorkflowDiaResult <- function(store, identity) {
    resolveSavedWorkflowEntity(store, "scientific_product", identity,
        factors = list(route = "iq", layer = "full_results"))
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
