# KIT-095: the transcriptomics binding of the shared omics declaration core (utils_workflow_omics_declarations.R), the
# step map `kit-001/step-map-trx.json`. KIT-093's kernels (R/func_trx_*.R) are called unchanged; only the DE kernel is
# taken as its two halves, the fit and the test, so that one saved fit serves one declaration per contrast (the step map's
# kernel_observations). States hold the count table and the design; logCPM, the grid, the RUVs display matrix, the DE
# results and the RUVIII matrix are products. Nothing here reads GUI inputs, the global environment or a Shiny session.

transcriptomicsOperationFields <- c("count_data", "design_matrix")

#' The merge rules of a partitioned transcriptomics declaration: the count table concatenates, the design is kept.
#' @noRd
transcriptomicsPartitionMerge <- c(count_data = "concatenate", design_matrix = "design")

#' The grid's declared parameters after its dependency `logcpm`.
#' @noRd
trxRuvGridParameterNames <- c("ruv_mode", "percentage_min", "percentage_max", "percentage_step", "k_penalty_weight", "manual_k",
    "manual_percentage")

#' The grid's parameters with the notebook's values as defaults (rmd:248-298; KIT-093's trxRuvGrid defaults).
#' @param ruv_mode `automatic` (the notebook) or `manual`; there is no skip mode.
#' @return The named list `trx_ruv_grid` takes after `logcpm`.
#' @noRd
trxRuvGridParameters <- function(ruv_mode = "automatic", percentage_min = 5, percentage_max = 30, percentage_step = 5,
                                 k_penalty_weight = 0.05, manual_k = NA_integer_, manual_percentage = NA_real_) {
    list(ruv_mode = ruv_mode, percentage_min = percentage_min, percentage_max = percentage_max, percentage_step = percentage_step,
        k_penalty_weight = k_penalty_weight, manual_k = manual_k, manual_percentage = manual_percentage)
}

# ---------------------------------------------------------------------------------------------------------------
# Tables, products and the template.
# ---------------------------------------------------------------------------------------------------------------

#' A genes-by-runs matrix as a saved table's frame: `GeneID`, then one column per run.
#' @noRd
trxMatrixTable <- function(m) data.frame(GeneID = rownames(m) %||% character(), m, check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)

#' A saved matrix table read back as the genes-by-runs matrix it was written from.
#' @param reference The saved table's reference.
#' @noRd
trxSavedMatrix <- function(manager, reference) {
    table <- materializeWorkflowTable(manager$.context, reference)
    m <- as.matrix(table[, setdiff(names(table), "GeneID"), drop = FALSE])
    rownames(m) <- table$GeneID
    m
}

#' The reference of a product the session roots, by its entity id (a declared dependency parameter).
#' @noRd
trxProductReference <- function(manager, entity_id) {
    reference <- if (is.character(entity_id) && length(entity_id) == 1L) manager$.context$snapshot$products[[entity_id]]
    if (is.null(reference)) cli::cli_abort("the transcriptomics product {.val {entity_id}} is not in the session", class = "workflow_operation_parameters")
    reference
}

#' The value of a product the session roots, by its entity id or its reference.
#' @noRd
trxProductValue <- function(manager, product) {
    reference <- if (isSavedWorkflowReference(product)) product else trxProductReference(manager, product)
    hydrateSavedWorkflowEntity(manager$.context$store, reference)$value
}

#' The RUVs step's product of a state (the corrected matrix, W, k and the logCPM and grid it used).
#' @noRd
trxStateRuvs <- function(manager, reference) {
    product <- readWorkflowStateMetadata(manager, reference)$data$products$trx_ruvs_science
    if (!isSavedWorkflowReference(product)) cli::cli_abort("the state has no RUVs result", class = "workflow_execution_invalid")
    hydrateSavedWorkflowEntity(manager$.context$store, product)$value
}

#' The template of a partitioned variant (KIT095-F05): the stored design and slots with each count table as its zero-row
#' prototype. The class's validity refuses an object with empty counts (a zero library), so the object is built without
#' it, and nothing here calls edgeR (DGEList refuses zero rows).
#' @noRd
trxReadTemplate <- function(manager, source) {
    context <- manager$.context
    spec <- readWorkflowStateMetadata(manager, source)$data
    prototype <- function(component) {
        scope <- newWorkflowTableReadScope(context$store, component, metadata_only = TRUE)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        readSavedWorkflowTablePrototype(scope)
    }
    object <- methods::new("TranscriptomicsAssayData")
    for (slot in intersect(names(spec$slots), methods::slotNames(object))) methods::slot(object, slot, check = FALSE) <- spec$slots[[slot]]
    methods::slot(object, "design_matrix", check = FALSE) <- materializeWorkflowTable(context, spec$components$design_matrix)
    methods::slot(object, "count_data", check = FALSE) <- lapply(spec$components$count_data, prototype)
    attr(object, "dia_component_source") <- source
    attr(object, "dia_execution_fields") <- "design_matrix"
    object
}

#' The design-ordered counts of a page of the count table.
#' @noRd
trxPageCounts <- function(page, runs, id_column) {
    m <- as.matrix(page[, runs, drop = FALSE])
    # as.matrix of a frame with no rows is logical; the template's output types are the counts'.
    if (!nrow(m)) storage.mode(m) <- "double"
    rownames(m) <- page[[id_column]]
    m
}

#' The whole read a partitioned variant's pre-pass makes: the validated inputs of the source state.
#' @noRd
trxSourceInputs <- function(binding, store) {
    object <- readWorkflowExecutionComponents(store$manager, store$source, fields = binding$fields)
    c(trxAssayInputs(object), list(id_column = object@gene_id_column))
}

# ---------------------------------------------------------------------------------------------------------------
# The expression filter (row, given the pre-pass statistics).
# ---------------------------------------------------------------------------------------------------------------

#' The filter's value: the gene counts at each stage and its global statistics.
#' @noRd
trxFilterValue <- function(before, after_expr, after_var, statistics) {
    list(counts = list(before = before, after_filter_by_expr = after_expr, after_zero_variance = after_var),
        effective_lib_size = statistics$effective_lib_size, min_sample_size = statistics$min_sample_size,
        median_lib_size = statistics$median_lib_size)
}

#' The expression filter (rmd:166-175): `filterByExpr` on TMM-scaled library sizes, then the zero-variance rows.
#' @noRd
trxExpressionFilterOperation <- function(binding) {
    label <- "trx-expression-filtered"
    declaration <- newOmicsOperation(binding, "trx_expression_filter", "row", parameters = character(),
        product_type = "trx.expression_filter_result",
        global_statistics = c(effective_lib_size = "lib.size x the TMM factors of the unfiltered counts (a whole read)",
            median_lib_size = "the median effective library size; the CPM cutoff is min.count / median x 1e6",
            min_sample_size = "the smallest group size with filterByExpr's large.n / min.prop adjustment",
            filter_defaults = "filterByExpr's min.count 10, min.total.count 15, large.n 10, min.prop 0.7"),
        kernel = function(object, parameters, store) {
            x <- trxAssayInputs(object)
            f <- trxExpressionFilter(x$counts, x$design$group)
            kept <- which(f$keep_expr)[f$keep_var]
            table <- object@count_data[[1L]][kept, , drop = FALSE]
            rownames(table) <- NULL
            object@count_data[[1L]] <- table
            c(list(final = store$checkpoint(object, label)), trxFilterValue(nrow(x$counts), sum(f$keep_expr), length(kept), f))
        })
    declaration$partitioned <- function(parameters, store) trxPartitionedFilter(binding, store, label)
    declaration
}

#' The filter's partitioned variant: the pre-pass reads the counts once for the effective library sizes (TMM is a whole
#' read); each page is then `filterByExpr` with those sizes fixed, and its zero-variance rows (step map: exact).
#' @noRd
trxPartitionedFilter <- function(binding, store, label) {
    requireOptionalPackage("edgeR", "the transcriptomics expression filter")
    x <- trxSourceInputs(binding, store)
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = x$counts, group = x$design$group), method = "TMM")
    effective <- dge$samples$lib.size * dge$samples$norm.factors
    statistics <- list(effective_lib_size = effective, min_sample_size = trxFilterMinSampleSize(x$design$group),
        median_lib_size = stats::median(effective))
    runs <- colnames(x$counts)
    group <- x$design$group
    x <- dge <- NULL
    kept_expr <- 0L
    page_filter <- function(page) {
        if (!nrow(page)) return(page)
        m <- as.matrix(page[, runs, drop = FALSE])
        keep <- edgeR::filterByExpr(m, group = group, lib.size = effective)
        kept_expr <<- kept_expr + sum(keep)
        variable <- if (any(keep)) apply(m[keep, , drop = FALSE], 1, stats::var) > 0 else logical()
        out <- page[which(keep)[variable], , drop = FALSE]
        rownames(out) <- NULL
        out
    }
    run <- omicsRunPartitioned(binding, store, function(template) template, label,
        write_assay = function(assay, component, prototype, template, partition) partition$write_rows(partition$prototype, page_filter))
    rows <- function(tables) as.integer(tables[[1L]]$metadata$shape[[1L]])
    c(list(final = run$final), trxFilterValue(rows(run$source_tables), kept_expr, rows(run$tables), statistics))
}

# ---------------------------------------------------------------------------------------------------------------
# TMM and logCPM, the RUV grid, RUVs (products and the RUVs state).
# ---------------------------------------------------------------------------------------------------------------

#' TMM factors and the TMM logCPM of the filtered counts (rmd:178-180), a product: the logCPM as a saved table and edgeR's
#' samples table. The raw logCPM is not stored (trx_qc_figures computes it for stage raw).
#' @noRd
trxTmmLogCpmOperation <- function(binding) {
    declaration <- newOmicsOperation(binding, "trx_tmm_logcpm", "whole_matrix", parameters = character(),
        product_type = "trx.tmm_logcpm", output_class = NULL, fixed = list(method = "TMM", prior_count = 2),
        global_statistics = c(lib_size = "the filtered library sizes x the filtered TMM factors (a whole read), prior.count 2"),
        kernel = function(object, parameters, store) {
            x <- trxAssayInputs(object)
            n <- trxTmmLogCpm(x$counts, x$design$group)
            list(log_cpm = store$table(trxMatrixTable(n$log_cpm_tmm)), samples = n$norm_factors_filtered)
        })
    declaration$partitioned <- function(parameters, store) trxPartitionedLogCpm(binding, store)
    declaration
}

#' The logCPM's partitioned variant: the TMM factors in the pre-pass, then the logCPM table written page by page with the
#' library sizes fixed (cpm.default is per cell given them; step map).
#' @noRd
trxPartitionedLogCpm <- function(binding, store) {
    manager <- store$manager
    context <- manager$.context
    x <- trxSourceInputs(binding, store)
    samples <- trxFilteredDge(x$counts, x$design$group)$samples
    lib_size <- samples$lib.size * samples$norm.factors
    runs <- colnames(x$counts)
    id_column <- x$id_column
    x <- NULL
    tables <- readWorkflowStateMetadata(manager, store$source)$data$components$count_data
    component <- tables[[1L]]
    scope <- newWorkflowTableReadScope(context$store, component, reader = if (component$metadata$shape[[2L]] > 64L) "duckdb" else "arrow")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    partition <- omicsAssayPartition(scope, context, component, "trx-logcpm", names(tables)[[1L]], store$worker_cap %||% (12 * 1024^3))
    page_cpm <- function(page) trxMatrixTable(edgeR::cpm(trxPageCounts(page, runs, id_column), lib.size = lib_size, log = TRUE, prior.count = 2))
    table <- store$remember(partition$write_rows(suppressWarnings(page_cpm(partition$prototype)), page_cpm))
    list(log_cpm = table, samples = samples)
}

#' The RUV grid (rmd:195-298): automatic, the notebook's grid and its first winner; manual, the controls ranked at
#' `manual_percentage` with `manual_k`. A product: the grid, its winner and the chosen k, percentage and controls.
#' @noRd
trxRuvGridOperation <- function(binding) {
    newOmicsOperation(binding, "trx_ruv_grid", "whole_matrix", parameters = c("logcpm", trxRuvGridParameterNames),
        product_type = "trx.ruv_grid", output_class = NULL,
        fixed = list(max_k = "3 below 15 samples, else 4 (the kernel's rule)",
            controls = "inferred: the least significant genes of the preliminary ~ group QL GLM, coefficient 2"),
        kernel = function(object, parameters, store) {
            x <- trxAssayInputs(object)
            trxChooseRuv(x, trxSavedMatrix(store$manager, trxProductValue(store$manager, parameters$logcpm)$log_cpm), parameters)
        })
}

#' Run the grid in the parameters' mode and choose k, the percentage and the controls.
#' @return trxRuvGrid's list with `mode`, `settings` and `chosen`.
#' @noRd
trxChooseRuv <- function(x, log_cpm, parameters) {
    mode <- parameters$ruv_mode
    if (!identical(mode, "automatic") && !identical(mode, "manual")) {
        trxAbort(sprintf("RUV mode '%s' is not supported (automatic or manual)", paste(mode, collapse = ",")), c("trx_ruv_mode_unsupported", "trx_error"))
    }
    manual <- identical(mode, "manual")
    if (manual && (!is.numeric(parameters$manual_k) || length(parameters$manual_k) != 1L || is.na(parameters$manual_k) || parameters$manual_k < 1 ||
        !is.numeric(parameters$manual_percentage) || length(parameters$manual_percentage) != 1L || is.na(parameters$manual_percentage))) {
        trxAbort("the manual mode needs one manual_k of at least 1 and one manual_percentage", c("trx_ruv_manual_invalid", "trx_error"))
    }
    range <- if (manual) as.numeric(parameters$manual_percentage) else
        seq(as.numeric(parameters$percentage_min), as.numeric(parameters$percentage_max), by = as.numeric(parameters$percentage_step))
    grid <- trxRuvGrid(x$counts, log_cpm, x$design, list(percentage_range = range, k_penalty_weight = parameters$k_penalty_weight))
    k <- if (manual) as.integer(parameters$manual_k) else grid$best$best_k
    if (k > grid$max_acceptable_k) {
        trxAbort(sprintf("manual k %d is above the cap of %d for %d samples", k, grid$max_acceptable_k, ncol(x$counts)),
            c("trx_ruv_k_above_cap", "trx_error"))
    }
    c(grid, list(mode = mode, settings = parameters[trxRuvGridParameterNames],
        chosen = list(k = k, percentage = grid$best$percentage, n_controls = grid$best$n_controls, control_genes = grid$control_genes)))
}

#' The record's RUV outcomes from the grid product (KIT095-F11): the names KIT-097's template and report read.
#' @noRd
trxRuvOutcomes <- function(grid) {
    list(ruv_mode = grid$mode, percentage_min = grid$settings$percentage_min, percentage_max = grid$settings$percentage_max,
        percentage_step = grid$settings$percentage_step, k_penalty_weight = grid$settings$k_penalty_weight,
        max_acceptable_k = grid$max_acceptable_k, k = grid$chosen$k, negative_control_percentage = grid$chosen$percentage,
        negative_control_count = grid$chosen$n_controls, controls = "inferred",
        control_rule = "least significant genes in the preliminary ~ group QL GLM, coefficient 2",
        separation_score = grid$best$separation_score, composite_score = grid$best$composite_score,
        failed_members = as.numeric(grid$grid_failures$percentage), ruv_applied = TRUE)
}

#' RUVs (rmd:406-426): the state keeps the counts and gains W_1..W_k in its design (the notebook's `ruvs_pdata`); the
#' removeBatchEffect display matrix is a saved table of the step's product.
#' @noRd
trxRuvsOperation <- function(binding) {
    newOmicsOperation(binding, "trx_ruvs", "whole_matrix", parameters = c("logcpm", "grid"), product_type = "trx.ruvs_result",
        kernel = function(object, parameters, store) {
            x <- trxAssayInputs(object)
            grid <- trxProductValue(store$manager, parameters$grid)
            log_cpm <- trxSavedMatrix(store$manager, trxProductValue(store$manager, parameters$logcpm)$log_cpm)
            s <- trxRuvs(x$counts, log_cpm, x$design, grid$chosen$control_genes, grid$chosen$k)
            for (w in names(s$W)) object@design_matrix[[w]] <- unname(s$W[[w]])
            final <- store$checkpoint(object, "trx-ruvs")
            # The products it used are kept as references: an entity id is run-specific, a reference is compared by content.
            list(final = final, corrected = store$table(trxMatrixTable(s$log_cpm_ruvs)), W = s$W, k = grid$chosen$k,
                control_source = s$control_source, logcpm = trxProductReference(store$manager, parameters$logcpm),
                grid = trxProductReference(store$manager, parameters$grid), outcomes = trxRuvOutcomes(grid))
        })
}

# ---------------------------------------------------------------------------------------------------------------
# Differential expression: the shared contrast declaration on KIT-093's fit and test halves.
# ---------------------------------------------------------------------------------------------------------------

#' The DE kernel's sample table: the validated design (Run, group, W_1..W_k) with the runs as row names, as `pData`.
#' @noRd
trxDaDesign <- function(object) {
    x <- trxAssayInputs(object)
    data.frame(x$design, row.names = x$design$Run, check.names = FALSE)
}

#' The realised DE formula of a design (rmd:529-533): `~ 0 + group + W_1 + ... + W_k`.
#' @param design A design (its column names are read).
#' @noRd
trxDaFormula <- function(design) {
    w_terms <- colnames(design)[startsWith(colnames(design), "W_")]
    paste("~ 0 + group", paste(w_terms, collapse = " + "), sep = " + ")
}

#' The realised formula of a state, from its design alone.
#' @noRd
trxRealisedFormula <- function(manager, source) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the transcriptomics DA needs the session's state owner", class = "saved_entity_scope")
    trxDaFormula(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager, source)$design_matrix)
}

#' The DE design (rmd:532-536) with its rank check.
#' @noRd
trxDaModelMatrix <- function(pdata, formula) {
    design_final <- stats::model.matrix(stats::as.formula(formula), data = pdata)
    colnames(design_final) <- make.names(colnames(design_final))
    if (qr(design_final)$rank < ncol(design_final)) {
        trxAbort("the final design is not of full column rank (W confounded with group)", c("trx_da_rank_deficient", "trx_error"))
    }
    design_final
}

#' The fit half of trxDaContrast (R/func_trx_da.R:50-61): TMM, the realised design, the robust dispersion and QL fit.
#' @noRd
trxDaFit <- function(counts, pdata, formula) {
    validateTrxCounts(counts)
    requireOptionalPackage("edgeR", "the transcriptomics differential expression")
    dge_final <- edgeR::DGEList(counts = counts, group = pdata$group)
    dge_final <- edgeR::calcNormFactors(dge_final, method = "TMM")
    design_final <- trxDaModelMatrix(pdata, formula)
    dge_final <- edgeR::estimateDisp(dge_final, design_final, robust = TRUE)
    edgeR::glmQLFit(dge_final, design_final, robust = TRUE)
}

#' The test half of trxDaContrast (R/func_trx_da.R:62-63) for one contrast entry, in the shape the shared declaration reads.
#' @noRd
trxDaApply <- function(fit, contrast_entry) {
    parsed <- trxParseContrastEntry(contrast_entry)
    contrasts_matrix <- trxContrastMatrix(parsed, fit$design)
    lrt <- edgeR::glmQLFTest(fit, contrast = contrasts_matrix[, parsed$name])
    top_genes <- edgeR::topTags(lrt, n = Inf, sort.by = "PValue")$table
    list(results = stats::setNames(list(list(top_table = top_genes, qlf = lrt, contrasts_matrix = contrasts_matrix)), contrast_entry),
        fit.eb = lrt, qvalue_warnings = list())
}

#' The preflight: the formula is the head's realised one, the design has full rank and the contrast is estimable.
#' @noRd
trxDaPreflight <- function(object, contrasts, formula) {
    pdata <- trxDaDesign(object)
    if (!identical(formula, trxDaFormula(pdata))) {
        trxAbort(sprintf("the formula '%s' is not the head's realised one", formula), c("trx_da_formula_mismatch", "trx_error"))
    }
    # KIT-093 revision 1.2 (any number of groups, pairwise contrasts): every contrast line is checked, not only the first.
    model <- trxDaModelMatrix(pdata, formula)
    for (entry in as.character(unlist(contrasts$contrasts))) trxContrastMatrix(trxParseContrastEntry(entry), model)
    invisible(NULL)
}

#' The contrast's long table: the annotated topTags table (rmd:558-562) with the comparison and the significance class.
#' @noRd
trxDaFormat <- function(results, contrast, assay, friendly_name, annotation_info, expr_matrix, design, sample_id, group_id, qvalue, lfc) {
    table <- trxAnnotateDe(results$top_table, annotation_info, friendly_name)$de_annotated
    table$comparison <- friendly_name
    table$significant <- ifelse(table$FDR < qvalue, ifelse(table$logFC > 0, "Up", "Down"), "NS")
    table
}

#' The significant counts of a long table.
#' @noRd
trxDaCount <- function(table) {
    if (is.null(table) || !nrow(table)) return(list(up = 0L, down = 0L, ns = 0L))
    list(up = sum(table$significant == "Up"), down = sum(table$significant == "Down"), ns = sum(table$significant == "NS"))
}

#' The model record of the manifest (decision 2): edgeR QL, as the kernel runs it.
#' @noRd
trxDaModel <- function(parameters) {
    list(model = "edgeR_QL", dispersion = "estimateDisp(robust = TRUE)", fit = "glmQLFit(robust = TRUE)", test = "QL F-test (glmQLFTest)",
        adjust = "BH", formula_used = parameters$formula,
        w_columns = regmatches(parameters$formula, gregexpr("W_[0-9]+", parameters$formula))[[1L]], treat_applied = FALSE)
}

#' The contrasts table of contrast entries (`name = expression`, trxParseContrastFile's).
#' @noRd
trxContrastsTable <- function(entries) {
    data.frame(contrasts = entries, friendly_names = vapply(entries, function(entry) trxParseContrastEntry(entry)$name, ""),
        stringsAsFactors = FALSE, row.names = NULL)
}

# ---------------------------------------------------------------------------------------------------------------
# The RUVIII export after DA (decision 4): a product, accepted on the head as the DA manifest is.
# ---------------------------------------------------------------------------------------------------------------

#' RUVIII on the TMM logCPM (rmd:704-712) with `M` from the group or the technical replicates.
#' @noRd
trxRuviiiExportOperation <- function(binding) {
    newOmicsOperation(binding, "trx_ruviii_export", "whole_matrix", parameters = c("logcpm", "grid", "ruviii_m_source"),
        product_type = "trx.ruviii", output_class = NULL,
        kernel = function(object, parameters, store) {
            x <- trxAssayInputs(object)
            grid <- trxProductValue(store$manager, parameters$grid)
            log_cpm <- trxSavedMatrix(store$manager, trxProductValue(store$manager, parameters$logcpm)$log_cpm)
            r <- trxRuviii(log_cpm, x$design, grid$chosen$control_genes, grid$chosen$k, m_source = parameters$ruviii_m_source,
                technical_replicate_id = object@technical_replicate_id)
            c(r, list(control_count = sum(r$ctl_mask)))
        })
}

#' Accept a RUVIII export on the head (the worker's or the wrapper's run): the product rooted, the step's provenance cleared
#' (it commits no state) and the export recorded under `snapshot$exports`, which the analysis record lists.
#' @param run The operation's run (`reference`, `value`).
#' @noRd
acceptTrxRuviiiExport <- function(manager, source, run) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the RUVIII export needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        clearWorkflowStepProvenance(manager)
        if (is.null(context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(manager, run$reference, list(source))
        snapshot <- context$snapshot
        snapshot$exports$trx_ruviii_export <- list(operation = "trx_ruviii_export", source = source, reference = run$reference,
            parameters = list(ruviii_m_source = run$value$m_source))
        publishWorkflowStateHead(context, snapshot)
        invisible(run)
    })
}

# ---------------------------------------------------------------------------------------------------------------
# The steps, the record and the replay.
# ---------------------------------------------------------------------------------------------------------------

#' The transcriptomics steps (the binding's `steps`, KIT-094), in registry order.
#' @noRd
trxRouteSteps <- function(binding) {
    list(trxExpressionFilterOperation(binding), omicsStoreOperation(binding, "trx_qc_finalize"), trxTmmLogCpmOperation(binding),
        trxRuvGridOperation(binding), trxRuvsOperation(binding), omicsDaContrastOperation(binding), trxRuviiiExportOperation(binding),
        omicsStoreOperation(binding, "trx_reset"))
}

#' The realised outcomes the record keeps for a transcriptomics step or export (the binding's `record_outcomes`): the
#' filter's gene counts, the RUV outcomes, and the RUVIII export's M source, k and control count.
#' @noRd
trxRouteRecordOutcomes <- function(binding, context, operation, record, outcomes) {
    product <- record$data$products[[paste0(operation, "_science")]]
    if (!isSavedWorkflowReference(product)) return(outcomes)
    value <- tryCatch(hydrateSavedWorkflowEntity(context$store, product)$value, error = function(error) NULL)
    if (is.null(value)) return(outcomes)
    if (identical(operation, "trx_expression_filter")) {
        outcomes[c("genes_before", "genes_after_filter_by_expr", "genes_after_zero_variance")] <- unname(value$counts)
    } else if (identical(operation, "trx_ruvs")) {
        outcomes[names(value$outcomes)] <- value$outcomes
    } else if (identical(operation, "trx_ruviii_export")) {
        outcomes[c("ruviii_m_source", "k", "negative_control_count")] <- list(value$m_source, value$k, value$control_count)
    }
    outcomes
}

#' The replay of the transcriptomics steps and export (the binding's `replay`), by the name after the prefix.
#' @noRd
trxRouteReplay <- local({
    store_step <- function(binding, manager, step, parameters, config, description) {
        commitWorkflowOmicsStoreStep(manager, step$operation, step$state, config, description)
    }
    list(
        expression_filter = function(binding, manager, step, parameters, config, description) {
            runTrxReferenceExpressionFilter(manager, config, description, state_name = step$state)
        },
        qc_finalize = store_step,
        reset = store_step,
        ruvs = function(binding, manager, step, parameters, config, description) {
            # The grid that ran last before the commit is the one the state used.
            previous <- Filter(function(p) identical(p$operation, "trx_ruv_grid"), step$preceded_by %||% list())
            if (!length(previous)) omicsReplayUnsupported(step)
            grid <- workflowReplayParameters(previous[[length(previous)]]$parameters, workflowOperationEntry("trx_ruv_grid")$parameters)
            runTrxReferenceRuvs(manager, grid[trxRuvGridParameterNames], config, description, state_name = step$state)
        },
        ruviii_export = function(binding, manager, step, parameters, config, description) {
            runTrxReferenceRuviiiExport(manager, parameters$ruviii_m_source)
        })
})

#' The replay's check of a recorded DA formula against the replayed head's (the binding's `da_replay_check`).
#' @noRd
trxDaReplayCheck <- function(manager, settings) {
    derived <- trxRealisedFormula(manager, getWorkflowStateReference(manager))
    if (!identical(settings$formula_string, derived)) {
        cli::cli_abort("the record's DA formula {.val {settings$formula_string}} is not the replayed head's {.val {derived}}",
            class = "workflow_replay_formula_mismatch")
    }
    invisible(NULL)
}

# ---------------------------------------------------------------------------------------------------------------
# The worker's publication operations: KIT-093's renderers on the products, the notebook's file names.
# ---------------------------------------------------------------------------------------------------------------

#' The design without the RUVs covariates: the figures join the notebook's design (Run, group and its own columns).
#' @noRd
trxFigureDesign <- function(x) x$design[, !startsWith(names(x$design), "W_"), drop = FALSE]

#' The RUVs display matrix from a RUVs state's inputs (rmd:421-426): the same removeBatchEffect call on the same TMM logCPM
#' and W, for the shared final PCA, which reads the state alone.
#' @noRd
trxStateCorrectedLogCpm <- function(x) {
    pdata <- data.frame(x$design, row.names = x$design$Run, check.names = FALSE)
    w_matrix <- pdata[, startsWith(names(pdata), "W_"), drop = FALSE]
    log_cpm_tmm <- trxTmmLogCpm(x$counts, x$design$group)$log_cpm_tmm
    limma::removeBatchEffect(log_cpm_tmm, covariates = w_matrix, design = stats::model.matrix(~ group, data = pdata))
}

#' One normalisation stage's matrix (genes by runs) for a RUVs state: raw and TMM logCPM, the RUVs display matrix, the
#' RUVIII matrix of the accepted export.
#' @noRd
trxStageMatrix <- function(manager, source, stage, x) {
    switch(stage,
        raw = trxTmmLogCpm(x$counts, x$design$group)$log_cpm_raw,
        tmm = trxSavedMatrix(manager, trxProductValue(manager, trxStateRuvs(manager, source)$logcpm)$log_cpm),
        ruvs = trxSavedMatrix(manager, trxStateRuvs(manager, source)$corrected),
        ruviii = {
            export <- manager$.context$snapshot$exports$trx_ruviii_export
            if (is.null(export)) cli::cli_abort("no accepted RUVIII export", class = "workflow_execution_invalid")
            t(hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value$newY)
        },
        cli::cli_abort("unknown normalisation stage {.val {stage}}", class = "workflow_execution_invalid"))
}

#' Save a figure as the notebook does (ggsave, PNG and PDF).
#' @noRd
trxSaveFigure <- function(plot, directory, name, width, height) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    for (extension in c(".png", ".pdf")) {
        ggplot2::ggsave(file.path(directory, paste0(name, extension)), plot, width = width, height = height, dpi = 300)
    }
    invisible(NULL)
}

#' A publication's completion value: the source's checkpoint, the display; a failure is the display's message.
#' @noRd
trxPublicationResult <- function(manager, request, produce) {
    source <- request$input_references$state
    display <- tryCatch(produce(source), error = function(error) list(saved = FALSE, message = substr(conditionMessage(error), 1L, 1024L)))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source, display = display)
}

#' The grid's publication (rmd:283-305, 373-393, 429-432): the CCA panels, the workbook, the RUVs display matrix TSV.
#' @noRd
computeWorkflowTrxRuvPublicationAction <- function(manager, request, variant) {
    trxPublicationResult(manager, request, function(source) {
        ruvs <- trxStateRuvs(manager, source)
        grid <- trxProductValue(manager, ruvs$grid)
        panels <- trxCcaPanels(grid$best, grid$grid_summary, length(grid$prelim_order))
        graphs <- omicsStagedPublicationDir(request, "publication_graphs_dir")
        trxSaveFigure(panels$grid$plot, graphs, "canonical_correlation_panels", 12, 5)
        wb <- openxlsx::createWorkbook()
        for (sheet in names(panels$workbook)) {
            openxlsx::addWorksheet(wb, sheet)
            openxlsx::writeData(wb, sheet, panels$workbook[[sheet]])
        }
        openxlsx::saveWorkbook(wb, file.path(graphs, "canonical_correlation_optimization_summary.xlsx"), overwrite = TRUE)
        vroom::vroom_write(as.data.frame(trxSavedMatrix(manager, ruvs$corrected)) |> tibble::rownames_to_column("GeneID"),
            file.path(omicsStagedPublicationDir(request, "gene_qc_dir"), "ruv_normalised_log_values.tsv"))
        list(saved = TRUE)
    })
}

#' One stage's PCA and RLE (rmd:441-498): the statistics saved for the tab, the two figures staged.
#' @noRd
computeWorkflowTrxQcFiguresAction <- function(manager, request, variant) {
    trxPublicationResult(manager, request, function(source) {
        stage <- request$parameters$stage
        x <- trxAssayInputs(readWorkflowOmicsState(transcriptomicsDeclarationBinding, manager, source))
        mat <- trxStageMatrix(manager, source, stage, x)
        design <- trxFigureDesign(x)
        stats <- list(pca = trxPcaStats(mat, design), rle = trxRleStats(mat, design))
        directory <- omicsStagedPublicationDir(request, "gene_qc_dir")
        ggplot2::ggsave(file.path(directory, paste0("trx_", stage, "_pca.png")), trxPcaPlot(stats$pca, paste("PCA:", stage)), width = 7, height = 6, dpi = 150)
        ggplot2::ggsave(file.path(directory, paste0("trx_", stage, "_rle.png")), trxRlePlot(stats$rle, paste("RLE:", stage)), width = 9, height = 6, dpi = 150)
        list(saved = TRUE, stats = saveWorkflowStateValue(manager$.context, stats, "trx-qc-stats", codec = "rds"))
    })
}

#' The notebook's composite figures: the normalisation-stage grid (rmd:490-511) or the RUVIII before and after (rmd:730-745).
#' @noRd
computeWorkflowTrxCompositeFigureAction <- function(manager, request, variant) {
    trxPublicationResult(manager, request, function(source) {
        x <- trxAssayInputs(readWorkflowOmicsState(transcriptomicsDeclarationBinding, manager, source))
        stage <- function(name) trxStageMatrix(manager, source, name, x)
        graphs <- omicsStagedPublicationDir(request, "publication_graphs_dir")
        if (identical(request$parameters$figure, "ruviii_before_and_after")) {
            grid <- trxRuviiiGrid(stage("tmm"), t(stage("ruviii")), trxFigureDesign(x))
            trxSaveFigure(grid$grid$plot, graphs, "ruviii_before_and_after", 11, 8)
        } else {
            grid <- trxNormStageGrid(stage("raw"), stage("tmm"), stage("ruvs"), trxFigureDesign(x), trxStateRuvs(manager, source)$k)
            trxSaveFigure(grid$grid$plot, graphs, "normalization_comparison_grid", 15, 9)
        }
        list(saved = TRUE)
    })
}

#' The DE publication (rmd:558-583, 597-664, 755-775) from the accepted manifest: per contrast the DE CSV, the volcano and
#' the significance bar; the top genes of every contrast; the summary text of the first contrast, the dataset line the
#' project's (`parameters$dataset`).
#' @noRd
computeWorkflowTrxDaPublicationAction <- function(manager, request, variant) {
    trxPublicationResult(manager, request, function(source) {
        context <- manager$.context
        if (is.null(context$snapshot$da_results)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
        manifest <- resolveWorkflowStateValue(context, context$snapshot$da_results)
        object <- readWorkflowOmicsState(transcriptomicsDeclarationBinding, manager, manifest$source)
        da_dir <- omicsStagedPublicationDir(request, "da_output_dir")
        graphs <- omicsStagedPublicationDir(request, "publication_graphs_dir")
        annotated <- list()
        for (contrast in names(manifest$results)) {
            entry <- manifest$results[[contrast]]
            if (is.null(entry$products$da_contrast)) next
            results <- hydrateSavedWorkflowEntity(context$store, entry$products$da_contrast)$value$assays[[1L]]$results
            name <- entry$comparison
            de <- trxAnnotateDe(results$top_table, object@count_data[[1L]], name, object@gene_id_column, object@annotation_id_column)
            utils::write.csv(de$de_annotated, file.path(da_dir, paste0("de_results_", name, ".csv")), row.names = FALSE)
            trxSaveFigure(trxVolcano(de$de_annotated, name)$plot, file.path(graphs, "Volcano_Plots"), paste0("volcano_plot_", name), 8, 7)
            trxSaveFigure(trxSignificanceBar(de$de_annotated, name)$plot, file.path(graphs, "NumSigDeMolecules"), paste0("summary_barplot_", name), 7, 5)
            annotated[[contrast]] <- de
        }
        if (!length(annotated)) cli::cli_abort("no contrast has a result", class = "workflow_execution_invalid")
        vroom::vroom_write(dplyr::bind_rows(lapply(annotated, `[[`, "top_de")), file.path(da_dir, "top_de_genes.tsv"))
        grid <- trxProductValue(manager, trxStateRuvs(manager, manifest$source)$grid)
        trxWriteSummaryText(trxSummaryText(grid$best, length(grid$prelim_order), manifest$parameters$formula, names(annotated)[[1L]],
            annotated[[1L]]$de_annotated, grid$grid_summary), request$parameters$dataset, file.path(da_dir, "ruv_optimization_summary.txt"))
        list(saved = TRUE)
    })
}

#' Write the summary text as the notebook does (rmd:756-775), the dataset line the project's.
#' @noRd
trxWriteSummaryText <- function(summary, dataset, path) {
    lines <- summary$lines
    lines[[4L]] <- paste("Dataset:", dataset %||% "not recorded")
    writeLines(lines, path)
    utils::write.table(summary$grid_table, file = path, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
    invisible(path)
}

#' The integration export's writer (the binding's `integration_export`): the accepted RUVIII matrix, samples by genes with
#' a `Run` column, as the notebook writes it for MOFA (rmd:714-724).
#' @noRd
trxWriteIntegrationExport <- function(manager, parameters, path) {
    export <- manager$.context$snapshot$exports$trx_ruviii_export
    if (is.null(export)) cli::cli_abort("no accepted RUVIII export", class = "workflow_execution_invalid")
    newY <- hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value$newY
    vroom::vroom_write(as.data.frame(newY) |> tibble::rownames_to_column("Run"), path)
    invisible(path)
}

#' The interactive volcano of one contrast (rmd:673-693; the binding's `volcano_widget`): Glimma on the contrast's DGELRT
#' with the annotation aligned to its order. A failure is the page action's recorded message.
#' @noRd
trxVolcanoWidget <- function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh = 0.05) {
    entry <- da_results_list$contrasts_results[[selected_assay]]$results[[selected_contrast]]
    if (is.null(entry)) return(NULL)
    name <- trxParseContrastEntry(selected_contrast)$name
    object <- da_results_list$theObject
    de <- trxAnnotateDe(entry$top_table, object@count_data[[1L]], name, object@gene_id_column, object@annotation_id_column)$de_annotated
    inputs <- trxGlimmaInputs(de, entry$qlf, name)
    Glimma::glimmaVolcano(entry$qlf, main = inputs$main, anno = inputs$anno, display.columns = inputs$display_columns)
}

#' The tab's final PCA (the binding's `final_pca`) on the RUVs display matrix.
#' @noRd
trxFinalPca <- function(object, color = "group", shape = NULL) {
    x <- trxAssayInputs(object)
    k <- sum(startsWith(names(x$design), "W_"))
    trxPcaPlot(trxPcaStats(trxStateCorrectedLogCpm(x), trxFigureDesign(x)), sprintf("PCA: After RUVs (k=%d)", k), group_col = color %||% "group")
}

#' The filtering metrics of a state (the binding's `filtering_metrics`): genes, and per run the library size and the genes
#' detected.
#' @noRd
trxFilteringMetrics <- function(object) {
    x <- trxAssayInputs(object)
    list(genes = nrow(x$counts), library_size = colSums(x$counts), detected_genes = colSums(x$counts > 0))
}

#' The DE kernel's inputs from the count table: the counts in the design's run order and the gene annotation.
#' @noRd
trxDaAssayInput <- function(table, assay, design, sample_id, id_column, annotation_id_column) {
    list(expr_matrix = trxPageCounts(table, design$Run, id_column),
        annotation_info = data.frame(GeneID = table[[id_column]], GeneName = table[[annotation_id_column]], stringsAsFactors = FALSE))
}

#' The worker's publication operations of the transcriptomics route (the binding's `publications`).
#' @noRd
trxRoutePublications <- list(c("filtering_metrics", "computeWorkflowOmicsFilteringMetricsAction"),
    c("ruv_publication", "computeWorkflowTrxRuvPublicationAction"), c("qc_figures", "computeWorkflowTrxQcFiguresAction"),
    c("composite_figure", "computeWorkflowTrxCompositeFigureAction"), c("da_publication", "computeWorkflowTrxDaPublicationAction"),
    c("da_tab_results", "computeWorkflowOmicsDaTabResultsAction"), c("final_pca", "computeWorkflowOmicsFinalPcaAction"),
    c("volcano_page", "computeWorkflowOmicsVolcanoPageAction"), c("session_export", "computeWorkflowOmicsSessionExportAction"),
    c("integration_export", "computeWorkflowOmicsIntegrationExportAction"))

# ---------------------------------------------------------------------------------------------------------------
# The binding.
# ---------------------------------------------------------------------------------------------------------------

#' The transcriptomics binding of the omics declaration core.
#' @noRd
transcriptomicsDeclarationBinding <- list(
    omic_type = "transcriptomics", prefix = "trx", class = "TranscriptomicsAssayData", data_slot = "count_data",
    id_slot = "gene_id_column", fields = transcriptomicsOperationFields, merge = transcriptomicsPartitionMerge,
    design_field = "design_matrix", feature_noun = "gene", long_table = "da_genes_long", design_input_prefix = "trx_counts_",
    step_scope = "the transcriptomics step", da_scope = "the transcriptomics DA",
    step_labels = c(trx_expression_filter = "2_Expression_Filtered", trx_qc_finalize = "3_QC_Complete", trx_ruvs = "4_RUVs",
        trx_reset = "0_Reset"),
    steps = trxRouteSteps,
    worker_steps = c("trx_expression_filter", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_da_contrast", "trx_ruviii_export"),
    # Every step defaults to `whole` until KIT-099 decides.
    default_variants = character(),
    publications = trxRoutePublications,
    filtering_steps = "trx_expression_filter",
    native_steps = "trx_da_contrast",
    dependency_parameters = c("logcpm", "grid"),
    replay = trxRouteReplay,
    record_outcomes = trxRouteRecordOutcomes,
    read_template = trxReadTemplate,
    # The shared DA declaration's edgeR QL settings, design, halves and model record (decision 2).
    da_fixed = list(model = "edgeR_QL", dispersion_robust = TRUE, fit_robust = TRUE, adjust = "BH", treat = FALSE),
    da_design = trxDaDesign,
    da_preflight = trxDaPreflight,
    da_assay_input = trxDaAssayInput,
    da_fit = function(data, design, formula_string, sample_id) trxDaFit(data, design, formula_string),
    da_apply = function(fit, contrast, ...) trxDaApply(fit, contrast),
    da_format = trxDaFormat,
    da_count = trxDaCount,
    da_model = trxDaModel,
    da_replay_check = trxDaReplayCheck,
    qc_dir_role = "gene_qc_dir",
    filtering_metrics = trxFilteringMetrics,
    final_pca = trxFinalPca,
    volcano_widget = trxVolcanoWidget,
    integration_export = trxWriteIntegrationExport
)

#' Identify the transcriptomics DA manifest.
#' @noRd
isTrxResultManifest <- function(value) is.list(value) && isTRUE(value[["schema", exact = TRUE]] %in% "trx.da_results.v2")

# ---------------------------------------------------------------------------------------------------------------
# The wrappers: what KIT-096/097's tabs and the record's replay call.
# ---------------------------------------------------------------------------------------------------------------

#' The transcriptomics state read for one action (KIT-075).
#' @noRd
readTrxWorkflowState <- function(manager, reference = NULL) readWorkflowOmicsState(transcriptomicsDeclarationBinding, manager, reference)

#' The design step on the owner: the `design_builder` provenance and the inputs (`trx_counts_<n>`, `design_file`).
#' @noRd
recordTrxDesignProvenance <- function(manager, object, source_dir = NULL, design_file = NULL) {
    recordWorkflowOmicsDesignProvenance(transcriptomicsDeclarationBinding, manager, object, source_dir, design_file)
}

#' The expression filter's tab step.
#' @noRd
trxExpressionFilterStep <- function(description = NULL, state_name = NULL) {
    list(operation = "trx_expression_filter", parameters = list(), state_name = state_name %||% "trx_expression_filtered",
        description = description %||% "Applied the expression filter (filterByExpr) and removed zero-variance genes")
}

#' Run the expression filter on the head and commit it.
#' @param variant `whole` (the default) or `partitioned`.
#' @noRd
runTrxReferenceExpressionFilter <- function(manager, config, description = NULL, state_name = NULL, variant = NULL) {
    step <- trxExpressionFilterStep(description, state_name)
    if (!is.null(variant)) step$parameters$variant <- variant
    runWorkflowOmicsReferenceStep(manager, step, config)
}

#' A product-only step on a source, rooted in the head (a hit from a store whose head does not root it).
#' @noRd
trxRunProduct <- function(manager, operation, source, parameters, dependencies = list()) {
    run <- runSavedOperation(manager, workflowOperationEntry(operation)$declaration, source, parameters, dependencies)
    if (is.null(manager$.context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(manager, run$reference, c(list(source), dependencies))
    run
}

#' TMM/logCPM, the grid and RUVs on the head, the RUVs state committed (the products before it are its `preceded_by`).
#' @param grid_parameters trxRuvGridParameters' list.
#' @return The RUVs step's run, reference and value, with `logcpm` and `grid`, the products' runs.
#' @noRd
runTrxReferenceRuvs <- function(manager, grid_parameters, config, description = NULL, state_name = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the transcriptomics RUVs needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        logcpm <- trxRunProduct(manager, "trx_tmm_logcpm", source, list())
        grid <- trxRunProduct(manager, "trx_ruv_grid", source, c(list(logcpm = logcpm$reference$entity_id), grid_parameters[trxRuvGridParameterNames]),
            list(logcpm$reference))
        step <- list(operation = "trx_ruvs", parameters = list(logcpm = logcpm$reference$entity_id, grid = grid$reference$entity_id),
            state_name = state_name %||% "trx_ruvs", description = description %||% sprintf("RUVs: W_1..W_%d added to the design", grid$value$chosen$k))
        result <- runWorkflowOmicsReferenceStep(manager, step, config, source = source, dependencies = list(logcpm$reference, grid$reference))
        c(result, list(logcpm = logcpm, grid = grid))
    })
}

#' The DA on a RUVs state, the formula realised from its design (one declaration per contrast on one saved fit).
#' @param settings `da_q_val_thresh` (FDR threshold) and `treat_lfc_cutoff` (the display threshold; no treat).
#' @noRd
runTrxReferenceDaAnalysis <- function(manager, source, settings, contrasts) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the transcriptomics DA needs the session's state owner", class = "saved_entity_scope")
    settings$formula_string <- withWorkflowStateScope(manager, function() trxRealisedFormula(manager, source))
    runWorkflowOmicsReferenceDaAnalysis(transcriptomicsDeclarationBinding, manager, source, settings, contrasts)
}

#' The RUVIII export on the head after DA, on the logCPM and grid its RUVs step used, accepted on the head.
#' @param m_source `group` (the notebook) or `technical_replicate`.
#' @noRd
runTrxReferenceRuviiiExport <- function(manager, m_source = "group") {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the RUVIII export needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        ruvs <- trxStateRuvs(manager, source)
        run <- runSavedOperation(manager, workflowOperationEntry("trx_ruviii_export")$declaration, source,
            list(logcpm = ruvs$logcpm$entity_id, grid = ruvs$grid$entity_id, ruviii_m_source = m_source), list(ruvs$logcpm, ruvs$grid))
        acceptTrxRuviiiExport(manager, source, run)
        run
    })
}
