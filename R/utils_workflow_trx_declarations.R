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

#' The declared parameters of a RUV design other than today's grid (KIT-129 W4), after the dependency `logcpm`.
#' @noRd
trxRuvDesignParameterNames <- c("ruv_design", "ruv_k", "ruv_k_rule", "ruv_control_ids", "ruv_k_seed")

#' The parameters of a declared RUV design (KIT-129 W4).
#' @param ruv_design One of `TRX_RUV_DESIGNS`.
#' @param ruv_k The number of factors, already resolved: a whole number, 0 when the rule found no factor. Not given for
#'   `none`.
#' @param ruv_k_rule The rule that gave `ruv_k`, recorded with it: `num_sv_be` (the estimated count), `sweep_bounded`
#'   (today's grid estimate) or `user`.
#' @param ruv_control_ids The supplied control gene ids, for `supplied_controls`; kept sorted and unique, so one list
#'   is one identity whatever its order.
#' @param ruv_k_seed KIT-129 item 3: the seed of the estimated count, read only when the kernel estimates k (`ruv_k_rule
#'   = "num_sv_be"` with no `ruv_k`); NA there is the default seed.
#' @return The named list the declared `trx_ruv_grid` takes after `logcpm`, in its canonical form
#'   (`trxCanonicalRuvDesignParameters`).
#' @noRd
trxRuvDesignParameters <- function(ruv_design, ruv_k = NA_integer_, ruv_k_rule = NA_character_, ruv_control_ids = character(),
                                   ruv_k_seed = NA_integer_) {
    trxCanonicalRuvDesignParameters(list(ruv_design = ruv_design, ruv_k = ruv_k, ruv_k_rule = ruv_k_rule, ruv_control_ids = ruv_control_ids,
        ruv_k_seed = ruv_k_seed))
}

#' The declared parameters in one form (KIT-129 audit F-07), so a product's identity and its recorded settings do not
#' depend on how the parameters were passed: `ruv_k` an integer (2 and 2L are one k; NA when absent, and for the design
#' `none`, which takes no k), `ruv_k_rule` NA when absent, `ruv_control_ids` `trxCanonicalControlIds`' list,
#' `ruv_k_seed` an integer when the kernel estimates k (the default seed when absent) and NA otherwise. The
#' design and the rule are checked by `trxChooseRuvDesign`, not here.
#' @param parameters A list with the names of `trxRuvDesignParameterNames`; other elements are not returned.
#' @return The five parameters. Refuses `trx_ruv_design_invalid` when `ruv_k` is given and is not one whole number of
#'   0 or more: a fractional k is never rounded.
#' @noRd
trxCanonicalRuvDesignParameters <- function(parameters) {
    k <- parameters$ruv_k
    absent <- is.null(k) || (length(k) == 1L && is.na(k))
    if (!absent && !checkmate::test_count(k)) {
        trxAbort(sprintf("ruv_k must be one whole number of 0 or more, not %s", substr(deparse1(k), 1L, 40L)), c("trx_ruv_design_invalid", "trx_error"))
    }
    none <- identical(parameters$ruv_design, "none")
    rule <- parameters$ruv_k_rule %||% NA_character_
    # KIT-129 item 3: the seed is part of the identity only when the kernel estimates k.
    estimated <- !none && absent && identical(rule, "num_sv_be")
    seed <- parameters$ruv_k_seed
    seed_absent <- is.null(seed) || (length(seed) == 1L && is.na(seed))
    if (!seed_absent && (!estimated || !checkmate::test_int(seed))) {
        trxAbort("ruv_k_seed is one whole number, given only when the kernel estimates k (ruv_k_rule num_sv_be with no ruv_k)",
            c("trx_ruv_design_invalid", "trx_error"))
    }
    list(ruv_design = parameters$ruv_design, ruv_k = if (absent || (none && k == 0)) NA_integer_ else as.integer(round(k)),
        ruv_k_rule = rule, ruv_control_ids = trxCanonicalControlIds(parameters$ruv_control_ids),
        ruv_k_seed = if (estimated) as.integer(if (seed_absent) RUV_FACTOR_COUNT_SEED else seed) else NA_integer_)
}

#' Refuse today's grid when its request carries a declared design's parameter (KIT-129 audit F-02). The declared
#' choice shares the operation's name, so a request with both parameter sets would otherwise run today's grid and drop
#' the declared design without a word. `ruv_k_seed` is the seed parameter of the shared k rule (KIT-129 W6).
#' The guard is the kernel's: a request whose today's product is already saved on the source is a hit of
#' `runSavedOperation` (the identity reads the declared names alone) and is returned without reaching it. Refusing
#' that too is the registration's parameter validation (the later wiring's). The class is the tab's
#' (`trx_ruv_grid_declared_parameters`, first, so its text is unchanged), then `trx_error` and the legacy refusal of
#' every omic, `workflow_execution_invalid` (KIT-129 item 1).
#' @noRd
trxRefuseDeclaredRuvParameters <- function(parameters) {
    declared <- intersect(c(trxRuvDesignParameterNames, "ruv_k_seed"), names(parameters))
    if (!length(declared)) return(invisible(NULL))
    trxAbort(sprintf("today's RUV grid was asked with the declared-design parameter(s) %s: it runs no declared design and does not ignore one",
        paste(declared, collapse = ", ")), c("trx_ruv_grid_declared_parameters", "trx_error", "workflow_execution_invalid"))
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
#' @param ruv_design KIT-129 W4. `"legacy"` (the default): today's declaration and kernel, unchanged. `"declared"`: the
#'   same product chosen by a declared design instead of the grid (`trxRuvDesignParameterNames`, `trxChooseRuvDesign`);
#'   the RUVs step reads either.
#' @noRd
trxRuvGridOperation <- function(binding, ruv_design = c("legacy", "declared")) {
    if (identical(match.arg(ruv_design), "declared")) return(trxRuvDesignOperation(binding))
    newOmicsOperation(binding, "trx_ruv_grid", "whole_matrix", parameters = c("logcpm", trxRuvGridParameterNames),
        product_type = "trx.ruv_grid", output_class = NULL,
        fixed = list(max_k = "3 below 15 samples, else 4 (the kernel's rule)",
            controls = "inferred: the least significant genes of the preliminary ~ group QL GLM, coefficient 2"),
        kernel = function(object, parameters, store) {
            trxRefuseDeclaredRuvParameters(parameters)
            x <- trxAssayInputs(object)
            trxChooseRuv(x, trxSavedMatrix(store$manager, trxProductValue(store$manager, parameters$logcpm)$log_cpm), parameters)
        })
}

#' The choice of a declared RUV design (KIT-129 W4), the product the RUVs step reads in place of the grid's: the
#' controls by the design's rule and the k it was given. `logcpm` stays a declared dependency: the shared rule for k
#' (KIT-129 W6) estimates the count on the TMM logCPM when it is called here.
#' @noRd
trxRuvDesignOperation <- function(binding) {
    newOmicsOperation(binding, "trx_ruv_grid", "whole_matrix", parameters = c("logcpm", trxRuvDesignParameterNames),
        product_type = "trx.ruv_grid", output_class = NULL,
        fixed = list(factors_estimated_on = "upper-quartile normalised counts (EDASeq::betweenLaneNormalization, which = upper)",
            replicate_sets = "the groups (RUVs)", first_pass = "the preliminary ~ group QL GLM, every group coefficient; its most significant 25% left out",
            test = "the filtered counts with W_1..W_k in the edgeR model"),
        kernel = function(object, parameters, store) {
            trxProductReference(store$manager, parameters$logcpm)
            # KIT-129 item 3: the estimated count (`ruv_k_rule = "num_sv_be"` with no `ruv_k`) is computed here, by the
            # shared resolver, on the TMM logCPM (as simulator S3 computes E). KIT-129 integration item 3: this is the run
            # path's one resolution of k, so a request that gives a k under the count is refused (the kernel would take
            # it as resolved); under the user's k the count is recorded beside it.
            canonical <- trxCanonicalRuvDesignParameters(parameters)
            none <- identical(canonical$ruv_design, "none")
            if (!none && !is.na(canonical$ruv_k) && identical(canonical$ruv_k_rule, "num_sv_be")) {
                trxAbort("a k is not given with ruv_k_rule num_sv_be: the estimated count gives k on the run path",
                    c("trx_ruv_design_invalid", "trx_error"))
            }
            estimated <- !none && is.na(canonical$ruv_k) && identical(canonical$ruv_k_rule, "num_sv_be")
            user <- !none && identical(canonical$ruv_k_rule, "user")
            # KIT-129 fifth pass item 4: the separation rule (`ruv_k_rule = "sweep_bounded"` with no `ruv_k`) is computed
            # here too, on the design's own controls (`trxSeparationRuleK`; arms 22 and 24 of simulator S3).
            separation <- !none && is.na(canonical$ruv_k) && identical(canonical$ruv_k_rule, "sweep_bounded")
            x <- trxAssayInputs(object)
            log_cpm <- if (estimated || user || separation) trxSavedMatrix(store$manager, trxProductValue(store$manager, parameters$logcpm)$log_cpm)
            separated <- if (separation) trxSeparationRuleK(x, log_cpm, canonical)
            choice <- trxChooseRuvDesign(x, parameters, k_record = separated$k_record, log_cpm = if (estimated) log_cpm)
            # the plot the rule read, as the classic sweep keeps its own (a grob, without its owners)
            if (separation && inherits(separated$plot, c("gg", "ggplot"))) {
                choice$separation_plot <- ggplot2::ggplotGrob(stripProtNormQcPlotOwners(separated$plot))
            }
            if (user) {
                beside <- recordRuvFactorCountBesideUserK(log_cpm, as.character(x$design$group))
                choice$k_record[c("k_estimate", "k_estimate_reason")] <- beside[c("k_estimate", "k_estimate_reason")]
                choice$k_record$k_seed <- beside$seed
                choice$k_estimate <- beside$estimate
            }
            choice
        })
}

#' The separation rule for k on a declared design's own controls (KIT-129 fifth pass item 4), as simulator S3's arms 22
#' (a supplied list) and 24 (the first pass's controls) compute it: `trxCancor` of the TMM logCPM against the group
#' columns of `~ group` with those controls, then `trxFindBestK` up to the grid's cap (3 below 15 runs, else 4), and the
#' shared resolver's limit (`resolveRuvFactorCount`; undefined, and refused, with all genes as controls).
#' @param x The validated inputs (`counts`, `design`). @param log_cpm Their TMM logCPM.
#' @param parameters `trxCanonicalRuvDesignParameters`' list.
#' @return `list(k_record, plot)`: the resolver's record (the `k_record` of `trxChooseRuvDesign`) and the
#'   canonical-correlation plot the rule read.
#' @noRd
trxSeparationRuleK <- function(x, log_cpm, parameters) {
    trxRuvCheckAlignment(x$counts, x$design, log_cpm)
    controls <- trxRuvDesignControls(x$counts, x$design, parameters$ruv_design, parameters$ruv_control_ids)
    groups <- as.character(x$design$group)
    limit <- ruvFactorCountLimit(ncol(x$counts), length(unique(groups)), length(unique(groups)), length(controls$control_genes))
    plot <- NULL
    k_record <- resolveRuvFactorCount("sweep_bounded", controls$control_rule, limit, separation_rule = function() {
        plot <<- trxCancor(Y = t(log_cpm), X = stats::model.matrix(~ group, data = x$design)[, -1, drop = FALSE],
            ctl = rownames(log_cpm) %in% controls$control_genes)
        trxFindBestK(plot, max_k = if (ncol(x$counts) < 15) 3L else 4L)
    })
    list(k_record = k_record, plot = plot)
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

#' The k rules a declared design records with its k (KIT-129 section 1e; the shared resolver's names).
#' @noRd
TRX_RUV_K_RULES <- c("num_sv_be", "sweep_bounded", "user")

#' The choice of a declared RUV design (KIT-129 W4): its controls, and the k it is given with the rule that gave it.
#' @param x The validated inputs (`counts`, `design`).
#' @param parameters `trxRuvDesignParameters`' list.
#' @param k_record The number of factors and its source: a list with `k` (a whole number; 0 when the rule found no
#'   factor) and `k_rule` (one of `TRX_RUV_K_RULES`), as the shared resolver returns it; its other fields (`reason`,
#'   `k_estimate`, `k_max`, `estimate`, ...) are kept with the product. Not read for `none`. NULL (the default): the
#'   parameters' `ruv_k` and `ruv_k_rule`, or, with `ruv_k_rule = "num_sv_be"` and no `ruv_k`, the shared resolver's
#'   estimated count on `log_cpm` under the limit `ruvFactorCountLimit` gives (KIT-129 item 3; a k above it is refused,
#'   `ruv_k_above_limit`).
#' @param log_cpm The TMM logCPM of the counts, for the estimated count.
#' @return The product the RUVs step reads: `chosen` (`k`, `percentage`, `n_controls`, `control_genes`, the grid
#'   product's names), `ruv_design` (`id`, `method`, `control_rule`, `unmatched`, `list_sha256`, `k`, `k_rule`, `applied`, `reason`,
#'   `test_form`), `k_record`, `k_estimate` (the count's record for `plotRuvFactorCount`, when the count ran), `control_genes`, `control_source`, `mode` (`"declared"`) and `settings`. A k of 0 and
#'   `none` give a choice without factors, with the reason. Refuses `trx_ruv_design_unsupported` and
#'   `trx_ruv_design_invalid` (parameters that do not fit the design), `trx_ruv_controls_invalid`.
#' @noRd
trxChooseRuvDesign <- function(x, parameters, k_record = NULL, log_cpm = NULL) {
    ruv_design <- parameters$ruv_design
    invalid <- function(message) trxAbort(message, c("trx_ruv_design_invalid", "trx_error"))
    if (!checkmate::test_choice(ruv_design, TRX_RUV_DESIGNS)) {
        trxAbort(sprintf("RUV design '%s' is not one of %s", paste(ruv_design, collapse = ","), paste(TRX_RUV_DESIGNS, collapse = ", ")),
            c("trx_ruv_design_unsupported", "trx_error"))
    }
    # Audit F-07: the canonical form, before `k_record`'s default reads the parameters.
    parameters <- trxCanonicalRuvDesignParameters(parameters)
    none <- identical(ruv_design, "none")
    supplied <- parameters$ruv_control_ids
    estimated <- !none && is.null(k_record) && is.na(parameters$ruv_k) && identical(parameters$ruv_k_rule, "num_sv_be")
    if (is.null(k_record)) k_record <- list(k = parameters$ruv_k, k_rule = parameters$ruv_k_rule)
    k <- k_record$k
    k_rule <- k_record$k_rule
    if (none) {
        if (!(is.null(k) || identical(length(k), 1L) && (is.na(k) || k == 0)) || !(is.null(k_rule) || all(is.na(k_rule)))) {
            invalid("the design none takes no k and no k rule")
        }
        k <- 0L
        k_rule <- NA_character_
    } else if (!estimated && (!checkmate::test_count(k) || !checkmate::test_choice(k_rule, TRX_RUV_K_RULES) || (identical(k_rule, "user") && k < 1))) {
        invalid(sprintf("the design %s needs k (a whole number; 1 or more when the user chose it) and its rule (%s)", ruv_design,
            paste(TRX_RUV_K_RULES, collapse = ", ")))
    }
    if (identical(ruv_design, "supplied_controls") != checkmate::test_character(supplied, min.len = 1L)) {
        invalid("ruv_control_ids names the control genes of the design supplied_controls, and is given for no other design")
    }
    controls <- trxRuvDesignControls(x$counts, x$design, ruv_design, supplied)
    if (estimated) {
        if (!is.matrix(log_cpm)) invalid("the estimated count (ruv_k_rule num_sv_be with no ruv_k) needs the TMM logCPM")
        trxRuvCheckAlignment(x$counts, x$design, log_cpm)
        groups <- as.character(x$design$group)
        limit <- ruvFactorCountLimit(ncol(x$counts), length(unique(groups)), length(unique(groups)), length(controls$control_genes))
        k_record <- resolveRuvFactorCount("num_sv_be", controls$control_rule, limit,
            estimate_count = \() estimateRuvFactorCount(log_cpm, groups, seed = parameters$ruv_k_seed))
        k <- k_record$k
        k_rule <- k_record$k_rule
    }
    k <- as.integer(round(k))
    applied <- k > 0L
    reason <- if (none) {
        "the design none: no correction for unwanted variation"
    } else if (!applied) {
        k_record$reason %||% "the k rule found no unwanted factor (k = 0)"
    } else {
        NA_character_
    }
    n_controls <- length(controls$control_genes)
    list(chosen = list(k = k, percentage = 100 * n_controls / nrow(x$counts), n_controls = n_controls, control_genes = controls$control_genes),
        ruv_design = list(id = ruv_design, method = controls$method, control_rule = controls$control_rule, unmatched = controls$unmatched,
            list_sha256 = controls$list_sha256, k = k, k_rule = k_rule, applied = applied, reason = reason, test_form = if (applied) "factors_in_model" else "none"),
        k_record = if (!none) k_record, k_estimate = if (estimated) k_record$estimate, control_genes = controls$control_genes, control_source = controls$control_rule, mode = "declared",
        settings = parameters[trxRuvDesignParameterNames])
}

#' The record's RUV outcomes of a declared design (KIT-129 W4), by the names of KIT-129 section 1e where the kernel
#' knows the value, and the names today's record reads (`k`, the control count and percentage, `ruv_applied`).
#' @param choice `trxChooseRuvDesign`'s value.
#' @noRd
trxRuvDesignOutcomes <- function(choice) {
    d <- choice$ruv_design
    list(ruv_design = d$id, ruv_method = d$method, ruv_replicate_structure = if (identical(d$method, "RUVs")) "group" else NA_character_,
        ruv_control_rule = d$control_rule, ruv_control_count = choice$chosen$n_controls, ruv_control_list_sha256 = d$list_sha256, ruv_control_list_unmatched = d$unmatched,
        ruv_k = d$k, ruv_k_rule = d$k_rule, ruv_test_form = d$test_form, ruv_applied = d$applied, ruv_not_applied_reason = d$reason,
        k = d$k, negative_control_percentage = choice$chosen$percentage, negative_control_count = choice$chosen$n_controls,
        # KIT-129 integration item 3: the estimate (the count's k, or the count beside the user's k) and its seed.
        ruv_k_estimate = choice$k_record$k_estimate %||% NA_integer_,
        ruv_k_seed = choice$k_record$k_seed %||% (if (identical(d$k_rule, "num_sv_be")) choice$settings$ruv_k_seed) %||% NA_integer_)
}

#' RUVs (rmd:406-426): the state keeps the counts and gains W_1..W_k in its design (the notebook's `ruvs_pdata`); the
#' removeBatchEffect display matrix is a saved table of the step's product. KIT-129 W4: a `grid` product that carries a
#' declared design (`trxChooseRuvDesign`) is run by `trxRuvFactors` (the factors on upper-quartile counts, by RUVs or
#' RUVg; no factor when k is 0); today's grid product takes today's path, unchanged.
#' @noRd
trxRuvsOperation <- function(binding) {
    newOmicsOperation(binding, "trx_ruvs", "whole_matrix", parameters = c("logcpm", "grid"), product_type = "trx.ruvs_result",
        kernel = function(object, parameters, store) {
            x <- trxAssayInputs(object)
            grid <- trxProductValue(store$manager, parameters$grid)
            log_cpm <- trxSavedMatrix(store$manager, trxProductValue(store$manager, parameters$logcpm)$log_cpm)
            declared <- grid[["ruv_design", exact = TRUE]]
            s <- if (is.null(declared)) {
                trxRuvs(x$counts, log_cpm, x$design, grid$chosen$control_genes, grid$chosen$k)
            } else {
                # Audit F-01: a source that already carries factors (a second RUV run on a head) gives them up, so the
                # state's design, and the formula the DA realises from it, hold this design's factors and no other.
                object@design_matrix <- object@design_matrix[, !startsWith(names(object@design_matrix), "W_"), drop = FALSE]
                c(trxRuvFactors(x$counts, log_cpm, x$design, grid$chosen$control_genes, grid$chosen$k, declared$method),
                    list(control_source = grid$control_source))
            }
            for (w in names(s$W)) object@design_matrix[[w]] <- unname(s$W[[w]])
            final <- store$checkpoint(object, "trx-ruvs")
            # The products it used are kept as references: an entity id is run-specific, a reference is compared by content.
            list(final = final, corrected = store$table(trxMatrixTable(s$log_cpm_ruvs)), W = s$W, k = grid$chosen$k,
                control_source = s$control_source, logcpm = trxProductReference(store$manager, parameters$logcpm),
                grid = trxProductReference(store$manager, parameters$grid),
                outcomes = if (is.null(declared)) trxRuvOutcomes(grid) else c(trxRuvDesignOutcomes(grid),
                    # KIT-129 fourth pass item 4: the factor-tied-to-groups check on this run's W (statistic, threshold, outcome).
                    ruvGroupAssociationCheck(if (length(s$W)) s$W, as.character(x$design$group), "transcriptomics",
                        declared$id, declared$k_rule, grid$chosen$k, grid$k_record$k_estimate %||% NA_integer_)[ruvGroupAssociationRecordFields]))
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
    # KIT-129 W4: a state without factors (the design none, or a k of 0) is tested on the groups alone.
    if (!length(w_terms)) return("~ 0 + group")
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
        # KIT-129 sixth pass item 4: a declared design's fields, one set in every omic (today's grid records none).
        o <- value$outcomes
        if (!is.null(o$ruv_design) && isSavedWorkflowReference(value$grid)) {
            choice <- tryCatch(hydrateSavedWorkflowEntity(context$store, value$grid)$value, error = function(error) NULL)
            listed <- if (identical(o$ruv_control_rule, "supplied_list")) length(unique(choice$settings$ruv_control_ids %||% character()))
            outcomes[RUV_DESIGN_RECORD_FIELDS] <- ruvDesignRecordFields(ruv_design = o$ruv_design,
                ruv_replicate_structure = o$ruv_replicate_structure,
                ruv_grouping_variable = if (identical(o$ruv_replicate_structure, "group")) "group",
                ruv_k_rule = o$ruv_k_rule, ruv_k_typed = if (identical(o$ruv_k_rule, "user")) choice$settings$ruv_k, ruv_k = o$ruv_k,
                ruv_k_estimate = o$ruv_k_estimate, ruv_k_seed = o$ruv_k_seed, ruv_k_permutations = choice$k_estimate$permutations,
                ruv_control_rule = o$ruv_control_rule, ruv_control_list_sha256 = o$ruv_control_list_sha256, ruv_control_list_count = listed,
                ruv_control_list_found = if (length(listed)) listed - length(o$ruv_control_list_unmatched %||% character()),
                ruv_control_count = o$ruv_control_count,
                ruv_first_pass_rule = if (identical(o$ruv_design, "first_pass_nulls")) "edger_qlf_f_test_groups",
                ruv_da_method = if (isTRUE(o$ruv_applied)) "edger_qlf_factors_in_model" else "edger_qlf_none")
        }
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
            # KIT-129 sixth pass item 5: a declared design's grid replays from its recorded declaration (the supplied
            # list is a declared parameter of the transcriptomics design, so the record holds it).
            recorded <- previous[[length(previous)]]$parameters
            if ("ruv_design" %in% vapply(recorded, function(parameter) parameter$name, character(1))) {
                design <- workflowReplayParameters(recorded, trxRuvDesignParameterNames)
                return(runTrxReferenceRuvs(manager, NULL, config, description, state_name = step$state, design_parameters = design))
            }
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
    # KIT-131: the PNG device on Windows writes nothing, and says nothing, from 260 characters.
    checkWorkflowFigureFiles(file.path(directory, paste0(name, c(".png", ".pdf"))), "a figure was drawn but its file was not written")
    invisible(NULL)
}

#' A publication's completion value: the source's checkpoint, the display; a failure is the display's message, except a
#' figure drawn but not written (`workflow_figure_unwritten`, KIT-131), which stops the step.
#' @noRd
trxPublicationResult <- function(manager, request, produce) {
    source <- request$input_references$state
    display <- tryCatch(produce(source), error = function(error) {
        if (inherits(error, "workflow_figure_unwritten")) stop(error)
        list(saved = FALSE, message = substr(conditionMessage(error), 1L, 1024L))
    })
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source, display = display)
}

#' The grid's publication (rmd:283-305, 373-393, 429-432): the CCA panels, the workbook, the RUVs display matrix TSV.
#' @noRd
computeWorkflowTrxRuvPublicationAction <- function(manager, request, variant) {
    trxPublicationResult(manager, request, function(source) {
        ruvs <- trxStateRuvs(manager, source)
        grid <- trxProductValue(manager, ruvs$grid)
        # KIT-129 fourth pass item 5: a declared design has no grid, so no canonical-correlation panels or grid workbook;
        # its corrected display matrix is published as today's.
        if (!is.null(grid[["ruv_design", exact = TRUE]])) {
            vroom::vroom_write(as.data.frame(trxSavedMatrix(manager, ruvs$corrected)) |> tibble::rownames_to_column("GeneID"),
                file.path(omicsStagedPublicationDir(request, "gene_qc_dir"), "ruv_normalised_log_values.tsv"))
            return(list(saved = TRUE, ruv_design = grid$ruv_design$id))
        }
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
        checkWorkflowFigureFiles(file.path(directory, paste0("trx_", stage, c("_pca.png", "_rle.png"))), "a figure was drawn but its file was not written")
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
#' project's (`parameters$dataset`): the grid summary only for a run of the legacy grid (KIT-192).
#' @noRd
computeWorkflowTrxDaPublicationAction <- function(manager, request, variant) {
    trxPublicationResult(manager, request, function(source) {
        context <- manager$.context
        if (is.null(context$snapshot$da_results)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
        manifest <- resolveWorkflowStateValue(context, context$snapshot$da_results)
        object <- readWorkflowOmicsState(transcriptomicsDeclarationBinding, manager, manifest$source)
        da_dir <- omicsStagedPublicationDir(request, "da_output_dir")
        graphs <- omicsStagedPublicationDir(request, "publication_graphs_dir")
        # KIT-220 (VAL-191 (1)): the figures, the top genes and the summary line at the run's threshold.
        qvalue <- manifest$parameters$qvalue
        annotated <- list()
        for (contrast in names(manifest$results)) {
            entry <- manifest$results[[contrast]]
            if (is.null(entry$products$da_contrast)) next
            results <- hydrateSavedWorkflowEntity(context$store, entry$products$da_contrast)$value$assays[[1L]]$results
            name <- entry$comparison
            de <- trxAnnotateDe(results$top_table, object@count_data[[1L]], name, object@gene_id_column, object@annotation_id_column,
                fdr_thresh = qvalue)
            utils::write.csv(de$de_annotated, file.path(da_dir, paste0("de_results_", name, ".csv")), row.names = FALSE)
            trxSaveFigure(trxVolcano(de$de_annotated, name, fdr_thresh = qvalue)$plot, file.path(graphs, "Volcano_Plots"),
                paste0("volcano_plot_", name), 8, 7)
            trxSaveFigure(trxSignificanceBar(de$de_annotated, name, fdr_thresh = qvalue)$plot, file.path(graphs, "NumSigDeMolecules"),
                paste0("summary_barplot_", name), 7, 5)
            annotated[[contrast]] <- de
        }
        if (!length(annotated)) cli::cli_abort("no contrast has a result", class = "workflow_execution_invalid")
        vroom::vroom_write(dplyr::bind_rows(lapply(annotated, `[[`, "top_de")), file.path(da_dir, "top_de_genes.tsv"))
        grid <- trxProductValue(manager, trxStateRuvs(manager, manifest$source)$grid)
        # KIT-192: a declared RUV design (none, a k of 0, or factors) ran no grid search, so there is no grid summary to
        # write (as the RUV publication, above); its choices are in the analysis record and the report.
        if (is.null(grid[["ruv_design", exact = TRUE]])) {
            trxWriteSummaryText(trxSummaryText(grid$best, length(grid$prelim_order), manifest$parameters$formula, names(annotated)[[1L]],
                annotated[[1L]]$de_annotated, grid$grid_summary, fdr_thresh = qvalue), request$parameters$dataset,
                file.path(da_dir, "ruv_optimization_summary.txt"))
        }
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
#' a `Run` column, as the notebook writes it for MOFA (rmd:714-724). KIT-136: when the head's RUV step applied no RUV
#' (design `none`, or k of 0; KIT-129), the TMM logCPM it read, in the same layout.
#' @noRd
trxWriteIntegrationExport <- function(manager, parameters, path) {
    state <- trxIntegrationRuvState(manager)
    if (!is.null(state) && !state$applied) {
        log_cpm <- trxSavedMatrix(manager, trxProductValue(manager, state$ruvs$logcpm)$log_cpm)
        vroom::vroom_write(as.data.frame(t(log_cpm)) |> tibble::rownames_to_column("Run"), path)
        return(invisible(path))
    }
    export <- manager$.context$snapshot$exports$trx_ruviii_export
    if (is.null(export)) cli::cli_abort("no accepted RUVIII export", class = "workflow_execution_invalid")
    # KIT-207 (SD-4 (a)): the export must have been made from the lineage's RUVs state, as `trxIntegrationInputSource`
    # requires; an earlier run's matrix is refused before the file is written.
    if (!is.null(state) && !identical(export$source$entity_id, state$reference$entity_id)) {
        abortIntegrationInput("the accepted RUVIII export was made from an earlier RUV step: run the export again, then save again")
    }
    newY <- hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value$newY
    vroom::vroom_write(as.data.frame(newY) |> tibble::rownames_to_column("Run"), path)
    invisible(path)
}

#' The head lineage's RUVs state (KIT-136): its reference, record and product, and whether it applied RUV (the
#' product's `ruv_applied` with k above 0; today's grid always does, a declared design `none` or a k of 0 does not).
#' @param manager The session's owner, inside its state scope.
#' @return `list(reference, record, ruvs, applied)`, or NULL when the head's lineage has no RUVs state.
#' @noRd
trxIntegrationRuvState <- function(manager) {
    lineage <- protIntegrationInputLineage(manager, until = "trx_ruvs")
    if (!"trx_ruvs" %in% names(lineage) || is.null(lineage$trx_ruvs)) return(NULL)
    reference <- getWorkflowStateReference(manager, "trx_ruvs")
    ruvs <- trxStateRuvs(manager, reference)
    list(reference = reference, record = lineage$trx_ruvs, ruvs = ruvs,
        applied = isTRUE(ruvs$outcomes$ruv_applied) && isTRUE(as.integer(ruvs$k) > 0L))
}

#' Write the transcriptomics integration input (KIT-136; KIT-061's contract): one view, `transcriptomics`, from the head
#' lineage's RUVs state (`trxIntegrationInputSource`); with no such state, or a manager that is not an owner, none is
#' written and an earlier one is removed (K061-3).
#' @param binding The binding (the shared writers' signature).
#' @param dir The omic's folder, `<project>/integration/transcriptomics/`.
#' @return TRUE when an input was written, else FALSE, invisibly.
#' @noRd
trxWriteIntegrationInput <- function(binding, manager, dir) {
    source <- if (inherits(manager, "DIAStateAdapter")) withWorkflowStateScope(manager, function() trxIntegrationInputSource(manager))
    if (is.null(source)) {
        removeIntegrationInput(dir)
        return(invisible(FALSE))
    }
    writeIntegrationInput(dir, source$manifest, source$samples, list(transcriptomics = source$matrix),
        list(transcriptomics = source$features),
        corrected = if (!is.null(source$corrected)) list(transcriptomics = source$corrected) else list())
    invisible(TRUE)
}

#' The transcriptomics integration input (KIT-136; KIT-060 audit rows C3a, C3b). Where the head's RUVs state applied RUV
#' the matrix is the accepted RUVIII export (INT-1, audit row C3b: RUV-III on the TMM logCPM with the RUV step's controls
#' and k, `newY`), genes by runs, role `corrected`; the export must have been made from that state. Where it applied
#' none (design `none` or k of 0) the matrix is the TMM logCPM the step read (no file held it before), role
#' `normalised`. Never the removeBatchEffect display matrix (C3a). The samples are the state's design without the RUVs
#' factors `W_1..W_k` (they are the DA model's, not the export's); the features are the gene id and name.
#' @param manager The session's owner, inside its state scope.
#' @return `list(manifest, samples, matrix, features)`, or NULL when the head's lineage has no RUVs state.
#' @noRd
trxIntegrationInputSource <- function(manager) {
    state <- trxIntegrationRuvState(manager)
    if (is.null(state)) return(NULL)
    reference <- state$reference
    spec <- readWorkflowStateMetadata(manager, reference)$data
    slots <- spec$slots
    design <- as.data.frame(materializeWorkflowTable(manager$.context, spec$components$design_matrix))
    design <- design[, !startsWith(names(design), "W_"), drop = FALSE]
    export <- NULL
    if (state$applied) {
        export <- manager$.context$snapshot$exports$trx_ruviii_export
        if (is.null(export)) {
            abortIntegrationInput("the corrected transcriptomics matrix is the RUVIII export: run the export, then save again")
        }
        if (!identical(export$source$entity_id, reference$entity_id)) {
            abortIntegrationInput("the accepted RUVIII export was made from an earlier RUV step: run the export again, then save again")
        }
        export$value <- hydrateSavedWorkflowEntity(manager$.context$store, export$reference)$value
        values <- t(export$value$newY)
        product <- list(type = "trx.ruviii", entity_id = export$reference$entity_id, sha256 = export$reference$payload$sha256)
    } else {
        values <- trxSavedMatrix(manager, trxProductValue(manager, state$ruvs$logcpm)$log_cpm)
        product <- list(type = "trx.tmm_logcpm", entity_id = state$ruvs$logcpm$entity_id, sha256 = state$ruvs$logcpm$payload$sha256)
    }
    sample_ids <- as.character(design[[slots$sample_id]])
    runs <- intersect(sample_ids, colnames(values))
    samples <- design[match(runs, sample_ids), , drop = FALSE]
    rownames(samples) <- NULL
    ids <- rownames(values)
    matrix <- data.frame(feature_id = ids, values[, runs, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE,
        row.names = NULL)
    values <- NULL
    # KIT-153 (INT-23): with RUV applied the matrix above is the RUVIII export (corrected); the data set's own matrix is
    # the TMM logCPM the export's RUV-III read (`trx_ruviii_export`'s `logcpm`, the RUVs state's), on its genes and runs.
    corrected <- NULL
    if (state$applied) {
        log_cpm <- trxSavedMatrix(manager, trxProductValue(manager, state$ruvs$logcpm)$log_cpm)
        uncorrected <- data.frame(feature_id = rownames(log_cpm), log_cpm, check.names = FALSE, stringsAsFactors = FALSE,
            row.names = NULL)
        log_cpm <- NULL
        corrected <- matrix
        matrix <- integrationInputUncorrectedOn(uncorrected, corrected, "transcriptomics")
        uncorrected <- NULL
    }
    gene_column <- slots$gene_id_column
    name_column <- slots$annotation_id_column
    genes <- materializeWorkflowTable(manager$.context, spec$components$count_data[[1L]], selected_columns = c(gene_column, name_column))
    gene_names <- as.character(genes[[name_column]])[match(ids, as.character(genes[[gene_column]]))]
    features <- data.frame(feature_id = ids, display_name = ifelse(is.na(gene_names) | !nzchar(gene_names), ids, gene_names),
        stringsAsFactors = FALSE)
    if (!name_column %in% names(features)) features[[name_column]] <- gene_names
    # KIT-167: with the gate on and a gene table saved at count import, the gene-level mapping columns and the gene
    # length (counts unchanged); without one the features are as before.
    if (referenceMappingEnabled()) features <- trxReferenceMappingFeatures(features, readReferenceGenes(manager))
    technical <- slots$technical_replicate_id
    if (!is.character(technical) || length(technical) != 1L || !isTRUE(technical %in% names(samples))) technical <- NA_character_
    manifest <- list(omic_type = "transcriptomics",
        workflow = manager$workflow_type %||% manager$.context$snapshot$workflow_type %||% NA_character_,
        matrix_role = "normalised", scale = "log2_cpm",
        normalisation = list(method = "TMM, log2 counts per million", source = list(name = "trx_tmm_logcpm",
            entity_id = state$ruvs$logcpm$entity_id, sha256 = state$ruvs$logcpm$payload$sha256)),
        correction = if (state$applied) trxIntegrationInputCorrection(manager, state, export, samples, technical) else NULL,
        source_state = list(name = "trx_ruvs", entity_id = reference$entity_id, sha256 = reference$payload$sha256,
            gene_id_column = gene_column, annotation_id_column = name_column, product = product),
        sample_id_column = slots$sample_id, group_column = slots$group_id, technical_replicate_column = technical,
        biological_sample_column = integrationInputBiologicalSampleColumn(samples, technical, slots$sample_id))
    list(manifest = manifest, samples = samples, matrix = matrix, features = features, corrected = corrected)
}

#' The correction behind a corrected transcriptomics input (KIT-136; audit row C3b): the RUVIII export's replicate
#' structure M, k and controls, and the RUV step that chose them (a new design's every field of
#' `RUV_DESIGN_RECORD_FIELDS`, as the analysis record reads them; today's grid names none).
#' @param state `trxIntegrationRuvState`'s value.
#' @param export The accepted export with its hydrated `value`.
#' @param samples The input's samples (for the replicate sets' levels).
#' @param technical The design's technical-replicate column, or NA.
#' @return The manifest's `correction`.
#' @noRd
trxIntegrationInputCorrection <- function(manager, state, export, samples, technical) {
    value <- export$value
    by_technical <- identical(value$m_source, "technical_replicate")
    grouping <- if (by_technical) technical else "group"
    correction <- list(id = "C3b", method = "RUVIII", m_source = value$m_source, k = value$k,
        control_count = value$control_count %||% sum(value$ctl_mask), control_source = value$control_source,
        replicate_structure = if (by_technical) "technical_replicates" else "group", grouping_variable = grouping,
        matrix = "ruviii_export_values", da_tests = "counts_with_ruvs_factors_in_model")
    outcomes <- workflowRecordStepOutcomes(manager$.context, "trx_ruvs", "trx_ruvs", state$record, list())
    if (all(RUV_DESIGN_RECORD_FIELDS %in% names(outcomes))) correction <- c(correction, outcomes[RUV_DESIGN_RECORD_FIELDS])
    correction <- lapply(correction, function(value) if (is.null(value)) NA else value)
    correction$grouping_levels <- I(if (is.character(grouping) && length(grouping) == 1L && isTRUE(grouping %in% names(samples)))
        unique(as.character(samples[[grouping]])) else character())
    c(list(state = "trx_ruvs"), correction)
}

#' The interactive volcano of one contrast (rmd:673-693; the binding's `volcano_widget`): Glimma on the contrast's DGELRT
#' with the annotation aligned to its order. A failure is the page action's recorded message.
#' KIT-220 (VAL-191 (1)): the status is Glimma's documented default (`edgeR::decideTests.DGEExact`, BH) at the given
#' threshold, not at its fixed 0.05.
#' @noRd
trxVolcanoWidget <- function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh = 0.05) {
    entry <- da_results_list$contrasts_results[[selected_assay]]$results[[selected_contrast]]
    if (is.null(entry)) return(NULL)
    name <- trxParseContrastEntry(selected_contrast)$name
    object <- da_results_list$theObject
    de <- trxAnnotateDe(entry$top_table, object@count_data[[1L]], name, object@gene_id_column, object@annotation_id_column,
        fdr_thresh = da_q_val_thresh)$de_annotated
    inputs <- trxGlimmaInputs(de, entry$qlf, name)
    Glimma::glimmaVolcano(entry$qlf, main = inputs$main, anno = inputs$anno, display.columns = inputs$display_columns,
        status = edgeR::decideTests.DGEExact(entry$qlf, adjust.method = "BH", p.value = da_q_val_thresh))
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
    # KIT-100: the filter and TMM/logCPM default to `partitioned` (user, 2026-11-07, on KIT-099's large pairs: exact, peak
    # .941 and .979 of `whole`); every other step to `whole`.
    default_variants = c(trx_expression_filter = "partitioned", trx_tmm_logcpm = "partitioned"),
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
    da_fit = function(data, design, formula_string, sample_id, ...) trxDaFit(data, design, formula_string),
    da_apply = function(fit, contrast, ...) trxDaApply(fit, contrast),
    da_format = trxDaFormat,
    da_count = trxDaCount,
    da_model = trxDaModel,
    da_replay_check = trxDaReplayCheck,
    qc_dir_role = "gene_qc_dir",
    filtering_metrics = trxFilteringMetrics,
    final_pca = trxFinalPca,
    volcano_widget = trxVolcanoWidget,
    integration_export = trxWriteIntegrationExport,
    # KIT-136: the declared integration input, written beside the integration export (`trxWriteIntegrationInput`).
    integration_input = trxWriteIntegrationInput
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
#' @param variant `whole` or `partitioned`; NULL runs the registry default (`partitioned`, KIT-100).
#' @noRd
runTrxReferenceExpressionFilter <- function(manager, config, description = NULL, state_name = NULL, variant = NULL) {
    step <- trxExpressionFilterStep(description, state_name)
    if (!is.null(variant)) step$parameters$variant <- variant
    runWorkflowOmicsReferenceStep(manager, step, config)
}

#' A product-only step on a source, rooted in the head (a hit from a store whose head does not root it).
#' @noRd
trxRunProduct <- function(manager, operation, source, parameters, dependencies = list()) {
    entry <- workflowOperationEntry(operation)
    # KIT-100: a product whose default is not `whole` runs that variant in process too, as in the worker (KIT-089's rule in
    # runWorkflowOmicsReferenceStep), so the provenance names what ran.
    if (is.null(parameters$variant) && !identical(entry$default_variant %||% "whole", "whole")) parameters$variant <- entry$default_variant
    # KIT-129 integration item 1: the declared design's choice when the request names it, as on the worker.
    chosen <- workflowRequestOperation(manager, entry, source, parameters)
    run <- runSavedOperation(manager, chosen$declaration, source, chosen$parameters, dependencies)
    if (is.null(manager$.context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(manager, run$reference, c(list(source), dependencies))
    run
}

#' The RUV step's request from the tab's settings (KIT-129 integration item 1), in process and on the worker: today's
#' grid parameters, or, with the gate on and a new design named (`ruv_design`), that design's canonical parameters.
#' @param settings The grid's parameters (`trxNormGridParameters`), with the design's when the tab names one.
#' @return `list(grid, design)`: the grid's parameters, or NULL; the design's, or NULL.
#' @noRd
trxNormRuvRequest <- function(settings) {
    design <- settings[["ruv_design", exact = TRUE]]
    if (!omicGuiRuvDesignEnabled() || !is.character(design) || identical(design, "legacy")) return(list(grid = settings, design = NULL))
    list(grid = NULL, design = trxCanonicalRuvDesignParameters(settings)[trxRuvDesignParameterNames])
}

#' The choice of a declared RUV design on a source (KIT-129 W4), rooted in the head as `trxRunProduct` roots the grid.
#' The declaration is not the registry's (`trxRuvGridOperation(binding, "declared")`): nothing selects it by default.
#' @noRd
trxRunRuvDesign <- function(manager, source, logcpm, design_parameters) {
    # Audit F-07: the identity is taken from the canonical parameters.
    design_parameters <- trxCanonicalRuvDesignParameters(design_parameters)
    run <- runSavedOperation(manager, trxRuvGridOperation(transcriptomicsDeclarationBinding, "declared"), source,
        c(list(logcpm = logcpm$entity_id), design_parameters[trxRuvDesignParameterNames]), list(logcpm))
    if (is.null(manager$.context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(manager, run$reference, list(source, logcpm))
    run
}

#' The RUVs state's description: today's words for today's grid; the design and what it did for a declared one.
#' @noRd
trxRuvsDescription <- function(choice) {
    declared <- choice[["ruv_design", exact = TRUE]]
    if (is.null(declared)) return(sprintf("RUVs: W_1..W_%d added to the design", choice$chosen$k))
    if (!declared$applied) return(sprintf("RUV design %s: no factor added to the design (%s)", declared$id, declared$reason))
    sprintf("RUV design %s (%s): W_1..W_%d added to the design", declared$id, declared$method, declared$k)
}

#' TMM/logCPM, the grid and RUVs on the head, the RUVs state committed (the products before it are its `preceded_by`).
#' @param grid_parameters trxRuvGridParameters' list.
#' @param design_parameters KIT-129 W4. NULL (the default): today's grid with `grid_parameters`. `trxRuvDesignParameters`'
#'   list: the declared design's choice in the grid's place; `grid_parameters` is not read.
#' @return The RUVs step's run, reference and value, with `logcpm` and `grid`, the products' runs.
#' @noRd
runTrxReferenceRuvs <- function(manager, grid_parameters, config, description = NULL, state_name = NULL, design_parameters = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the transcriptomics RUVs needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        logcpm <- trxRunProduct(manager, "trx_tmm_logcpm", source, list())
        grid <- if (is.null(design_parameters)) {
            trxRunProduct(manager, "trx_ruv_grid", source, c(list(logcpm = logcpm$reference$entity_id), grid_parameters[trxRuvGridParameterNames]),
                list(logcpm$reference))
        } else {
            trxRunRuvDesign(manager, source, logcpm$reference, design_parameters)
        }
        step <- list(operation = "trx_ruvs", parameters = list(logcpm = logcpm$reference$entity_id, grid = grid$reference$entity_id),
            state_name = state_name %||% "trx_ruvs", description = description %||% trxRuvsDescription(grid$value))
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
