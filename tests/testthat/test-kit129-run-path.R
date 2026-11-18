# KIT-129 integration: the new-design RUV kernels on the application's run path, behind the gate
# (`MultiScholaR.ruv_design_panel`, off by default). Item 1: each omic's new-design RUV step and its DA step are requested
# through the compute halves the tab uses in process and the managed worker uses (`workflowRequestOperation`), selected
# by the unified parameters; saved products name the design and the DA method; today's declarations are untouched.

kit129Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

# Two assays, two groups of six, one sample-wise factor uncorrelated with the groups, values under a detection limit
# missing (KIT-127's fixture shape, smaller).
kit129OmicsObject <- function(omic = "lipidomics", seed = 129L, p = 50L, technical = FALSE) {
    set.seed(seed)
    design <- data.frame(Sample_ID = sprintf("S%02d", 1:12), group = rep(c("A", "B"), each = 6L), stringsAsFactors = FALSE)
    # (third pass item 3) three biological units per group, each injected twice: the technical-replicate declaration
    if (technical) design$extract <- paste(design$group, rep(rep(1:3, each = 2L), 2L), sep = "_")
    f <- stats::residuals(stats::lm(stats::rnorm(12L) ~ factor(design$group)))
    f <- f / sqrt(mean(f^2))
    tables <- lapply(c(pos = "pos", neg = "neg"), function(assay) {
        effect <- c(rep(1, 5L), rep(-1, 5L), rep(0, p - 10L))
        y <- stats::rnorm(p, 20, 2) + outer(effect, as.numeric(design$group == "B")) + outer(stats::rnorm(p, sd = 0.6), f) +
            matrix(stats::rnorm(p * 12L, sd = 0.5), p)
        seen <- y >= stats::quantile(y, 0.08)
        # every feature observed at least twice in each group (the intensity filter's rule; KIT-127's fixtures)
        kept <- rowSums(seen[, design$group == "A"]) >= 2L & rowSums(seen[, design$group == "B"]) >= 2L
        seen[!kept, ] <- TRUE
        y[!seen] <- NA_real_
        dimnames(y) <- list(sprintf("%s_f%03d", assay, seq_len(p)), design$Sample_ID)
        table <- tibble::as_tibble(as.data.frame(y, check.names = FALSE))
        table <- tibble::add_column(table, database_identifier = rownames(y), .before = 1L)
        tibble::add_column(table, name = rownames(y), .after = 1L)
    })
    object <- if (identical(omic, "lipidomics")) {
        methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = "database_identifier",
            annotation_id_column = "name", design_matrix = design, sample_id = "Sample_ID", group_id = "group", args = list())
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = tables, metabolite_id_column = "database_identifier",
            annotation_id_column = "name", design_matrix = design, sample_id = "Sample_ID", group_id = "group", args = list())
    }
    if (technical) object@technical_replicate_id <- "extract"
    object
}

kit129OmicsOwner <- function(omic, env = parent.frame(), technical = FALSE) {
    root <- withr::local_tempdir(.local_envir = env)
    object <- kit129OmicsObject(omic, technical = technical)
    fixture <- if (identical(omic, "lipidomics")) lipidDeclarationsOwner(root, object = object) else metabDeclarationsOwner(root, object = object)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    fixture$object <- object
    fixture$contrasts <- data.frame(contrasts = "groupB-groupA", friendly_names = "B vs A", full_format = "B_vs_A=groupB-groupA",
        stringsAsFactors = FALSE)
    fixture
}

kit129OmicsRuvParameters <- function(ruv_design = "groups_all_features", ruv_k = 1L, ruv_k_rule = "user", grouping = "group") {
    list(ruv_design = ruv_design, ruv_grouping_variable = grouping, ruv_k = ruv_k, ruv_k_rule = ruv_k_rule, ruv_k_seed = NULL,
        ruv_control_ids = NULL)
}

kit129ProtStore <- function(env = parent.frame(), object = makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = env), "kit-129-run-path", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    manager
}

kit129ProtParameters <- function(ruv_mode = "manual", ruv_k = 1) {
    list(norm_method = "quantile", ruv_mode = ruv_mode, auto_percentage_min = NULL, auto_percentage_max = NULL, separation_metric = NULL,
        k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = ruv_k, grouping = "group")
}

test_that("(item 1) a request runs the new design's declaration only when it names the design, and today's otherwise", {
    entry <- function(name) workflowOperationEntry(name)
    pick <- function(name, parameters) workflowRequestOperation(NULL, entry(name), NULL, parameters)$declaration
    # today's declarations are the registry's, unchanged
    expect_identical(entry("normalization-ruv")$parameters, diaNormalizationOperation()$parameters)
    expect_identical(entry("lipid_ruv")$parameters, omicsRuvParameterNames)
    expect_identical(entry("lipid_da_contrast")$parameters, c("formula", "qvalue", "lfc", "contrast", "friendly_name"))
    # classic: the selector ruv_design
    classic <- kit129ProtParameters()
    expect_identical(pick("normalization-ruv", classic)$parameters, entry("normalization-ruv")$parameters)
    expect_identical(pick("normalization-ruv", c(classic, list(ruv_design = "groups_all_features", ruv_k_rule = "user",
        ruv_k_seed = NULL)))$parameters, diaNormalizationOperation("declared")$parameters)
    # DIA with limpa: the selector ruv_control_rule
    dia <- resolveDIARuvCovariatesParameters(list(), "group")
    expect_identical(pick("protein_ruv_covariates", dia)$parameters, entry("protein_ruv_covariates")$parameters)
    expect_identical(pick("protein_ruv_covariates", resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features"))$parameters,
        diaRuvCovariatesOperation("declared")$parameters)
    # transcriptomics: the declared design alone; a request with today's grid parameters stays today's (and is refused there)
    design <- c(list(logcpm = "x"), trxRuvDesignParameters("groups_all_features", 2L, "user"))
    expect_identical(pick("trx_ruv_grid", design)$parameters, c("logcpm", trxRuvDesignParameterNames))
    expect_identical(pick("trx_ruv_grid", c(design, trxRuvGridParameters()))$parameters, entry("trx_ruv_grid")$parameters)
    # lipidomics and metabolomics: the declared design alone
    for (name in c("lipid_ruv", "metab_ruv")) {
        expect_identical(pick(name, kit129OmicsRuvParameters())$parameters, omicsRuvDesignParameterNames, info = name)
        expect_identical(pick(name, c(kit129OmicsRuvParameters(), lipidDeclarationsRuvParameters()))$parameters, omicsRuvParameterNames,
            info = name)
    }
})

test_that("(item 1) today's declarations refuse the new design's parameters, the RUV step's and the DA method", {
    for (name in c("lipid_ruv", "metab_ruv")) {
        refuse <- workflowLegacyDesignRefusal(workflowOperationEntry(name)$declaration)
        expect_identical(refuse, assertOmicsLegacyRuvParameters, info = name)
        expect_error(refuse(list(ruv_k_rule = "user")), class = "workflow_execution_invalid")
        expect_silent(refuse(lipidDeclarationsRuvParameters()))
    }
    expect_null(workflowLegacyDesignRefusal(omicsRuvOperation(lipidomicsDeclarationBinding, "declared")))
    for (declaration in list(diaContrastOperation(), lipidDaContrastOperation(), metabDaContrastOperation())) {
        refuse <- workflowLegacyDesignRefusal(declaration)
        expect_error(refuse(list(ruv_da_method = "limma_ruv_factors_exact_se")), class = "workflow_execution_invalid")
    }
    expect_null(workflowLegacyDesignRefusal(diaContrastOperation("iq", "ruv_factors_exact_se_trend_robust")))
    expect_null(workflowLegacyDesignRefusal(omicsDaContrastOperation(lipidomicsDeclarationBinding, "declared")))
    # the declared DA declarations add the method to the identity's parameters and nothing else
    legacy <- diaContrastOperation("limpa")
    declared <- diaContrastOperation("limpa", "dpcde_ruv_factors_exact_se")
    expect_identical(declared$parameters, c(legacy$parameters, "ruv_da_method"))
    expect_identical(declared[setdiff(names(declared), c("kernel", "parameters"))], legacy[setdiff(names(legacy), c("kernel", "parameters"))])
})

test_that("(item 1) with the gate off the tab's requests are today's; with it on they carry the design the tab names", {
    settings <- c(kit129ProtParameters(), list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = 7L))
    norm_data <- new.env()
    norm_data$normalization_input_reference <- list(entity_id = "source")
    request <- function() diaNormalizationRequest(NULL, settings, norm_data, "group")$parameters
    off <- withr::with_options(list(MultiScholaR.ruv_design_panel = NULL), request())
    expect_identical(names(off), c(setdiff(names(kit129ProtParameters()), "grouping"), "grouping"))
    on <- withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), request())
    expect_identical(on[c("ruv_design", "ruv_k_rule", "ruv_k_seed")], list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = 7L))
    dia_settings <- list(ruv_mode = "covariates", ruv_control_rule = "all_features")
    expect_identical(withr::with_options(list(MultiScholaR.ruv_design_panel = NULL), diaRuvCovariatesRequestParameters(list(), "group", dia_settings)),
        resolveDIARuvCovariatesParameters(list(), "group"))
    expect_identical(withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), diaRuvCovariatesRequestParameters(list(), "group", dia_settings)),
        resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features"))
    omics <- c(lipidDeclarationsRuvParameters(), kit129OmicsRuvParameters())
    expect_identical(withr::with_options(list(MultiScholaR.ruv_design_panel = NULL), omicsRuvStep(lipidomicsDeclarationBinding, omics))$parameters,
        lipidDeclarationsRuvParameters())
    expect_identical(withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), omicsRuvStep(lipidomicsDeclarationBinding, omics))$parameters,
        kit129OmicsRuvParameters())
    grid <- c(trxRuvGridParameters(), list(ruv_design = "groups_all_features", ruv_k = 2L, ruv_k_rule = "user"))
    expect_identical(withr::with_options(list(MultiScholaR.ruv_design_panel = NULL), trxNormRuvRequest(grid)), list(grid = grid, design = NULL))
    expect_identical(withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), trxNormRuvRequest(grid))$design,
        trxRuvDesignParameters("groups_all_features", 2L, "user"))
})

test_that("(item 1) classic proteomics: the compute half runs the new design and its DA names the method in the identity", {
    manager <- kit129ProtStore()
    source <- getDIAStateReference(manager)
    declared <- c(kit129ProtParameters(), list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = NULL))
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        today <- computeDIAReferenceNormalization(manager, source, kit129ProtParameters())
        new <- computeDIAReferenceNormalization(manager, source, declared)
        spec <- readDIAStateMetadata(manager, new$value$final)$data
        step <- list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA, contrast = "a=grouptreated-groupcontrol",
            annotations = spec$components$protein_id_table$entity_id, name = "a", comparison = "a")
        list(today = today, new = new, method = workflowSourceRuvDaMethod(manager, new$value$final),
            today_method = workflowSourceRuvDaMethod(manager, today$value$final),
            da = computeDIAReferenceDaContrast(manager, new$value$final, step, "iq")$run)
    }))
    expect_identical(out$new$status, "computed")
    expect_identical(out$new$value$parameters$ruv_design, "groups_all_features")
    expect_false(identical(out$new$reference$entity_id, out$today$reference$entity_id))
    expect_identical(out$method, "ruv_factors_exact_se_trend_robust")
    expect_null(out$today_method)
    expect_identical(out$da$status, "computed")
    expect_identical(out$da$value$model$tested_matrix, "uncorrected")
    # the product's identity was taken over the method: the same request without it is another identity
    expect_true(grepl("^[a-f0-9]{64}$", out$da$reference$entity_id))
})

test_that("(item 1) lipidomics and metabolomics: the new design and its DA in process and on the worker", {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) {
            owner$job <- NULL
            list(status = status)
        }, .package = "MultiScholaR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    for (omic in c("lipidomics", "metabolomics")) {
        prefix <- if (identical(omic, "lipidomics")) "lipid" else "metab"
        fixture <- kit129OmicsOwner(omic)
        owner <- fixture$owner
        binding <- omicsDeclarationBindings()[[omic]]
        # the worker first, on the head: it computes and saves the new design's product without moving the head
        worker <- omicsWorkerRun(fixture, paste0(prefix, "_ruv"), kit129OmicsRuvParameters())
        expect_identical(worker$status, "complete", info = paste(omic, worker$error))
        expect_identical(worker$result$checkpoint$slots$args$ruvFactorsDeclaration$ruv_design, "groups_all_features", info = omic)
        # in process, the same request is the worker's product: one identity (the controller's handle indexes the worker's
        # products only when its action is accepted, so it computes it again here)
        ruv <- kit129Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit129OmicsRuvParameters(), list()))
        expect_identical(ruv$run$status, "computed", info = omic)
        expect_identical(ruv$run$reference$identity, worker$result$reference$identity, info = omic)
        expect_identical(ruv$value$ruv_design, "groups_all_features", info = omic)
        expect_identical(unname(unlist(ruv$value$best_k)), c(1L, 1L), info = omic)
        expect_identical(workflowSourceRuvDaMethod(owner, ruv$reference), "limma_ruv_factors_exact_se", info = omic)
        state <- withWorkflowStateScope(owner, function() readWorkflowExecutionComponents(owner, ruv$reference, fields = binding$fields))
        expect_identical(state@args$ruvFactorsDeclaration$ruv_design, "groups_all_features", info = omic)
        expect_identical(names(state@args$ruv_factors), c("pos", "neg"), info = omic)
        expect_identical(state@args$ruv_uncorrected_assays, methods::slot(fixture$object, binding$data_slot), info = omic)
        da <- kit129Quiet(runWorkflowOmicsReferenceDaAnalysis(binding, owner, ruv$reference,
            list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$contrasts))
        result <- da$manifest$results[[1L]]
        expect_identical(result$model$ruv_da_method, "limma_ruv_factors_exact_se", info = omic)
        expect_identical(names(result$assays), c("pos", "neg"), info = omic)
        # the worker's DA on the new design's state names the method too
        step <- omicsDaContrastSteps(binding, list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            fixture$contrasts)[[1L]]
        worker_da <- omicsWorkerRun(fixture, step$operation, step$parameters)
        expect_identical(worker_da$status, "complete", info = paste(omic, worker_da$error))
        expect_true(isSavedWorkflowReference(worker_da$result$reference), info = omic)
    }
})

# ---- Item 2: the DA of every new-design object is the new design's --------------------------------------------------

# The DA step on a classic state: one contrast through the compute half, in the reference store.
kit129ProtDa <- function(manager, state) {
    spec <- readDIAStateMetadata(manager, state)$data
    step <- list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA, contrast = "a=grouptreated-groupcontrol",
        annotations = spec$components$protein_id_table$entity_id, name = "a", comparison = "a")
    computed <- computeDIAReferenceDaContrast(manager, state, step, "iq")
    table <- as.data.frame(materializeDIATable(manager$.context, computed$run$value$tables$da_proteins_long))
    list(run = computed$run, base_fit = computed$base_fit, table = table[order(table$Protein.Group %||% table[[1L]]), , drop = FALSE])
}

# The new designs' test without W, by hand: limma on the normalised matrix, eBayes and treat() with the settings.
kit129NoFactorsLimma <- function(object) {
    x <- as.matrix(object@protein_quant_table[, object@design_matrix$Run])
    rownames(x) <- object@protein_quant_table[[object@protein_id_column]]
    design <- stats::model.matrix(~ 0 + group, object@design_matrix)
    technical <- resolveProtTechnicalReplicateBlock(object@design_matrix, object@group_id, object@technical_replicate_id)
    fit <- if (isTRUE(technical$technical)) {
        correlation <- limma::duplicateCorrelation(x, design, block = technical$block)$consensus.correlation
        limma::lmFit(x, design, block = technical$block, correlation = correlation)
    } else limma::lmFit(x, design)
    fit <- limma::contrasts.fit(fit, limma::makeContrasts(contrasts = "grouptreated-groupcontrol", levels = colnames(design)))
    tested <- limma::treat(limma::eBayes(fit, trend = TRUE, robust = TRUE), lfc = 0, trend = TRUE, robust = TRUE)
    stats::setNames(tested$p.value[, 1L], rownames(tested$coefficients))
}

test_that("(item 2) classic proteomics: the option set's None declares itself and is tested with the trend and robust settings (M30)", {
    manager <- kit129ProtStore()
    source <- getDIAStateReference(manager)
    none <- c(kit129ProtParameters("skip"), list(ruv_design = "none", ruv_k_rule = NULL, ruv_k_seed = NULL))
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        today <- computeDIAReferenceNormalization(manager, source, kit129ProtParameters("skip"))
        new <- computeDIAReferenceNormalization(manager, source, none)
        list(today = today, new = new, state = readDIAExecutionComponents(manager, new$value$final),
            normalised = readDIAExecutionComponents(manager, new$value$normalized),
            method = workflowSourceRuvDaMethod(manager, new$value$final),
            today_da = kit129ProtDa(manager, today$value$final), new_da = kit129ProtDa(manager, new$value$final),
            refused = tryCatch(computeDIAReferenceNormalization(manager, source, c(kit129ProtParameters(), list(ruv_design = "none",
                ruv_k_rule = "user", ruv_k_seed = NULL))), error = function(e) e))
    }))
    expect_s3_class(out$refused, "workflow_execution_invalid")
    expect_identical(out$new$value$parameters$ruv_design, "none")
    expect_false(identical(out$new$value$final, out$new$value$normalized))
    expect_identical(out$state@args$ruvCovariates[c("ruv_mode", "ruv_design", "k", "ruv_applied")],
        list(ruv_mode = "none", ruv_design = "none", k = 0L, ruv_applied = FALSE))
    expect_identical(out$state@protein_quant_table, out$normalised@protein_quant_table)
    expect_identical(out$method, "ruv_none_trend_robust")
    # the DA: the new design's test (trend and robust at TREAT), another product and another base fit than today's skip
    expect_identical(out$new_da$run$value$model[c("ruv_mode", "tested_matrix", "ruv_design", "exact_standard_errors")],
        list(ruv_mode = "none", tested_matrix = "normalised", ruv_design = "none", exact_standard_errors = FALSE))
    expect_false(identical(out$new_da$run$reference$entity_id, out$today_da$run$reference$entity_id))
    expect_false(identical(out$new_da$base_fit$entity_id, out$today_da$base_fit$entity_id))
    expected <- kit129NoFactorsLimma(out$normalised)
    got <- out$new_da$table
    ids <- got[[out$normalised@protein_id_column]]
    expect_equal(got$raw_pvalue, unname(expected[ids]), tolerance = 1e-12)
})

test_that("(item 2) classic proteomics: a count of 0 declares the design without factors and is tested as None", {
    manager <- kit129ProtStore()
    source <- getDIAStateReference(manager)
    testthat::local_mocked_bindings(computeProtRuvModelFactors = function(theObject, ...) list(status = "none", k = 0L,
        reason = "the count found no unwanted factor (k = 0)", ruv_k_rule = "num_sv_be", ruv_k_estimate = 0L, k_estimate = NULL,
        W = NULL, controls = NULL, ruv_controls = "all"), .package = "MultiScholaR")
    automatic <- kit129ProtParameters("automatic")
    automatic[c("auto_percentage_min", "auto_percentage_max", "separation_metric", "k_penalty_weight", "adaptive_k_penalty")] <-
        list(5, 10, "max_difference", 0.5, TRUE)
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, source, c(automatic, list(ruv_design = "groups_all_features",
            ruv_k_rule = "num_sv_be", ruv_k_seed = NULL)))
        list(run = run, state = readDIAExecutionComponents(manager, run$value$final),
            normalised = readDIAExecutionComponents(manager, run$value$normalized), da = kit129ProtDa(manager, run$value$final))
    }))
    expect_false(out$run$value$ruv_result$ruv_applied)
    expect_identical(out$state@args$ruvCovariates[c("ruv_mode", "ruv_design", "k", "ruv_k_rule", "ruv_k_estimate")],
        list(ruv_mode = "none", ruv_design = "groups_all_features", k = 0L, ruv_k_rule = "num_sv_be", ruv_k_estimate = 0L))
    expect_identical(out$da$run$value$model$ruv_mode, "none")
    expected <- kit129NoFactorsLimma(out$normalised)
    ids <- out$da$table[[out$normalised@protein_id_column]]
    expect_equal(out$da$table$raw_pvalue, unname(expected[ids]), tolerance = 1e-12)
})

test_that("(item 2) classic proteomics: the factors-in-the-model test is saved under identities that name it", {
    manager <- kit129ProtStore()
    source <- getDIAStateReference(manager)
    declared <- c(kit129ProtParameters(), list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = NULL))
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        new <- computeDIAReferenceNormalization(manager, source, declared)
        first <- kit129ProtDa(manager, new$value$final)
        state <- readDIAExecutionComponents(manager, new$value$final)
        list(first = first, state = state)
    }))
    # the base fit is in the store, its identity naming the test; the result table equals the kernel called directly
    expect_true(isSavedWorkflowReference(out$first$base_fit))
    direct <- kit129Quiet(differentialAbundanceAnalysis(out$state, contrasts_tbl = data.frame(contrasts = "a=grouptreated-groupcontrol"),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = out$state@protein_id_column, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
    table <- direct$contrasts_results$results[[1L]]
    ids <- out$first$table[[out$state@protein_id_column]]
    expect_equal(out$first$table$raw_pvalue, table[ids, "raw_pvalue"], tolerance = 1e-12)
    expect_equal(out$first$table$log2FC %||% out$first$table$logFC, table[ids, "logFC"], tolerance = 1e-12)
})

test_that("(item 2) DIA with limpa: the DA declaration of a declared state refuses a source whose declaration names another method", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    object@args$ruvCovariates <- list(k = 1L, ruv_mode = "covariates", ruv_k_rule = "num_sv_be")
    expect_identical(diaObjectRuvDaMethod(object), "dpcde_ruv_factors_exact_se")
    expect_error(diaContrastOperation("limpa", "dpcde_ruv_none")$kernel(object, list(ruv_da_method = "dpcde_ruv_none"), NULL),
        class = "workflow_execution_invalid")
    object@args$ruvCovariates <- list(k = 0L, ruv_mode = "none", ruv_k_rule = "num_sv_be")
    expect_identical(diaObjectRuvDaMethod(object), "dpcde_ruv_none")
})

test_that("(item 2) lipidomics: the option set's None is tested as today's test, under its own identity", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    fixture <- kit129OmicsOwner("lipidomics")
    owner <- fixture$owner
    binding <- lipidomicsDeclarationBinding
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    none <- kit129OmicsRuvParameters("none", NULL, NULL)
    ruv <- kit129Quiet(runWorkflowOmicsReferenceRuv(binding, owner, none, list()))
    expect_identical(ruv$value$ruv_design, "none")
    expect_identical(unname(unlist(ruv$value$best_k)), c(0L, 0L))
    da <- kit129Quiet(runWorkflowOmicsReferenceDaAnalysis(binding, owner, ruv$reference, settings, fixture$contrasts))
    expect_identical(da$manifest$results[[1L]]$model$ruv_da_method, "limma_ruv_factors_exact_se")
    today <- kit129Quiet(runLipidsDA(fixture$object, fixture$contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE))
    key <- function(t) t[order(t$assay, t$lipid_id), c("assay", "lipid_id", "logFC", "P.Value", "fdr_qvalue")]
    expect_equal(as.data.frame(key(da$results$da_lipids_long)), as.data.frame(key(today$da_lipids_long)), tolerance = 1e-12,
        ignore_attr = TRUE)
})

# ---- Item 3: k is resolved once on the run path, and the step records the rule, the estimate and the seed -----------

kit129TrxOwner <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    reference <- trxDeclReference("trx_edge_12", env = env)
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    kit129Quiet(runTrxReferenceExpressionFilter(fixture$owner, list()))
    commitWorkflowOmicsStoreStep(fixture$owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
    fixture
}

test_that("(item 3) transcriptomics: a k given under the estimated count is refused on the run path; the user's k records the count beside it", {
    skip_if_not_installed("sva")
    fixture <- kit129TrxOwner()
    owner <- fixture$owner
    out <- kit129Quiet(withWorkflowStateScope(owner, function() {
        source <- getWorkflowStateReference(owner)
        logcpm <- trxRunProduct(owner, "trx_tmm_logcpm", source, list())
        design <- function(...) c(list(logcpm = logcpm$reference$entity_id), trxRuvDesignParameters("groups_all_features", ...))
        list(given = tryCatch(trxRunProduct(owner, "trx_ruv_grid", source, design(2L, "num_sv_be"), list(logcpm$reference)),
                error = function(e) e),
            user = trxRunProduct(owner, "trx_ruv_grid", source, design(2L, "user"), list(logcpm$reference)),
            log_cpm = trxSavedMatrix(owner, logcpm$value$log_cpm), x = trxAssayInputs(readTrxWorkflowState(owner, source)))
    }))
    expect_s3_class(out$given, "trx_ruv_design_invalid")
    expected <- estimateRuvFactorCount(out$log_cpm, as.character(out$x$design$group))
    expect_identical(out$user$value$ruv_design[c("k", "k_rule")], list(k = 2L, k_rule = "user"))
    expect_identical(out$user$value$k_record[c("k_estimate", "k_seed")], list(k_estimate = expected$k, k_seed = RUV_FACTOR_COUNT_SEED))
    outcomes <- trxRuvDesignOutcomes(out$user$value)
    expect_identical(outcomes[c("ruv_k", "ruv_k_rule", "ruv_k_estimate", "ruv_k_seed")],
        list(ruv_k = 2L, ruv_k_rule = "user", ruv_k_estimate = expected$k, ruv_k_seed = RUV_FACTOR_COUNT_SEED))
})

test_that("(item 3) classic proteomics: the step's record names the design, the rule, the estimate and the seed", {
    skip_if_not_installed("sva")
    manager <- kit129ProtStore()
    source <- getDIAStateReference(manager)
    declared <- c(kit129ProtParameters(), list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = 11L))
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, source, declared)
        record <- list(config = list(ruv_k = run$value$best_k), data = list(products = list(normalization_science = run$reference)))
        list(run = run, normalised = readDIAExecutionComponents(manager, run$value$normalized),
            outcomes = workflowRecordStepOutcomes(manager$.context, "ruv_corrected", "normalization-ruv", record, list()),
            today = workflowRecordStepOutcomes(manager$.context, "ruv_corrected", "normalization-ruv",
                list(config = list(ruv_k = 1), data = list(products = list(normalization_science = computeDIAReferenceNormalization(manager,
                    source, kit129ProtParameters())$reference))), list()))
    }))
    n <- out$normalised
    technical <- resolveProtTechnicalReplicateBlock(n@design_matrix, n@group_id, n@technical_replicate_id)
    expected <- estimateRuvFactorCount(as.matrix(n@protein_quant_table[, n@design_matrix$Run]), n@design_matrix$group, seed = 11L,
        technical_replicates = isTRUE(technical$technical), units = if (isTRUE(technical$technical)) technical$block)
    expect_identical(out$run$value$ruv_result[c("ruv_k_rule", "ruv_k_estimate", "ruv_k_seed")],
        list(ruv_k_rule = "user", ruv_k_estimate = expected$k, ruv_k_seed = 11L))
    expect_identical(out$outcomes[c("ruv_design", "ruv_k_rule", "ruv_k_estimate", "ruv_k_seed", "ruv_applied")],
        list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_estimate = expected$k, ruv_k_seed = 11L, ruv_applied = TRUE))
    expect_null(out$today$ruv_design)
    expect_null(out$today$ruv_k_rule)
})

test_that("(item 3) lipidomics: the step's record names the design, the rule, the estimate per assay and the seed", {
    skip_if_not_installed("sva")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    fixture <- kit129OmicsOwner("lipidomics")
    owner <- fixture$owner
    binding <- lipidomicsDeclarationBinding
    counted <- kit129Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit129OmicsRuvParameters(ruv_k = NULL, ruv_k_rule = "num_sv_be"), list()))
    refused <- tryCatch(kit129Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit129OmicsRuvParameters(ruv_k = 1L, ruv_k_rule = "num_sv_be"),
        list())), error = function(e) e)
    expect_s3_class(refused, "workflow_execution_invalid")
    # the user's k on its own session (each run reads the normalised state, the head)
    other <- kit129OmicsOwner("lipidomics")
    owner <- other$owner
    user <- kit129Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit129OmicsRuvParameters(), list()))
    direct <- kit129Quiet(runOmicsRuvFactors(fixture$object, "groups_all_features", ruv_k_rule = "num_sv_be"))
    expect_identical(counted$value$ruv_k_estimate, lapply(direct$factors, `[[`, "ruv_k_estimate"))
    expect_identical(user$value$ruv_k_estimate, lapply(direct$factors, `[[`, "ruv_k_estimate"))
    expect_identical(user$value$ruv_k_seed, RUV_FACTOR_COUNT_SEED)
    record <- buildWorkflowAnalysisRecord(owner)
    step <- Filter(function(step) identical(step$operation, "lipid_ruv"), record$steps)
    expect_length(step, 1L)
    last <- step[[1L]]
    expect_identical(last$outcomes[c("ruv_design", "ruv_k_rule", "ruv_k_seed")],
        list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = RUV_FACTOR_COUNT_SEED))
    expect_identical(last$outcomes$ruv_k_estimate_per_assay, lapply(direct$factors, `[[`, "ruv_k_estimate"))
    expect_true(all(c("ruv_design", "ruv_k", "ruv_k_rule", "ruv_k_seed") %in% vapply(last$parameters, `[[`, "", "name")))
})

# ---- Item 5: end to end through the run path, in process and by the worker-facing entry, equal to the kernel -----------

# DIA with limpa's fixture: test-dia-limpa-ruv-rules.R's simulation and DPC quantification (copied; each test file runs
# on its own).
kit129LimpaSimulate <- function(seed, proteins = 200L, per_group = 6L, groups = 2L, factors = 1L, strength = 3, units = 0L) {
    set.seed(seed)
    sim <- limpa::simProteinDataSet(n.peptides = proteins * 4L, n.groups = groups, samples.per.group = per_group,
        peptides.per.protein = 4L, mu.range = c(2, 10), prop.de = 0, dpc.slope = 0.7, prop.missing = 0.25)
    y <- sim$other$E.complete
    samples <- ncol(y)
    ids <- unique(sim$genes$Protein)
    protein_of <- match(sim$genes$Protein, ids)
    group_index <- sim$targets$Group
    contrast <- c(-0.5, 0.5, 0)[group_index]
    n_effects <- round(0.2 * proteins)
    delta <- numeric(proteins)
    delta[sample.int(proteins, n_effects)] <- rep(c(1, -1), c(ceiling(0.5 * n_effects), n_effects - ceiling(0.5 * n_effects)))
    term <- outer(delta, contrast)
    if (factors > 0L) {
        scores <- matrix(stats::rnorm(samples * factors), samples)
        loadings <- sqrt(strength * 0.09) * matrix(stats::rnorm(proteins * factors), proteins)
        term <- term + loadings %*% t(scores)
    }
    within_group <- stats::ave(group_index, group_index, FUN = seq_along)
    replicates <- within_group
    if (units > 0L) {
        replicates <- ceiling(within_group / (per_group / units))
        unit <- (group_index - 1L) * units + replicates
        term <- term + matrix(stats::rnorm(proteins * max(unit), 0, 0.3), proteins)[, unit]
    }
    y <- y + term[protein_of, , drop = FALSE]
    probability <- stats::plogis(sim$dpc[[1L]] + y * sim$dpc[[2L]], lower.tail = FALSE)
    y[stats::rbinom(length(y), prob = probability, size = 1L) == 1L] <- NA
    order <- sample.int(nrow(y))
    runs <- sprintf("run_%03d", seq_len(samples))
    design <- data.frame(Run = runs, group = c("A", "B", "C")[group_index], replicates = replicates, stringsAsFactors = FALSE)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    list(y = y[order, , drop = FALSE], protein = sim$genes$Protein[order], precursor = rownames(y)[order], design = design)
}

kit129LimpaQuantify <- function(sim, engine) {
    y <- sim$y
    dimnames(y) <- list(paste(sim$protein, sim$precursor, sep = "%"), sim$design$Run)
    object <- methods::new("PeptideQuantitativeData")
    object@peptide_data <- data.frame(Protein.Ids = sim$protein, Stripped.Sequence = sim$precursor, stringsAsFactors = FALSE)
    object@peptide_matrix <- y
    object@protein_id_column <- "Protein.Ids"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@is_logged_data <- TRUE
    object@design_matrix <- sim$design
    object@sample_id <- "Run"
    object@group_id <- "group"
    object@technical_replicate_id <- "tech_rep_group"
    object@args <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))
    counts <- tibble::tibble(Protein.Ids = unique(sim$protein), peptide_count = 4L, peptidoform_count = 4L)
    quantified <- NULL
    utils::capture.output(quantified <- suppressWarnings(suppressMessages(proteinMissingValueImputationLimpa(object,
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = engine, peptide_summary = counts))))
    quantified
}


kit129LimpaManager <- function(quantified, env = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = env)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$saveState("correlation_filtered", quantified, list(), "fixture")
    manager$setWorkflowType("DIA_LIMPA")
    manager
}

# The contrast's step as the DA request half builds it, and its table read back from the product.
kit129DaStep <- function(manager, state, contrast) {
    spec <- readDIAStateMetadata(manager, state)$data
    list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA, contrast = contrast,
        annotations = spec$components$protein_id_table$entity_id, name = sub("=.*$", "", contrast), comparison = sub("=.*$", "", contrast))
}
kit129DaTable <- function(manager, value, id_column) {
    table <- as.data.frame(materializeDIATable(manager$.context, value$tables$contrasts_results_table))
    table <- table[order(table[[id_column]]), c(id_column, "logFC", "raw_pvalue", "fdr_qvalue")]
    rownames(table) <- NULL
    table
}
kit129KernelTable <- function(result, id_column) {
    table <- result$contrasts_results$results[[1L]]
    table <- data.frame(id = rownames(table), logFC = table$logFC, raw_pvalue = table$raw_pvalue, fdr_qvalue = table$fdr_qvalue,
        stringsAsFactors = FALSE)
    names(table)[[1L]] <- id_column
    table <- table[order(table[[id_column]]), ]
    rownames(table) <- NULL
    table
}

test_that("(item 5) classic proteomics: in process, by the worker's handlers and by the kernel, one result", {
    skip_if_not_installed("sva")
    declared <- c(kit129ProtParameters("automatic", NULL), list(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = 5L))
    declared[c("auto_percentage_min", "auto_percentage_max", "separation_metric", "k_penalty_weight", "adaptive_k_penalty")] <-
        list(5, 10, "max_difference", 0.5, TRUE)
    contrast <- "a=grouptreated-groupcontrol"
    run_in <- function(path) {
        manager <- kit129ProtStore()
        source <- getDIAStateReference(manager)
        kit129Quiet(withDIAStateScope(manager, function() {
            if (identical(path, "process")) {
                run <- computeDIAReferenceNormalization(manager, source, declared)
                final <- run$value$final
                da <- computeDIAReferenceDaContrast(manager, final, kit129DaStep(manager, final, contrast), "iq")$run$value
            } else {
                worker <- computeWorkflowDiaNormalizationAction(manager, list(input_references = list(state = source), parameters = declared), "whole")
                run <- list(value = hydrateSavedWorkflowEntity(manager$.context$store, worker$reference)$value)
                final <- run$value$final
                handoff <- computeWorkflowDiaContrastAction(manager, list(input_references = list(state = final),
                    parameters = c(kit129DaStep(manager, final, contrast), list(route = "iq"))), "whole")
                da <- hydrateSavedWorkflowEntity(manager$.context$store, handoff$reference)$value
            }
            list(k = run$value$best_k, state = readDIAExecutionComponents(manager, final),
                normalised = readDIAExecutionComponents(manager, run$value$normalized), table = kit129DaTable(manager, da, "Protein.Group"))
        }))
    }
    process <- run_in("process")
    worker <- run_in("worker")
    expect_identical(worker$k, process$k)
    expect_identical(worker$table, process$table)
    # the kernel called directly on the normalised object: the count, the factors, the DA entry point
    normalised <- process$normalised
    factors <- kit129Quiet(computeProtRuvModelFactors(normalised, ruv_k_rule = "num_sv_be", ruv_grouping_variable = "group", ruv_k_seed = 5L))
    expect_identical(factors$k, process$k)
    skip_if(factors$k == 0L, "the count found no factor on this fixture")
    corrected <- kit129Quiet(removeRowsWithMissingValuesPercent(theObject = ruvIII_C_Varying(protRuvGroupReplicateObject(normalised),
        ruv_grouping_variable = "group", ruv_number_k = factors$k, ctrl = factors$controls)))
    state <- attachProtRuvModelFactors(corrected, normalised, factors)
    kernel <- kit129Quiet(differentialAbundanceAnalysis(state, contrasts_tbl = data.frame(contrasts = contrast), formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Group",
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
    expected <- kit129KernelTable(kernel, "Protein.Group")
    expect_identical(process$table$Protein.Group, expected$Protein.Group)
    expect_equal(process$table[-1L], expected[-1L], tolerance = 1e-12)
})

test_that("(item 5) DIA with limpa: in process, by the worker's handlers and by the kernel, one result", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    skip_if_not_installed("fastlimpadpc")
    sim <- kit129LimpaSimulate(129L, proteins = 120L, per_group = 5L)
    quantified <- kit129Quiet(kit129LimpaQuantify(sim, "fastlimpadpc"))
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 77L)
    contrast <- "B_vs_A=groupB-groupA"
    run_in <- function(path) {
        manager <- kit129LimpaManager(quantified)
        source <- getDIAStateReference(manager)
        kit129Quiet(withDIAStateScope(manager, function() {
            if (identical(path, "process")) {
                run <- computeDIAReferenceRuvCovariates(manager, source, parameters)
                final <- run$value$final
                da <- computeDIAReferenceDaContrast(manager, final, kit129DaStep(manager, final, contrast), "limpa")$run$value
            } else {
                worker <- computeWorkflowDiaRuvCovariatesAction(manager, list(input_references = list(state = source), parameters = parameters), "whole")
                run <- list(value = hydrateSavedWorkflowEntity(manager$.context$store, worker$reference)$value)
                final <- run$value$final
                handoff <- computeWorkflowDiaContrastAction(manager, list(input_references = list(state = final),
                    parameters = c(kit129DaStep(manager, final, contrast), list(route = "limpa"))), "whole")
                da <- hydrateSavedWorkflowEntity(manager$.context$store, handoff$reference)$value
            }
            list(k = run$value$best_k, model = da$model, state = readDIAExecutionComponents(manager, final),
                table = kit129DaTable(manager, da, "Protein.Ids"))
        }))
    }
    process <- run_in("process")
    worker <- run_in("worker")
    expect_identical(worker$k, process$k)
    expect_identical(worker$table, process$table)
    expect_identical(process$model$model, "dpcDE")
    # the kernel called directly: the declared kernel in memory, then the DA entry point on its state
    saved <- list()
    store <- list(checkpoint = function(value, label, from = NULL) {
        saved[[label]] <<- value
        label
    })
    value <- kit129Quiet(diaRuvCovariatesOperation("declared")$kernel(quantified, parameters, store))
    expect_identical(value$best_k, process$k)
    kernel <- kit129Quiet(differentialAbundanceAnalysis(saved[[value$final]], contrasts_tbl = data.frame(contrasts = contrast),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
    expected <- kit129KernelTable(kernel, "Protein.Ids")
    expect_identical(process$table$Protein.Ids, expected$Protein.Ids)
    expect_equal(process$table[-1L], expected[-1L], tolerance = 1e-12)
})

test_that("(item 5) transcriptomics: the declared design in process, on the worker and by the kernel, one result", {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
    fixture <- kit129TrxOwner()
    owner <- fixture$owner
    design <- trxRuvDesignParameters("groups_all_features", 2L, "user")
    # the worker first, on the head (the RUVs step's source), with the logCPM product the tab's sequence delivers first
    logcpm <- kit129Quiet(withWorkflowStateScope(owner, function() trxRunProduct(owner, "trx_tmm_logcpm", getWorkflowStateReference(owner), list())))
    handoff <- omicsWorkerRun(fixture, "trx_ruv_grid", c(list(logcpm = logcpm$reference$entity_id), design))
    expect_identical(handoff$status, "complete", info = handoff$error)
    worker <- withWorkflowStateScope(owner, function() hydrateSavedWorkflowEntity(owner$.context$store, handoff$result$reference)$value)
    # in process, through the tab's wrapper
    run <- kit129Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = design))
    expect_identical(worker$ruv_design, run$grid$value$ruv_design)
    expect_identical(worker$chosen, run$grid$value$chosen)
    # the kernel directly, on the RUVs step's source: the choice, the factors, the fit and test halves
    filtered <- trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, "trx_qc_complete")))
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, run$logcpm$value$log_cpm))
    choice <- kit129Quiet(trxChooseRuvDesign(filtered, design))
    expect_identical(choice$chosen, run$grid$value$chosen)
    s <- kit129Quiet(trxRuvFactors(filtered$counts, log_cpm, filtered$design, choice$chosen$control_genes, choice$chosen$k, choice$ruv_design$method))
    da <- kit129Quiet(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner), list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1),
        fixture$workflow$contrasts_tbl))
    entry <- names(da$manifest$results)[[1L]]
    top <- trxDeclValue(owner, da$manifest$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table
    direct <- kit129Quiet(trxDaApply(trxDaFit(s$ruvs_counts, s$pdata, trxDaFormula(s$pdata)), entry))$results[[1L]]$top_table
    expect_identical(top, direct)
})

test_that("(item 5) lipidomics and metabolomics: in process, on the worker and by the kernel, one result", {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) {
            owner$job <- NULL
            list(status = status)
        }, .package = "MultiScholaR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    for (omic in c("lipidomics", "metabolomics")) {
        prefix <- if (identical(omic, "lipidomics")) "lipid" else "metab"
        long <- if (identical(omic, "lipidomics")) "da_lipids_long" else "da_metabolites_long"
        id <- if (identical(omic, "lipidomics")) "lipid_id" else "metabolite_id"
        binding <- omicsDeclarationBindings()[[omic]]
        key <- function(t) {
            t <- as.data.frame(t)[order(t$assay, t[[id]]), c("assay", id, "logFC", "P.Value", "fdr_qvalue")]
            rownames(t) <- NULL
            t
        }
        # in process
        a <- kit129OmicsOwner(omic)
        ruv <- kit129Quiet(runWorkflowOmicsReferenceRuv(binding, a$owner, kit129OmicsRuvParameters(), list()))
        process <- kit129Quiet(runWorkflowOmicsReferenceDaAnalysis(binding, a$owner, ruv$reference, settings, a$contrasts))
        # on the worker: the RUV step accepted on a controller, then the contrast's action on its state
        b <- kit129OmicsOwner(omic)
        controller <- omicsWorkerController(b)
        accepted <- kit129Quiet(omicsWorkerAction(b, controller, paste0(prefix, "_ruv"), kit129OmicsRuvParameters(),
            context = list(state_name = paste0(prefix, "_ruv_corrected"), config = list(), description = "RUV")))
        expect_identical(accepted$acceptance, "complete", info = omic)
        step <- omicsDaContrastSteps(binding, settings, b$contrasts)[[1L]]
        handoff <- omicsWorkerRun(b, step$operation, step$parameters)
        expect_identical(handoff$status, "complete", info = paste(omic, handoff$error))
        worker <- withWorkflowStateScope(b$owner, function() {
            value <- hydrateSavedWorkflowEntity(b$owner$.context$store, handoff$result$reference)$value
            do.call(rbind, lapply(value$assays, function(entry) as.data.frame(materializeWorkflowTable(b$owner$.context, entry$tables[[long]]))))
        })
        # the kernel directly
        direct <- kit129Quiet(runOmicsRuvFactors(a$object, "groups_all_features", 1L, "user"))
        kernel <- kit129Quiet(if (identical(omic, "lipidomics")) runLipidsDA(a$object, a$contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE,
            ruv_factors = direct$factors) else runMetabolitesDA(a$object, a$contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = direct$factors))
        expect_equal(key(process$results[[long]]), key(kernel[[long]]), tolerance = 1e-12, info = omic)
        expect_equal(key(worker), key(kernel[[long]]), tolerance = 1e-12, info = omic)
    }
})

# ---- Third pass item 3: the technical-replicate design and first-pass nulls on the run path --------------------------

# The classic kernel called directly on the normalised object with the design's replicate sets, then the DA entry point.
kit129ProtKernelDa <- function(normalised, factors, replicate_sets, contrast) {
    corrected <- kit129Quiet(removeRowsWithMissingValuesPercent(theObject = ruvIII_C_Varying(protRuvGroupReplicateObject(normalised),
        ruv_grouping_variable = replicate_sets, ruv_number_k = factors$k, ctrl = factors$controls)))
    state <- attachProtRuvModelFactors(corrected, normalised, factors)
    kernel <- kit129Quiet(differentialAbundanceAnalysis(state, contrasts_tbl = data.frame(contrasts = contrast), formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Group",
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
    kit129KernelTable(kernel, "Protein.Group")
}

# One declared classic run in process: the normalisation's compute half, then the contrast's.
kit129ProtRunPath <- function(parameters, contrast, object = makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))) {
    manager <- kit129ProtStore(object = object)
    source <- getDIAStateReference(manager)
    kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, source, parameters)
        final <- run$value$final
        da <- computeDIAReferenceDaContrast(manager, final, kit129DaStep(manager, final, contrast), "iq")$run$value
        list(value = run$value, state = readDIAExecutionComponents(manager, final),
            normalised = readDIAExecutionComponents(manager, run$value$normalized), table = kit129DaTable(manager, da, "Protein.Group"),
            method = workflowSourceRuvDaMethod(manager, final))
    }))
}

test_that("(third pass item 3) classic proteomics: the technical-replicate design on the run path equals the kernel; refused without technical replicates", {
    contrast <- "a=grouptreated-groupcontrol"
    declared <- c(kit129ProtParameters(), list(ruv_design = "technical_replicates", ruv_k_rule = "user", ruv_k_seed = NULL))
    out <- kit129ProtRunPath(declared, contrast)
    expect_identical(out$value$parameters$ruv_design, "technical_replicates")
    expect_identical(out$value$ruv_model[c("ruv_design", "replicate_structure", "ruv_grouping_variable", "k")],
        list(ruv_design = "technical_replicates", replicate_structure = "technical_replicates", ruv_grouping_variable = "tech_rep_group", k = 1L))
    expect_identical(out$state@args$ruvCovariates[c("ruv_design", "replicate_structure", "ruv_grouping_variable")],
        list(ruv_design = "technical_replicates", replicate_structure = "technical_replicates", ruv_grouping_variable = "tech_rep_group"))
    expect_identical(out$method, "ruv_factors_exact_se_trend_robust")
    factors <- kit129Quiet(computeProtRuvModelFactors(out$normalised, ruv_k = 1L, ruv_k_rule = "user", ruv_grouping_variable = "tech_rep_group"))
    expected <- kit129ProtKernelDa(out$normalised, factors, "tech_rep_group", contrast)
    expect_identical(out$table$Protein.Group, expected$Protein.Group)
    expect_equal(out$table[-1L], expected[-1L], tolerance = 1e-12)
    # the identity names the design: the groups' design on the same request is another product
    groups <- kit129ProtRunPath(utils::modifyList(declared, list(ruv_design = "groups_all_features")), contrast)
    expect_false(identical(groups$value$final, out$value$final))
    # the k rules the kernel defines: the separation rule is undefined with every protein a control
    automatic <- kit129ProtParameters("automatic", NULL)
    automatic[c("auto_percentage_min", "auto_percentage_max", "separation_metric", "k_penalty_weight", "adaptive_k_penalty")] <-
        list(5, 10, "max_difference", 0.5, TRUE)
    separation <- tryCatch(kit129ProtRunPath(c(automatic, list(ruv_design = "technical_replicates", ruv_k_rule = "sweep_bounded",
        ruv_k_seed = NULL)), contrast), error = function(e) e)
    expect_s3_class(separation, "ruv_k_separation_undefined")
    # a design without technical replicates is refused, not run as the groups' design
    plain <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    plain@technical_replicate_id <- NA_character_
    refused <- tryCatch(kit129ProtRunPath(declared, contrast, object = plain), error = function(e) e)
    expect_s3_class(refused, "workflow_execution_invalid")
    expect_match(conditionMessage(refused), "needs declared technical replicates")
})

test_that("(third pass item 3) classic proteomics: first-pass nulls on the run path equal the kernel, with the user's k and the separation rule", {
    contrast <- "a=grouptreated-groupcontrol"
    declared <- c(kit129ProtParameters(), list(ruv_design = "first_pass_nulls", ruv_k_rule = "user", ruv_k_seed = NULL))
    out <- kit129ProtRunPath(declared, contrast)
    # the fixture declares technical replicates: the replicate sets are the units (arm 23's structure)
    expect_identical(out$value$ruv_model[c("ruv_design", "ruv_controls", "ruv_grouping_variable", "k")],
        list(ruv_design = "first_pass_nulls", ruv_controls = "first_pass", ruv_grouping_variable = "tech_rep_group", k = 1L))
    expect_identical(out$state@args$ruvCovariates$ruv_design, "first_pass_nulls")
    factors <- kit129Quiet(computeProtRuvModelFactors(out$normalised, ruv_k = 1L, ruv_k_rule = "user", ruv_grouping_variable = "tech_rep_group",
        ruv_design = "first_pass_nulls"))
    expect_identical(sum(out$value$ruv_result$best_control_genes_index), sum(factors$controls))
    expect_lt(sum(factors$controls), length(factors$controls))
    expected <- kit129ProtKernelDa(out$normalised, factors, "tech_rep_group", contrast)
    expect_identical(out$table$Protein.Group, expected$Protein.Group)
    expect_equal(out$table[-1L], expected[-1L], tolerance = 1e-12)
    # the separation rule is defined for the first pass's controls: k is the sweep's estimate the run path hands the kernel
    automatic <- kit129ProtParameters("automatic", NULL)
    automatic[c("auto_percentage_min", "auto_percentage_max", "separation_metric", "k_penalty_weight", "adaptive_k_penalty")] <-
        list(5, 10, "max_difference", 0.5, TRUE)
    swept <- tryCatch(kit129ProtRunPath(c(automatic, list(ruv_design = "first_pass_nulls", ruv_k_rule = "sweep_bounded", ruv_k_seed = NULL)),
        contrast), ruv_k_above_limit = function(e) e)
    if (inherits(swept, "ruv_k_above_limit")) {
        # the sweep's k is above the units' limit on this fixture (decision M28): refused, never reduced
        expect_match(conditionMessage(swept), "limit")
    } else {
        expect_identical(swept$value$ruv_model$ruv_k_rule, "sweep_bounded")
        expect_identical(as.integer(swept$value$best_k), as.integer(swept$value$ruv_result$sweep_diagnostic$best_k))
    }
})

# One lipid/metab RUV step and its DA in process, and the kernel called directly, as tables keyed by assay and feature.
kit129OmicsRunPath <- function(omic, parameters, technical = FALSE, env = parent.frame()) {
    binding <- omicsDeclarationBindings()[[omic]]
    fixture <- kit129OmicsOwner(omic, env = env, technical = technical)
    ruv <- kit129Quiet(runWorkflowOmicsReferenceRuv(binding, fixture$owner, parameters, list()))
    da <- kit129Quiet(runWorkflowOmicsReferenceDaAnalysis(binding, fixture$owner, ruv$reference,
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$contrasts))
    list(fixture = fixture, ruv = ruv, da = da, binding = binding)
}
kit129OmicsKey <- function(omic, t) {
    id <- if (identical(omic, "lipidomics")) "lipid_id" else "metabolite_id"
    t <- as.data.frame(t)[order(t$assay, t[[id]]), c("assay", id, "logFC", "P.Value", "fdr_qvalue")]
    rownames(t) <- NULL
    t
}
kit129OmicsKernelDa <- function(omic, fixture, factors) {
    kit129Quiet(if (identical(omic, "lipidomics")) runLipidsDA(fixture$object, fixture$contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = factors) else runMetabolitesDA(fixture$object, fixture$contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = factors))
}

test_that("(third pass item 3) lipidomics and metabolomics: the technical-replicate design on the run path equals the kernel and records its structure", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    for (omic in c("lipidomics", "metabolomics")) {
        long <- if (identical(omic, "lipidomics")) "da_lipids_long" else "da_metabolites_long"
        prefix <- if (identical(omic, "lipidomics")) "lipid" else "metab"
        run <- kit129OmicsRunPath(omic, kit129OmicsRuvParameters("technical_replicates", grouping = NULL), technical = TRUE)
        value <- run$ruv$value
        expect_identical(value[c("ruv_design", "ruv_replicate_structure", "ruv_grouping_variable")],
            list(ruv_design = "technical_replicates", ruv_replicate_structure = "technical_replicates", ruv_grouping_variable = "extract"), info = omic)
        expect_identical(unname(unlist(value$best_k)), c(1L, 1L), info = omic)
        state <- withWorkflowStateScope(run$fixture$owner, function() readWorkflowExecutionComponents(run$fixture$owner, run$ruv$reference,
            fields = run$binding$fields))
        expect_identical(state@args$ruvFactorsDeclaration[c("ruv_design", "ruv_replicate_structure", "ruv_grouping_variable")],
            list(ruv_design = "technical_replicates", ruv_replicate_structure = "technical_replicates", ruv_grouping_variable = "extract"), info = omic)
        expect_identical(unname(vapply(state@args$ruv_factors, `[[`, "", "ruv_replicate_structure")), rep("technical_replicates", 2L), info = omic)
        expect_identical(unname(vapply(state@args$ruv_factors, `[[`, "", "ruv_design")), rep("technical_replicates", 2L), info = omic)
        expect_identical(run$da$manifest$results[[1L]]$model$ruv_da_method, "limma_ruv_factors_exact_se", info = omic)
        # the kernel directly: all features as controls, the declared units as the replicate sets
        direct <- kit129Quiet(runOmicsRuvFactors(run$fixture$object, "groups_all_features", 1L, "user", ruv_grouping_variable = "extract"))
        expect_equal(lapply(state@args$ruv_factors, `[[`, "W"), lapply(direct$factors, `[[`, "W"), tolerance = 0, info = omic)
        kernel <- kit129OmicsKernelDa(omic, run$fixture, direct$factors)
        expect_equal(kit129OmicsKey(omic, run$da$results[[long]]), kit129OmicsKey(omic, kernel[[long]]), tolerance = 1e-12, info = omic)
        # the record names the structure and its column
        record <- buildWorkflowAnalysisRecord(run$fixture$owner)
        step <- Filter(function(step) identical(step$operation, paste0(prefix, "_ruv")), record$steps)
        expect_length(step, 1L)
        expect_identical(step[[1L]]$outcomes[c("ruv_design", "ruv_replicate_structure", "ruv_grouping_variable")],
            list(ruv_design = "technical_replicates", ruv_replicate_structure = "technical_replicates", ruv_grouping_variable = "extract"), info = omic)
    }
    # refused without technical replicates, with another replicate-set column, and under the separation rule
    binding <- lipidomicsDeclarationBinding
    plain <- kit129OmicsOwner("lipidomics")
    expect_error(kit129Quiet(runWorkflowOmicsReferenceRuv(binding, plain$owner, kit129OmicsRuvParameters("technical_replicates", grouping = NULL),
        list())), "needs declared technical replicates", class = "workflow_execution_invalid")
    technical <- kit129OmicsOwner("lipidomics", technical = TRUE)
    expect_error(kit129Quiet(runWorkflowOmicsReferenceRuv(binding, technical$owner, kit129OmicsRuvParameters("technical_replicates"), list())),
        "not `group`", class = "workflow_execution_invalid")
    expect_error(kit129Quiet(runWorkflowOmicsReferenceRuv(binding, technical$owner, kit129OmicsRuvParameters("technical_replicates", 1L,
        "sweep_bounded", grouping = NULL), list())), class = "omics_ruv_factors_invalid")
})

test_that("(third pass item 3) lipidomics and metabolomics: first-pass nulls on the run path equal the kernel", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    for (omic in c("lipidomics", "metabolomics")) {
        long <- if (identical(omic, "lipidomics")) "da_lipids_long" else "da_metabolites_long"
        run <- kit129OmicsRunPath(omic, kit129OmicsRuvParameters("first_pass_nulls"))
        expect_identical(run$ruv$value[c("ruv_design", "ruv_replicate_structure", "ruv_grouping_variable")],
            list(ruv_design = "first_pass_nulls", ruv_replicate_structure = "group", ruv_grouping_variable = "group"), info = omic)
        direct <- kit129Quiet(runOmicsRuvFactors(run$fixture$object, "first_pass_nulls", 1L, "user"))
        expect_identical(run$ruv$value$control_count, lapply(direct$factors, `[[`, "control_count"), info = omic)
        expect_true(all(unlist(run$ruv$value$control_count) < unlist(lapply(direct$factors, `[[`, "complete_features"))), info = omic)
        kernel <- kit129OmicsKernelDa(omic, run$fixture, direct$factors)
        expect_equal(kit129OmicsKey(omic, run$da$results[[long]]), kit129OmicsKey(omic, kernel[[long]]), tolerance = 1e-12, info = omic)
    }
})

test_that("(third pass item 3) DIA with limpa and transcriptomics: first-pass nulls reach the run path's declared step", {
    # the tab's requests with the gate on carry the first pass's rule, and the request selects the declared declaration
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    dia <- diaRuvCovariatesRequestParameters(list(), "group", list(ruv_mode = "covariates", ruv_control_rule = "first_pass_drop_top_25",
        ruv_k_rule = "user", ruv_k = 1L))
    expect_identical(dia[c("ruv_control_rule", "ruv_k_rule", "ruv_k")], list(ruv_control_rule = "first_pass_drop_top_25", ruv_k_rule = "user", ruv_k = 1L))
    expect_identical(workflowRequestOperation(NULL, workflowOperationEntry("protein_ruv_covariates"), NULL, dia)$declaration$parameters,
        diaRuvCovariatesOperation("declared")$parameters)
    trx <- trxNormRuvRequest(c(trxRuvGridParameters(), list(ruv_design = "first_pass_nulls", ruv_k = 1L, ruv_k_rule = "user")))
    expect_null(trx$grid)
    expect_identical(trx$design$ruv_design, "first_pass_nulls")
    expect_identical(workflowRequestOperation(NULL, workflowOperationEntry("trx_ruv_grid"), NULL, c(list(logcpm = "x"), trx$design))$declaration$parameters,
        c("logcpm", trxRuvDesignParameterNames))
    # DIA with limpa end to end: the compute half equals the declared kernel called directly
    skip_if_not_installed("limpa")
    skip_if_not_installed("fastlimpadpc")
    sim <- kit129LimpaSimulate(129L, proteins = 120L, per_group = 5L)
    quantified <- kit129Quiet(kit129LimpaQuantify(sim, "fastlimpadpc"))
    manager <- kit129LimpaManager(quantified)
    contrast <- "B_vs_A=groupB-groupA"
    process <- kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceRuvCovariates(manager, getDIAStateReference(manager), dia)
        final <- run$value$final
        da <- computeDIAReferenceDaContrast(manager, final, kit129DaStep(manager, final, contrast), "limpa")$run$value
        list(value = run$value, table = kit129DaTable(manager, da, "Protein.Ids"))
    }))
    expect_identical(process$value$ruv_control_rule, "first_pass_drop_top_25")
    saved <- list()
    store <- list(checkpoint = function(value, label, from = NULL) {
        saved[[label]] <<- value
        label
    })
    value <- kit129Quiet(diaRuvCovariatesOperation("declared")$kernel(quantified, dia, store))
    expect_identical(value$best_k, process$value$best_k)
    kernel <- kit129Quiet(differentialAbundanceAnalysis(saved[[value$final]], contrasts_tbl = data.frame(contrasts = contrast),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
    expected <- kit129KernelTable(kernel, "Protein.Ids")
    expect_identical(process$table$Protein.Ids, expected$Protein.Ids)
    expect_equal(process$table[-1L], expected[-1L], tolerance = 1e-12)
})

test_that("(fourth pass item 1) lipidomics and metabolomics: the declared RUV step's identity states the replicate structure by design", {
    for (binding in list(lipidomicsDeclarationBinding, metabolomicsDeclarationBinding)) {
        fixed <- omicsRuvOperation(binding, "declared")$fixed
        expect_match(fixed$replicate_sets, "the groups (RUV-III)", fixed = TRUE, info = binding$prefix)
        expect_match(fixed$replicate_sets, "under the design technical_replicates the declared technical replicates", fixed = TRUE,
            info = binding$prefix)
        # today's declaration has no fixed replicate text
        expect_null(omicsRuvOperation(binding)$fixed$replicate_sets, info = binding$prefix)
    }
})

test_that("(fourth pass item 2) the panel's choice runs in process and on the worker as one product; the head's technical replicates are read", {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) {
            owner$job <- NULL
            list(status = status)
        }, .package = "MultiScholaR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    plain <- kit129OmicsOwner("lipidomics")
    technical <- kit129OmicsOwner("lipidomics", technical = TRUE)
    expect_false(omicNormRuvHeadTechnicalReplicates(plain$owner))
    expect_true(omicNormRuvHeadTechnicalReplicates(technical$owner))
    expect_false(omicNormRuvHeadTechnicalReplicates(NULL))
    input <- list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_user = 1, ruv_k_seed = RUV_FACTOR_COUNT_SEED,
        ruv_control_list = "", ruv_mode = "automatic", ruv_grouping_variable = "group")
    settings <- omicNormRunSettings(input, list(color_var = "group", shape_var = "group"), list(assay_names = NULL), list(),
        lipidomicsDeclarationBinding, ruvDesign = omicNormRuvDesignChoice(input, omicNormRuvRunTable("lipidomics")))
    worker <- omicsWorkerRun(plain, "lipid_ruv", omicsRuvStep(lipidomicsDeclarationBinding, settings$ruv_parameters)$parameters)
    expect_identical(worker$status, "complete", info = worker$error)
    ruv <- kit129Quiet(runWorkflowOmicsReferenceRuv(lipidomicsDeclarationBinding, plain$owner, settings$ruv_parameters, list()))
    expect_identical(ruv$run$reference$identity, worker$result$reference$identity)
    expect_identical(ruv$value$ruv_design, "groups_all_features")
    expect_identical(unname(unlist(ruv$value$best_k)), c(1L, 1L))
})
