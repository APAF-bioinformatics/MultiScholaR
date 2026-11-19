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
    # (fourth pass item 3) what the panel shows of this run: k per assay, the count beside the typed k, the sets, the seed
    result <- omicsRuvDesignResult(ruv$value)
    expect_identical(result[c("design", "k_rule", "seed", "replicate_structure")], list(design = "groups_all_features", k_rule = "user",
        seed = RUV_FACTOR_COUNT_SEED, replicate_structure = "group"))
    expect_identical(lapply(result$assays, `[[`, "k"), list(pos = 1L, neg = 1L))
    expect_identical(lapply(result$assays, `[[`, "k_estimate"), ruv$value$ruv_k_estimate)
})

# ---- Fourth pass item 4: the factor-tied-to-groups check on the run's W, its threshold and the record ------------------

test_that("(fourth pass item 4) the thresholds table is the calibration's: validated designs only, with their sources", {
    thresholds <- ruvGroupAssociationThresholds()
    key <- paste(thresholds$omic, thresholds$ruv_design, thresholds$ruv_k_rule)
    # decision M27: the recommended design in classic proteomics, transcriptomics and lipid/metab; classic supplied
    # controls and first-pass nulls (arm 5, the separation rule); lipid/metab supplied controls. None in DIA with limpa.
    expect_setequal(key, c("proteomics groups_all_features num_sv_be", "proteomics supplied_controls num_sv_be",
        "proteomics first_pass_nulls sweep_bounded", "transcriptomics groups_all_features num_sv_be",
        "lipidomics groups_all_features num_sv_be", "metabolomics groups_all_features num_sv_be",
        "lipidomics supplied_controls num_sv_be", "metabolomics supplied_controls num_sv_be"))
    expect_false(any(thresholds$omic == "dia_limpa"))
    # the shares of the validation runs as the results files report them (results-20261118.md, followup-results-20261118.md)
    expect_identical(thresholds$false_alarms, c(0.051, 0.043, 0.047, 0.091, 0.095, 0.095, 0.060, 0.060))
    expect_identical(thresholds$detection, c(1.000, 1.000, 0.989, 1.000, 0.929, 0.929, 0.915, 0.915))
    expect_true(all(thresholds$false_alarms <= 0.10 & thresholds$detection >= 0.80))
    # each threshold is its source file's row, to the digits written there
    root <- testthat::test_path("..", "..")
    for (i in seq_len(nrow(thresholds))) {
        file <- file.path(root, thresholds$source[[i]])
        if (!file.exists(file)) next
        table <- utils::read.delim(file, stringsAsFactors = FALSE)
        expect_identical(thresholds$threshold[[i]], table$threshold[table$arm == thresholds$arm[[i]]], info = key[[i]])
    }
})

test_that("(fourth pass item 4) the statistic is the calibration's; the check gives a verdict where validated and none elsewhere", {
    set.seed(4)
    group <- rep(c("A", "B", "C"), each = 4L)
    W <- cbind(W1 = stats::rnorm(12L), W2 = stats::rnorm(12L) + (group == "B"))
    # the simulators' association mode (association.R, `factorAssociations`)
    expected <- vapply(1:2, function(j) sqrt(summary(stats::lm(W[, j] ~ factor(group)))$r.squared), numeric(1L))
    expect_identical(ruvFactorGroupAssociations(W, group), expected)
    largest <- max(expected)
    check <- function(...) ruvGroupAssociationCheck(W, group, ...)
    flagged <- check("lipidomics", "groups_all_features", "num_sv_be")
    expect_identical(flagged[c("ruv_group_association", "ruv_group_association_threshold", "ruv_group_association_flag",
        "ruv_group_association_outcome")], list(ruv_group_association = largest, ruv_group_association_threshold = 0.088611623,
        ruv_group_association_flag = largest > 0.088611623, ruv_group_association_outcome = if (largest > 0.088611623) "detectable" else "not_detectable"))
    expect_match(flagged$ruv_group_association_source, "threshold-by-design-S4.tsv (simulator S4, arm 7; set on run 44203265622, validated on run 44227129836)",
        fixed = TRUE)
    low <- ruvGroupAssociationCheck(W * 0 + rep(c(1, -1), 6L), group, "lipidomics", "groups_all_features", "num_sv_be")
    expect_identical(low$ruv_group_association_outcome, "not_detectable")
    expect_false(low$ruv_group_association_flag)
    for (args in list(list("dia_limpa", NA_character_, "num_sv_be"), list("proteomics", "groups_all_features", "user"),
        list("proteomics", "technical_replicates", "num_sv_be"), list("proteomics", "first_pass_nulls", "num_sv_be"),
        list("transcriptomics", "supplied_controls", "num_sv_be"))) {
        uncalibrated <- do.call(check, args)
        expect_identical(uncalibrated[c("ruv_group_association", "ruv_group_association_flag", "ruv_group_association_outcome")],
            list(ruv_group_association = largest, ruv_group_association_flag = NA, ruv_group_association_outcome = "no_limit"),
            info = paste(unlist(args), collapse = " "))
    }
    expect_identical(ruvGroupAssociationCheck(NULL, group, "lipidomics", "groups_all_features", "num_sv_be")$ruv_group_association_outcome,
        "no_factors")
    # factors whose samples the design does not name are not checked
    expect_identical(ruvGroupAssociationCheck(W, ruvFactorGroups(W, data.frame(s = sprintf("s%d", 12:1), g = rev(group)), "s", "g"),
        "lipidomics", "groups_all_features", "num_sv_be")$ruv_group_association_outcome, "not_computed")
    rownames(W) <- sprintf("s%d", 1:12)
    expect_identical(ruvFactorGroups(W, data.frame(s = sprintf("s%d", 12:1), g = rev(group)), "s", "g"), group)
})

test_that("(sixth pass M35) with a typed k the check gives a verdict only when the typed k equals the estimated number", {
    set.seed(4)
    group <- rep(c("A", "B"), each = 6L)
    W <- cbind(W1 = stats::rnorm(12L) + (group == "B"))
    largest <- max(ruvFactorGroupAssociations(W, group))
    estimated <- ruvGroupAssociationCheck(W, group, "proteomics", "groups_all_features", "num_sv_be", 1L, 1L)
    equal <- ruvGroupAssociationCheck(W, group, "proteomics", "groups_all_features", "user", 1L, 1L)
    # the calibrated calculation: the estimated number's limit and figures, and the record says so
    expect_identical(equal[setdiff(names(equal), "ruv_group_association_basis")], estimated[setdiff(names(estimated), "ruv_group_association_basis")])
    expect_identical(equal$ruv_group_association_basis, "typed_k_equals_estimated_number")
    expect_identical(estimated$ruv_group_association_basis, "estimated_number")
    expect_identical(equal$ruv_group_association_outcome, if (largest > 0.058077853) "detectable" else "not_detectable")
    for (case in list(list(k = 2L, estimate = 1L, basis = "typed_k_differs_from_estimated_number"),
        list(k = 1L, estimate = NA_integer_, basis = "typed_k_without_estimated_number"))) {
        other <- ruvGroupAssociationCheck(W, group, "proteomics", "groups_all_features", "user", case$k, case$estimate)
        expect_identical(other[c("ruv_group_association", "ruv_group_association_threshold", "ruv_group_association_flag",
            "ruv_group_association_outcome", "ruv_group_association_basis", "ruv_group_association_fdp_06")],
            list(ruv_group_association = largest, ruv_group_association_threshold = NA_real_, ruv_group_association_flag = NA,
                ruv_group_association_outcome = "no_limit", ruv_group_association_basis = case$basis, ruv_group_association_fdp_06 = NA_real_),
            info = case$basis)
    }
    # first-pass nulls has a limit with the separation rule only: a typed k equal to the estimate has none
    expect_identical(ruvGroupAssociationCheck(W, group, "proteomics", "first_pass_nulls", "user", 1L, 1L)$ruv_group_association_outcome,
        "no_limit")
})

test_that("(sixth pass M36) each limit's measured consequences are its approach's own arm, read from the scored tables", {
    consequences <- ruvGroupAssociationConsequences()
    thresholds <- ruvGroupAssociationThresholds()
    key <- function(x) paste(x$omic, x$ruv_design, x$ruv_k_rule)
    # every validated limit has its approach's figures
    expect_true(all(key(thresholds) %in% key(consequences)))
    expect_identical(anyDuplicated(key(consequences)), 0L)
    # the arm each shipped combination reproduces
    expect_identical(consequences$arm[match(key(thresholds), key(consequences))], thresholds$arm)
    root <- testthat::test_path("..", "..")
    read <- function(path) {
        file <- file.path(root, "docs/qa/workflow-execution/kit/kit-124/scored", path)
        if (file.exists(file)) utils::read.delim(file, stringsAsFactors = FALSE)
    }
    figure <- function(table, simulator, cell, arm) {
        hit <- table$mean_fdp[table$simulator == simulator & table$cell == cell & table$arm == arm & table$subset == "full" & !table$void]
        if (length(hit)) hit else NA_real_
    }
    for (i in seq_len(nrow(consequences))) {
        row <- consequences[i, ]
        final <- read(sprintf("final-20261118/%s/k124-cells.tsv", tolower(row$simulator)))
        own <- if (row$arm == 24L) read("newarms-gha/k124-cells.tsv") else final
        if (is.null(final) || is.null(own)) next
        expect_equal(c(row$fdp_0, row$fdp_03, row$fdp_06), c(figure(own, row$simulator, "C0", row$arm), figure(own, row$simulator, "C5", row$arm),
            figure(own, row$simulator, "S3", row$arm)), tolerance = 1e-6, info = key(row))
        expect_equal(row$uncorrected_06, figure(final, row$simulator, "S3", 1L), tolerance = 1e-6, info = key(row))
    }
    # M36's quoted figures for the recommended design
    recommended <- consequences[consequences$ruv_design == "groups_all_features", ]
    expect_identical(sprintf("%.3f", unlist(recommended[recommended$omic == "proteomics", c("fdp_0", "fdp_03", "fdp_06", "uncorrected_06")])),
        c("0.054", "0.056", "0.062", "0.593"))
    expect_identical(sprintf("%.3f", unlist(recommended[recommended$omic == "lipidomics", c("fdp_0", "fdp_03", "fdp_06", "uncorrected_06")])),
        c("0.057", "0.072", "0.130", "0.642"))
    # a check against a limit records its approach's figures
    group <- rep(c("A", "B"), each = 6L)
    check <- ruvGroupAssociationCheck(cbind(W1 = rep(c(1, -1), 6L)), group, "transcriptomics", "groups_all_features", "num_sv_be")
    expect_identical(unlist(check[c("ruv_group_association_fdp_0", "ruv_group_association_fdp_03", "ruv_group_association_fdp_06",
        "ruv_group_association_uncorrected_06")], use.names = FALSE), c(0.048062, 0.0480989, 0.0524544, 0.450032))
})

test_that("(fourth pass item 4) each omic's step writes the check on its W into the record; the panel shows it", {
    skip_if_not_installed("sva")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    # lipidomics: the estimated count on the recommended design, a validated combination, per assay
    fixture <- kit129OmicsOwner("lipidomics")
    ruv <- kit129Quiet(runWorkflowOmicsReferenceRuv(lipidomicsDeclarationBinding, fixture$owner,
        kit129OmicsRuvParameters(ruv_k = NULL, ruv_k_rule = "num_sv_be"), list()))
    state <- withWorkflowStateScope(fixture$owner, function() readWorkflowExecutionComponents(fixture$owner, ruv$reference,
        fields = lipidomicsDeclarationBinding$fields))
    for (assay in names(ruv$value$factors)) {
        record <- ruv$value$factors[[assay]]
        if (is.null(record$W)) next
        groups <- as.character(fixture$object@design_matrix$group[match(rownames(record$W), fixture$object@design_matrix$Sample_ID)])
        expect_identical(record$ruv_group_association$ruv_group_association, max(ruvFactorGroupAssociations(record$W, groups)), info = assay)
        expect_identical(record$ruv_group_association$ruv_group_association_threshold, 0.088611623, info = assay)
        expect_true(record$ruv_group_association$ruv_group_association_outcome %in% c("detectable", "not_detectable"), info = assay)
    }
    step <- Filter(function(step) identical(step$operation, "lipid_ruv"), buildWorkflowAnalysisRecord(fixture$owner)$steps)[[1L]]
    expect_identical(step$outcomes$ruv_group_association_outcome_per_assay,
        lapply(ruv$value$factors, function(record) record$ruv_group_association$ruv_group_association_outcome))
    expect_identical(step$outcomes$ruv_group_association_threshold_per_assay,
        lapply(ruv$value$factors, function(record) record$ruv_group_association$ruv_group_association_threshold))
    shown <- as.character(omicNormRuvResultTags(omicsRuvDesignResult(ruv$value), omicNormRuvPresets(omicGuiConfig("lipidomics")$norm$ruv_design),
        omicNormRuvRunTable("lipidomics")$texts))
    if (any(vapply(ruv$value$factors, function(record) !is.null(record$W), NA))) expect_match(shown, "values above 0.089|detection limit 0.089")
    # classic proteomics: the user's k is not a validated combination: the value without a verdict, recorded
    manager <- kit129ProtStore()
    source <- getDIAStateReference(manager)
    declared <- c(kit129ProtParameters(), list(ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_seed = NULL))
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, source, declared)
        record <- list(config = list(ruv_k = run$value$best_k), data = list(products = list(normalization_science = run$reference)))
        list(run = run, outcomes = workflowRecordStepOutcomes(manager$.context, "ruv_corrected", "normalization-ruv", record, list()))
    }))
    association <- out$run$value$ruv_result$ruv_group_association
    expect_identical(association$ruv_group_association_outcome, "no_limit")
    expect_true(is.finite(association$ruv_group_association))
    expect_identical(out$outcomes[ruvGroupAssociationRecordFields], association[ruvGroupAssociationRecordFields])
    shown <- as.character(omicNormRuvResultTags(protNormRuvDesignResult(out$run$value),
        omicNormRuvPresets(omicGuiConfig("proteomics")$norm$ruv_design), omicNormRuvRunTable("proteomics", "LFQ")$texts))
    expect_match(shown, sprintf("Largest association of an estimated unwanted factor with your groups: %.3f", association$ruv_group_association),
        fixed = TRUE)
    # today's design records nothing of the check
    today <- kit129Quiet(withDIAStateScope(manager, function() computeDIAReferenceNormalization(manager, source, kit129ProtParameters())))
    expect_null(today$value$ruv_result$ruv_group_association)
    expect_true(isTRUE(protNormRuvDesignResult(today$value)$legacy))
})

test_that("(fourth pass item 4) transcriptomics: the RUVs step's outcomes carry the check of its W", {
    fixture <- kit129TrxOwner()
    owner <- fixture$owner
    run <- kit129Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("groups_all_features", 2L, "user")))
    outcomes <- run$value$outcomes
    expect_identical(outcomes$ruv_group_association_outcome, "no_limit")
    design <- withWorkflowStateScope(owner, function() readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, owner)$design_matrix)
    W <- design[, startsWith(names(design), "W_"), drop = FALSE]
    expect_identical(outcomes$ruv_group_association, max(ruvFactorGroupAssociations(W, as.character(design$group))))
    result <- trxRuvDesignResult(run$grid$value, outcomes)
    expect_identical(result$assays$all$association$ruv_group_association_outcome, "no_limit")
})

# ---- Fourth pass item 5: the transcriptomics tab under a declared design ---------------------------------------------

test_that("(fourth pass item 5) transcriptomics: a declared design's run publishes what exists for it and shows no grid", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    fixture <- kit129TrxOwner()
    owner <- fixture$owner
    root <- withr::local_tempdir()
    paths <- list(gene_qc_dir = file.path(root, "qc"), publication_graphs_dir = file.path(root, "graphs"))
    norm_data <- new.env()
    norm_data$ruv_optimization_results <- list(counts = list(success = TRUE, best_k = 3L))
    lines <- character()
    outcome <- new.env()
    input <- list(ruv_mode = "automatic", ruv_design = "groups_all_features", ruv_k_rule = "user", ruv_k_user = 2, ruv_k_seed = 7)
    kit129Quiet(runTrxNormInProcess(transcriptomicsDeclarationBinding, owner, trxNormGridParameters(input), norm_data, list(),
        function(line) lines <<- c(lines, line), paths, finish = function() outcome$finished <- TRUE,
        failed = function(error) outcome$error <- conditionMessage(error)))
    expect_null(outcome$error)
    expect_true(isTRUE(outcome$finished))
    # the grid of the earlier run is not shown; the panel's result is this run's
    expect_identical(norm_data$ruv_optimization_results, list())
    expect_identical(norm_data$ruv_design_result$design, "groups_all_features")
    expect_identical(norm_data$ruv_design_result$assays$all$k, 2L)
    # the RUV publication writes the corrected display matrix, and no canonical-correlation panels or grid workbook
    display <- kit129Quiet(runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, owner, "trx_ruv_publication", list(), paths))
    expect_true(isTRUE(display$saved), info = display$message)
    expect_true(file.exists(file.path(paths$gene_qc_dir, "ruv_normalised_log_values.tsv")))
    expect_false(any(grepl("canonical_correlation", list.files(paths$publication_graphs_dir))))
    # the stage figures of the declared state are there
    expect_true(file.exists(file.path(paths$gene_qc_dir, "trx_ruvs_pca.png")))
})

test_that("(fourth pass item 5) the RUV QC tab shows the legacy grid's frame only for a legacy run", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    bus <- shiny::reactiveValues(state_manager = NULL, design_matrix = NULL, config_list = list(),
        tab_status = list(normalization = "pending"), state_update_trigger = NULL)
    suppressWarnings(suppressMessages(shiny::testServer(mod_omic_norm_server, args = list(id = "norm", workflow_data = bus,
        experiment_paths = list(), omic_type = "transcriptomics", experiment_label = "kit129",
        config = omicGuiConfig("transcriptomics")$norm), {
        html <- function() paste(as.character(output$ruv_qc_ui$html), collapse = "")
        session$flushReact()
        expect_match(html(), "Assay: counts", fixed = TRUE)
        omicNormSetRuvDesignResult(session$returned, omicNormRuvDesignResult("groups_all_features", "user", 7L, "group",
            list(all = list(k = 1L, applied = TRUE))))
        session$flushReact()
        expect_false(grepl("Assay: counts", html(), fixed = TRUE))
        omicNormSetRuvDesignResult(session$returned, omicNormRuvLegacyResult())
        session$flushReact()
        expect_match(html(), "Assay: counts", fixed = TRUE)
    })))
})

# ---- Fifth pass item 3: supplied control lists on the classic and lipid/metab run paths -------------------------------

# The run's `ruv_control_ids_unmatched` warning (its fields), the other output and conditions muffled; NULL without one.
kit129UnmatchedWarning <- function(expr) {
    tryCatch({
        utils::capture.output(suppressMessages(withCallingHandlers(expr, warning = function(w) {
            if (!inherits(w, "ruv_control_ids_unmatched")) invokeRestart("muffleWarning")
        })))
        NULL
    }, ruv_control_ids_unmatched = function(w) list(matched = w$matched, listed = w$listed, unmatched = w$unmatched))
}

test_that("(fifth pass item 3) classic proteomics: a supplied list on the run path equals the kernel, named by its hash and count", {
    contrast <- "a=grouptreated-groupcontrol"
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    ids <- as.character(object@protein_quant_table[[object@protein_id_column]])
    listed <- c(ids[seq(1L, length(ids), by = 3L)], "NOT_A_PROTEIN")
    declared <- c(kit129ProtParameters(), list(ruv_design = "supplied_controls", ruv_k_rule = "user", ruv_k_seed = NULL,
        ruv_control_ids = rev(listed)))
    out <- kit129ProtRunPath(declared, contrast, object = object)
    expect_identical(out$value$ruv_model[c("ruv_design", "ruv_controls", "ruv_grouping_variable", "k")],
        list(ruv_design = "supplied_controls", ruv_controls = "supplied", ruv_grouping_variable = "group", k = 1L))
    expect_identical(out$state@args$ruvCovariates$ruv_design, "supplied_controls")
    expect_identical(out$method, "ruv_factors_exact_se_trend_robust")
    # the list by section 1e's names; the request's identity is its hash and the number listed, not the ids
    expect_identical(out$value$ruv_result[c("ruv_control_rule", "ruv_control_list_sha256", "ruv_control_list_count",
        "ruv_control_list_unmatched")], list(ruv_control_rule = "supplied_list", ruv_control_list_sha256 = trxControlListSha256(listed),
        ruv_control_list_count = length(listed), ruv_control_list_unmatched = "NOT_A_PROTEIN"))
    expect_identical(out$value$parameters[RUV_SUPPLIED_LIST_PARAMETERS], list(ruv_control_list_sha256 = trxControlListSha256(listed),
        ruv_control_count = length(listed)))
    expect_false("ruv_control_ids" %in% names(out$value$parameters))
    # the kernel directly: the listed proteins of the table as controls, the groups as replicate sets (arm 8)
    factors <- kit129Quiet(computeProtRuvModelFactors(out$normalised, ruv_k = 1L, ruv_k_rule = "user", ruv_grouping_variable = "group",
        controls = ids %in% listed, ruv_design = "supplied_controls"))
    expect_identical(out$value$ruv_model$control_count, factors$control_count)
    expected <- kit129ProtKernelDa(out$normalised, factors, "group", contrast)
    expect_identical(out$table$Protein.Group, expected$Protein.Group)
    expect_equal(out$table[-1L], expected[-1L], tolerance = 1e-12)
    # what the identity reads: one list in another order (an id repeated) is one identity, another list another one
    identity <- function(ids) {
        chosen <- workflowRequestOperation(NULL, workflowOperationEntry("normalization-ruv"), NULL,
            utils::modifyList(declared, list(ruv_control_ids = ids)))
        chosen$parameters[chosen$declaration$parameters]
    }
    expect_identical(identity(c(listed, listed[[1L]])), identity(rev(listed)))
    expect_false(identical(identity(listed[-1L]), identity(listed)))
    expect_false("ruv_control_ids" %in% names(identity(listed)))
    # the unknown identifiers are said with their counts
    manager <- kit129ProtStore(object = object)
    warned <- kit129UnmatchedWarning(withDIAStateScope(manager, function() computeDIAReferenceNormalization(manager,
        getDIAStateReference(manager), declared)))
    expect_identical(warned[c("matched", "listed", "unmatched")], list(matched = length(listed) - 1L, listed = length(listed),
        unmatched = "NOT_A_PROTEIN"))
    # refused: a list that names no protein; too few controls (the kernel's rule); a list under another design; today's
    # declaration handed a list
    expect_error(kit129ProtRunPath(utils::modifyList(declared, list(ruv_control_ids = c("X1", "X2"))), contrast, object = object),
        "none of the 2 listed", class = "workflow_execution_invalid")
    expect_error(kit129ProtRunPath(utils::modifyList(declared, list(ruv_control_ids = ids[[1L]])), contrast, object = object),
        class = "prot_ruv_factors_unavailable")
    expect_error(kit129ProtRunPath(utils::modifyList(declared, list(ruv_design = "groups_all_features")), contrast, object = object),
        "ruv_control_ids names the controls", class = "workflow_execution_invalid")
    expect_error(kit129ProtRunPath(c(kit129ProtParameters(), list(ruv_control_ids = listed)), contrast, object = object),
        class = "workflow_execution_invalid")
})

test_that("(fifth pass item 3) lipidomics and metabolomics: one supplied list on the run path, per assay, equals the kernel", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    controls <- list(pos = sprintf("pos_f%03d", 11:30), neg = sprintf("neg_f%03d", 21:40))
    listed <- c(controls$neg, controls$pos, "absent_f999")
    parameters <- utils::modifyList(kit129OmicsRuvParameters("supplied_controls"), list(ruv_control_ids = listed))
    for (omic in c("lipidomics", "metabolomics")) {
        long <- if (identical(omic, "lipidomics")) "da_lipids_long" else "da_metabolites_long"
        run <- kit129OmicsRunPath(omic, parameters)
        value <- run$ruv$value
        expect_identical(value[c("ruv_design", "ruv_control_rule", "ruv_control_list_sha256", "ruv_control_list_count",
            "ruv_control_list_unmatched")], list(ruv_design = "supplied_controls", ruv_control_rule = "supplied_list",
            ruv_control_list_sha256 = trxControlListSha256(listed), ruv_control_list_count = length(listed),
            ruv_control_list_unmatched = "absent_f999"), info = omic)
        expect_identical(value$parameters[RUV_SUPPLIED_LIST_PARAMETERS], list(ruv_control_list_sha256 = trxControlListSha256(listed),
            ruv_control_count = length(listed)), info = omic)
        expect_false("ruv_control_ids" %in% names(value$parameters), info = omic)
        # the kernel directly, each assay its own listed features
        direct <- kit129Quiet(runOmicsRuvFactors(run$fixture$object, "supplied_controls", 1L, "user", ruv_control_ids = controls))
        expect_identical(value$control_count, lapply(direct$factors, `[[`, "control_count"), info = omic)
        expect_equal(lapply(value$factors, `[[`, "W"), lapply(direct$factors, `[[`, "W"), tolerance = 0, info = omic)
        kernel <- kit129OmicsKernelDa(omic, run$fixture, direct$factors)
        expect_equal(kit129OmicsKey(omic, run$da$results[[long]]), kit129OmicsKey(omic, kernel[[long]]), tolerance = 1e-12, info = omic)
        # what the panel shows: the identifiers listed and found, the controls used per assay
        shown <- omicsRuvDesignResult(value)
        expect_identical(shown$control_list, list(listed = length(listed), unmatched = 1L), info = omic)
        expect_identical(lapply(shown$assays, `[[`, "control_count"), value$control_count, info = omic)
    }
    binding <- lipidomicsDeclarationBinding
    # an assay that holds none of the listed features is refused; too few controls in an assay are the kernel's refusal
    only_pos <- kit129OmicsOwner("lipidomics")
    expect_error(kit129Quiet(runWorkflowOmicsReferenceRuv(binding, only_pos$owner, utils::modifyList(parameters,
        list(ruv_control_ids = controls$pos)), list())), "assay 'neg': none of the 20 listed", class = "workflow_execution_invalid")
    one_each <- kit129OmicsOwner("lipidomics")
    expect_error(kit129Quiet(runWorkflowOmicsReferenceRuv(binding, one_each$owner, utils::modifyList(parameters,
        list(ruv_control_ids = c("pos_f011", "neg_f021"))), list())), class = "ruv_k_above_limit")
    # the unknown identifier is said with the counts
    warned <- kit129OmicsOwner("lipidomics")
    warning <- kit129UnmatchedWarning(runWorkflowOmicsReferenceRuv(binding, warned$owner, parameters, list()))
    expect_identical(warning[c("matched", "listed", "unmatched")], list(matched = length(listed) - 1L, listed = length(listed),
        unmatched = "absent_f999"))
})

# ---- Fifth pass item 4: the separation rule on the run path where it is defined ---------------------------------------

test_that("(fifth pass item 4) transcriptomics: the separation rule on the declared run path reads the design's own controls", {
    fixture <- kit129TrxOwner()
    owner <- fixture$owner
    filtered <- trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, "trx_qc_complete")))
    genes <- rownames(filtered$counts)
    for (design in list(trxRuvDesignParameters("first_pass_nulls", ruv_k_rule = "sweep_bounded"),
        trxRuvDesignParameters("supplied_controls", ruv_k_rule = "sweep_bounded", ruv_control_ids = genes[seq(1L, length(genes), by = 7L)]))) {
        run <- kit129Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = design))
        choice <- run$grid$value
        expect_identical(choice$ruv_design[c("id", "k_rule")], list(id = design$ruv_design, k_rule = "sweep_bounded"))
        expect_s3_class(choice$separation_plot, "grob")
        # arms 22 and 24 of simulator S3: trxCancor on the TMM logCPM with the design's controls, trxFindBestK to the cap
        log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, run$logcpm$value$log_cpm))
        controls <- kit129Quiet(trxRuvDesignControls(filtered$counts, filtered$design, design$ruv_design, design$ruv_control_ids))
        k <- trxFindBestK(trxCancor(Y = t(log_cpm), X = stats::model.matrix(~ group, data = filtered$design)[, -1, drop = FALSE],
            ctl = rownames(log_cpm) %in% controls$control_genes), max_k = if (ncol(filtered$counts) < 15) 3L else 4L)
        expect_identical(choice$chosen$k, k, info = design$ruv_design)
        expect_identical(choice$chosen$control_genes, controls$control_genes, info = design$ruv_design)
        # the panel shows the rule's plot
        expect_s3_class(trxRuvDesignResult(choice, run$value$outcomes)$separation_plot, "grob")
    }
    # undefined with all genes as controls
    expect_error(kit129Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("groups_all_features",
        ruv_k_rule = "sweep_bounded"))), class = "ruv_k_separation_undefined")
})

test_that("(fifth pass item 4) classic proteomics: with a supplied list the separation rule reads the list (arm 22)", {
    contrast <- "a=grouptreated-groupcontrol"
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    ids <- as.character(object@protein_quant_table[[object@protein_id_column]])
    listed <- ids[seq(1L, length(ids), by = 2L)]
    automatic <- kit129ProtParameters("automatic", NULL)
    automatic[c("auto_percentage_min", "auto_percentage_max", "separation_metric", "k_penalty_weight", "adaptive_k_penalty")] <-
        list(5, 10, "max_difference", 0.5, TRUE)
    declared <- c(automatic, list(ruv_design = "supplied_controls", ruv_k_rule = "sweep_bounded", ruv_k_seed = NULL, ruv_control_ids = listed))
    out <- tryCatch(kit129ProtRunPath(declared, contrast, object = object), ruv_k_refused = function(e) e)
    separated <- kit129Quiet(tryCatch(protRuvSeparationRuleK(kit129ProtRunPath(c(kit129ProtParameters("skip", NULL)), contrast,
        object = object)$normalised, ids %in% listed, "group"), ruv_k_refused = function(e) e))
    if (inherits(out, "ruv_k_refused")) {
        # the rule's k is refused on this fixture (no k, or above the limit): refused, never reduced
        expect_true(inherits(separated, "ruv_k_refused") || separated$k > 0L)
    } else {
        expect_identical(out$value$ruv_model$ruv_k_rule, "sweep_bounded")
        expect_identical(as.integer(out$value$best_k), separated$k)
        expect_s3_class(out$value$ruv_result$best_cancor_plot, "grob")
        expect_s3_class(protNormRuvDesignResult(out$value)$separation_plot, "grob")
        factors <- kit129Quiet(computeProtRuvModelFactors(out$normalised, ruv_k = separated$k, ruv_k_rule = "sweep_bounded",
            ruv_grouping_variable = "group", controls = ids %in% listed, ruv_design = "supplied_controls"))
        expected <- kit129ProtKernelDa(out$normalised, factors, "group", contrast)
        expect_identical(out$table$Protein.Group, expected$Protein.Group)
        expect_equal(out$table[-1L], expected[-1L], tolerance = 1e-12)
    }
})

test_that("(fifth pass item 4) lipidomics and metabolomics: the separation rule on the run path, per assay, with its plots", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    controls <- list(pos = sprintf("pos_f%03d", 11:40), neg = sprintf("neg_f%03d", 11:40))
    for (omic in c("lipidomics", "metabolomics")) {
        long <- if (identical(omic, "lipidomics")) "da_lipids_long" else "da_metabolites_long"
        for (design in c("first_pass_nulls", "supplied_controls")) {
            parameters <- utils::modifyList(kit129OmicsRuvParameters(design, ruv_k = NULL, ruv_k_rule = "sweep_bounded"),
                list(ruv_control_ids = if (identical(design, "supplied_controls")) unlist(controls, use.names = FALSE)))
            run <- kit129OmicsRunPath(omic, parameters)
            value <- run$ruv$value
            separated <- kit129Quiet(omicsRuvSeparationRuleK(run$fixture$object, design,
                if (identical(design, "supplied_controls")) controls, "group"))
            expect_identical(value$ruv_k_rule, "sweep_bounded", info = paste(omic, design))
            kernel_k <- lapply(separated$k, function(k) k)
            direct <- kit129Quiet(runOmicsRuvFactors(run$fixture$object, design, kernel_k, "sweep_bounded",
                ruv_control_ids = if (identical(design, "supplied_controls")) controls))
            expect_identical(value$best_k, lapply(direct$factors, `[[`, "k"), info = paste(omic, design))
            kernel <- kit129OmicsKernelDa(omic, run$fixture, direct$factors)
            expect_equal(kit129OmicsKey(omic, run$da$results[[long]]), kit129OmicsKey(omic, kernel[[long]]), tolerance = 1e-12,
                info = paste(omic, design))
            # the plot per assay, for the panel's assay rows
            shown <- omicsRuvDesignResult(value)
            for (assay in names(shown$assays)) expect_s3_class(shown$assays[[assay]]$separation_plot, "grob")
        }
    }
})

# ---- Sixth pass items 4 to 6: the record's fields, the study parameters, save -> reopen -> replay ------------------------

# A replay project: a copy of the store under `source_dir` cut to `preconditions`, as the omics replay tests make it.
kit129ReplayProject <- function(source_dir, preconditions, env = parent.frame()) {
    harness <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    copy <- withr::local_tempdir(.local_envir = env)
    stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
    store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
    copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
    list(prepared = harness$kitPrepareReplayProject(copy_root, preconditions, api), origin = store_root, harness = harness, api = api)
}

# The record's input files, by the basenames the record names, under `source_dir`.
kit129ReplayInputs <- function(record, source_dir) {
    files <- list.files(source_dir, recursive = TRUE, full.names = TRUE)
    inputs <- Filter(function(input) is.character(input$basename) && any(basename(files) == input$basename), record$inputs)
    stats::setNames(vapply(inputs, function(input) files[basename(files) == input$basename][[1L]], ""),
        vapply(inputs, function(input) input$role, ""))
}

# Save the record as JSON, reopen it, replay it on a copy cut to `preconditions`: the states and their fingerprints.
kit129SaveReopenReplay <- function(owner, source_dir, preconditions, paths = list(), env = parent.frame()) {
    json <- file.path(withr::local_tempdir(.local_envir = env), "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    project <- kit129ReplayProject(source_dir, preconditions, env = env)
    replay <- kit129Quiet(project$harness$kitReplayJourney(record, project$prepared, project$origin,
        kit129ReplayInputs(record, source_dir), paths, withr::local_tempdir(.local_envir = env), project$api))
    list(record = record, replay = replay, project = project)
}

# A record's parameter of a step, decoded.
kit129StepParameter <- function(step, name) {
    hit <- Filter(function(parameter) identical(parameter$name, name), step$parameters)
    if (length(hit)) decodeWorkflowRecordValue(hit[[1L]]$value)
}

# The same record with one recorded parameter of the named operation's step replaced.
kit129TamperParameter <- function(record, operation, name, value) {
    for (i in seq_along(record$steps)) {
        if (!identical(record$steps[[i]]$operation, operation)) next
        for (j in seq_along(record$steps[[i]]$parameters)) {
            if (identical(record$steps[[i]]$parameters[[j]]$name, name)) record$steps[[i]]$parameters[[j]]$value <- encodeWorkflowRecordValue(value)
        }
    }
    record
}

test_that("(sixth pass item 4) the record fields are one set; the methods follow from the DA method", {
    fields <- ruvDesignRecordFields(ruv_design = "groups_all_features", ruv_k = 2L, ruv_da_method = "ruv_factors_exact_se_trend_robust",
        ruv_k_typed = NULL, ruv_control_list_count = integer())
    expect_identical(names(fields), RUV_DESIGN_RECORD_FIELDS)
    expect_identical(fields[c("ruv_design", "ruv_k", "ruv_se_method", "ruv_test_method", "ruv_k_typed", "ruv_control_list_count")],
        list(ruv_design = "groups_all_features", ruv_k = 2L, ruv_se_method = "exact_per_comparison", ruv_test_method = "limma_treat_trend_robust",
            ruv_k_typed = NA, ruv_control_list_count = NA))
    expect_identical(ruvDesignDaMethods("limma_ruv_factors_exact_se"), c(se = "exact_per_comparison", test = "limma_da_settings"))
    expect_identical(ruvDesignDaMethods(NA_character_), c(se = NA_character_, test = NA_character_))
    expect_error(ruvDesignRecordFields(not_a_field = 1))
    # every code has words
    texts <- omicGuiRuvRecordTexts()
    for (method in c("ruv_factors_exact_se_trend_robust", "ruv_none_trend_robust", "dpcde_ruv_factors_exact_se", "dpcde_ruv_none",
        "limma_ruv_factors_exact_se", "edger_qlf_factors_in_model", "edger_qlf_none")) {
        codes <- ruvDesignDaMethods(method)
        expect_true(all(paste0(c("se_", "test_"), codes) %in% names(texts)), info = method)
    }
    expect_false(any(grepl("=", texts, fixed = TRUE)))
})

test_that("(sixth pass items 4 to 6) classic proteomics: the record and the study parameters carry the design; save, reopen and replay give the same states", {
    skip_if_not_installed("sva")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    ids <- as.character(object@protein_quant_table[[object@protein_id_column]])
    listed <- c(ids[seq(1L, length(ids), by = 3L)], "NOT_A_PROTEIN")
    run_on <- function(parameters) {
        root <- withr::local_tempdir(.local_envir = parent.frame())
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit129-replay", "DIA"), WorkflowState$new(audit_enabled = FALSE))
        manager$.context$reference_states <- TRUE
        manager$setWorkflowType("DIA")
        setDIAStepProvenance(manager, "design_builder")
        manager$saveState("protein_replicate_filtered", object, list(), "fixture")
        run <- kit129Quiet(withDIAStateScope(manager, function() {
            run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
            commit <- diaNormalizationCommit(manager, parameters, run)
            commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            run
        }))
        list(manager = manager, root = root, run = run)
    }
    declared <- c(kit129ProtParameters(), list(ruv_design = "supplied_controls", ruv_k_rule = "user", ruv_k_seed = 11L,
        ruv_control_ids = rev(listed)))
    origin <- run_on(declared)
    withr::defer(closeDIAStateAdapter(origin$manager))
    record <- buildWorkflowAnalysisRecord(origin$manager)
    step <- Filter(function(step) identical(step$operation, "normalization-ruv"), record$steps)[[1L]]
    value <- origin$run$value
    expect_identical(step$outcomes[c("ruv_design", "ruv_replicate_structure", "ruv_grouping_variable", "ruv_k_rule", "ruv_control_rule",
        "ruv_control_list_sha256", "ruv_da_method", "ruv_se_method", "ruv_test_method", "ruv_first_pass_rule")],
        list(ruv_design = "supplied_controls", ruv_replicate_structure = value$ruv_model$replicate_structure, ruv_grouping_variable = "group",
            ruv_k_rule = "user", ruv_control_rule = "supplied_list", ruv_control_list_sha256 = trxControlListSha256(listed),
            ruv_da_method = "ruv_factors_exact_se_trend_robust", ruv_se_method = "exact_per_comparison",
            ruv_test_method = "limma_treat_trend_robust", ruv_first_pass_rule = NA))
    expect_equal(unlist(step$outcomes[c("ruv_k_typed", "ruv_k", "ruv_k_seed", "ruv_k_permutations", "ruv_control_list_count",
        "ruv_control_list_found", "ruv_control_count")]), c(ruv_k_typed = 1, ruv_k = 1, ruv_k_seed = 11, ruv_k_permutations = 200,
        ruv_control_list_count = length(listed), ruv_control_list_found = length(listed) - 1L, ruv_control_count = value$ruv_model$control_count))
    expect_identical(step$outcomes$ruv_k_estimate, value$ruv_result$ruv_k_estimate)
    expect_true(step$outcomes$ruv_group_association_outcome %in% c("detectable", "not_detectable", "no_limit"))
    # the list itself is kept in the step's record (its identity holds its digest and count only)
    expect_identical(kit129StepParameter(step, "ruv_control_ids"), rev(listed))
    expect_identical(kit129StepParameter(step, "ruv_control_list_sha256"), trxControlListSha256(listed))
    # the study parameters' section, in words, once
    lines <- workflowRuvDesignStudyLines(record)
    expect_identical(lines[[1L]], "Unwanted variation (RUV):")
    expect_identical(sum(lines == "Unwanted variation (RUV):"), 1L)
    expect_false(any(grepl("=", lines, fixed = TRUE)))
    texts <- omicGuiRuvRecordTexts()
    expect_true(any(lines == paste0("* ", texts[["approach"]], ": ", texts[["design_supplied_controls"]])))
    expect_true(any(lines == paste0("* ", texts[["list"]], ": ", sprintf(texts[["list_value"]], length(listed), length(listed) - 1L,
        trxControlListSha256(listed)))))
    expect_true(any(startsWith(lines, paste0("* ", texts[["methods"]], ": Unwanted variation was estimated with RUV-III"))))
    # the record's workbook has the same rows
    expect_identical(workflowAnalysisRecordTables(record)$ruv_design, workflowRuvDesignSummary(record)$rows)
    # save -> reopen -> replay: every state equal, the list from the record
    saved <- kit129SaveReopenReplay(origin$manager, origin$root, c("initial", "protein_replicate_filtered"))
    expect_identical(saved$replay$states, unlist(origin$manager$.context$snapshot$state_history))
    expect_identical(as.character(saved$replay$mismatches), character())
    # a record whose list does not hash to its recorded digest is refused before any step runs
    tampered <- kit129TamperParameter(saved$record, "normalization-ruv", "ruv_control_ids", listed[-1L])
    project <- kit129ReplayProject(origin$root, c("initial", "protein_replicate_filtered"))
    expect_error(kit129Quiet(replayWorkflowAnalysisRecord(tampered, project$prepared, character(), list(), withr::local_tempdir())),
        class = "workflow_replay_invalid")
    # today's design: no field, no section, the record's parameters are today's declaration's, and it replays as before
    legacy <- run_on(kit129ProtParameters())
    withr::defer(closeDIAStateAdapter(legacy$manager))
    legacy_record <- buildWorkflowAnalysisRecord(legacy$manager)
    legacy_step <- Filter(function(step) identical(step$operation, "normalization-ruv"), legacy_record$steps)[[1L]]
    expect_false(any(RUV_DESIGN_RECORD_FIELDS %in% names(legacy_step$outcomes)))
    recorded_names <- vapply(legacy_step$parameters, function(p) p$name, "")
    expect_true(all(diaNormalizationOperation()$parameters %in% recorded_names))
    expect_false(any(setdiff(RUV_DESIGN_REQUEST_PARAMETERS, "ruv_k") %in% recorded_names))
    expect_null(workflowRuvDesignSummary(legacy_record))
    expect_identical(workflowRuvDesignStudyLines(legacy_record), character())
    expect_null(workflowAnalysisRecordTables(legacy_record)$ruv_design)
    replayed <- kit129SaveReopenReplay(legacy$manager, legacy$root, c("initial", "protein_replicate_filtered"))
    expect_identical(as.character(replayed$replay$mismatches), character())
})

test_that("(sixth pass items 4 and 5) lipidomics and metabolomics: the record carries the design per assay; save, reopen and replay give the same states", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    for (omic in c("lipidomics", "metabolomics")) {
        fixture <- kit129OmicsOwner(omic)
        binding <- omicsDeclarationBindings()[[omic]]
        # the features of kit129OmicsObject's two assays (the fifth pass's lists)
        listed <- c(sprintf("neg_f%03d", 21:40), sprintf("pos_f%03d", 11:30), "NOT_A_FEATURE")
        parameters <- utils::modifyList(kit129OmicsRuvParameters("supplied_controls", ruv_k = 1L, ruv_k_rule = "user"),
            list(ruv_control_ids = rev(listed)))
        # the module's configuration, as the tab passes it (the replay commits with the head's)
        kit129Quiet(runWorkflowOmicsReferenceRuv(binding, fixture$owner, parameters, fixture$workflow$config_list))
        record <- buildWorkflowAnalysisRecord(fixture$owner)
        step <- Filter(function(step) identical(step$operation, paste0(binding$prefix, "_ruv")), record$steps)[[1L]]
        o <- step$outcomes
        expect_identical(o[c("ruv_design", "ruv_k_rule", "ruv_control_rule", "ruv_control_list_sha256", "ruv_da_method", "ruv_se_method")],
            list(ruv_design = "supplied_controls", ruv_k_rule = "user", ruv_control_rule = "supplied_list",
                ruv_control_list_sha256 = trxControlListSha256(listed), ruv_da_method = "limma_ruv_factors_exact_se",
                ruv_se_method = "exact_per_comparison"), info = omic)
        expect_equal(unlist(o[c("ruv_k_typed", "ruv_control_list_count", "ruv_control_list_found")]),
            c(ruv_k_typed = 1, ruv_control_list_count = length(listed), ruv_control_list_found = length(listed) - 1L), info = omic)
        expect_setequal(names(o$ruv_k_per_assay), c("pos", "neg"))
        expect_setequal(names(o$ruv_control_count_per_assay), c("pos", "neg"))
        expect_identical(kit129StepParameter(step, "ruv_control_ids"), rev(listed), info = omic)
        lines <- workflowRuvDesignStudyLines(record)
        expect_identical(lines[[1L]], "Unwanted variation (RUV):", info = omic)
        expect_false(any(grepl("=", lines, fixed = TRUE)), info = omic)
        preconditions <- unlist(fixture$owner$.context$snapshot$state_history)[1:2]
        saved <- kit129SaveReopenReplay(fixture$owner, fixture$source_dir, preconditions, fixture$paths)
        expect_identical(saved$replay$states, unlist(fixture$owner$.context$snapshot$state_history), info = omic)
        expect_identical(as.character(saved$replay$mismatches), character(), info = omic)
    }
})

test_that("(sixth pass items 4 and 5) transcriptomics: the record carries the design; save, reopen and replay give the same states", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    fixture <- kit129TrxOwner()
    owner <- fixture$owner
    genes <- rownames(trxAssayInputs(readTrxWorkflowState(owner, getWorkflowStateReference(owner)))$counts)
    listed <- c(genes[seq(1L, length(genes), by = 7L)], "NOT_A_GENE")
    kit129Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("supplied_controls", 2L, "user",
        ruv_control_ids = rev(listed))))
    record <- buildWorkflowAnalysisRecord(owner)
    step <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)[[1L]]
    o <- step$outcomes
    expect_identical(o[c("ruv_design", "ruv_k_rule", "ruv_control_rule", "ruv_control_list_sha256", "ruv_da_method", "ruv_se_method",
        "ruv_test_method")], list(ruv_design = "supplied_controls", ruv_k_rule = "user", ruv_control_rule = "supplied_list",
        ruv_control_list_sha256 = trxControlListSha256(listed), ruv_da_method = "edger_qlf_factors_in_model", ruv_se_method = "edger_contrast",
        ruv_test_method = "edger_qlf"))
    expect_equal(unlist(o[c("ruv_k_typed", "ruv_k", "ruv_control_list_count", "ruv_control_list_found")]),
        c(ruv_k_typed = 2, ruv_k = 2, ruv_control_list_count = length(listed), ruv_control_list_found = length(listed) - 1L))
    lines <- workflowRuvDesignStudyLines(record)
    expect_true(any(startsWith(lines, paste0("* ", omicGuiRuvRecordTexts()[["methods"]], ": Unwanted variation was estimated with RUVg"))))
    saved <- kit129SaveReopenReplay(owner, fixture$source_dir, c("initial", "trx_raw_data_s4"), fixture$paths)
    expect_identical(saved$replay$states, unlist(owner$.context$snapshot$state_history))
    expect_identical(as.character(saved$replay$mismatches), character())
})

test_that("(sixth pass items 4 and 5) DIA with limpa: the record carries the declared rules; save, reopen and replay give the same states", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    skip_if_not_installed("fastlimpadpc")
    sim <- kit129LimpaSimulate(129L, proteins = 120L, per_group = 5L)
    quantified <- kit129Quiet(kit129LimpaQuantify(sim, "fastlimpadpc"))
    ids <- as.character(quantified@protein_quant_table[[quantified@protein_id_column]])
    listed <- ids[seq(2L, length(ids), by = 3L)]
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit129-dia-replay", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", quantified, list(), "fixture")
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "supplied_list", ruv_control_ids = rev(listed),
        ruv_k_rule = "user", ruv_k = 1L, ruv_k_seed = 5L)
    run <- kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceRuvCovariates(manager, getDIAStateReference(manager), parameters)
        commit <- diaRuvCovariatesCommit(manager, parameters, run)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        run
    }))
    record <- buildWorkflowAnalysisRecord(manager)
    step <- Filter(function(step) identical(step$operation, "protein_ruv_covariates"), record$steps)[[1L]]
    o <- step$outcomes
    expect_identical(o[c("ruv_design", "ruv_k_rule", "ruv_control_rule", "ruv_control_list_sha256", "ruv_da_method", "ruv_test_method")],
        list(ruv_design = "supplied_controls", ruv_k_rule = "user", ruv_control_rule = "supplied_list",
            ruv_control_list_sha256 = trxControlListSha256(listed), ruv_da_method = "dpcde_ruv_factors_exact_se", ruv_test_method = "limpa_dpcde"))
    expect_equal(unlist(o[c("ruv_k_typed", "ruv_k", "ruv_k_seed", "ruv_control_list_count", "ruv_control_list_found", "ruv_control_count")]),
        c(ruv_k_typed = 1, ruv_k = 1, ruv_k_seed = 5, ruv_control_list_count = length(listed), ruv_control_list_found = length(listed),
            ruv_control_count = length(listed)))
    expect_true(isTRUE(o$ruv_applied))
    # the list as the declaration holds it: its canonical form (each id once, sorted)
    expect_identical(kit129StepParameter(step, "ruv_control_ids"), trxCanonicalControlIds(listed))
    saved <- kit129SaveReopenReplay(manager, root, c("initial", "correlation_filtered"))
    expect_identical(saved$replay$states, unlist(manager$.context$snapshot$state_history))
    expect_identical(as.character(saved$replay$mismatches), character())
})

# Seventh pass items 2 and 3: the proteomics integration input (KIT-061) of a new design's run, committed as the tab commits
# it (`diaNormalizationCommit`) and written by the Summary's writer (`protSummaryIntegrationWrite`), read back by the contract.
kit129IntegrationInput <- function(parameters, env = parent.frame()) {
    manager <- kit129ProtStore(env)
    manager$setWorkflowType("DIA")
    root <- withr::local_tempdir(.local_envir = env)
    path <- file.path(root, "integration", "proteomics", INTEGRATION_INPUT_MANIFEST)
    out <- kit129Quiet(withDIAStateScope(manager, function() {
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), parameters)
        commit <- diaNormalizationCommit(manager, parameters, run)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        list(run = run, commit = commit, normalised = readDIAExecutionComponents(manager, run$value$normalized))
    }))
    out$written <- protSummaryIntegrationWrite(list2env(list(state_manager = manager)), path)
    out$input <- readIntegrationInput(dirname(path))
    out$path <- path
    out
}

test_that("(seventh pass item 2) classic proteomics: a new design whose count gave no factor is written to the integration input as not corrected", {
    testthat::local_mocked_bindings(computeProtRuvModelFactors = function(theObject, ...) list(status = "none", k = 0L,
        reason = "the count found no unwanted factor (k = 0)", ruv_k_rule = "num_sv_be", ruv_k_estimate = 0L, k_estimate = NULL,
        W = NULL, controls = NULL, ruv_controls = "all"), .package = "MultiScholaR")
    automatic <- kit129ProtParameters("automatic")
    automatic[c("auto_percentage_min", "auto_percentage_max", "separation_metric", "k_penalty_weight", "adaptive_k_penalty")] <-
        list(5, 10, "max_difference", 0.5, TRUE)
    out <- kit129IntegrationInput(c(automatic, list(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = NULL)))
    expect_false(out$run$value$ruv_result$ruv_applied)
    expect_identical(out$written, "Saved Integration input")
    # the state keeps its name and its matrix; its record says that no RUV was applied
    expect_identical(out$commit$state_name, "ruv_corrected")
    expect_false(out$commit$config$ruv_applied)
    manifest <- out$input$manifest
    expect_identical(manifest$matrix_role, "normalised")
    expect_null(manifest$correction)
    expect_identical(manifest$source_state$name, "ruv_corrected")
    expected <- out$normalised@protein_quant_table
    names(expected)[[1L]] <- "feature_id"
    expect_equal(out$input$views$proteomics$matrix, as.data.frame(expected), tolerance = 0, ignore_attr = TRUE)
})
