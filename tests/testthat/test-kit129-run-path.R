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
kit129OmicsObject <- function(omic = "lipidomics", seed = 129L, p = 50L) {
    set.seed(seed)
    design <- data.frame(Sample_ID = sprintf("S%02d", 1:12), group = rep(c("A", "B"), each = 6L), stringsAsFactors = FALSE)
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
    if (identical(omic, "lipidomics")) {
        return(methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = "database_identifier",
            annotation_id_column = "name", design_matrix = design, sample_id = "Sample_ID", group_id = "group", args = list()))
    }
    methods::new("MetaboliteAssayData", metabolite_data = tables, metabolite_id_column = "database_identifier",
        annotation_id_column = "name", design_matrix = design, sample_id = "Sample_ID", group_id = "group", args = list())
}

kit129OmicsOwner <- function(omic, env = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = env)
    object <- kit129OmicsObject(omic)
    fixture <- if (identical(omic, "lipidomics")) lipidDeclarationsOwner(root, object = object) else metabDeclarationsOwner(root, object = object)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    fixture$object <- object
    fixture$contrasts <- data.frame(contrasts = "groupB-groupA", friendly_names = "B vs A", full_format = "B_vs_A=groupB-groupA",
        stringsAsFactors = FALSE)
    fixture
}

kit129OmicsRuvParameters <- function(ruv_design = "groups_all_features", ruv_k = 1L, ruv_k_rule = "user") {
    list(ruv_design = ruv_design, ruv_grouping_variable = "group", ruv_k = ruv_k, ruv_k_rule = ruv_k_rule, ruv_k_seed = NULL,
        ruv_control_ids = NULL)
}

kit129ProtStore <- function(env = parent.frame()) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = env), "kit-129-run-path", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_replicate_filtered", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), list(), "source")
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
