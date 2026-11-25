# KIT-200 (VAL-116 (4)): a protein with no value in any analysed sample has no mean (limma's Amean NaN), so the trended
# empirical-Bayes step of the proteomics DA stops (limma 3.68.4: "prior.weights contain NA values" with unequal residual
# df, "covariate contains NA or infinite values" with equal df). Such proteins are removed before eBayes and kept in every
# contrast table as rows without statistics, as KIT-196 does for lipidomics and metabolomics. Each test first asserts
# that its call completed (`expect_null` on the caught message), so the tests-only commit fails by assertion at the
# stated limma message. Synthetic data only; builders are copied from the files the ticket names, never called.

kit200Strings <- c("groupB-groupA", "groupC-groupA", "groupC-groupB")
kit200Columns <- c("logFC", "AveExpr", "t", "raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment")
kit200Scatter <- list(c("P2", "S1"), c("P3", "S5"), c("P3", "S9"), c("P4", "S12"))
kit200SkipParameters <- list(norm_method = "quantile", ruv_mode = "skip", auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = NULL, ruv_k = NULL, grouping = "group")

kit200Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

kit200WithSeed <- function(seed, code) withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# The fixture (section 4): proteins P1..P60, levels N(22, 2^2), residual sd 0.3, P11..P16 with a group-B effect of +1.
kit200Matrix <- function(sizes = c(A = 4L, B = 4L, C = 4L), seed = 200001L, proteins = 60L) {
    groups <- rep(names(sizes), sizes)
    x <- kit200WithSeed(seed, {
        level <- stats::rnorm(proteins, 22, 2)
        effect <- as.numeric(seq_len(proteins) %in% 11:16)
        level + outer(effect, as.numeric(groups == "B")) +
            matrix(stats::rnorm(proteins * length(groups), sd = 0.3), proteins)
    })
    dimnames(x) <- list(paste0("P", seq_len(proteins)), paste0("S", seq_along(groups)))
    x
}

kit200Design <- function(sizes = c(A = 4L, B = 4L, C = 4L)) {
    runs <- paste0("S", seq_len(sum(sizes)))
    data.frame(Run = runs, group = rep(names(sizes), sizes), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

kit200SetScatter <- function(x, cells = kit200Scatter) {
    for (cell in cells) x[cell[[1L]], cell[[2L]]] <- NA_real_
    x
}

kit200AddRow <- function(x, id, values = numeric()) {
    row <- matrix(NA_real_, 1L, ncol(x), dimnames = list(id, colnames(x)))
    row[1L, names(values)] <- values
    rbind(x, row)
}

# The chain fixtures' samples that the correlation filter must remove: each column's values permuted across the
# proteins, so that the sample's within-group Pearson correlations are near 0, far below the threshold (0.5).
kit200Decorrelate <- function(x, samples, seed = 200002L) {
    kit200WithSeed(seed, for (sample in samples) x[, sample] <- x[sample.int(nrow(x)), sample])
    x
}

kit200Run <- function(x, design, lfc = NA, trend = TRUE, robust = TRUE, defaults = FALSE) {
    tryCatch(list(value = kit200Quiet(if (defaults) runTestsContrasts(x, kit200Strings, design, "~ 0 + group") else
        runTestsContrasts(x, kit200Strings, design, "~ 0 + group", treat_lfc_cutoff = lfc, eBayes_trend = trend,
            eBayes_robust = robust))), error = function(e) list(error = conditionMessage(e)))
}

# The base-side guard: the call completed. On the tests-only commit this fails and prints limma's message.
kit200Completed <- function(run, label) {
    expect_null(run$error, label = label)
    is.null(run$error)
}

kit200ExpectUntested <- function(tables, id) {
    for (contrast in names(tables)) {
        table <- tables[[contrast]]
        expect_identical(sum(rownames(table) == id), 1L, info = contrast)
        present <- intersect(c(kit200Columns, "B"), names(table))
        expect_true(all(is.na(unlist(table[id, present]))), info = contrast)
        expect_identical(table[id, "AveExpr"], NA_real_, info = contrast)
        expect_identical(sum(is.na(table$AveExpr)), 1L, info = contrast)
    }
}

# "Equal" (section 9): identical, or with a complete reference to rounding at 1e-12 (outcome 3), the largest relative
# difference printed. AveExpr is a row mean, identical on both paths.
kit200ExpectEqual <- function(got, ref, ids, tolerance = NULL) {
    for (contrast in names(ref)) {
        expect_true(all(ids %in% rownames(got[[contrast]])), info = contrast)
        for (column in kit200Columns) {
            a <- got[[contrast]][ids, column]
            b <- ref[[contrast]][ids, column]
            if (is.null(tolerance) || identical(column, "AveExpr")) {
                expect_identical(a, b, info = paste(contrast, column))
            } else {
                largest <- suppressWarnings(max(abs(a - b) / pmax(abs(b), .Machine$double.xmin), na.rm = TRUE))
                expect_equal(a, b, tolerance = tolerance,
                    info = sprintf("%s %s: largest relative difference %.3g", contrast, column, largest))
            }
        }
    }
}

# ---- In-app chains (T6, T7, T10): the declared protein QC steps, the correlation filter, the DA tab ------------------

kit200Contrasts <- function() data.frame(contrasts = c("B.vs.A", "C.vs.A"),
    full_format = c("B.vs.A=groupB-groupA", "C.vs.A=groupC-groupA"), friendly_names = c("B.vs.A", "C.vs.A"),
    stringsAsFactors = FALSE)

kit200Object <- function(x, design) {
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = list())
}

# A state read back from one store, made savable in another: the record of the store it came from is dropped.
kit200Detach <- function(object) {
    attr(object, "dia_component_source") <- NULL
    object
}

# An owner with the paths the steps and the DA publication write to (savedDaW1Owner's shape, copied), holding `object`
# as `state`. The owner is closed with the calling frame.
kit200Owner <- function(object, type = "DIA", state = "protein_s4_initial", env = parent.frame()) {
    root <- normalizePath(withr::local_tempdir(.local_envir = env))
    project <- paste0("kit200-", tolower(type))
    paths <- list(source_dir = file.path(root, "source"), protein_qc_dir = file.path(root, "protein_qc"),
        time_dir = file.path(root, "time"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, project, type), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(type)
    manager$.context$experiment_paths <- paths
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState(state, object, list(), "fixture")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list(),
        tab_status = list()))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = env)
    list(root = root, project = project, paths = paths, workflow = workflow, da = new.env(parent = emptyenv()))
}

kit200Head <- function(fixture) {
    manager <- fixture$workflow$state_manager
    withDIAStateScope(manager, function() readDIAExecutionComponents(manager, getDIAStateReference(manager)))
}

kit200Ids <- function(object) as.character(object@protein_quant_table$Protein.Ids)

# The protein QC tab's steps in process (their apply-step functions), then normalisation with RUV skipped (T6 (d)) and
# the correlation filter's declarations (computeDIAReferenceCorrelation, as the normalisation tab commits it).
kit200Qc <- function(fixture, m, g, default_chain = FALSE, normalise = FALSE, threshold = 0.5) {
    workflow <- fixture$workflow
    manager <- workflow$state_manager
    out <- list()
    kit200Quiet(runProteinIntensityFilterApplyStep(workflow, FALSE, m, g, 1))
    out$intensity <- kit200Head(fixture)
    if (default_chain) {
        kit200Quiet(runProteinDuplicateRemovalStep(workflow, "mean"))
        kit200Quiet(runProteinReplicateFilterApplyStep(workflow, fixture$paths, "group", 1L))
    }
    out$qc <- kit200Head(fixture)
    kit200Quiet(withDIAStateScope(manager, function() {
        source <- getDIAStateReference(manager)
        if (normalise) {
            run <- computeDIAReferenceNormalization(manager, source, kit200SkipParameters)
            commit <- diaNormalizationCommit(manager, kit200SkipParameters, run)
            commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
            source <- run$value$final
        }
        correlation <- computeDIAReferenceCorrelation(manager, source, "group", threshold)
        commit <- diaCorrelationCommit(manager, correlation$final, threshold, FALSE)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
    }))
    out$input <- kit200Head(fixture)
    out
}

# ---- The managed DA click (savedDaW2Dispatch, savedDaW2Front, savedDaW2Managed, savedDaW2Dispatcher, savedDaW2Mocks of
# test-workflow-saved-da.R, copied) ----------------------------------------------------------------------------------

kit200Dispatch <- function(fixture, controller, operation, parameters, input) {
    store <- fixture$workflow$state_manager$.context$store
    id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
    paths <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation,
        expected_head = workflowExecutionHead(store)$entity_id, generation = controller$control$generation,
        controller_id = controller$control$controller_id)
    writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"),
        list(experiment_label = "core", experiment_paths = fixture$paths, omic_type = "proteomics"))
    writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id,
        status = "admitting", expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
        input_references = list(state = input), parameters = parameters, token = token, owner = store$identity,
        project_root = dirname(store$directory), directory = paths$request_dir, config = list(mode = "adapters"),
        worker_cap = 12 * 1024^3)
    path <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
        workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
    message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
    controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
        parameters = parameters), paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, NULL))
    receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
    accepted <- resumeWorkflowAction(controller, id)
    c(readManagedOmicsResult(fixture$workflow$state_manager, accepted, NULL),
        list(acceptance = receipt$status, accepted = accepted))
}

kit200Dispatcher <- function(fixture) {
    owner <- NULL
    controller <- NULL
    function(operation, parameters, input) {
        manager <- fixture$workflow$state_manager
        if (!identical(manager, owner)) {
            owner <<- manager
            controller <<- omicsWorkerController(list(owner = manager))
        }
        kit200Dispatch(fixture, controller, operation, parameters, input)
    }
}

kit200Managed <- function(fixture, dispatch, settings, contrasts) {
    manager <- fixture$workflow$state_manager
    annotations <- fixture$workflow$uniprot_dat_cln
    request <- withDIAStateScope(manager, function() diaDaRequest(manager, settings, fixture$da, fixture$workflow,
        contrasts, annotations))
    route <- withDIAStateScope(manager, function() diaActiveRoute(manager))
    results <- list()
    for (step in request$steps) {
        delivered <- dispatch("DA-contrast-output", c(step, list(route = route)), request$source)
        results[[step$name]] <- delivered$value
        if (!is.null(delivered$display$base_fit)) manager$.context$da_base_fit_reference <- delivered$display$base_fit
    }
    manifest <- withDIAStateScope(manager, function() acceptDIADaManifest(manager, request$manifest, results))
    published <- dispatch("da_publication", list(annotations = saveDIADisplayTable(manager,
        fixture$workflow$uniprot_dat_cln)), getDIAStateReference(manager))
    delivered <- withDIAStateScope(manager, function() readManagedDaResult(manager, published$accepted, NULL))
    completeDIADaDisplay(fixture$da, fixture$workflow, delivered$manifest)
    manifest
}

kit200ManagedMocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

kit200ManifestTables <- function(manager, manifest) {
    withDIAStateScope(manager, function() lapply(stats::setNames(nm = names(manifest$results)), function(name)
        lapply(c(results = "contrasts_results_table", long = "da_proteins_long"), function(field)
            as.data.frame(materializeDIATable(manager$.context, manifest$results[[name]]$tables[[field]])))))
}

# The DA tab's click on the chain's head (in process, `runDIAReferenceDaAnalysis`; or managed), with the rows of every
# public fit.eb and of every base fit recorded as the kernel returns them.
kit200Da <- function(fixture, contrasts = kit200Contrasts(), publish = TRUE, managed = FALSE) {
    manager <- fixture$workflow$state_manager
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    fixture$da$current_s4_object <- getDIAStateReference(manager)
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    seen <- new.env(parent = emptyenv())
    seen$fit <- list()
    seen$base <- list()
    run_tests <- get("runTestsContrasts", envir = asNamespace("MultiScholaR"))
    base_fit <- get("fitDIAExecutionBaseModel", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(runTestsContrasts = function(...) {
        result <- run_tests(...)
        seen$fit[[length(seen$fit) + 1L]] <- rownames(result$fit.eb$coefficients)
        result
    }, fitDIAExecutionBaseModel = function(...) {
        artifact <- base_fit(...)
        seen$base[[length(seen$base) + 1L]] <- rownames(artifact$fit$coefficients)
        artifact
    }, .package = "MultiScholaR")
    run <- tryCatch(list(value = kit200Quiet(if (managed)
        kit200Managed(fixture, kit200Dispatcher(fixture), settings, contrasts) else
        runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts, NULL,
            publish = publish))), error = function(e) list(error = conditionMessage(e)))
    run$fit <- seen$fit
    run$base <- seen$base
    if (is.null(run$error)) run$tables <- kit200ManifestTables(manager, run$value)
    run
}

# One chain: the owner, the QC steps, the DA. `x` holds the fixture's proteins and the special protein.
kit200Chain <- function(x, sizes, m, g, type = "DIA", default_chain = FALSE, normalise = FALSE, managed = FALSE,
    publish = TRUE, env = parent.frame()) {
    fixture <- kit200Owner(kit200Object(x, kit200Design(sizes)), type, env = env)
    qc <- kit200Qc(fixture, m, g, default_chain, normalise)
    list(fixture = fixture, qc = qc, da = kit200Da(fixture, publish = publish, managed = managed))
}

kit200InputRow <- function(object, id) {
    table <- object@protein_quant_table
    unlist(table[table$Protein.Ids == id, setdiff(names(table), "Protein.Ids")])
}

# T6's head assertions on one chain's tables: one untested row per contrast table, with its statistics and AveExpr NA
# in the results table and its fold change, p, q and every log2norm column NA in the long table.
kit200ExpectChainUntested <- function(tables, id) {
    for (contrast in names(tables)) {
        results <- tables[[contrast]]$results
        row <- results[results$Protein.Ids == id, , drop = FALSE]
        expect_identical(nrow(row), 1L, info = contrast)
        expect_true(all(is.na(unlist(row[intersect(c(kit200Columns, "B"), names(row))]))), info = contrast)
        expect_identical(row$AveExpr, NA_real_, info = contrast)
        expect_identical(sum(is.na(results$AveExpr)), 1L, info = contrast)
        long <- tables[[contrast]]$long
        row <- long[long$Protein.Ids == id, , drop = FALSE]
        expect_identical(nrow(row), 1L, info = contrast)
        norm <- grep("^log2norm\\.", names(row), value = TRUE)
        expect_gt(length(norm), 0L)
        expect_true(all(is.na(unlist(row[c("log2FC", "raw_pvalue", "fdr_qvalue", norm)]))), info = contrast)
    }
}

kit200ExpectChainEqual <- function(tables, reference, id) {
    for (contrast in names(reference)) {
        got <- tables[[contrast]]$results
        ref <- reference[[contrast]]$results
        ids <- ref$Protein.Ids
        expect_false(id %in% ids)
        rownames(got) <- got$Protein.Ids
        rownames(ref) <- ref$Protein.Ids
        expect_true(all(ids %in% rownames(got)), info = contrast)
        for (column in intersect(kit200Columns, names(ref))) {
            expect_identical(got[ids, column], ref[ids, column], info = paste(contrast, column))
        }
    }
}

kit200ExpectFits <- function(da, id) {
    expect_gt(length(da$fit), 0L)
    expect_false(any(vapply(da$fit, function(rows) id %in% rows, logical(1))))
    expect_true(any(vapply(da$base, function(rows) id %in% rows, logical(1))))
}

# T6's fixture: three groups of four, the scatter, S1 decorrelated, P_one with its one value (22) in S1.
kit200T6Matrix <- function(with_row = TRUE) {
    x <- kit200SetScatter(kit200Decorrelate(kit200Matrix(), "S1"))
    if (with_row) x <- kit200AddRow(x, "P_one", c(S1 = 22))
    x
}

kit200ExpectT6Preconditions <- function(chain, id, removed) {
    expect_true(id %in% kit200Ids(chain$qc$intensity))
    expect_true(id %in% kit200Ids(chain$qc$qc))
    input <- chain$qc$input
    expect_false(any(removed %in% input@design_matrix$Run))
    expect_false(any(removed %in% names(input@protein_quant_table)))
    expect_identical(sum(kit200Ids(input) == id), 1L)
    expect_true(all(is.na(kit200InputRow(input, id))))
}

# ---- T1 to T5: the kernel ----------------------------------------------------------------------------------------

test_that("KIT-200 T1: equal df, an empty protein is reported untested and the others equal the analysis without it", {
    design <- kit200Design()
    reference_x <- kit200Matrix()
    x <- kit200AddRow(reference_x, "P_empty")
    for (lfc in c(NA, 0)) {
        run <- kit200Run(x, design, lfc = lfc)
        if (!kit200Completed(run, sprintf("T1 error (treat_lfc_cutoff %s)", lfc))) next
        reference <- kit200Run(reference_x, design, lfc = lfc)$value
        expect_identical(run$value$features_without_values, "P_empty")
        kit200ExpectUntested(run$value$results, "P_empty")
        kit200ExpectEqual(run$value$results, reference$results, rownames(reference_x), tolerance = 1e-12)
        expect_false("P_empty" %in% rownames(run$value$fit.eb$coefficients))
        expect_null(reference$features_without_values)
    }
})

test_that("KIT-200 T2: unequal df, an empty protein is reported untested and the others are identical", {
    design <- kit200Design()
    reference_x <- kit200SetScatter(kit200Matrix())
    x <- kit200AddRow(reference_x, "P_empty")
    for (lfc in c(NA, 0)) {
        run <- kit200Run(x, design, lfc = lfc)
        if (!kit200Completed(run, sprintf("T2 error (treat_lfc_cutoff %s)", lfc))) next
        reference <- kit200Run(reference_x, design, lfc = lfc)$value
        expect_identical(run$value$features_without_values, "P_empty")
        kit200ExpectUntested(run$value$results, "P_empty")
        kit200ExpectEqual(run$value$results, reference$results, rownames(reference_x))
        expect_false("P_empty" %in% rownames(run$value$fit.eb$coefficients))
    }
})

test_that("KIT-200 T3 (guard): runTestsContrasts' defaults (trend off) complete and the tested rows are unchanged", {
    design <- kit200Design()
    reference_x <- kit200SetScatter(kit200Matrix())
    run <- kit200Run(kit200AddRow(reference_x, "P_empty"), design, defaults = TRUE)
    if (!kit200Completed(run, "T3 error")) return(invisible())
    reference <- kit200Run(reference_x, design, defaults = TRUE)$value
    kit200ExpectEqual(run$value$results, reference$results, rownames(reference_x))
})

test_that("KIT-200 T4 (guards): a group without values and zero residual df complete with the trend on", {
    design <- kit200Design()
    x <- kit200SetScatter(kit200Matrix())
    # (a) P5 without a value in group C: no C fold change, a finite B-A fold change
    a <- x
    a["P5", paste0("S", 9:12)] <- NA_real_
    run <- kit200Run(a, design, lfc = 0)
    if (kit200Completed(run, "T4 (a) error")) {
        tables <- run$value$results
        expect_true(is.na(tables[["groupC-groupA"]]["P5", "logFC"]))
        expect_true(is.na(tables[["groupC-groupB"]]["P5", "logFC"]))
        expect_true(is.finite(tables[["groupB-groupA"]]["P5", "logFC"]))
    }
    # (b) P6 with one value per group: zero residual df, finite fold changes and p-values
    b <- x
    b["P6", setdiff(colnames(b), c("S1", "S5", "S9"))] <- NA_real_
    run <- kit200Run(b, design, lfc = 0)
    if (kit200Completed(run, "T4 (b) error")) {
        for (contrast in kit200Strings) {
            expect_true(is.finite(run$value$results[[contrast]]["P6", "logFC"]), info = contrast)
            expect_true(is.finite(run$value$results[[contrast]]["P6", "raw_pvalue"]), info = contrast)
        }
    }
})

test_that("KIT-200 T5: input without any value is refused with KIT-196's class and the lipidomics kernel's text", {
    design <- kit200Design()
    x <- kit200Matrix()
    x[] <- NA_real_
    refuse <- function(data) kit200Quiet(runTestsContrasts(data, kit200Strings, design, "~ 0 + group",
        treat_lfc_cutoff = NA, eBayes_trend = TRUE, eBayes_robust = TRUE))
    expect_error(refuse(x), class = "omics_da_no_observed_features")
    unnamed <- x
    rownames(unnamed) <- NULL
    expect_error(refuse(unnamed), class = "omics_da_no_observed_features")
    # D-C: one wording across omics (judged at R2, after KIT-196's merge)
    lipid_fit <- suppressWarnings(limma::lmFit(x, stats::model.matrix(~ 0 + group, design)))
    lipid <- tryCatch(kit200Quiet(applyLipidDaContrasts(lipid_fit, "groupB-groupA")), error = conditionMessage)
    proteomics <- tryCatch(refuse(x), error = conditionMessage)
    expect_identical(proteomics, lipid)
})

# ---- T6 to T10: the routes in the app ----------------------------------------------------------------------------

test_that("KIT-200 T6: in process, replicate filter skipped, the correlation filter empties P_one and the DA completes", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    chain <- kit200Chain(kit200T6Matrix(), c(A = 4L, B = 4L, C = 4L), 1, 1)
    kit200ExpectT6Preconditions(chain, "P_one", "S1")
    if (!kit200Completed(chain$da, "T6 DA error")) return(invisible())
    kit200ExpectChainUntested(chain$da$tables, "P_one")
    reference <- kit200Chain(kit200T6Matrix(FALSE), c(A = 4L, B = 4L, C = 4L), 1, 1)
    kit200ExpectChainEqual(chain$da$tables, reference$da$tables, "P_one")
    kit200ExpectFits(chain$da, "P_one")
})

test_that("KIT-200 T6 (b): replicate filter skipped, GUI defaults, groups of 3, 3 and 6, two removed samples", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    sizes <- c(A = 3L, B = 3L, C = 6L)
    build <- function(with_row) {
        x <- kit200SetScatter(kit200Decorrelate(kit200Matrix(sizes), c("S1", "S4")), kit200Scatter[-1L])
        if (with_row) kit200AddRow(x, "P_two", c(S1 = 22, S4 = 22)) else x
    }
    chain <- kit200Chain(build(TRUE), sizes, 2, 2)
    kit200ExpectT6Preconditions(chain, "P_two", c("S1", "S4"))
    if (!kit200Completed(chain$da, "T6 (b) DA error")) return(invisible())
    kit200ExpectChainUntested(chain$da$tables, "P_two")
    kit200ExpectChainEqual(chain$da$tables, kit200Chain(build(FALSE), sizes, 2, 2)$da$tables, "P_two")
    kit200ExpectFits(chain$da, "P_two")
})

test_that("KIT-200 T6 (c): label-free and TMT owners through T6's steps", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    for (type in c("LFQ", "TMT")) {
        chain <- kit200Chain(kit200T6Matrix(), c(A = 4L, B = 4L, C = 4L), 1, 1, type = type)
        kit200ExpectT6Preconditions(chain, "P_one", "S1")
        if (!kit200Completed(chain$da, paste("T6 (c)", type, "DA error"))) next
        kit200ExpectChainUntested(chain$da$tables, "P_one")
        reference <- kit200Chain(kit200T6Matrix(FALSE), c(A = 4L, B = 4L, C = 4L), 1, 1, type = type)
        kit200ExpectChainEqual(chain$da$tables, reference$da$tables, "P_one")
    }
})

test_that("KIT-200 T6 (d): the default QC chain, normalisation with RUV skipped, four samples removed", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    sizes <- c(A = 6L, B = 6L, C = 6L)
    removed <- c("S1", "S2", "S7", "S8")
    x <- kit200SetScatter(kit200Decorrelate(kit200Matrix(sizes), removed), kit200Scatter[-1L])
    x <- kit200AddRow(x, "P_four", c(S1 = 22, S2 = 22, S7 = 22, S8 = 22))
    chain <- kit200Chain(x, sizes, 2, 2, default_chain = TRUE, normalise = TRUE)
    kit200ExpectT6Preconditions(chain, "P_four", removed)
    if (!kit200Completed(chain$da, "T6 (d) DA error")) return(invisible())
    kit200ExpectChainUntested(chain$da$tables, "P_four")
    # The same DA declaration on the chain's DA input without P_four's row (normalisation reads every row).
    reduced <- chain$qc$input
    reduced@protein_quant_table <- reduced@protein_quant_table[reduced@protein_quant_table$Protein.Ids != "P_four", ,
        drop = FALSE]
    reference <- kit200Da(kit200Owner(kit200Detach(reduced), state = "correlation_filtered"))
    expect_null(reference$error)
    kit200ExpectChainEqual(chain$da$tables, reference$tables, "P_four")
    kit200ExpectFits(chain$da, "P_four")
})

test_that("KIT-200 T7: the managed DA click on T6's chain completes on the worker", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    kit200ManagedMocks()
    chain <- kit200Chain(kit200T6Matrix(), c(A = 4L, B = 4L, C = 4L), 1, 1, managed = TRUE)
    kit200ExpectT6Preconditions(chain, "P_one", "S1")
    if (!kit200Completed(chain$da, "T7 managed DA error")) return(invisible())
    kit200ExpectChainUntested(chain$da$tables, "P_one")
    reference <- kit200Chain(kit200T6Matrix(FALSE), c(A = 4L, B = 4L, C = 4L), 1, 1)
    kit200ExpectChainEqual(chain$da$tables, reference$da$tables, "P_one")
    kit200ExpectFits(chain$da, "P_one")
})

kit200ProtStore <- function(object, env = parent.frame()) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = env), "kit-200", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    manager
}

kit200ProtDa <- function(manager, state) {
    spec <- readDIAStateMetadata(manager, state)$data
    step <- list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA, contrast = "a=grouptreated-groupcontrol",
        annotations = spec$components$protein_id_table$entity_id, name = "a", comparison = "a")
    run <- computeDIAReferenceDaContrast(manager, state, step, "iq")$run
    as.data.frame(materializeDIATable(manager$.context, run$value$tables$contrasts_results_table))
}

kit200ExpectDeclaredEqual <- function(got, ref, id_column, id) {
    rownames(got) <- got[[id_column]]
    rownames(ref) <- ref[[id_column]]
    ids <- rownames(ref)
    expect_false(id %in% ids)
    row <- got[id, , drop = FALSE]
    expect_true(all(is.na(unlist(row[intersect(kit200Columns, names(row))]))))
    expect_identical(row$AveExpr, NA_real_)
    for (column in intersect(kit200Columns, names(ref))) {
        expect_identical(got[ids, column], ref[ids, column], info = column)
    }
}

# test-kit129-run-path.R's kit129ProtParameters, copied.
kit200ProtParameters <- function(ruv_mode = "manual", ruv_k = 1) {
    list(norm_method = "quantile", ruv_mode = ruv_mode, auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = ruv_k,
        grouping = "group")
}

# test-prot-ruv-factors-in-model.R's ruvModelFixture, ruvModelObject, ruvModelPackageState and ruvModelDa, copied.
kit200RuvModelFixture <- function(seed = 20261118L, proteins = 300L, per_group = 6L, missing = 0.25) {
    set.seed(seed)
    groups <- rep(c("control", "treated"), each = per_group)
    n_samples <- length(groups)
    mu <- stats::rnorm(proteins, 22, 2)
    s <- sqrt(4 * 0.25^2 / stats::rchisq(proteins, 4))
    effect <- numeric(proteins)
    where <- sample.int(proteins, round(0.2 * proteins))
    effect[where] <- sample(rep(c(1, -1), length.out = length(where)))
    scale1 <- function(v) { v <- v - mean(v); v / stats::sd(v) }
    w <- scale1(stats::lm.fit(stats::model.matrix(~ factor(groups)), stats::rnorm(n_samples))$residuals)
    tau <- sqrt(3 * 4 * 0.25^2 / stats::qchisq(0.5, 4))
    x <- matrix(mu, proteins, n_samples) + outer(effect, as.numeric(groups == "treated")) +
        outer(stats::rnorm(proteins, 0, tau), w) + matrix(stats::rnorm(proteins * n_samples), proteins, n_samples) * s
    if (missing > 0) {
        lod <- stats::uniroot(function(l) mean(stats::plogis(l - x)) - missing, c(0, 60))$root
        x[stats::runif(length(x)) < stats::plogis(lod - x)] <- NA_real_
    }
    dimnames(x) <- list(sprintf("PG_%05d", seq_len(proteins)), sprintf("run_%03d", seq_len(n_samples)))
    keep <- rowSums(is.finite(x[, groups == "control"])) >= per_group / 2 &
        rowSums(is.finite(x[, groups == "treated"])) >= per_group / 2
    design <- data.frame(Run = colnames(x), group = groups, replicates = seq_len(n_samples),
        tech_rep_group = paste0("u", seq_len(n_samples)), stringsAsFactors = FALSE)
    list(x = x[keep, , drop = FALSE], design = design)
}

kit200RuvModelObject <- function(x, design, technical = "tech_rep_group") {
    quant <- as.data.frame(c(list(Protein.Group = rownames(x)), as.data.frame(x)), check.names = FALSE,
        stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(x), stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = technical,
        args = list(globalParameters = list(workflow_type = "LFQ")))
}

kit200RuvModelPackageState <- function(fixture, k, ruv_k_rule = "user", technical = "tech_rep_group") {
    normalised <- kit200RuvModelObject(fixture$x, fixture$design, technical)
    factors <- computeProtRuvModelFactors(normalised, ruv_k = k, ruv_k_rule = ruv_k_rule, ruv_grouping_variable = "group")
    corrected <- suppressMessages(ruvIII_C_Varying(protRuvGroupReplicateObject(normalised), ruv_grouping_variable = "group",
        ruv_number_k = factors$k, ctrl = factors$controls))
    list(state = attachProtRuvModelFactors(corrected, normalised, factors), factors = factors, corrected = corrected)
}

kit200RuvModelDa <- function(object, formula = "~ 0 + group") {
    result <- NULL
    utils::capture.output(result <- suppressWarnings(suppressMessages(differentialAbundanceAnalysis(object,
        contrasts_tbl = data.frame(contrasts = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE),
        formula_string = formula, da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))))
    result
}

test_that("KIT-200 T8: the declared designs (factors in the model, and none) complete with an empty protein", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    # (a) W in the model (test-prot-ruv-factors-in-model.R's builders, copied): the empty protein and the scatter in the
    # stored normalised matrix, which is the DA's input (decision M14).
    fixture <- kit200RuvModelFixture()
    fixture$x[2L, 1L] <- NA_real_
    fixture$x[3L, c(5L, 9L)] <- NA_real_
    fixture$x[4L, 12L] <- NA_real_
    reference_state <- kit200RuvModelPackageState(fixture, 1L)$state
    state <- reference_state
    stored <- state@args$ruv_uncorrected_quant
    state@args$ruv_uncorrected_quant <- rbind(stored, P_empty = rep(NA_real_, ncol(stored)))
    analysed <- state@args$ruv_uncorrected_quant[, as.character(state@design_matrix$Run), drop = FALSE]
    expect_true(all(is.na(analysed["P_empty", ])))
    expect_true(all(is.na(analysed[2L, 1L]), is.na(analysed[3L, c(5L, 9L)]), is.na(analysed[4L, 12L])))
    run <- tryCatch(list(value = kit200RuvModelDa(state)), error = function(e) list(error = conditionMessage(e)))
    if (kit200Completed(run, "T8 (a) error")) {
        got <- run$value$contrasts_results$results[[1L]]
        ref <- kit200RuvModelDa(reference_state)$contrasts_results$results[[1L]]
        got$id <- rownames(got)
        ref$id <- rownames(ref)
        kit200ExpectDeclaredEqual(got, ref, "id", "P_empty")
    }
    # (b) declared none (test-kit129-run-path.R:245-275's builder, copied) with the empty protein and the scatter
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    table <- object@protein_quant_table
    scatter <- list(c(2L, 1L), c(3L, 5L), c(3L, 9L), c(4L, 12L))
    scatter <- lapply(scatter, function(cell) c(table$Protein.Group[[cell[[1L]]]], sprintf("run_%04d", cell[[2L]])))
    for (cell in scatter) table[table$Protein.Group == cell[[1L]], cell[[2L]]] <- NA_real_
    empty <- table[1L, , drop = FALSE]
    empty$Protein.Group <- "P_empty"
    empty[setdiff(names(empty), "Protein.Group")] <- NA_real_
    object@protein_quant_table <- rbind(table, empty)
    identifiers <- object@protein_id_table[1L, , drop = FALSE]
    identifiers[] <- list("P_empty", "P_EMPTY", "SYNTHETIC_EMPTY", "")
    object@protein_id_table <- rbind(object@protein_id_table, identifiers)
    manager <- kit200ProtStore(object)
    none <- c(kit200ProtParameters("skip"), list(ruv_design = "none", ruv_k_rule = NULL, ruv_k_seed = NULL))
    declared <- kit200Quiet(withDIAStateScope(manager, function() {
        new <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), none)
        list(final = new$value$final, state = readDIAExecutionComponents(manager, new$value$final))
    }))
    input <- declared$state@protein_quant_table
    samples <- setdiff(names(input), "Protein.Group")
    expect_identical(sum(input$Protein.Group == "P_empty"), 1L)
    expect_true(all(is.na(unlist(input[input$Protein.Group == "P_empty", samples]))))
    for (cell in scatter) expect_true(is.na(input[input$Protein.Group == cell[[1L]], cell[[2L]]]), info = cell[[1L]])
    run <- tryCatch(list(value = kit200Quiet(withDIAStateScope(manager, function() kit200ProtDa(manager,
        declared$final)))), error = function(e) list(error = conditionMessage(e)))
    if (kit200Completed(run, "T8 (b) error")) {
        reduced <- declared$state
        reduced@protein_quant_table <- input[input$Protein.Group != "P_empty", , drop = FALSE]
        reference_manager <- kit200ProtStore(kit200Detach(reduced))
        ref <- kit200Quiet(withDIAStateScope(reference_manager, function() kit200ProtDa(reference_manager,
            getDIAStateReference(reference_manager))))
        kit200ExpectDeclaredEqual(run$value, ref, "Protein.Group", "P_empty")
    }
})

test_that("KIT-200 T9: the scripted wrapper daAnalysisWrapperFunction's DA step completes on T2's object", {
    # The owner's decision VAL-180 (option (a)): T9 asserts only that the wrapper's differential-abundance step
    # (runTestsContrasts) completes and reports the protein without values as not tested. The wrapper's later steps stop
    # on any data (R/func_prot_da_model_wrapper.R:272-275, a list nested one level too deep for getSignificantData); that
    # defect has its own ticket. The wrapper's RLE, PCA and value-count plots are replaced: they precede the DA step.
    run_tests <- get("runTestsContrasts", envir = asNamespace("MultiScholaR"))
    seen <- new.env(parent = emptyenv())
    testthat::local_mocked_bindings(plotRle = function(...) ggplot2::ggplot(), plotPca = function(...) ggplot2::ggplot(),
        plotNumOfValuesNoLog = function(...) ggplot2::ggplot(), runTestsContrasts = function(...) {
            seen$value <- tryCatch(run_tests(...), error = function(e) {
                seen$error <- conditionMessage(e)
                stop(e)
            })
            seen$value
        }, .package = "MultiScholaR")
    x <- kit200AddRow(kit200SetScatter(kit200Matrix()), "P_empty")
    object <- kit200Object(x, kit200Design())
    try(kit200Quiet(daAnalysisWrapperFunction(object,
        contrasts_tbl = data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids")), silent = TRUE)
    if (!kit200Completed(list(error = seen$error), "T9 DA step error")) return(invisible())
    expect_false(is.null(seen$value))
    expect_identical(seen$value$features_without_values, "P_empty")
    kit200ExpectUntested(seen$value$results, "P_empty")
})

test_that("KIT-200 T10: T6's chain is recorded and replayed to the same DA tables", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    chain <- kit200Chain(kit200T6Matrix(), c(A = 4L, B = 4L, C = 4L), 1, 1, publish = FALSE)
    kit200ExpectT6Preconditions(chain, "P_one", "S1")
    if (!kit200Completed(chain$da, "T10 DA error")) return(invisible())
    fixture <- chain$fixture
    manager <- fixture$workflow$state_manager
    json <- file.path(fixture$root, "analysis_record.json")
    writeWorkflowAnalysisRecord(manager, json)
    record <- jsonlite::read_json(json, simplifyVector = FALSE)
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(fixture$root, "."), copy)) == 0L)
    h <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), h)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- h$kitPrepareReplayProject(copy, c("initial", "protein_s4_initial"), api)
    replayed <- NULL
    kit200Quiet(replayWorkflowAnalysisRecord(record, prepared, character(), fixture$paths, withr::local_tempdir(),
        on_step = function(step, replay_manager) {
            if (identical(step$operation, "DA-contrast-output")) {
                reference <- replay_manager$.context$snapshot$da_results
                manifest <- withDIAStateScope(replay_manager, function()
                    resolveWorkflowStateValue(replay_manager$.context, reference))
                replayed <<- kit200ManifestTables(replay_manager, manifest)
            }
        }))
    expect_false(is.null(replayed))
    expect_identical(replayed, chain$da$tables)
})

# ---- T11 and T12 ---------------------------------------------------------------------------------------------------

test_that("KIT-200 T11: the lipidomics and metabolomics exported differentialAbundanceAnalysis complete", {
    design <- kit200Design()
    rownames(design) <- NULL
    x <- kit200AddRow(kit200SetScatter(kit200Matrix()), "F_empty")
    assay <- data.frame(Name = rownames(x), x, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(assay) <- NULL
    objects <- list(
        lipidomics = createLipidomicsAssayData(lipid_data = list(Assay1 = assay), design_matrix = design,
            lipid_id_column = "Name", sample_id = "Run", group_id = "group"),
        metabolomics = createMetaboliteAssayData(metabolite_data = list(Assay1 = assay), design_matrix = design,
            metabolite_id_column = "Name", sample_id = "Run", group_id = "group"))
    for (omic in names(objects)) {
        run <- tryCatch(list(value = kit200Quiet(differentialAbundanceAnalysis(list(Assay1 = objects[[omic]]),
            contrasts_tbl = data.frame(contrasts = kit200Strings[1:2], stringsAsFactors = FALSE),
            formula_string = "~ 0 + group", treat_lfc_cutoff = 0))), error = function(e) list(error = conditionMessage(e)))
        if (!kit200Completed(run, paste("T11", omic, "error"))) next
        for (table in run$value$Assay1@contrasts_results_table) {
            row <- table[table$Name == "F_empty", , drop = FALSE]
            expect_identical(nrow(row), 1L, info = omic)
            expect_true(all(is.na(unlist(row[intersect(kit200Columns, names(row))]))), info = omic)
        }
    }
})

test_that("KIT-200 T12: with the trend off the removal applies too (D2)", {
    design <- kit200Design()
    reference_x <- kit200SetScatter(kit200Matrix())
    x <- kit200AddRow(reference_x, "P_empty")
    # (a) robust off: a guard, equal on base and head
    run <- kit200Run(x, design, lfc = 0, trend = FALSE, robust = FALSE)
    if (kit200Completed(run, "T12 (a) error")) {
        reference <- kit200Run(reference_x, design, lfc = 0, trend = FALSE, robust = FALSE)$value
        kit200ExpectEqual(run$value$results, reference$results, rownames(reference_x))
    }
    # (b) robust on: head equals the analysis without the empty protein; on base a measurement, printed
    run <- kit200Run(x, design, lfc = 0, trend = FALSE, robust = TRUE)
    reference <- kit200Run(reference_x, design, lfc = 0, trend = FALSE, robust = TRUE)$value
    if (is.null(run$error)) {
        ids <- rownames(reference_x)
        largest <- max(vapply(names(reference$results), function(contrast) max(abs(
            run$value$results[[contrast]][ids, "raw_pvalue"] - reference$results[[contrast]][ids, "raw_pvalue"])), 0))
        same <- identical(lapply(run$value$results, function(table) table[ids, kit200Columns]),
            lapply(reference$results, function(table) table[ids, kit200Columns]))
        cat(sprintf("KIT-200 T12 (b) measurement: completed; equal to the reference: %s; largest raw_pvalue difference: %.6g\n",
            same, largest))
    } else {
        cat(sprintf("KIT-200 T12 (b) measurement: stopped: %s\n", run$error))
    }
    if (!kit200Completed(run, "T12 (b) error")) return(invisible())
    kit200ExpectEqual(run$value$results, reference$results, rownames(reference_x))
})

test_that("KIT-200 audit I05: without row names, every row keeps its input position as its label", {
    design <- kit200Design()
    reference_x <- kit200SetScatter(kit200Matrix())
    x <- rbind(reference_x[1:2, ], NA_real_, reference_x[-(1:2), ])
    rownames(x) <- NULL
    rownames(reference_x) <- NULL
    for (trend in c(TRUE, FALSE)) {
        run <- kit200Run(x, design, lfc = 0, trend = trend, robust = trend)
        if (!kit200Completed(run, sprintf("I05 error (trend %s)", trend))) next
        reference <- kit200Run(reference_x, design, lfc = 0, trend = trend, robust = trend)$value
        expect_identical(run$value$features_without_values, "3")
        for (contrast in kit200Strings) {
            got <- run$value$results[[contrast]]
            ref <- reference$results[[contrast]]
            expect_setequal(rownames(got), as.character(1:61))
            expect_true(all(is.na(unlist(got["3", kit200Columns]))), info = contrast)
            expect_identical(as.list(got[as.character(c(1:2, 4:61)), kit200Columns]),
                as.list(ref[as.character(1:60), kit200Columns]), info = contrast)
        }
    }
})
