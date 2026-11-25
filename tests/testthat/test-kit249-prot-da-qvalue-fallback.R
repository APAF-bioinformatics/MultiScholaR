# KIT-249 (DAS-058, DAS-064): when qvalue() fails for a comparison, runTestsContrasts' error handler computed the
# Benjamini-Hochberg fallback into a local copy of the q-value vector, so fdr_qvalue stayed NA for every protein of that
# comparison while qvalue_warnings, the DA tab's warning and the analysis record said BH was used. These tests assert
# the decided rule on a constructed comparison with known truth (ticket section 8): on B_vs_A no p-value reaches 0.95,
# so qvalue() fails (qvalue 2.44.0, pi0est.R:106) and fdr_qvalue must hold BH; C_vs_A is the success path. Each call is
# wrapped so that the tests-only commit fails by assertion. Synthetic data only; builders copied, never called.

kit249Full <- c(B_vs_A = "B_vs_A=groupB-groupA", C_vs_A = "C_vs_A=groupC-groupA")
kit249Plain <- c(B_vs_A = "groupB-groupA", C_vs_A = "groupC-groupA")
kit249Planted <- sprintf("P%02d", 1:6)

kit249Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

# The matrix of section 8 (no random numbers): groups A, B, C of 4 samples; P01..P06 planted (+2 log2 in B and C);
# P07..P30 true nulls whose B shift is a fixed realisation c_j sigma_j, c_j = 0.30 + 0.05 (j - 7).
kit249Matrix <- function() {
    j <- 1:30
    mu <- 20 + 0.25 * (j - 1)
    sigma <- 0.20 * (1 + 0.05 * ((j - 1) %% 5))
    beta <- ifelse(j <= 6, 2, (0.30 + 0.05 * (j - 7)) * sigma)
    gamma <- ifelse(j <= 6, 2, 0)
    e <- list(A = c(-1.5, -0.5, 0.5, 1.5), B = c(0.5, 1.5, -1.5, -0.5), C = c(1.5, -1.5, -0.5, 0.5))
    x <- cbind(mu + outer(sigma, e$A), mu + beta + outer(sigma, e$B), mu + gamma + outer(sigma, e$C))
    dimnames(x) <- list(sprintf("P%02d", j), paste0("S", 1:12))
    x
}

kit249Design <- function() {
    runs <- paste0("S", 1:12)
    data.frame(Run = runs, group = rep(c("A", "B", "C"), each = 4L), replicates = runs, row.names = runs,
        stringsAsFactors = FALSE)
}

kit249Object <- function(args = list(globalParameters = list(workflow_type = "LFQ"))) {
    x <- kit249Matrix()
    design <- kit249Design()
    rownames(design) <- NULL
    ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
            stringsAsFactors = FALSE), args = args)
}

kit249Run <- function(lfc) {
    tryCatch(list(value = kit249Quiet(runTestsContrasts(kit249Matrix(), unname(kit249Full), kit249Design(),
        "~ 0 + group", treat_lfc_cutoff = lfc, eBayes_trend = TRUE, eBayes_robust = TRUE))),
        error = function(e) list(error = conditionMessage(e)))
}

# The base-side guard: the call completed.
kit249Completed <- function(run, label) {
    expect_null(run$error, label = label)
    is.null(run$error)
}

kit249Ids <- function(table, id) if (!is.null(id) && id %in% names(table)) as.character(table[[id]]) else rownames(table)

# E1 (ii) to (vi) on one B_vs_A table (`bh`: the table's own BH column; NULL skips (iii)).
kit249ExpectFallback <- function(table, label, id = NULL, bh = "fdr_value_bh_adjustment", checks = c("ii", "iii", "iv", "v", "vi")) {
    table <- as.data.frame(table)
    p <- table$raw_pvalue
    q <- table$fdr_qvalue
    planted <- kit249Ids(table, id) %in% kit249Planted
    expect_identical(sum(planted), 6L, label = paste(label, "planted rows"))
    if ("ii" %in% checks) expect_identical(sum(is.na(q[is.finite(p)])), 0L, label = paste(label, "(ii) missing fdr_qvalue"))
    if ("iii" %in% checks && !is.null(bh)) expect_identical(q, table[[bh]], label = paste(label, "(iii) fdr_qvalue"))
    if ("iv" %in% checks) expect_identical(q, stats::p.adjust(p, "BH"), label = paste(label, "(iv) fdr_qvalue"))
    if ("v" %in% checks) expect_identical(sum(!is.na(q[planted]) & q[planted] < 0.05), 6L,
        label = paste(label, "(v) planted significant"))
    # (vi) a guard on both commits: NA compared is not a call
    if ("vi" %in% checks) expect_false(any(q < 0.05 & !planted, na.rm = TRUE), label = paste(label, "(vi) a null called"))
}

# G1: the success path (C_vs_A) keeps Storey's q-values, unflagged, and calls the planted proteins.
kit249ExpectStorey <- function(table, warnings, label, id = NULL) {
    table <- as.data.frame(table)
    p <- table$raw_pvalue
    expect_identical(table$fdr_qvalue, qvalue::qvalue(p)$qvalues, label = paste(label, "G1 fdr_qvalue"))
    expect_identical(length(warnings), 0L, label = paste(label, "G1 qvalue_warnings"))
    planted <- kit249Ids(table, id) %in% kit249Planted
    expect_true(all(table$fdr_qvalue[planted] < 0.05), label = paste(label, "G1 planted significant"))
}

# ---- E1, E1b, E2: the kernel --------------------------------------------------------------------------------------

test_that("KIT-249 E1: topTable branch, qvalue() fails on B_vs_A and fdr_qvalue holds the Benjamini-Hochberg values", {
    run <- kit249Run(NA)
    if (!kit249Completed(run, "E1 error")) return(invisible())
    # (i) a guard on both commits: the comparison is flagged
    expect_true(isTRUE(run$value$qvalue_warnings[[kit249Full[["B_vs_A"]]]]))
    kit249ExpectFallback(run$value$results[[kit249Full[["B_vs_A"]]]], "E1 B_vs_A")
    kit249ExpectStorey(run$value$results[[kit249Full[["C_vs_A"]]]], run$value$qvalue_warnings[[kit249Full[["C_vs_A"]]]],
        "E1 C_vs_A")
})

test_that("KIT-249 E1b: any qvalue() error, not only the pi0 one, keeps the Benjamini-Hochberg values", {
    run <- module_ci_prot_da_with_namespace_bindings(list(qvalue = function(p) stop("forced")), kit249Run(NA))
    if (!kit249Completed(run, "E1b error")) return(invisible())
    for (cmp in names(kit249Full)) {
        expect_true(isTRUE(run$value$qvalue_warnings[[kit249Full[[cmp]]]]), info = cmp)
        kit249ExpectFallback(run$value$results[[kit249Full[[cmp]]]], paste("E1b", cmp), checks = c("ii", "iii", "iv"))
    }
})

test_that("KIT-249 E2: topTreat branch (the application's TREAT cutoff 0), the same", {
    run <- kit249Run(0)
    if (!kit249Completed(run, "E2 error")) return(invisible())
    expect_true(isTRUE(run$value$qvalue_warnings[[kit249Full[["B_vs_A"]]]]))
    kit249ExpectFallback(run$value$results[[kit249Full[["B_vs_A"]]]], "E2 B_vs_A")
    kit249ExpectStorey(run$value$results[[kit249Full[["C_vs_A"]]]], run$value$qvalue_warnings[[kit249Full[["C_vs_A"]]]],
        "E2 C_vs_A")
})

# ---- E3: the S4 method --------------------------------------------------------------------------------------------

test_that("KIT-249 E3: S4 differentialAbundanceAnalysis on a label-free object, one comparison per call", {
    run <- tryCatch(list(value = kit249Quiet(differentialAbundanceAnalysis(kit249Object(),
        contrasts_tbl = data.frame(contrasts = kit249Full[["B_vs_A"]], stringsAsFactors = FALSE),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
        eBayes_robust = TRUE, args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))),
        error = function(e) list(error = conditionMessage(e)))
    if (!kit249Completed(run, "E3 error")) return(invisible())
    kit249ExpectFallback(run$value$contrasts_results_table, "E3 B_vs_A", id = "Protein.Ids", checks = c("ii", "iii", "v", "vi"))
    counts <- as.data.frame(run$value$num_sig_da_molecules)
    up <- counts$counts[counts$status == "Significant and Up"]
    expect_identical(if (length(up)) as.integer(up) else 0L, 6L, label = "E3 count table, Significant and Up")
})

# ---- E4 to E6: the application's DA tab (in process and managed) and the analysis record --------------------------

# An owner with the paths the DA publication writes to, holding the fixture as the correlation filter's state
# (test-kit200-prot-da-no-value-proteins.R's kit200Owner, copied; T6 (d)'s reference owner). Closed with `env`.
kit249Owner <- function(type = "DIA", env = parent.frame()) {
    root <- normalizePath(withr::local_tempdir(.local_envir = env))
    paths <- list(source_dir = file.path(root, "source"), protein_qc_dir = file.path(root, "protein_qc"),
        time_dir = file.path(root, "time"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, paste0("kit249-", tolower(type)), type),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(type)
    manager$.context$experiment_paths <- paths
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("correlation_filtered", kit249Object(list()), list(), "fixture")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list(),
        tab_status = list()))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = env)
    list(root = root, paths = paths, workflow = workflow, da = new.env(parent = emptyenv()))
}

kit249Contrasts <- function() data.frame(contrasts = names(kit249Full), full_format = unname(kit249Full),
    friendly_names = names(kit249Full), stringsAsFactors = FALSE)

# The managed DA click (test-kit200-prot-da-no-value-proteins.R's kit200Dispatch, kit200Dispatcher, kit200Managed and
# kit200ManagedMocks, which copy savedDaW2Dispatch, savedDaW2Managed, savedDaW2Dispatcher and savedDaW2Mocks of
# test-workflow-saved-da.R; copied) ---------------------------------------------------------------------------------
kit249Dispatch <- function(fixture, controller, operation, parameters, input) {
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

kit249Dispatcher <- function(fixture) {
    owner <- NULL
    controller <- NULL
    function(operation, parameters, input) {
        manager <- fixture$workflow$state_manager
        if (!identical(manager, owner)) {
            owner <<- manager
            controller <<- omicsWorkerController(list(owner = manager))
        }
        kit249Dispatch(fixture, controller, operation, parameters, input)
    }
}

kit249Managed <- function(fixture, dispatch, settings, contrasts) {
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

kit249ManagedMocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

kit249ManifestTables <- function(manager, manifest) {
    withDIAStateScope(manager, function() lapply(stats::setNames(nm = names(manifest$results)), function(name)
        lapply(c(results = "contrasts_results_table", long = "da_proteins_long"), function(field)
            as.data.frame(materializeDIATable(manager$.context, manifest$results[[name]]$tables[[field]])))))
}

# The DA tab's click on the owner's head: in process (`runDIAReferenceDaAnalysis`, the publication included) or managed.
kit249Da <- function(fixture, managed = FALSE) {
    manager <- fixture$workflow$state_manager
    contrasts <- kit249Contrasts()
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    fixture$da$current_s4_object <- getDIAStateReference(manager)
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    run <- tryCatch(list(value = kit249Quiet(if (managed)
        kit249Managed(fixture, kit249Dispatcher(fixture), settings, contrasts) else
        runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts, NULL))),
        error = function(e) list(error = conditionMessage(e)))
    if (is.null(run$error)) run$tables <- kit249ManifestTables(manager, run$value)
    run
}

# E4's assertions on one DA tab's manifest: B_vs_A's results table (ii), (iii), (iv), (v), (vi); its long table (ii),
# (v) and the same fdr_qvalue per protein; the manifest's qvalue_warnings; C_vs_A's G1.
kit249ExpectTab <- function(run, label) {
    b <- run$tables[["B_vs_A"]]
    expect_false(is.null(b), label = paste(label, "B_vs_A tables"))
    if (is.null(b)) return(invisible())
    expect_true(isTRUE(run$value$results[["B_vs_A"]]$qvalue_warnings[[kit249Full[["B_vs_A"]]]]),
        label = paste(label, "manifest qvalue_warnings names B_vs_A"))
    kit249ExpectFallback(b$results, paste(label, "B_vs_A results"), id = "Protein.Ids")
    kit249ExpectFallback(b$long, paste(label, "B_vs_A da_proteins_long"), id = "Protein.Ids", checks = c("ii", "v", "vi"))
    expect_identical(b$long$fdr_qvalue[match(b$results$Protein.Ids, b$long$Protein.Ids)], b$results$fdr_qvalue,
        label = paste(label, "B_vs_A long and results fdr_qvalue"))
    kit249ExpectStorey(run$tables[["C_vs_A"]]$results, run$value$results[["C_vs_A"]]$qvalue_warnings,
        paste(label, "C_vs_A"), id = "Protein.Ids")
}

test_that("KIT-249 E4 and E6: the DA tab in process on DIA, label-free and TMT owners, and the analysis record", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    for (type in c("DIA", "LFQ", "TMT")) {
        fixture <- kit249Owner(type)
        run <- kit249Da(fixture)
        if (!kit249Completed(run, paste("E4", type, "DA error"))) next
        kit249ExpectTab(run, paste("E4", type))
        if (!identical(type, "DIA")) next
        # E6: the record names BH for B_vs_A (a guard on both commits) and its saved fdr_qvalue holds BH values
        record <- buildWorkflowAnalysisRecord(fixture$workflow$state_manager)
        results <- record$analyses$differential_abundance$results
        methods <- stats::setNames(vapply(results, function(entry) entry$fdr_method %||% NA_character_, ""),
            vapply(results, function(entry) entry$contrast %||% NA_character_, ""))
        expect_identical(methods[["B_vs_A"]], "BH")
        expect_identical(methods[["C_vs_A"]], "qvalue")
        b <- run$tables[["B_vs_A"]]$results
        expect_identical(b$fdr_qvalue, b$fdr_value_bh_adjustment, label = "E6 B_vs_A saved fdr_qvalue")
    }
})

test_that("KIT-249 E5: the managed DA click on the worker", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    kit249ManagedMocks()
    run <- kit249Da(kit249Owner("DIA"), managed = TRUE)
    if (!kit249Completed(run, "E5 managed DA error")) return(invisible())
    kit249ExpectTab(run, "E5")
})

# ---- E7, E8: the scripted wrapper and the lipidomics and metabolomics exported S4 method --------------------------

# The kernel's value as runTestsContrasts returns it inside a route (KIT-200 T9's capture); the wrapper's RLE, PCA and
# value-count plots, which precede the DA step, are replaced as in T9.
kit249Capture <- function(code, plots = FALSE, env = parent.frame()) {
    run_tests <- get("runTestsContrasts", envir = asNamespace("MultiScholaR"))
    seen <- new.env(parent = emptyenv())
    bindings <- list(runTestsContrasts = function(...) {
        seen$value <- tryCatch(run_tests(...), error = function(e) {
            seen$error <- conditionMessage(e)
            stop(e)
        })
        seen$value
    })
    if (plots) bindings <- c(bindings, list(plotRle = function(...) ggplot2::ggplot(),
        plotPca = function(...) ggplot2::ggplot(), plotNumOfValuesNoLog = function(...) ggplot2::ggplot()))
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "MultiScholaR", .env = env)))
    result <- tryCatch(list(value = kit249Quiet(code())), error = function(e) list(error = conditionMessage(e)))
    list(result = result, kernel = seen$value, error = seen$error)
}

test_that("KIT-249 E7: the scripted wrapper protein_deAnalysisWrapperFunction", {
    captured <- kit249Capture(function() protein_deAnalysisWrapperFunction(kit249Object(),
        contrasts_tbl = tibble::tibble(contrasts = unname(kit249Full)), formula_string = "~ 0 + group",
        de_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids"), plots = TRUE)
    if (!kit249Completed(list(error = captured$error), "E7 DA step error")) return(invisible())
    expect_false(is.null(captured$kernel))
    kit249ExpectFallback(captured$kernel$results[[kit249Full[["B_vs_A"]]]], "E7 B_vs_A", checks = c("ii", "v", "vi"))
})

test_that("KIT-249 E8: the lipidomics and metabolomics exported differentialAbundanceAnalysis", {
    design <- kit249Design()
    rownames(design) <- NULL
    x <- kit249Matrix()
    assay <- data.frame(Name = rownames(x), x, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(assay) <- NULL
    objects <- list(
        lipidomics = createLipidomicsAssayData(lipid_data = list(Assay1 = assay), design_matrix = design,
            lipid_id_column = "Name", sample_id = "Run", group_id = "group"),
        metabolomics = createMetaboliteAssayData(metabolite_data = list(Assay1 = assay), design_matrix = design,
            metabolite_id_column = "Name", sample_id = "Run", group_id = "group"))
    for (omic in names(objects)) {
        run <- tryCatch(list(value = kit249Quiet(differentialAbundanceAnalysis(list(Assay1 = objects[[omic]]),
            contrasts_tbl = data.frame(contrasts = unname(kit249Plain), stringsAsFactors = FALSE),
            formula_string = "~ 0 + group", treat_lfc_cutoff = 0))), error = function(e) list(error = conditionMessage(e)))
        if (!kit249Completed(run, paste("E8", omic, "error"))) next
        kit249ExpectFallback(run$value$Assay1@contrasts_results_table[[kit249Plain[["B_vs_A"]]]], paste("E8", omic, "B_vs_A"),
            id = "Name", checks = c("ii", "iii"))
    }
})

test_that("KIT-249 E7 (b): the scripted wrapper daAnalysisWrapperFunction (KIT-218 merged), its returned tables", {
    contrasts <- data.frame(contrasts = unname(kit249Plain), full_format = unname(kit249Full), stringsAsFactors = FALSE)
    captured <- kit249Capture(function() daAnalysisWrapperFunction(kit249Object(), contrasts_tbl = contrasts,
        formula_string = "~ 0 + group", treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids"), plots = TRUE)
    if (!kit249Completed(list(error = captured$error), "E7 (b) DA step error")) return(invisible())
    if (!kit249Completed(captured$result, "E7 (b) wrapper error")) return(invisible())
    kit249ExpectFallback(captured$kernel$results[[kit249Full[["B_vs_A"]]]], "E7 (b) kernel B_vs_A", checks = c("ii", "v", "vi"))
    kit249ExpectFallback(captured$result$value$contrasts_results_table[[kit249Full[["B_vs_A"]]]], "E7 (b) returned B_vs_A",
        checks = c("ii", "iii", "v", "vi"))
})

# G1b (AUDIT-KIT-249-impl I03; a guard on both commits): qvalue() replaced by a double that succeeds with q = min(1, 1.5 p),
# so Storey's path is told from BH by a factor, not by rounding. A fix that assigned BH on the success path would put
# p.adjust(raw_pvalue, "BH") here (for the planted proteins about 30/6 p, not 1.5 p), and qvalue_warnings would stay empty:
# this expectation would then fail on every comparison of both branches.
test_that("KIT-249 G1b: a q-value that qvalue() returns is assigned as returned, on both branches", {
    for (lfc in c(NA, 0)) {
        run <- module_ci_prot_da_with_namespace_bindings(list(qvalue = \(p) list(q = pmin(1, 1.5 * p))), kit249Run(lfc))
        if (!kit249Completed(run, paste("G1b error", lfc))) next
        expect_identical(length(run$value$qvalue_warnings), 0L, label = paste("G1b qvalue_warnings", lfc))
        for (cmp in names(kit249Full)) {
            table <- run$value$results[[kit249Full[[cmp]]]]
            expect_identical(table$fdr_qvalue, pmin(1, 1.5 * table$raw_pvalue), label = paste("G1b", cmp, lfc))
        }
    }
})
