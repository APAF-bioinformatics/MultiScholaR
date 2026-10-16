# KIT-016: the enrichment declaration on the test double and the DA display products on the fixture.

diaProductsFixtureManager <- function(object, root, state = "protein_replicate_filtered") {
    dir.create(root, showWarnings = FALSE)
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState(state, object, list(), "fixture")
    manager
}

withDIAProductsSession <- function(manager, expr) {
    expr <- substitute(expr)
    env <- parent.frame()
    out <- NULL
    shiny::testServer(function(input, output, session) {}, {
        session$userData$dia_capacity_manager <- manager
        out <<- withDIAStateScope(manager, function() eval(expr, env))
    })
    out
}

diaProductsManifest <- function(manager, contrast = "groupB-groupA", friendly = "groupB_vs_groupA") {
    object <- manager$getState()
    source <- getDIAStateReference(manager)
    ids <- object@protein_quant_table[[object@protein_id_column]]
    table <- withr::with_seed(16L, data.frame(check.names = FALSE, stringsAsFactors = FALSE,
        Protein.Group = ids, comparison = friendly, log2FC = stats::rnorm(length(ids)),
        raw_pvalue = stats::runif(length(ids)), fdr_qvalue = stats::runif(length(ids)),
        gene_name = paste0("G", seq_along(ids))))
    reference <- saveDIATableComponent(manager$.context, table)
    # The contrast's NumSig product as the DA declaration saves it (the publication kernel aggregates them).
    numsig <- data.frame(comparison = friendly, status = c("Up", "Down", "Not significant"),
        counts = c(3L, 2L, length(ids) - 5L), stringsAsFactors = FALSE)
    numsig_reference <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowScientificProduct(list(table = numsig), "dia-num_sig_da_molecules_first_go"),
        basename(tempfile("dia-numsig-")), diaComponentScientificIdentity(manager, source, "DA-numsig-fixture", list(contrast = contrast)),
        "scientific_product")
    list(schema = "dia.da_results.v1", source = source, protein_id_column = object@protein_id_column,
        sample_id = object@sample_id,
        contrasts = data.frame(contrasts = contrast, friendly_names = friendly, full_format = paste0(contrast, "=1"),
            stringsAsFactors = FALSE),
        results = stats::setNames(list(list(schema = "dia.contrast_reference.v1", source = source, contrast = contrast,
            comparison = friendly, tables = list(da_proteins_long = reference),
            products = list(num_sig_da_molecules_first_go = numsig_reference),
            columns = names(table), rows = nrow(table))), contrast),
        parameters = list(formula = "~ group", qvalue = 0.05, lfc = 0, weights = NA), table = table)
}

test_that("the enrichment declaration runs the kernel on the DA manifest's source state, records the service identity and is not cacheable (KIT-016)", {
    declaration <- diaEnrichmentOperation()
    expect_identical(declaration$unit, "whole_matrix")
    expect_identical(declaration$output_class, "EnrichmentResults")
    expect_false(declaration$factors$cacheable)
    expect_true("results" %in% declaration$parameters)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    manifest <- diaProductsManifest(manager)
    da_results <- module_ci_prot_enrich_da_results()
    pathway <- withr::local_tempdir()
    arguments <- list(da_results = da_results, taxon_id = 9606, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05,
        pathway_dir = pathway, go_annotations = NULL, exclude_iea = FALSE, protein_id_column = "uniprot_acc",
        contrast_names = "groupB-groupA", correction_method = "gSCS")
    workflow <- list2env(list(state_manager = manager))
    declared <- withDIAProductsSession(manager, runDIAReferenceEnrichment(workflow, manifest, arguments,
        process_fn = runProtEnrichDeterministicProcessEnrichments))
    expect_s4_class(declared$results, "EnrichmentResults")
    expect_identical(declared$service, "gprofiler2:gprofiler2_test_double")
    expect_true(isSavedWorkflowReference(declared$reference))
    # The pre-ticket kernel path on the same inputs gives the same enrichment data.
    direct <- do.call(runProtEnrichDeterministicProcessEnrichments, c(arguments, list(pathway_dir = withr::local_tempdir()))[names(arguments)])
    expect_identical(digest::digest(declared$results@enrichment_data), digest::digest(direct@enrichment_data))
    expect_identical(digest::digest(declared$results@enrichment_summaries), digest::digest(direct@enrichment_summaries))
    # The product carries the service identity and the result tables; a second run saves a new product (no cache hit).
    product <- hydrateSavedWorkflowEntity(manager$.context$store, declared$reference)$value
    expect_identical(product$service, declared$service)
    expect_identical(unname(product$results), manifest$results[[1]]$tables$da_proteins_long$entity_id)
    expect_identical(product$provenance$cacheable %||% FALSE, FALSE)
    products_before <- length(manager$.context$snapshot$products)
    again <- withDIAProductsSession(manager, runDIAReferenceEnrichment(workflow, manifest, arguments,
        process_fn = runProtEnrichDeterministicProcessEnrichments))
    expect_false(identical(again$reference$entity_id, declared$reference$entity_id))
    expect_identical(length(manager$.context$snapshot$products), products_before + 1L)
    # The results checkpoint as the KIT-003 result class on the adapter, as the module's completed state does.
    manager$saveState("enrichment_completed", declared$results, list(q_cutoff = 0.05), "fixture enrichment")
    expect_identical(manager$current_state, "enrichment_completed")
    restored <- manager$getState("enrichment_completed")
    expect_identical(digest::digest(restored@enrichment_data), digest::digest(declared$results@enrichment_data))
    expect_error(runDIAReferenceEnrichment(workflow, list(schema = "other"), arguments), class = "saved_entity_invalid")
})

test_that("the enrichment service identity reads a live gost result, the double and reports unavailable otherwise (KIT-016 audit PH-02)", {
    results_of <- function(data) methods::new("EnrichmentResults", contrasts = tibble::tibble(), enrichment_data = data,
        enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())
    gost <- list(result = data.frame(term_id = "GO:0008150", p_value = 0.01), meta = list(version = "e111_eg58_p18_f463989"))
    double <- list(result = data.frame(term_id = "GO:0008150", deterministic_backend = "gprofiler2_test_double"))
    expect_identical(diaEnrichmentServiceIdentity(results_of(list(c1 = list(up = gost, down = NULL)))),
        "gprofiler2:e111_eg58_p18_f463989")
    expect_identical(diaEnrichmentServiceIdentity(results_of(list(c1 = list(up = double)))), "gprofiler2:gprofiler2_test_double")
    expect_identical(diaEnrichmentServiceIdentity(results_of(list(c1 = list(up = gost, down = double), c2 = list(up = gost)))),
        c("gprofiler2:e111_eg58_p18_f463989", "gprofiler2:gprofiler2_test_double"))
    expect_identical(diaEnrichmentServiceIdentity(results_of(list())), "unavailable")
    expect_identical(diaEnrichmentServiceIdentity(results_of(list(c1 = list(up = NULL, down = list(result = NULL))))), "unavailable")
})

test_that("the static volcano display is a saved product of the contrast's result table, rebuilt on demand with unchanged content (KIT-016)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    manifest <- diaProductsManifest(manager)
    withDIAProductsSession(manager, {
        before <- length(manager$.context$snapshot$products)
        product <- generateDIAReferenceVolcanoProduct(manifest, "groupB_vs_groupA")
        expect_identical(product$schema, "dia_volcano_display_v1")
        expect_identical(product$rows, nrow(manifest$table))
        expect_setequal(product$columns, c("Protein.Group", "comparison", "log2FC", "raw_pvalue", "fdr_qvalue", "gene_name"))
        expect_equal(as.data.frame(product$data)[product$columns], manifest$table[product$columns], ignore_attr = TRUE)
        expect_identical(length(manager$.context$snapshot$products), before + 1L)
        # A second request is a cache hit (no new product); the rebuilt plot is the kernel's plot on the same data.
        again <- generateDIAReferenceVolcanoProduct(manifest, "groupB-groupA")
        expect_identical(length(manager$.context$snapshot$products), before + 1L)
        expect_identical(digest::digest(again$data), digest::digest(product$data))
        rebuilt <- generateProtDAVolcanoStatic(manifest, "groupB_vs_groupA", da_q_val_thresh = 0.05, lfc_threshold = 1)
        direct <- generateProtDAVolcanoStatic(list(da_proteins_long = manifest$table), "groupB_vs_groupA",
            da_q_val_thresh = 0.05, lfc_threshold = 1)
        expect_s3_class(rebuilt, "ggplot")
        expect_equal(as.data.frame(rebuilt$data), as.data.frame(direct$data), ignore_attr = TRUE)
        # The heatmap was already a saved product (KIT-004): its identity resolves the same way.
        expect_true(is.function(generateDIAReferenceHeatmap))
    })
})

test_that("the interactive volcano page is a saved rendered asset the module shows without loading Glimma (KIT-016)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    manifest <- diaProductsManifest(manager)
    page <- "<html><body><div id='glimma'>volcano</div></body></html>"
    withDIAProductsSession(manager, {
        expect_null(resolveDIAInteractiveVolcanoAsset(manifest, "groupB_vs_groupA", 0.05))
        reference <- saveDIAInteractiveVolcanoAsset(manager, manifest, "groupB-groupA", 0.05, charToRaw(page))
        expect_true(isSavedWorkflowReference(reference))
        expect_identical(resolveDIAInteractiveVolcanoAsset(manifest, "groupB_vs_groupA", 0.05), page)
        # Another threshold has no asset: the builder would rebuild the widget (not exercised here).
        expect_null(resolveDIAInteractiveVolcanoAsset(manifest, "groupB_vs_groupA", 0.01))
        loaded_before <- "Glimma" %in% loadedNamespaces()
        shown <- generateProtDAVolcanoPlotGlimma(manifest, "groupB_vs_groupA", da_q_val_thresh = 0.05, args_row_id = "Protein.Group")
        expect_s3_class(shown, "shiny.tag")
        expect_identical(shown$name, "iframe")
        expect_identical(shown$attribs$srcdoc, page)
        expect_identical("Glimma" %in% loadedNamespaces(), loaded_before)
    })
})

test_that("the DA publication saves the interactive pages as assets when it gets the owner (KIT-016)", {
    skip_if_not_installed("Glimma")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    manifest <- diaProductsManifest(manager)
    directory <- withr::local_tempdir()
    files <- withDIAProductsSession(manager, writeDIAInteractiveVolcanoOutputs(manifest, directory,
        data.frame(), 0.05, manager = manager))
    expect_true(any(grepl("groupB_vs_groupA[.]html$", files)))
    withDIAProductsSession(manager, {
        page <- resolveDIAInteractiveVolcanoAsset(manifest, "groupB_vs_groupA", 0.05)
        expect_true(is.character(page) && nchar(page) > 1000)
        expect_identical(page, rawToChar(readBin(grep("groupB_vs_groupA[.]html$", files, value = TRUE), "raw",
            n = file.size(grep("groupB_vs_groupA[.]html$", files, value = TRUE)))))
    })
})

test_that("a display table is saved once and paged from the store, with a column filter (KIT-018)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    table <- data.frame(term = paste0("GO:", seq_len(25)), p_value = seq_len(25) / 100,
        directionality = rep(c("positive", "negative"), length.out = 25), stringsAsFactors = FALSE)
    withDIAProductsSession(manager, {
        expect_null(saveDIADisplayTable(manager, NULL))
        expect_null(saveDIADisplayTable(manager, table[0L, ]))
        reference <- saveDIADisplayTable(manager, table)
        expect_true(isWorkflowTableReference(reference))
        expect_identical(saveDIADisplayTable(manager, table)$entity_id, reference$entity_id)
        page <- readDIADisplayTablePage(reference, limit = 10L)
        expect_identical(page$rows, 25L)
        expect_identical(nrow(page$data), 10L)
        expect_identical(names(page$data), names(table))
        expect_identical(page$data$term, table$term[1:10])
        up <- readDIADisplayTablePage(reference, filter = list(column = "directionality", value = "positive"))
        expect_identical(up$rows, 13L)
        expect_identical(up$data$term, table$term[table$directionality == "positive"])
        second <- readDIADisplayTablePage(reference, columns = "term", filter = list(column = "directionality", value = "negative"),
            offset = 10L, limit = 5L)
        expect_identical(second$rows, 12L)
        expect_identical(second$data$term, table$term[table$directionality == "negative"][11:12])
        none <- readDIADisplayTablePage(reference, filter = list(column = "directionality", value = "other"))
        expect_identical(none$rows, 0L)
        expect_identical(nrow(none$data), 0L)
        expect_error(readDIADisplayTablePage(list(), limit = 1L), class = "saved_entity_invalid")
    })
})

test_that("the RUV optimisation result is read from the accepted normalisation product (KIT-018)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    expect_null(resolveDIARuvOptimizationResult(NULL))
    expect_null(resolveDIARuvOptimizationResult(WorkflowState$new(audit_enabled = FALSE)))
    withDIAProductsSession(manager, {
        expect_null(resolveDIARuvOptimizationResult(manager))
        ruv <- list(optimization_results = data.frame(percentage = c(0.1, 0.2), k = c(1L, 2L), separation_score = c(0.5, 0.7)),
            best_percentage = 0.2, best_k = 2L)
        source <- getDIAStateReference(manager)
        identity <- diaComponentScientificIdentity(manager, source, "normalization-fixture", list())
        product <- saveWorkflowConsumerEntity(manager$.context$store,
            newSavedWorkflowScientificProduct(list(ruv_result = ruv), "dia-normalization-ruv"),
            basename(tempfile("dia-norm-")), identity, "scientific_product")
        final <- readDIAStateMetadata(manager, source)$data
        final$products$normalization_science <- product
        commitDIAReferenceState(manager, "ruv_corrected", final, list(ruv_k = 2L), "fixture normalisation")
        expect_identical(manager$current_state, "ruv_corrected")
        resolved <- resolveDIARuvOptimizationResult(manager)
        expect_identical(resolved$optimization_results, ruv$optimization_results)
        expect_identical(resolved$best_percentage, 0.2)
    })
})

diaProductsCommitManifest <- function(manager, manifest) {
    context <- manager$.context
    reference <- saveDIAStateValue(context, manifest, "da-result-references", codec = "rds")
    snapshot <- context$snapshot
    snapshot$da_results <- reference
    publishDIAStateHead(context, snapshot)
    reference
}

test_that("the DA publication worker handler stages the exports, figures and pages under the request directory without saving assets (KIT-019)", {
    skip_if_not_installed("Glimma")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    manifest <- diaProductsManifest(manager)
    manifest_reference <- diaProductsCommitManifest(manager, manifest)
    head <- getDIAStateReference(manager)
    annotations <- data.frame(Entry = manifest$table$Protein.Group[1:3], gene_names = c("A", "B", "C"), stringsAsFactors = FALSE)
    annotations_reference <- saveDIADisplayTable(manager, annotations)
    directory <- withr::local_tempdir()
    request <- list(operation = "da_publication", input_references = list(state = head),
        parameters = list(annotations = annotations_reference), directory = directory, worker_cap = 12 * 1024^3)
    products_before <- length(manager$.context$snapshot$products)
    value <- suppressMessages(computeWorkflowDaPublicationAction(manager, request, "whole"))
    expect_identical(value$reference$entity_id, head$entity_id)
    expect_identical(value$state_reference$entity_id, head$entity_id)
    expect_true(isDIAComponentCheckpoint(value$checkpoint))
    expect_identical(value$manifest_reference$entity_id, manifest_reference$entity_id)
    expect_true(any(grepl("^da_output_dir/da_proteins_groupB-groupA_long_annot[.]tsv$", value$files)))
    expect_true(any(grepl("^publication_graphs_dir/Volcano_Plots/groupB-groupA[.]png$", value$files)))
    expect_true(any(grepl("^publication_graphs_dir/Interactive_Volcano_Plots/groupB_vs_groupA[.]html$", value$files)))
    expect_true(any(grepl("NumSigDaMolecules/da_proteins_num_sig_da_molecules[.]tab$", value$files)))
    # The annotations joined: the export carries the gene names column.
    export <- vroom::vroom(file.path(directory, "da-publication", "da_output_dir", "da_proteins_groupB-groupA_long_annot.tsv"),
        show_col_types = FALSE)
    expect_true("gene_names" %in% names(export))
    # No asset and no product saved by the worker; the head is unchanged.
    expect_identical(length(manager$.context$snapshot$products), products_before)
    expect_identical(getDIAStateReference(manager)$entity_id, head$entity_id)
    expect_lt(as.numeric(object.size(value)), 8388608)
    # Not the head: refused before any work.
    stale <- request
    stale$input_references$state <- manifest$results[[1]]$tables$da_proteins_long
    expect_error(computeWorkflowDaPublicationAction(manager, stale, "whole"), class = "workflow_execution_stale")
    # The completion republishes the head unchanged and names the head state.
    workflow <- list2env(list(state_manager = manager))
    completed <- completeDIADaPublication(workflow, value)
    expect_identical(completed$stateReference$entity_id, head$entity_id)
    expect_identical(getDIAStateReference(manager)$entity_id, head$entity_id)
    # The controller's publication moves the staged files into the project and saves the pages as assets.
    root <- withr::local_tempdir()
    request_dir <- file.path(root, "request")
    dir.create(request_dir)
    file.copy(file.path(directory, "da-publication"), request_dir, recursive = TRUE)
    owner <- list(manager = manager, control = list(generation = 1L))
    result <- list(operation = "da_publication", state_reference = head)
    withDIAProductsSession(manager, suppressMessages(publishWorkflowDaResult(owner, list(request_dir = request_dir, root = root), result)))
    published <- list.files(file.path(request_dir, "required-output"), recursive = TRUE)
    expect_true(any(grepl("da_proteins_groupB-groupA_long_annot[.]tsv$", published)))
    expect_true(any(grepl("Volcano_Plots/all_volcano_plots_combined[.]pdf$", published)))
    expect_true(any(grepl("Interactive_Volcano_Plots/groupB_vs_groupA[.]html$", published)))
    withDIAProductsSession(manager, {
        page <- resolveDIAInteractiveVolcanoAsset(manifest, "groupB_vs_groupA", 0.05)
        expect_true(is.character(page) && nchar(page) > 1000)
        expect_identical(page, rawToChar(readBin(file.path(directory, "da-publication", "publication_graphs_dir",
            "Interactive_Volcano_Plots", "groupB_vs_groupA.html"), "raw", n = 1e8)))
        read <- readManagedDaResult(manager, result, NULL)
        expect_identical(read$manifest$source$entity_id, head$entity_id)
        expect_identical(read$resultText, diaDaStatusText(1L))
        expect_identical(read$stateReference$entity_id, head$entity_id)
    })
})

test_that("the worker entry point runs the DA publication from a frozen request (KIT-019; README pre-launch chain)", {
    skip_if_not_installed("Glimma")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    manifest <- diaProductsManifest(manager)
    diaProductsCommitManifest(manager, manifest)
    head <- getDIAStateReference(manager)
    store <- manager$.context$store
    directory <- withr::local_tempdir()
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = "da_publication",
        input_references = list(state = head), parameters = list(annotations = NULL),
        token = list(request_id = "wfx-test"), owner = store$identity, project_root = dirname(store$directory),
        directory = directory, config = list(mode = "adapters"), worker_cap = 12 * 1024^3)
    path <- file.path(directory, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::local_envvar(MULTISCHOLAR_WORKFLOW_WORKER = "1")
    closeDIAStateAdapter(manager)
    expect_null(suppressMessages(workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256")))))
    result <- readRDS(file.path(directory, "worker-result.rds"))
    expect_identical(result$status, "complete")
    expect_identical(result$result$operation, "da_publication")
    expect_identical(result$result$reference$entity_id, head$entity_id)
    expect_true(file.exists(file.path(directory, "da-publication", "publication_graphs_dir", "Interactive_Volcano_Plots", "groupB_vs_groupA.html")))
})

test_that("the active owner resolves from the open state scope before the session (KIT-022)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    other <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(other))
    expect_null(getDIAActiveManager())
    seen <- withDIAStateScope(manager, function() {
        inner <- withDIAStateScope(other, function() getDIAActiveManager())
        list(inner = inner, outer = getDIAActiveManager())
    })
    expect_identical(seen$inner, other)
    expect_identical(seen$outer, manager)
    expect_null(getDIAActiveManager())
    # An error inside the scope still pops it.
    expect_error(withDIAStateScope(manager, function() stop("boom")), "boom")
    expect_null(getDIAActiveManager())
    # A kit-level reader works inside the scope with no session at all.
    manifest <- diaProductsManifest(manager)
    page <- withDIAStateScope(manager, function() readDIAResultPage(manifest, "groupB_vs_groupA", limit = 5L))
    expect_identical(nrow(page$data), 5L)
    # The session's owner still serves the modules outside any scope.
    withDIAProductsSession(manager, expect_identical(getDIAActiveManager(), manager))
})

test_that("the server-side page filter serves store pages to the table widget and leaves in-memory tables to the widget's filter (KIT-024)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    table <- data.frame(term = paste0("GO:", seq_len(23)), p_value = seq_len(23) / 100, stringsAsFactors = FALSE)
    withDIAProductsSession(manager, {
        reference <- saveDIADisplayTable(manager, table)
        pager <- new.env(parent = emptyenv())
        filter <- diaDisplayTableFilter(pager)
        params <- function(start, length, draw = 3L) list(start = start, length = length, draw = draw, escape = "true",
            columns = lapply(names(table), function(name) list(name = name, searchable = "false",
                search = list(value = "", regex = "false"))),
            order = list(), search = list(value = "", regex = "false", smart = "true", caseInsensitive = "true"))
        # In memory: the widget's own filter (first page of the data it was given).
        plain <- filter(table, params(0L, 10L))
        expect_identical(plain$recordsTotal, 23L)
        expect_identical(nrow(plain$data), 10L)
        # Installed pager: pages from the store, with the table's row count.
        pager$rows <- reference$metadata$shape[[1L]]
        pager$read <- function(offset, limit) readDIADisplayTablePage(reference, offset = offset, limit = limit)$data
        page <- filter(data.frame(), params(10L, 10L))
        expect_identical(page$draw, 3L)
        expect_identical(page$recordsTotal, 23L)
        expect_identical(page$recordsFiltered, 23L)
        expect_identical(page$data[[1L]], table$term[11:20])
        expect_identical(page$DT_rows_current, 11:20)
        last <- filter(data.frame(), params(20L, 10L))
        expect_identical(last$data[[1L]], table$term[21:23])
        beyond <- filter(data.frame(), params(30L, 10L))
        expect_identical(nrow(beyond$data), 0L)
        all_rows <- filter(data.frame(), params(0L, -1L))
        expect_identical(nrow(all_rows$data), 23L)
        expect_true(all(c("pageLength", "serverSide", "ordering", "searching") %in% names(diaDisplayTableOptions(25L))))
    })
})

test_that("the enrichment method is inferred from the results object (KIT-024)", {
    results_of <- function(data) methods::new("EnrichmentResults", contrasts = tibble::tibble(), enrichment_data = data,
        enrichment_plots = list(), enrichment_plotly = list(), enrichment_summaries = list())
    gost <- list(result = data.frame(term_id = "GO:1"), meta = list(version = "e111"))
    expect_identical(diaEnrichmentMethodOf(results_of(list(c1 = list(up = gost)))), "gprofiler2")
    expect_identical(diaEnrichmentMethodOf(results_of(list())), "gprofiler2")
    fake <- structure(list(), class = "enrichResult")
    expect_identical(diaEnrichmentMethodOf(results_of(list(c1 = list(up = NULL, down = gost)))), "gprofiler2")
})

test_that("after an Open the RUV grid and the enrichment tables resolve from the accepted state without rerunning their steps (KIT-024)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    withDIAProductsSession(manager, {
        # The accepted normalisation product on a reopened owner: the module holds no table reference.
        ruv <- list(optimization_results = data.frame(percentage = c(0.1, 0.2, 0.3), k = 1:3, separation_score = c(0.5, 0.7, 0.6)),
            best_percentage = 0.2, best_k = 2L)
        source <- getDIAStateReference(manager)
        product <- saveWorkflowConsumerEntity(manager$.context$store,
            newSavedWorkflowScientificProduct(list(ruv_result = ruv), "dia-normalization-ruv"),
            basename(tempfile("dia-norm-")), diaComponentScientificIdentity(manager, source, "normalization-fixture", list()), "scientific_product")
        final <- readDIAStateMetadata(manager, source)$data
        final$products$normalization_science <- product
        commitDIAReferenceState(manager, "ruv_corrected", final, list(ruv_k = 2L), "fixture normalisation")
        norm_data <- shiny::reactiveValues(ruv_optimization_table = NULL, ruv_optimization_result = NULL)
        reference <- resolveProtNormRuvOptimizationTable(norm_data)
        expect_true(isWorkflowTableReference(reference))
        expect_identical(readDIADisplayTablePage(reference)$data$k, 1:3)
        expect_identical(shiny::isolate(norm_data$ruv_optimization_table)$entity_id, reference$entity_id)
        # The enrichment results from the completed state: the module holds nothing after an Open.
        da_results <- module_ci_prot_enrich_da_results()
        results <- runProtEnrichDeterministicProcessEnrichments(da_results = da_results, taxon_id = 9606, up_cutoff = 0,
            down_cutoff = 0, q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), go_annotations = NULL, exclude_iea = FALSE,
            protein_id_column = "uniprot_acc", contrast_names = "groupB-groupA", correction_method = "gSCS")
        manager$saveState("enrichment_completed", results, list(q_cutoff = 0.05), "fixture enrichment")
        enrichment_data <- shiny::reactiveValues(all_enrichment_results = NULL, enrichment_results_full = NULL, analysis_complete = FALSE)
        workflow <- shiny::reactiveValues(state_manager = manager)
        expect_true(restoreProtEnrichResultsFromState(enrichment_data, workflow))
        restored <- shiny::isolate(enrichment_data$all_enrichment_results)
        expect_identical(names(restored), names(results@enrichment_data))
        expect_true(isWorkflowTableReference(restored[[1]]$gprofiler_results_table))
        expect_identical(readDIADisplayTablePage(restored[[1]]$gprofiler_results_table)$rows, nrow(restored[[1]]$gprofiler_results))
        expect_true(shiny::isolate(enrichment_data$analysis_complete))
        expect_false(restoreProtEnrichResultsFromState(enrichment_data, workflow))
    })
})

test_that("the page readers take the owner from the render, because the widget's page requests run outside any reactive domain (KIT-017)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    table <- data.frame(term = paste0("GO:", seq_len(12)), stringsAsFactors = FALSE)
    reference <- withDIAProductsSession(manager, saveDIADisplayTable(manager, table))
    # Shiny serves the widget's page requests with no reactive domain and no open scope (the KIT-017 screen run:
    # the RUV grid's first page request failed and the browser's error dialog blocked the driver).
    expect_null(getDIAActiveManager())
    expect_error(readDIADisplayTablePage(reference, limit = 1L), class = "dia_owner_missing")
    pager <- new.env(parent = emptyenv())
    pager$rows <- 12L
    pager$read <- function(offset, limit) readDIADisplayTablePage(reference, offset = offset, limit = limit, manager = manager)$data
    page <- diaDisplayTableFilter(pager)(data.frame(), list(start = 10L, length = 10L, draw = 1L, escape = "true"))
    expect_identical(page$data[[1L]], table$term[11:12])
    expect_identical(page$recordsTotal, 12L)
    # The widget reads each row as an array (as DT's own filter serves it), not as an object keyed by column.
    expect_identical(as.character(jsonlite::toJSON(page$data, dataframe = "rows")), '[["GO:11"],["GO:12"]]')
})

test_that("a self-contained widget page is written in R with the same data and dependencies as the pandoc page, without pandoc (KIT-036)", {
    skip_if_not_installed("Glimma")
    skip_if_not(rmarkdown::pandoc_available(), "the pandoc reference page needs pandoc")
    set.seed(36)
    widget <- Glimma::glimmaXY(x = stats::rnorm(20), y = stats::rnorm(20), html = NULL)
    inline <- file.path(withr::local_tempdir(), "inline.html")
    reference <- file.path(withr::local_tempdir(), "pandoc.html")
    htmlwidgets::saveWidget(widget, reference, selfcontained = TRUE)
    testthat::local_mocked_bindings(pandoc_convert = function(...) stop("pandoc must not run"), .package = "rmarkdown")
    written <- saveSelfContainedWidget(widget, inline)
    expect_identical(attr(written, "writer"), "inline")
    read <- function(path) readChar(path, file.size(path), useBytes = TRUE)
    page <- read(inline)
    # The widget's data payload is the same JSON in both pages.
    # The widget's element id is drawn at render time; the JSON payload itself must be identical.
    payload <- function(text) sub("data-for=\"[^\"]+\"", "data-for=\"\"",
        regmatches(text, gregexpr("<script type=\"application/json\" data-for=\"[^\"]+\">.*?</script>", text, perl = TRUE))[[1L]])
    expect_identical(payload(page), payload(read(reference)))
    # Every dependency's script and stylesheet is inlined verbatim; nothing is referenced by path.
    deps <- htmltools::resolveDependencies(htmltools::renderTags(htmltools::as.tags(widget, standalone = TRUE))$dependencies)
    for (dependency in deps) {
        root <- dependency$src$file
        if (!is.null(dependency$package)) root <- system.file(root, package = dependency$package)
        for (name in c(dependency$script, dependency$stylesheet))
            expect_true(grepl(readChar(file.path(root, name), file.size(file.path(root, name)), useBytes = TRUE), page, fixed = TRUE, useBytes = TRUE), info = name)
    }
    expect_false(grepl("<script src=", page, fixed = TRUE))
    expect_false(grepl("<link href=", page, fixed = TRUE))
    expect_match(page, "^<!DOCTYPE html>")
    # The pandoc writer stays available behind the option.
    withr::local_options(multischolar.widget_pandoc = TRUE)
    testthat::local_mocked_bindings(saveWidget = function(widget, file, ...) { writeLines("pandoc", file); invisible(file) }, .package = "htmlwidgets")
    expect_identical(attr(saveSelfContainedWidget(widget, inline), "writer"), "pandoc")
})

test_that("the head records who produced each state and under what software, and nothing on the state record (KIT-025)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- diaProductsFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    snapshot <- manager$.context$snapshot
    expect_identical(snapshot$step_provenance$initial$operation, "import")
    software <- workflowSoftwareIdentity()
    expect_match(software$build, "^[0-9a-f]{64}$")
    installed <- nzchar(system.file("R", "MultiScholaR.rdb", package = "MultiScholaR"))
    if (installed) expect_identical(software$build, workflowExecutionBuildIdentity())
    expect_true(is.character(software$source_commit) && length(software$source_commit) == 1L)
    expect_true(all(c("MultiScholaR", "limma") %in% names(software$packages)))
    expect_match(software$r_version, "^R version")
    expect_identical(names(snapshot$software_environments), software$build)
    # A commit with no pending provenance is legacy; one after an operation set it carries the operation.
    withDIAProductsSession(manager, {
        head <- getDIAStateReference(manager)
        record <- readDIAStateMetadata(manager, head)
        checkpoint <- record$data
        before <- digest::digest(list(checkpoint$slots, checkpoint$dimensions, record$config), algo = "sha256")
        commitDIAReferenceState(manager, "legacy_state", checkpoint, record$config, "no operation")
        setDIAStepProvenance(manager, "protein_intensity_filter", "whole", c("strict_mode", "min_reps_per_group"),
            fixed = list(normalization = "none"))
        commitDIAReferenceState(manager, "declared_state", checkpoint, record$config, "declared operation")
        commitDIAReferenceState(manager, "after_state", checkpoint, record$config, "pending consumed")
    })
    provenance <- manager$.context$snapshot$step_provenance
    expect_setequal(names(provenance), c("initial", manager$.context$snapshot$state_history[-1]))
    expect_identical(provenance$legacy_state$operation, "legacy")
    expect_identical(provenance$declared_state$operation, "protein_intensity_filter")
    expect_identical(provenance$declared_state$variant, "whole")
    expect_identical(provenance$declared_state$parameter_names, c("strict_mode", "min_reps_per_group"))
    expect_identical(provenance$declared_state$fixed, list(normalization = "none"))
    expect_identical(provenance$declared_state$build, software$build)
    expect_match(provenance$declared_state$defaults_sha256, "^[0-9a-f]{64}$")
    expect_true(is.na(provenance$declared_state$seed))
    expect_identical(provenance$after_state$operation, "legacy")
    # The state record is untouched: slots, dimensions and config hash as before.
    after <- withDIAProductsSession(manager, {
        r <- readDIAStateMetadata(manager, getDIAStateReference(manager, "declared_state"))
        digest::digest(list(r$data$slots, r$data$dimensions, r$config), algo = "sha256")
    })
    expect_identical(after, before)
    # A revert drops the provenance of the states it drops; a reopened store keeps its entries.
    withDIAProductsSession(manager, selectDIAStateReference(manager, "legacy_state"))
    expect_setequal(names(manager$.context$snapshot$step_provenance), c("initial", manager$.context$snapshot$state_history[-1]))
    expect_false("declared_state" %in% names(manager$.context$snapshot$step_provenance))
})
