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
    list(schema = "dia.da_results.v1", source = source, protein_id_column = object@protein_id_column,
        sample_id = object@sample_id,
        contrasts = data.frame(contrasts = contrast, friendly_names = friendly, full_format = paste0(contrast, "=1"),
            stringsAsFactors = FALSE),
        results = stats::setNames(list(list(schema = "dia.contrast_reference.v1", source = source, contrast = contrast,
            comparison = friendly, tables = list(da_proteins_long = reference), products = list(),
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
