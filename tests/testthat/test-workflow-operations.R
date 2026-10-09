# KIT-004: the operation runner (identity, cache hit, projected hydrate, checkpoint through the
# declaration, refusals) and the three DIA whole-matrix declarations run with no Shiny session.

newOperationTestManager <- function(label) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = parent.frame()), label, "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = parent.frame())
    manager$.context$reference_states <- TRUE
    manager
}

proteinFields <- c("protein_quant_table", "design_matrix", "protein_id_table")

newCountingOperation <- function(counter, fields = proteinFields, columns = list(), output_class = NULL,
                                 kernel = NULL) {
    newWorkflowOperation("test-operation", "ProteinQuantitativeData", fields = fields, columns = columns,
        unit = "whole_matrix", output_class = output_class, product_type = "test-product", parameters = "alpha",
        kernel = kernel %||% function(object, parameters, store) {
            counter$calls <- counter$calls + 1L
            counter$object <- object
            list(rows = nrow(object@protein_quant_table), columns = names(object@protein_quant_table),
                samples = nrow(object@design_matrix), alpha = parameters$alpha, extra = parameters$extra)
        })
}

test_that("operation declarations are validated against the class registry", {
    kernel <- function(object, parameters, store) NULL
    declare <- function(...) newWorkflowOperation("op", "ProteinQuantitativeData", fields = proteinFields,
        unit = "whole_matrix", kernel = kernel, parameters = "a", product_type = "p", ...)
    expect_identical(declare()$unit, "whole_matrix")
    expect_error(newWorkflowOperation("op", "data.frame", fields = "x", unit = "whole_matrix", kernel = kernel,
        parameters = "a", product_type = "p"), class = "workflow_class_undeclared")
    expect_error(declare(output_class = "list"), class = "workflow_class_undeclared")
    expect_error(newWorkflowOperation("op", "ProteinQuantitativeData", fields = proteinFields, unit = "per_plex",
        kernel = kernel, parameters = "a", product_type = "p"), class = "workflow_operation_invalid")
    expect_error(newWorkflowOperation("op", "ProteinQuantitativeData", fields = "no_such_slot", unit = "whole_matrix",
        kernel = kernel, parameters = "a", product_type = "p"), class = "workflow_operation_invalid")
    expect_error(newWorkflowOperation("op", "ProteinQuantitativeData", fields = proteinFields, unit = "whole_matrix",
        kernel = function(object) NULL, parameters = "a", product_type = "p"), class = "workflow_operation_invalid")
    expect_error(declare(merge = c(design_matrix = "design")), class = "workflow_operation_invalid")
    expect_error(newWorkflowOperation("op", "ProteinQuantitativeData", fields = proteinFields, unit = "feature_group",
        kernel = kernel, parameters = "a", product_type = "p"), class = "workflow_operation_invalid")
    partitioned <- newWorkflowOperation("op", "ProteinQuantitativeData", fields = proteinFields, unit = "feature_group",
        kernel = kernel, parameters = "a", product_type = "p", merge = c(design_matrix = "design"),
        global_statistics = c(threshold = "dataset-wide quantile"))
    expect_identical(partitioned$merge, c(design_matrix = "design"))
    expect_error(newWorkflowOperation("op", "ProteinQuantitativeData", fields = proteinFields, unit = "feature_group",
        kernel = kernel, parameters = "a", product_type = "p", merge = c(design_matrix = "average")),
        class = "workflow_operation_invalid")
    for (declaration in list(diaNormalizationOperation(), diaPearsonPairsOperation(), diaPearsonSelectionOperation(),
        diaContrastOperation())) {
        expect_identical(validateWorkflowOperation(declaration), declaration)
        expect_identical(declaration$unit, "whole_matrix")
    }
    expect_identical(diaNormalizationOperation()$name, "normalization-ruv")
    expect_identical(diaContrastOperation()$factors, list(route = "iq", layer = "contrast_references"))
    expect_identical(diaContrastOperation()$parameters, c("formula", "qvalue", "lfc", "weights", "contrast", "annotations"))
})

test_that("the runner computes once under the existing identity and reuses the saved product on a hit", {
    manager <- newOperationTestManager("operations-hit")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    source <- getDIAStateReference(manager)
    counter <- new.env()
    counter$calls <- 0L
    declaration <- newCountingOperation(counter)
    withDIAStateScope(manager, function() {
        first <- runSavedOperation(manager, declaration, source, list(alpha = 1, extra = "kept"))
        expect_identical(first$status, "computed")
        expect_identical(counter$calls, 1L)
        expect_identical(first$value$rows, nrow(object@protein_quant_table))
        expect_identical(first$value$extra, "kept")
        identity <- diaComponentScientificIdentity(manager, source, "test-operation", list(alpha = 1))
        expect_identical(resolveSavedWorkflowEntity(manager$.context$store, "scientific_product", identity)$reference,
            first$reference)
        expect_true(first$reference$entity_id %in% names(manager$.context$snapshot$products))
        second <- runSavedOperation(manager, declaration, source, list(alpha = 1, extra = "ignored by the identity"))
        expect_identical(second$status, "hit")
        expect_identical(counter$calls, 1L)
        expect_identical(second$reference, first$reference)
        expect_identical(second$value, first$value)
        third <- runSavedOperation(manager, declaration, source, list(alpha = 2))
        expect_identical(third$status, "computed")
        expect_identical(counter$calls, 2L)
        expect_error(runSavedOperation(manager, declaration, source, list(beta = 1)), class = "workflow_operation_parameters")
        expect_error(runSavedOperation(manager, declaration, object, list(alpha = 1)), class = "workflow_operation_source")
    })
    peptide <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager$saveState("initial", peptide, list(), "peptide")
    withDIAStateScope(manager, function() {
        expect_error(runSavedOperation(manager, declaration, getDIAStateReference(manager), list(alpha = 1)),
            class = "workflow_operation_input")
    })
    expect_identical(manager$.context$full_s4_reads, 0L)
})

test_that("the runner hydrates only the declared fields and columns", {
    manager <- newOperationTestManager("operations-projection")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    source <- getDIAStateReference(manager)
    counter <- new.env()
    counter$calls <- 0L
    withDIAStateScope(manager, function() {
        design_only <- runSavedOperation(manager, newCountingOperation(counter, fields = "design_matrix"), source,
            list(alpha = "design"))
        expect_identical(design_only$value$rows, 0L)
        expect_identical(design_only$value$samples, nrow(object@design_matrix))
        expect_identical(design_only$value$columns, names(object@protein_quant_table))
        expect_identical(attr(counter$object, "dia_component_source", exact = TRUE), source)
        expect_identical(attr(counter$object, "dia_execution_fields", exact = TRUE), "design_matrix")
        # A column projection must name a hydrated field; the class validity keeps the sample columns together.
        expect_error(newCountingOperation(counter, fields = "design_matrix", columns = list(protein_quant_table = "x")),
            class = "workflow_operation_invalid")
    })
})

test_that("the store handle checkpoints through the output declaration and refuses other classes", {
    manager <- newOperationTestManager("operations-store")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    source <- getDIAStateReference(manager)
    peptide <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    kernel <- function(object, parameters, store) {
        if (identical(parameters$alpha, "wrong-class")) return(store$checkpoint(peptide, "wrong"))
        if (identical(parameters$alpha, "product-only")) return(store$checkpoint(object, "wrong"))
        subset <- object
        subset@protein_quant_table <- subset@protein_quant_table[1:5, ]
        reference <- store$checkpoint(subset, "subset-science")
        list(reference = reference, table = store$table(data.frame(id = 1:3)),
            derived = store$product(list(count = 5L), "count", "test-count"))
    }
    checkpointing <- newCountingOperation(NULL, output_class = "ProteinQuantitativeData", kernel = kernel)
    product_only <- newCountingOperation(NULL, kernel = kernel)
    withDIAStateScope(manager, function() {
        expect_error(runSavedOperation(manager, checkpointing, source, list(alpha = "wrong-class")),
            class = "workflow_operation_output")
        expect_error(runSavedOperation(manager, product_only, source, list(alpha = "product-only")),
            class = "workflow_operation_output")
        run <- runSavedOperation(manager, checkpointing, source, list(alpha = "ok"))
        expect_identical(run$status, "computed")
        expect_true(isSavedWorkflowReference(run$value$reference))
        record <- readDIAStateMetadata(manager, run$value$reference)$data
        expect_identical(record$dimensions$rows, 5L)
        expect_true(source$entity_id %in% names(record$dependencies))
        expect_identical(manager$.context$snapshot$gui_references[["subset-science"]], run$value$reference)
        expect_true(isWorkflowTableReference(run$value$table))
        identity <- diaComponentScientificIdentity(manager, source, "test-operation", list(alpha = "ok"))
        identity$selectors <- digest::digest(list(identity$selectors, "count"), algo = "sha256")
        derived <- resolveSavedWorkflowEntity(manager$.context$store, "scientific_product", identity)
        expect_identical(derived$status, "hit")
        expect_identical(derived$value$value, list(count = 5L))
        expect_identical(derived$value$product_type, "test-count")
        roots <- manager$.context$snapshot$product_dependencies[[run$reference$entity_id]]
        expect_true(all(c(source$entity_id, run$value$reference$entity_id, run$value$table$entity_id,
            run$value$derived$entity_id) %in% names(roots)))
    })
})

test_that("the three DIA declarations run in a fresh R session with no Shiny session and leave the global environment unchanged", {
    skip_if_not_installed("callr")
    root <- normalizePath(test_path("..", ".."))
    result <- callr::r(function(root, libraries) {
        .libPaths(unique(c(libraries, .libPaths())))
        # A source tree is loaded with pkgload; an installed suite or a check directory has no DESCRIPTION
        # at the root, so the installed package is used there.
        if (file.exists(file.path(root, "DESCRIPTION"))) {
            suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE, attach_testthat = FALSE))
        } else {
            suppressPackageStartupMessages(library(MultiScholaR))
        }
        ns <- asNamespace("MultiScholaR")
        helpers <- new.env(parent = ns)
        sys.source(file.path(root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), helpers)
        attach(ns, name = "MultiScholaR-namespace", warn.conflicts = FALSE)
        on.exit(detach("MultiScholaR-namespace"), add = TRUE)
        snapshot <- function() {
            names <- sort(ls(.GlobalEnv, all.names = TRUE))
            stats::setNames(lapply(names, function(name) digest::digest(get(name, envir = .GlobalEnv), algo = "sha256")), names)
        }
        before <- snapshot()
        root_dir <- tempfile("operations-")
        dir.create(root_dir)
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root_dir, "fresh-session", "DIA"),
            WorkflowState$new(audit_enabled = FALSE))
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        manager$.context$reference_states <- TRUE
        object <- helpers$makeWorkflowCheckpointProtein(helpers$readWorkflowCheckpointRecipe("small"))
        manager$saveState("protein_replicate_filtered", object, list(), "source")
        source <- getDIAStateReference(manager)
        groups <- as.character(unique(object@design_matrix$group))
        contrasts <- normaliseProtDaContrastsTable(data.frame(contrasts = paste0("group", groups[[2L]], "-group", groups[[1L]])))
        outcome <- withDIAStateScope(manager, function() {
            parameters <- list(norm_method = "quantile", ruv_mode = "manual", auto_percentage_min = NULL,
                auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
                adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1, grouping = "group")
            normalization <- runSavedOperation(manager, diaNormalizationOperation(), source, parameters)
            final <- normalization$value$final
            pairs <- runSavedOperation(manager, diaPearsonPairsOperation(), final,
                list(group = "group", tech_rep_remove_regex = "pool"))
            selection <- runSavedOperation(manager, diaPearsonSelectionOperation(), final,
                list(pairs = pairs$reference$entity_id, threshold = 0.5, correlations = pairs$value),
                dependencies = list(pairs$reference))
            spec <- readDIAStateMetadata(manager, selection$value)$data
            contrast <- runSavedOperation(manager, diaContrastOperation(), selection$value,
                list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA,
                    contrast = contrasts$full_format[[1L]], annotations = spec$components$protein_id_table$entity_id,
                    name = contrasts$contrasts[[1L]], comparison = contrasts$friendly_names[[1L]]))
            list(statuses = c(normalization$status, pairs$status, selection$status, contrast$status),
                best_k = normalization$value$best_k,
                final_rows = readDIAStateMetadata(manager, final)$data$dimensions$rows,
                selected_samples = readDIAStateMetadata(manager, selection$value)$data$dimensions$samples,
                contrast_rows = contrast$value$rows, contrast_tables = names(contrast$value$tables),
                contrast_products = names(contrast$value$products),
                open_resources = manager$.context$store$open_resources, full_s4_reads = manager$.context$full_s4_reads)
        })
        list(before = before, after = snapshot(), outcome = outcome, shiny = !is.null(shiny::getDefaultReactiveDomain()))
    }, args = list(root, .libPaths()), timeout = 600)
    expect_identical(result$after, result$before)
    expect_false(result$shiny)
    expect_identical(result$outcome$statuses, rep("computed", 4L))
    expect_identical(result$outcome$best_k, 1L)
    expect_true(result$outcome$final_rows > 0L)
    expect_true(result$outcome$selected_samples > 0L)
    expect_identical(result$outcome$contrast_rows, result$outcome$final_rows)
    expect_true(all(c("da_proteins_long", "da_proteins_wide") %in% result$outcome$contrast_tables))
    expect_true("num_sig_da_molecules_first_go" %in% result$outcome$contrast_products)
    expect_identical(result$outcome$open_resources, 0L)
    expect_identical(result$outcome$full_s4_reads, 0L)
})
