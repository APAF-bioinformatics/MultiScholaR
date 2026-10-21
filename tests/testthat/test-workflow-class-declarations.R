# KIT-002: class declarations and declaration-driven component checkpoints.
# Both DIA classes round-trip exactly; records without the additive `declaration` field read as before;
# the scientific identity digest ignores the field.

test_that("the two DIA classes are declared and unknown classes are refused", {
    peptide <- workflowClassDeclaration("PeptideQuantitativeData")
    expect_identical(peptide$assay, "peptide_data")
    expect_identical(peptide$layout, "long")
    expect_identical(peptide$design, "design_matrix")
    expect_null(peptide$describe)
    protein <- workflowClassDeclaration("ProteinQuantitativeData")
    expect_identical(protein$assay, "protein_quant_table")
    expect_identical(protein$layout, "wide")
    expect_setequal(protein$tables, c("protein_quant_table", "design_matrix", "protein_id_table"))
    expect_identical(protein$describe, describeDIAProteinState)
    expect_identical(workflowClassRequiredComponents(peptide), c("peptide_data", "design_matrix"))
    expect_identical(workflowDeclarationField(protein), list(class = "ProteinQuantitativeData", version = 1L))
    expect_error(workflowClassDeclaration("data.frame"), class = "workflow_class_undeclared")
    expect_error(workflowClassDeclaration(c("a", "b")), class = "workflow_class_undeclared")
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    expect_identical(workflowClassDeclaration(object)$class, "PeptideQuantitativeData")
    broken <- peptide
    broken$assay <- "no_such_slot"
    expect_error(validateWorkflowClassDeclaration(broken), class = "workflow_class_undeclared")
})

newDeclarationTestManager <- function(label) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = parent.frame()), label, "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = parent.frame())
    manager$.context$reference_states <- TRUE
    manager
}

test_that("both declared classes round-trip exactly through the declaration-driven checkpoint", {
    recipe <- readWorkflowCheckpointRecipe("small")
    for (class in c("PeptideQuantitativeData", "ProteinQuantitativeData")) {
        object <- if (class == "PeptideQuantitativeData") makeWorkflowCheckpointPeptide(recipe) else
            makeWorkflowCheckpointProtein(recipe)
        declaration <- workflowClassDeclaration(class)
        manager <- newDeclarationTestManager(paste0("roundtrip-", declaration$layout))
        manager$saveState("state", object, list(), "declared")
        record <- readDIAStateMetadata(manager)
        expect_identical(record$data$class, class, info = class)
        expect_identical(record$data$declaration, list(class = class, version = 1L), info = class)
        expect_true(all(workflowClassRequiredComponents(declaration) %in% names(record$data$components)), info = class)
        expect_false("declaration" %in% names(record$data$slots), info = class)
        expect_identical(record$data$dimensions$rows, nrow(methods::slot(object, declaration$assay)), info = class)
        if (is.function(declaration$describe)) {
            expect_identical(record$data$metadata, declaration$describe(object), info = class)
        } else expect_null(record$data$metadata, info = class)
        expect_equal(manager$getState(), object, tolerance = 1e-12, info = class)
        full <- readDIAExecutionComponents(manager, fields = declaration$tables)
        for (field in declaration$tables) expect_equal(methods::slot(full, field), methods::slot(object, field),
            tolerance = 1e-12, info = paste(class, field))
        projected <- readDIAExecutionComponents(manager, fields = declaration$design)
        expect_identical(nrow(methods::slot(projected, declaration$assay)), 0L, info = class)
        expect_identical(names(methods::slot(projected, declaration$assay)), names(methods::slot(object, declaration$assay)), info = class)
        expect_equal(methods::slot(projected, declaration$design), methods::slot(object, declaration$design), tolerance = 1e-12, info = class)
        expect_identical(attr(projected, "dia_execution_fields"), declaration$design, info = class)
    }
})

test_that("a record without the declaration field reads as before and identities ignore the field", {
    recipe <- readWorkflowCheckpointRecipe("small")
    object <- makeWorkflowCheckpointProtein(recipe)
    manager <- newDeclarationTestManager("pre-ticket")
    manager$saveState("protein_replicate_filtered", object, list(), "input")
    with_field <- readDIAStateMetadata(manager)
    expect_identical(with_field$data$declaration$class, "ProteinQuantitativeData")
    # A checkpoint exactly as the pre-ticket build wrote it: the same record without the field.
    legacy <- with_field$data
    legacy$declaration <- NULL
    commitDIAReferenceState(manager, "legacy_record", legacy, list(), "pre-ticket shape")
    read_back <- readDIAStateMetadata(manager)
    expect_null(read_back$data$declaration)
    expect_identical(read_back$data$components, with_field$data$components)
    expect_identical(read_back$data$slots, with_field$data$slots)
    expect_equal(manager$getState(), object, tolerance = 1e-12)
    identity_with <- diaComponentScientificIdentity(manager, getDIAStateReference(manager, "protein_replicate_filtered"),
        "declaration-test", list(a = 1))
    identity_without <- diaComponentScientificIdentity(manager, getDIAStateReference(manager, "legacy_record"),
        "declaration-test", list(a = 1))
    expect_identical(identity_with, identity_without)
    # Unknown classes refuse at the metadata boundary with the existing class.
    foreign <- with_field$data
    foreign$class <- "LipidomicsAssayData"
    expect_error(commitDIAReferenceState(manager, "foreign", foreign, list(), "unknown") |> invisible() |>
        (\(reference) readDIAStateMetadata(manager, reference))(), class = "saved_entity_corrupt")
})

test_that("an undeclared S4 result is refused with a typed class before any component is written", {
    methods::setClass("DeclarationAuditForeign", representation(args = "list", x = "numeric"),
        where = environment())
    manager <- newDeclarationTestManager("undeclared")
    foreign <- methods::new("DeclarationAuditForeign", args = list(), x = 1)
    writes <- manager$.context$writes
    expect_error(checkpointDIAExecutionResult(manager, foreign), class = "workflow_class_undeclared")
    expect_identical(manager$.context$writes, writes)
})

# ---- KIT-003: list-of-tables slots, result classes, nested and opaque slots --------------------

# Opaque exactness: plots are compared by their built layer data (their environments are re-created
# on read); everything else by a serialize hash.
opaqueDigest <- function(value) {
    portable <- function(x) {
        if (inherits(x, "ggplot")) return(ggplot2::ggplot_build(x)$data)
        if (is.list(x) && !is.data.frame(x) && !isS4(x)) return(lapply(x, portable))
        x
    }
    digest::digest(portable(value), algo = "sha256")
}

# Plot objects capture the calling environment; building them in a function with forced arguments
# keeps that environment free of lazy bindings, which the entity store refuses.
newDeclarationPlot <- function(table, x, y = NULL) {
    force(table); force(x); force(y)
    if (is.null(y)) ggplot2::ggplot(table, ggplot2::aes(x = .data[[x]])) + ggplot2::geom_histogram(bins = 5L) else
        ggplot2::ggplot(table, ggplot2::aes(x = .data[[x]], y = .data[[y]])) + ggplot2::geom_point()
}

newDeclarationDaResults <- function(class, object, assay_slot) {
    tables <- methods::slot(object, assay_slot)
    first <- tables[[1L]]
    sample_columns <- intersect(names(first), as.character(object@design_matrix[[object@sample_id]]))
    matrix <- as.matrix(first[, sample_columns, drop = FALSE])
    rownames(matrix) <- as.character(seq_len(nrow(matrix)))
    groups <- object@design_matrix[[object@group_id]][match(sample_columns, object@design_matrix[[object@sample_id]])]
    design <- stats::model.matrix(~ 0 + factor(groups))
    colnames(design) <- levels(factor(groups))
    fit <- limma::eBayes(limma::lmFit(matrix, design))
    long <- data.frame(feature = rownames(matrix), log2FC = fit$coefficients[, 1L] - fit$coefficients[, 2L],
        stringsAsFactors = FALSE, row.names = NULL)
    wide <- data.frame(feature = rownames(matrix), a = fit$coefficients[, 1L], b = fit$coefficients[, 2L],
        stringsAsFactors = FALSE, row.names = NULL)
    plot <- newDeclarationPlot(long, "log2FC")
    methods::new(class, theObject = object, fit.eb = fit,
        contrasts_results_table = list(a_vs_b = long, b_vs_a = transform(long, log2FC = -log2FC)),
        num_sig_diff_exp_bar_plot = list(a_vs_b = plot), num_sig_diff_table = data.frame(contrast = "a_vs_b", n = 1L),
        volcano_plot = list(a_vs_b = plot), interactive_volcano_plot = list(), p_value_dist_plot = list(a_vs_b = plot),
        results_table_long = long, results_table_wide = wide)
}

test_that("list-assay classes and result classes are declared with slot kinds; the DA manifest is a declared output", {
    for (class in c("LipidomicsAssayData", "MetaboliteAssayData")) {
        declaration <- workflowClassDeclaration(class)
        expect_identical(declaration$kind, "assay", info = class)
        expect_true(declaration$assay %in% declaration$lists, info = class)
        expect_identical(declaration$layout, "wide", info = class)
        expect_identical(workflowClassRequiredComponents(declaration), c(declaration$assay, "design_matrix"), info = class)
    }
    for (class in c("LipidomicsDifferentialAbundanceResults", "MetabolomicsDifferentialAbundanceResults")) {
        declaration <- workflowClassDeclaration(class)
        expect_identical(declaration$kind, "result", info = class)
        expect_identical(declaration$nested, "theObject", info = class)
        expect_true("fit.eb" %in% declaration$opaque, info = class)
        expect_identical(workflowClassRequiredComponents(declaration), character(), info = class)
    }
    expect_identical(workflowClassDeclaration("da_results_for_enrichment")$lists, "da_data")
    expect_setequal(workflowClassDeclaration("EnrichmentResults")$opaque,
        c("enrichment_plots", "enrichment_plotly", "enrichment_summaries"))
    manifest <- workflowOutputDeclaration("dia.da_results.v1")
    expect_identical(manifest$kind, "manifest")
    expect_true("results" %in% manifest$fields)
    expect_error(workflowOutputDeclaration("dia.unknown.v1"), class = "workflow_class_undeclared")
    broken <- workflowClassDeclaration("LipidomicsAssayData")
    broken$lists <- c(broken$lists, "design_matrix")
    expect_error(validateWorkflowClassDeclaration(broken), class = "workflow_class_undeclared")
})

test_that("list-of-tables assays round-trip exactly with one component per element and projected elements", {
    for (class in c("LipidomicsAssayData", "MetaboliteAssayData")) {
        object <- if (class == "LipidomicsAssayData") module_ci_lipid_norm_object() else module_ci_metab_norm_object()
        declaration <- workflowClassDeclaration(class)
        assay <- declaration$assay
        elements <- methods::slot(object, assay)
        expect_true(length(elements) >= 2L, info = class)
        manager <- newDeclarationTestManager(paste0("list-", tolower(substr(class, 1, 5))))
        manager$saveState("state", object, list(), "declared")
        record <- readComponentMetadata(manager)
        component <- record$data$components[[assay]]
        expect_identical(names(component), names(elements), info = class)
        expect_true(all(vapply(component, isWorkflowTableReference, logical(1))), info = class)
        expect_false(assay %in% names(record$data$slots), info = class)
        expect_identical(record$data$dimensions$rows, sum(vapply(elements, nrow, integer(1))), info = class)
        expect_equal(manager$getState(), object, tolerance = 1e-12, info = class)
        full <- readDIAExecutionComponents(manager, fields = c(assay, "design_matrix"))
        expect_equal(methods::slot(full, assay), elements, tolerance = 1e-12, info = class)
        first <- names(elements)[[1L]]
        projected <- readDIAExecutionComponents(manager, fields = c(paste(assay, first, sep = "/"), "design_matrix"))
        got <- methods::slot(projected, assay)
        expect_identical(names(got), names(elements), info = class)
        expect_equal(got[[first]], elements[[first]], tolerance = 1e-12, info = class)
        for (other in setdiff(names(elements), first)) {
            expect_identical(nrow(got[[other]]), 0L, info = paste(class, other))
            expect_identical(names(got[[other]]), names(elements[[other]]), info = paste(class, other))
        }
    }
})

test_that("differential abundance result classes round-trip with nested, opaque, list and table slots", {
    for (class in c("LipidomicsDifferentialAbundanceResults", "MetabolomicsDifferentialAbundanceResults")) {
        lipid <- class == "LipidomicsDifferentialAbundanceResults"
        object <- if (lipid) module_ci_lipid_norm_object() else module_ci_metab_norm_object()
        assay <- if (lipid) "lipid_data" else "metabolite_data"
        results <- newDeclarationDaResults(class, object, assay)
        manager <- newDeclarationTestManager(paste0("da-", if (lipid) "lipid" else "metab"))
        # The assay object is saved first; the nested copy inside the result must reuse its tables.
        manager$saveState("assay", object, list(), "input")
        assay_components <- readComponentMetadata(manager)$data$components
        manager$saveState("da", results, list(), "results")
        record <- readComponentMetadata(manager)
        components <- record$data$components
        expect_true(isSavedWorkflowReference(components$theObject), info = class)
        expect_identical(components$theObject$payload$codec, "rds", info = class)
        nested <- resolveDIAStateValue(manager$.context, components$theObject)
        expect_identical(nested$schema, "dia.nested_checkpoint.v1", info = class)
        expect_identical(nested$data$class, class(object)[[1L]], info = class)
        expect_identical(lapply(nested$data$components[[assay]], `[[`, "entity_id"),
            lapply(assay_components[[assay]], `[[`, "entity_id"), info = paste(class, "nested tables stored once"))
        expect_identical(nested$data$components$design_matrix$entity_id, assay_components$design_matrix$entity_id, info = class)
        expect_true(all(c("fit.eb", "volcano_plot", "p_value_dist_plot", "num_sig_diff_exp_bar_plot") %in% names(components)), info = class)
        expect_false("interactive_volcano_plot" %in% names(components), info = paste(class, "empty opaque stays a prototype"))
        expect_identical(names(components$contrasts_results_table), c("a_vs_b", "b_vs_a"), info = class)
        restored <- manager$getState()
        expect_s4_class(restored, class)
        for (field in c("num_sig_diff_table", "results_table_long", "results_table_wide")) {
            expect_equal(methods::slot(restored, field), methods::slot(results, field), tolerance = 1e-12, info = paste(class, field))
        }
        expect_equal(restored@contrasts_results_table, results@contrasts_results_table, tolerance = 1e-12, info = class)
        expect_equal(restored@theObject, object, tolerance = 1e-12, info = class)
        for (field in c("fit.eb", "volcano_plot", "p_value_dist_plot", "num_sig_diff_exp_bar_plot")) {
            expect_identical(opaqueDigest(methods::slot(restored, field)), opaqueDigest(methods::slot(results, field)),
                info = paste(class, field))
        }
        expect_identical(restored@interactive_volcano_plot, list(), info = class)
        projected <- readDIAExecutionComponents(manager, fields = "results_table_long")
        expect_equal(projected@results_table_long, results@results_table_long, tolerance = 1e-12, info = class)
        expect_identical(nrow(projected@results_table_wide), 0L, info = class)
        expect_null(projected@fit.eb, info = paste(class, "unrequested opaque keeps the prototype"))
        # An unrequested nested object keeps its design (its own validity needs it) and zero-row assay elements.
        expect_equal(projected@theObject@design_matrix, object@design_matrix, tolerance = 1e-12, info = class)
        nested_assay <- methods::slot(projected@theObject, assay)
        expect_identical(names(nested_assay), names(methods::slot(object, assay)), info = class)
        expect_true(all(vapply(nested_assay, nrow, integer(1)) == 0L), info = paste(class, "unrequested nested assay is a prototype"))
    }
})

test_that("enrichment input and result classes round-trip", {
    object <- module_ci_lipid_norm_object()
    design <- object@design_matrix
    contrasts <- tibble::tibble(contrasts = c("a_vs_b", "b_vs_a"))
    da_table <- data.frame(feature = c("f1", "f2"), log2FC = c(1, -1), fdr_qvalue = c(0.01, 0.2), stringsAsFactors = FALSE)
    input <- methods::new("da_results_for_enrichment", contrasts = contrasts, design_matrix = design,
        da_data = list(a_vs_b = da_table, b_vs_a = NULL))
    enrichment <- createEnrichmentResults(contrasts)
    enrichment@enrichment_data <- list(a_vs_b = list(up = da_table, down = da_table[0, ]), b_vs_a = list(up = NULL, down = NULL))
    plot <- newDeclarationPlot(da_table, "log2FC", "fdr_qvalue")
    enrichment@enrichment_plots <- list(a_vs_b = list(up = plot, down = NULL))
    enrichment@enrichment_summaries <- list(a_vs_b = list(terms = 3L))
    manager <- newDeclarationTestManager("enrichment")
    manager$saveState("input", input, list(), "da input")
    record <- readComponentMetadata(manager)$data
    expect_true(isWorkflowTableReference(record$components$da_data$a_vs_b))
    expect_true(isSavedWorkflowReference(record$components$da_data$b_vs_a))
    expect_true(isWorkflowTableReference(record$components$contrasts))
    restored <- manager$getState()
    expect_s4_class(restored, "da_results_for_enrichment")
    expect_equal(restored@contrasts, contrasts)
    expect_equal(restored@design_matrix, design)
    expect_equal(restored@da_data$a_vs_b, da_table)
    expect_null(restored@da_data$b_vs_a)
    manager$saveState("enrichment", enrichment, list(), "enrichment")
    record <- readComponentMetadata(manager)$data
    expect_setequal(names(record$components), c("contrasts", "enrichment_data", "enrichment_plots", "enrichment_summaries"))
    restored <- manager$getState()
    expect_s4_class(restored, "EnrichmentResults")
    expect_equal(restored@contrasts, contrasts)
    expect_equal(restored@enrichment_data, enrichment@enrichment_data)
    expect_identical(opaqueDigest(restored@enrichment_plots), opaqueDigest(enrichment@enrichment_plots))
    expect_identical(restored@enrichment_summaries, enrichment@enrichment_summaries)
    expect_identical(restored@enrichment_plotly, list())
})

# ---- KIT-047: a declared args entry stored out of line as one owned value ------------------------

test_that("a declared owned args entry is saved out of line, reattached on hydrate and reused when unchanged (KIT-047)", {
    protein <- workflowClassDeclaration("ProteinQuantitativeData")
    expect_identical(protein$args_as_owned, "limpa_dpc_quant_results")
    expect_identical(workflowClassDeclaration("PeptideQuantitativeData")$args_as_owned, character())
    broken <- protein
    broken$args_as_owned <- c("a", "a")
    expect_error(validateWorkflowClassDeclaration(broken), class = "workflow_class_undeclared")
    broken$args_as_owned <- 1
    expect_error(validateWorkflowClassDeclaration(broken), class = "workflow_class_undeclared")

    recipe <- readWorkflowCheckpointRecipe("small")
    object <- makeWorkflowCheckpointProtein(recipe)
    # A synthetic entry of about 20 MiB, above the checkpoint's inline cap (8 MiB) and the field spill (64 KiB).
    entry <- list(big = seq_len(2.7e6) + 0.5, tag = "synthetic", nested = list(a = 1L, b = "x"))
    expect_gt(as.numeric(object.size(entry)), 20 * 1024^2)
    object@args$limpa_dpc_quant_results <- entry
    manager <- newDeclarationTestManager("owned-args")
    manager$saveState("protein_dpc_quant", object, list(), "owned")
    record <- readDIAStateMetadata(manager)
    stored <- record$data$slots$args$limpa_dpc_quant_results
    expect_true(isSavedWorkflowReference(stored))
    expect_lt(as.numeric(object.size(record$data$slots)), 1024^2)
    expect_false("limpa_dpc_quant_results" %in% names(record$data$support))
    expect_true(stored$entity_id %in% vapply(record$data$dependencies, `[[`, character(1), "entity_id"))
    # Hydrated for a kernel: the entry is the value again; a projected read keeps the reference.
    hydrated <- readDIAExecutionComponents(manager)
    expect_identical(hydrated@args$limpa_dpc_quant_results, entry)
    expect_identical(manager$getState()@args$limpa_dpc_quant_results, entry)
    projected <- readDIAExecutionComponents(manager, fields = "design_matrix")
    expect_true(isSavedWorkflowReference(projected@args$limpa_dpc_quant_results))
    # A later state of the same object (the shared protein steps) reuses the owned value; a changed entry is saved anew.
    again <- checkpointDIAExecutionResult(manager, inheritDIAExecutionSource(hydrated, manager, getDIAStateReference(manager)))
    expect_identical(again$slots$args$limpa_dpc_quant_results$entity_id, stored$entity_id)
    from_projection <- checkpointDIAExecutionResult(manager, inheritDIAExecutionSource(projected, manager, getDIAStateReference(manager)))
    expect_identical(from_projection$slots$args$limpa_dpc_quant_results$entity_id, stored$entity_id)
    changed <- hydrated
    changed@args$limpa_dpc_quant_results$tag <- "changed"
    renewed <- checkpointDIAExecutionResult(manager, inheritDIAExecutionSource(changed, manager, getDIAStateReference(manager)))
    expect_true(isSavedWorkflowReference(renewed$slots$args$limpa_dpc_quant_results))
    expect_false(identical(renewed$slots$args$limpa_dpc_quant_results$entity_id, stored$entity_id))
    commitDIAReferenceState(manager, "later", renewed, list(), "changed entry")
    expect_identical(manager$getState()@args$limpa_dpc_quant_results$tag, "changed")
    # An object without the entry (every DIA route state) checkpoints as before: no reference in its args.
    plain <- makeWorkflowCheckpointProtein(recipe)
    manager2 <- newDeclarationTestManager("no-owned-args")
    manager2$saveState("protein_replicate_filtered", plain, list(), "plain")
    args <- readDIAStateMetadata(manager2)$data$slots$args
    expect_false("limpa_dpc_quant_results" %in% names(args))
    expect_false(any(vapply(args, isSavedWorkflowReference, logical(1))))
    expect_equal(manager2$getState(), plain, tolerance = 1e-12)
})
