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
    expect_error(workflowClassDeclaration("LipidomicsAssayData"), class = "workflow_class_undeclared")
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
