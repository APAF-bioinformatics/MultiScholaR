# Class declarations for the workflow kit (KIT-002): the storage facts of a scientific S4 class
# that the component checkpoint needs. A declaration binds storage only; kernels and execution
# are untouched. Unknown classes are refused rather than guessed.

#' Declared workflow classes.
#' @return Named list of declarations: `class`, `version`, `assay` (slot name), `layout` (`long` or
#'   `wide`), `design` (slot name), `tables` (data.frame slots saved as table components), `matrices`
#'   (matrix slots saved as values), `roles` (identity column slots) and `describe` (NULL or a
#'   function producing the bounded state metadata).
#' @noRd
workflowClassDeclarations <- function() {
    list(
        PeptideQuantitativeData = list(class = "PeptideQuantitativeData", version = 1L,
            assay = "peptide_data", layout = "long", design = "design_matrix",
            tables = c("peptide_data", "design_matrix"), matrices = "peptide_matrix",
            roles = c("protein_id_column", "peptide_sequence_column", "sample_id", "group_id",
                "technical_replicate_id"),
            describe = NULL),
        ProteinQuantitativeData = list(class = "ProteinQuantitativeData", version = 1L,
            assay = "protein_quant_table", layout = "wide", design = "design_matrix",
            tables = c("protein_quant_table", "design_matrix", "protein_id_table"), matrices = character(),
            roles = c("protein_id_column", "sample_id", "group_id", "technical_replicate_id"),
            describe = describeDIAProteinState))
}

#' Resolve one class declaration.
#' @param class Class name (a character scalar or an object whose first class is declared).
#' @return Validated declaration.
#' @noRd
workflowClassDeclaration <- function(class) {
    if (!is.character(class)) class <- class(class)[[1L]]
    if (length(class) != 1L || is.na(class) || !nzchar(class)) {
        cli::cli_abort("a single class name is required", class = "workflow_class_undeclared")
    }
    declaration <- workflowClassDeclarations()[[class]]
    if (is.null(declaration)) {
        cli::cli_abort("class {.cls {class}} has no workflow declaration", class = "workflow_class_undeclared")
    }
    validateWorkflowClassDeclaration(declaration)
}

#' Check a declaration against the S4 definition it describes.
#' @return The declaration, invisibly unchanged.
#' @noRd
validateWorkflowClassDeclaration <- function(declaration) {
    scalar <- function(value) is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
    slots <- names(methods::getSlots(methods::getClass(declaration$class)))
    valid <- is.list(declaration) && scalar(declaration$class) && scalar(declaration$assay) &&
        scalar(declaration$design) && declaration$layout %in% c("long", "wide") &&
        is.integer(declaration$version) && length(declaration$version) == 1L && declaration$version >= 1L &&
        is.character(declaration$tables) && length(declaration$tables) >= 2L &&
        is.character(declaration$matrices) && is.character(declaration$roles) &&
        declaration$assay %in% declaration$tables && declaration$design %in% declaration$tables &&
        all(c(declaration$tables, declaration$matrices, declaration$roles) %in% slots) &&
        (is.null(declaration$describe) || is.function(declaration$describe))
    if (!isTRUE(valid)) {
        cli::cli_abort("workflow declaration for {.cls {declaration$class}} is invalid",
            class = "workflow_class_undeclared")
    }
    declaration
}

#' Components a declared class must carry in a component checkpoint.
#' @noRd
workflowClassRequiredComponents <- function(declaration) unique(c(declaration$assay, declaration$design))

#' The additive `declaration` field written into a component checkpoint record.
#' @noRd
workflowDeclarationField <- function(declaration) list(class = declaration$class, version = declaration$version)
