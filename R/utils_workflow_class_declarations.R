# Class declarations for the workflow kit (KIT-002, KIT-003): the storage facts of a scientific S4
# class that the component checkpoint needs. A declaration binds storage only; kernels and
# execution are untouched. Unknown classes are refused rather than guessed.
#
# Slot kinds: `tables` (data.frame slots saved as table components), `lists` (named lists whose
# data.frame elements are saved as one table component each, other elements as values),
# `matrices` (matrix slots saved as values), `nested` (S4 slots saved as a reference to their own
# component checkpoint, never a second copy), `opaque` (slots saved as one serialized value: a
# limma fit, lists of plot objects), `roles` (identity column slots kept inline).

#' Declared workflow classes.
#' @return Named list of declarations: `class`, `version`, `kind` (`assay` or `result`), `assay`
#'   (slot name or NULL), `layout` (`long`, `wide` or NULL), `design` (slot name or NULL),
#'   `tables`, `lists`, `matrices`, `nested`, `opaque`, `roles` and `describe` (NULL or a function
#'   producing the bounded state metadata).
#' @noRd
workflowClassDeclarations <- function() {
    # KIT-047: `args_as_owned` names the entries of the `args` slot a checkpoint stores out of line as one owned
    # value each (a saved reference in `slots$args`, the value reattached on hydrate); the entry is absent from
    # every state that never carried it, so a class with no such entry checkpoints exactly as before.
    declare <- function(class, kind, assay, layout, design, tables, roles, lists = character(),
                        matrices = character(), nested = character(), opaque = character(), describe = NULL,
                        args_as_owned = character()) {
        list(class = class, version = 1L, kind = kind, assay = assay, layout = layout, design = design,
            tables = tables, lists = lists, matrices = matrices, nested = nested, opaque = opaque,
            roles = roles, describe = describe, args_as_owned = args_as_owned)
    }
    plots <- c("num_sig_diff_exp_bar_plot", "volcano_plot", "interactive_volcano_plot", "p_value_dist_plot")
    list(
        PeptideQuantitativeData = declare("PeptideQuantitativeData", "assay", "peptide_data", "long", "design_matrix",
            tables = c("peptide_data", "design_matrix"), matrices = "peptide_matrix",
            roles = c("protein_id_column", "peptide_sequence_column", "sample_id", "group_id", "technical_replicate_id")),
        ProteinQuantitativeData = declare("ProteinQuantitativeData", "assay", "protein_quant_table", "wide", "design_matrix",
            tables = c("protein_quant_table", "design_matrix", "protein_id_table"),
            roles = c("protein_id_column", "sample_id", "group_id", "technical_replicate_id"),
            describe = describeDIAProteinState, args_as_owned = "limpa_dpc_quant_results"),
        LipidomicsAssayData = declare("LipidomicsAssayData", "assay", "lipid_data", "wide", "design_matrix",
            tables = "design_matrix", lists = "lipid_data",
            roles = c("lipid_id_column", "annotation_id_column", "sample_id", "group_id", "technical_replicate_id")),
        MetaboliteAssayData = declare("MetaboliteAssayData", "assay", "metabolite_data", "wide", "design_matrix",
            tables = "design_matrix", lists = "metabolite_data",
            roles = c("metabolite_id_column", "annotation_id_column", "sample_id", "group_id", "technical_replicate_id")),
        LipidomicsDifferentialAbundanceResults = declare("LipidomicsDifferentialAbundanceResults", "result", NULL, NULL, NULL,
            tables = c("num_sig_diff_table", "results_table_long", "results_table_wide"),
            lists = "contrasts_results_table", nested = "theObject", opaque = c("fit.eb", plots), roles = character()),
        MetabolomicsDifferentialAbundanceResults = declare("MetabolomicsDifferentialAbundanceResults", "result", NULL, NULL, NULL,
            tables = c("num_sig_diff_table", "results_table_long", "results_table_wide"),
            lists = "contrasts_results_table", nested = "theObject", opaque = c("fit.eb", plots), roles = character()),
        da_results_for_enrichment = declare("da_results_for_enrichment", "result", NULL, NULL, "design_matrix",
            tables = c("contrasts", "design_matrix"), lists = "da_data", roles = character()),
        EnrichmentResults = declare("EnrichmentResults", "result", NULL, NULL, NULL,
            tables = "contrasts", lists = "enrichment_data",
            opaque = c("enrichment_plots", "enrichment_plotly", "enrichment_summaries"), roles = character()))
}

#' The declared class names, optionally of one kind (KIT-070: the owner's readers take their class list here).
#' @param kind `NULL` for every declaration, else `"assay"` or `"result"`.
#' @return Character vector of class names.
#' @noRd
workflowDeclaredClassNames <- function(kind = NULL) {
    declarations <- workflowClassDeclarations()
    if (!is.null(kind)) declarations <- Filter(function(d) identical(d$kind, kind), declarations)
    names(declarations)
}

#' Declared non-class outputs: bounded manifests a step delivers instead of an S4 object.
#' @return Named list by schema. The proteomics DA output is the existing `dia.da_results.v1`
#'   manifest (audit B11): a reference-only record whose contrast products are saved entities.
#' @noRd
workflowOutputDeclarations <- function() {
    # KIT-031: `operation` names the registry entry that produces the manifest and `long_table` the long result table
    # the analysis record lists per result (the record reads both from the declaration, not from the omic).
    list(`dia.da_results.v1` = list(schema = "dia.da_results.v1", version = 1L, kind = "manifest",
        producer = "runDIAReferenceDaAnalysis", operation = "DA-contrast-output", long_table = "da_proteins_long",
        fields = c("source", "protein_id_column", "sample_id", "contrasts", "results", "parameters"),
        products = c("da_proteins_long", "da_proteins_wide", "contrasts_results_table", "significant_rows",
            "norm_counts", "num_sig_da_molecules_first_go", "num_sig_da_molecules", "qvalue_warnings"),
        identifies = "isDIAResultManifest"),
        # KIT-031: the lipidomics DA output (the GUI never builds the S4 result class): one kernel call over every
        # contrast and assay; `results` is keyed by assay, each with its long table; `products$da_results` is the
        # kernel's list for the module's display and files.
        `lipid.da_results.v1` = list(schema = "lipid.da_results.v1", version = 1L, kind = "manifest",
        producer = "runLipidReferenceDaAnalysis", operation = "lipid_da_output", long_table = "da_lipids_long",
        fields = c("source", "lipid_id_column", "sample_id", "contrasts", "results", "tables", "products", "parameters", "model"),
        products = c("da_lipids_long", "da_results", "significant_counts", "qvalue_warnings"),
        identifies = "isLipidResultManifest"))
}

#' Resolve one declared output manifest by schema.
#' @noRd
workflowOutputDeclaration <- function(schema) {
    declaration <- workflowOutputDeclarations()[[schema]]
    if (is.null(declaration)) {
        cli::cli_abort("output schema {.val {schema}} has no workflow declaration", class = "workflow_class_undeclared")
    }
    declaration
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
    optional <- function(value) is.null(value) || scalar(value)
    slots <- names(methods::getSlots(methods::getClass(declaration$class)))
    kinds <- c("tables", "lists", "matrices", "nested", "opaque", "roles")
    named <- unlist(declaration[kinds], use.names = FALSE)
    valid <- is.list(declaration) && scalar(declaration$class) && declaration$kind %in% c("assay", "result") &&
        is.integer(declaration$version) && length(declaration$version) == 1L && declaration$version >= 1L &&
        optional(declaration$assay) && optional(declaration$design) &&
        (is.null(declaration$layout) || declaration$layout %in% c("long", "wide")) &&
        all(vapply(declaration[kinds], is.character, logical(1))) &&
        !anyDuplicated(named) && all(named %in% slots) &&
        (declaration$kind == "result" || (scalar(declaration$assay) && scalar(declaration$design) &&
            scalar(declaration$layout) && declaration$design %in% declaration$tables &&
            (declaration$assay %in% declaration$tables || declaration$assay %in% declaration$lists))) &&
        (is.null(declaration$design) || declaration$design %in% declaration$tables) &&
        (is.null(declaration$describe) || is.function(declaration$describe)) &&
        (is.null(declaration$args_as_owned) || (is.character(declaration$args_as_owned) &&
            !anyDuplicated(declaration$args_as_owned) && all(nzchar(declaration$args_as_owned)) &&
            (!length(declaration$args_as_owned) || "args" %in% slots)))
    if (!isTRUE(valid)) {
        cli::cli_abort("workflow declaration for {.cls {declaration$class}} is invalid",
            class = "workflow_class_undeclared")
    }
    declaration
}

#' Components a declared class must carry in a component checkpoint.
#' @return Slot names (an assay class needs its assay and design; a result class needs none).
#' @noRd
workflowClassRequiredComponents <- function(declaration) as.character(unique(c(declaration$assay, declaration$design)))

#' The additive `declaration` field written into a component checkpoint record.
#' @noRd
workflowDeclarationField <- function(declaration) list(class = declaration$class, version = declaration$version)

#' Bounded dimensions of a declared object's assay (rows, columns, samples).
#' A list assay reports the summed element rows and the union of element columns.
#' @noRd
workflowAssayDimensions <- function(object, declaration) {
    rows <- NA_integer_
    columns <- character()
    if (!is.null(declaration$assay)) {
        assay <- methods::slot(object, declaration$assay)
        if (is.data.frame(assay)) {
            rows <- nrow(assay)
            columns <- names(assay)
        } else if (is.list(assay)) {
            frames <- Filter(is.data.frame, assay)
            rows <- sum(vapply(frames, nrow, integer(1)))
            columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
        }
    }
    samples <- if (!is.null(declaration$design)) nrow(methods::slot(object, declaration$design)) else NA_integer_
    list(rows = rows, columns = columns, samples = samples)
}
