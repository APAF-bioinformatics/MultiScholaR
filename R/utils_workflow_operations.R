# Operation runner for the workflow kit (KIT-004): identity, cache, projected hydrate, unchanged
# kernel, checkpoint through the class declaration, release. A declaration binds storage and
# placement only; the scientific kernels are called unchanged inside its `kernel` function. Nothing
# here reads GUI inputs, the global environment or a Shiny session.

#' Declare one workflow operation.
#' @param name Operation name used verbatim by the scientific identity (`normalization-ruv`,
#'   `pearson-pairs`, `pearson-selection`, `DA-contrast-output`), so saved products keep hitting.
#' @param input_class Declared class of the source checkpoint (KIT-002 registry).
#' @param fields,columns Components hydrated for the kernel and an optional column projection per table.
#' @param unit Independence unit: `row`, `complete_sample`, `feature_group`, `whole_matrix` or `store`.
#' @param kernel `function(object, parameters, store)` composing the unchanged scientific kernels.
#'   `store$source` is the source reference; `store$checkpoint(value, label, from)` saves an S4 of
#'   `output_class` as an owned reference (`from` defaults to the operation source);
#'   `store$table(table)` saves a table component; `store$product(value, name, product_type)` saves a
#'   derived product under the operation identity. The kernel's return value is the product saved
#'   under the identity (references in place of any S4 it checkpointed).
#' @param output_class Declared class of S4 values the kernel checkpoints, or NULL for a product-only
#'   operation.
#' @param parameters Names of the population-affecting parameters, in identity order.
#' @param product_type Product type of the value saved under the identity.
#' @param global_statistics Named character vector (statistic = reason) of pre-pass statistics a
#'   partitioned run computes once over the whole input before any batch; NULL when none.
#' @param merge Named character vector (component = rule) for the non-assay components a partitioned
#'   kernel writes: `design`, `concatenate` or `sum`. Required for a partitioned unit, NULL otherwise.
#' @param factors Store factors applied to the product identity.
#' @return Validated declaration.
#' @noRd
newWorkflowOperation <- function(name, input_class, fields, columns = list(), unit, kernel, output_class = NULL,
                                 parameters, product_type, global_statistics = NULL, merge = NULL, factors = list()) {
    validateWorkflowOperation(list(name = name, input_class = input_class, fields = fields, columns = columns,
        unit = unit, kernel = kernel, output_class = output_class, parameters = parameters,
        product_type = product_type, global_statistics = global_statistics, merge = merge, factors = factors))
}

#' Check an operation declaration against the class registry.
#' @return The declaration, unchanged.
#' @noRd
validateWorkflowOperation <- function(declaration) {
    scalar <- function(value) is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
    named <- function(value, allowed = NULL) is.null(value) || (is.character(value) && !is.null(names(value)) &&
        all(nzchar(names(value))) && (is.null(allowed) || all(value %in% allowed)))
    units <- c("row", "complete_sample", "feature_group", "whole_matrix", "store")
    valid <- is.list(declaration) && scalar(declaration$name) && scalar(declaration$input_class) &&
        scalar(declaration$unit) && declaration$unit %in% units &&
        is.character(declaration$fields) && length(declaration$fields) > 0L &&
        is.list(declaration$columns) && is.function(declaration$kernel) && length(formals(declaration$kernel)) == 3L &&
        is.character(declaration$parameters) && !anyDuplicated(declaration$parameters) &&
        scalar(declaration$product_type) && is.list(declaration$factors) &&
        (is.null(declaration$output_class) || scalar(declaration$output_class)) &&
        named(declaration$global_statistics) && named(declaration$merge, c("design", "concatenate", "sum"))
    if (!isTRUE(valid)) {
        cli::cli_abort("workflow operation declaration is invalid", class = "workflow_operation_invalid")
    }
    partitioned <- declaration$unit %in% c("row", "complete_sample", "feature_group")
    if (partitioned != !is.null(declaration$merge)) {
        cli::cli_abort("operation {.val {declaration$name}} needs merge rules exactly when its unit is partitioned",
            class = "workflow_operation_invalid")
    }
    input <- workflowClassDeclaration(declaration$input_class)
    storage <- unlist(input[c("tables", "lists", "matrices", "nested", "opaque")], use.names = FALSE)
    if (!all(declaration$fields %in% storage) || !all(names(declaration$columns) %in% declaration$fields)) {
        cli::cli_abort("operation {.val {declaration$name}} hydrates fields that {.cls {declaration$input_class}} does not store",
            class = "workflow_operation_invalid")
    }
    if (!is.null(declaration$output_class)) workflowClassDeclaration(declaration$output_class)
    declaration
}

#' Run a declared operation on a saved source: identity, cache, projected hydrate, kernel, checkpoint.
#' @param manager Current DIA state owner.
#' @param declaration Declaration from `newWorkflowOperation`.
#' @param source Saved reference of the input checkpoint.
#' @param parameters Named list holding every declared parameter (identity) and any extra values
#'   the kernel needs; the identity uses the declared names in declared order.
#' @param dependencies Further references the product depends on (rooted with the source).
#' @return `list(status, value, reference)`: `hit` reuses the saved product without computing;
#'   `computed` ran the kernel on one projected object and saved the product.
#' @noRd
runSavedOperation <- function(manager, declaration, source, parameters, dependencies = list()) {
    validateWorkflowOperation(declaration)
    if (!isSavedWorkflowReference(source)) {
        cli::cli_abort("operation {.val {declaration$name}} needs a saved source reference",
            class = "workflow_operation_source")
    }
    missing <- setdiff(declaration$parameters, names(parameters))
    if (!is.list(parameters) || length(missing)) {
        cli::cli_abort("operation {.val {declaration$name}} is missing parameter{?s} {.val {missing}}",
            class = "workflow_operation_parameters")
    }
    spec <- readComponentMetadata(manager, source)$data
    if (!identical(spec$class, declaration$input_class)) {
        cli::cli_abort("operation {.val {declaration$name}} expects {.cls {declaration$input_class}}, the source is {.cls {spec$class}}",
            class = "workflow_operation_input")
    }
    identity <- diaComponentScientificIdentity(manager, source, declaration$name, parameters[declaration$parameters])
    hit <- resolveSavedWorkflowEntity(manager$.context$store, "scientific_product", identity, factors = declaration$factors)
    if (identical(hit$status, "hit")) {
        return(list(status = "hit", value = hit$value$value, reference = hit$reference))
    }
    store <- newWorkflowOperationStore(manager, declaration, source, identity)
    value <- runProjected(manager, declaration, source, parameters, store)
    reference <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowScientificProduct(value, declaration$product_type),
        basename(tempfile(paste0(declaration$name, "-"))), identity, "scientific_product", factors = declaration$factors)
    rootDIAProduct(manager, reference, c(list(source), dependencies, store$created()))
    list(status = "computed", value = value, reference = reference)
}

#' Whole-matrix placement: one projected object, the kernel, nothing else resident.
#' @inheritParams runSavedOperation
#' @param store Operation store handle, or NULL for a kernel that saves nothing.
#' @return The kernel's value.
#' @noRd
runProjected <- function(manager, declaration, source, parameters, store = NULL) {
    if (is.null(store)) store <- newWorkflowOperationStore(manager, declaration, source, NULL)
    object <- readDIAExecutionComponents(manager, source, fields = declaration$fields, columns = declaration$columns)
    declaration$kernel(object, parameters, store)
}

#' Store handle handed to a kernel: checkpoints through the output class declaration and records
#' every reference it creates so the product can be rooted on them.
#' @noRd
newWorkflowOperationStore <- function(manager, declaration, source, identity) {
    created <- list()
    remember <- function(reference) {
        created[[length(created) + 1L]] <<- reference
        reference
    }
    list(
        source = source,
        checkpoint = function(value, label, from = NULL) {
            if (is.null(declaration$output_class) || !methods::is(value, declaration$output_class)) {
                cli::cli_abort("operation {.val {declaration$name}} may only checkpoint {.cls {declaration$output_class}}",
                    class = "workflow_operation_output")
            }
            value <- inheritDIAExecutionSource(value, manager, from %||% source)
            remember(saveDIAOwnedValue(value, label, manager))
        },
        table = function(table) remember(saveDIATableComponent(manager$.context, table)),
        product = function(value, name, product_type) {
            if (is.null(identity)) {
                cli::cli_abort("derived products need the operation identity", class = "workflow_operation_output")
            }
            product_identity <- identity
            product_identity$selectors <- digest::digest(list(identity$selectors, name), algo = "sha256")
            remember(saveWorkflowConsumerEntity(manager$.context$store,
                newSavedWorkflowScientificProduct(value, product_type),
                basename(tempfile(paste0(declaration$name, "-", name, "-"))), product_identity, "scientific_product"))
        },
        created = function() created)
}
