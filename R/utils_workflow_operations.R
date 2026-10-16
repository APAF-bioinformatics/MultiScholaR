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
#'   `store$source` is the source reference and `store$manager` the state owner (for kernels that run a
#'   reference entry point, KIT-006); `store$checkpoint(value, label, from)` saves an S4 of
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
                                 parameters, product_type, global_statistics = NULL, merge = NULL, factors = list(),
                                 fixed = NULL) {
    validateWorkflowOperation(list(name = name, input_class = input_class, fields = fields, columns = columns,
        unit = unit, kernel = kernel, output_class = output_class, parameters = parameters,
        product_type = product_type, global_statistics = global_statistics, merge = merge, factors = factors,
        fixed = fixed))
}

#' The registered operation behind a kernel section (KIT-025 audit PH-02): the in-process finishing steps name
#' their state's operation from the section the compute recorded, so a registered step run without a managed
#' application is attributed to its declaration rather than recorded as `legacy`.
#' @noRd
workflowSectionOperations <- c(srlQvalueProteotypicPeptideClean = "peptide_qvalue_filter",
    filterMinNumPeptidesPerProtein = "peptide_protein_evidence_filter", rollUpPrecursorToPeptide = "precursor_rollup",
    peptideIntensityFiltering = "peptide_intensity_filter", filterMinNumPeptidesPerSample = "minimum_peptides_per_sample",
    removePeptidesWithOnlyOneReplicate = "peptide_replicate_filter", peptideMissingValueImputation = "peptide_imputation",
    protein_rollup = "protein_iq_rollup", chooseBestProteinAccession = "protein_accession_cleanup",
    removeRowsWithMissingValuesPercent = "protein_intensity_filter", resolveProteinDuplicateRows = "protein_duplicate_removal",
    removeProteinsWithOnlyOneReplicate = "protein_replicate_filter")

#' @noRd
workflowOperationForSection <- function(section) {
    if (!is.character(section) || length(section) != 1L || is.na(section)) return(NULL)
    unname(workflowSectionOperations[section])[!is.na(workflowSectionOperations[section])]
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
        # KIT-025: a declaration may fix parameter values the kernel always uses; they are recorded on the head.
        (is.null(declaration$fixed) || (is.list(declaration$fixed) && !is.null(names(declaration$fixed)) &&
            all(nzchar(names(declaration$fixed))) && !anyDuplicated(names(declaration$fixed)))) &&
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
    # KIT-025: the state the wrapper commits after this operation carries the declaration's provenance.
    setDIAStepProvenance(manager, declaration$name, parameters$variant %||% "whole", declaration$parameters, declaration$fixed)
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
        manager = manager,
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

# Operation registry (KIT-007). One process-level registry, keyed by operation name, carries every
# per-operation fact the execution layer used to name by hand: the input role and class, the
# parameter validator, resolver and byte limit, the compute handler, the completion and
# publication roles, the step label, the allowed qualification variants and the diagnostic-only
# flag. Entries hold handler *names*, never closures: handlers never travel in requests, and the
# execution layer resolves a name in its own environment at the call site. The registry is filled
# at namespace load (`.onLoad`) and on first use, so the controller and every fresh worker
# `Rscript` resolve the same entries from `request$operation`.

.workflowOperationRegistry <- new.env(parent = emptyenv())

#' Validate one registry entry.
#' @noRd
validateWorkflowOperationEntry <- function(entry) {
    scalar <- function(value) is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
    flag <- function(value) is.logical(value) && length(value) == 1L && !is.na(value)
    valid <- is.list(entry) && scalar(entry$name) && scalar(entry$input) && entry$input %in% c("checkpoint", "peptide") &&
        scalar(entry$input_class) && scalar(entry$compute) && entry$compute %in% c("fixture", "peptide", "protein", "in_process") &&
        # KIT-024: an in-process declaration names no worker handler; every other role must.
        (if (identical(entry$compute, "in_process")) is.null(entry$handler) && !is.null(entry$declaration) else scalar(entry$handler)) &&
        scalar(entry$validate) && (is.null(entry$resolve) || scalar(entry$resolve)) &&
        is.numeric(entry$parameter_limit) && length(entry$parameter_limit) == 1L && entry$parameter_limit > 0 &&
        is.character(entry$parameters) && !anyDuplicated(entry$parameters) &&
        scalar(entry$completion) && entry$completion %in% c("fixture", "peptide", "protein", "da", "in_process") &&
        scalar(entry$publication) && entry$publication %in% c("fixture", "peptide", "protein", "da", "in_process") &&
        (is.null(entry$step_label) || scalar(entry$step_label)) && scalar(entry$omic_type) &&
        is.character(entry$variants) &&
        (is.null(entry$default_variant) || (scalar(entry$default_variant) && entry$default_variant %in% entry$variants)) &&
        flag(entry$diagnostic_only) && flag(entry$inclusion_absence) &&
        (is.null(entry$declaration) || is.list(entry$declaration))
    if (!isTRUE(valid)) cli::cli_abort("workflow operation registry entry is invalid", class = "workflow_operation_invalid")
    if (!is.null(entry$declaration)) {
        validateWorkflowOperation(entry$declaration)
        if (!identical(entry$declaration$name, entry$name) || !identical(entry$declaration$input_class, entry$input_class))
            cli::cli_abort("registry entry {.val {entry$name}} does not match its declaration", class = "workflow_operation_invalid")
    }
    entry
}

#' Register one operation.
#' @param entry Entry list (see `validateWorkflowOperationEntry`); an existing name is replaced only
#'   when `replace = TRUE`.
#' @return The entry, invisibly.
#' @noRd
registerWorkflowOperation <- function(entry, replace = FALSE) {
    entry <- validateWorkflowOperationEntry(entry)
    if (!replace && exists(entry$name, envir = .workflowOperationRegistry, inherits = FALSE))
        cli::cli_abort("operation {.val {entry$name}} is already registered", class = "workflow_operation_invalid")
    assign(entry$name, entry, envir = .workflowOperationRegistry)
    invisible(entry)
}

#' Names of the registered operations (built-ins registered on first use).
#' @noRd
workflowOperationNames <- function() {
    if (!length(ls(.workflowOperationRegistry))) registerBuiltinWorkflowOperations()
    sort(ls(.workflowOperationRegistry))
}

#' Look up one registered operation; an unregistered or malformed name refuses with the
#' execution layer's class.
#' @noRd
workflowOperationEntry <- function(operation) {
    if (!is.character(operation) || length(operation) != 1L || is.na(operation))
        cli::cli_abort("operation adapter is not qualified", class = "workflow_execution_unsupported")
    if (!length(ls(.workflowOperationRegistry))) registerBuiltinWorkflowOperations()
    entry <- get0(operation, envir = .workflowOperationRegistry, inherits = FALSE)
    if (is.null(entry)) cli::cli_abort("operation adapter is not qualified", class = "workflow_execution_unsupported")
    entry
}

#' Register the built-in operations: the diagnostic fixture, the two DIA declarations of KIT-006, the
#' imputation declaration of KIT-009, the protein replicate filter of KIT-010, the four remaining peptide
#' QC declarations of KIT-014 and the four protein QC declarations of KIT-015.
#' `default_variant` (KIT-008) is the variant an unqualified launch runs, in the worker and in process.
#' Called from `.onLoad` and from the first lookup in a process that bypassed it.
#' @noRd
registerBuiltinWorkflowOperations <- function() {
    rm(list = ls(.workflowOperationRegistry), envir = .workflowOperationRegistry)
    registerWorkflowOperation(list(name = "wfx_fixture", input = "checkpoint", input_class = "PeptideQuantitativeData",
        compute = "fixture", handler = "computeWorkflowFixture", validate = "validateWorkflowFixtureParameters", resolve = NULL,
        parameter_limit = 65536, parameters = c("mode", "seconds", "bytes"), completion = "fixture", publication = "fixture",
        step_label = NULL, omic_type = "proteomics", variants = character(), default_variant = NULL,
        diagnostic_only = TRUE, inclusion_absence = FALSE, declaration = NULL))
    registerWorkflowOperation(list(name = "precursor_rollup", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowPrecursorRollupAction", validate = "validateWorkflowPrecursorRollupParameters",
        resolve = NULL, parameter_limit = 65536, parameters = character(), completion = "peptide", publication = "peptide",
        step_label = "4_precursor_rollup", omic_type = "proteomics", variants = "p14", default_variant = "p14",
        diagnostic_only = FALSE, inclusion_absence = FALSE, declaration = diaPrecursorRollupOperation()))
    # The frozen 4,000-sample filter includes duplicate/whitespace-bearing inclusion arguments: its
    # parameter bytes keep the larger finite cap.
    registerWorkflowOperation(list(name = "minimum_peptides_per_sample", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowSampleFilterAction", validate = "validateWorkflowSampleFilterParameters",
        resolve = "resolveWorkflowSampleFilterParameters", parameter_limit = 262144,
        parameters = c("peptides_per_sample_cutoff", "inclusion_list"), completion = "peptide", publication = "peptide",
        step_label = "6_sample_filtered", omic_type = "proteomics", variants = c("native", "physical", "selection"),
        # KIT-008: production runs the KIT-006 mechanism (adaptive complete-sample batches); the native R
        # control and the selection arm stay reachable through an explicit adapters qualification.
        default_variant = "physical", diagnostic_only = FALSE, inclusion_absence = TRUE,
        declaration = diaSampleFilterOperation()))
    # KIT-009: the imputation step on the feature-group partitioner; `grouped` is the per-group loop as it
    # stood before KIT-009 (the control arm of the measurement), `partitioned` the production variant.
    registerWorkflowOperation(list(name = "peptide_imputation", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowPeptideImputationAction", validate = "validateWorkflowPeptideImputationParameters",
        resolve = NULL, parameter_limit = 65536, parameters = c("proportion_missing_values", "exclusion_column"),
        completion = "peptide", publication = "peptide", step_label = "8_imputed", omic_type = "proteomics",
        variants = c("grouped", "partitioned"), default_variant = "partitioned", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaPeptideImputationOperation()))
    # KIT-014: the remaining peptide QC steps. The three feature-group sections share the family's variants
    # (`grouped`: the per-group loop before KIT-009, the adapters-only control arm; `partitioned`: production);
    # the Q-value step has one path, `streamed` (the evidence pre-pass over whole protein groups, then the
    # ordered row stream).
    registerWorkflowOperation(list(name = "peptide_qvalue_filter", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowPeptideQvalueAction", validate = "validateWorkflowPeptideQvalueParameters",
        resolve = NULL, parameter_limit = 65536,
        parameters = c("qvalue_threshold", "global_qvalue_threshold", "global_pg_qvalue_threshold", "proteotypic_only",
            "confidence_provenance_mode"),
        completion = "peptide", publication = "peptide", step_label = "2_qval_filtered", omic_type = "proteomics",
        variants = "streamed", default_variant = "streamed", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaQvalueFilterOperation()))
    registerWorkflowOperation(list(name = "peptide_protein_evidence_filter", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowPeptideProteinEvidenceAction",
        validate = "validateWorkflowPeptideProteinEvidenceParameters", resolve = NULL, parameter_limit = 65536,
        parameters = c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"),
        completion = "peptide", publication = "peptide", step_label = "3_protein_peptide_filtered", omic_type = "proteomics",
        variants = c("grouped", "partitioned"), default_variant = "partitioned", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaPeptideProteinEvidenceOperation()))
    registerWorkflowOperation(list(name = "peptide_intensity_filter", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowPeptideIntensityAction", validate = "validateWorkflowPeptideIntensityParameters",
        resolve = NULL, parameter_limit = 65536,
        parameters = c("strict_mode", "min_reps_per_group", "min_groups", "peptides_intensity_cutoff_percentile", "grouping_variable"),
        completion = "peptide", publication = "peptide", step_label = "5_intensity_filtered", omic_type = "proteomics",
        variants = c("grouped", "partitioned"), default_variant = "partitioned", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaPeptideIntensityFilterOperation()))
    registerWorkflowOperation(list(name = "peptide_replicate_filter", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "peptide", handler = "computeWorkflowPeptideReplicateAction", validate = "validateWorkflowPeptideReplicateParameters",
        resolve = NULL, parameter_limit = 65536, parameters = "replicate_group_column",
        completion = "peptide", publication = "peptide", step_label = "7_replicate_filtered", omic_type = "proteomics",
        variants = c("grouped", "partitioned"), default_variant = "partitioned", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaPeptideReplicateFilterOperation()))
    # KIT-015: the protein QC steps, each placed whole (2,500 proteins wide); `whole` is their one variant.
    registerWorkflowOperation(list(name = "protein_iq_rollup", input = "peptide", input_class = "PeptideQuantitativeData",
        compute = "protein", handler = "computeWorkflowProteinIQRollupAction", validate = "validateWorkflowProteinIQRollupParameters",
        resolve = NULL, parameter_limit = 65536, parameters = character(),
        completion = "protein", publication = "protein", step_label = "9_protein_s4_created", omic_type = "proteomics",
        variants = "whole", default_variant = "whole", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaProteinIQRollupOperation()))
    registerWorkflowOperation(list(name = "protein_accession_cleanup", input = "checkpoint", input_class = "ProteinQuantitativeData",
        compute = "protein", handler = "computeWorkflowProteinAccessionAction", validate = "validateWorkflowProteinAccessionParameters",
        resolve = NULL, parameter_limit = 65536, parameters = c("delimiter", "aggregation_method"),
        completion = "protein", publication = "protein", step_label = "10_protein_accession_cleaned", omic_type = "proteomics",
        variants = "whole", default_variant = "whole", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaProteinAccessionCleanupOperation()))
    registerWorkflowOperation(list(name = "protein_intensity_filter", input = "checkpoint", input_class = "ProteinQuantitativeData",
        compute = "protein", handler = "computeWorkflowProteinIntensityAction", validate = "validateWorkflowProteinIntensityParameters",
        resolve = NULL, parameter_limit = 65536,
        parameters = c("strict_mode", "min_reps_per_group", "min_groups", "proteins_intensity_cutoff_percentile"),
        completion = "protein", publication = "protein", step_label = "11_protein_intensity_filtered", omic_type = "proteomics",
        variants = "whole", default_variant = "whole", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaProteinIntensityFilterOperation()))
    registerWorkflowOperation(list(name = "protein_duplicate_removal", input = "checkpoint", input_class = "ProteinQuantitativeData",
        compute = "protein", handler = "computeWorkflowProteinDuplicateAction", validate = "validateWorkflowProteinDuplicateParameters",
        resolve = NULL, parameter_limit = 65536, parameters = "aggregation_method",
        completion = "protein", publication = "protein", step_label = "12_duplicates_removed", omic_type = "proteomics",
        variants = "whole", default_variant = "whole", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaProteinDuplicateRemovalOperation()))
    # KIT-010: the protein replicate filter on the row partitioner; `whole` is the GUI's in-process path (the
    # object hydrated, the kernel once; the control arm of the measurement), `rows` the production variant.
    registerWorkflowOperation(list(name = "protein_replicate_filter", input = "checkpoint", input_class = "ProteinQuantitativeData",
        compute = "protein", handler = "computeWorkflowProteinReplicateAction", validate = "validateWorkflowProteinReplicateParameters",
        resolve = NULL, parameter_limit = 65536, parameters = c("grouping_variable", "parallel_cores"),
        completion = "protein", publication = "protein", step_label = "13_protein_replicate_filtered", omic_type = "proteomics",
        variants = c("whole", "rows"), default_variant = "rows", diagnostic_only = FALSE,
        inclusion_absence = FALSE, declaration = diaProteinReplicateFilterOperation()))
    # KIT-019: the DA publication phase (the annotated exports, the static volcano and NumSig figures and the
    # interactive Glimma pages) as a worker operation on the accepted DA manifest; the DA kernel (limma) stays in
    # process. Completion records nothing new (the manifest is already committed); publication moves the worker's
    # staged files into the project and saves the pages as the KIT-016 rendered assets.
    registerWorkflowOperation(list(name = "da_publication", input = "checkpoint", input_class = "ProteinQuantitativeData",
        compute = "protein", handler = "computeWorkflowDaPublicationAction", validate = "validateWorkflowDaPublicationParameters",
        resolve = NULL, parameter_limit = 65536, parameters = "annotations",
        completion = "da", publication = "da", step_label = "14_da_published", omic_type = "proteomics",
        variants = "whole", default_variant = "whole", diagnostic_only = FALSE, inclusion_absence = FALSE))
    # KIT-024: the whole-matrix declarations that run in process only (KIT-004, KIT-016) are registered too, so the
    # registry is the one place an operation is named (KIT-007); a dispatch of one to the worker is refused.
    for (declaration in list(diaNormalizationOperation(), diaPearsonPairsOperation(), diaPearsonSelectionOperation(),
        diaContrastOperation(), diaEnrichmentOperation())) {
        registerWorkflowOperation(list(name = declaration$name, input = "checkpoint", input_class = declaration$input_class,
            compute = "in_process", handler = NULL, validate = "validateWorkflowInProcessParameters", resolve = NULL,
            parameter_limit = 65536, parameters = declaration$parameters %||% character(),
            completion = "in_process", publication = "in_process", step_label = NULL, omic_type = "proteomics",
            variants = "whole", default_variant = "whole", diagnostic_only = FALSE, inclusion_absence = FALSE,
            declaration = declaration))
    }
    invisible(workflowOperationNames())
}

#' Refuse to dispatch an in-process declaration to the worker (KIT-024).
#' @noRd
assertWorkflowOperationDispatchable <- function(operation) {
    entry <- workflowActionSpec(operation)
    if (identical(entry$compute, "in_process")) cli::cli_abort(
        "operation {.val {operation}} runs in process only (whole-matrix declaration); it cannot be dispatched to a worker",
        class = "workflow_execution_unsupported")
    invisible(entry)
}
