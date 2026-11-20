# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# KIT-062: the integration owner. An integration project has its own on-disk state owner (store workflow id
# `INTEGRATION`); its first state is the omics' integration inputs as read by `readIntegrationInput()` (KIT-061), held
# in one declared class, `IntegrationViewsData`, with each input's manifest sha256 and path. The import is a replay
# precondition (as an omic's import is), each manifest an input role `integration_input_<omic>` checked by sha256. The
# model (MOFA2) is a product-only declaration on that state, accepted on the head as the transcriptomics RUVIII export
# is and listed in the analysis record under `analyses`.

INTEGRATION_WORKFLOW_ID <- "INTEGRATION"
INTEGRATION_INPUTS_STATE <- "integration_inputs"
INTEGRATION_IMPORT_OPERATION <- "integration_import"
INTEGRATION_INPUT_ROLE_PREFIX <- "integration_input_"

#' The integration owner's inputs state (KIT-062): every view of every imported omic, its feature table, one design
#' row per run (`omic`, `run`, `group`, `biological_sample`), and per omic the manifest, its sha256 and path, the
#' input's samples table and the column that names the biological sample.
#' @slot views Named list by view of tables: `feature_id`, then one numeric column per run of the view's omic.
#' @slot features Named list by view of the input's feature tables (`feature_id`, `display_name`, annotation).
#' @slot design_matrix One row per run of every omic.
#' @slot inputs Named list by omic: `omic`, `path`, `manifest_sha256`, `views`, `manifest`, `samples`,
#'   `biological_sample_column`.
#' @slot sample_id,group_id,biological_sample_id The design's run, group and biological-sample columns.
#' @slot args Unused; the declared classes' argument slot.
#' @noRd
setClass("IntegrationViewsData",
    slots = c(views = "list", features = "list", design_matrix = "data.frame", inputs = "list",
        sample_id = "character", group_id = "character", biological_sample_id = "character", args = "list"),
    prototype = list(views = list(), features = list(), design_matrix = data.frame(), inputs = list(),
        sample_id = "run", group_id = "group", biological_sample_id = "biological_sample", args = list()),
    validity = function(object) {
        if (!identical(names(object@views), names(object@features))) return("views and features must name the same views")
        roles <- list(object@sample_id, object@group_id, object@biological_sample_id)
        if (!all(vapply(roles, function(role) is.character(role) && length(role) == 1L && !is.na(role), logical(1)))) {
            return("the run, group and biological-sample columns are one name each")
        }
        TRUE
    })

#' Refuse an integration step with a classed condition in plain words (KIT-062).
#' @param message What failed and why, lower case.
#' @param class The refusal's own class; every refusal is also `integration_refused`.
#' @noRd
integrationRefuse <- function(message, class) {
    cli::cli_abort(message, class = c(class, "integration_refused"), call = NULL)
}

#' A new integration owner in a project folder (KIT-062): the on-disk state owner with the store workflow id
#' `INTEGRATION`, reference states on, its omic `integration` (the record's and the replay's).
#' @param root The integration project's folder (created when absent); the store is `<root>/.workflow-entities`.
#' @param project_id The store's project id.
#' @return The owner (`WorkflowStateOwner`).
#' @noRd
newIntegrationStateOwner <- function(root, project_id = "integration") {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    manager <- newWorkflowStateOwner(newWorkflowEntityStore(root, project_id, INTEGRATION_WORKFLOW_ID),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$.context$omic_type <- "integration"
    manager
}

#' Read and check the omics' integration inputs (KIT-062; KIT-147: no matrix is loaded). Each folder is read and
#' validated by `readIntegrationInput()`; the omics are ordered by name, their views in manifest order. Refused in plain
#' words: fewer than two omics (D7), two inputs of one omic, a view name used by two omics (INT-4 names views
#' `<omic>_<assay>`), a biological-sample column the samples table lacks or leaves empty.
#' @param dirs The omics' folders, `<project>/integration/<omic>/`.
#' @return `list(inputs, views, features, design, record, manifests)`: the readers' values by omic, per view its
#'   folder and omic, the feature tables by view, one design row per run, the `inputs` slot's entries by omic, and the
#'   manifests' paths by omic.
#' @noRd
integrationReadInputs <- function(dirs) {
    dirs <- unique(as.character(dirs))
    if (length(dirs) < 2L) {
        integrationRefuse(sprintf("integration needs the inputs of two or more omics; %d given", length(dirs)),
            "integration_too_few_omics")
    }
    inputs <- lapply(dirs, readIntegrationInput)
    omics <- vapply(inputs, function(input) as.character(input$manifest$omic_type), character(1))
    if (anyDuplicated(omics)) {
        integrationRefuse(sprintf("two integration inputs are of the same omic (%s); give one per omic",
            paste(unique(omics[duplicated(omics)]), collapse = ", ")), "integration_input_duplicate")
    }
    order <- order(omics, method = "radix")
    inputs <- inputs[order]
    dirs <- dirs[order]
    omics <- omics[order]
    views <- list()
    features <- list()
    design <- list()
    record <- list()
    for (i in seq_along(inputs)) {
        input <- inputs[[i]]
        omic <- omics[[i]]
        manifest <- input$manifest
        samples <- input$samples
        sample_column <- manifest$sample_id_column
        technical <- manifest$technical_replicate_column
        if (!is.character(technical) || length(technical) != 1L) technical <- NA_character_
        biological <- manifest$biological_sample_column %||%
            integrationInputBiologicalSampleColumn(samples, technical, sample_column)
        if (!is.character(biological) || length(biological) != 1L || !biological %in% names(samples)) {
            integrationRefuse(sprintf("the integration input of %s names no biological-sample column of its samples table",
                omic), "integration_input_invalid")
        }
        key <- as.character(samples[[biological]])
        if (anyNA(key) || any(!nzchar(key))) {
            integrationRefuse(sprintf("the integration input of %s leaves a run without its biological sample in column %s",
                omic, biological), "integration_input_invalid")
        }
        for (view in names(input$views)) {
            if (view %in% names(views)) {
                integrationRefuse(sprintf("the view name %s is used by two omics; views are named <omic>_<assay>", view),
                    "integration_input_duplicate")
            }
            views[[view]] <- list(dir = dirs[[i]], omic = omic)
            features[[view]] <- input$views[[view]]$features
        }
        design[[i]] <- data.frame(omic = omic, run = as.character(samples[[sample_column]]),
            group = as.character(samples[[manifest$group_column]]), biological_sample = key, stringsAsFactors = FALSE)
        path <- file.path(dirs[[i]], INTEGRATION_INPUT_MANIFEST)
        record[[omic]] <- list(omic = omic, path = normalizePath(path, mustWork = TRUE),
            manifest_sha256 = digest::digest(file = path, algo = "sha256"), views = names(input$views),
            manifest = manifest, samples = samples, biological_sample_column = biological)
    }
    design <- do.call(rbind, design)
    rownames(design) <- NULL
    list(inputs = stats::setNames(inputs, omics), views = views, features = features, design = design, record = record,
        manifests = stats::setNames(vapply(record, function(entry) entry$path, character(1)), omics))
}

#' The reader of each view of the read inputs (`openIntegrationInputView`), in the inputs' view order.
#' @param read `integrationReadInputs()`'s value.
#' @noRd
integrationInputViewReaders <- function(read) {
    readers <- lapply(names(read$views), function(view) {
        entry <- read$views[[view]]
        openIntegrationInputView(entry$dir, view, read$inputs[[entry$omic]])
    })
    stats::setNames(readers, names(read$views))
}

#' Read the omics' integration inputs into one `IntegrationViewsData` in memory (KIT-062): `integrationReadInputs()`
#' with every view's matrix read whole. For small inputs and tests; the import streams the matrices into the store
#' instead (`importIntegrationInputs`).
#' @param dirs The omics' folders, `<project>/integration/<omic>/`.
#' @return `list(object, manifests)`: the object, and the manifests' paths named by omic.
#' @noRd
integrationViewsFromInputs <- function(dirs) {
    read <- integrationReadInputs(dirs)
    matrices <- lapply(integrationInputViewReaders(read), function(reader) reader$read())
    object <- methods::new("IntegrationViewsData", views = matrices, features = read$features, design_matrix = read$design,
        inputs = read$record)
    list(object = object, manifests = read$manifests)
}

#' Save one data set's matrix in the integration store, streamed row group by row group from its Parquet file (KIT-147):
#' the saved table of the same file (by its sha256) is reused; otherwise each row group is appended in slices of at
#' most 32 MiB of values (the table writer's batch bound), so no more than one row group is held at a time.
#' @param context The integration owner's context.
#' @param reader `openIntegrationInputView()`'s value.
#' @return The saved table's reference.
#' @noRd
integrationImportViewTable <- function(context, reader) {
    identity <- newSavedWorkflowScientificIdentity(reader$sha256, "integration-input-view",
        digest::digest(list(view = reader$view, columns = reader$columns), algo = "sha256"), INTEGRATION_INPUT_SCHEMA,
        context$software)
    hit <- resolveSavedWorkflowEntity(context$store, "result_table", identity, hydrate = FALSE)
    if (identical(hit$status, "hit")) return(hit$reference)
    slice <- max(1L, as.integer(floor(33554432 / (8 * length(reader$columns)))))
    saveWorkflowPartitionedTable(context, reader$prototype, identity, function(append) {
        for (index in seq_len(reader$row_groups)) {
            page <- reader$read_row_group(index)
            rows <- nrow(page)
            if (rows <= slice) {
                append(page)
            } else {
                for (start in seq.int(1L, rows, by = slice)) append(page[seq.int(start, min(rows, start + slice - 1L)), , drop = FALSE])
            }
            page <- NULL
        }
        invisible(NULL)
    })
}

#' Import the omics' integration inputs into the integration owner (KIT-062): the inputs state committed with the
#' provenance `integration_import` (a replay precondition), and each manifest recorded as the input role
#' `integration_input_<omic>` with its sha256. KIT-147: each data set is streamed from its Parquet file into the store
#' as a table (`integrationImportViewTable`); the state is checkpointed from zero-row prototypes of the data sets, the
#' feature tables and the design (`omicsRunPartitioned`'s pattern), its `views` components replaced by the saved tables,
#' and committed. No matrix is held whole.
#' @param manager The integration owner (`newIntegrationStateOwner`).
#' @param dirs The omics' folders, `<project>/integration/<omic>/`.
#' @param description The state's description.
#' @return The inputs state's reference, invisibly.
#' @noRd
importIntegrationInputs <- function(manager, dirs, description = "Integration inputs imported") {
    if (!usesWorkflowReferenceStates(manager)) {
        cli::cli_abort("the integration import needs the integration owner", class = "saved_entity_scope")
    }
    read <- integrationReadInputs(dirs)
    config <- list(outcomes = list(omics = names(read$record), views = names(read$views), runs = nrow(read$design),
        biological_samples = length(unique(read$design$biological_sample))))
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        readers <- integrationInputViewReaders(read)
        tables <- lapply(readers, function(reader) integrationImportViewTable(context, reader))
        object <- methods::new("IntegrationViewsData", views = lapply(readers, function(reader) reader$prototype),
            features = read$features, design_matrix = read$design, inputs = read$record)
        setWorkflowStepProvenance(manager, INTEGRATION_IMPORT_OPERATION)
        checkpoint <- checkpointWorkflowExecutionResult(manager, object)
        checkpoint$components$views <- tables
        checkpoint$dimensions$rows <- sum(vapply(tables, function(table) as.integer(table$metadata$shape[[1L]]), integer(1)))
        commitWorkflowReferenceState(manager, INTEGRATION_INPUTS_STATE, checkpoint, config, description)
        recordWorkflowInputs(manager, paths = stats::setNames(read$manifests,
            paste0(INTEGRATION_INPUT_ROLE_PREFIX, names(read$manifests))))
    })
    invisible(getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE))
}

# ---------------------------------------------------------------------------------------------------------------
# The model: MOFA2 on the inputs state, a product-only declaration (KIT-062 item 3).
# ---------------------------------------------------------------------------------------------------------------

INTEGRATION_MODEL_OPERATION <- "integration_mofa_model"
# Decision D7 (KIT-060): below this many biological samples measured in every view the model warns; it does not refuse.
INTEGRATION_FEW_SAMPLES <- 16L

#' The model's method parameters with their provisional defaults: the one table the user's integration workbook
#' settles (`settled_by`). The provisional values are the ones the MOFA2 smoke lane trained with on the tier-1 corpus
#' (run 44423980807: MOFA2 1.22.0, mofapy2 0.7.3; every other MOFA2 option at its default), and the decisions of
#' 2026-11-19 (INT-2, INT-3) and KIT-060 D4.
#' @return Named list by parameter: `default`, `choices` (NULL for a number or a flag), `settled_by`.
#' @noRd
integrationMofaParameterSpec <- function() {
    workbook <- "workbook (provisional)"
    list(
        technical_replicates = list(default = "average", choices = "average",
            settled_by = "INT-2 (averaged per biological sample); the workbook may override"),
        internal_standards = list(default = "drop", choices = c("drop", "keep"),
            settled_by = "INT-3 (dropped); the workbook may override"),
        missing_samples = list(default = "keep_missing", choices = c("keep_missing", "drop_incomplete"),
            settled_by = "KIT-060 D4 (kept as missing in that view); the workbook may override"),
        feature_selection = list(default = "none", choices = c("none", "top_variance"), settled_by = workbook),
        features_per_view = list(default = NA_integer_, choices = NULL, settled_by = workbook),
        scaling = list(default = "none", choices = c("none", "feature_zscore"), settled_by = workbook),
        scale_views = list(default = FALSE, choices = NULL, settled_by = workbook),
        center_groups = list(default = TRUE, choices = NULL, settled_by = workbook),
        num_factors = list(default = 5L, choices = NULL, settled_by = workbook),
        convergence_mode = list(default = "fast", choices = c("fast", "medium", "slow"), settled_by = workbook),
        maxiter = list(default = 1000L, choices = NULL, settled_by = workbook),
        seed = list(default = 20261119L, choices = NULL, settled_by = workbook),
        # INT-6 (user, 2026-11-19): decided, not provisional. Double precision gives the same model on every computer
        # (KIT-062 diagnosis run 44462198417: MOFA2's float32 default followed one of two training paths depending on
        # the CPU model and thread count); single precision uses less memory for very large studies, and its results
        # can differ slightly between computers.
        precision = list(default = "double", choices = c("double", "single"),
            settled_by = "INT-6 (decided 2026-11-19: double, the same model on every computer)"))
}

#' The provisional-parameter table (one row per method parameter), for the ticket, the ledger and the GUI.
#' @return `data.frame(parameter, default, choices, settled_by)`.
#' @noRd
integrationMofaParameterTable <- function() {
    spec <- integrationMofaParameterSpec()
    data.frame(parameter = names(spec),
        default = vapply(spec, function(entry) paste(format(entry$default), collapse = ""), character(1)),
        choices = vapply(spec, function(entry) if (is.null(entry$choices)) "" else paste(entry$choices, collapse = ", "),
            character(1)),
        settled_by = vapply(spec, function(entry) entry$settled_by, character(1)),
        stringsAsFactors = FALSE, row.names = NULL)
}

#' The model's declared parameters in identity order: the inputs (each omic's manifest sha256, the views), then the
#' method parameters.
#' @noRd
integrationMofaParameterNames <- function() c("input_manifests", "views", names(integrationMofaParameterSpec()))

#' The manifests' sha256 the inputs state holds, named by omic in the state's order.
#' @param object An `IntegrationViewsData` with its `inputs` hydrated.
#' @noRd
integrationStateManifests <- function(object) {
    vapply(object@inputs, function(entry) as.character(entry$manifest_sha256), character(1))
}

#' Resolve and check the model's parameters against the inputs state (KIT-062): each absent method parameter takes its
#' provisional default, the inputs default to the state's (every view, every manifest), values are put in one form
#' (whole numbers as integers) so that the identity does not depend on how a caller wrote them. Refused in plain
#' words: a value outside its choices, manifests that are not the state's (`integration_input_stale`), an unknown
#' view, fewer than two views (`integration_too_few_views`).
#' @param object An `IntegrationViewsData` with its `inputs` hydrated (the views' names are enough).
#' @param parameters The caller's parameters (named list).
#' @return The declared parameters, named in identity order.
#' @noRd
integrationResolveMofaParameters <- function(object, parameters = list()) {
    spec <- integrationMofaParameterSpec()
    invalid <- function(message) integrationRefuse(message, "integration_parameters_invalid")
    given <- function(name) parameters[[name, exact = TRUE]]
    out <- list(input_manifests = given("input_manifests") %||% integrationStateManifests(object),
        views = given("views") %||% names(object@views))
    for (name in names(spec)) {
        value <- given(name)
        out[name] <- list(if (is.null(value)) spec[[name]]$default else value)
    }
    for (name in names(spec)) {
        value <- out[[name]]
        choices <- spec[[name]]$choices
        if (!is.null(choices) && !(is.character(value) && length(value) == 1L && isTRUE(value %in% choices))) {
            invalid(sprintf("%s must be one of %s", name, paste(choices, collapse = ", ")))
        }
    }
    for (name in c("scale_views", "center_groups")) {
        if (!(is.logical(out[[name]]) && length(out[[name]]) == 1L && !is.na(out[[name]]))) {
            invalid(sprintf("%s must be TRUE or FALSE", name))
        }
    }
    whole <- function(name, minimum) {
        value <- out[[name]]
        if (!(is.numeric(value) && length(value) == 1L && !is.na(value) && value == round(value) && value >= minimum)) {
            invalid(sprintf("%s must be a whole number of at least %d", name, minimum))
        }
        as.integer(value)
    }
    out$num_factors <- whole("num_factors", 1L)
    out$maxiter <- whole("maxiter", 1L)
    out$seed <- whole("seed", 0L)
    if (identical(out$feature_selection, "top_variance")) {
        out$features_per_view <- whole("features_per_view", 1L)
    } else if (!(length(out$features_per_view) == 1L && is.na(out$features_per_view))) {
        invalid("features_per_view is set only with feature_selection = top_variance")
    } else {
        out$features_per_view <- NA_integer_
    }
    manifests <- out$input_manifests
    recorded <- integrationStateManifests(object)
    if (!is.character(manifests) || is.null(names(manifests)) || !setequal(names(manifests), names(recorded)) ||
        !identical(unname(manifests[names(recorded)]), unname(recorded))) {
        integrationRefuse("the integration inputs named for the model are not the ones imported (a manifest's sha256 differs); import the inputs again",
            "integration_input_stale")
    }
    out$input_manifests <- recorded
    views <- out$views
    if (!is.character(views) || anyNA(views) || anyDuplicated(views) || !all(views %in% names(object@views))) {
        invalid(sprintf("the views must be distinct views of the inputs state (%s)", paste(names(object@views), collapse = ", ")))
    }
    if (length(views) < 2L) {
        integrationRefuse(sprintf("the model needs two or more views; %d selected", length(views)), "integration_too_few_views")
    }
    out <- out[integrationMofaParameterNames()]
    # KIT-063 item A (KIT-139 follow-up): an optional user-edited sample match. Absent, or giving no run another
    # biological sample, it adds nothing, so the parameters and the model's identity are those of before.
    table <- given("sample_match_table")
    declared <- given("sample_match")
    if (!is.null(table)) {
        match <- integrationResolveSampleMatch(object, table)
        if (!is.null(match$table)) {
            if (!is.null(declared) && !identical(declared, match$sha256)) {
                invalid("the sample match table does not have the checksum recorded with it")
            }
            out$sample_match <- match$sha256
            out$sample_match_table <- match$table
        }
    } else if (!is.null(declared)) {
        invalid("a sample match checksum was given without its table")
    }
    out
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-063 item A (KIT-139 follow-up): the user-edited sample match. Samples are matched across omics by the
# biological-sample column each omic saved; a user whose omics name one biological sample differently uploads a match
# table (per omic, a run or a saved biological-sample name, and the biological sample it belongs to). Only the runs it
# gives another biological sample enter the model, as the parameter `sample_match` (the sha256 of that canonical table,
# identity) with the table beside it (`sample_match_table`, not identity, checked against the checksum).
# ---------------------------------------------------------------------------------------------------------------

#' The sample match as it stands (KIT-063 item A): one row per run of every imported omic, the biological sample its
#' omic saved and the one it is matched to (the saved one, or the match's).
#' @param object An `IntegrationViewsData` with its `inputs` hydrated.
#' @param match NULL, or a canonical match (`integrationResolveSampleMatch()$table`).
#' @return `data.frame(omic, run, saved_biological_sample, biological_sample)`.
#' @noRd
integrationSampleMatchTemplate <- function(object, match = NULL) {
    rows <- lapply(object@inputs, function(entry) {
        samples <- entry$samples
        data.frame(omic = entry$omic, run = as.character(samples[[entry$manifest$sample_id_column]]),
            saved_biological_sample = as.character(samples[[entry$biological_sample_column]]), stringsAsFactors = FALSE)
    })
    table <- do.call(rbind, rows)
    rownames(table) <- NULL
    table$biological_sample <- table$saved_biological_sample
    if (!is.null(match)) {
        hit <- match(paste(table$omic, table$run, sep = "\r"), paste(match$omic, match$run, sep = "\r"))
        table$biological_sample[!is.na(hit)] <- match$biological_sample[hit[!is.na(hit)]]
    }
    table
}

#' Check a user's sample match table against the imported omics and reduce it to its canonical form (KIT-063 item A).
#' Each row names an omic, then a run (`run`) or, with `run` empty, a saved biological-sample name
#' (`saved_biological_sample`, every run of that name), and the biological sample (`biological_sample`). Refused
#' (`integration_sample_match_invalid`): a table without those columns, a known run or name left without a biological
#' sample, a run mapped more than once. Reported, not refused: rows naming an omic, run or saved name the import does
#' not have (counted, ignored). A biological sample given two groups is refused later by the alignment, as without a
#' table.
#' @param object An `IntegrationViewsData` with its `inputs` hydrated.
#' @param table A data frame, or a list of equal-length columns (the canonical form is accepted as it is).
#' @return `list(table, sha256, report)`: the runs given another biological sample than their omic saved, as
#'   `list(omic, run, biological_sample)` sorted by omic and run (NULL when none), its sha256 (NULL when none), and the
#'   counts (`rows`, `runs_mapped`, `runs_changed`, `unknown_omics`, `unknown_runs`, `unknown_names`, `rows_empty`).
#' @noRd
integrationResolveSampleMatch <- function(object, table) {
    invalid <- function(message) integrationRefuse(message, "integration_sample_match_invalid")
    if (is.list(table) && !is.data.frame(table)) {
        table <- tryCatch(as.data.frame(table, stringsAsFactors = FALSE, optional = TRUE),
            error = function(error) invalid("the sample match table is not a table of equal-length columns"))
    }
    if (!is.data.frame(table)) invalid("the sample match table is not a table")
    names(table) <- trimws(names(table))
    if (!all(c("omic", "biological_sample") %in% names(table)) || !any(c("run", "saved_biological_sample") %in% names(table))) {
        invalid("the sample match table needs the columns omic, run (or saved_biological_sample) and biological_sample")
    }
    column <- function(name) {
        value <- table[[name, exact = TRUE]]
        if (is.null(value)) return(rep("", nrow(table)))
        value <- trimws(as.character(value))
        value[is.na(value)] <- ""
        value
    }
    omic <- column("omic")
    run <- column("run")
    saved <- column("saved_biological_sample")
    target <- column("biological_sample")
    template <- integrationSampleMatchTemplate(object)
    run_key <- paste(template$omic, template$run, sep = "\r")
    saved_key <- paste(template$omic, template$saved_biological_sample, sep = "\r")
    known <- omic %in% template$omic
    by_run <- nzchar(run)
    by_name <- !by_run & nzchar(saved)
    unknown_run <- known & by_run & !paste(omic, run, sep = "\r") %in% run_key
    unknown_name <- known & by_name & !paste(omic, saved, sep = "\r") %in% saved_key
    usable <- known & ((by_run & !unknown_run) | (by_name & !unknown_name))
    empty <- usable & !nzchar(target)
    if (any(empty)) {
        # Worded without "sample (": the registry test's scan for random-number calls reads message text too.
        invalid(sprintf("%d row(s) of the sample match table give no biological sample; the first is %s, %s", sum(empty),
            omic[empty][[1L]], if (nzchar(run[empty][[1L]])) run[empty][[1L]] else saved[empty][[1L]]))
    }
    rows <- which(usable)
    hits <- lapply(rows, function(i) if (by_run[[i]]) which(run_key == paste(omic[[i]], run[[i]], sep = "\r")) else
        which(saved_key == paste(omic[[i]], saved[[i]], sep = "\r")))
    mapped <- unlist(hits)
    mapped_target <- rep(target[rows], lengths(hits))
    twice <- unique(mapped[duplicated(mapped)])
    if (length(twice)) {
        invalid(sprintf("%d run(s) are mapped more than once in the sample match table (first: %s, %s); map each run once",
            length(twice), template$omic[[twice[[1L]]]], template$run[[twice[[1L]]]]))
    }
    changed <- mapped_target != template$saved_biological_sample[mapped]
    sorted <- order(template$omic[mapped][changed], template$run[mapped][changed], method = "radix")
    canonical <- if (any(changed)) list(omic = template$omic[mapped][changed][sorted],
        run = template$run[mapped][changed][sorted], biological_sample = mapped_target[changed][sorted])
    list(table = canonical, sha256 = if (!is.null(canonical)) digest::digest(canonical, algo = "sha256"),
        report = list(rows = nrow(table), runs_mapped = length(mapped), runs_changed = sum(changed),
            unknown_omics = sum(!known & (nzchar(omic) | by_run | by_name)), unknown_runs = sum(unknown_run),
            unknown_names = sum(unknown_name), rows_empty = sum(known & !by_run & !by_name)))
}

#' The design with the sample match applied (KIT-063 item A): each matched run's biological sample replaced.
#' @param design The inputs state's design (`omic`, `run`, `group`, `biological_sample`).
#' @param match NULL (the design as it is), or the canonical match.
#' @noRd
integrationApplySampleMatch <- function(design, match) {
    if (is.null(match)) return(design)
    hit <- match(paste(design$omic, design$run, sep = "\r"), paste(match$omic, match$run, sep = "\r"))
    design$biological_sample[!is.na(hit)] <- as.character(match$biological_sample[hit[!is.na(hit)]])
    design
}

#' Keep the user's sample match as a table of the integration store (KIT-148): its runs given another biological sample
#' (`omic`, `run`, `biological_sample`), the canonical match whose sha256 is the model's parameter `sample_match`. A
#' request to the managed worker carries this table's reference, not the table (a match of thousands of runs is beyond
#' the request's 64 KiB); the model's product depends on it. The same match is the same table.
#' @param manager The integration owner.
#' @param match The canonical match (`integrationResolveSampleMatch()$table`).
#' @return The table's reference.
#' @noRd
integrationStoreSampleMatch <- function(manager, match) {
    withWorkflowStateScope(manager, function() saveWorkflowOwnedTable(manager$.context, data.frame(omic = as.character(match$omic),
        run = as.character(match$run), biological_sample = as.character(match$biological_sample), stringsAsFactors = FALSE)))
}

#' The canonical sample match kept in the store (KIT-148; `integrationStoreSampleMatch`).
#' @param manager The owner of the store.
#' @param reference The table's reference.
#' @return `list(omic, run, biological_sample)`; its sha256 is checked against `sample_match` when the parameters resolve.
#' @noRd
integrationStoredSampleMatch <- function(manager, reference) {
    table <- withWorkflowStateScope(manager, function() materializeWorkflowTable(manager$.context, reference))
    list(omic = as.character(table$omic), run = as.character(table$run), biological_sample = as.character(table$biological_sample))
}

#' Refuse the model when an imported input has changed on disk since the import (KIT-062): each recorded manifest
#' must still be at its path with the recorded sha256.
#' @param object An `IntegrationViewsData` with its `inputs` hydrated.
#' @return TRUE invisibly, or an `integration_input_stale` refusal.
#' @noRd
integrationAssertInputsCurrent <- function(object) {
    for (entry in object@inputs) {
        if (!file.exists(entry$path)) {
            integrationRefuse(sprintf("the integration input of %s is no longer at %s; save it again at that omic's Summary and import again",
                entry$omic, entry$path), "integration_input_stale")
        }
        if (!identical(digest::digest(file = entry$path, algo = "sha256"), entry$manifest_sha256)) {
            integrationRefuse(sprintf("the integration input of %s has changed since it was imported (its manifest's sha256 no longer matches); import the inputs again",
                entry$omic), "integration_input_stale")
        }
    }
    invisible(TRUE)
}

#' Which features of a view are internal standards (INT-3): the writer's `is_internal_standard` column (lipidomics,
#' metabolomics; read back as text); none when the column is absent.
#' @param features The view's feature table.
#' @noRd
integrationInternalStandards <- function(features) {
    flag <- features[["is_internal_standard", exact = TRUE]]
    if (is.null(flag)) return(rep(FALSE, nrow(features)))
    toupper(trimws(as.character(flag))) %in% "TRUE"
}

#' Average a view's runs per biological sample (INT-2): the mean of the observed values of the sample's runs, missing
#' when none is observed.
#' @param values Features x runs.
#' @param key The biological sample of each run (in column order).
#' @return Features x biological samples, in C-locale order of the keys.
#' @noRd
integrationAverageRuns <- function(values, key) {
    levels <- sort(unique(key), method = "radix")
    out <- vapply(levels, function(level) {
        means <- rowMeans(values[, key == level, drop = FALSE], na.rm = TRUE)
        means[is.nan(means)] <- NA_real_
        means
    }, numeric(nrow(values)))
    out <- matrix(out, nrow = nrow(values), dimnames = list(rownames(values), levels))
    out
}

#' Align the selected views on their biological samples (KIT-062): per view, internal standards dropped or kept
#' (`internal_standards`), runs averaged per biological sample (`technical_replicates`, INT-2), the samples of every
#' view kept with missing values where a view lacks them, or only the samples every view has (`missing_samples`);
#' features observed in fewer than two samples dropped (and, under per-feature scaling, features without variance);
#' then the `features_per_view` most variable features (`feature_selection = top_variance`), then the per-feature
#' Z-score (`scaling = feature_zscore`). Refused: no biological sample shared by every view, a biological sample with
#' two groups, a view left without features.
#' @param object An `IntegrationViewsData`, hydrated.
#' @param parameters Resolved parameters (`integrationResolveMofaParameters`).
#' @param prefix_rows KIT-148: name each view's rows `<view>:<feature_id>` when the view is read (the model's hand-off
#'   names, unique across views), so that the hand-off makes no renamed copy; FALSE names them `<feature_id>`.
#' @return `list(views, samples, view_summary, shared_samples, notes)`: the views as features x samples matrices in
#'   the selected order with one column set, the aligned sample table, a summary row per view.
#' @noRd
integrationAlignViews <- function(object, parameters, prefix_rows = FALSE) {
    # KIT-063 item A: the user's sample match, when one is given (NULL: the design as saved).
    design <- integrationApplySampleMatch(object@design_matrix, parameters$sample_match_table)
    view_omic <- integrationViewOmics(object)
    views <- parameters$views
    averaged <- list()
    runs <- list()
    standards <- integer()
    for (view in views) {
        omic <- view_omic[[view]]
        table <- object@views[[view]]
        values <- as.matrix(table[-1L])
        storage.mode(values) <- "double"
        rownames(values) <- if (isTRUE(prefix_rows)) paste(view, table[[1L]], sep = ":") else as.character(table[[1L]])
        dropped <- 0L
        if (identical(parameters$internal_standards, "drop")) {
            standard <- integrationInternalStandards(object@features[[view]])
            values <- values[!standard, , drop = FALSE]
            dropped <- sum(standard)
        }
        standards[[view]] <- dropped
        rows <- design[design$omic == omic, , drop = FALSE]
        key <- rows$biological_sample[match(colnames(values), rows$run)]
        if (anyNA(key)) {
            integrationRefuse(sprintf("view %s has a run that its omic's design does not list", view), "integration_input_invalid")
        }
        averaged[[view]] <- integrationAverageRuns(values, key)
        runs[[view]] <- vapply(colnames(averaged[[view]]), function(level) {
            paste(sort(colnames(values)[key == level], method = "radix"), collapse = ";")
        }, character(1))
    }
    present <- lapply(averaged, colnames)
    plan <- integrationAlignSamples(design, view_omic, views, present, parameters$missing_samples)
    zscore <- identical(parameters$scaling, "feature_zscore")
    summary <- list()
    final <- list()
    for (view in views) {
        block <- integrationAlignedRows(averaged[[view]], plan$kept, zscore)
        m <- block$values
        if (identical(parameters$feature_selection, "top_variance") && nrow(m) > parameters$features_per_view) {
            m <- m[integrationSelectTopFeatures(integrationFeatureSelectionStatistic(m), parameters$features_per_view), ,
                drop = FALSE]
        }
        if (zscore) m <- integrationZscoreRows(m)
        if (!nrow(m)) integrationRefuseEmptyView(view)
        final[[view]] <- m
        summary[[view]] <- integrationViewSummaryRow(view, view_omic[[view]], nrow(averaged[[view]]) + standards[[view]],
            standards[[view]], block$sparse, block$constant, nrow(m), block$present)
    }
    integrationAlignedValue(final, plan, present, runs, views, summary)
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-149: the alignment's rules, shared by the whole-matrix alignment (`integrationAlignViews`) and the partitioned
# one (`integrationAlignViewsPartitioned`), so that the two apply one definition. Every rule on values is row by row
# (each feature's average, its count of observed samples, its spread, its selection statistic and Z-score depend on
# that feature's values alone); only the most-variable selection compares features, and it is taken once over the
# whole data set from one statistic per feature.
# ---------------------------------------------------------------------------------------------------------------

#' The omic of each data set of the inputs state, named by data set (the inputs' order).
#' @param object An `IntegrationViewsData` with its `inputs` hydrated.
#' @noRd
integrationViewOmics <- function(object) {
    unlist(unname(lapply(object@inputs, function(entry) stats::setNames(rep(entry$omic, length(entry$views)), entry$views))))
}

#' The biological samples the alignment keeps (KIT-062): those every data set measures (`drop_incomplete`) or those any
#' data set measures (`keep_missing`), each with its group. Refused: no biological sample shared by every data set, a
#' biological sample with two groups. Fewer shared samples than `INTEGRATION_FEW_SAMPLES` is warned
#' (`integration_few_samples`) and kept as a note.
#' @param design The design with the sample match applied.
#' @param view_omic `integrationViewOmics()`.
#' @param views The selected data sets.
#' @param present Named list by data set: the biological samples it measures, in C-locale order.
#' @param missing_samples The parameter's value.
#' @return `list(shared, kept, groups, notes)`.
#' @noRd
integrationAlignSamples <- function(design, view_omic, views, present, missing_samples) {
    shared <- sort(Reduce(intersect, present), method = "radix")
    if (!length(shared)) {
        integrationRefuse(sprintf("the views %s share no biological sample; check the biological-sample column of each omic",
            paste(views, collapse = ", ")), "integration_no_shared_samples")
    }
    kept <- if (identical(missing_samples, "drop_incomplete")) shared else sort(unique(unlist(present)), method = "radix")
    selected <- design[design$omic %in% unique(view_omic[views]) & design$biological_sample %in% kept, , drop = FALSE]
    groups <- tapply(selected$group, selected$biological_sample, function(group) unique(group[!is.na(group)]),
        simplify = FALSE)
    conflicting <- names(groups)[lengths(groups) > 1L]
    if (length(conflicting)) {
        integrationRefuse(sprintf("biological sample %s has more than one group across the omics", conflicting[[1L]]),
            "integration_group_conflict")
    }
    notes <- character()
    if (length(shared) < INTEGRATION_FEW_SAMPLES) {
        notes <- sprintf("only %d biological samples are measured in every view (fewer than %d); factor models need more samples",
            length(shared), INTEGRATION_FEW_SAMPLES)
        cli::cli_warn(notes, class = "integration_few_samples")
    }
    list(shared = shared, kept = kept, groups = groups, notes = notes)
}

#' Rows of a data set on the kept biological samples (KIT-062's rules): placed on the kept samples (missing where the
#' data set lacks one), the features observed in fewer than two samples dropped and, under per-feature scaling, those
#' without variance. Row by row: a block of rows gives the rows the whole data set gives.
#' @param averaged Features x biological samples (`integrationAverageRuns`).
#' @param kept The kept biological samples.
#' @param zscore Whether features without variance are dropped (per-feature scaling).
#' @return `list(values, sparse, constant, present)`: the rows left, the counts dropped as sparse and as constant, and
#'   how many kept samples the data set measures.
#' @noRd
integrationAlignedRows <- function(averaged, kept, zscore) {
    m <- matrix(NA_real_, nrow(averaged), length(kept), dimnames = list(rownames(averaged), kept))
    here <- intersect(kept, colnames(averaged))
    m[, here] <- averaged[, here, drop = FALSE]
    sparse <- rowSums(!is.na(m)) < 2L
    m <- m[!sparse, , drop = FALSE]
    constant <- 0L
    if (isTRUE(zscore)) {
        spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
        flat <- !(spread > 0)
        constant <- sum(flat)
        m <- m[!flat, , drop = FALSE]
    }
    list(values = m, sparse = sum(sparse), constant = constant, present = length(here))
}

#' The statistic the most-variable selection ranks features by (`feature_selection = top_variance`): each feature's
#' variance over the samples it is observed in. The one definition both alignments use; the selection rule is under
#' scientific review (KIT-151), and a change made here changes both alike.
#' @param m Features x samples.
#' @return One number per feature, in row order.
#' @noRd
integrationFeatureSelectionStatistic <- function(m) apply(m, 1L, stats::var, na.rm = TRUE)

#' The features the most-variable selection keeps: the `n` largest statistics, a tie broken by position (the earlier
#' feature first), returned in position order.
#' @param statistic `integrationFeatureSelectionStatistic()` of every feature of a data set, in its row order.
#' @param n The number kept (`features_per_view`).
#' @return Positions in `statistic`, increasing.
#' @noRd
integrationSelectTopFeatures <- function(statistic, n) sort(order(-statistic, seq_along(statistic))[seq_len(n)])

#' The per-feature Z-score (`scaling = feature_zscore`): each feature centred on its mean and divided by its standard
#' deviation over the samples it is observed in.
#' @param m Features x samples.
#' @noRd
integrationZscoreRows <- function(m) {
    centre <- rowMeans(m, na.rm = TRUE)
    spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
    (m - centre) / spread
}

#' Refuse a data set left without features.
#' @noRd
integrationRefuseEmptyView <- function(view) {
    integrationRefuse(sprintf("view %s has no feature left for the model (none is observed in two or more samples)", view),
        "integration_view_empty")
}

#' One row of the alignment's data set summary.
#' @noRd
integrationViewSummaryRow <- function(view, omic, features_input, standards, sparse, constant, features_model, present) {
    data.frame(view = view, omic = omic, features_input = features_input, internal_standards_dropped = standards,
        sparse_dropped = sparse, constant_dropped = constant, features_model = features_model, samples_present = present,
        stringsAsFactors = FALSE)
}

#' The alignment's value from its parts (`integrationAlignViews`' form).
#' @param final The aligned data sets (NULL when not built).
#' @param plan `integrationAlignSamples()`.
#' @param present,runs Named lists by data set: the biological samples it measures, and per sample its runs.
#' @param views The selected data sets.
#' @param summary List of `integrationViewSummaryRow()` by data set.
#' @noRd
integrationAlignedValue <- function(final, plan, present, runs, views, summary) {
    kept <- plan$kept
    groups <- plan$groups
    samples <- data.frame(biological_sample = kept,
        group = unname(vapply(kept, function(name) if (length(groups[[name]])) groups[[name]][[1L]] else NA_character_,
            character(1))),
        views_present = unname(vapply(kept, function(name) sum(vapply(present, function(p) name %in% p, logical(1))),
            integer(1))),
        stringsAsFactors = FALSE, row.names = NULL)
    for (view in views) samples[[paste0(view, "_runs")]] <- unname(ifelse(kept %in% names(runs[[view]]),
        runs[[view]][kept], ""))
    summary <- do.call(rbind, summary)
    rownames(summary) <- NULL
    list(views = final, samples = samples, view_summary = summary, shared_samples = length(plan$shared), notes = plan$notes)
}

#' Align the selected data sets as `integrationAlignViews` does, reading each data set in row pages (KIT-149): no data
#' set is held whole. The biological samples, runs and internal standards come from the design, the feature tables
#' and the data sets' column names, without values. Then per data set, pass 1 reads its pages and keeps, per feature
#' left after the sparse and constant filters, its selection statistic (one number per feature); the most-variable
#' selection is taken once over the whole data set from those numbers (`integrationSelectTopFeatures`, ties by
#' position); pass 2 reads the pages again and keeps the selected rows, Z-scored under per-feature scaling, in a matrix
#' of the selected features only. Every rule is the whole-matrix path's (`integrationAlignedRows` and the functions
#' beside it), applied row by row, so the value is `identical()` to `integrationAlignViews()`'s.
#' @param object An `IntegrationViewsData` with `features`, `design_matrix` and `inputs` hydrated; its `views` may be
#'   zero-row prototypes (their column names are read).
#' @param parameters Resolved parameters (`integrationResolveMofaParameters`).
#' @param pages `function(view, columns, kernel)` calling `kernel(page, offset)` on each page of the data set's
#'   `columns` in row order (`offset`: the page's first row's zero-based position): `integrationStateViewPages()` or
#'   `integrationObjectViewPages()`.
#' @param prefix_rows As `integrationAlignViews()`.
#' @param values FALSE: pass 2 is skipped and `views` is NULL (the counts and the sample table only).
#' @return `integrationAlignViews()`'s value.
#' @noRd
integrationAlignViewsPartitioned <- function(object, parameters, pages, prefix_rows = FALSE, values = TRUE) {
    design <- integrationApplySampleMatch(object@design_matrix, parameters$sample_match_table)
    view_omic <- integrationViewOmics(object)
    views <- parameters$views
    drop_standards <- identical(parameters$internal_standards, "drop")
    columns <- list()
    keys <- list()
    present <- list()
    runs <- list()
    standards <- list()
    for (view in views) {
        columns[[view]] <- names(object@views[[view]])
        run_names <- columns[[view]][-1L]
        standards[view] <- list(if (drop_standards) integrationInternalStandards(object@features[[view]]))
        rows <- design[design$omic == view_omic[[view]], , drop = FALSE]
        key <- rows$biological_sample[match(run_names, rows$run)]
        if (anyNA(key)) {
            integrationRefuse(sprintf("view %s has a run that its omic's design does not list", view), "integration_input_invalid")
        }
        keys[[view]] <- key
        present[[view]] <- sort(unique(key), method = "radix")
        runs[[view]] <- vapply(present[[view]], function(level) {
            paste(sort(run_names[key == level], method = "radix"), collapse = ";")
        }, character(1))
    }
    plan <- integrationAlignSamples(design, view_omic, views, present, parameters$missing_samples)
    zscore <- identical(parameters$scaling, "feature_zscore")
    cap <- if (identical(parameters$feature_selection, "top_variance")) parameters$features_per_view else NA_integer_
    summary <- list()
    final <- if (isTRUE(values)) list()
    for (view in views) {
        # One page of the data set on the kept samples: the whole path's read, standards and averaging, row by row.
        block_of <- function(page, offset) {
            m <- as.matrix(page[-1L])
            storage.mode(m) <- "double"
            rownames(m) <- if (isTRUE(prefix_rows)) paste(view, page[[1L]], sep = ":") else as.character(page[[1L]])
            if (drop_standards) m <- m[!standards[[view]][offset + seq_len(nrow(m))], , drop = FALSE]
            integrationAlignedRows(integrationAverageRuns(m, keys[[view]]), plan$kept, zscore)
        }
        # Pass 1: the counts, and the selection statistic of every feature left.
        rows_read <- 0L
        sparse <- 0L
        constant <- 0L
        left <- 0L
        statistic <- list()
        pages(view, columns[[view]], function(page, offset) {
            if (!nrow(page)) return(invisible(NULL))
            block <- block_of(page, offset)
            rows_read <<- rows_read + nrow(page)
            sparse <<- sparse + block$sparse
            constant <<- constant + block$constant
            left <<- left + nrow(block$values)
            if (!is.na(cap)) statistic[[length(statistic) + 1L]] <<- unname(integrationFeatureSelectionStatistic(block$values))
            invisible(NULL)
        })
        # The global statistic: the cap's selection over the whole data set (NULL: every feature left is kept).
        take <- NULL
        if (!is.na(cap) && left > cap) {
            take <- logical(left)
            take[integrationSelectTopFeatures(unlist(statistic), cap)] <- TRUE
        }
        statistic <- NULL
        features_model <- if (is.null(take)) left else sum(take)
        if (!features_model) integrationRefuseEmptyView(view)
        # Pass 2: the selected rows, Z-scored under per-feature scaling, into a matrix of the selected features only.
        if (isTRUE(values)) {
            out <- matrix(NA_real_, features_model, length(plan$kept))
            ids <- character(features_model)
            seen <- 0L
            filled <- 0L
            pages(view, columns[[view]], function(page, offset) {
                if (!nrow(page)) return(invisible(NULL))
                m <- block_of(page, offset)$values
                rows <- nrow(m)
                if (!is.null(take)) m <- m[take[seen + seq_len(rows)], , drop = FALSE]
                seen <<- seen + rows
                if (nrow(m)) {
                    if (zscore) m <- integrationZscoreRows(m)
                    at <- filled + seq_len(nrow(m))
                    out[at, ] <<- m
                    ids[at] <<- rownames(m)
                    filled <<- filled + nrow(m)
                }
                invisible(NULL)
            })
            dimnames(out) <- list(ids, plan$kept)
            final[[view]] <- out
            out <- NULL
        }
        dropped <- if (drop_standards) sum(standards[[view]]) else 0L
        summary[[view]] <- integrationViewSummaryRow(view, view_omic[[view]], rows_read, dropped, sparse, constant,
            features_model, length(intersect(plan$kept, present[[view]])))
    }
    integrationAlignedValue(final, plan, present, runs, views, summary)
}

#' Pages of the inputs state's data sets from the integration store (KIT-149; `integrationAlignViewsPartitioned`'s
#' `pages`): each data set's table read by `partitionByRows`, at most 10,000 rows per page, a table of more than 64
#' columns read in blocks of 64 columns bound in column order.
#' @param manager The integration owner.
#' @param source The inputs state's reference.
#' @param worker_cap The admitted worker memory (NULL: 12 GiB, the omics' partitioned steps' default).
#' @noRd
integrationStateViewPages <- function(manager, source, worker_cap = NULL) {
    function(view, columns, kernel) withWorkflowStateScope(manager, function() {
        context <- manager$.context
        table <- readWorkflowStateMetadata(manager, source)$data$components$views[[view]]
        scope <- newWorkflowTableReadScope(context$store, table, reader = if (length(columns) > 64L) "duckdb" else "arrow")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        partitionByRows(scope, columns, kernel = function(page, index, offset) kernel(page, offset),
            consume = function(value, index) invisible(NULL), worker_cap = worker_cap %||% (12 * 1024^3),
            working_bytes_per_row = omicsPartitionedWorkingBytesPerColumn * length(columns),
            overhead_bytes = 2 * 67108864 + 2 * 8388608, policy = "integration-row-pages-v1")
        invisible(NULL)
    })
}

#' Pages of a hydrated inputs state's data sets (KIT-149): each data set as one page.
#' @param object An `IntegrationViewsData` with its `views` hydrated.
#' @noRd
integrationObjectViewPages <- function(object) {
    function(view, columns, kernel) {
        kernel(object@views[[view]][columns], 0L)
        invisible(NULL)
    }
}

#' The MOFA2 options the declared parameters set (every other MOFA2 option stays at its default).
#' @noRd
integrationMofaOptions <- function(parameters) {
    parameters[c("scale_views", "center_groups", "num_factors", "convergence_mode", "maxiter", "seed", "precision")]
}

#' Whether MOFA2 can train here (KIT-062; D6: Python optional and probed): the R packages MOFA2 and reticulate, and a
#' Python reticulate reaches with the module `mofapy2` (the bundle sets RETICULATE_PYTHON). The package loads without
#' either (both are in Suggests).
#' @return `list(available, reason)`.
#' @noRd
integrationMofaCapability <- function() {
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        return(list(available = FALSE, reason = "the R package MOFA2 is not installed"))
    }
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        return(list(available = FALSE, reason = "the R package reticulate is not installed"))
    }
    found <- tryCatch(reticulate::py_module_available("mofapy2"), error = function(error) FALSE)
    if (!isTRUE(found)) {
        python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
        return(list(available = FALSE, reason = sprintf("no Python with the module mofapy2 is available to reticulate%s",
            if (nzchar(python)) sprintf(" (RETICULATE_PYTHON is %s)", python) else " (RETICULATE_PYTHON is not set)")))
    }
    list(available = TRUE, reason = NA_character_)
}

#' The versions a model was trained with: R, MOFA2, reticulate, Python and its packages.
#' @noRd
integrationMofaVersions <- function() {
    metadata <- tryCatch(reticulate::import("importlib.metadata"), error = function(error) NULL)
    python_package <- function(name) tryCatch(as.character(metadata$version(name)), error = function(error) NA_character_)
    list(R = as.character(getRversion()), MOFA2 = as.character(utils::packageVersion("MOFA2")),
        reticulate = as.character(utils::packageVersion("reticulate")),
        python = tryCatch(as.character(reticulate::py_config()$version_string), error = function(error) NA_character_),
        mofapy2 = python_package("mofapy2"), numpy = python_package("numpy"), scipy = python_package("scipy"),
        h5py = python_package("h5py"))
}

#' MOFA2's check that a trained factor follows each sample's overall level (MOFA2 1.22.0, `.quality_control` of a
#' loaded model), on the data the model was trained with (KIT-148). MOFA2 runs the check only when the model file holds
#' the training data; the model file no longer does (`save_data = FALSE`), so the same rule is applied here: per data
#' set, the features centred on their mean when MOFA2 centred them (`center_groups`), the mean of each sample over the
#' observed features, and its absolute correlation with each factor over the samples (pairwise complete); a factor
#' above 0.75 in any data set is named. Scaling a data set by one number (`scale_views`) does not change a correlation.
#' One data set is held at a time.
#' @param views Named list of features x samples matrices as handed to MOFA2.
#' @param factors Samples x factors (the trained model's, named by sample).
#' @param center Whether MOFA2 centred the features (`center_groups`).
#' @return MOFA2's warning text, or NULL when no factor follows the overall level.
#' @noRd
integrationMofaOverallLevelWarning <- function(views, factors, center) {
    r <- vapply(views, function(m) {
        if (isTRUE(center)) m <- m - rowMeans(m, na.rm = TRUE)
        means <- colMeans(m, na.rm = TRUE)
        abs(suppressWarnings(stats::cor(means, factors[colnames(m), , drop = FALSE], use = "pairwise.complete.obs")))[1L, ]
    }, numeric(ncol(factors)))
    r <- matrix(r, nrow = ncol(factors))
    flagged <- which(rowSums(r > 0.75) > 0)
    if (!length(flagged)) return(NULL)
    # MOFA2 1.22.0's wording (QC.R), so that the Results tab reads it as it read MOFA2's own warning.
    sprintf(paste("Factor(s) %s are strongly correlated with the average expression of features for at least one of",
        "your omics. Such factors appear when there are differences in the total 'levels' between your samples,",
        "*sometimes* because of poor normalisation in the preprocessing steps."), paste(flagged, collapse = ", "))
}

#' Whether the direct hand-off is used (KIT-148 item 2): asked for (`handoff = "direct"`, the internal option
#' `MultiScholaR.integration_mofa_handoff`; the default is `run_mofa`), and the versions it is written for are the ones
#' installed (MOFA2 1.22.0, mofapy2 0.7.3); otherwise MOFA2's own `run_mofa` hands the data over.
#' @param handoff `run_mofa` or `direct`.
#' @noRd
integrationMofaDirectHandoff <- function(handoff) {
    if (!identical(handoff, "direct") || !requireNamespace("MOFA2", quietly = TRUE) ||
        !identical(as.character(utils::packageVersion("MOFA2")), "1.22.0")) return(FALSE)
    mofapy2 <- tryCatch(as.character(reticulate::import("mofapy2")$version$`__version__`), error = function(error) NA_character_)
    identical(mofapy2, "0.7.3")
}

#' Hand a prepared MOFA2 object to mofapy2 one data set at a time and train it (KIT-148 item 2; internal, off by
#' default). MOFA2 1.22.0's `run_mofa` (its `.run_mofa_reticulate`) converts every data set to Python before it hands
#' the list over, while its R copy of every data set is still held; here each data set is converted and its R copy
#' released before the next, through the same mofapy2 0.7.3 entry point with the same arguments in the same order, so
#' the model is the same. Written for those two versions only (`integrationMofaDirectHandoff`); no covariates (MEFISTO).
#' @param holder An environment holding the prepared object as `object`; the object is released from it.
#' @param outfile The model file; written without the training data.
#' @return `outfile`, invisibly.
#' @noRd
integrationRunMofapy2 <- function(holder, outfile) {
    if (!is.null(holder$object@covariates)) cli::cli_abort("the direct hand-off takes no covariates", class = "integration_refused")
    mofa <- reticulate::import("mofapy2")
    entry <- mofa$run$entry_point$entry_point()
    data_options <- holder$object@data_options
    entry$set_data_options(scale_views = data_options$scale_views, scale_groups = data_options$scale_groups,
        center_groups = data_options$center_groups, use_float32 = data_options$use_float32)
    samples <- holder$object@samples_metadata
    entry$data_opts$samples_metadata <- reticulate::r_to_py(lapply(data_options$groups, function(g) samples[samples$group == g, ]))
    features <- holder$object@features_metadata
    entry$data_opts$features_metadata <- reticulate::r_to_py(unname(lapply(data_options$views,
        function(m) features[features$view == m, ])))
    one <- function(xs) if (length(xs) > 1) xs else list(xs)
    samples_names <- reticulate::r_to_py(lapply(unname(lapply(holder$object@data[[1L]], colnames)), one))
    features_names <- reticulate::r_to_py(lapply(unname(lapply(holder$object@data, function(x) rownames(x[[1L]]))), one))
    data <- reticulate::r_to_py(list())
    for (m in seq_along(holder$object@data)) {
        groups <- reticulate::r_to_py(list())
        for (g in seq_along(holder$object@data[[m]])) groups$append(reticulate::r_to_py(t(holder$object@data[[m]][[g]])))
        data$append(groups)
        # This data set is now mofapy2's; its R copy is released before the next one is converted.
        holder$object@data[m] <- list(NULL)
    }
    entry$set_data_matrix(data = data, likelihoods = unname(holder$object@model_options$likelihoods),
        views_names = reticulate::r_to_py(as.list(data_options$views)), groups_names = reticulate::r_to_py(as.list(data_options$groups)),
        samples_names = samples_names, features_names = features_names)
    model_options <- holder$object@model_options
    entry$set_model_options(factors = model_options$num_factors, spikeslab_factors = model_options$spikeslab_factors,
        spikeslab_weights = model_options$spikeslab_weights, ard_factors = model_options$ard_factors,
        ard_weights = model_options$ard_weights)
    training <- holder$object@training_options
    entry$set_train_options(iter = training$maxiter, convergence_mode = training$convergence_mode,
        dropR2 = training$drop_factor_threshold, startELBO = training$startELBO, freqELBO = training$freqELBO,
        seed = training$seed, gpu_mode = training$gpu_mode, gpu_device = training$gpu_device, verbose = training$verbose,
        outfile = training$outfile, weight_views = training$weight_views, save_interrupted = training$save_interrupted)
    if (isTRUE(training$stochastic)) {
        stochastic <- holder$object@stochastic_options
        entry$set_stochastic_options(learning_rate = stochastic$learning_rate, forgetting_rate = stochastic$forgetting_rate,
            batch_size = stochastic$batch_size, start_stochastic = stochastic$start_stochastic)
    }
    holder$object <- NULL
    data <- NULL
    reclaimWorkflowMemory()
    entry$build()
    entry$run()
    entry$save(outfile, save_data = FALSE)
    invisible(outfile)
}

#' Train MOFA2 as the proven lane does (KIT-060 status 2026-11-19, run 44423980807): `create_mofa` on the views,
#' the default data, model and training options with the declared ones set, `prepare_mofa`, `run_mofa(use_basilisk =
#' FALSE)` (a fall-back to basilisk is refused: the bundle's Python must be the one used). KIT-148: the model file is
#' written without the training data (`save_data = FALSE`; the aligned values are kept in the project's store instead),
#' the model `run_mofa` loads is the one used (one load, not two), the R copy MOFA2 made for the hand-off is released
#' after training, and MOFA2's overall-level check, which needs the training data, is applied to the views here
#' (`integrationMofaOverallLevelWarning`), its warning kept where MOFA2 would have given it.
#' @param views Named list of features x samples matrices (feature names unique across views).
#' @param options `integrationMofaOptions()`.
#' @param handoff KIT-148 item 2 (internal): `run_mofa` (the default) hands the data over through MOFA2's `run_mofa`;
#'   `direct` one data set at a time through mofapy2's entry point (`integrationRunMofapy2`; MOFA2 1.22.0 and mofapy2
#'   0.7.3 only, `run_mofa` otherwise), the same model; then the model is loaded from its file as `run_mofa` loads it.
#' @return `list(factors, weights, r2_per_factor (factors x views), r2_total, elbo, iterations, options, versions,
#'   messages, model_file)`.
#' @noRd
integrationTrainMofa <- function(views, options, handoff = getOption("MultiScholaR.integration_mofa_handoff", "run_mofa")) {
    object <- MOFA2::create_mofa(views)
    data_options <- MOFA2::get_default_data_options(object)
    data_options$scale_views <- options$scale_views
    data_options$center_groups <- options$center_groups
    # INT-6: MOFA2 trains in single precision unless told otherwise; double is the shipped setting.
    data_options$use_float32 <- identical(options$precision, "single")
    model_options <- MOFA2::get_default_model_options(object)
    model_options$num_factors <- options$num_factors
    training_options <- MOFA2::get_default_training_options(object)
    training_options$seed <- options$seed
    training_options$maxiter <- options$maxiter
    training_options$convergence_mode <- options$convergence_mode
    training_options$verbose <- FALSE
    object <- MOFA2::prepare_mofa(object, data_options = data_options, model_options = model_options,
        training_options = training_options)
    reclaimWorkflowMemory()
    outfile <- tempfile("integration-mofa-", fileext = ".hdf5")
    on.exit(unlink(outfile), add = TRUE)
    # Kept, not dropped: MOFA2's messages and warnings are returned with the model. KIT-063: they include the checks
    # MOFA2 makes when `run_mofa` loads the trained model, for the Results tab.
    messages <- character()
    direct <- integrationMofaDirectHandoff(handoff)
    model <- withCallingHandlers(
        if (direct) {
            holder <- new.env(parent = emptyenv())
            holder$object <- object
            object <- NULL
            integrationRunMofapy2(holder, outfile)
            MOFA2::load_model(outfile)
        } else {
            MOFA2::run_mofa(object, outfile = outfile, save_data = FALSE, use_basilisk = FALSE)
        },
        message = function(m) {
            messages <<- c(messages, trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            messages <<- c(messages, paste("WARNING:", trimws(conditionMessage(w))))
            invokeRestart("muffleWarning")
        })
    object <- NULL
    reclaimWorkflowMemory()
    if (any(grepl("Setting use_basilisk = TRUE", messages, fixed = TRUE))) {
        integrationRefuse("MOFA2 could not use the configured Python with mofapy2 and fell back to basilisk; the model is not trained",
            "integration_mofa_unavailable")
    }
    # KIT-148: MOFA2's overall-level check, placed where its loading gave it (before the correlated-factors check, the
    # last of its checks).
    overall <- integrationMofaOverallLevelWarning(views, MOFA2::get_factors(model)[[1L]], options$center_groups)
    if (!is.null(overall) && !paste("WARNING:", overall) %in% messages) {
        at <- grep("WARNING: The model contains highly correlated factors", messages, fixed = TRUE)
        messages <- if (length(at)) append(messages, paste("WARNING:", overall), after = at[[1L]] - 1L) else
            c(messages, paste("WARNING:", overall))
    }
    variance <- MOFA2::get_variance_explained(model)
    r2 <- variance$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
    stats_slot <- methods::slot(model, "training_stats")
    elbo <- as.numeric(stats_slot$elbo)
    list(factors = MOFA2::get_factors(model)[[1L]], weights = MOFA2::get_weights(model), r2_per_factor = r2,
        r2_total = variance$r2_total[[1L]], elbo = elbo, iterations = length(stats_slot$time),
        options = list(data = unclass(data_options), model = unclass(model_options), training = unclass(training_options)),
        versions = integrationMofaVersions(), messages = messages,
        model_file = list(bytes = readBin(outfile, "raw", file.size(outfile)),
            sha256 = digest::digest(file = outfile, algo = "sha256"), size = file.size(outfile)))
}

# KIT-148: the model product v2. The small tables are kept in the product itself; each data set's weights and the values
# handed to MOFA2 are tables of the integration store (one per data set, read by page); the model file is an entity of
# its own. A reopened store refuses data frames among the head's roots, so the head keeps only the product's reference.
INTEGRATION_MODEL_SCHEMA <- "multischolar.integration_mofa_model.v2"

#' The model's saved product (KIT-062; KIT-148 v2): `schema`; inline the factors (`biological_sample`, then one column
#' per factor), the variance explained per factor and data set and in all per data set, the ELBO trace summary, the
#' aligned sample table and data set summary, the shared samples, notes, MOFA2's options, versions and messages; per data
#' set (named by data set, in the model's order) a table reference in `weights` (`feature_id`, then one column per
#' factor) and in `aligned_values` (`feature_id`, then one double column per biological sample in the factors' order:
#' the values handed to MOFA2, before MOFA2 centres and scales them), both with the data set's model features in one
#' order; and `model_file` (`reference` to the model file's own entity, `sha256`, `size`, `media_type`). Each data set's
#' tables are written one data set at a time.
#' @param aligned `integrationAlignViews()`'s value, its rows named `<view>:<feature_id>`.
#' @param fit The trainer's value (`integrationTrainMofa()`'s form), on those views.
#' @param store The operation's store handle (`store$table`, `store$product`).
#' @noRd
integrationModelProduct <- function(aligned, fit, store) {
    factors <- fit$factors
    factor_table <- data.frame(biological_sample = rownames(factors), factors, check.names = FALSE,
        stringsAsFactors = FALSE, row.names = NULL)
    table <- function(m, view) data.frame(feature_id = substring(rownames(m), nchar(view) + 2L), m, check.names = FALSE,
        stringsAsFactors = FALSE, row.names = NULL)
    weights <- list()
    values <- list()
    for (view in names(aligned$views)) {
        weights[[view]] <- store$table(table(fit$weights[[view]], view))
        values[[view]] <- store$table(table(aligned$views[[view]], view))
    }
    r2 <- fit$r2_per_factor
    variance <- data.frame(factor = rep(rownames(r2), times = ncol(r2)), view = rep(colnames(r2), each = nrow(r2)),
        r2_percent = as.vector(r2), stringsAsFactors = FALSE)
    elbo <- fit$elbo[is.finite(fit$elbo)]
    file <- fit$model_file
    file_reference <- store$product(list(bytes = file$bytes %||% raw(0), sha256 = file$sha256 %||% NA_character_,
        size = file$size %||% NA_real_), "model_file", "integration.mofa_model_file")
    list(schema = INTEGRATION_MODEL_SCHEMA, factors = factor_table, weights = weights, aligned_values = values,
        variance_explained = variance,
        variance_total = data.frame(view = names(fit$r2_total), r2_percent = unname(as.numeric(fit$r2_total)),
            stringsAsFactors = FALSE),
        elbo = list(iterations = as.integer(fit$iterations), evaluations = length(elbo),
            first = if (length(elbo)) elbo[[1L]] else NA_real_, last = if (length(elbo)) elbo[[length(elbo)]] else NA_real_,
            trace = elbo),
        samples = aligned$samples, view_summary = aligned$view_summary, shared_samples = aligned$shared_samples,
        notes = aligned$notes, options = fit$options, versions = fit$versions, messages = fit$messages,
        model_file = list(reference = file_reference, sha256 = file$sha256 %||% NA_character_, size = file$size %||% NA_real_,
            media_type = "application/x-hdf5"))
}

#' Check that a value is a model product v2 (KIT-148), refused in words otherwise (`integration_model_outdated`).
#' @noRd
integrationAssertModelProduct <- function(value) {
    if (!is.list(value) || !identical(value$schema, INTEGRATION_MODEL_SCHEMA)) {
        integrationRefuse("the integration model was saved in an earlier form by a development build; train it again",
            "integration_model_outdated")
    }
    invisible(value)
}

#' The table reference of one data set's weights or aligned values in a model product (KIT-148).
#' @param value The model product (v2).
#' @param kind `weights` or `aligned_values`.
#' @param view The data set.
#' @noRd
integrationModelTableReference <- function(value, kind = c("weights", "aligned_values"), view) {
    kind <- match.arg(kind)
    integrationAssertModelProduct(value)
    reference <- value[[kind]][[view, exact = TRUE]]
    if (!isWorkflowTableReference(reference)) {
        integrationRefuse(sprintf("the integration model has no %s table for the data set %s", kind, view),
            "integration_model_outdated")
    }
    reference
}

#' One data set's weights or aligned values of a model product, whole (KIT-148): one data set is read.
#' @param manager The owner of the store that holds the model.
#' @param value The model product (v2).
#' @param kind `weights` or `aligned_values`.
#' @param view The data set.
#' @return `data.frame(feature_id, ...)`: one column per factor (weights) or per biological sample (aligned values).
#' @noRd
integrationModelViewTable <- function(manager, value, kind, view) {
    reference <- integrationModelTableReference(value, kind, view)
    withWorkflowStateScope(manager, function() materializeWorkflowTable(manager$.context, reference))
}

#' One page of rows of a data set's weights or aligned values (KIT-148): rows `offset + 1` to `offset + limit` in the
#' model's feature order, every column. The tables' layout: one row per model feature of the data set, `feature_id`
#' first, then one double column per factor (weights) or per biological sample (aligned values; a study may have
#' hundreds). The store's read rule is at most 10,000 rows and 64 columns per read (`readSavedWorkflowTable`), so a page
#' is at most 10,000 rows and is read in blocks of at most 64 columns, bound in column order.
#' @param manager,value,kind,view As `integrationModelViewTable()`.
#' @param offset Rows skipped (zero-based position of the page's first row).
#' @param limit Rows read (at most 10,000).
#' @return The page, a data frame (automatic row names).
#' @noRd
integrationModelTablePage <- function(manager, value, kind, view, offset = 0L, limit = 10000L) {
    reference <- integrationModelTableReference(value, kind, view)
    withWorkflowStateScope(manager, function() {
        store <- manager$.context$store
        scope <- newWorkflowTableReadScope(store, reference)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        columns <- scope$metadata$columns
        blocks <- split(columns, ceiling(seq_along(columns) / 64L))
        parts <- lapply(blocks, function(block) {
            part <- readSavedWorkflowTable(store, reference, block, offset = offset, limit = limit, scope = scope)$data
            as.list(part)
        })
        page <- as.data.frame(stats::setNames(do.call(c, unname(parts)), columns), stringsAsFactors = FALSE,
            check.names = FALSE, optional = TRUE)
        page
    })
}

#' Every data set's weights of a model product as one table (KIT-148): `view`, `feature_id`, one column per factor,
#' the data sets in the model's order, read one at a time (the form of the product v1's `weights`).
#' @noRd
integrationModelWeights <- function(manager, value) {
    integrationAssertModelProduct(value)
    tables <- lapply(names(value$weights), function(view) {
        w <- integrationModelViewTable(manager, value, "weights", view)
        data.frame(view = rep(view, nrow(w)), w, check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
    })
    out <- do.call(rbind, tables)
    rownames(out) <- NULL
    out
}

#' The model file's bytes (KIT-148): the HDF5 MOFA2 wrote, without the training data, from its own entity.
#' @noRd
integrationModelFileBytes <- function(manager, value) {
    integrationAssertModelProduct(value)
    reference <- value$model_file$reference
    if (!isSavedWorkflowReference(reference)) {
        integrationRefuse("the integration model has no saved model file", "integration_model_outdated")
    }
    withWorkflowStateScope(manager, function() hydrateSavedWorkflowEntity(manager$.context$store, reference)$value$bytes)
}

#' A model product as the comparison and the replay hold it outside its store (KIT-148): the inline tables as they
#' are, `weights` as one table (`integrationModelWeights()`), `model_file` without its reference; no table reference is
#' left in it, so it stays readable once the store is closed.
#' @noRd
integrationModelComparable <- function(manager, value) {
    value$weights <- integrationModelWeights(manager, value)
    value$aligned_values <- NULL
    value$model_file$reference <- NULL
    value
}

#' The MOFA2 model declaration (KIT-062): whole matrix, on the inputs state, product only. Its identity is the state's
#' (every view's table, the design), the manifests' sha256, the views and every method parameter. A caller may hand
#' the kernel a `trainer` beside the declared parameters (not identity; the tests' stand-in for MOFA2).
#' @param match FALSE: the registry's declaration. TRUE (KIT-063 item A): the same declaration with the user's sample
#'   match as one more parameter (`sample_match`, the match's sha256), the match itself passed beside it
#'   (`sample_match_table`). A model without a match keeps the registry's declaration, so its identity is unchanged.
#' KIT-149: the declaration's `partitioned` variant aligns the selected data sets from the store in row pages
#' (`integrationAlignViewsPartitioned`), never hydrating a data set whole, and gives the product of the whole-matrix
#' kernel (the alignment is `identical()`); the variant is not part of the identity.
#' @noRd
integrationMofaOperation <- function(match = FALSE) {
    declaration <- newWorkflowOperation(INTEGRATION_MODEL_OPERATION, "IntegrationViewsData",
        # KIT-149 (audit ENG-05): the data sets are not among the fields read before the kernel; the kernel reads the
        # selected ones only.
        fields = c("features", "design_matrix", "inputs"), unit = "whole_matrix", output_class = NULL,
        parameters = c(integrationMofaParameterNames(), if (isTRUE(match)) "sample_match"),
        product_type = "integration.mofa_model",
        # KIT-148: the product's form is part of what a saved model is looked up by, so a model saved in the earlier
        # form (v1) is not reused.
        factors = list(product_schema = INTEGRATION_MODEL_SCHEMA),
        kernel = function(object, parameters, store) {
            resolved <- integrationResolveMofaParameters(object, parameters)
            # KIT-149 (audit ENG-05): only the selected data sets are read whole; the others stay zero-row prototypes.
            object@views[resolved$views] <- readWorkflowExecutionComponents(store$manager, store$source,
                fields = paste0("views/", resolved$views))@views[resolved$views]
            # KIT-148: the rows carry the hand-off names (`<view>:<feature_id>`) from the read, so the aligned views are
            # handed to the trainer as they are, with no renamed copy.
            integrationModelFromAligned(integrationAlignViews(object, resolved, prefix_rows = TRUE), resolved, parameters, store)
        })
    declaration$partitioned <- function(parameters, store) {
        object <- readWorkflowExecutionComponents(store$manager, store$source, fields = c("features", "design_matrix", "inputs"))
        resolved <- integrationResolveMofaParameters(object, parameters)
        aligned <- integrationAlignViewsPartitioned(object, resolved,
            integrationStateViewPages(store$manager, store$source, store$worker_cap), prefix_rows = TRUE)
        object <- NULL
        integrationModelFromAligned(aligned, resolved, parameters, store)
    }
    declaration
}

#' Train on the aligned data sets and save the product (both variants of the model declaration).
#' @param aligned The alignment, its rows named `<view>:<feature_id>`.
#' @param resolved The resolved parameters.
#' @param parameters The run's parameters (a stand-in `trainer` may be among them).
#' @param store The operation's store handle.
#' @noRd
integrationModelFromAligned <- function(aligned, resolved, parameters, store) {
    trainer <- parameters[["trainer", exact = TRUE]] %||% integrationTrainMofa
    fit <- trainer(aligned$views, integrationMofaOptions(resolved))
    reclaimWorkflowMemory()
    integrationModelProduct(aligned, fit, store)
}

#' The declaration a model run takes (KIT-063 item A): the registry's, or with a sample match its declaration with
#' `sample_match`. In process and on the worker alike.
#' @param parameters The resolved parameters.
#' @noRd
integrationMofaDeclaration <- function(parameters) {
    if (is.null(parameters[["sample_match", exact = TRUE]])) return(workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$declaration)
    integrationMofaOperation(match = TRUE)
}

#' The parameters an accepted model records: the declared ones, and the sample match with its table when one was given.
#' @noRd
integrationMofaRecordedNames <- function(parameters) {
    c(integrationMofaParameterNames(), intersect(c("sample_match", "sample_match_table"), names(parameters)))
}

#' Accept a model on the head (the RUVIII export's pattern): the product rooted, the step's provenance cleared (it
#' commits no state) and the model recorded under `snapshot$exports`, which the analysis record lists.
#' @param manager The integration owner, inside its state scope.
#' @param source The inputs state's reference.
#' @param run The operation's run (`reference`, `value`).
#' @param parameters The resolved declared parameters.
#' @param confirmation KIT-142: the sample table's confirmation as the record keeps it (`integrationConfirmationRecord`),
#'   or NULL; kept with the model. The same model accepted again keeps the results saved from it.
#' @noRd
acceptIntegrationModel <- function(manager, source, run, parameters, confirmation = NULL) {
    context <- manager$.context
    clearWorkflowStepProvenance(manager)
    if (is.null(context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(manager, run$reference, list(source))
    snapshot <- context$snapshot
    previous <- snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
    export <- list(operation = INTEGRATION_MODEL_OPERATION, source = source,
        reference = run$reference, parameters = parameters[integrationMofaRecordedNames(parameters)])
    if (!is.null(confirmation)) export$confirmation <- confirmation
    if (identical(previous$reference$entity_id, run$reference$entity_id) && !is.null(previous$results)) {
        export$results <- previous$results
    }
    # KIT-141: the enrichments of the same model are kept; a new model starts without any.
    if (identical(previous$reference$entity_id, run$reference$entity_id) && !is.null(previous$enrichment)) {
        export$enrichment <- previous$enrichment
    }
    snapshot$exports[[INTEGRATION_MODEL_OPERATION]] <- export
    publishWorkflowStateHead(context, snapshot)
    invisible(run)
}

#' Train the integration model on the owner's inputs state and accept it (KIT-062). Refused in plain words before any
#' training: no inputs state at the head, an input changed on disk since the import (`check_inputs`), parameters
#' that do not fit the state, MOFA2 or its Python not available (`integration_mofa_unavailable`, with the probe's
#' reason).
#' @param manager The integration owner.
#' @param parameters The declared parameters to set (absent ones take the provisional defaults). KIT-149: `variant`
#'   (`whole` or `partitioned`; absent, the registry's default variant) chooses how the data sets are aligned; it is
#'   not a parameter of the model and not part of its identity.
#' @param check_inputs Check each input's manifest on disk (the replay checks its inputs itself and passes FALSE).
#' @param probe The capability probe (`integrationMofaCapability`).
#' @param trainer NULL for MOFA2, or a stand-in with `integrationTrainMofa`'s form (tests; then nothing is probed).
#' @param confirmation KIT-142: the sample table's confirmation, kept with the model (`acceptIntegrationModel`).
#' @return The run (`status`, `value`, `reference`).
#' @noRd
runIntegrationModel <- function(manager, parameters = list(), check_inputs = TRUE, probe = integrationMofaCapability,
    trainer = NULL, confirmation = NULL) {
    if (!usesWorkflowReferenceStates(manager)) {
        cli::cli_abort("the integration model needs the integration owner", class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager, INTEGRATION_INPUTS_STATE)
        if (is.null(source) || !identical(manager$current_state, INTEGRATION_INPUTS_STATE)) {
            integrationRefuse("the integration inputs are not imported; import them first", "integration_inputs_missing")
        }
        object <- readWorkflowExecutionComponents(manager, source, fields = "inputs")
        if (isTRUE(check_inputs)) integrationAssertInputsCurrent(object)
        resolved <- integrationResolveMofaParameters(object, parameters)
        entry <- workflowOperationEntry(INTEGRATION_MODEL_OPERATION)
        variant <- parameters[["variant", exact = TRUE]] %||% entry$default_variant
        if (!(is.character(variant) && length(variant) == 1L && isTRUE(variant %in% entry$variants))) {
            integrationRefuse(sprintf("variant must be one of %s", paste(entry$variants, collapse = ", ")),
                "integration_parameters_invalid")
        }
        if (is.null(trainer)) {
            capability <- probe()
            if (!isTRUE(capability$available)) {
                integrationRefuse(sprintf("the integration model cannot be trained here: %s", capability$reason),
                    "integration_mofa_unavailable")
            }
        }
        declaration <- integrationMofaDeclaration(resolved)
        # KIT-148: a sample match is kept as a table of the store, and the model depends on it.
        match <- if (!is.null(resolved$sample_match_table)) list(integrationStoreSampleMatch(manager, resolved$sample_match_table))
        # A refused or failed training leaves no pending provenance for a later commit to take.
        run <- withCallingHandlers(runSavedOperation(manager, declaration, source,
            c(resolved, list(trainer = trainer, variant = variant)),
            dependencies = match %||% list()), error = function(error) clearWorkflowStepProvenance(manager))
        acceptIntegrationModel(manager, source, run, resolved, confirmation)
        run
    })
}

#' The model on the managed worker (KIT-140; the registry's `omics` handler for `integration_mofa_model`): the same
#' declaration `runIntegrationModel` runs, on the request's inputs state with its resolved parameters, MOFA2 probed in
#' the worker's own process (the worker may not reach the Python the app reaches). The product commits no state; the
#' app accepts it on the head (`acceptIntegrationModel`) when the action is delivered, as in process.
#' @param manager The worker's reopened integration store, inside its state scope.
#' @param request The worker request (`input_references$state`, `parameters`; a sample match as its table's reference,
#'   `integrationModelPrecheck()$request`).
#' @param variant KIT-149: the variant the request runs (`whole` or `partitioned`; the registry's default when the
#'   request names none), as in process.
#' @return The completion value: the source's checkpoint, the product's reference and the run's status (computed or hit).
#' @noRd
computeWorkflowIntegrationModelAction <- function(manager, request, variant) {
    source <- request$input_references$state
    object <- readWorkflowExecutionComponents(manager, source, fields = "inputs")
    # KIT-148: the request carries the sample match as the reference of its table in the store; it is read here and
    # checked against its sha256 (`sample_match`) as any match is.
    parameters <- request$parameters
    match <- parameters[["sample_match_table", exact = TRUE]]
    if (isSavedWorkflowReference(match)) parameters$sample_match_table <- integrationStoredSampleMatch(manager, match)
    resolved <- integrationResolveMofaParameters(object, parameters)
    capability <- integrationMofaCapability()
    if (!isTRUE(capability$available)) {
        integrationRefuse(sprintf("the integration model cannot be trained here: %s", capability$reason),
            "integration_mofa_unavailable")
    }
    run <- runSavedOperation(manager, integrationMofaDeclaration(resolved), source, c(resolved, list(variant = variant)),
        dependencies = if (isSavedWorkflowReference(match)) list(match) else list())
    # A product an earlier worker computed, whose action never completed, has no root yet (the omics handler's rule).
    if (identical(run$status, "hit")) rootWorkflowProduct(manager, run$reference, list(source))
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = run$reference,
        display = list(status = run$status))
}

#' The realised outcomes the analysis record lists for an accepted model: the views and their feature counts, the
#' samples, the factors, the training's length and last ELBO, the versions and the model file's sha256.
#' @param context The owner's context.
#' @param export The accepted export (`snapshot$exports$integration_mofa_model`).
#' @noRd
integrationModelRecordOutcomes <- function(context, export) {
    value <- tryCatch(hydrateSavedWorkflowEntity(context$store, export$reference)$value, error = function(error) NULL)
    if (!is.list(value)) return(list())
    list(views = value$view_summary$view, features_model = value$view_summary$features_model,
        internal_standards_dropped = value$view_summary$internal_standards_dropped, samples = nrow(value$samples),
        shared_samples = value$shared_samples, factors = ncol(value$factors) - 1L,
        iterations = value$elbo$iterations, elbo_last = value$elbo$last, mofa2 = value$versions$MOFA2 %||% NA_character_,
        mofapy2 = value$versions$mofapy2 %||% NA_character_, python = value$versions$python %||% NA_character_,
        model_sha256 = value$model_file$sha256 %||% NA_character_, notes = value$notes)
}

#' Replay an accepted model from the record's parameters on the replayed inputs state (the replay's analyses loop).
#' @param manager The replay's owner.
#' @param analysis The record's analysis entry.
#' @noRd
replayIntegrationModel <- function(manager, analysis) {
    # KIT-063 item A: a model trained on a user's sample match replays with that match.
    parameters <- workflowReplayParameters(analysis$parameters, c(workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$parameters,
        "sample_match", "sample_match_table"))
    runIntegrationModel(manager, parameters, check_inputs = FALSE)
}
