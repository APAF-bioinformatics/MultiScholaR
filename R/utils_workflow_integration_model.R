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

#' Read the omics' integration inputs into one `IntegrationViewsData` (KIT-062). Each folder is read and validated by
#' `readIntegrationInput()`; the omics are ordered by name, their views in manifest order. Refused in plain words: fewer
#' than two omics (D7), two inputs of one omic, a view name used by two omics (INT-4 names views `<omic>_<assay>`), a
#' biological-sample column the samples table lacks or leaves empty.
#' @param dirs The omics' folders, `<project>/integration/<omic>/`.
#' @return `list(object, manifests)`: the object, and the manifests' paths named by omic.
#' @noRd
integrationViewsFromInputs <- function(dirs) {
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
            views[[view]] <- input$views[[view]]$matrix
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
    object <- methods::new("IntegrationViewsData", views = views, features = features, design_matrix = design,
        inputs = record)
    list(object = object, manifests = stats::setNames(vapply(record, function(entry) entry$path, character(1)), omics))
}

#' Import the omics' integration inputs into the integration owner (KIT-062): the inputs state committed with the
#' provenance `integration_import` (a replay precondition), and each manifest recorded as the input role
#' `integration_input_<omic>` with its sha256.
#' @param manager The integration owner (`newIntegrationStateOwner`).
#' @param dirs The omics' folders, `<project>/integration/<omic>/`.
#' @param description The state's description.
#' @return The inputs state's reference, invisibly.
#' @noRd
importIntegrationInputs <- function(manager, dirs, description = "Integration inputs imported") {
    if (!usesWorkflowReferenceStates(manager)) {
        cli::cli_abort("the integration import needs the integration owner", class = "saved_entity_scope")
    }
    read <- integrationViewsFromInputs(dirs)
    object <- read$object
    config <- list(outcomes = list(omics = names(object@inputs), views = names(object@views),
        runs = nrow(object@design_matrix),
        biological_samples = length(unique(object@design_matrix$biological_sample))))
    withWorkflowStateScope(manager, function() {
        setWorkflowStepProvenance(manager, INTEGRATION_IMPORT_OPERATION)
        manager$saveState(INTEGRATION_INPUTS_STATE, object, config, description)
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
        seed = list(default = 20261119L, choices = NULL, settled_by = workbook))
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
    out[integrationMofaParameterNames()]
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
#' @return `list(views, samples, view_summary, shared_samples, notes)`: the views as features x samples matrices in
#'   the selected order with one column set, the aligned sample table, a summary row per view.
#' @noRd
integrationAlignViews <- function(object, parameters) {
    design <- object@design_matrix
    view_omic <- unlist(unname(lapply(object@inputs, function(entry) stats::setNames(rep(entry$omic, length(entry$views)),
        entry$views))))
    views <- parameters$views
    averaged <- list()
    runs <- list()
    standards <- integer()
    for (view in views) {
        omic <- view_omic[[view]]
        table <- object@views[[view]]
        values <- as.matrix(table[-1L])
        storage.mode(values) <- "double"
        rownames(values) <- as.character(table[[1L]])
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
    shared <- sort(Reduce(intersect, present), method = "radix")
    if (!length(shared)) {
        integrationRefuse(sprintf("the views %s share no biological sample; check the biological-sample column of each omic",
            paste(views, collapse = ", ")), "integration_no_shared_samples")
    }
    kept <- if (identical(parameters$missing_samples, "drop_incomplete")) shared else
        sort(unique(unlist(present)), method = "radix")
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
    summary <- list()
    final <- list()
    for (view in views) {
        m <- matrix(NA_real_, nrow(averaged[[view]]), length(kept), dimnames = list(rownames(averaged[[view]]), kept))
        here <- intersect(kept, colnames(averaged[[view]]))
        m[, here] <- averaged[[view]][, here, drop = FALSE]
        features_input <- nrow(m) + standards[[view]]
        sparse <- rowSums(!is.na(m)) < 2L
        m <- m[!sparse, , drop = FALSE]
        constant <- 0L
        if (identical(parameters$scaling, "feature_zscore")) {
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            flat <- !(spread > 0)
            constant <- sum(flat)
            m <- m[!flat, , drop = FALSE]
        }
        if (identical(parameters$feature_selection, "top_variance") && nrow(m) > parameters$features_per_view) {
            variance <- apply(m, 1L, stats::var, na.rm = TRUE)
            keep <- sort(order(-variance, seq_along(variance))[seq_len(parameters$features_per_view)])
            m <- m[keep, , drop = FALSE]
        }
        if (identical(parameters$scaling, "feature_zscore")) {
            centre <- rowMeans(m, na.rm = TRUE)
            spread <- apply(m, 1L, stats::sd, na.rm = TRUE)
            m <- (m - centre) / spread
        }
        if (!nrow(m)) {
            integrationRefuse(sprintf("view %s has no feature left for the model (none is observed in two or more samples)", view),
                "integration_view_empty")
        }
        final[[view]] <- m
        summary[[view]] <- data.frame(view = view, omic = view_omic[[view]], features_input = features_input,
            internal_standards_dropped = standards[[view]], sparse_dropped = sum(sparse), constant_dropped = constant,
            features_model = nrow(m), samples_present = length(here), stringsAsFactors = FALSE)
    }
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
    list(views = final, samples = samples, view_summary = summary, shared_samples = length(shared), notes = notes)
}

#' The MOFA2 options the declared parameters set (every other MOFA2 option stays at its default).
#' @noRd
integrationMofaOptions <- function(parameters) {
    parameters[c("scale_views", "center_groups", "num_factors", "convergence_mode", "maxiter", "seed")]
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

#' Train MOFA2 as the proven lane does (KIT-060 status 2026-11-19, run 44423980807): `create_mofa` on the views,
#' the default data, model and training options with the declared ones set, `prepare_mofa`, `run_mofa(use_basilisk =
#' FALSE)` (a fall-back to basilisk is refused: the bundle's Python must be the one used), then the trained model
#' loaded from its file.
#' @param views Named list of features x samples matrices (feature names unique across views).
#' @param options `integrationMofaOptions()`.
#' @return `list(factors, weights, r2_per_factor (factors x views), r2_total, elbo, iterations, options, versions,
#'   messages, model_file)`.
#' @noRd
integrationTrainMofa <- function(views, options) {
    object <- MOFA2::create_mofa(views)
    data_options <- MOFA2::get_default_data_options(object)
    data_options$scale_views <- options$scale_views
    data_options$center_groups <- options$center_groups
    model_options <- MOFA2::get_default_model_options(object)
    model_options$num_factors <- options$num_factors
    training_options <- MOFA2::get_default_training_options(object)
    training_options$seed <- options$seed
    training_options$maxiter <- options$maxiter
    training_options$convergence_mode <- options$convergence_mode
    training_options$verbose <- FALSE
    object <- MOFA2::prepare_mofa(object, data_options = data_options, model_options = model_options,
        training_options = training_options)
    outfile <- tempfile("integration-mofa-", fileext = ".hdf5")
    on.exit(unlink(outfile), add = TRUE)
    # Kept, not dropped: MOFA2's messages and warnings are returned with the model.
    messages <- character()
    withCallingHandlers(
        MOFA2::run_mofa(object, outfile = outfile, save_data = TRUE, use_basilisk = FALSE),
        message = function(m) {
            messages <<- c(messages, trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            messages <<- c(messages, paste("WARNING:", trimws(conditionMessage(w))))
            invokeRestart("muffleWarning")
        })
    if (any(grepl("Setting use_basilisk = TRUE", messages, fixed = TRUE))) {
        integrationRefuse("MOFA2 could not use the configured Python with mofapy2 and fell back to basilisk; the model is not trained",
            "integration_mofa_unavailable")
    }
    model <- MOFA2::load_model(outfile)
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

#' The model's saved product: the factors, the weights, the variance explained per view and factor, the ELBO trace
#' summary, the versions and options, the aligned sample table and the view summary (KIT-062).
#' @param aligned `integrationAlignViews()`'s value.
#' @param fit The trainer's value (`integrationTrainMofa()`'s form), on views whose features are `<view>:<feature_id>`.
#' @noRd
integrationModelProduct <- function(aligned, fit) {
    factors <- fit$factors
    factor_table <- data.frame(biological_sample = rownames(factors), factors, check.names = FALSE,
        stringsAsFactors = FALSE, row.names = NULL)
    weights <- do.call(rbind, lapply(names(aligned$views), function(view) {
        w <- fit$weights[[view]]
        data.frame(view = view, feature_id = substring(rownames(w), nchar(view) + 2L), w, check.names = FALSE,
            stringsAsFactors = FALSE, row.names = NULL)
    }))
    r2 <- fit$r2_per_factor
    variance <- data.frame(factor = rep(rownames(r2), times = ncol(r2)), view = rep(colnames(r2), each = nrow(r2)),
        r2_percent = as.vector(r2), stringsAsFactors = FALSE)
    elbo <- fit$elbo[is.finite(fit$elbo)]
    list(schema = "multischolar.integration_mofa_model.v1", factors = factor_table, weights = weights,
        variance_explained = variance,
        variance_total = data.frame(view = names(fit$r2_total), r2_percent = unname(as.numeric(fit$r2_total)),
            stringsAsFactors = FALSE),
        elbo = list(iterations = as.integer(fit$iterations), evaluations = length(elbo),
            first = if (length(elbo)) elbo[[1L]] else NA_real_, last = if (length(elbo)) elbo[[length(elbo)]] else NA_real_,
            trace = elbo),
        samples = aligned$samples, view_summary = aligned$view_summary, shared_samples = aligned$shared_samples,
        notes = aligned$notes, options = fit$options, versions = fit$versions, messages = fit$messages,
        model_file = fit$model_file)
}

#' The MOFA2 model declaration (KIT-062): whole matrix, on the inputs state, product only. Its identity is the state's
#' (every view's table, the design), the manifests' sha256, the views and every method parameter. A caller may hand
#' the kernel a `trainer` beside the declared parameters (not identity; the tests' stand-in for MOFA2).
#' @noRd
integrationMofaOperation <- function() {
    newWorkflowOperation(INTEGRATION_MODEL_OPERATION, "IntegrationViewsData",
        fields = c("views", "features", "design_matrix", "inputs"), unit = "whole_matrix", output_class = NULL,
        parameters = integrationMofaParameterNames(), product_type = "integration.mofa_model",
        kernel = function(object, parameters, store) {
            resolved <- integrationResolveMofaParameters(object, parameters)
            aligned <- integrationAlignViews(object, resolved)
            trainer <- parameters[["trainer", exact = TRUE]] %||% integrationTrainMofa
            prefixed <- lapply(names(aligned$views), function(view) {
                m <- aligned$views[[view]]
                rownames(m) <- paste(view, rownames(m), sep = ":")
                m
            })
            names(prefixed) <- names(aligned$views)
            integrationModelProduct(aligned, trainer(prefixed, integrationMofaOptions(resolved)))
        })
}

#' Accept a model on the head (the RUVIII export's pattern): the product rooted, the step's provenance cleared (it
#' commits no state) and the model recorded under `snapshot$exports`, which the analysis record lists.
#' @param manager The integration owner, inside its state scope.
#' @param source The inputs state's reference.
#' @param run The operation's run (`reference`, `value`).
#' @param parameters The resolved declared parameters.
#' @noRd
acceptIntegrationModel <- function(manager, source, run, parameters) {
    context <- manager$.context
    clearWorkflowStepProvenance(manager)
    if (is.null(context$snapshot$products[[run$reference$entity_id]])) rootWorkflowProduct(manager, run$reference, list(source))
    snapshot <- context$snapshot
    snapshot$exports[[INTEGRATION_MODEL_OPERATION]] <- list(operation = INTEGRATION_MODEL_OPERATION, source = source,
        reference = run$reference, parameters = parameters[integrationMofaParameterNames()])
    publishWorkflowStateHead(context, snapshot)
    invisible(run)
}

#' Train the integration model on the owner's inputs state and accept it (KIT-062). Refused in plain words before any
#' training: no inputs state at the head, an input changed on disk since the import (`check_inputs`), parameters
#' that do not fit the state, MOFA2 or its Python not available (`integration_mofa_unavailable`, with the probe's
#' reason).
#' @param manager The integration owner.
#' @param parameters The declared parameters to set (absent ones take the provisional defaults).
#' @param check_inputs Check each input's manifest on disk (the replay checks its inputs itself and passes FALSE).
#' @param probe The capability probe (`integrationMofaCapability`).
#' @param trainer NULL for MOFA2, or a stand-in with `integrationTrainMofa`'s form (tests; then nothing is probed).
#' @return The run (`status`, `value`, `reference`).
#' @noRd
runIntegrationModel <- function(manager, parameters = list(), check_inputs = TRUE, probe = integrationMofaCapability,
    trainer = NULL) {
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
        if (is.null(trainer)) {
            capability <- probe()
            if (!isTRUE(capability$available)) {
                integrationRefuse(sprintf("the integration model cannot be trained here: %s", capability$reason),
                    "integration_mofa_unavailable")
            }
        }
        declaration <- workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$declaration
        # A refused or failed training leaves no pending provenance for a later commit to take.
        run <- withCallingHandlers(runSavedOperation(manager, declaration, source, c(resolved, list(trainer = trainer))),
            error = function(error) clearWorkflowStepProvenance(manager))
        acceptIntegrationModel(manager, source, run, resolved)
        run
    })
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
    parameters <- workflowReplayParameters(analysis$parameters, workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$parameters)
    runIntegrationModel(manager, parameters, check_inputs = FALSE)
}
