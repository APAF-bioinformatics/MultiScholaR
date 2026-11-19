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
            integrationRefuse(sprintf("the integration input of %s leaves a run without its biological sample (column %s)",
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
