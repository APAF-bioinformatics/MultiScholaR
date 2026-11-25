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
#'   `biological_sample_column`; KIT-184: and, only for an omic whose input holds them, `removed_factors`, named by data
#'   set: the manifest's record of the unwanted factors its correction removed, the file's `sha256` and `factors` (runs
#'   by factors).
#' @slot sample_id,group_id,biological_sample_id The design's run, group and biological-sample columns.
#' @slot args Unused; the declared classes' argument slot.
#' @slot corrected_views KIT-153: named list by view of the RUV-corrected tables of the views whose omic saved one (the
#'   same form as `views`); `views` holds each view's normalised matrix before RUV.
#' @noRd
setClass("IntegrationViewsData",
    slots = c(views = "list", features = "list", design_matrix = "data.frame", inputs = "list",
        sample_id = "character", group_id = "character", biological_sample_id = "character", args = "list",
        corrected_views = "list"),
    prototype = list(views = list(), features = list(), design_matrix = data.frame(), inputs = list(),
        sample_id = "run", group_id = "group", biological_sample_id = "biological_sample", args = list(),
        corrected_views = list()),
    validity = function(object) {
        if (!identical(names(object@views), names(object@features))) return("views and features must name the same views")
        if (!all(names(object@corrected_views) %in% names(object@views))) return("every corrected view must be a view")
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
            views[[view]] <- list(dir = dirs[[i]], omic = omic, corrected = isTRUE(input$views[[view]]$corrected))
            features[[view]] <- input$views[[view]]$features
        }
        design[[i]] <- data.frame(omic = omic, run = as.character(samples[[sample_column]]),
            group = as.character(samples[[manifest$group_column]]), biological_sample = key, stringsAsFactors = FALSE)
        path <- file.path(dirs[[i]], INTEGRATION_INPUT_MANIFEST)
        record[[omic]] <- list(omic = omic, path = normalizePath(path, mustWork = TRUE),
            manifest_sha256 = digest::digest(file = path, algo = "sha256"), views = names(input$views),
            manifest = manifest, samples = samples, biological_sample_column = biological)
        # KIT-184: the removed unwanted factors of each data set that has them (record, sha256, runs by factors). The
        # element is added only then, so an input without them gives the entry it gave before.
        for (view in names(input$views)) {
            removed <- input$views[[view]]$removed
            if (is.null(removed)) next
            record[[omic]]$removed_factors[[view]] <- c(removed,
                list(factors = readIntegrationInputRemovedFactors(dirs[[i]], view, input)))
        }
    }
    design <- do.call(rbind, design)
    rownames(design) <- NULL
    list(inputs = stats::setNames(inputs, omics), views = views, features = features, design = design, record = record,
        manifests = stats::setNames(vapply(record, function(entry) entry$path, character(1)), omics))
}

#' The reader of each view of the read inputs (`openIntegrationInputView`), in the inputs' view order.
#' @param read `integrationReadInputs()`'s value.
#' @param matrix KIT-153: `uncorrected` (every view's normalised matrix before RUV) or `ruv_corrected` (the corrected
#'   matrix of each view that has one).
#' @noRd
integrationInputViewReaders <- function(read, matrix = "uncorrected") {
    views <- names(read$views)
    if (identical(matrix, "ruv_corrected")) views <- views[vapply(read$views[views], function(entry) isTRUE(entry$corrected), logical(1))]
    readers <- lapply(views, function(view) {
        entry <- read$views[[view]]
        openIntegrationInputView(entry$dir, view, read$inputs[[entry$omic]], matrix)
    })
    stats::setNames(readers, views)
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
    corrected <- lapply(integrationInputViewReaders(read, "ruv_corrected"), function(reader) reader$read())
    object <- methods::new("IntegrationViewsData", views = matrices, features = read$features, design_matrix = read$design,
        inputs = read$record, corrected_views = corrected)
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
        # KIT-153: the RUV-corrected matrices are imported too (where an omic saved one), so that the model can take
        # either without a second import; each is a table of the store as the normalised ones are.
        corrected_readers <- integrationInputViewReaders(read, "ruv_corrected")
        corrected_tables <- lapply(corrected_readers, function(reader) integrationImportViewTable(context, reader))
        object <- methods::new("IntegrationViewsData", views = lapply(readers, function(reader) reader$prototype),
            features = read$features, design_matrix = read$design, inputs = read$record,
            corrected_views = lapply(corrected_readers, function(reader) reader$prototype))
        setWorkflowStepProvenance(manager, INTEGRATION_IMPORT_OPERATION)
        checkpoint <- checkpointWorkflowExecutionResult(manager, object)
        checkpoint$components$views <- tables
        checkpoint$components$corrected_views <- corrected_tables
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

#' The model's method parameters with their defaults, each with the decision that settled it (`settled_by`: the
#' methods-ledger row). The first values were the ones the MOFA2 smoke lane trained with on the tier-1 corpus (run
#' 44423980807: MOFA2 1.22.0, mofapy2 0.7.3; every other MOFA2 option at its default); KIT-240: those not decided
#' earlier were kept by the owner's decision VAL-221 (2026-11-25), the values unchanged.
#' @return Named list by parameter: `default`, `choices` (NULL for a number or a flag), `settled_by`.
#' @noRd
integrationMofaParameterSpec <- function() {
    # KIT-240 (VAL-221 (4)): the values kept, the decision named.
    kept <- "VAL-221 (decided 2026-11-25: %s)"
    list(
        # KIT-153 (INT-23, user 2026-11-21): each data set's normalised matrix before RUV by default; the RUV-corrected
        # matrix selectable where every data set in the model has one.
        input_matrix = list(default = "uncorrected", choices = c("uncorrected", "ruv_corrected"),
            settled_by = "INT-23 (decided 2026-11-21: the normalised matrix before RUV; the RUV-corrected matrix selectable)"),
        technical_replicates = list(default = "average", choices = "average",
            settled_by = "INT-2 (decided 2026-11-19: averaged per biological sample)"),
        internal_standards = list(default = "drop", choices = c("drop", "keep"),
            settled_by = "INT-3 (decided 2026-11-19: dropped)"),
        missing_samples = list(default = "keep_missing", choices = c("keep_missing", "drop_incomplete"),
            settled_by = "INT-010 (KIT-060 D4, decided 2026-11-17: kept as missing in that data set)"),
        # INT-12 (user, 2026-11-21): decided. The features whose variance lies furthest above the data set's fitted
        # mean-variance trend (`integrationFeatureSelectionScore`), 1,000 per data set (a smaller data set keeps all),
        # among the features observed in at least half the biological samples; raw variance and none stay selectable.
        feature_selection = list(default = "variance_above_trend", choices = c("variance_above_trend", "top_variance", "none"),
            settled_by = "INT-12 (decided 2026-11-21: variance above the fitted mean-variance trend)"),
        features_per_view = list(default = 1000L, choices = NULL,
            settled_by = "INT-12 (decided 2026-11-21: about 1,000 per data set; a smaller data set keeps all)"),
        min_observed_fraction = list(default = 0.5, choices = NULL, type = "number",
            settled_by = "INT-12 (decided 2026-11-21: observed in at least half the biological samples)"),
        # INT-24 (user, 2026-11-21): decided. Each selected feature scaled to mean 0 and standard deviation 1 (after the
        # selection, which is computed on the unscaled values), and each data set divided by its overall standard
        # deviation; the z-score stays selectable off.
        scaling = list(default = "feature_zscore", choices = c("none", "feature_zscore"),
            settled_by = "INT-24 (decided 2026-11-21: each selected feature scaled to mean 0 and standard deviation 1)"),
        scale_views = list(default = TRUE, choices = NULL,
            settled_by = "INT-24 (decided 2026-11-21: each data set divided by its overall standard deviation)"),
        # MOFA2's group centring: the app passes no sample groups (`create_mofa(views)`), so it centres each feature over
        # all samples.
        center_groups = list(default = TRUE, choices = NULL,
            settled_by = sprintf(kept, "TRUE kept; with one sample group, each feature is centred over all samples")),
        num_factors = list(default = 5L, choices = NULL, settled_by = sprintf(kept, "5 factors")),
        # INT-15 (user, 2026-11-21): decided. Training stops by its convergence rule well within 10,000 rounds (none of
        # the 1,020 slow fits of KIT-151 run b passed 3,486); whether it stopped by the rule or at the limit is kept.
        convergence_mode = list(default = "slow", choices = c("fast", "medium", "slow"),
            settled_by = "INT-15 (decided 2026-11-21: slow)"),
        maxiter = list(default = 10000L, choices = NULL,
            settled_by = "INT-15 (decided 2026-11-21: a limit training stops well within)"),
        seed = list(default = 20261119L, choices = NULL, settled_by = sprintf(kept, "20261119, the seed of the reported model")),
        # INT-14 (user, 2026-11-21): the model is trained from several seeds and each factor's agreement across them
        # reported (`integrationSeedStability`); the reported model is the declared seed's. VAL-221 (2): five seeds.
        training_seeds = list(default = 5L, choices = NULL,
            settled_by = paste("INT-14 (several seeds, 2026-11-21);", sprintf(kept, "5 seeds, the declared seed and the next four"))),
        # INT-6 (user, 2026-11-19): decided, not provisional. Double precision (KIT-062 diagnosis run 44462198417:
        # MOFA2's float32 default followed one of two training paths depending on the CPU model and thread count);
        # single precision uses less memory for very large studies, and its results differ more between computers.
        # INT-17 (2026-11-21) corrects INT-6: double precision gives the identical model on the same computer and
        # software, not on every computer; with the maths library pinned (`integrationPinMathsLibrary`) models were
        # bit-identical across the AMD processors tested; without it they agree to about 1e-12 and a factor near the
        # detection limit can differ between processor models.
        precision = list(default = "double", choices = c("double", "single"),
            settled_by = paste("INT-6 (decided 2026-11-19: double), corrected by INT-17 (2026-11-21: identical on the same",
                "computer and software, not guaranteed across processors)")))
}

#' The method-parameter table (one row per method parameter, with the decision behind each default), for the ticket,
#' the ledger and the GUI.
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

#' The data sets of the inputs state whose omic saved a RUV-corrected matrix (KIT-153), from the imported manifests.
#' @param object An `IntegrationViewsData` with its `inputs` hydrated.
#' @return A character vector of view names.
#' @noRd
integrationCorrectedViewNames <- function(object) {
    unlist(lapply(unname(object@inputs), function(entry) {
        named <- vapply(entry$manifest$views, function(view) if (!is.null(view$corrected_matrix_file))
            as.character(view$name) else NA_character_, character(1))
        named[!is.na(named)]
    }), use.names = FALSE) %||% character()
}

#' The slot of the inputs state the model reads its data sets from (KIT-153): `views` (the normalised matrices before
#' RUV) or, with `input_matrix = ruv_corrected`, `corrected_views`.
#' @param parameters The resolved parameters.
#' @noRd
integrationModelViewsSlot <- function(parameters) {
    if (identical(parameters$input_matrix, "ruv_corrected")) "corrected_views" else "views"
}

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
    # KIT-151 (plan item E): centring is required unless every feature is already centred by its z-score.
    if (!isTRUE(out$center_groups) && !identical(out$scaling, "feature_zscore")) {
        invalid("center_groups can be FALSE only with scaling = feature_zscore (the features are otherwise not centred)")
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
    out$training_seeds <- whole("training_seeds", 1L)
    # The other seeds are the next whole numbers after the declared one (`integrationTrainingSeeds`).
    if (out$seed > .Machine$integer.max - (out$training_seeds - 1L)) {
        invalid(sprintf("seed plus training_seeds - 1 must be at most %d", .Machine$integer.max))
    }
    # INT-12: the number kept and the observed share apply only when features are selected; with no selection a value
    # the caller gives is refused and an absent one is NA (so the identity does not carry an unused default).
    if (!identical(out$feature_selection, "none")) {
        out$features_per_view <- whole("features_per_view", 1L)
        fraction <- out$min_observed_fraction
        if (!(is.numeric(fraction) && length(fraction) == 1L && !is.na(fraction) && fraction > 0 && fraction <= 1)) {
            invalid("min_observed_fraction must be a number above 0 and at most 1")
        }
        out$min_observed_fraction <- as.numeric(fraction)
    } else {
        for (name in c("features_per_view", "min_observed_fraction")) {
            value <- given(name)
            if (!is.null(value) && !(length(value) == 1L && is.na(value))) {
                invalid(sprintf("%s is set only with a feature selection (feature_selection other than none)", name))
            }
        }
        out$features_per_view <- NA_integer_
        out$min_observed_fraction <- NA_real_
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
    # KIT-153: the RUV-corrected matrix only where every selected data set has one.
    if (identical(out$input_matrix, "ruv_corrected")) {
        lacking <- setdiff(views, integrationCorrectedViewNames(object))
        if (length(lacking)) {
            integrationRefuse(sprintf(paste("the RUV-corrected matrix was chosen, but data set(s) %s have none (their omic",
                "applied no correction for unwanted variation); choose the uncorrected matrix, or leave them out"),
                paste(lacking, collapse = ", ")), c("integration_input_matrix_unavailable", "integration_parameters_invalid"))
        }
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
#' features observed in fewer than two samples dropped, or, when features are selected, in fewer than
#' `min_observed_fraction` of the kept samples (INT-12; `integrationMinObserved`), and, under per-feature scaling,
#' features without variance; then, in a view with more than `features_per_view` features, the `features_per_view`
#' with the largest selection score (computed on the values before any per-feature scaling; KIT-151, INT-24) (`feature_selection`: `variance_above_trend` or `top_variance`;
#' `integrationFeatureSelectionScore`), then the per-feature Z-score (`scaling = feature_zscore`). Refused: no
#' biological sample shared by every view, a biological sample with two groups, a view left without features.
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
    min_observed <- integrationMinObserved(parameters, length(plan$kept))
    selection <- parameters$feature_selection
    summary <- list()
    final <- list()
    for (view in views) {
        block <- integrationAlignedRows(averaged[[view]], plan$kept, zscore, min_observed)
        m <- block$values
        if (!identical(selection, "none") && nrow(m) > parameters$features_per_view) {
            score <- integrationFeatureSelectionScore(integrationFeatureSelectionStatistic(m), selection)
            m <- m[integrationSelectTopFeatures(score, parameters$features_per_view), , drop = FALSE]
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
# (each feature's average, its count of observed samples, its spread, its mean and variance and Z-score depend on
# that feature's values alone); only the feature selection compares features (KIT-151: the mean-variance trend is
# fitted over them), and it is taken once over the whole data set from the per-feature mean and variance.
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

#' The number of kept biological samples a feature must be observed in (KIT-062; INT-12): two without a feature
#' selection; with one, `min_observed_fraction` of the kept samples rounded up (the product first rounded to 8 decimals,
#' so that 0.3 x 10 is 3, not 4), and never fewer than two. Half of 20 samples is 10, half of 35 is 18.
#' @param parameters Resolved parameters.
#' @param samples The number of kept biological samples.
#' @noRd
integrationMinObserved <- function(parameters, samples) {
    if (identical(parameters$feature_selection, "none")) return(2L)
    max(2L, as.integer(ceiling(round(parameters$min_observed_fraction * samples, 8L))))
}

#' Rows of a data set on the kept biological samples (KIT-062's rules): placed on the kept samples (missing where the
#' data set lacks one), the features observed in fewer than `min_observed` samples dropped and, under per-feature
#' scaling, those without variance. Row by row: a block of rows gives the rows the whole data set gives.
#' @param averaged Features x biological samples (`integrationAverageRuns`).
#' @param kept The kept biological samples.
#' @param zscore Whether features without variance are dropped (per-feature scaling).
#' @param min_observed The fewest kept samples a feature is observed in (`integrationMinObserved`).
#' @return `list(values, sparse, constant, present)`: the rows left, the counts dropped as sparse and as constant, and
#'   how many kept samples the data set measures.
#' @noRd
integrationAlignedRows <- function(averaged, kept, zscore, min_observed = 2L) {
    m <- matrix(NA_real_, nrow(averaged), length(kept), dimnames = list(rownames(averaged), kept))
    here <- intersect(kept, colnames(averaged))
    m[, here] <- averaged[, here, drop = FALSE]
    sparse <- rowSums(!is.na(m)) < min_observed
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

#' The per-feature statistics the feature selection is computed from (KIT-149; INT-12): each feature's mean and
#' variance over the samples it is observed in. Row by row, so the partitioned alignment gathers them page by page in
#' its first pass and needs no further pass over the values; the one definition both alignments use.
#' @param m Features x samples.
#' @return `list(mean, variance)`: one number per feature each, in row order.
#' @noRd
integrationFeatureSelectionStatistic <- function(m) {
    list(mean = as.numeric(rowMeans(m, na.rm = TRUE)), variance = as.numeric(apply(m, 1L, stats::var, na.rm = TRUE)))
}

#' The span of the mean-variance trend's fit (INT-12): the evidence run's (KIT-151 run b, `k151b.R`, TREND_SPAN).
#' @noRd
INTEGRATION_TREND_SPAN <- 0.3

#' The fitted mean-variance trend (INT-12), as in the evidence run (KIT-151 run b, `.github/mofa-smoke/k151b.R`,
#' `fitted_trend`): a lowess fit of `y` on the features' mean (`stats::lowess`, span 0.3, its default three robustness
#' iterations and `delta`), read at each feature's mean by linear interpolation (`stats::approx`, equal means averaged,
#' the ends held constant). With fewer than two distinct means there is no trend: 0 for every feature.
#' @param abundance Each feature's mean.
#' @param y Each feature's variance.
#' @return The fitted value per feature.
#' @noRd
integrationFittedTrend <- function(abundance, y) {
    if (length(unique(abundance)) < 2L) return(rep(0, length(y)))
    fit <- stats::lowess(abundance, y, f = INTEGRATION_TREND_SPAN)
    stats::approx(fit$x, fit$y, xout = abundance, ties = base::mean, rule = 2)$y
}

#' The score the feature selection ranks a data set's features by (INT-12): `variance_above_trend`, each feature's
#' variance minus the fitted mean-variance trend at its mean (`integrationFittedTrend`; the "difference" form of the
#' evidence run); `top_variance`, its variance.
#' @param statistic `integrationFeatureSelectionStatistic()` of every feature of the data set, in its row order.
#' @param selection The `feature_selection` parameter.
#' @return One number per feature, in row order.
#' @noRd
integrationFeatureSelectionScore <- function(statistic, selection) {
    if (identical(selection, "top_variance")) return(statistic$variance)
    statistic$variance - integrationFittedTrend(statistic$mean, statistic$variance)
}

#' The features the selection keeps: the `n` largest scores, a tie broken by position (the earlier feature first),
#' returned in position order.
#' @param statistic `integrationFeatureSelectionScore()` of every feature of a data set, in its row order.
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
    integrationRefuse(sprintf("view %s has no feature left for the model (none is observed in enough biological samples)", view),
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
#' left after the sparse and constant filters, its selection statistics (its mean and variance; KIT-151); the selection
#' is taken once over the whole data set from those numbers (`integrationFeatureSelectionScore`, the trend fitted once,
#' then `integrationSelectTopFeatures`, ties by position); pass 2 reads the pages again and keeps the selected rows, Z-scored under per-feature scaling, in a matrix
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
    selection <- parameters$feature_selection
    cap <- if (!identical(selection, "none")) parameters$features_per_view else NA_integer_
    min_observed <- integrationMinObserved(parameters, length(plan$kept))
    summary <- list()
    final <- if (isTRUE(values)) list()
    for (view in views) {
        # One page of the data set on the kept samples: the whole path's read, standards and averaging, row by row.
        block_of <- function(page, offset) {
            m <- as.matrix(page[-1L])
            storage.mode(m) <- "double"
            rownames(m) <- if (isTRUE(prefix_rows)) paste(view, page[[1L]], sep = ":") else as.character(page[[1L]])
            if (drop_standards) m <- m[!standards[[view]][offset + seq_len(nrow(m))], , drop = FALSE]
            # KIT-155: a page of internal standards only leaves no row (the averaging needs at least one).
            if (!nrow(m)) {
                return(list(values = matrix(NA_real_, 0L, length(plan$kept), dimnames = list(character(), plan$kept)),
                    sparse = 0L, constant = 0L))
            }
            integrationAlignedRows(integrationAverageRuns(m, keys[[view]]), plan$kept, zscore, min_observed)
        }
        # Pass 1: the counts, and the selection statistics (mean, variance) of every feature left.
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
            if (!is.na(cap)) statistic[[length(statistic) + 1L]] <<- integrationFeatureSelectionStatistic(block$values)
            invisible(NULL)
        })
        # The global statistic: the selection over the whole data set from the per-feature statistics in row order
        # (the trend is fitted once, on every feature left; NULL: every feature left is kept).
        take <- NULL
        if (!is.na(cap) && left > cap) {
            gathered <- list(mean = unlist(lapply(statistic, `[[`, "mean")),
                variance = unlist(lapply(statistic, `[[`, "variance")))
            take <- logical(left)
            take[integrationSelectTopFeatures(integrationFeatureSelectionScore(gathered, selection), cap)] <- TRUE
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

#' The byte budget of one page of a data set read in row pages (KIT-155; INT-29): 8 MiB of values. KIT-150 (INT-27)
#' measured the training stage's R heap growing by about 13 MiB per MiB of a page (pages of 10,000 rows: 15 MiB at
#' 200 samples, 38 MiB at 500), the page's copies in the alignment; at 8 MiB a page adds at most about 105 MiB
#' whatever the number of samples, which keeps the growth above the stage's fixed overhead under claim 1's bound at
#' KIT-150's sizes. A table of at most 104 columns keeps its pages of 10,000 rows (the read ceiling).
#' @noRd
INTEGRATION_PAGE_BYTES <- 8388608

#' Rows per page of a data set read in row pages (KIT-155): the budget over 8 bytes per value of every column read,
#' at least 1 and at most 10,000 (the store's read ceiling).
#' @param columns The number of columns read (the feature identifier and the runs).
#' @param page_bytes The byte budget (`INTEGRATION_PAGE_BYTES`).
#' @noRd
integrationPageRows <- function(columns, page_bytes = INTEGRATION_PAGE_BYTES) {
    as.integer(min(10000, max(1, floor(page_bytes / (8 * columns)))))
}

#' Pages of the inputs state's data sets from the integration store (KIT-149; `integrationAlignViewsPartitioned`'s
#' `pages`): each data set's table read by `partitionByRows` in pages sized by a byte budget (KIT-155;
#' `integrationPageRows`), a table of more than 64 columns read in blocks of 64 columns bound in column order.
#' @param manager The integration owner.
#' @param source The inputs state's reference.
#' @param worker_cap The admitted worker memory (NULL: 12 GiB, the omics' partitioned steps' default).
#' @param slot KIT-153: `views` (the normalised matrices) or `corrected_views` (`integrationModelViewsSlot`).
#' @param page_bytes KIT-155: the byte budget of a page (`INTEGRATION_PAGE_BYTES`; tests pass others).
#' @noRd
integrationStateViewPages <- function(manager, source, worker_cap = NULL, slot = "views",
    page_bytes = INTEGRATION_PAGE_BYTES) {
    function(view, columns, kernel) withWorkflowStateScope(manager, function() {
        context <- manager$.context
        table <- readWorkflowStateMetadata(manager, source)$data$components[[slot]][[view]]
        scope <- newWorkflowTableReadScope(context$store, table, reader = if (length(columns) > 64L) "duckdb" else "arrow")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        partitionByRows(scope, columns, kernel = function(page, index, offset) kernel(page, offset),
            consume = function(value, index) invisible(NULL), worker_cap = worker_cap %||% (12 * 1024^3),
            working_bytes_per_row = omicsPartitionedWorkingBytesPerColumn * length(columns),
            overhead_bytes = 2 * 67108864 + 2 * 8388608, policy = "integration-row-pages-v1",
            page_rows = integrationPageRows(length(columns), page_bytes))
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
#' @return `list(available, reason)`; KIT-240: when not available, also `missing` (`"MOFA2"`, `"reticulate"` or
#'   `"mofapy2"`), the piece the Model tab names.
#' @noRd
integrationMofaCapability <- function() {
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
        return(list(available = FALSE, reason = "the R package MOFA2 is not installed", missing = "MOFA2"))
    }
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        return(list(available = FALSE, reason = "the R package reticulate is not installed", missing = "reticulate"))
    }
    # KIT-151 (INT-17): the maths library is pinned before anything starts Python in this process.
    integrationPinMathsLibrary()
    found <- tryCatch(reticulate::py_module_available("mofapy2"), error = function(error) FALSE)
    if (!isTRUE(found)) {
        python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
        return(list(available = FALSE, reason = sprintf("no Python with the module mofapy2 is available to reticulate%s",
            if (nzchar(python)) sprintf(" (RETICULATE_PYTHON is %s)", python) else " (RETICULATE_PYTHON is not set)"),
            missing = "mofapy2"))
    }
    list(available = TRUE, reason = NA_character_)
}

#' The maths-library pin for training (KIT-151; INT-17, audit SCI-05): the variables the KIT-151 c run set before
#' Python started, under which every model was bit-identical across the processors seen (run 44678183859: AMD EPYC
#' 7763, 9V74, 9V45; an Intel processor not yet tested): the OpenBLAS copies of numpy and scipy on one thread and one
#' code path (Haswell: AVX2 and FMA, which every hosted-runner processor seen supports), OpenMP and MKL (if loaded) on
#' one thread. In the evidence run R's own OpenBLAS was pinned too (set before R started); in the app R has loaded its
#' maths library before any of this runs, so only Python's libraries are pinned, where MOFA2 trains.
#' @noRd
INTEGRATION_MATHS_PIN <- c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    OPENBLAS_CORETYPE = "Haswell")

#' Pin the maths library of the training's Python (KIT-151; INT-17), once per process, in the app's process and on the
#' worker alike: when Python has not loaded numpy yet, `INTEGRATION_MATHS_PIN` is set, Python started and mofapy2's
#' entry point (numpy, scipy, scikit-learn) loaded, then the variables are restored, so that only those libraries,
#' which read them when they load, run pinned and the rest of the R session is unchanged. When numpy was already loaded
#' in this process (Python started earlier by other code), its libraries keep the settings they loaded with: the pin
#' cannot be applied, training goes ahead unpinned, and the model's record says so (`applied = FALSE` with the reason);
#' such a model is still identical on the same computer and software, agrees to about 1e-12 across processors, and a
#' factor near the detection limit can differ between processor models.
#' @return `list(applied, variables, reason)`.
#' @noRd
integrationPinMathsLibrary <- local({
    state <- NULL
    function() {
        if (!is.null(state)) return(state)
        variables <- as.list(INTEGRATION_MATHS_PIN)
        if (isTRUE(reticulate::py_available(initialize = FALSE))) {
            loaded <- tryCatch(isTRUE(reticulate::py_eval("'numpy' in __import__('sys').modules")), error = function(error) TRUE)
            if (loaded) {
                state <<- list(applied = FALSE, variables = variables, reason = paste("numpy was already loaded in this",
                    "process before the integration model was first trained, so its maths library could not be pinned"))
                return(state)
            }
        }
        previous <- Sys.getenv(names(INTEGRATION_MATHS_PIN), unset = NA_character_, names = TRUE)
        on.exit({
            for (name in names(previous)) {
                if (is.na(previous[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(previous[[name]]), name))
            }
        }, add = TRUE)
        do.call(Sys.setenv, variables)
        loaded <- tryCatch({
            reticulate::import("mofapy2.run.entry_point")
            TRUE
        }, error = function(error) conditionMessage(error))
        state <<- list(applied = isTRUE(loaded), variables = variables,
            reason = if (isTRUE(loaded)) NA_character_ else paste("mofapy2 could not be loaded:", loaded))
        state
    }
})

#' The maths libraries Python has loaded, as threadpoolctl reports them (KIT-151; INT-17): per library its interface,
#' version, threads and code path; an empty list when threadpoolctl cannot be asked.
#' @noRd
integrationMathsLibraries <- function() {
    info <- tryCatch(reticulate::import("threadpoolctl")$threadpool_info(), error = function(error) list())
    lapply(info, function(entry) list(api = as.character(entry$internal_api %||% NA_character_),
        version = as.character(entry$version %||% NA_character_), threads = as.integer(entry$num_threads %||% NA_integer_),
        architecture = as.character(entry$architecture %||% NA_character_)))
}

#' The versions a model was trained with: R, MOFA2, reticulate, Python and its packages; KIT-151 (INT-17): the maths
#' library's pin and the maths libraries loaded.
#' @noRd
integrationMofaVersions <- function() {
    metadata <- tryCatch(reticulate::import("importlib.metadata"), error = function(error) NULL)
    python_package <- function(name) tryCatch(as.character(metadata$version(name)), error = function(error) NA_character_)
    list(R = as.character(getRversion()), MOFA2 = as.character(utils::packageVersion("MOFA2")),
        reticulate = as.character(utils::packageVersion("reticulate")),
        python = tryCatch(as.character(reticulate::py_config()$version_string), error = function(error) NA_character_),
        mofapy2 = python_package("mofapy2"), numpy = python_package("numpy"), scipy = python_package("scipy"),
        h5py = python_package("h5py"), maths_pin = integrationPinMathsLibrary(), maths_libraries = integrationMathsLibraries())
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

#' Train MOFA2 as the proven lane does (KIT-060 status 2026-11-19, run 44423980807): `create_mofa` on the views,
#' the default data, model and training options with the declared ones set, `prepare_mofa`, `run_mofa(use_basilisk =
#' FALSE)` (a fall-back to basilisk is refused: the bundle's Python must be the one used). KIT-148: the model file is
#' written without the training data (`save_data = FALSE`; the aligned values are kept in the project's store instead),
#' the model `run_mofa` loads is the one used (one load, not two), the R copy MOFA2 made for the hand-off is released
#' after training, and MOFA2's overall-level check, which needs the training data, is applied to the views here
#' (`integrationMofaOverallLevelWarning`), its warning kept where MOFA2 would have given it.
#' @param views Named list of features x samples matrices (feature names unique across views).
#' @param options `integrationMofaOptions()`.
#' @return `list(factors, weights, r2_per_factor (factors x views), r2_total, elbo, iterations, options, versions,
#'   messages, model_file)`.
#' @noRd
integrationTrainMofa <- function(views, options) {
    # KIT-151 (INT-17): pinned before MOFA2 starts Python, when nothing has started it yet.
    integrationPinMathsLibrary()
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
    model <- withCallingHandlers(
        MOFA2::run_mofa(object, outfile = outfile, save_data = FALSE, use_basilisk = FALSE),
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
    stats_slot <- methods::slot(model, "training_stats")
    # KIT-151 (INT-14): a further seed's model keeps only its factors and its training length; KIT-240 (VAL-221 (1)): and
    # its final ELBO.
    if (isTRUE(options$factors_only)) {
        return(list(factors = MOFA2::get_factors(model)[[1L]], iterations = length(stats_slot$time),
            elbo_last = integrationLastFinite(stats_slot$elbo)))
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
#' @param maxiter The round limit training ran under (KIT-151: `elbo$maxiter` and `elbo$stopped_by`,
#'   `integrationStoppedBy`); NA when not known.
#' @param seeds KIT-151 (INT-14): `integrationSeedStability()`'s value, kept as `seeds`; NULL when not trained so.
#' @noRd
integrationModelProduct <- function(aligned, fit, store, maxiter = NA_integer_, seeds = NULL) {
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
            trace = elbo, maxiter = as.integer(maxiter),
            stopped_by = if (is.na(maxiter)) NA_character_ else integrationStoppedBy(fit$iterations, maxiter)),
        samples = aligned$samples, view_summary = aligned$view_summary, shared_samples = aligned$shared_samples,
        notes = aligned$notes, options = fit$options, versions = fit$versions, messages = fit$messages, seeds = seeds,
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
            # KIT-153: from the matrices `input_matrix` names (the normalised ones, or the RUV-corrected ones).
            slot <- integrationModelViewsSlot(resolved)
            object@views[resolved$views] <- methods::slot(readWorkflowExecutionComponents(store$manager, store$source,
                fields = paste0(slot, "/", resolved$views)), slot)[resolved$views]
            # KIT-148: the rows carry the hand-off names (`<view>:<feature_id>`) from the read, so the aligned views are
            # handed to the trainer as they are, with no renamed copy.
            integrationModelFromAligned(integrationAlignViews(object, resolved, prefix_rows = TRUE), resolved, parameters, store)
        })
    declaration$partitioned <- function(parameters, store) {
        object <- readWorkflowExecutionComponents(store$manager, store$source, fields = c("features", "design_matrix", "inputs"))
        resolved <- integrationResolveMofaParameters(object, parameters)
        aligned <- integrationAlignViewsPartitioned(object, resolved,
            integrationStateViewPages(store$manager, store$source, store$worker_cap, integrationModelViewsSlot(resolved)),
            prefix_rows = TRUE)
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
    options <- integrationMofaOptions(resolved)
    fit <- trainer(aligned$views, options)
    reclaimWorkflowMemory()
    # KIT-151 (INT-14): the other seeds, one model at a time; of each only its factors and its training length are kept,
    # KIT-240 (VAL-221 (1)): and its final ELBO (NA from a trainer that gives none).
    others <- list()
    for (seed in integrationTrainingSeeds(resolved)[-1L]) {
        other <- trainer(aligned$views, utils::modifyList(options, list(seed = seed, factors_only = TRUE)))
        others[[length(others) + 1L]] <- list(seed = seed, factors = other$factors, iterations = as.integer(other$iterations),
            elbo_last = as.numeric(other$elbo_last %||% NA_real_))
        other <- NULL
        reclaimWorkflowMemory()
    }
    seeds <- integrationSeedStability(fit, others, resolved)
    integrationModelProduct(aligned, fit, store, maxiter = resolved$maxiter, seeds = seeds)
}

#' The seeds a model is trained from (KIT-151; INT-14): the declared seed, then the next `training_seeds - 1` whole
#' numbers (seed 20261119 with 5 seeds: 20261119 to 20261123). The first is the reported model's.
#' @param parameters Resolved parameters.
#' @noRd
integrationTrainingSeeds <- function(parameters) as.integer(parameters$seed) + seq_len(as.integer(parameters$training_seeds)) - 1L

#' Pair the factors of two models one to one (KIT-151; INT-14): every pair's absolute correlation of the factor values
#' over the biological samples both have (pairwise complete; a factor without spread correlates 0); the pair with the
#' largest correlation first, then the next largest among the factors left, and so on (equal values in the order of the
#' reported model's factors, then the other's). A factor left without a partner (the other model kept fewer factors)
#' agrees 0.
#' @param reported,other Samples x factors matrices, named by sample and factor.
#' @return `data.frame(factor, matched_factor, abs_r)`, one row per factor of `reported`, in its order.
#' @noRd
integrationMatchFactors <- function(reported, other) {
    samples <- intersect(rownames(reported), rownames(other))
    r <- suppressWarnings(abs(stats::cor(reported[samples, , drop = FALSE], other[samples, , drop = FALSE],
        use = "pairwise.complete.obs")))
    r[is.na(r)] <- 0
    pairs <- expand.grid(i = seq_len(nrow(r)), j = seq_len(ncol(r)))
    pairs$abs_r <- r[cbind(pairs$i, pairs$j)]
    pairs <- pairs[order(-pairs$abs_r, pairs$i, pairs$j), , drop = FALSE]
    out <- data.frame(factor = colnames(reported), matched_factor = NA_character_, abs_r = 0, stringsAsFactors = FALSE)
    used <- integer()
    for (k in seq_len(nrow(pairs))) {
        i <- pairs$i[[k]]
        j <- pairs$j[[k]]
        if (!is.na(out$matched_factor[[i]]) || j %in% used) next
        out$matched_factor[[i]] <- colnames(other)[[j]]
        out$abs_r[[i]] <- pairs$abs_r[[k]]
        used <- c(used, j)
    }
    out
}

#' Each factor's agreement across the training seeds (KIT-151; INT-14, audit SCI-06, SCI-20): the reported model's
#' factors paired with each other seed's (`integrationMatchFactors`), and per factor the median and the lowest absolute
#' correlation over the other seeds.
#' @param fit The reported model (`integrationTrainMofa()`'s form).
#' @param others Per other seed `list(seed, factors, iterations, elbo_last)`.
#' @param parameters Resolved parameters.
#' @return `list(rule, seeds, models (seed, iterations, factors, stopped_by, elbo_last), agreement (factor, seed,
#'   matched_factor, abs_r), stability (factor, other_models, agreement_median, agreement_lowest), elbo_rank)`; KIT-240
#'   (VAL-221 (1)): `elbo_last` each model's final ELBO, `elbo_rank` the reported model's rank by it
#'   (`integrationElboRank`). The stability table describes only the reported model's factors.
#' @noRd
integrationSeedStability <- function(fit, others, parameters) {
    reported <- fit$factors
    seeds <- integrationTrainingSeeds(parameters)
    models <- data.frame(seed = seeds,
        iterations = as.integer(c(fit$iterations, vapply(others, function(o) o$iterations, integer(1)))),
        factors = as.integer(c(ncol(reported), vapply(others, function(o) ncol(o$factors), integer(1)))),
        stringsAsFactors = FALSE)
    models$stopped_by <- vapply(models$iterations, integrationStoppedBy, character(1), maxiter = parameters$maxiter)
    models$elbo_last <- c(integrationLastFinite(fit$elbo),
        vapply(others, function(o) as.numeric(o$elbo_last %||% NA_real_), numeric(1)))
    agreement <- do.call(rbind, c(list(data.frame(factor = character(), seed = integer(), matched_factor = character(),
        abs_r = numeric(), stringsAsFactors = FALSE)), lapply(others, function(o) {
        m <- integrationMatchFactors(reported, o$factors)
        data.frame(factor = m$factor, seed = rep(as.integer(o$seed), nrow(m)), matched_factor = m$matched_factor,
            abs_r = m$abs_r, stringsAsFactors = FALSE)
    })))
    rownames(agreement) <- NULL
    factors <- colnames(reported)
    stability <- data.frame(factor = factors, other_models = length(others),
        agreement_median = vapply(factors, function(f) if (length(others)) stats::median(agreement$abs_r[agreement$factor == f])
            else NA_real_, numeric(1), USE.NAMES = FALSE),
        agreement_lowest = vapply(factors, function(f) if (length(others)) min(agreement$abs_r[agreement$factor == f])
            else NA_real_, numeric(1), USE.NAMES = FALSE),
        stringsAsFactors = FALSE)
    out <- list(rule = "the declared seed, then the next whole numbers", seeds = seeds, models = models, agreement = agreement,
        stability = stability)
    # Absent (not NULL-valued) when the rank is not known, as in a product trained before KIT-240.
    out$elbo_rank <- integrationElboRank(models$elbo_last)
    out
}

#' The last finite value of an ELBO trace (KIT-240); NA when there is none.
#' @param trace The trace (MOFA2's `training_stats$elbo`).
#' @noRd
integrationLastFinite <- function(trace) {
    trace <- as.numeric(trace)
    trace <- trace[is.finite(trace)]
    if (length(trace)) trace[[length(trace)]] else NA_real_
}

#' The reported model's rank among the seeds by final ELBO (KIT-240; VAL-221 (1)). The ELBO is MOFA2's training
#' objective, a lower bound on the model evidence; higher is better. All seeds train one model on the same data with the
#' same options, so their values are comparable.
#' @param elbo_last The final ELBO of each model, the reported (declared seed's) model first.
#' @return `list(rank, of, below_best, range, relative_spread)`: the rank (1 = highest; ties take the best rank), the
#'   number of models, the highest value less the reported model's and the highest less the lowest (ELBO units), and the
#'   range over the absolute highest value (kept in the product and the record, not shown); NULL with fewer than two
#'   models or a missing value.
#' @noRd
integrationElboRank <- function(elbo_last) {
    if (length(elbo_last) < 2L || anyNA(elbo_last)) return(NULL)
    best <- max(elbo_last)
    list(rank = as.integer(sum(elbo_last > elbo_last[[1L]]) + 1L), of = length(elbo_last), below_best = best - elbo_last[[1L]],
        range = best - min(elbo_last), relative_spread = (best - min(elbo_last)) / abs(best))
}

#' An ELBO difference as the Results tab and the methods paragraph write it (KIT-240): two decimals, in ELBO units; a
#' positive distance below 0.005 reads "less than 0.01", never "0.00" (AUDIT-KIT-240-impl I01).
#' @param x A difference of final ELBOs.
#' @noRd
integrationElboText <- function(x) {
    x <- as.numeric(x)
    if (is.finite(x) && x > 0 && x < 0.005) "less than 0.01" else sprintf("%.2f", x)
}

#' How training stopped (KIT-151; INT-15): `convergence_rule` when it took fewer rounds than the limit (MOFA2's rule
#' stopped it), `round_limit` when it reached the limit.
#' @param iterations The rounds trained.
#' @param maxiter The limit.
#' @noRd
integrationStoppedBy <- function(iterations, maxiter) {
    if (isTRUE(as.integer(iterations) < as.integer(maxiter))) "convergence_rule" else "round_limit"
}

#' The declaration a model run takes (KIT-063 item A): the registry's, or with a sample match its declaration with
#' `sample_match`. In process and on the worker alike. KIT-151 (audit ENG-24): with MOFA2 training, the MOFA2, mofapy2
#' and numpy versions are among the factors a saved model is looked up by, so a model saved before an upgrade of those
#' is not reused.
#' @param parameters The resolved parameters.
#' @param software NULL (a stand-in trainer), or `integrationMofaSoftware()`.
#' @noRd
integrationMofaDeclaration <- function(parameters, software = NULL) {
    declaration <- if (is.null(parameters[["sample_match", exact = TRUE]])) {
        workflowOperationEntry(INTEGRATION_MODEL_OPERATION)$declaration
    } else {
        integrationMofaOperation(match = TRUE)
    }
    # Nothing known (MOFA2 not installed, Python not running): nothing to key on.
    if (!is.null(software) && !all(unlist(software) == "unknown")) declaration$factors <- c(declaration$factors, software)
    declaration
}

#' The software a MOFA2 model depends on beyond the package build (KIT-151; audit ENG-24): the versions of MOFA2,
#' mofapy2 and numpy, as lookup factors (`mofa2`, `mofapy2`, `numpy`; `unknown` when one cannot be read). Python is
#' asked only when it is running (the capability probe has started it before a model is trained); nothing here starts it.
#' @noRd
integrationMofaSoftware <- function() {
    running <- requireNamespace("reticulate", quietly = TRUE) &&
        isTRUE(tryCatch(reticulate::py_available(initialize = FALSE), error = function(error) FALSE))
    metadata <- if (running) tryCatch(reticulate::import("importlib.metadata"), error = function(error) NULL)
    python_package <- function(name) {
        if (is.null(metadata)) return("unknown")
        tryCatch(as.character(metadata$version(name)), error = function(error) "unknown")
    }
    list(mofa2 = tryCatch(as.character(utils::packageVersion("MOFA2")), error = function(error) "unknown"),
        mofapy2 = python_package("mofapy2"), numpy = python_package("numpy"))
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
        declaration <- integrationMofaDeclaration(resolved, software = if (is.null(trainer)) integrationMofaSoftware())
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
    run <- runSavedOperation(manager, integrationMofaDeclaration(resolved, software = integrationMofaSoftware()), source,
        c(resolved, list(variant = variant)),
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
