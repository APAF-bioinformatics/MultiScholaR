# ------------------------------------------------------------------------------------------------
# The integration workbook's functions (KIT-154): one exported function per workbook step, each a thin wrapper over
# the module's internal functions (the same calls the workbook made through `:::` until KIT-154, in the same order and
# with the same arguments). No behaviour of the module changes here; the wrappers check their arguments in plain words
# and return lists of tables, figures and file paths.
# ------------------------------------------------------------------------------------------------

#' A refusal of a workbook function's argument, in plain words (KIT-154).
#' @param message The sentence.
#' @noRd
integrationApiRefuse <- function(message) {
    stop(structure(class = c("integration_api_invalid", "integration_refused", "error", "condition"),
        list(message = message, call = NULL)))
}

#' Whether a value is one non-empty string.
#' @noRd
integrationApiString <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)

#' Check that a value is a list holding the named elements (a value returned by an earlier workbook function).
#' @noRd
integrationApiAssertList <- function(x, elements, argument, source) {
    if (!is.list(x) || !all(elements %in% names(x))) {
        integrationApiRefuse(sprintf("%s must be the value returned by %s", argument, source))
    }
    invisible(x)
}

#' Check settings overrides: a named list whose names are settings of the table.
#' @noRd
integrationApiAssertOverrides <- function(overrides, spec, argument, excluded = character()) {
    if (!is.list(overrides)) integrationApiRefuse(sprintf("%s must be a list of settings and values", argument))
    if (length(overrides) == 0L) return(invisible(overrides))
    if (is.null(names(overrides)) || any(!nzchar(names(overrides))) || anyDuplicated(names(overrides))) {
        integrationApiRefuse(sprintf("every value in %s must be named once by its setting", argument))
    }
    unknown <- setdiff(names(overrides), names(spec))
    if (length(unknown)) {
        integrationApiRefuse(sprintf("these are not settings of the table: %s (see integrationSettings())",
            paste(unknown, collapse = ", ")))
    }
    fixed <- intersect(names(overrides), excluded)
    if (length(fixed)) {
        integrationApiRefuse(sprintf("these settings are fixed by the confirmed sample match and cannot be changed here: %s",
            paste(fixed, collapse = ", ")))
    }
    invisible(overrides)
}

#' Check a project handle (`integrationOpenProject()`'s value).
#' @noRd
integrationApiAssertProject <- function(project) {
    integrationApiAssertList(project, c("manager", "store_root", "paths", "views"), "project", "integrationOpenProject()")
}

#' The integration module's settings tables
#'
#' @description The settings of the model, of the results and of the enrichment, each with its default and the
#'   decision that settled it, as the Integration tab and the workbook show them. A setting is changed by naming it
#'   in the `overrides` of [integrationTrainModel()], [integrationGetResults()] or [integrationEnrichFactor()].
#' @return A list: `model` and `results`, data frames with columns `parameter`, `default`, `choices` and
#'   `settled_by`; `enrichment`, a data frame with columns `setting`, `default` and `settled_by`; and
#'   `input_matrix_help`, the sentence explaining which values each data set enters the model with.
#' @examples
#' \dontrun{
#' settings <- integrationSettings()
#' settings$model
#' }
#' @export
integrationSettings <- function() {
    enrichment_spec <- integrationEnrichmentParameterSpec()
    list(model = integrationMofaParameterTable(), results = integrationResultsParameterTable(),
        enrichment = data.frame(setting = names(enrichment_spec),
            default = vapply(enrichment_spec, function(entry) entry$default_text %||% paste(format(entry$default),
                collapse = ""), character(1)),
            settled_by = vapply(enrichment_spec, `[[`, character(1), "settled_by"), row.names = NULL),
        input_matrix_help = integrationGuiTexts()[["help_input_matrix"]])
}

#' Find the omics' tables saved for integration
#'
#' @description Lays out a project's integration folders and lists the tables each omic's workflow saved at its
#'   Summary step under `<project>/integration/<omic>/`, saying for each omic whether it can be used and, if not, why.
#'   The folder for the workbook's scratch files (`<project>/scripts/integration`) is created when absent.
#' @param project_dir The project's folder.
#' @return A list: `project_dir`, the project's folder (normalised); `paths`, the integration folders (a named list
#'   of paths: `base_dir`, `integration_dir`, `source_dir`, `results_dir`, `multiomic_inputs_dir`, `mofa_tables_dir`,
#'   `mofa_plots_dir`, `study_report_dir`); `omics`, a data frame with one row per omic found and its status; `views`,
#'   a data frame with one row per data set found; and `scan`, the scan itself, passed on to
#'   [integrationOpenProject()].
#' @examples
#' \dontrun{
#' # A project holding the synthetic corpus under integration/<omic>/.
#' inputs <- integrationFindInputs("path/to/project")
#' inputs$omics
#' }
#' @export
integrationFindInputs <- function(project_dir) {
    if (!integrationApiString(project_dir) || !dir.exists(project_dir)) {
        integrationApiRefuse("project_dir must be the path of an existing project folder")
    }
    project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
    results_dir <- file.path(project_dir, "results", "integration")
    paths <- list(
        base_dir = project_dir,
        integration_dir = file.path(project_dir, "integration"),
        source_dir = file.path(project_dir, "scripts", "integration"),
        results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"),
        mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(project_dir, "results_summary", "integration", "Study_report")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    scan <- integrationScanProject(paths$integration_dir)
    list(project_dir = project_dir, paths = paths, omics = integrationScanDisplay(scan), views = scan$views, scan = scan)
}

#' Import the omics' tables into a new integration store
#'
#' @description Opens a new store for the analysis and imports the data sets of every omic that can be used. The
#'   store holds every later step (the model, the results, the enrichments) until [integrationCloseProject()].
#' @param inputs The value of [integrationFindInputs()].
#' @param store_root The folder of the new store, or `NULL` for a new folder under the project's
#'   `scripts/integration`.
#' @return A list: `manager`, the open store (passed on to the later steps); `store_root`, its folder; `paths`, the
#'   project's integration folders; `views`, the names of the data sets imported.
#' @examples
#' \dontrun{
#' inputs <- integrationFindInputs("path/to/project")
#' project <- integrationOpenProject(inputs)
#' project$views
#' integrationCloseProject(project)
#' }
#' @export
integrationOpenProject <- function(inputs, store_root = NULL) {
    integrationApiAssertList(inputs, c("paths", "scan"), "inputs", "integrationFindInputs()")
    if (!is.null(store_root) && !integrationApiString(store_root)) {
        integrationApiRefuse("store_root must be NULL or the path of a folder for the new store")
    }
    request <- integrationImportRequest(inputs$scan, integrationScanReadyViews(inputs$scan)$view)
    if (is.null(store_root)) store_root <- tempfile("integration-store-", tmpdir = inputs$paths$source_dir)
    manager <- newIntegrationStateOwner(store_root)
    importIntegrationInputs(manager, request$dirs)
    list(manager = manager, store_root = store_root, paths = inputs$paths, views = request$views)
}

#' Match the samples across the omics
#'
#' @description Matches every run to the biological sample its omic saved (never by order), with an edited sample
#'   match table when one is given, and shows the result before it is confirmed.
#' @param project The value of [integrationOpenProject()].
#' @param sample_match_file `NULL` (or `""`) to match by the saved biological samples, or the path of an edited
#'   sample match table (tab-separated, as written by [integrationConfirmMatch()]).
#' @return A list: `facts`, sentences describing the match; `data_sets`, a data frame with one row per data set
#'   (runs, biological samples, internal standards dropped); `samples`, a data frame with one row per biological
#'   sample (its group and the data sets it is in); `preview` and `object`, the match and the imported inputs,
#'   passed on to [integrationConfirmMatch()].
#' @examples
#' \dontrun{
#' match <- integrationMatchSamples(project)
#' writeLines(match$facts)
#' match$samples
#' }
#' @export
integrationMatchSamples <- function(project, sample_match_file = NULL) {
    integrationApiAssertProject(project)
    if (identical(sample_match_file, "")) sample_match_file <- NULL
    if (!is.null(sample_match_file) && (!integrationApiString(sample_match_file) || !file.exists(sample_match_file))) {
        integrationApiRefuse("sample_match_file must be NULL, empty, or the path of an existing sample match table")
    }
    object <- integrationImportedObject(project$manager, INTEGRATION_PREVIEW_FIELDS)
    match_table <- if (!is.null(sample_match_file)) integrationReadSampleMatch(sample_match_file)
    preview <- integrationAlignmentPreview(object, project$views, match_table, manager = project$manager)
    list(facts = integrationAlignmentFacts(preview), data_sets = integrationAlignmentViewTable(preview),
        samples = integrationAlignmentSampleTable(preview), preview = preview, object = object)
}

#' Confirm the sample match
#'
#' @description Confirms a sample match: writes the sample match template (`integration_sample_match_template.tsv`,
#'   to edit when the omics name one biological sample differently), the sample design and, when a match table was
#'   given, the match, all to the project's `results/integration/multiomic_inputs` folder.
#' @param project The value of [integrationOpenProject()].
#' @param match The value of [integrationMatchSamples()].
#' @return The confirmation, a list: `samples`, the matched samples (a data frame); `view_summary`, a data frame with
#'   one row per data set; `parameters`, the matching settings; `manifests`, each input's manifest checksum; `file`,
#'   the path of the saved design; `confirmed_at`, the time; `sample_match_file`, the path of the saved match (only
#'   when a match table was given); and `template_file`, the path of the sample match template. Passed on to
#'   [integrationTrainModel()].
#' @examples
#' \dontrun{
#' confirmed <- integrationConfirmMatch(project, integrationMatchSamples(project))
#' confirmed$template_file
#' }
#' @export
integrationConfirmMatch <- function(project, match) {
    integrationApiAssertProject(project)
    integrationApiAssertList(match, c("preview", "object"), "match", "integrationMatchSamples()")
    preview <- match$preview
    dir <- project$paths$multiomic_inputs_dir
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    template <- integrationSampleMatchTemplate(match$object, preview$parameters$sample_match_table)
    template_file <- file.path(dir, "integration_sample_match_template.tsv")
    utils::write.table(template, template_file, sep = "\t", quote = FALSE, row.names = FALSE)
    confirmed <- list(
        samples = preview$aligned$samples,
        view_summary = preview$aligned$view_summary,
        parameters = preview$parameters,
        manifests = integrationStateManifests(match$object),
        file = integrationSaveDesign(preview$aligned$samples, dir),
        confirmed_at = Sys.time()
    )
    if (!is.null(preview$parameters$sample_match_table)) {
        confirmed$sample_match_file <- integrationSaveSampleMatch(preview$parameters$sample_match_table, dir)
    }
    confirmed$template_file <- template_file
    confirmed
}

#' Train the integration model
#'
#' @description Trains the MOFA2 model on the confirmed match with the model settings' defaults (see
#'   [integrationSettings()]), changed by `input_matrix` and `overrides`. A model already trained on the same inputs
#'   with the same settings and software is reused. How samples are matched (technical replicates, internal
#'   standards, missing samples) is fixed by the confirmed match.
#' @param project The value of [integrationOpenProject()].
#' @param confirmed The value of [integrationConfirmMatch()].
#' @param input_matrix `NULL` for the default, or which values each data set enters the model with (one of the
#'   `input_matrix` choices in `integrationSettings()$model`).
#' @param overrides A named list of model settings and their values, for example `list(num_factors = 8L)`.
#' @return A list: `status`, `"computed"` or `"hit"` (a saved model reused); `value`, the trained model;
#'   `reference`, its reference in the store; `parameters`, the settings it was trained with; `summary`, sentences
#'   on training (convergence, factors, samples, the saved model file); and `variance`, a data frame of the variance
#'   explained per data set. Passed on to [integrationGetResults()].
#' @examples
#' \dontrun{
#' run <- integrationTrainModel(project, confirmed, overrides = list(num_factors = 5L))
#' writeLines(run$summary)
#' run$variance
#' }
#' @export
integrationTrainModel <- function(project, confirmed, input_matrix = NULL, overrides = list()) {
    integrationApiAssertProject(project)
    integrationApiAssertList(confirmed, c("parameters", "samples"), "confirmed", "integrationConfirmMatch()")
    model_spec <- integrationMofaParameterSpec()
    if (is.null(input_matrix)) input_matrix <- model_spec$input_matrix$default
    if (!integrationApiString(input_matrix) || !input_matrix %in% model_spec$input_matrix$choices) {
        integrationApiRefuse(sprintf("input_matrix must be one of %s", paste(model_spec$input_matrix$choices, collapse = ", ")))
    }
    integrationApiAssertOverrides(overrides, model_spec, "overrides", INTEGRATION_GUI_MATCHING_PARAMETERS)
    settings <- lapply(model_spec, `[[`, "default")
    settings$input_matrix <- input_matrix
    settings[names(overrides)] <- overrides
    settings <- settings[setdiff(names(settings), INTEGRATION_GUI_MATCHING_PARAMETERS)]
    # Without a feature selection there is no number of features to keep and no observed share.
    if (identical(settings$feature_selection, "none")) {
        settings$features_per_view <- NA_integer_
        settings$min_observed_fraction <- NA_real_
    }
    model_parameters <- integrationModelParameters(confirmed$parameters, settings)
    run <- runIntegrationModel(project$manager, model_parameters, confirmation = integrationConfirmationRecord(confirmed))
    run_parameters <- project$manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]$parameters
    model_summary <- integrationModelSummary(run$value, run_parameters)
    list(status = run$status, value = run$value, reference = run$reference, parameters = run_parameters,
        summary = c(model_summary$convergence, model_summary$factors, model_summary$samples, model_summary$saved),
        variance = model_summary$variance)
}

#' Get the results of a trained model
#'
#' @description Reads the standard results of a trained model with the results settings' defaults (see
#'   [integrationSettings()]), changed by `overrides`: the variance explained, each factor's stability across the
#'   training seeds, the factors against the groups, the weights, the factors' correlation, the top features'
#'   values and, when chosen, UMAP and t-SNE maps of the factors.
#' @param project The value of [integrationOpenProject()].
#' @param model The value of [integrationTrainModel()].
#' @param overrides A named list of results settings and their values, for example `list(umap = TRUE)`.
#' @return A list: `messages`, MOFA2's warnings in words and the model's notes; `claim_limit`, the statement to read
#'   before the factors' association with the groups; `maps_text`, a map's refusal or the caution beside the maps;
#'   `factor`, the first factor shown; `tables`, a list of data frames (`variance`, `stability`, `association`,
#'   `correlation`); `figures`, a list of ggplot objects (`variance`, `variance_total`, `factors_by_group`,
#'   `factor_correlation`, `factor_pairs`; `weights`, `top_feature_scatter` and `top_feature_heatmap`, each a list by
#'   data set for `factor`; `maps`, a list by method of the maps drawn); and `results`, the results themselves,
#'   passed on to [integrationExportResults()] and [integrationEnrichFactor()].
#' @examples
#' \dontrun{
#' results <- integrationGetResults(project, run)
#' print(results$figures$variance)
#' results$tables$association
#' }
#' @export
integrationGetResults <- function(project, model, overrides = list()) {
    integrationApiAssertProject(project)
    integrationApiAssertList(model, c("value", "reference", "parameters"), "model", "integrationTrainModel()")
    results_spec <- integrationResultsParameterSpec()
    integrationApiAssertOverrides(overrides, results_spec, "overrides")
    results_parameters <- lapply(results_spec, `[[`, "default")
    results_parameters[names(overrides)] <- overrides
    results <- integrationResults(model$value, integrationImportedObject(project$manager), results_parameters,
        manager = project$manager, model_parameters = model$parameters)
    factor <- results$shown[[1L]]
    views <- stats::setNames(results$views, results$views)
    drawn <- Filter(function(map) is.null(map$refusal), results$maps)
    refused <- Filter(function(map) !is.null(map$refusal), results$maps)
    maps_text <- unlist(lapply(refused, `[[`, "refusal"), use.names = FALSE)
    if (length(drawn)) {
        maps_text <- c(maps_text, integrationMapCaution(results, drawn[[1L]]$settings$samples, names(drawn)))
    }
    list(
        messages = c(vapply(results$advisories, function(advisory) advisory$sentence, character(1)), results$notes),
        claim_limit = integrationClaimLimitText(results),
        maps_text = as.character(maps_text),
        factor = factor,
        tables = list(variance = integrationResultsVarianceTable(results),
            stability = integrationResultsStabilityTable(results),
            association = integrationResultsAssociationTable(results),
            correlation = integrationResultsCorrelationTable(results)),
        figures = list(variance = integrationPlotVariance(results),
            variance_total = integrationPlotVarianceTotal(results),
            factors_by_group = integrationPlotFactorsByGroup(results),
            weights = lapply(views, function(view) integrationPlotWeights(results, factor, view)),
            factor_correlation = integrationPlotFactorCorrelation(results),
            factor_pairs = integrationPlotFactorPairs(results),
            top_feature_scatter = lapply(views, function(view) integrationPlotTopFeatureScatter(results, factor, view)),
            top_feature_heatmap = lapply(views, function(view) integrationPlotTopFeatureHeatmap(results, factor, view)),
            maps = lapply(stats::setNames(names(drawn), names(drawn)), function(method) integrationPlotFactorMap(results,
                method))),
        results = results)
}

#' Save the results of a trained model
#'
#' @description Writes every results table and figure to the project's `results/integration/mofa/tables` and
#'   `results/integration/mofa/plots` folders, and keeps them with the model for the analysis record. A map chosen
#'   but refused for its settings stops the save with its sentence.
#' @param project The value of [integrationOpenProject()].
#' @param model The value of [integrationTrainModel()].
#' @param results The value of [integrationGetResults()] for that model.
#' @return A list: `tables` and `figures`, the paths of the files written; `tables_dir` and `plots_dir`, the folders.
#' @examples
#' \dontrun{
#' saved <- integrationExportResults(project, run, results)
#' basename(c(saved$tables, saved$figures))
#' }
#' @export
integrationExportResults <- function(project, model, results) {
    integrationApiAssertProject(project)
    integrationApiAssertList(model, c("value", "reference", "parameters"), "model", "integrationTrainModel()")
    integrationApiAssertList(results, "results", "results", "integrationGetResults()")
    saved <- integrationSaveResults(project$manager, model$reference, results$results, project$paths)
    list(tables = file.path(saved$tables_dir, saved$tables), figures = file.path(saved$plots_dir, saved$figures),
        tables_dir = saved$tables_dir, plots_dir = saved$plots_dir)
}

#' Enrich a factor of the trained model
#'
#' @description Enriches each data set of the model for one factor, ranking every measured feature by its
#'   correlation with the factor: protein and gene data sets against STRING, metabolite and lipid data sets against
#'   KEGG and Reactome. STRING is skipped without a key; KEGG and Reactome are skipped for a data set whose
#'   identifier column is not named, and their pathway tables are downloaded with their release recorded. The key is
#'   used for the query only and is never written to the project, the record or the value returned.
#' @param project The value of [integrationOpenProject()].
#' @param results The value of [integrationGetResults()].
#' @param factor The factor to enrich, or `NULL` for the first factor shown (`results$factor`).
#' @param api_key Your STRING key, or `NULL` (or `""`) to skip STRING. Read it from the environment (for example
#'   `Sys.getenv("STRING_API_KEY")`); never write it in a script.
#' @param id_columns A named list giving, for a metabolite or lipid data set, the column of its feature table that
#'   holds KEGG compound or ChEBI identifiers, for example `list(metabolomics_LCMS_Pos = "kegg_id")`.
#' @param overrides A named list of enrichment settings and their values (see `integrationSettings()$enrichment`).
#'   `test` chooses the test for KEGG and Reactome: `"camera"`, the test of record and the default (limma's camera
#'   on each feature's regression on the factor, allowing for the correlation measured between a pathway's features,
#'   on the features measured in every biological sample); `"camera_floor"`, the same with negative estimated
#'   correlations counted as zero (it never calls a pathway `"camera"` does not call); `"permutation"`, PCGSE's
#'   permutation test of the factor's values (each feature's correlation with the factor against `permutations`
#'   random reorderings of the factor's values across the biological samples, 1000 or 10000, with `seed`; a message
#'   gives the predicted time first); or `"gsea"`, the pre-ranked test used before (each feature's correlation with
#'   the factor, tested as if features were independent). STRING ignores it. `adjust_ruv_factors` (`TRUE` or `FALSE`),
#'   "Adjust for removed RUV factors", applies to `"camera"` and `"camera_floor"`: with `TRUE` the test puts the unwanted
#'   factors that the omic's RUV correction removed from the data set into its model, so a feature counts for what it
#'   shares with the factor beyond them. Left out, it is `TRUE` for a model trained on RUV-corrected values
#'   (`input_matrix = "ruv_corrected"`) and `FALSE` otherwise. `FALSE` on a model trained on RUV-corrected values is a
#'   sensitivity analysis (the values before RUV, without those factors); `TRUE` needs a data set to which the omic
#'   applied RUV, saved by this version or later. Every result says which form it is. The pre-ranked test, the
#'   permutation test and STRING do not use the setting.
#' @return A list: `enrichments`, a list named `"<data set> <database>"` of the enrichments made, each a list with
#'   `summary` (sentences) and `terms` (a data frame of the terms found); and `skipped`, sentences naming each data
#'   set and database skipped or refused, and why.
#' @examples
#' \dontrun{
#' enriched <- integrationEnrichFactor(project, results, api_key = Sys.getenv("STRING_API_KEY"))
#' writeLines(enriched$skipped)
#' # The pre-ranked test used before KIT-183, for KEGG and Reactome:
#' earlier <- integrationEnrichFactor(project, results, id_columns = list(metabolomics_GCMS = "kegg_id"),
#'     overrides = list(test = "gsea"))
#' # The permutation test of the factor's values (PCGSE), with its seed:
#' permuted <- integrationEnrichFactor(project, results, id_columns = list(metabolomics_GCMS = "kegg_id"),
#'     overrides = list(test = "permutation", seed = 1L))
#' # For a model trained on RUV-corrected values, the sensitivity form without the removed RUV factors:
#' unadjusted <- integrationEnrichFactor(project, results, id_columns = list(metabolomics_GCMS = "kegg_id"),
#'     overrides = list(adjust_ruv_factors = FALSE))
#' names(enriched$enrichments)
#' }
#' @export
integrationEnrichFactor <- function(project, results, factor = NULL, api_key = NULL, id_columns = list(),
    overrides = list()) {
    integrationApiAssertProject(project)
    integrationApiAssertList(results, c("results", "factor"), "results", "integrationGetResults()")
    if (is.null(factor)) factor <- results$factor
    if (!integrationApiString(factor)) integrationApiRefuse("factor must be NULL or the name of one factor")
    if (is.null(api_key)) api_key <- ""
    if (!is.character(api_key) || length(api_key) != 1L || is.na(api_key)) {
        integrationApiRefuse("api_key must be NULL or one string")
    }
    if (!is.list(id_columns) || (length(id_columns) && (is.null(names(id_columns)) ||
        !all(vapply(id_columns, integrationApiString, logical(1)))))) {
        integrationApiRefuse("id_columns must be a list naming, per data set, one identifier column")
    }
    enrichment_spec <- integrationEnrichmentParameterSpec()
    integrationApiAssertOverrides(overrides, enrichment_spec, "overrides")
    enrichment_settings <- lapply(enrichment_spec, `[[`, "default")
    enrichment_settings[names(overrides)] <- overrides
    # KIT-183: the test is sent only when named in `overrides` (absent, the test of record for the data set), and
    # never for STRING.
    if (is.null(overrides[["test"]])) enrichment_settings$test <- NULL
    # KIT-184 (VAL-62): the adjustment for removed RUV factors is sent only when named in `overrides` (absent, it
    # follows the model's input), only with a test that has a design (an absent test is the test of record, which has
    # one) and never for STRING.
    if (is.null(overrides[["adjust_ruv_factors"]]) || !(is.null(overrides[["test"]]) ||
        isTRUE(overrides[["test"]] %in% INTEGRATION_ENRICHMENT_DESIGN_TESTS))) {
        enrichment_settings$adjust_ruv_factors <- NULL
    }
    manager <- project$manager
    view_omics <- integrationViewOmics(integrationImportedObject(manager, "inputs"))
    enrichments <- list()
    skipped <- character()
    for (view in results$results$views) {
        omic <- view_omics[[view]]
        databases <- INTEGRATION_ENRICHMENT_DATABASES[[omic]]
        for (database in databases) {
            if (database == "string" && !nzchar(api_key)) {
                skipped <- c(skipped, sprintf("%s: STRING skipped (no key).", view))
                next
            }
            if (database != "string" && is.null(id_columns[[view]])) {
                skipped <- c(skipped, sprintf("%s: %s skipped (no identifier column named).", view, toupper(database)))
                next
            }
            request <- c(list(factor = factor, view = view, database = database), enrichment_settings)
            if (database == "string") {
                request$test <- NULL
                request$adjust_ruv_factors <- NULL
            }
            mapping <- NULL
            if (database != "string") {
                request$id_column <- id_columns[[view]]
                mapping <- if (database == "kegg") {
                    integrationDownloadKeggMapping(enrichment_settings$kegg_organism)
                } else {
                    integrationDownloadReactomeMapping()
                }
            }
            enrichment <- tryCatch(runIntegrationEnrichment(manager, request,
                api_key = if (database == "string") api_key, mapping = mapping),
                integration_refused = function(condition) {
                    skipped <<- c(skipped, sprintf("%s, %s: %s", view, database, conditionMessage(condition)))
                    NULL
                })
            if (is.null(enrichment)) next
            enrichments[[paste(view, database)]] <- list(summary = integrationEnrichmentSummary(enrichment$run$value,
                entry = enrichment$entry),
                terms = as.data.frame(enrichment$run$value$terms))
        }
    }
    list(enrichments = enrichments, skipped = skipped)
}

#' Write the analysis record and report
#'
#' @description Writes the analysis record (`analysis_record.json` and `analysis_record.xlsx`), the study
#'   parameters, the report, the model file and the sample match to the project's
#'   `results_summary/integration/Study_report` folder, and gives the methods paragraph written from the record.
#' @param project The value of [integrationOpenProject()], with a trained model.
#' @param format The report's format, `"docx"` or `"html"`.
#' @return A list: `directory`, the folder; `files`, the paths of the files written; `methods`, the methods
#'   paragraph.
#' @examples
#' \dontrun{
#' bundle <- integrationWriteRecord(project, format = "html")
#' basename(bundle$files)
#' writeLines(bundle$methods)
#' }
#' @export
integrationWriteRecord <- function(project, format = "docx") {
    integrationApiAssertProject(project)
    if (!integrationApiString(format) || !format %in% c("docx", "html")) {
        integrationApiRefuse("format must be \"docx\" or \"html\"")
    }
    integrationWriteRecordBundle(project$manager, project$paths$study_report_dir, format = format)
}

#' Close an integration store
#'
#' @description Closes the store opened by [integrationOpenProject()]. Its folder is kept, so the analysis can be
#'   repeated from its record with [integrationReplayProject()].
#' @param project The value of [integrationOpenProject()].
#' @return The store's folder (`project$store_root`), invisibly.
#' @examples
#' \dontrun{
#' integrationCloseProject(project)
#' }
#' @export
integrationCloseProject <- function(project) {
    integrationApiAssertProject(project)
    closeWorkflowStateOwner(project$manager)
    invisible(project$store_root)
}

#' Repeat an integration analysis from its record
#'
#' @description Reopens a closed integration store, repeats its analysis from its record alone and compares the
#'   repeated model, results and enrichments with the saved ones. The enrichments are repeated from their recorded
#'   answers and tables, without querying the services again.
#' @param store_root The store's folder (`project$store_root`, or the value of [integrationCloseProject()]).
#' @param scratch A folder for the repeat's own store, or `NULL` for a new folder beside `store_root`.
#' @return A list: `comparison`, the repeated model against the saved one (`equal`, `identical`,
#'   `largest_relative` and a data frame `tables`, one row per table compared); `outputs_comparison`, the repeated
#'   results against the saved ones (`equal`, `identical`, `tables`; `NULL` when none were saved);
#'   `enrichment_comparison`, a data frame with one row per enrichment repeated (`NULL` when there was none).
#' @examples
#' \dontrun{
#' store_root <- integrationCloseProject(project)
#' replay <- integrationReplayProject(store_root)
#' replay$comparison$tables
#' }
#' @export
integrationReplayProject <- function(store_root, scratch = NULL) {
    if (!integrationApiString(store_root) || !dir.exists(store_root)) {
        integrationApiRefuse("store_root must be the folder of an existing integration store")
    }
    if (is.null(scratch)) scratch <- tempfile("integration-replay-", tmpdir = dirname(store_root))
    if (!integrationApiString(scratch)) integrationApiRefuse("scratch must be NULL or the path of a folder")
    replay <- replayIntegrationProject(store_root, scratch)
    list(comparison = replay$comparison, outputs_comparison = replay$outputs_comparison,
        enrichment_comparison = replay$enrichment_comparison)
}
