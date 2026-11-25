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

# KIT-141: the integration workflow's Enrichment tab, behind the integration gate like every integration tab. For the
# model the Model tab trained last, it ranks every measured feature of a chosen data set by its correlation with a
# chosen factor and enriches the ranked list: STRING for proteins and genes (the user's key, used for this query only
# and never kept), KEGG or Reactome for metabolites and lipids (tables downloaded in the app, their release and
# SHA-256 kept). Everything runs in the app process (`runIntegrationEnrichment`); the key never reaches the managed
# worker. The step is optional: the tab opens the Record tab as soon as the Results tab is complete. Recorded divergence
# from conventions/R-shiny.md (module function names, fluidRow layout, no ExtendedTask): the integration tabs' existing
# pattern is kept, under the existing-code rule of conventions/R.md.

#' A setting's label on the Enrichment tab (KIT-240, R1-04): the label, and, for a default not yet confirmed
#' (`integrationParameterProvisional`: its `settled_by` says it awaits confirmation), the mark in brackets after it.
#' @param texts `integrationGuiTexts()`.
#' @param key The label's key in the wording table.
#' @param entry The setting's entry of `integrationEnrichmentParameterSpec()`.
#' @return One string.
#' @noRd
integrationEnrichmentLabel <- function(texts, key, entry) {
    if (integrationParameterProvisional(entry)) sprintf("%s (%s)", texts[[key]], texts[["mark_provisional"]]) else texts[[key]]
}

#' The Enrichment tab's UI.
#' @param id Module id.
#' @noRd
mod_integration_enrichment_ui <- function(id) {
    ns <- shiny::NS(id)
    texts <- integrationGuiTexts()
    spec <- integrationEnrichmentParameterSpec()
    # KIT-240 (R1-04): a default not yet confirmed carries the mark in its label.
    label <- function(key, name) integrationEnrichmentLabel(texts, key, spec[[name]])
    shiny::fluidRow(shiny::column(12,
        shiny::h3(texts[["enrich_heading"]]),
        shiny::p(texts[["enrich_intro"]]),
        shiny::p(shiny::em(texts[["enrich_optional"]])),
        testid(shiny::uiOutput(ns("status")), "integration-enrichment-status"),
        shiny::fluidRow(
            shiny::column(4, shiny::selectInput(ns("factor"), texts[["enrich_factor"]], choices = character())),
            shiny::column(4, shiny::selectInput(ns("view"), texts[["enrich_view"]], choices = character())),
            shiny::column(4, shiny::selectInput(ns("database"), texts[["enrich_database"]], choices = character()))),
        testid(shiny::uiOutput(ns("id_column_ui")), "integration-enrichment-id-column"),
        shiny::h4(texts[["enrich_settings_heading"]]),
        shiny::fluidRow(
            shiny::column(4, shiny::selectInput(ns("correlation_method"), label("enrich_correlation_method", "correlation_method"),
                choices = stats::setNames(c("pearson", "spearman"), unname(texts[c("enrich_correlation_pearson",
                    "enrich_correlation_spearman")])), selected = spec$correlation_method$default)),
            shiny::column(4, shiny::numericInput(ns("min_samples"), label("enrich_min_samples", "min_samples"), spec$min_samples$default,
                min = 3L, step = 1L)),
            # KIT-191: the seed is shown for the tests that use one (the pre-ranked and permutation tests).
            shiny::column(4, testid(shiny::conditionalPanel(sprintf(paste("(input['%s'] == 'kegg' || input['%s'] ==",
                "'reactome') && (input['%s'] == 'gsea' || input['%s'] == 'permutation')"), ns("database"), ns("database"),
                ns("test"), ns("test")),
                shiny::numericInput(ns("seed"), label("enrich_seed", "seed"), spec$seed$default, step = 1L)),
                "integration-enrichment-seed"))),
        testid(shiny::div(
            shiny::fluidRow(
                shiny::column(4, shiny::numericInput(ns("string_species"), label("enrich_string_species", "string_species"),
                    spec$string_species$default, step = 1L)),
                shiny::column(4, shiny::numericInput(ns("string_fdr"), texts[["enrich_string_fdr"]], spec$string_fdr$default,
                    min = 0, max = 1, step = 0.01)),
                shiny::column(4, shiny::passwordInput(ns("string_key"), texts[["enrich_key"]], value = ""))),
            shiny::p(shiny::em(texts[["enrich_key_note"]]))), "integration-enrichment-string"),
        testid(shiny::div(
            shiny::fluidRow(
                shiny::column(3, shiny::textInput(ns("kegg_organism"), label("enrich_kegg_organism", "kegg_organism"),
                    spec$kegg_organism$default)),
                shiny::column(3, shiny::textInput(ns("reactome_organism"), label("enrich_reactome_organism", "reactome_organism"),
                    spec$reactome_organism$default)),
                shiny::column(3, shiny::numericInput(ns("min_set_size"), label("enrich_min_set", "min_set_size"), spec$min_set_size$default,
                    min = 5L, step = 1L)),
                shiny::column(3, shiny::numericInput(ns("max_set_size"), label("enrich_max_set", "max_set_size"), spec$max_set_size$default,
                    min = 5L, step = 1L))),
            # KIT-183: the test, pre-selected with the test of record for the chosen data set, and its sentence.
            testid(shiny::uiOutput(ns("test_ui")), "integration-enrichment-test"),
            # KIT-191 (D1): the number of reorderings, for the permutation test only.
            testid(shiny::conditionalPanel(sprintf("input['%s'] == 'permutation'", ns("test")),
                shiny::selectInput(ns("permutations"), texts[["enrich_permutations"]],
                    choices = stats::setNames(as.character(INTEGRATION_ENRICHMENT_PERMUTATIONS),
                        unname(texts[paste0("enrich_permutations_", INTEGRATION_ENRICHMENT_PERMUTATIONS)])),
                    selected = as.character(spec$permutations$default))), "integration-enrichment-permutations"),
            # KIT-184 (VAL-62): "Adjust for removed RUV factors", for the tests that have a design.
            testid(shiny::uiOutput(ns("adjust_ui")), "integration-enrichment-adjust"),
            testid(shiny::uiOutput(ns("test_note")), "integration-enrichment-test-note"),
            # KIT-191 (D6 (a)): how each test behaved, for the data set's kind.
            testid(shiny::uiOutput(ns("evidence")), "integration-enrichment-evidence"),
            shiny::p(shiny::em(texts[["enrich_kegg_terms"]])),
            testid(shiny::actionButton(ns("download"), texts[["enrich_download"]], icon = shiny::icon("download")),
                "integration-enrichment-download"),
            testid(shiny::uiOutput(ns("download_status")), "integration-enrichment-download-status")),
            "integration-enrichment-pathways"),
        shiny::br(),
        testid(shiny::actionButton(ns("run"), texts[["enrich_run"]], class = "btn-primary", icon = shiny::icon("play")),
            "integration-enrichment-run"),
        testid(shiny::uiOutput(ns("result")), "integration-enrichment-result"),
        shiny::h4(texts[["enrich_terms_heading"]]),
        shiny::div(`data-testid` = "integration-enrichment-terms", DT::DTOutput(ns("terms")))))
}

#' The Enrichment tab's settings as the request `runIntegrationEnrichment` takes (KIT-141): the chosen factor, data set
#' and database, the identifier column for metabolites and lipids, and every setting. The key is not among them.
#' @param input The tab's inputs (a list or Shiny's input).
#' @return Named list.
#' @noRd
integrationEnrichmentRequest <- function(input) {
    out <- list(factor = input$factor, view = input$view, database = input$database,
        correlation_method = input$correlation_method, min_samples = input$min_samples, seed = input$seed,
        string_species = input$string_species, string_fdr = input$string_fdr, kegg_organism = input$kegg_organism,
        reactome_organism = input$reactome_organism, min_set_size = input$min_set_size, max_set_size = input$max_set_size)
    if (isTRUE(input$database %in% c("kegg", "reactome"))) {
        out$id_column <- input$id_column
        # KIT-183: the test for KEGG and Reactome only; STRING ignores it.
        out$test <- input$test
        # KIT-191: the number of reorderings, for the permutation test only.
        if (identical(input$test, "permutation") && !is.null(input$permutations)) {
            out$permutations <- as.integer(input$permutations)
        }
        # KIT-184 (VAL-62): the adjustment for removed RUV factors, for a test that has a design only (the control is
        # hidden for the others, so a hidden control never causes a refusal).
        if (isTRUE(input$test %in% INTEGRATION_ENRICHMENT_DESIGN_TESTS) && !is.null(input$adjust_ruv_factors)) {
            out$adjust_ruv_factors <- isTRUE(input$adjust_ruv_factors)
        }
    }
    Filter(Negate(is.null), out)
}

#' The Enrichment tab's server (KIT-141).
#' @param id Module id.
#' @param workflow_data The integration bus (`state_manager`, `integration_model`, `integration_enrichment`, `tab_status`).
#' @param paths `project_dirs$integration` (unused; the results are kept in the integration store).
#' @param texts `integrationGuiTexts()`.
#' @param query The STRING client (`integrationStringRankedEnrichment`; tests pass a stand-in).
#' @param download The downloader of the pathway tables (`integrationEnrichmentDownload`; tests pass a stand-in).
#' @param withProgressFn The progress display around the download and the run (`shiny::withProgress`).
#' @return The latest run's display, as a reactive.
#' @noRd
mod_integration_enrichment_server <- function(id, workflow_data, paths, texts = integrationGuiTexts(),
    query = integrationStringRankedEnrichment, download = integrationEnrichmentDownload, withProgressFn = shiny::withProgress) {
    shiny::moduleServer(id, function(input, output, session) {
        state <- shiny::reactiveValues(model = NULL, run = NULL, error = NULL, mapping = NULL, download_error = NULL)
        # The Record tab does not wait for this optional step.
        shiny::observeEvent(workflow_data$tab_status$integration_results, {
            status <- workflow_data$tab_status
            if (identical(status$integration_results, "complete") && identical(status$integration_record, "disabled")) {
                status$integration_record <- "pending"
                workflow_data$tab_status <- status
            }
        }, ignoreNULL = TRUE)
        # A new model (or none): its factors and data sets; no earlier result is shown.
        shiny::observeEvent(workflow_data$integration_model, {
            state$run <- NULL
            state$error <- NULL
            workflow_data$integration_enrichment <- NULL
            model <- workflow_data$integration_model
            manager <- workflow_data$state_manager
            state$model <- NULL
            if (is.null(model) || is.null(manager)) return()
            facts <- tryCatch(withWorkflowStateScope(manager, function() {
                value <- hydrateSavedWorkflowEntity(manager$.context$store, model$reference)$value
                object <- integrationImportedObject(manager, fields = c("features", "design_matrix", "inputs"))
                parameters <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]$parameters
                list(factors = setdiff(names(value$factors), "biological_sample"), views = names(value$weights),
                    omics = integrationViewOmics(object)[names(value$weights)],
                    columns = lapply(object@features[names(value$weights)], function(f) setdiff(names(f), "feature_id")),
                    # KIT-183: each data set's biological samples (the test of record and the few-samples sentence)
                    # and the values the model was trained on (the RUV note).
                    samples = vapply(names(value$weights), function(view) length(integrationEnrichmentViewSamples(object,
                        value, parameters, view)), integer(1)),
                    input_matrix = as.character(parameters$input_matrix %||% "uncorrected"))
            }), error = function(error) NULL)
            state$model <- facts
            if (is.null(facts)) return()
            shiny::updateSelectInput(session, "factor", choices = facts$factors, selected = facts$factors[[1L]])
            shiny::updateSelectInput(session, "view", choices = facts$views, selected = facts$views[[1L]])
        }, ignoreNULL = FALSE)
        # The databases of the chosen data set's omic.
        shiny::observe({
            facts <- state$model
            view <- input$view
            if (is.null(facts) || is.null(view) || !view %in% facts$views) return()
            databases <- INTEGRATION_ENRICHMENT_DATABASES[[facts$omics[[view]]]]
            shiny::updateSelectInput(session, "database", choices = stats::setNames(databases,
                unname(texts[paste0("enrich_database_", databases)])), selected = databases[[1L]])
        })
        output$id_column_ui <- shiny::renderUI({
            facts <- state$model
            if (is.null(facts) || !isTRUE(input$database %in% c("kegg", "reactome")) || !isTRUE(input$view %in% facts$views)) {
                return(NULL)
            }
            shiny::tagList(shiny::selectInput(session$ns("id_column"), texts[["enrich_id_column"]],
                choices = facts$columns[[input$view]]), shiny::p(shiny::em(texts[[paste0("enrich_id_column_", input$database)]])))
        })
        # KIT-183: the test selector, for KEGG and Reactome only (STRING has its own test).
        output$test_ui <- shiny::renderUI({
            facts <- state$model
            if (is.null(facts) || !isTRUE(input$database %in% c("kegg", "reactome")) || !isTRUE(input$view %in% facts$views)) {
                return(NULL)
            }
            integrationEnrichmentTestInput(session$ns("test"), facts$samples[[input$view]], texts)
        })
        # KIT-184 (VAL-62, VAL-160): the adjustment's control, shown for a test that has a design; off unless the user
        # switches it on, whatever the model was trained on.
        output$adjust_ui <- shiny::renderUI({
            facts <- state$model
            if (is.null(facts) || !isTRUE(input$database %in% c("kegg", "reactome")) || !isTRUE(input$view %in% facts$views) ||
                !isTRUE(input$test %in% INTEGRATION_ENRICHMENT_DESIGN_TESTS)) {
                return(NULL)
            }
            shiny::tagList(shiny::checkboxInput(session$ns("adjust_ruv_factors"), texts[["enrich_adjust_ruv"]], value = FALSE),
                shiny::p(shiny::em(texts[["enrich_adjust_ruv_help"]])))
        })
        output$test_note <- shiny::renderUI({
            facts <- state$model
            if (is.null(facts) || !isTRUE(input$database %in% c("kegg", "reactome")) || !isTRUE(input$view %in% facts$views) ||
                is.null(input$test)) {
                return(NULL)
            }
            lapply(integrationEnrichmentTestNotes(input$test, facts$samples[[input$view]], facts$input_matrix,
                facts$omics[[input$view]], texts, permutations = as.integer(input$permutations %||% "1000")),
                function(sentence) shiny::p(shiny::em(sentence)))
        })
        output$evidence <- shiny::renderUI({
            facts <- state$model
            if (is.null(facts) || !isTRUE(input$database %in% c("kegg", "reactome")) || !isTRUE(input$view %in% facts$views)) {
                return(NULL)
            }
            integrationEnrichmentEvidenceUi(facts$omics[[input$view]], texts)
        })
        shiny::observeEvent(input$download, {
            state$download_error <- NULL
            database <- input$database
            if (!isTRUE(database %in% c("kegg", "reactome"))) return()
            mapping <- tryCatch(withProgressFn(message = texts[["enrich_downloading"]], value = 0.5,
                if (database == "kegg") integrationDownloadKeggMapping(input$kegg_organism, download) else
                    integrationDownloadReactomeMapping(download)),
                error = function(error) {
                    state$download_error <- integrationRefusalText(error, texts)
                    NULL
                })
            if (!is.null(mapping)) state$mapping[[database]] <- mapping
        })
        output$download_status <- shiny::renderUI({
            if (!is.null(state$download_error)) {
                return(integrationRefusalUi(state$download_error$sentence, state$download_error$detail, texts = texts))
            }
            mapping <- state$mapping[[input$database %||% ""]]
            if (is.null(mapping)) return(NULL)
            shiny::div(class = "alert alert-info", sprintf(texts[["enrich_downloaded"]], mapping$database, mapping$release,
                nrow(mapping$table)))
        })
        shiny::observeEvent(input$run, {
            state$error <- NULL
            manager <- workflow_data$state_manager
            if (is.null(manager) || is.null(workflow_data$integration_model)) {
                state$error <- list(sentence = texts[["enrich_not_trained"]], detail = NULL)
                return()
            }
            request <- integrationEnrichmentRequest(shiny::reactiveValuesToList(input))
            # The key is read here and handed to this call only; the field is emptied straight after.
            key <- input$string_key
            shiny::updateTextInput(session, "string_key", value = "")
            run <- tryCatch(withProgressFn(message = texts[["enrich_running"]], value = 0.5,
                runIntegrationEnrichment(manager, request, api_key = if (identical(request$database, "string")) key,
                    mapping = state$mapping[[request$database %||% ""]], query = query,
                    time_limit = INTEGRATION_ENRICHMENT_TAB_SECONDS)),
                error = function(error) {
                    state$error <- integrationRefusalText(error, texts)
                    NULL
                })
            key <- NULL
            if (is.null(run)) return()
            state$run <- run
            workflow_data$integration_enrichment <- list(model = workflow_data$integration_model$reference$entity_id,
                entry = integrationEnrichmentKey(run$entry$factor, run$entry$view, run$entry$database, run$entry$settings$test,
                    run$entry$settings$adjust_ruv_factors))
            status <- shiny::isolate(workflow_data$tab_status)
            status$integration_enrichment <- "complete"
            workflow_data$tab_status <- status
        })
        output$status <- shiny::renderUI({
            if (is.null(workflow_data$integration_model)) return(shiny::div(class = "alert alert-info", texts[["enrich_not_trained"]]))
            NULL
        })
        output$result <- shiny::renderUI({
            if (!is.null(state$error)) {
                return(integrationRefusalUi(sprintf(texts[["enrich_refused"]], state$error$sentence), state$error$detail,
                    texts = texts))
            }
            run <- state$run
            if (is.null(run)) return(NULL)
            shiny::div(class = "alert alert-success", lapply(integrationEnrichmentSummary(run$run$value, texts, run$entry),
                shiny::p))
        })
        output$terms <- DT::renderDT({
            run <- state$run
            shiny::req(!is.null(run))
            DT::datatable(as.data.frame(run$run$value$terms), rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
        })
        shiny::reactive(state$run)
    })
}

#' KIT-183: the test selector, pre-selected with the test of record for a data set of `n_samples` biological samples.
#' @param id The input's (namespaced) id.
#' @param n_samples The data set's biological samples.
#' @noRd
integrationEnrichmentTestInput <- function(id, n_samples, texts = integrationGuiTexts()) {
    shiny::selectInput(id, texts[["enrich_test"]], choices = stats::setNames(INTEGRATION_ENRICHMENT_TESTS_OFFERED,
        unname(texts[paste0("enrich_test_", INTEGRATION_ENRICHMENT_TESTS_OFFERED)])),
        selected = integrationEnrichmentTestOfRecord(n_samples))
}

#' KIT-183: the sentences beside the test selector: the selected test's note; for a data set with fewer biological
#' samples than the check covered, the few-samples sentence; for a model trained on RUV-corrected values, the RUV note
#' (KIT-184: for the tests that have no design, the pre-ranked and permutation tests; the others have the setting
#' "Adjust for removed RUV factors" and their own sentences with each result).
#' KIT-191: the note and the few-samples sentence of the test and of the data set's kind (lipids: the metabolites').
#' @param test The selected test (one of `INTEGRATION_ENRICHMENT_TESTS_OFFERED`).
#' @param n_samples The data set's biological samples.
#' @param input_matrix The values the model was trained on (`uncorrected` or `ruv_corrected`).
#' @param omic The data set's omic.
#' @param permutations VAL-115 (4): the number of reorderings chosen, named in the permutation test's note.
#' @noRd
integrationEnrichmentTestNotes <- function(test, n_samples, input_matrix, omic, texts = integrationGuiTexts(),
    permutations = INTEGRATION_ENRICHMENT_PERMUTATIONS[[1L]]) {
    if (!isTRUE(test %in% INTEGRATION_ENRICHMENT_TESTS_OFFERED)) return(character())
    kind <- integrationEnrichmentEvidenceKind(omic)
    note <- texts[[paste0("enrich_test_", test, "_note_", kind)]]
    if (identical(test, "permutation")) {
        note <- sprintf(note, format(as.numeric(permutations), big.mark = ",", scientific = FALSE, trim = TRUE))
    }
    c(note,
        if (isTRUE(n_samples < INTEGRATION_ENRICHMENT_CHECKED_SAMPLES))
            sprintf(texts[[paste0("enrich_test_few_samples_", test, "_", kind)]], as.integer(n_samples)),
        if (identical(input_matrix, "ruv_corrected") && !test %in% INTEGRATION_ENRICHMENT_DESIGN_TESTS)
            texts[["enrich_values_ruv_note"]])
}

#' KIT-191 (5.4, D6 (a)): the evidence table beside the selector: the four tests' figures at 12 biological samples for
#' the data set's kind (metabolite and lipid data sets: rows A, B, D, E, C3, I; gene and protein data sets: A, B, D, E,
#' G, H), with, for a lipid data set, the line saying no lipid data were simulated.
#' @param omic The data set's omic.
#' @noRd
integrationEnrichmentEvidenceUi <- function(omic, texts = integrationGuiTexts()) {
    kind <- integrationEnrichmentEvidenceKind(omic)
    tests <- INTEGRATION_ENRICHMENT_TESTS_OFFERED
    header <- shiny::tags$tr(shiny::tags$th(texts[["enrich_evidence_measure"]]),
        lapply(tests, function(test) shiny::tags$th(texts[[paste0("enrich_test_", test)]])))
    rows <- lapply(integrationEnrichmentEvidenceRows(kind), function(measure) {
        label <- if (measure == "E") paste0("enrich_evidence_E_", kind) else paste0("enrich_evidence_", measure)
        shiny::tags$tr(shiny::tags$td(texts[[label]]), lapply(tests, function(test)
            shiny::tags$td(integrationEnrichmentEvidenceValue(kind, test, measure))))
    })
    shiny::div(shiny::h5(sprintf(texts[["enrich_evidence_heading"]], texts[[paste0("enrich_evidence_kind_", kind)]])),
        if (identical(omic, "lipidomics")) shiny::p(shiny::em(texts[["enrich_evidence_lipids"]])),
        shiny::tags$table(class = "table table-condensed table-bordered", shiny::tags$thead(header), shiny::tags$tbody(rows)))
}

#' An enrichment's facts in words (KIT-141): what was enriched against which release, the counts, the mapping, the
#' direction and the number of tests. KIT-183: for a correlated test the counts it used and the relative direction;
#' for KEGG and Reactome, whatever the test, the caution and, with the accepted entry, the factor's strength and the
#' RUV note for a model trained on RUV-corrected values.
#' KIT-184 (5.5): for a test that has a design, the label of its form (from the product) and that form's sentences
#' (`integrationEnrichmentAdjustmentSentences()`); the RUV note stays with the pre-ranked and permutation tests.
#' @param value The enrichment product.
#' @param entry The entry accepted on the head (`runIntegrationEnrichment()$entry`), or NULL: the factor-strength
#'   sentence and the sentences that depend on the values the model was trained on are left out.
#' @return A character vector of sentences.
#' @noRd
integrationEnrichmentSummary <- function(value, texts = integrationGuiTexts(), entry = NULL) {
    service <- value$service
    release <- if (identical(value$database, "string")) sprintf("STRING %s", service$version$string_version %||% "?") else
        as.character(service$release)
    r <- value$ranking_counts
    m <- value$mapping_counts
    d <- value$direction
    # AUDIT-KIT-191-impl I04: exactly `test` (a pre-ranked product has `tests`, the count, and no `test`).
    test <- value[["test", exact = TRUE]]
    correlated <- isTRUE(test %in% INTEGRATION_ENRICHMENT_CORRELATED_TESTS)
    permutation <- identical(test, "permutation")
    pathways <- !identical(value$database, "string")
    counts <- value$test_counts
    c(sprintf(texts[["enrich_done"]], value$factor, value$view, integrationEnrichmentDatabaseName(value$database), release),
        sprintf(texts[["enrich_counts"]], as.integer(r$features), as.integer(r$ranked), as.integer(r$too_few_samples),
            as.integer(r$no_variation), as.integer(r$internal_standard)),
        sprintf(texts[["enrich_mapping"]], as.integer(m$mapped), as.integer(m$unmapped), as.integer(m$duplicated)),
        if (correlated) sprintf(texts[[if (permutation) "enrich_test_counts_permutation" else "enrich_test_counts"]],
            as.integer(counts$universe), as.integer(counts$samples), as.integer(counts$incomplete),
            as.integer(counts$collapsed)),
        if (correlated) paste(sprintf(texts[["enrich_direction_oriented"]], as.character(d$feature_id), as.character(d$view),
            as.integer(d$sign)), texts[["enrich_direction_camera"]]) else
            sprintf(texts[["enrich_direction"]], as.character(d$feature_id), as.character(d$view), as.integer(d$sign)),
        if (as.integer(value$tests) > 0L) sprintf(texts[["enrich_tests"]], as.integer(value$tests)) else texts[["enrich_none"]],
        # KIT-191 (5.2 step 8): the permutation test's floor of p-values and of BH.
        if (permutation) integrationEnrichmentPermutationFloorText(value$permutation, as.integer(value$tests), texts),
        if (pathways) texts[["enrich_caution"]],
        if (pathways && !is.null(entry)) {
            strength <- integrationEnrichmentStrengthText(entry, texts)
            if (nzchar(strength)) strength
        },
        if (pathways && isTRUE(test %in% INTEGRATION_ENRICHMENT_DESIGN_TESTS))
            integrationEnrichmentAdjustmentSentences(value, entry, texts),
        if (pathways && !isTRUE(test %in% INTEGRATION_ENRICHMENT_DESIGN_TESTS) &&
            identical(entry$input_matrix, "ruv_corrected")) texts[["enrich_values_ruv_note"]])
}

#' KIT-184 (5.5; VAL-58 (5), (6), VAL-62): the sentences of a result of a test that has a design, by its form. The form
#' is read from the product (`integrationEnrichmentEntryForm()`), the values the model was trained on from the entry.
#' Always the label. With the removed RUV factors in the design: the note for the model's input; the design, the features
#' measured in every run and the share of the factor the removed factors explain beside the share expected by chance; for
#' an approximate set of factors the share of what was removed that lies within them, and no sentence claiming an exact
#' match; the overlap sentence when the adjusted share is at least 0.5; the few-degrees sentence below the smallest
#' count checked; for `camera_floor` its own sentence. Without them (the default) on a model trained on RUV-corrected
#' values: the note of that form.
#' @param value The enrichment product of a test that has a design.
#' @param entry The entry accepted on the head (`input_matrix`), or NULL.
#' @return A character vector of sentences.
#' @noRd
integrationEnrichmentAdjustmentSentences <- function(value, entry = NULL, texts = integrationGuiTexts()) {
    form <- integrationEnrichmentEntryForm(value, texts)
    trained_corrected <- identical(entry$input_matrix, "ruv_corrected")
    if (!form$adjusted) {
        return(c(form$label, if (trained_corrected) texts[["enrich_values_ruv_unadjusted_note"]]))
    }
    detail <- value$design_detail
    q <- form$factors
    n <- as.integer(value$test_counts$samples)
    residual_df <- as.integer(detail$residual_df)
    exact <- identical(as.character(detail$status), "exact")
    percent <- function(x) integrationEnrichmentPercent(100 * as.numeric(x))
    share <- sprintf(texts[["enrich_test_covariates_share"]], percent(detail$factor_r2), percent(q / (n - 1L)))
    design <- sprintf(texts[["enrich_test_covariates_design"]], q, n, residual_df)
    c(form$label,
        if (is.null(entry)) NULL else if (!trained_corrected) sprintf(texts[["enrich_values_adjusted_default_note"]], q) else
            if (exact) sprintf(texts[["enrich_values_ruv_covariates_note"]], q, q, residual_df) else
                sprintf(texts[["enrich_values_ruv_covariates_note_approximate"]], q, residual_df),
        if (!exact) sprintf(texts[["enrich_removed_approximate"]], as.integer(detail$k_declared), percent(detail$share)),
        if (exact) paste(design, sprintf(texts[["enrich_test_covariates_runs"]], as.integer(detail$run_complete),
            as.integer(value$test_counts$universe)), share) else paste(design, share),
        if (isTRUE(integrationEnrichmentAdjustedShare(as.numeric(detail$factor_r2), n, q) >= 0.5))
            texts[["enrich_factor_removed_overlap"]],
        if (isTRUE(as.integer(detail$n_eq) < INTEGRATION_ENRICHMENT_CHECKED_N_EQ)) sprintf(texts[["enrich_test_few_df"]], q,
            residual_df, INTEGRATION_ENRICHMENT_CHECKED_N_EQ - 2L),
        if (identical(value[["test", exact = TRUE]], "camera_floor")) texts[["enrich_floor_covariates_note"]])
}

#' KIT-191 (5.2 steps 7 and 8, corrected by AUDIT-KIT-193 A05): the permutation test's floor. With d reorderings that
#' left the factor unchanged no p-value is below p_min = (1 + d) / (1 + B), printed exactly as a fraction
#' (`integrationEnrichmentFloorFraction()`; AUDIT-KIT-191-impl I01: a decimal rounded to nearest can exceed it); with m
#' terms tested, BH at 5 % calls none unless it calls at least k_min = ceiling(m p_min / 0.05) together (it calls k*
#' terms with p_(k*) <= 0.05 k* / m, and p_(k*) >= p_min), k_min computed in whole numbers as
#' ceiling(20 m (1 + d) / (1 + B)).
#' @param permutation The product's `permutation` (`permutations`, `identity_draws`, `samples`).
#' @param terms The number of terms tested (m).
#' @return `enrich_test_permutation_floor`, or `enrich_test_permutation_floor_none` when p_min > 0.05 (no term can be
#'   called); NULL without terms.
#' @noRd
integrationEnrichmentPermutationFloorText <- function(permutation, terms, texts = integrationGuiTexts()) {
    if (is.null(permutation) || !isTRUE(terms > 0L) || is.na(permutation$identity_draws)) return(NULL)
    reorderings <- as.numeric(permutation$permutations)
    unchanged <- as.numeric(permutation$identity_draws)
    floor_text <- integrationEnrichmentFloorFraction(unchanged, reorderings)
    if (20 * (1 + unchanged) > 1 + reorderings) {
        return(sprintf(texts[["enrich_test_permutation_floor_none"]], as.integer(reorderings),
            as.integer(permutation$samples), floor_text))
    }
    k_min <- integrationEnrichmentFloorMinimum(unchanged, reorderings, terms)
    sprintf(texts[["enrich_test_permutation_floor"]], as.integer(reorderings), as.integer(permutation$samples), floor_text,
        as.integer(terms), k_min)
}
