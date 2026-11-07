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

# KIT-103 W4: the shared design builder (kit-103 design 4.5), one module for every omic. The UI is
# the lipid builder's (mod_lipid_design_builder_ui.R at 616e1147), its 38 ids unchanged, with the
# omic's texts from `config$design`; the server is the lipid builder's decomposition as one module
# (mod_lipid_design_builder_{server,state_helpers,action_helpers,display_helpers,helpers}.R at
# 616e1147), its observers and outputs created in the same order, over the design core
# (func_omics_design_shared.R: C3 initial state, C4 edits, C5 payload). The data reach it through
# the configured adapter (`config$design$data_adapter`): `wide_assay_list` (lipid, metab, trx) or,
# from KIT-103 W6, `dia_long` (proteomics: its DIA references, initial state, owner sync, renames,
# texts and save kept by name; its builder UI stays its own plugin). testServer masking (KIT-102
# W4): the module function keeps its reactive named `design_matrix` (MCI-027 reads it) and holds no
# other local name that a test's own variables would meet; the configuration is read outside it.
# The `DEBUG66` console messages are not kept (KIT-103 D5). Recorded divergence from
# conventions/R-shiny.md: the fluidPage/wellPanel layout is kept (the byte-identical contract),
# under the existing-code scope rule of conventions/R.md.

# ------------------------------------------------------------------------------------------------
# The UI.
# ------------------------------------------------------------------------------------------------

#' The shared design builder's UI: the rename, removal, factor, metadata, technical-replicate,
#' contrast, formula and settings tabs, the information panel and the design table.
#' @param id Module id (the design tab's `builder`).
#' @param config An omic GUI configuration (`omicGuiConfig()`); its `design$multi_assay_note` and
#'   `design$matrix_heading` are the omic's texts.
#' @noRd
mod_omic_design_builder_ui <- function(id, config) {
    ns <- shiny::NS(id)
    design <- config$design
    shiny::fluidPage(
        shiny::titlePanel("Design Matrix Builder"),
        shiny::fluidRow(
            shiny::column(4,
                shiny::wellPanel(omicsDesignBuilderTabs(ns)),
                shiny::wellPanel(
                    shiny::h4("Information & Assignments"),
                    shiny::helpText(shiny::HTML(design$multi_assay_note,
                        "(e.g., LCMS_Pos, LCMS_Neg, GCMS). All assays share the SAME design matrix - ",
                        "sample assignments, factors, and contrasts apply uniformly across all assays.")),
                    shiny::hr(),
                    shiny::helpText(shiny::HTML("<b>Factors:</b> Experimental variables distinguishing your samples ",
                        "(e.g., Treatment, Timepoint, Genotype).")),
                    shiny::helpText(shiny::HTML(
                        "<b>Technical Replicates:</b> Repeated measurements of the same biological sample.")),
                    shiny::helpText(shiny::HTML(
                        "<b>Contrasts:</b> Statistical comparisons between groups for differential analysis.")),
                    shiny::hr(),
                    shiny::h5("Available Factors"),
                    shiny::uiOutput(ns("available_factors_display")),
                    shiny::hr(),
                    shiny::h5("Defined Contrasts"),
                    shiny::uiOutput(ns("defined_contrasts_display")))),
            shiny::column(8,
                shiny::wellPanel(
                    shiny::h4(design$matrix_heading),
                    DT::DTOutput(ns("data_table"))),
                shiny::actionButton(ns("save_results"), "Save Design", class = "btn-primary",
                    icon = shiny::icon("save"), style = "float: right;"))))
}

#' The builder's tabset (`main_tabs`).
#' @noRd
omicsDesignBuilderTabs <- function(ns) {
    shiny::tabsetPanel(id = ns("main_tabs"),
        shiny::tabPanel("Rename Samples",
            shiny::h4("Individual Rename"),
            shiny::fluidRow(
                shiny::column(6, shiny::selectizeInput(ns("sample_to_rename"), "Select Sample:", choices = NULL)),
                shiny::column(6, shiny::textInput(ns("new_sample_name"), "New Name:"))),
            shiny::actionButton(ns("rename_sample"), "Rename"),
            shiny::hr(),
            shiny::h4("Bulk Rename"),
            shiny::selectizeInput(ns("samples_to_transform"), "Select Samples:", choices = NULL, multiple = TRUE),
            shiny::radioButtons(ns("transform_mode"), "Transformation Mode:", choices = c(
                "Keep Text BEFORE First Underscore" = "before_underscore",
                "Keep Text AFTER Last Underscore" = "after_underscore",
                "Extract by Position" = "range")),
            shiny::helpText(shiny::HTML(
                "Select how to shorten the selected sample names using underscores as delimiters.",
                "<br><b>Keep Text BEFORE First Underscore:</b> Retains only the part before the very first '_'. Example: 'SampleA_T1_Rep1' &rarr; 'SampleA'.",
                "<br><b>Keep Text AFTER Last Underscore:</b> Retains only the part after the very last '_'. Example: 'SampleA_T1_Rep1' &rarr; 'Rep1'.",
                "<br><b>Extract by Underscore Position:</b> Keeps text segments based on underscore indices (1-based).")),
            shiny::conditionalPanel(condition = paste0("input['", ns("transform_mode"), "'] == 'range'"),
                shiny::numericInput(ns("range_start"), "Start Underscore Index:", 1, min = 1),
                shiny::numericInput(ns("range_end"), "End Underscore Index:", 1, min = 1),
                shiny::h5("Preview (First Selected Sample):"),
                shiny::verbatimTextOutput(ns("range_preview"))),
            shiny::actionButton(ns("bulk_rename"), "Apply Transformation")),
        shiny::tabPanel("Remove Samples",
            shiny::h4("Remove Samples from Analysis"),
            shiny::helpText(shiny::HTML("Select samples you want to exclude from the analysis. ",
                "Removed samples will not appear in the design matrix or be included in downstream analysis.")),
            shiny::selectizeInput(ns("samples_to_remove"), "Select Samples to Remove:", choices = NULL,
                multiple = TRUE),
            shiny::actionButton(ns("remove_samples"), "Remove Selected Samples", class = "btn-danger",
                icon = shiny::icon("trash")),
            shiny::hr(),
            shiny::h5("Currently Removed Samples"),
            shiny::verbatimTextOutput(ns("removed_samples_display"))),
        shiny::tabPanel("Factors",
            shiny::h4("Add New Factor"),
            shiny::textInput(ns("new_factor"), "New Factor Name:"),
            shiny::actionButton(ns("add_factor"), "Add Factor")),
        shiny::tabPanel("Assign Metadata",
            shiny::h4("Assign Metadata"),
            shiny::selectizeInput(ns("selected_runs"), "Select Samples:", choices = NULL, multiple = TRUE),
            shiny::selectInput(ns("factor1_select"), "Select Factor 1:", choices = c("")),
            shiny::selectInput(ns("factor2_select"), "Select Factor 2:", choices = c("")),
            shiny::selectInput(ns("factor3_select"), "Select Factor 3:", choices = c("")),
            shiny::uiOutput(ns("replicate_inputs")),
            shiny::actionButton(ns("assign_metadata"), "Assign")),
        shiny::tabPanel("Technical Replicates",
            shiny::h4("Assign Technical Replicates"),
            shiny::helpText(shiny::HTML("Select samples that are technical replicates of each other. ",
                "These samples will be grouped under the same biological replicate number.")),
            shiny::selectizeInput(ns("tech_rep_samples"), "Select Technical Replicates:", choices = NULL,
                multiple = TRUE),
            shiny::radioButtons(ns("tech_rep_assignment_mode"), "Assignment Mode:", choices = c(
                "Use lowest replicate number" = "lowest",
                "Use first selected sample's replicate" = "first",
                "Specify replicate number" = "manual")),
            shiny::conditionalPanel(condition = paste0("input['", ns("tech_rep_assignment_mode"), "'] == 'manual'"),
                shiny::numericInput(ns("manual_replicate_number"), "Replicate Number:", value = 1, min = 1)),
            shiny::actionButton(ns("assign_tech_reps"), "Assign Technical Replicates", class = "btn-primary"),
            shiny::hr(),
            shiny::h5("Current Technical Replicate Assignments"),
            shiny::verbatimTextOutput(ns("tech_rep_summary"))),
        shiny::tabPanel("Contrasts",
            shiny::h4("Define Contrasts"),
            shiny::selectInput(ns("contrast_group1"), "Group 1:", choices = c("")),
            shiny::selectInput(ns("contrast_group2"), "Group 2:", choices = c("")),
            shiny::verbatimTextOutput(ns("contrast_factors_info")),
            shiny::actionButton(ns("add_contrast"), "Add Contrast")),
        shiny::tabPanel("Formula",
            shiny::h4("Model Formula"),
            shiny::textInput(ns("formula_string"), "Formula:", value = "~ 0 + group"),
            shiny::helpText("Default: '~ 0 + group' creates groupX-groupY contrast format for limma.")),
        shiny::tabPanel("Settings",
            shiny::h4("Session Management"),
            shiny::p("Reset options allow you to revert changes made during this session."),
            shiny::div(style = "margin-top: 15px;",
                shiny::selectInput(ns("reset_scope"), "Reset Scope:", choices = c(
                    "All Changes" = "all", "Sample Names Only" = "sample_names",
                    "Removed Samples Only" = "removed_samples", "Factors Only" = "factors",
                    "Contrasts Only" = "contrasts", "Formula Only" = "formula"), selected = "all")),
            shiny::div(style = "margin-top: 15px;",
                shiny::actionButton(ns("reset_changes"), "Reset to Initial State", class = "btn-warning",
                    icon = shiny::icon("rotate-left"), width = "100%")),
            shiny::hr(),
            shiny::h4("Information"),
            shiny::verbatimTextOutput(ns("session_info"))))
}

# ------------------------------------------------------------------------------------------------
# The data adapter.
# ------------------------------------------------------------------------------------------------

#' The builder's data calls by the configured adapter: `wide_assay_list` (lipid, metab, trx) or
#' `dia_long` (proteomics, `omicsDesignDiaLongAdapter()`).
#' @return A list of functions:
#'   - `holders(currentOwner)`: the saved result's and the working tables' holders (`result`,
#'     `data_cln`, each with a `reactiveVal`'s get/set contract);
#'   - `initial_state(dataTbl, configList, columnMapping, existingDesignMatrix, existingContrasts,
#'     currentOwner)`: the initial state's reactive;
#'   - `on_start(input, session, dataTbl, initialState, state)`: the observer that loads the initial
#'     state into the builder's values and sample inputs;
#'   - `scope(operation)`: runs an edit or the save in the data owner's scope;
#'   - `read_data(holder)`: the working tables as the edits and the save take them;
#'   - `rename(designMatrix, data, sampleNames, renameMap)`, `assign_metadata(...)` and
#'     `tech_reps(...)` (C4's signatures);
#'   - `reset_modal(resetScope, ns)`, `reset(scope, initialState, state, updateFormula)`;
#'   - `save(designMatrix, removedSamples, data, colMap, contrastData, configList, formulaString)`:
#'     the payload, or NULL when no run is assigned;
#'   - `keep_result(value, owner)`: the payload as the result holder keeps it.
#' @noRd
omicsDesignBuilderAdapter <- function(config) {
    switch(config$design$data_adapter,
        wide_assay_list = omicsDesignWideAdapter(config),
        dia_long = omicsDesignDiaLongAdapter(),
        cli::cli_abort("the {.val {config$design$data_adapter}} design adapter is unknown.",
            class = "omic_gui_config_invalid"))
}

#' `wide_assay_list` (lipid, metab, trx): the core's initial state, edits and payload over a list
#' of wide assay tables, held by reference in the owner when there is one (KIT-075).
#' @noRd
omicsDesignWideAdapter <- function(config) {
    list(
        holders = function(currentOwner) {
            list(result = shiny::reactiveVal(NULL), data_cln = omicsBuilderTableVal("design_builder_data_cln", currentOwner))
        },
        initial_state = function(dataTbl, configList, columnMapping, existingDesignMatrix, existingContrasts,
                                 currentOwner) {
            shiny::reactive({
                shiny::req(dataTbl())
                shiny::req(configList())
                state <- buildOmicsDesignInitialState(dataTbl(), configList(), columnMapping(), config,
                    existingDesignMatrix = existingDesignMatrix, existingContrasts = existingContrasts)
                # KIT-075: a reactive keeps its value; it keeps the tables' reference.
                state$data_cln <- saveWorkflowOwnedValue(state$data_cln, "design_builder_initial_data_cln", currentOwner())
                state
            })
        },
        on_start = registerOmicsDesignInitialState,
        scope = function(operation) operation(),
        read_data = function(holder) holder(),
        rename = applyOmicsDesignSampleRenameMap,
        assign_metadata = applyOmicsDesignMetadataAssignment,
        tech_reps = applyOmicsDesignTechRepAssignment,
        reset_modal = omicsDesignResetModal,
        reset = applyOmicsDesignResetState,
        save = function(designMatrix, removedSamples, dataCln, colMap, contrastData, configList, formulaString) {
            buildOmicsDesignSavePayload(designMatrix, removedSamples, dataCln, colMap, contrastData, configList,
                formulaString, idKey = config$design$id_mapping_key)
        },
        keep_result = function(value, owner) {
            # KIT-075: the returned result carries the tables' reference.
            value$data_cln <- saveWorkflowOwnedValue(value$data_cln, "design_builder_result_data_cln", owner)
            value
        })
}

#' `dia_long` (proteomics; the DIA adaptor's design entries): the long peptide or protein table
#' held as a DIA reference and read only by the owner's transforms (KIT-075/082). Today's
#' proteomics builder pieces, names kept: the saved holders (`createProtDesignBuilderState`), the
#' two initial-state paths (`createProtDesignInitialStateReactive`: a DIA reference through
#' `buildDIAReferenceDesignState`, a plain table through `buildProtDesignInitialState`), the owner
#' replacement and restore (`registerProtDesignInitialStateSyncObserver`) with the checkpoint's
#' design reader, the renames (`applyProtDesign*RenameUpdate*`), the integer replicate numbers and
#' the technical-replicate and reset texts, and the save's DIA handling
#' (`runProtDesignSaveResultsObserverShell`). The edits and the save run in the owner's scope.
#' @noRd
omicsDesignDiaLongAdapter <- function() {
    list(
        holders = function(currentOwner) {
            builder <- createProtDesignBuilderState()
            list(result = builder$resultRv, data_cln = builder$dataClnReactive)
        },
        initial_state = function(dataTbl, configList, columnMapping, existingDesignMatrix, existingContrasts,
                                 currentOwner) {
            createProtDesignInitialStateReactive(dataTbl, configList, columnMapping)
        },
        on_start = function(input, session, dataTbl, initialState, state) {
            registerProtDesignInitialStateSyncObserver(input, session, initialState, dataTbl,
                designMatrix = state$design_matrix, dataClnReactive = state$data_cln, groups = state$groups,
                factors = state$factors, contrasts = state$contrasts, removedSamples = state$removed_samples,
                resultRv = state$result)
            registerProtDesignCheckpointReader(input, session, state)
        },
        scope = function(operation) withDIAStateScope(getDIAActiveManager(), operation),
        read_data = getDIAReactiveReference,
        rename = function(designMatrix, dataCln, sampleNames, renameMap) {
            renamed <- applyProtDesignBulkRenameUpdates(designMatrix, dataCln, names(renameMap), unname(renameMap))
            list(designMatrix = renamed$designMatrix, assayList = renamed$dataCln, sampleNames = sampleNames)
        },
        assign_metadata = applyProtDesignMetadataAssignment,
        tech_reps = applyProtDesignTechRepAssignment,
        reset_modal = protDesignResetModal,
        reset = function(scope, initialState, state, updateFormula) {
            applyProtDesignResetState(scope, initialState, state$design_matrix, state$design_matrix, state$data_cln,
                state$removed_samples(), state$removed_samples, state$factors, state$groups,
                state$contrasts, updateFormula)
        },
        save = function(designMatrix, removedSamples, dataCln, colMap, contrastData, configList, formulaString) {
            # The builder shows the texts and keeps the result.
            runProtDesignSaveResultsObserverShell(designMatrix, removedSamples, dataCln, contrastData, configList,
                formulaString, resultSetter = function(value) invisible(NULL),
                showNotification = function(...) invisible(NULL))
        },
        keep_result = function(value, owner) value)
}

# ------------------------------------------------------------------------------------------------
# The display formatters.
# ------------------------------------------------------------------------------------------------

#' The technical replicate assignments, by group and biological replicate.
#' @noRd
formatOmicsDesignTechRepSummary <- function(designMatrix) {
    techRepGroups <- designMatrix |>
        dplyr::filter(!is.na(tech_reps)) |>
        dplyr::group_by(group, replicates) |>
        dplyr::summarise(samples = paste(Run, collapse = ", "), techRepNumbers = paste(tech_reps, collapse = ", "),
            .groups = "drop")
    if (nrow(techRepGroups) == 0) return("No technical replicates assigned yet.")
    outputText <- techRepGroups |>
        dplyr::mutate(formatted = sprintf(
            "Group: %s, Biological Replicate: %s\n  Samples: %s\n  Technical Replicates: %s",
            group, replicates, samples, techRepNumbers)) |>
        dplyr::pull(formatted)
    paste(outputText, collapse = "\n\n")
}

#' The removed runs, natural-sorted.
#' @noRd
formatOmicsDesignRemovedSamples <- function(removedSamples) {
    if (length(removedSamples) == 0) return("No samples have been removed.")
    paste0("Removed ", length(removedSamples), " sample(s):\n", paste(gtools::mixedsort(removedSamples), collapse = "\n"))
}

#' The contrast note: the omic's text for a `~ 0 + group` formula, else its as-is text. A NULL
#' formula is an error, as today (the output shows it until the formula input arrives).
#' @noRd
formatOmicsDesignContrastFactorsInfo <- function(formulaString, builderTexts) {
    if (grepl("~ *0 *\\+ *group", formulaString)) {
        builderTexts[["contrast_group_prefix"]]
    } else {
        builderTexts[["contrast_as_is"]]
    }
}

#' The range rename's preview of the first selected sample.
#' @noRd
formatOmicsDesignRangePreview <- function(samplesToTransform, rangeStart, rangeEnd) {
    firstSample <- samplesToTransform[1]
    tryCatch({
        previewResult <- extract_experiment(firstSample, mode = "range", start = rangeStart, end = rangeEnd)
        paste0("\"", firstSample, "\" -> \"", previewResult, "\"")
    }, error = function(error) {
        paste("Error:", error$message)
    })
}

# ------------------------------------------------------------------------------------------------
# The reset.
# ------------------------------------------------------------------------------------------------

#' The reset's confirmation modal.
#' @noRd
omicsDesignResetModal <- function(resetScope, ns) {
    shiny::modalDialog(
        title = "Confirm Reset",
        shiny::HTML(paste0("<p>This will revert <strong>", resetScope, "</strong> to their initial state.</p>")),
        footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_reset"), "Reset", class = "btn-danger")),
        easyClose = TRUE)
}

#' Revert the builder's state, by scope, to the initial state.
#' @param state The builder's reactive values (`design_matrix`, `data_cln`, `sample_names`,
#'   `removed_samples`, `factors`, `groups`, `contrasts`).
#' @noRd
applyOmicsDesignResetState <- function(scope, initialState, state, updateFormula) {
    if (scope == "all" || scope == "sample_names") {
        state$design_matrix(initialState$design_matrix)
        state$data_cln(initialState$data_cln)
        state$sample_names(initialState$sample_names)
    }
    if (scope == "all" || scope == "removed_samples") state$removed_samples(character(0))
    if (scope == "all" || scope == "factors") {
        state$factors(initialState$factors)
        currentMatrix <- state$design_matrix()
        currentMatrix$factor1 <- NA_character_
        currentMatrix$factor2 <- NA_character_
        currentMatrix$factor3 <- NA_character_
        currentMatrix$group <- NA_character_
        currentMatrix$tech_reps <- NA_integer_
        state$design_matrix(currentMatrix)
        state$groups(initialState$groups)
    }
    if (scope == "all" || scope == "contrasts") state$contrasts(initialState$contrasts)
    if (scope == "all" || scope == "formula") updateFormula(initialState$formula)
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# The server.
# ------------------------------------------------------------------------------------------------

#' The shared design builder's server.
#' @param id Module id.
#' @param data_tbl A function (or reactive) returning the assay tables, read per use (KIT-075). For
#'   `dia_long`, the table's provider with its `dia_reference` attribute (`registerProtDesignBuilderModule`),
#'   or a plain reactive (the standalone applet).
#' @param config_list,column_mapping Reactives of the workflow configuration and column mapping.
#' @param existing_design_matrix,existing_contrasts Optional reactives of an existing design (an
#'   imported one) and its contrasts: a design with an assigned group is the initial state
#'   (`wide_assay_list` only).
#' @param state_owner Optional function returning the session's state owner; the builder then holds
#'   its tables by reference (KIT-075). Without it they stay plain values. (`dia_long` reads the
#'   session's DIA owner itself.)
#' @param config An omic GUI configuration (`omicGuiConfig()`): the run order, the id key, the texts,
#'   the level order and excluded runs (`record_excluded_runs`, `group_levels`).
#' @return A reactiveVal holding the saved result (`design_matrix`, `data_cln` held by reference,
#'   `contrasts_tbl`, `config_list`; with `group_levels` and `excluded_runs` when configured). For
#'   `dia_long`, proteomics' saved holder, whose result names the table by `data_reference` under a
#'   DIA owner.
#' @noRd
mod_omic_design_builder_server <- function(id, data_tbl, config_list, column_mapping, existing_design_matrix = NULL,
                                           existing_contrasts = NULL, state_owner = NULL, config) {
    adapter <- omicsDesignBuilderAdapter(config)
    builder_texts <- config$design$builder_texts
    replicate_testid <- unname(config$design$testids["replicate_start"])
    record_excluded_runs <- isTRUE(config$design$record_excluded_runs)
    group_levels <- config$design$group_levels
    shiny::moduleServer(id, function(input, output, session) {
        # KIT-075: the builder's table lists are held by reference in the owner (plain without one).
        currentOwner <- function() if (is.function(state_owner)) state_owner()

        # The final result's and the working tables' holders (the adapter's).
        holders <- adapter$holders(currentOwner)
        result_rv <- holders$result

        initial_state <- adapter$initial_state(data_tbl, config_list, column_mapping, existing_design_matrix,
            existing_contrasts, currentOwner)

        # Reactive values for the current, mutable state (MCI-027 reads `design_matrix` by name).
        design_matrix <- shiny::reactiveVal()
        data_cln_reactive <- holders$data_cln
        sample_names_reactive <- shiny::reactiveVal()
        groups <- shiny::reactiveVal()
        factors <- shiny::reactiveVal()
        contrasts <- shiny::reactiveVal()
        removed_samples <- shiny::reactiveVal(character(0))

        registerOmicsDesignBuilderObservers(input, output, session, data_tbl, initial_state, column_mapping,
            config_list, list(design_matrix = design_matrix, data_cln = data_cln_reactive,
                sample_names = sample_names_reactive, groups = groups, factors = factors, contrasts = contrasts,
                removed_samples = removed_samples, result = result_rv),
            adapter, builder_texts, replicate_testid, function(value) {
                # The level order the user entered and the runs the user removed (trx).
                if (record_excluded_runs) {
                    value$group_levels <- group_levels(value$design_matrix, factors())
                    value$excluded_runs <- removed_samples()
                }
                result_rv(adapter$keep_result(value, currentOwner()))
            })

        result_rv
    })
}

#' The builder's observers and outputs, created in today's order: the initial state, the sample
#' and factor dropdowns, the design table, the summaries, the adjacent outputs, the rename,
#' metadata, technical-replicate, contrast and removal edits, the reset and the save.
#' @param state The builder's reactive values (as `applyOmicsDesignResetState`) and the result's
#'   holder (`result`).
#' @param setResult The function taking the saved payload.
#' @noRd
registerOmicsDesignBuilderObservers <- function(input, output, session, dataTbl, initialState, columnMapping,
                                                configList, state, adapter, builderTexts, replicateTestid,
                                                setResult) {
    adapter$on_start(input, session, dataTbl, initialState, state)
    proxyDataTable <- DT::dataTableProxy("data_table")
    registerOmicsDesignSelectionInputs(input, session, state)
    registerOmicsDesignOutputs(input, output, session, proxyDataTable, state, builderTexts, replicateTestid)
    registerOmicsDesignEdits(input, session, state, adapter)
    shiny::observeEvent(input$reset_changes, {
        shiny::showModal(adapter$reset_modal(input$reset_scope, session$ns))
    })
    shiny::observeEvent(input$confirm_reset, {
        scope <- input$reset_scope
        adapter$reset(scope, initialState(), state, function(value) {
            shiny::updateTextInput(session, "formula_string", value = value)
        })
        shiny::removeModal()
        shiny::showNotification(paste("Reset of", scope, "completed."), type = "message")
    })
    shiny::observeEvent(input$save_results, adapter$scope(function() {
        finalResult <- adapter$save(state$design_matrix(), state$removed_samples(), adapter$read_data(state$data_cln),
            columnMapping(), state$contrasts(), configList(), input$formula_string)
        if (is.null(finalResult)) {
            shiny::showNotification(builderTexts[["no_assignment"]], type = "warning")
            return(NULL)
        }
        setResult(finalResult)
        shiny::showNotification(builderTexts[["save_success"]], type = "message", duration = 5)
    }))
    invisible(NULL)
}

#' The initial state: the reactive values and the sample inputs, on the data and on each reset.
#' @noRd
registerOmicsDesignInitialState <- function(input, session, dataTbl, initialState, state) {
    initializationObserver <- shiny::observe({
        initial <- initialState()
        shiny::req(initial)
        state$design_matrix(initial$design_matrix)
        state$data_cln(initial$data_cln)
        state$sample_names(initial$sample_names)
        state$groups(initial$groups)
        state$factors(initial$factors)
        state$contrasts(initial$contrasts)
        state$removed_samples(character(0))
        sortedRuns <- initial$design_matrix$Run
        for (inputId in c("sample_to_rename", "selected_runs", "samples_to_transform", "tech_rep_samples",
            "samples_to_remove")) {
            shiny::updateSelectizeInput(session, inputId, choices = sortedRuns, selected = "")
        }
        shiny::updateTextInput(session, "formula_string", value = initial$formula)
    })
    shiny::bindEvent(initializationObserver, dataTbl(), input$reset_changes)
    invisible(NULL)
}

#' The sample selections (the runs not removed, natural-sorted, each selection kept where it is
#' still offered) and the factor and group dropdowns (the selections kept).
#' @noRd
registerOmicsDesignSelectionInputs <- function(input, session, state) {
    shiny::observe({
        currentDesignMatrix <- state$design_matrix()
        shiny::req(currentDesignMatrix)
        allRuns <- gtools::mixedsort(currentDesignMatrix$Run)
        availableRuns <- allRuns[!allRuns %in% state$removed_samples()]
        inputIds <- c("sample_to_rename", "selected_runs", "samples_to_transform", "tech_rep_samples",
            "samples_to_remove")
        # Each selection kept where it is still offered (a NULL selection stays NULL); the single
        # rename choice is "" otherwise.
        selected <- shiny::isolate(lapply(stats::setNames(inputIds, inputIds), function(inputId) input[[inputId]]))
        rename <- selected$sample_to_rename
        selected$sample_to_rename <- if (!is.null(rename) && rename %in% availableRuns) rename else ""
        for (inputId in inputIds) {
            choice <- if (identical(inputId, "sample_to_rename")) {
                selected[[inputId]]
            } else {
                selected[[inputId]][selected[[inputId]] %in% availableRuns]
            }
            shiny::updateSelectizeInput(session, inputId, choices = availableRuns, selected = choice)
        }
    })
    shiny::observe({
        selected <- shiny::isolate(list(factor1_select = input$factor1_select, factor2_select = input$factor2_select,
            factor3_select = input$factor3_select, contrast_group1 = input$contrast_group1,
            contrast_group2 = input$contrast_group2))
        currentFactors <- state$factors()
        currentGroups <- state$groups()
        for (inputId in names(selected)) {
            choices <- c("", if (startsWith(inputId, "factor")) currentFactors else currentGroups)
            shiny::updateSelectInput(session, inputId, choices = choices, selected = selected[[inputId]])
        }
    })
    invisible(NULL)
}

#' The builder's outputs: the design table (rendered, then replaced in place), the technical
#' replicate, removal and contrast-note texts, the factor and contrast lists, the range preview and
#' the replicate-start input.
#' @noRd
registerOmicsDesignOutputs <- function(input, output, session, proxyDataTable, state, builderTexts, replicateTestid) {
    activeTable <- function(designMatrix) designMatrix |> dplyr::filter(!Run %in% state$removed_samples())
    output$data_table <- DT::renderDT({
        currentDesignMatrix <- state$design_matrix()
        shiny::req(currentDesignMatrix)
        activeTable(currentDesignMatrix)
    }, selection = "none", options = list(pageLength = 10, scrollX = TRUE, server = FALSE))
    shiny::observe({
        shiny::req(proxyDataTable)
        currentDesignMatrix <- state$design_matrix()
        shiny::req(currentDesignMatrix)
        DT::replaceData(proxyDataTable, activeTable(currentDesignMatrix), resetPaging = FALSE)
    })
    output$tech_rep_summary <- shiny::renderText({
        currentDesignMatrix <- state$design_matrix()
        shiny::req(currentDesignMatrix)
        formatOmicsDesignTechRepSummary(currentDesignMatrix)
    })
    output$removed_samples_display <- shiny::renderText({
        formatOmicsDesignRemovedSamples(state$removed_samples())
    })
    output$contrast_factors_info <- shiny::renderText({
        formatOmicsDesignContrastFactorsInfo(input$formula_string, builderTexts)
    })
    output$available_factors_display <- shiny::renderUI({
        currentFactors <- state$factors()
        shiny::p(if (length(currentFactors) == 0) {
            "No factors defined yet (use the 'Factors' tab)."
        } else {
            paste(currentFactors, collapse = ", ")
        })
    })
    output$defined_contrasts_display <- shiny::renderUI({
        contrastData <- state$contrasts()
        if (is.null(contrastData) || nrow(contrastData) == 0) return(shiny::p("No contrasts defined yet."))
        contrastLines <- buildOmicsDesignSaveContrastsTable(contrastData, input$formula_string)$full_format
        shiny::tagList(lapply(contrastLines, function(fullFormat) shiny::p(shiny::tags$code(fullFormat))))
    })
    output$range_preview <- shiny::renderText({
        currentSamplesToTransform <- input$samples_to_transform
        currentRangeStart <- input$range_start
        currentRangeEnd <- input$range_end
        shiny::req(currentSamplesToTransform)
        shiny::req(currentRangeStart, currentRangeEnd)
        formatOmicsDesignRangePreview(currentSamplesToTransform, currentRangeStart, currentRangeEnd)
    })
    output$replicate_inputs <- shiny::renderUI({
        currentSelectedRuns <- input$selected_runs
        shiny::req(currentSelectedRuns)
        replicateInput <- shiny::numericInput(session$ns("replicate_start"),
            paste("Starting replicate number for", length(currentSelectedRuns), builderTexts[["replicate_suffix"]]),
            value = 1, min = 1)
        if (is.na(replicateTestid)) replicateInput else testid(replicateInput, replicateTestid)
    })
    invisible(NULL)
}

#' The design edits: the single and bulk rename, the factor, the metadata assignment, the
#' technical replicates, the contrast and the run removal.
#' @noRd
registerOmicsDesignEdits <- function(input, session, state, adapter) {
    applyRename <- function(renameMap) {
        renamed <- adapter$rename(state$design_matrix(), adapter$read_data(state$data_cln), state$sample_names(),
            renameMap)
        state$design_matrix(renamed$designMatrix)
        state$data_cln(renamed$assayList)
        state$sample_names(renamed$sampleNames)
    }
    shiny::observeEvent(input$rename_sample, adapter$scope(function() {
        shiny::req(input$sample_to_rename, input$new_sample_name)
        applyRename(stats::setNames(input$new_sample_name, input$sample_to_rename))
        shiny::updateTextInput(session, "new_sample_name", value = "")
    }))
    shiny::observeEvent(input$bulk_rename, adapter$scope(function() {
        shiny::req(input$samples_to_transform)
        applyRename(buildOmicsDesignBulkRenameMap(input$samples_to_transform, input$transform_mode,
            input$range_start, input$range_end))
    }))
    shiny::observeEvent(input$add_factor, {
        shiny::req(input$new_factor)
        state$factors(appendOmicsDesignFactorName(state$factors(), input$new_factor))
        shiny::updateTextInput(session, "new_factor", value = "")
    })
    shiny::observeEvent(input$assign_metadata, {
        shiny::req(input$selected_runs, input$factor1_select)
        updatedDesignMatrix <- adapter$assign_metadata(state$design_matrix(), input$selected_runs,
            input$factor1_select, input$factor2_select, input$factor3_select, input$replicate_start)
        state$design_matrix(updatedDesignMatrix)
        state$groups(listOmicsDesignAssignedGroups(updatedDesignMatrix))
    })
    shiny::observeEvent(input$assign_tech_reps, {
        shiny::req(input$tech_rep_samples)
        assignment <- adapter$tech_reps(state$design_matrix(), input$tech_rep_samples,
            input$tech_rep_assignment_mode, input$manual_replicate_number)
        if (isTRUE(assignment$ok)) state$design_matrix(assignment$designMatrix)
        shiny::showNotification(assignment$message, type = assignment$type)
    })
    shiny::observeEvent(input$add_contrast, {
        shiny::req(input$contrast_group1, input$contrast_group2)
        state$contrasts(appendOmicsDesignContrast(state$contrasts(), input$contrast_group1, input$contrast_group2))
    })
    shiny::observeEvent(input$remove_samples, {
        shiny::req(input$samples_to_remove)
        samplesToRemove <- input$samples_to_remove
        state$removed_samples(appendOmicsDesignRemovedSamples(state$removed_samples(), samplesToRemove))
        shiny::updateSelectizeInput(session, "samples_to_remove", selected = "")
        shiny::showNotification(paste("Removed", length(samplesToRemove), "sample(s) from analysis."), type = "message")
    })
    invisible(NULL)
}
