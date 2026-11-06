#' Lipidomics Design Matrix Builder Server Module
#'
#' @param id Module ID.
#' @param data_tbl A reactive expression that returns the data table LIST.
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#' @param existing_design_matrix Optional reactive expression that returns an existing
#'   design matrix (e.g., from "Import Existing Design"). If provided and has valid
#'   group assignments, will be used instead of creating a fresh matrix.
#' @param existing_contrasts Optional reactive expression that returns an existing
#'   contrasts table (e.g., from "Import Existing Design").
#' @param state_owner Optional function returning the session's state owner; the builder
#'   then holds its tables by reference. Without it they stay plain values.
#' @param profile An omic profile (KIT-096): its run order for a fresh design, and the level order and removed runs on
#'   the saved result. `NULL`, the default, is lipidomics.
#'
#' @rdname lipidomicsDesignMatrixBuilderModule
#' @export
mod_lipid_design_builder_server <- function(id, data_tbl, config_list, column_mapping
                                            , existing_design_matrix = NULL
                                            , existing_contrasts = NULL
                                            , state_owner = NULL
                                            , profile = NULL) {
    shiny::moduleServer(id, function(input, output, session) {

        # Reactive value to store the final results
        result_rv <- shiny::reactiveVal(NULL)

        # KIT-075: the builder's table lists are held by reference in the owner (plain without one).
        currentOwner <- function() if (is.function(state_owner)) state_owner()

        # == Initial State Setup =================================================
        initial_state <- shiny::reactive({
            shiny::req(data_tbl())
            shiny::req(config_list())

            state <- buildLipidDesignInitialState(
                assayList = data_tbl()
                , configList = config_list()
                , colMap = column_mapping()
                , existingDesignMatrix = existing_design_matrix
                , existingContrasts = existing_contrasts
                , mixedsortFn = profile$run_order %||% gtools::mixedsort
            )
            # KIT-075: a reactive keeps its value; it keeps the tables' reference.
            state$data_cln <- saveWorkflowOwnedValue(state$data_cln, "design_builder_initial_data_cln", currentOwner())
            state
        })

        # Reactive values for the current, mutable state
        design_matrix <- shiny::reactiveVal()
        data_cln_reactive <- omicsBuilderTableVal("design_builder_data_cln", currentOwner)
        sample_names_reactive <- shiny::reactiveVal()
        groups <- shiny::reactiveVal()
        factors <- shiny::reactiveVal()
        contrasts <- shiny::reactiveVal()
        removed_samples <- shiny::reactiveVal(character(0))

        registerLipidDesignInitialStateShells(
            dataTbl = data_tbl
            , input = input
            , initialState = initial_state
            , session = session
            , designMatrix = design_matrix
            , dataClnReactive = data_cln_reactive
            , sampleNamesReactive = sample_names_reactive
            , groups = groups
            , factors = factors
            , contrasts = contrasts
            , removedSamples = removed_samples
        )


        # Create a proxy for the main data table
        proxy_data_table <- DT::dataTableProxy("data_table")

        # == UI Rendering and Updates =============================================

        # Update sample selection inputs when names change or samples are removed
        registerLipidDesignSampleSelectionInputShells(
            input = input
            , session = session
            , designMatrix = design_matrix
            , removedSamples = removed_samples
        )

        registerLipidDesignFactorGroupDropdownShells(
            input = input
            , session = session
            , factors = factors
            , groups = groups
        )

        registerLipidDesignDataTableShells(
            output = output
            , proxyDataTable = proxy_data_table
            , designMatrix = design_matrix
            , removedSamples = removed_samples
        )

        registerLipidDesignSummaryOutputShells(
            output = output
            , designMatrix = design_matrix
            , removedSamples = removed_samples
            , formulaString = function() input$formula_string
        )

        registerLipidDesignAdjacentOutputShells(
            output = output
            , session = session
            , factors = factors
            , contrasts = contrasts
            , formulaString = function() input$formula_string
            , samplesToTransform = function() input$samples_to_transform
            , rangeStart = function() input$range_start
            , rangeEnd = function() input$range_end
            , selectedRuns = function() input$selected_runs
        )

        registerLipidDesignSampleRenameShells(
            input = input
            , session = session
            , designMatrix = design_matrix
            , dataClnReactive = data_cln_reactive
            , sampleNamesReactive = sample_names_reactive
        )

        registerLipidDesignMetadataAssignmentShells(
            input = input
            , session = session
            , designMatrix = design_matrix
            , factors = factors
            , groups = groups
        )

        registerLipidDesignTechRepAssignmentShells(
            input = input
            , designMatrix = design_matrix
        )

        registerLipidDesignContrastManagementShells(
            input = input
            , contrasts = contrasts
        )

        registerLipidDesignSampleRemovalShells(
            input = input
            , removedSamples = removed_samples
            , session = session
        )

        registerLipidDesignResetRequestShells(
            input = input
            , session = session
        )

        registerLipidDesignResetConfirmationShells(
            input = input
            , initialState = initial_state
            , designMatrix = design_matrix
            , dataClnReactive = data_cln_reactive
            , sampleNamesReactive = sample_names_reactive
            , removedSamples = removed_samples
            , factors = factors
            , groups = groups
            , contrasts = contrasts
            , session = session
        )

        registerLipidDesignSaveResultsShells(
            input = input
            , designMatrix = design_matrix
            , removedSamples = removed_samples
            , dataClnReactive = data_cln_reactive
            , columnMapping = column_mapping
            , contrasts = contrasts
            , configList = config_list
            , resultRv = function(value) {
                # KIT-096: under a profile, the level order the user entered and the runs the user removed.
                if (!is.null(profile)) {
                    value$group_levels <- profile$group_levels(value$design_matrix, factors())
                    value$excluded_runs <- removed_samples()
                }
                # KIT-075: the returned result carries the tables' reference.
                value$data_cln <- saveWorkflowOwnedValue(value$data_cln, "design_builder_result_data_cln", currentOwner())
                result_rv(value)
            }
        )

        # Return the reactive containing results
        return(result_rv)
    })
}

