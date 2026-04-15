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
#'
#' @rdname lipidomicsDesignMatrixBuilderModule
#' @export
mod_lipid_design_builder_server <- function(id, data_tbl, config_list, column_mapping
                                            , existing_design_matrix = NULL
                                            , existing_contrasts = NULL) {
    shiny::moduleServer(id, function(input, output, session) {

        # Reactive value to store the final results
        result_rv <- shiny::reactiveVal(NULL)

        # == Initial State Setup =================================================
        initial_state <- shiny::reactive({
            shiny::req(data_tbl())
            shiny::req(config_list())

            buildLipidDesignInitialState(
                assayList = data_tbl()
                , configList = config_list()
                , colMap = column_mapping()
                , existingDesignMatrix = existing_design_matrix
                , existingContrasts = existing_contrasts
            )
        })

        # Reactive values for the current, mutable state
        design_matrix <- shiny::reactiveVal()
        data_cln_reactive <- shiny::reactiveVal()
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
            , resultRv = result_rv
        )

        # Return the reactive containing results
        return(result_rv)
    })
}

