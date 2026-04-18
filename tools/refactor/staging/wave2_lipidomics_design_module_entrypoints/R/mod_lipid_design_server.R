#' @rdname lipidomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification outputOptions renderText
#' @importFrom logger log_info log_error log_warn
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
#' @importFrom jsonlite write_json read_json
mod_lipid_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL) {
    message(sprintf("--- Entering mod_lipid_design_server ---"))
    message(sprintf("   mod_lipid_design_server Arg: id = %s", id))

    shiny::moduleServer(id, function(input, output, session) {
        message(sprintf("   mod_lipid_design_server Step: Inside moduleServer function"))

        # == Setup shinyFiles =======================================================
        import_bootstrap <- initializeLipidDesignImportBootstrap(
            input = input,
            session = session,
            experimentPaths = experiment_paths,
            volumes = volumes
        )
        resolved_volumes <- import_bootstrap$resolvedVolumes

        # == Modal Logic for Import =================================================

        registerLipidDesignImportModalShell(
            input = input,
            output = output,
            session = session,
            resolvedVolumes = resolved_volumes
        )

        # == Handle Import Confirmation =============================================

        registerLipidDesignImportConfirmationObserver(
            input = input,
            resolvedVolumes = resolved_volumes,
            workflowData = workflow_data,
            experimentPaths = experiment_paths,
            qcTrigger = qc_trigger
        )

        # == Reactivity Checks ======================================================

        output$data_available <- shiny::reactive({
            !is.null(workflow_data$data_tbl) && !is.null(workflow_data$config_list)
        })
        shiny::outputOptions(output, "data_available", suspendWhenHidden = FALSE)

        output$design_matrix_exists <- shiny::reactive({
            !is.null(workflow_data$design_matrix)
        })
        shiny::outputOptions(output, "design_matrix_exists", suspendWhenHidden = FALSE)

        # == Module Integration =====================================================

        builder_results_rv <- registerLipidDesignBuilderModule(
            workflowData = workflow_data
        )

        # == Handle Builder Results =================================================

        registerLipidDesignBuilderResultsObserver(
            builderResultsRv = builder_results_rv,
            workflowData = workflow_data,
            experimentPaths = experiment_paths,
            qcTrigger = qc_trigger
        )

        # == Previews of Saved Data =================================================

        registerLipidDesignPreviewOutputs(
            output = output,
            workflowData = workflow_data
        )
    })
}

