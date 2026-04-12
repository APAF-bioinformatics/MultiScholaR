#' @rdname normalizationAppletModule 
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req showNotification renderPlot renderText withProgress incProgress observe updateSelectInput
#' @importFrom readr write_tsv
#' @importFrom grid grid.draw grid.newpage pushViewport popViewport viewport grid.layout rasterGrob
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom vroom vroom_write
#' @importFrom dplyr filter select mutate distinct
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
mod_prot_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    message("=== NORMALIZATION MODULE SERVER STARTED ===")
    message(sprintf("Module ID: %s", id))
    message(sprintf("workflow_data is NULL: %s", is.null(workflow_data)))
    if (!is.null(workflow_data$state_manager)) {
      message(sprintf("Current state at module start: %s", workflow_data$state_manager$current_state))
    }
    norm_data <- createProtNormReactiveState()

    getPlotAesthetics <- function() {
      getProtNormPlotAesthetics(input$color_variable, input$shape_variable)
    }

    getRuvGroupingVariable <- function() {
      getProtNormRuvGroupingVariable(input$ruv_grouping_variable)
    }

    generatePreNormalizationQc <- function() {
      norm_data$qc_plot_paths <- generateProtNormPreNormalizationQcArtifacts(
        stateManager = workflow_data$state_manager,
        qcDir = experiment_paths$protein_qc_dir,
        aesthetics = getPlotAesthetics(),
        qcPlotPaths = norm_data$qc_plot_paths
      )
    }

    generatePostNormalizationQc <- function(normalized_s4) {
      norm_data$qc_plot_paths <- generateProtNormPostNormalizationQcArtifacts(
        normalizedS4 = normalized_s4,
        qcDir = experiment_paths$protein_qc_dir,
        aesthetics = getPlotAesthetics(),
        qcPlotPaths = norm_data$qc_plot_paths
      )

      norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1
    }

    generateRuvCorrectedQc <- function(ruv_corrected_s4) {
      norm_data$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(
        ruvCorrectedS4 = ruv_corrected_s4,
        qcDir = experiment_paths$protein_qc_dir,
        aesthetics = getPlotAesthetics(),
        qcPlotPaths = norm_data$qc_plot_paths
      )

      norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1
    }

    registerProtNormServerObservers(
      input = input,
      output = output,
      session = session,
      selectedTab = selected_tab,
      workflowData = workflow_data,
      normData = norm_data,
      experimentPaths = experiment_paths,
      omicType = omic_type,
      experimentLabel = experiment_label,
      generatePreNormalizationQcFn = generatePreNormalizationQc,
      generatePostNormalizationQcFn = generatePostNormalizationQc,
      generateRuvCorrectedQcFn = generateRuvCorrectedQc,
      getPlotAestheticsFn = getPlotAesthetics,
      getRuvGroupingVariableFn = getRuvGroupingVariable,
      checkMemoryUsageFn = checkProtNormMemoryUsage
    )

    registerProtNormQcImageOutputs(
      output = output,
      normData = norm_data,
      proteinQcDir = experiment_paths$protein_qc_dir
    )
    registerProtNormRenderOutputs(
      output = output,
      normData = norm_data,
      ruvMode = input$ruv_mode,
      groupingVariable = getRuvGroupingVariable()
    )

    return(norm_data)
  })
}
