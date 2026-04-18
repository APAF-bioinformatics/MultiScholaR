#' @rdname differentialAbundanceAppletModule
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req showNotification renderPlot renderText withProgress incProgress
#' @importFrom purrr detect map set_names
#' @importFrom dplyr filter select mutate bind_rows slice
#' @importFrom stringr str_detect str_extract
#' @importFrom DT renderDT datatable formatRound
mod_prot_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns # Define namespace function for module

    cat("--- Entering mod_prot_da_server ---\n")
    cat(sprintf("   mod_prot_da_server Arg: id = %s\n", id))
    cat(sprintf("   mod_prot_da_server Arg: workflow_data is NULL = %s\n", is.null(workflow_data)))
    cat(sprintf("   mod_prot_da_server Arg: selected_tab is NULL = %s\n", is.null(selected_tab)))

    cat("=== DIFFERENTIAL ABUNDANCE MODULE SERVER STARTED ===\n")
    cat(sprintf("Module ID: %s\n", id))

    # Initialize reactive values for DA state
    da_data <- shiny::reactiveValues(
      da_results_list = NULL,
      contrasts_available = NULL,
      analysis_complete = FALSE,
      current_s4_object = NULL,
      formula_from_s4 = NULL,
      current_row_clusters = NULL,
      current_col_clusters = NULL
    )

    cat("   mod_prot_da_server Step: Reactive values initialized\n")

    # --- Call Top-Level Server Handlers ---

    # Handler 1: Initialization (Tab selection, state updates, contrast detection)
    da_server_init_handlers(input, output, session, da_data, workflow_data, selected_tab)

    # Handler 2: Session Loading
    da_server_load_session_handler(input, output, session, da_data, workflow_data, experiment_paths)

    # Handler 3: DA Analysis Execution
    da_server_run_analysis_handler(input, output, session, ns, da_data, workflow_data, experiment_paths)

    # Handler 4: Volcano Plot Rendering (Restored)
    da_server_volcano_render_handler(input, output, session, ns, da_data, experiment_paths)

    # Handler 5: Heatmap Rendering
    da_server_heatmap_render_handler(input, output, session, ns, da_data, experiment_paths)

    # Handler 6: Results Table Rendering
    da_server_table_render_handler(input, output, session, da_data)

    # Return DE data for potential use by parent module
    return(da_data)
  })
}

