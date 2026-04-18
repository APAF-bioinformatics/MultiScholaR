# Handler 4: Volcano Plot Rendering
da_server_volcano_render_handler <- function(input, output, session, ns, da_data, experiment_paths) {
  # --- TESTTHAT CHECKPOINT CP08 (see test-prot-08-volcano.R) ---
  # Capture in a separate observer so we don't block the UI rendering
  shiny::observeEvent(list(input$volcano_contrast, da_data$da_results_list), {
  shiny::req(input$volcano_contrast, da_data$da_results_list)

  .capture_checkpoint(list(
    da_results_list = da_data$da_results_list,
    selected_contrast = input$volcano_contrast,
    da_q_val_thresh = input$da_q_val_thresh,
    args_row_id = da_data$current_s4_object@protein_id_column
  ), "cp08", "volcano_input")
  }, ignoreInit = TRUE)
  # --- END CP08 ---
  # Render Interactive Volcano Plot (Glimma)
  output$volcano_glimma <- shiny::renderUI({
    shiny::req(input$volcano_contrast, da_data$da_results_list)
    
    # Call the dedicated Glimma generation function from func_prot_da.R
    generateProtDAVolcanoPlotGlimma(
      da_results_list = da_data$da_results_list,
      selected_contrast = input$volcano_contrast,
      da_q_val_thresh = input$da_q_val_thresh,
      args_row_id = da_data$current_s4_object@protein_id_column,
      output_dir = experiment_paths$da_output_dir
    )
  })

  # Render Static Volcano Plot
  output$volcano_plot_static <- shiny::renderPlot({
    shiny::req(input$volcano_contrast, da_data$da_results_list)

    # Call the dedicated static volcano generation function from func_prot_da.R
    generateProtDAVolcanoStatic(
      da_results_list = da_data$da_results_list,
      selected_contrast = input$volcano_contrast,
      da_q_val_thresh = input$da_q_val_thresh,
      lfc_threshold = input$treat_lfc_cutoff,
      show_labels = input$volcano_show_labels,
      n_labels = input$volcano_label_top_n
    )
  })
}

