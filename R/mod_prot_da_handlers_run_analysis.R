# Handler 3: DA Analysis Execution
da_server_run_analysis_handler <- function(input, output, session, ns, da_data, workflow_data, experiment_paths,
                                           omic_type = "proteomics", experiment_label = NULL) {
  # KIT-019: under a managed application the DA publication (exports, figures, Glimma pages) runs in the worker
  # as the registered `da_publication` operation; the limma kernel runs in process before the request.
  update_contrast_choices <- function(manifest) {
    choices <- stats::setNames(manifest$contrasts$friendly_names, manifest$contrasts$friendly_names)
    for (id in c("volcano_contrast", "heatmap_contrast", "table_contrast"))
      shiny::updateSelectInput(session, id, choices = choices)
  }
  managed_step <- installWorkflowManagedStep(session, input, workflow_data, "da_publication",
    publication = list(experiment_label = experiment_label %||% "dia_baseline_smoke", experiment_paths = experiment_paths,
      omic_type = omic_type),
    # KIT-037: the step is acknowledged on the status text it shows at completion; the static volcano renders
    # only when its sub-tab is shown, so its identity is an observation (the alt text), not the acknowledgement.
    ui = list(text = "da_status", plot = NULL, ack = "workflow_ack", notification = "da_publication",
      label = "DA publication", step_label = "14_da_published",
      plot_identity = function(id) da_data$volcano_plot_identity <- id),
    deliver = function(result) {
      completeDIADaDisplay(da_data, workflow_data, result$manifest)
      update_contrast_choices(result$manifest)
      NULL
    },
    read_result = readManagedDaResult)
  shiny::observeEvent(input$run_da_analysis, {
    withDIAStateScope(getDIAActiveManager(), function() {
    cat("=== STARTING DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

    if (is.null(da_data$current_s4_object)) {
      shiny::showNotification(
        "Error in DA analysis: no filtered protein data is loaded.",
        type = "error",
        duration = 10
      )
      return(invisible(NULL))
    }

    shiny::showNotification("Running differential abundance analysis...", id = "da_working", duration = NULL)

    tryCatch(
      {
        contrasts <- workflow_data$contrasts_tbl
        annotations <- workflow_data$uniprot_dat_cln
        settings <- list(formula_string = input$formula_string, da_q_val_thresh = input$da_q_val_thresh,
          treat_lfc_cutoff = input$treat_lfc_cutoff)
        managed <- !is.null(workflowSessionApplication(session))
        manifest <- runDIAReferenceDaAnalysis(settings, da_data, workflow_data, experiment_paths,
          contrasts, annotations, publish = !managed)
        if (managed) {
          shiny::removeNotification("da_working")
          managed_step$apply(list(annotations = if (is.data.frame(annotations) && nrow(annotations))
            saveDIADisplayTable(workflow_data$state_manager, annotations) else NULL))
          return(invisible(NULL))
        }
        update_contrast_choices(manifest)
        warnings <- unlist(lapply(manifest$results, function(result) result$qvalue_warnings), recursive = FALSE)
        if (length(warnings)) shiny::showNotification(
          "Some contrasts used the existing BH fallback after q-value estimation failed.", type = "warning", duration = NULL)
        shiny::showNotification("Differential abundance analysis completed successfully!", type = "message", duration = 5)
        return(manifest)
      },
      error = function(e) {
        cat(paste("*** ERROR in DA analysis:", e$message, "\n"))
        shiny::showNotification(
          paste("Error in DA analysis:", e$message),
          type = "error",
          duration = 10
        )
      }
    )

    shiny::removeNotification("da_working")
    })
  })
}
