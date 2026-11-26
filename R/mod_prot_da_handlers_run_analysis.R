# The notice for the contrasts whose q-values fell back to Benjamini-Hochberg (KIT-268, VAL-312): it names them (GUI audit
# A112-03 had shown one sentence without names). `manifest$results` is named by contrast; `manifest$contrasts` maps each
# to its friendly name.
warnDaQvalueFallback <- function(manifest) {
  fell_back <- names(Filter(function(result) length(result$qvalue_warnings) > 0L, manifest$results))
  if (length(fell_back)) {
    notifyDaQvalueFallback(daQvalueFallbackNames(fell_back, manifest$contrasts), NA_character_)
  }
  invisible(NULL)
}

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
  # GUI audit A112-03 (F112-01): one notice for both paths (warnDaQvalueFallback), from the manifest's contrast records:
  # in process when the run returns, on a managed session when the publication step delivers the manifest.
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
      warnDaQvalueFallback(result$manifest)
      NULL
    },
    read_result = readManagedDaResult)
  # KIT-112 D3: the proteomics omics runner the contrasts dispatch through on a managed session: the one the
  # normalisation tab installed (KIT-111, one per session), else this tab's own (the plugin mounted alone). An
  # unmanaged session has none, and the DA runs in process as before.
  runner <- if (inherits(session, c("ShinySession", "session_proxy", "MockShinySession")))
    session$userData$workflow_omics_actions[["proteomics"]] %||% installWorkflowOmicsActions(session, workflow_data,
      "proteomics", experiment_label %||% "dia_baseline_smoke", experiment_paths)
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

    report_error <- function(e) {
      cat(paste("*** ERROR in DA analysis:", e$message, "\n"))
      shiny::showNotification(
        paste("Error in DA analysis:", e$message),
        type = "error",
        duration = 10
      )
    }
    tryCatch(
      {
        contrasts <- workflow_data$contrasts_tbl
        annotations <- workflow_data$uniprot_dat_cln
        settings <- list(formula_string = input$formula_string, da_q_val_thresh = input$da_q_val_thresh,
          treat_lfc_cutoff = input$treat_lfc_cutoff)
        managed <- !is.null(workflowSessionApplication(session))
        publication <- function() managed_step$apply(list(annotations =
          if (is.data.frame(annotations) && nrow(annotations))
            saveDIADisplayTable(workflow_data$state_manager, annotations) else NULL))
        if (managed && !is.null(runner)) {
          # KIT-112 F112-02: before the request half, which may commit the DA arguments (moving the head), a click while
          # the contrasts or the publication run is refused, and a publication waiting for its retry is retried without
          # recomputing the contrasts.
          running <- if (!is.null(shiny::isolate(managed_step$active()))) "DA publication" else
            if (!is.null(shiny::isolate(runner$step$active()))) "proteomics step"
          if (!is.null(running)) {
            shiny::showNotification(paste("A", running, "action is already running"), type = "warning")
            shiny::removeNotification("da_working")
            return(invisible(NULL))
          }
          # GUI audit A112-05 / P-39 (a): the runner is shared by the proteomics tabs. While it holds another tab's
          # publication for a retry, a click here commits nothing (the runner would resume that tab's retry on a moved head).
          held <- shiny::isolate(runner$step$pending())
          if (!is.null(held) && !identical(held$reservation$operation, proteomicsDiaBinding$da$operation)) {
            shiny::showNotification(paste0("An earlier proteomics step", if (is.character(held$reservation$operation))
              paste0(" (", held$reservation$operation, ")"), " could not write its files. Press that step's button again, ",
              "then run the differential abundance analysis."), type = "warning", duration = 10)
            shiny::removeNotification("da_working")
            return(invisible(NULL))
          }
          pending <- shiny::isolate(managed_step$pending())
          if (!is.null(held) || !is.null(pending)) {
            # (b): the retry keeps the earlier run's request, and says so. (c): the publication is retried with its own
            # parameters, so annotations changed since the failure do not refuse it; a pending contrast (F112-06) is the
            # runner's retry, finished by the first request's callbacks.
            shiny::showNotification(paste("Retrying the earlier differential abundance run, whose files could not be written,",
              "with that run's settings. Run again afterwards to apply changed settings."), type = "warning", duration = 10)
            if (is.null(pending)) return(runner$run(list()))
            shiny::removeNotification("da_working")
            managed_step$apply(pending$parameters)
            return(invisible(NULL))
          }
          # KIT-112 D3: the request half here, one DA-contrast-output per contrast on the worker (its source named, its
          # route in the request), then the accept half and the publication.
          manager <- workflow_data$state_manager
          entry <- proteomicsDiaBinding$da
          request <- entry$request(manager, settings, da_data, workflow_data, contrasts, annotations)
          route <- diaActiveRoute(manager)
          results <- list()
          deliver <- function(name) function(result) {
            results[[name]] <<- result$value
            # R1: the controller holds the base fit an in-process run would, for the manifest's model.
            if (!is.null(result$display$base_fit)) manager$.context$da_base_fit_reference <- result$display$base_fit
          }
          runner$run(lapply(request$steps, function(step) protNormManagedAction(entry$operation, c(step, route = route),
            request$source, NULL, deliver(step$name))),
            # The runner calls `done` outside the owner's scope.
            done = function() withDIAStateScope(manager, function() {
              entry$accept(manager, request$manifest, results)
              shiny::removeNotification("da_working")
              # KIT-243 (S2 (b)): the earlier files under this DA's names removed at the accept, so a worker-side
              # publication failure (which never reaches `publishWorkflowDaResult`) leaves none; the marker refuses Copy.
              removed <- tryCatch({
                removeWorkflowDaOwnFigures(manager, experiment_paths, omic_type)
                TRUE
              }, error = function(e) {
                report_error(e)
                FALSE
              })
              if (removed) publication()
            }),
            failed = function(e) {
              report_error(e)
              shiny::removeNotification("da_working")
            })
          return(invisible(NULL))
        }
        manifest <- runDIAReferenceDaAnalysis(settings, da_data, workflow_data, experiment_paths,
          contrasts, annotations, publish = !managed)
        if (managed) {
          shiny::removeNotification("da_working")
          publication()
          return(invisible(NULL))
        }
        update_contrast_choices(manifest)
        warnDaQvalueFallback(manifest)
        shiny::showNotification("Differential abundance analysis completed successfully!", type = "message", duration = 5)
        # KIT-111 F106-02: the in-process run returns here, so it removes its working notice first.
        shiny::removeNotification("da_working")
        return(manifest)
      },
      error = report_error
    )

    shiny::removeNotification("da_working")
    })
  })
}
