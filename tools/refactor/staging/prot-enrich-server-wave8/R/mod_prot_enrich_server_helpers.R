registerProtEnrichAnalysisMethodDisplayOutput <- function(output,
                                                          currentAnalysisMethodFn,
                                                          renderTextFn = shiny::renderText,
                                                          formatAnalysisMethodTextFn = formatProtEnrichAnalysisMethodText) {
  analysisMethodDisplay <- renderTextFn({
    formatAnalysisMethodTextFn(currentAnalysisMethodFn())
  })

  output$analysis_method_display <- analysisMethodDisplay
  analysisMethodDisplay
}

registerProtEnrichTaxonIdObserver <- function(workflowData,
                                              session,
                                              observeEventFn = shiny::observeEvent,
                                              updateTextInputFn = shiny::updateTextInput) {
  observeEventFn(
    eventExpr = workflowData$taxon_id,
    handlerExpr = {
    taxonId <- workflowData$taxon_id
    observerState <- list(
      taxonId = taxonId,
      updated = FALSE,
      reason = NULL
    )

    if (!is.null(taxonId)) {
      updateTextInputFn(session, "organism_taxid", value = taxonId)
      observerState$updated <- TRUE
      observerState$reason <- "updated"
    } else {
      observerState$reason <- "missing_taxon_id"
    }

    observerState
  },
    ignoreInit = TRUE,
    ignoreNULL = TRUE
  )
}

registerProtEnrichMixedSpeciesObserver <- function(workflowData,
                                                   session,
                                                   observeEventFn = shiny::observeEvent,
                                                   updateCheckboxInputFn = shiny::updateCheckboxInput,
                                                   showNotificationFn = shiny::showNotification,
                                                   catFn = cat) {
  observeEventFn(
    eventExpr = workflowData$mixed_species_analysis,
    handlerExpr = {
    mixedSpeciesAnalysis <- workflowData$mixed_species_analysis
    observerState <- list(
      mixedSpeciesAnalysis = mixedSpeciesAnalysis,
      filterEnabled = FALSE,
      notificationShown = FALSE,
      reason = NULL
    )

    if (!is.null(mixedSpeciesAnalysis) &&
        isTRUE(mixedSpeciesAnalysis$enabled)) {
      updateCheckboxInputFn(session, "enable_organism_filter", value = TRUE)
      catFn("*** ENRICHMENT: Auto-enabled organism filter (mixed species FASTA detected at import) ***\n")

      showNotificationFn(
        sprintf(
          "Multi-species FASTA detected. Filtering to %s enabled.",
          mixedSpeciesAnalysis$selected_organism
        ),
        type = "message",
        duration = 5
      )

      observerState$filterEnabled <- TRUE
      observerState$notificationShown <- TRUE
      observerState$reason <- "enabled_filter"
    } else {
      observerState$reason <- "mixed_species_disabled"
    }

    observerState
  },
    ignoreInit = TRUE,
    ignoreNULL = TRUE
  )
}

registerProtEnrichContrastsDisplayOutput <- function(output,
                                                     enrichmentData,
                                                     renderTextFn = shiny::renderText,
                                                     formatContrastsTextFn = formatProtEnrichContrastsText) {
  contrastsDisplay <- renderTextFn({
    formatContrastsTextFn(enrichmentData$contrasts_available)
  })

  output$contrasts_display <- contrastsDisplay
  contrastsDisplay
}

