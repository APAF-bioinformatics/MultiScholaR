#' @rdname enrichmentAnalysisAppletModule
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req showNotification
#' @importFrom purrr detect
#' @importFrom readr write_tsv
#' @importFrom dplyr filter select mutate case_when slice_head all_of pull
#' @importFrom stringr str_split str_detect str_extract
mod_prot_enrich_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    cat("--- Entering mod_prot_enrich_server ---\n")
    cat(sprintf("   mod_prot_enrich_server Arg: id = %s\n", id))
    cat(sprintf("   mod_prot_enrich_server Arg: workflow_data is NULL = %s\n", is.null(workflow_data)))
    cat(sprintf("   mod_prot_enrich_server Arg: selected_tab is NULL = %s\n", is.null(selected_tab)))

    cat("=== ENRICHMENT ANALYSIS MODULE SERVER STARTED ===\n")
    cat(sprintf("Module ID: %s\n", id))

    enrichment_data_setup <- setupProtEnrichReactiveValues()
    enrichment_data <- enrichment_data_setup$enrichmentData

    analysis_method_bootstrap <- setupProtEnrichAnalysisMethodBootstrap(
      input = input,
      enrichmentData = enrichment_data
    )
    supported_organisms <- analysis_method_bootstrap$supportedOrganisms
    current_analysis_method <- analysis_method_bootstrap$currentAnalysisMethod

    setupProtEnrichObserverRegistrationBootstrap(
      selectedTabFn = selected_tab,
      input = input,
      workflowData = workflow_data,
      enrichmentData = enrichment_data,
      session = session,
      setupTaxonIdObserverRegistrationFn = setupProtEnrichTaxonIdObserverRegistration,
      setupMixedSpeciesObserverRegistrationFn = setupProtEnrichMixedSpeciesObserverRegistration,
      setupSelectedContrastObserverRegistrationFn = setupProtEnrichSelectedContrastObserverRegistration,
      setupSelectedTabObserverRegistrationFn = setupProtEnrichSelectedTabObserverRegistration,
      setupDaResultsObserverRegistrationFn = setupProtEnrichDaResultsObserverRegistration
    )

    setupProtEnrichDisplayStatusOutputBootstrap(
      output = output,
      input = input,
      enrichmentData = enrichment_data,
      currentAnalysisMethodFn = current_analysis_method
    )

    setupProtEnrichRunOutputDownloadBootstrap(
      output = output,
      input = input,
      enrichmentData = enrichment_data,
      workflowData = workflow_data,
      experimentPaths = experiment_paths,
      currentAnalysisMethodFn = current_analysis_method,
      runObserverPreflightFn = runProtEnrichObserverPreflight
    )

    # Return enrichment data for potential use by parent module
    return(enrichment_data)
  })
}

