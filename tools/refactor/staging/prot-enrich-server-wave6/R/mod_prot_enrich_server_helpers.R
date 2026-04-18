setupProtEnrichSelectedTabObserverRegistration <- function(selectedTabFn,
                                                           workflowData,
                                                           enrichmentData,
                                                           session,
                                                           registerSelectedTabObserverFn = registerProtEnrichSelectedTabObserver,
                                                           catFn = cat) {
  registrationState <- list(
    selectedTabProvided = !is.null(selectedTabFn),
    registration = NULL,
    reason = NULL
  )

  if (!is.null(selectedTabFn)) {
    catFn("   mod_prot_enrich_server Step: Setting up tab selection observer\n")
    registrationState$registration <- registerSelectedTabObserverFn(
      selectedTabFn = selectedTabFn,
      workflowData = workflowData,
      enrichmentData = enrichmentData,
      session = session
    )
    registrationState$reason <- "selected_tab_provided"
  } else {
    catFn("   mod_prot_enrich_server Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
    registrationState$reason <- "selected_tab_missing"
  }

  registrationState
}

setupProtEnrichTaxonIdObserverRegistration <- function(workflowData,
                                                       session,
                                                       registerTaxonIdObserverFn = registerProtEnrichTaxonIdObserver) {
  registration <- registerTaxonIdObserverFn(
    workflowData = workflowData,
    session = session
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichMixedSpeciesObserverRegistration <- function(workflowData,
                                                            session,
                                                            registerMixedSpeciesObserverFn = registerProtEnrichMixedSpeciesObserver) {
  registration <- registerMixedSpeciesObserverFn(
    workflowData = workflowData,
    session = session
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

