runProtEnrichObserverShell <- function(selectedContrast,
                                       runAnalysisBody,
                                       withProgressFn = shiny::withProgress,
                                       finalizeFailureFn = finalizeProtEnrichObserverFailure,
                                       finalizeObserverRunFn = finalizeProtEnrichObserverRun,
                                       progressMessage = "Running enrichment analysis...",
                                       progressValue = 0) {
  analysisCompleted <- FALSE
  failureResult <- NULL
  successResult <- NULL

  tryCatch({
    withProgressFn(message = progressMessage, value = progressValue, {
      runAnalysisBody()
    })
    analysisCompleted <- TRUE
  }, error = function(e) {
    failureResult <<- finalizeFailureFn(errorMessage = e$message)
  })

  if (isTRUE(analysisCompleted)) {
    successResult <- finalizeObserverRunFn(
      completed = TRUE,
      selectedContrast = selectedContrast
    )
  }

  list(
    completed = analysisCompleted,
    selectedContrast = selectedContrast,
    failureResult = failureResult,
    successResult = successResult
  )
}

finalizeProtEnrichObserverFailure <- function(errorMessage,
                                              reportAnalysisErrorFn = reportProtEnrichAnalysisError,
                                              removeWorkingNotificationFn = removeProtEnrichWorkingNotification) {
  reportResult <- reportAnalysisErrorFn(errorMessage = errorMessage)
  cleanupResult <- removeWorkingNotificationFn()

  list(
    errorMessage = errorMessage,
    reportResult = reportResult,
    cleanupResult = cleanupResult
  )
}

reportProtEnrichCompletion <- function(selectedContrast,
                                       notifyCompletionFn = notifyProtEnrichCompletion,
                                       logCompletionFn = logProtEnrichCompletion) {
  notificationResult <- notifyCompletionFn(selectedContrast = selectedContrast)
  logResult <- logCompletionFn()

  list(
    selectedContrast = selectedContrast,
    notificationResult = notificationResult,
    logResult = logResult
  )
}

finalizeProtEnrichObserverRun <- function(completed,
                                          selectedContrast = NULL,
                                          reportCompletionFn = reportProtEnrichCompletion,
                                          removeWorkingNotificationFn = removeProtEnrichWorkingNotification) {
  reportResult <- NULL

  if (isTRUE(completed)) {
    reportResult <- reportCompletionFn(selectedContrast = selectedContrast)
  }

  cleanupResult <- removeWorkingNotificationFn()

  list(
    completed = isTRUE(completed),
    selectedContrast = selectedContrast,
    reportResult = reportResult,
    cleanupResult = cleanupResult
  )
}

logProtEnrichCompletion <- function(message = "=== ENRICHMENT ANALYSIS COMPLETED ===\n",
                                    catFn = cat) {
  catFn(message)

  list(message = message)
}

removeProtEnrichWorkingNotification <- function(notificationId = "enrichment_working",
                                                removeNotificationFn = shiny::removeNotification) {
  removeNotificationFn(notificationId)

  list(notificationId = notificationId)
}

handoffProtEnrichObserverRun <- function(input,
                                         enrichmentData,
                                         workflowData,
                                         experimentPaths,
                                         currentAnalysisMethodFn,
                                         showNotificationFn = shiny::showNotification,
                                         runObserverShellFn = runProtEnrichObserverShell,
                                         runAnalysisBodyFn = runProtEnrichAnalysisBody) {
  selectedContrast <- input$selected_contrast
  notificationMessage <- "Running enrichment analysis..."
  notificationId <- "enrichment_working"
  notificationDuration <- NULL

  showNotificationFn(
    notificationMessage,
    id = notificationId,
    duration = notificationDuration
  )

  shellResult <- runObserverShellFn(
    selectedContrast = selectedContrast,
    runAnalysisBody = function() {
      runAnalysisBodyFn(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
    }
  )

  list(
    selectedContrast = selectedContrast,
    notificationMessage = notificationMessage,
    notificationId = notificationId,
    notificationDuration = notificationDuration,
    shellResult = shellResult
  )
}

runProtEnrichObserverPreflight <- function(input,
                                           enrichmentData,
                                           workflowData,
                                           experimentPaths,
                                           currentAnalysisMethodFn,
                                           reqFn = shiny::req,
                                           handoffObserverRunFn = handoffProtEnrichObserverRun) {
  selectedContrast <- input$selected_contrast
  daResultsData <- enrichmentData$da_results_data

  reqFn(selectedContrast, daResultsData)

  handoffResult <- handoffObserverRunFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    selectedContrast = selectedContrast,
    daResultsData = daResultsData,
    handoffResult = handoffResult
  )
}

buildProtEnrichSupportedOrganisms <- function() {
  tibble::tribble(
    ~taxid,     ~id,             ~name,
    "9606",     "hsapiens",      "Homo sapiens",
    "10090",    "mmusculus",     "Mus musculus",
    "10116",    "rnorvegicus",   "Rattus norvegicus",
    "7227",     "dmelanogaster", "Drosophila melanogaster",
    "6239",     "celegans",      "Caenorhabditis elegans",
    "4932",     "scerevisiae",   "Saccharomyces cerevisiae",
    "3702",     "athaliana",     "Arabidopsis thaliana",
    "7955",     "drerio",        "Danio rerio",
    "9031",     "ggallus",       "Gallus gallus",
    "9823",     "sscrofa",       "Sus scrofa",
    "9913",     "btaurus",       "Bos taurus",
    "9544",     "mmulatta",      "Macaca mulatta",
    "9598",     "ptroglodytes",  "Pan troglodytes"
  )
}

createProtEnrichSupportedOrganismsReactive <- function(reactiveFn = shiny::reactive,
                                                       buildSupportedOrganismsFn = buildProtEnrichSupportedOrganisms) {
  reactiveFn({
    buildSupportedOrganismsFn()
  })
}

createProtEnrichReactiveValues <- function(reactiveValuesFn = shiny::reactiveValues) {
  reactiveValuesFn(
    enrichment_results = NULL,
    contrasts_available = NULL,
    analysis_complete = FALSE,
    current_s4_object = NULL,
    da_results_data = NULL,
    gprofiler_results = NULL,
    clusterprofiler_results = NULL,
    stringdb_results = NULL,
    analysis_method = NULL,
    organism_supported = NULL,
    all_enrichment_results = list(),
    current_contrast_results = list(),
    enrichment_plots = list()
  )
}

normaliseProtEnrichTaxonId <- function(organismTaxid) {
  if (is.null(organismTaxid) || length(organismTaxid) != 1L || is.na(organismTaxid)) {
    stop("organism_taxid must be a single positive integer NCBI taxonomy ID.", call. = FALSE)
  }

  taxonId <- trimws(as.character(organismTaxid))
  if (!nzchar(taxonId) || !grepl("^[0-9]+$", taxonId)) {
    stop("organism_taxid must be a single positive integer NCBI taxonomy ID.", call. = FALSE)
  }

  taxonId <- sub("^0+", "", taxonId)
  if (!nzchar(taxonId) || identical(taxonId, "0")) {
    stop("organism_taxid must be a single positive integer NCBI taxonomy ID.", call. = FALSE)
  }

  taxonId
}

validateProtEnrichProcessParameters <- function(organismTaxid,
                                                upCutoff,
                                                downCutoff,
                                                qCutoff,
                                                correctionMethod = NULL) {
  taxonId <- normaliseProtEnrichTaxonId(organismTaxid)
  numericParams <- list(
    up_cutoff = upCutoff,
    down_cutoff = downCutoff,
    q_cutoff = qCutoff
  )
  numericParams <- lapply(names(numericParams), function(name) {
    value <- numericParams[[name]]
    if (is.null(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("%s must be a single finite numeric value.", name), call. = FALSE)
    }
    value <- suppressWarnings(as.numeric(value))
    if (!is.finite(value)) {
      stop(sprintf("%s must be a single finite numeric value.", name), call. = FALSE)
    }
    value
  }) |>
    stats::setNames(names(numericParams))

  if (numericParams$q_cutoff < 0 || numericParams$q_cutoff > 1) {
    stop("q_cutoff must be between 0 and 1.", call. = FALSE)
  }

  if (!is.null(correctionMethod) &&
      (length(correctionMethod) != 1L || is.na(correctionMethod) || !nzchar(trimws(as.character(correctionMethod))))) {
    stop("correction_method must be a single non-empty string.", call. = FALSE)
  }

  list(
    organism_taxid = taxonId,
    up_cutoff = numericParams$up_cutoff,
    down_cutoff = numericParams$down_cutoff,
    q_cutoff = numericParams$q_cutoff,
    correction_method = if (is.null(correctionMethod)) NULL else trimws(as.character(correctionMethod))
  )
}
