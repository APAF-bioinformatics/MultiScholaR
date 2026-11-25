# KIT-250: the proteomics import's organism dialog ('Mixed species FASTA'). With the option ticked the import
# completes only when an organism is confirmed; Cancel clears the import as a failed import is cleared, so no design
# can be saved until an organism is confirmed or the option is unticked. Ids E1a to E7 and G1 to G3 are the ticket's
# (section 10). Calls drop arguments a function does not have (`k250Call`), so that on the code before the fix each
# E test fails by assertion at the value the ticket states, not by an "unused argument" error.

k250W1 <- paste(
  "Import not finished: no organism was chosen. Press 'Process Imported Data' again and choose the organism in the",
  "window that follows, or untick 'Mixed species FASTA' and enter the organism yourself."
)
k250W2 <- "Choose the organism in the window to finish the import."

k250Call <- function(fun, args) do.call(fun, args[names(args) %in% names(formals(fun))])

k250ImportResult <- function() list(data = data.frame(Protein.Group = c("P1", "P2")), import_summary = list(n_rows = 2L))

k250Table <- function() {
  data.frame(Protein.Group = c("P1;P3", "P3", "P2-2"), value = c(1, 2, 3), stringsAsFactors = FALSE)
}

k250WorkflowData <- function(status = list(setup_import = "pending", design_matrix = "disabled",
                                           quality_control = "disabled")) {
  list2env(list(
    data_tbl = k250Table(), data_cln = k250Table(), data_format = "diann", data_type = "peptide",
    column_mapping = list(protein_col = "Protein.Group"), fasta_file_path = "proteins.fasta",
    aa_seq_tbl_final = data.frame(x = 1), config_list = list(a = 1),
    processing_log = list(setup_import = list(taxon_id = 10090L, organism = "Mus musculus", n_proteins = 3L)),
    mixed_species_analysis = NULL, taxon_id = 10090L, organism_name = "Mus musculus", tab_status = status
  ))
}

k250LocalData <- function(waiting = TRUE) {
  list2env(list(
    processing = FALSE,
    organism_distribution = data.frame(organism_name = c("Homo sapiens", "Mus musculus"), taxon_id = c(9606L, 10090L),
      protein_count = c(2L, 1L), percentage = c(66.67, 33.33), stringsAsFactors = FALSE),
    organism_mapping = data.frame(uniprot_acc = c("P1", "P2", "P3"),
      organism_name = c("Homo sapiens", "Homo sapiens", "Mus musculus"), taxon_id = c(9606L, 9606L, 10090L),
      stringsAsFactors = FALSE),
    waiting_for_organism_selection = waiting,
    organism_unconfirmed = FALSE
  ))
}

# The success path with recording seams; `complete` writes the statuses as the frame's helper does.
k250CompleteSuccess <- function(workflowData, localData) {
  events <- character()
  result <- completeProtImportSuccessState(
    workflowData = workflowData, localData = localData, dataImportResult = k250ImportResult(), format = "diann",
    complete = function(notify = TRUE) {
      events <<- c(events, paste0("complete:notify=", notify))
      completeOmicImport(workflowData, omicGuiConfig("proteomics"), "pending", notify,
        showNotification = function(...) invisible(NULL))
    },
    removeModal = function() events <<- c(events, "removeModal"),
    showNotification = function(message, type = "default", duration = NULL) {
      events <<- c(events, paste0("notification:", message, ":", type))
    },
    messageFn = function(...) invisible(NULL)
  )
  list(result = result, events = events)
}

k250Confirm <- function(workflowData, localData, filter) {
  events <- character()
  seen <- list()
  confirmed <- k250Call(confirmProtImportOrganismSelection, list(
    selectedTaxon = 9606L, filterToOrganism = filter, workflowData = workflowData, localData = localData,
    session = structure(list(), class = "MockShinySession"),
    updateNumericInput = function(...) NULL, updateTextInput = function(...) NULL,
    showNotification = function(message, ...) events <<- c(events, paste0("notification:", message)),
    removeModal = function() events <<- c(events, "removeModal"),
    logInfo = function(...) NULL,
    complete = function(notify = TRUE) {
      events <<- c(events, paste0("complete:notify=", notify))
      seen <<- list(taxon_id = workflowData$taxon_id, rows = nrow(workflowData$data_tbl))
      completeOmicImport(workflowData, omicGuiConfig("proteomics"), "pending", notify,
        showNotification = function(...) invisible(NULL))
    }
  ))
  list(confirmed = confirmed, events = events, seen = seen)
}

k250Cancel <- function(workflowData, localData) {
  notes <- list()
  removed <- 0L
  if (!exists("cancelProtImportOrganismSelection", mode = "function")) return(NULL)
  cancelProtImportOrganismSelection(
    workflowData = workflowData, localData = localData,
    removeModal = function() removed <<- removed + 1L,
    showNotification = function(message, type = "default", duration = 5) {
      notes[[length(notes) + 1L]] <<- list(message = message, type = type, duration = duration)
    },
    logInfo = function(...) NULL
  )
  list(notes = notes, removed = removed)
}

test_that("E1a: with the organism dialog open the import is not complete (first import)", {
  workflowData <- k250WorkflowData()
  run <- k250CompleteSuccess(workflowData, k250LocalData(waiting = TRUE))
  expect_identical(run$events, character())
  expect_identical(workflowData$tab_status$setup_import, "incomplete")
  expect_identical(workflowData$tab_status$design_matrix, "disabled")
})

test_that("E1b: with the organism dialog open a second import is not complete either", {
  workflowData <- k250WorkflowData(list(setup_import = "complete", design_matrix = "pending",
    quality_control = "disabled"))
  run <- k250CompleteSuccess(workflowData, k250LocalData(waiting = TRUE))
  expect_false(identical(workflowData$tab_status$setup_import, "complete"))
  expect_false(any(grepl("^complete", run$events)))
  expect_identical(workflowData$tab_status$design_matrix, "pending")
})

test_that("E2: Confirm completes the import once, after the organism and the optional filter, before its notices", {
  for (filter in c(TRUE, FALSE)) {
    workflowData <- k250WorkflowData()
    localData <- k250LocalData(waiting = TRUE)
    k250CompleteSuccess(workflowData, localData)
    localData$organism_unconfirmed <- TRUE
    run <- k250Confirm(workflowData, localData, filter)
    expect_true(run$confirmed)
    expect_identical(sum(grepl("^complete", run$events)), 1L)
    expected <- c(if (filter) "notification:Filtered to Homo sapiens: kept 2 rows, removed 1 rows",
      "complete:notify=FALSE", "removeModal", "notification:Primary organism set to: Homo sapiens (Taxon: 9606)",
      "notification:Data import successful!")
    expect_identical(run$events, expected)
    expect_identical(run$seen$taxon_id, 9606L)
    expect_identical(run$seen$rows, if (filter) 2L else 3L)
    expect_identical(workflowData$tab_status$setup_import, "complete")
    expect_identical(workflowData$tab_status$design_matrix, "pending")
    expect_false(isTRUE(localData$organism_unconfirmed))
    expect_false(isTRUE(localData$waiting_for_organism_selection))
  }
})

test_that("E3: the dialog's Cancel is the module's own button, not a dismiss button", {
  modal <- buildProtImportOrganismSelectionModal(ns = shiny::NS("import"),
    organismDistribution = k250LocalData()$organism_distribution)
  footer <- htmltools::renderTags(modal)$html
  footer <- sub("^.*<div class=\"modal-footer\">", "", footer)
  expect_match(footer, "id=\"import-cancel_organism\"", fixed = TRUE)
  expect_false(grepl("data-dismiss", footer, fixed = TRUE))
  expect_false(grepl("data-bs-dismiss", footer, fixed = TRUE))
  expect_match(footer, "id=\"import-confirm_organism\"", fixed = TRUE)
})

test_that("E4: Cancel leaves the import unfinished, closes the dialog and says so until closed", {
  expect_true(exists("cancelProtImportOrganismSelection", mode = "function"))
  workflowData <- k250WorkflowData()
  localData <- k250LocalData(waiting = TRUE)
  run <- k250Cancel(workflowData, localData)
  if (is.null(run)) return(invisible(NULL))
  expect_false(isTRUE(localData$waiting_for_organism_selection))
  expect_true(isTRUE(localData$organism_unconfirmed))
  expect_identical(run$removed, 1L)
  expect_length(run$notes, 1L)
  expect_identical(run$notes[[1L]]$message, k250W1)
  expect_identical(run$notes[[1L]]$type, "warning")
  expect_null(run$notes[[1L]]$duration)
  expect_false(identical(workflowData$tab_status$setup_import, "complete"))
  expect_identical(workflowData$taxon_id, 10090L)
})

test_that("E5: a stale waiting flag is reset when an import starts and when an import fails", {
  localData <- k250LocalData(waiting = TRUE)
  localData$organism_unconfirmed <- TRUE
  startProtImportProcessing(searchResultsPath = "/tmp/report.tsv", fastaPath = "/tmp/proteins.fasta",
    format = "diann", localData = localData, ns = function(id) id, requireValue = function(value) value,
    announceStart = function(...) NULL, buildProcessingModal = function(ns) NULL, showModal = function(...) NULL)
  expect_false(isTRUE(localData$waiting_for_organism_selection))
  expect_false(isTRUE(localData$organism_unconfirmed))
  run <- k250CompleteSuccess(k250WorkflowData(), localData)
  expect_identical(run$events,
    c("complete:notify=FALSE", "removeModal", "notification:Data import successful!:message"))

  failed <- k250LocalData(waiting = TRUE)
  failed$organism_unconfirmed <- TRUE
  resetProtImportWorkflowStateOnError(workflowData = k250WorkflowData(), localData = failed, errorMessage = "boom",
    logError = function(...) NULL, removeModal = function() NULL, showNotification = function(...) NULL,
    messageFn = function(...) NULL)
  expect_false(isTRUE(failed$waiting_for_organism_selection))
  expect_false(isTRUE(failed$organism_unconfirmed))
})

test_that("E6: the import page shows no success and no organism before one is confirmed", {
  text <- function(tag) if (is.null(tag)) "" else as.character(htmltools::renderTags(tag)$html)
  status <- function(...) text(k250Call(buildProtImportStatusUi, list(isProcessing = FALSE,
    dataTbl = data.frame(x = 1), dataFormat = "diann", dataType = "peptide", ...)))
  waiting <- status(waitingForOrganism = TRUE)
  expect_match(waiting, k250W2, fixed = TRUE)
  expect_false(grepl("Data imported successfully", waiting, fixed = TRUE))
  cancelled <- status(organismUnconfirmed = TRUE)
  expect_match(cancelled, "Import not finished: no organism was chosen.", fixed = TRUE)
  expect_false(grepl("Data imported successfully", cancelled, fixed = TRUE))

  summary <- text(k250Call(buildProtImportDataSummaryUi, list(dataTbl = data.frame(x = 1),
    setupImportLog = list(detected_format = "diann", data_type = "peptide", n_rows = 8L, n_runs = 2L,
      n_proteins = 3L, n_peptides = NA, organism = "Mus musculus", taxon_id = 10090L),
    organismPending = TRUE)))
  expect_match(summary, "Organism: not confirmed", fixed = TRUE)
  expect_false(grepl("Mus musculus", summary, fixed = TRUE))
  expect_false(grepl("10090", summary, fixed = TRUE))
})

test_that("E7: after Cancel there is no table for the design builder and the import is incomplete", {
  expect_true(exists("cancelProtImportOrganismSelection", mode = "function"))
  workflowData <- k250WorkflowData(list(setup_import = "complete", design_matrix = "pending",
    quality_control = "disabled"))
  run <- k250Cancel(workflowData, k250LocalData(waiting = TRUE))
  if (is.null(run)) return(invisible(NULL))
  expect_null(workflowData$data_tbl)
  expect_null(workflowData$config_list)
  expect_null(workflowData$processing_log$setup_import)
  expect_identical(workflowData$tab_status$setup_import, "incomplete")
  expect_false(any(vapply(run$notes, function(n) identical(n$type, "error"), logical(1))))
})

test_that("G1: without the organism dialog the import completes as before", {
  workflowData <- k250WorkflowData()
  run <- k250CompleteSuccess(workflowData, k250LocalData(waiting = FALSE))
  expect_identical(run$events,
    c("complete:notify=FALSE", "removeModal", "notification:Data import successful!:message"))
  expect_identical(workflowData$tab_status$setup_import, "complete")
  expect_identical(workflowData$tab_status$design_matrix, "pending")
})

test_that("G2: Confirm writes the chosen organism, the mixed-species record and its two notices as before", {
  workflowData <- k250WorkflowData()
  run <- k250Confirm(workflowData, k250LocalData(waiting = TRUE), filter = FALSE)
  expect_true(run$confirmed)
  expect_identical(workflowData$taxon_id, 9606L)
  expect_identical(workflowData$organism_name, "Homo sapiens")
  expect_true(workflowData$mixed_species_analysis$enabled)
  expect_identical(workflowData$mixed_species_analysis$selected_taxon_id, 9606L)
  expect_identical(workflowData$mixed_species_analysis$selected_organism, "Homo sapiens")
  expect_false(workflowData$mixed_species_analysis$filter_applied_at_import)
  expect_identical(workflowData$processing_log$setup_import$taxon_id, 9606L)
  expect_identical(workflowData$processing_log$setup_import$organism, "Homo sapiens")
  expect_identical(workflowData$processing_log$setup_import$mixed_species_selection$selected_taxon_id, 9606L)
  expect_true("notification:Primary organism set to: Homo sapiens (Taxon: 9606)" %in% run$events)
  expect_true("notification:Data import successful!" %in% run$events)
})

test_that("G3: FASTA headers without organisms: no dialog, the import completes with the form's values", {
  workflowData <- list2env(list(organism_distribution = NULL))
  localData <- k250LocalData(waiting = FALSE)
  notes <- character()
  shown <- FALSE
  analysed <- analyzeProtImportMixedSpeciesData(
    workflowData = workflowData, localData = localData,
    dataImportResult = list(data = data.frame(Protein.Group = c("P1", "P2")),
      column_mapping = list(protein_col = "Protein.Group")),
    fastaPath = "/tmp/proteins.fasta", session = shiny::MockShinySession$new(),
    updateStatus = function(...) NULL, messageFn = function(...) NULL, logInfo = function(...) NULL,
    logWarn = function(...) NULL,
    showNotification = function(message, ...) notes <<- c(notes, message),
    extractOrganisms = function(...) data.frame(uniprot_acc = character(), taxon_id = integer(),
      organism_name = character()),
    analyzeDistribution = function(...) data.frame(taxon_id = integer(), organism_name = character(),
      protein_count = integer(), percentage = numeric()),
    buildSelectionModal = function(...) NULL,
    showModal = function(...) shown <<- TRUE
  )
  expect_false(analysed)
  expect_false(shown)
  expect_identical(notes, "Could not extract organism information from FASTA headers. Using default organism settings.")
  expect_false(isTRUE(localData$waiting_for_organism_selection))
  workflowData <- k250WorkflowData()
  run <- k250CompleteSuccess(workflowData, localData)
  expect_identical(run$events,
    c("complete:notify=FALSE", "removeModal", "notification:Data import successful!:message"))
  expect_identical(workflowData$taxon_id, 10090L)
})
