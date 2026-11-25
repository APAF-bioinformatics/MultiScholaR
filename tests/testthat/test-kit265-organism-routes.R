# KIT-265: the proteomics import with 'Mixed species FASTA' ticked, when no organism window can be shown. Route (i): the
# import continues and a persistent notice names the organism it will use; (i)(c) (the FASTA names organisms and none of
# the data's proteins is in it), (ii) an error counting organisms, (iii) no protein column and (iv) FASTA processing
# failed: the import is refused with a classed error, which the import's error path shows and uses to clear the import.
# Row ids are the ticket's (section 10.2). Calls drop the arguments a function does not have (`k265Call`), so that on
# the code before the fix each E test fails by assertion at the value the ticket states, not by an "unused argument".
# The message strings are the ticket's drafts (VAL-41; the owner's wording pass): `docs/qa/workflow-execution/kit/
# kit-265/wording.md`.

k265Head <- paste0(
  "No organism with a taxonomy ID could be matched to the proteins in the search results (the FASTA headers name no ",
  "organism or taxonomy ID for them, or the proteins are not in the FASTA file), so no organism window is shown. "
)
k265Named <- function(name, taxon) {
  sprintf("This analysis will use the organism entered in the form: %s (taxon %s). ", name, taxon)
}
k265Blank <- paste0(
  "The form has no organism (its name or Taxonomy ID is empty), so UniProt annotation and enrichment will have no ",
  "organism to use. "
)
k265Tail <- "To use another organism, untick 'Mixed species FASTA', enter it and press 'Process Imported Data' again."
k265Again <- "untick 'Mixed species FASTA', enter the organism yourself and press 'Process Imported Data' again."

k265Counting <- function(text) {
  paste0(
    "The import was stopped: the organisms in the FASTA file could not be counted (", text, "). The organism is ",
    "chosen from that count, so nothing was imported. Check the FASTA file, or ", k265Again
  )
}
k265NoColumn <- function(what) {
  paste0(
    "The import was stopped: the search results have no ", what, ", so their proteins cannot be matched to the ",
    "organisms in the FASTA file. Nothing was imported. Export the search results with that column and import them ",
    "again."
  )
}
k265FastaFailed <- paste0(
  "The import was stopped: the FASTA file could not be read or processed, so the organisms in it cannot be counted. ",
  "Nothing was imported. Check that the file is a FASTA file of protein sequences, or ", k265Again
)
k265NoneFound <- paste0(
  "The import was stopped: none of the proteins in the search results could be found in the FASTA file, so their ",
  "organisms cannot be counted. Nothing was imported. Check that this is the FASTA file the search used, or ",
  k265Again
)
k265WindowSentence <- paste0(
  "Only the organism you choose here is looked up in UniProt and tested for enrichment; the other organisms' ",
  "proteins are not, and they stay in the normalisation, RUV and the differential analysis (and in its correction ",
  "for multiple testing) unless you tick 'Filter data to keep only proteins from selected organism'."
)

k265Call <- function(fun, args) do.call(fun, args[names(args) %in% names(formals(fun))])

# The condition an expression signals, or NULL.
k265Err <- function(expr) tryCatch({ expr; NULL }, error = function(e) e)
k265IsRefusal <- function(err) !is.null(err) && inherits(err, "prot_import_organism_unavailable")
k265Message <- function(err) if (is.null(err)) "" else err$message

k265DataTable <- function(accessions = c("P1", "P2")) {
  data.frame(Protein.Group = accessions, value = seq_along(accessions), stringsAsFactors = FALSE)
}

# The non-owner shape (a table with the protein column named in the mapping).
k265Import <- function(accessions = c("P1", "P2")) {
  list(data = k265DataTable(accessions), column_mapping = list(protein_col = "Protein.Group"))
}

k265Fasta <- function(acc, organism, taxon) {
  data.frame(uniprot_acc = acc, organism_name = organism, taxon_id = as.integer(taxon), stringsAsFactors = FALSE)
}

k265EmptyFasta <- function() k265Fasta(character(), character(), integer())
k265EmptyDistribution <- function() {
  data.frame(taxon_id = integer(), organism_name = character(), protein_count = integer(), percentage = numeric())
}

# One call of `analyzeProtImportMixedSpeciesData`: the arguments the function has are passed, a notification and a
# modal are recorded, and the condition it signals (if any) is returned.
k265Run <- function(dataImportResult = k265Import(), extractOrganisms = function(...) k265EmptyFasta(),
                    analyzeDistribution = NULL, buildSelectionModal = NULL, form = list()) {
  workflowData <- list2env(list(organism_distribution = NULL))
  localData <- list2env(list(organism_mapping = NULL, organism_distribution = NULL,
    waiting_for_organism_selection = FALSE, organism_unconfirmed = FALSE))
  notes <- list()
  shown <- FALSE
  args <- list(
    workflowData = workflowData, localData = localData, dataImportResult = dataImportResult,
    fastaPath = "/tmp/proteins.fasta", session = shiny::MockShinySession$new(),
    updateStatus = function(...) NULL, messageFn = function(...) NULL, logInfo = function(...) NULL,
    logWarn = function(...) NULL,
    showNotification = function(message, type = "default", duration = 5, id = NULL) {
      notes[[length(notes) + 1L]] <<- list(message = message, type = type, duration = duration, id = id)
    },
    extractOrganisms = extractOrganisms,
    showModal = function(...) shown <<- TRUE
  )
  if (!is.null(analyzeDistribution)) args$analyzeDistribution <- analyzeDistribution
  if (!is.null(buildSelectionModal)) args$buildSelectionModal <- buildSelectionModal
  result <- NULL
  err <- k265Err(result <- k265Call(analyzeProtImportMixedSpeciesData, c(args, form)))
  list(result = result, err = err, notes = notes, shown = shown, localData = localData)
}

# The notice of route (i): one warning, persistent, under the fixed id, with the stated text.
k265ExpectNotice <- function(run, text) {
  expect_null(run$err)
  expect_false(isTRUE(run$result))
  expect_false(run$shown)
  expect_length(run$notes, 1L)
  note <- run$notes[[1L]]
  expect_identical(note$message, text)
  expect_identical(note$type, "warning")
  expect_null(note$duration)
  expect_identical(note$id, "prot_import_organism_cancelled")
  expect_false(isTRUE(run$localData$waiting_for_organism_selection))
}

k265ExpectRefusal <- function(run, text) {
  expect_true(k265IsRefusal(run$err))
  expect_identical(k265Message(run$err), text)
  expect_length(run$notes, 0L)
  expect_false(isTRUE(run$localData$waiting_for_organism_selection))
}

test_that("E1: route (i) names the form's organism in a persistent notice and the import continues", {
  run <- k265Run(form = list(organismName = "Mus musculus", taxonId = 10090))
  k265ExpectNotice(run, paste0(k265Head, k265Named("Mus musculus", "10090"), k265Tail))
})

test_that("E1b: route (i) with no form values passed says the form has no organism", {
  run <- k265Run()
  k265ExpectNotice(run, paste0(k265Head, k265Blank, k265Tail))
})

test_that("E1c: route (i) with a cleared form says the form has no organism and never prints 'taxon NA'", {
  run <- k265Run(form = list(organismName = "", taxonId = NA_real_))
  k265ExpectNotice(run, paste0(k265Head, k265Blank, k265Tail))
  expect_false(any(grepl("taxon NA", vapply(run$notes, `[[`, character(1), "message"), fixed = TRUE)))
})

test_that("E2: route (ii), an error counting organisms, refuses the import with the cause", {
  run <- k265Run(extractOrganisms = function(...) stop("boom"))
  k265ExpectRefusal(run, k265Counting("boom"))
})

test_that("E2c: the FASTA names organisms and none of the data's proteins is in it: the import is refused", {
  run <- k265Run(
    dataImportResult = k265Import(c("P1", "P2;P3")),
    extractOrganisms = function(...) {
      k265Fasta(c("Q1", "Q2"), c("Homo sapiens", "Mus musculus"), c(9606, 10090))
    },
    analyzeDistribution = analyzeOrganismDistribution,
    buildSelectionModal = buildProtImportOrganismSelectionModal,
    form = list(organismName = "Mus musculus", taxonId = 10090)
  )
  k265ExpectRefusal(run, k265NoneFound)
  expect_false(grepl("could not be counted", k265Message(run$err), fixed = TRUE))
  expect_false(run$shown)
})

test_that("E2d: proteins found under headers without organism fields are cause (a): the import continues", {
  run <- k265Run(
    dataImportResult = k265Import(c("P1", "P2;P3")),
    extractOrganisms = function(...) {
      k265Fasta(c("Q1", "P1", "P2", "P3"), c("Homo sapiens", NA, NA, NA), c(9606, NA, NA, NA))
    },
    analyzeDistribution = analyzeOrganismDistribution,
    buildSelectionModal = buildProtImportOrganismSelectionModal,
    form = list(organismName = "Mus musculus", taxonId = 10090)
  )
  k265ExpectNotice(run, paste0(k265Head, k265Named("Mus musculus", "10090"), k265Tail))
})

test_that("E3: route (iii), the column named in the mapping is absent from the table (non-owner shape)", {
  run <- k265Run(dataImportResult = list(
    data = data.frame(Other = c("P1", "P2")), column_mapping = list(protein_col = "PG.ProteinGroups")
  ))
  k265ExpectRefusal(run, k265NoColumn("'PG.ProteinGroups' column"))
  expect_false(grepl("could not be counted", k265Message(run$err), fixed = TRUE))
})

test_that("E3a: route (iii) with no protein column in the mapping says 'no protein column'", {
  run <- k265Run(dataImportResult = list(data = data.frame(Other = "P1"), column_mapping = list()))
  k265ExpectRefusal(run, k265NoColumn("protein column"))
})

test_that("E3b: route (iii), the owner shape: protein ids absent, the real distribution and window builder", {
  run <- k265Run(
    dataImportResult = list(column_mapping = list(protein_col = "PG.ProteinGroups"), import_summary = list(n_rows = 2L)),
    extractOrganisms = function(...) k265Fasta("Q1", "Homo sapiens", 9606),
    analyzeDistribution = analyzeOrganismDistribution,
    buildSelectionModal = buildProtImportOrganismSelectionModal
  )
  k265ExpectRefusal(run, k265NoColumn("'PG.ProteinGroups' column"))
  expect_false(grepl("could not be counted", k265Message(run$err), fixed = TRUE))
})

test_that("E4: route (iv), ticked with no FASTA table, refuses the import and does not call the analysis", {
  calls <- 0L
  err <- k265Err(runProtImportMixedSpeciesAnalysisIfNeeded(
    workflowData = list2env(list(aa_seq_tbl_final = NULL)), localData = list2env(list()),
    dataImportResult = k265Import(), fastaPath = "/tmp/proteins.fasta",
    session = structure(list(), class = "MockShinySession"), mixedSpeciesFasta = TRUE,
    analyzeMixedSpeciesData = function(...) {
      calls <<- calls + 1L
      TRUE
    }
  ))
  expect_true(k265IsRefusal(err))
  expect_identical(k265Message(err), k265FastaFailed)
  expect_identical(calls, 0L)
})

test_that("G4: ticked with a FASTA table the analysis is called once and its value returned; unticked it is not", {
  calls <- 0L
  run <- function(ticked, table) {
    runProtImportMixedSpeciesAnalysisIfNeeded(
      workflowData = list2env(list(aa_seq_tbl_final = table)), localData = list2env(list()),
      dataImportResult = k265Import(), fastaPath = "/tmp/proteins.fasta",
      session = structure(list(), class = "MockShinySession"), mixedSpeciesFasta = ticked,
      analyzeMixedSpeciesData = function(...) {
        calls <<- calls + 1L
        "analyzed"
      }
    )
  }
  expect_identical(run(TRUE, data.frame(id = "P1")), "analyzed")
  expect_identical(calls, 1L)
  expect_false(run(FALSE, NULL))
  expect_false(run(FALSE, data.frame(id = "P1")))
  expect_identical(calls, 1L)
})

test_that("E7: the organism window says that only the selected organism is looked up and tested", {
  distribution <- data.frame(organism_name = c("Homo sapiens", "Mus musculus"), taxon_id = c(9606L, 10090L),
    protein_count = c(2L, 1L), percentage = c(66.67, 33.33), stringsAsFactors = FALSE)
  html <- htmltools::renderTags(buildProtImportOrganismSelectionModal(ns = shiny::NS("import"),
    organismDistribution = distribution))$html
  html <- gsub("&#39;|&#x27;|&apos;", "'", html)
  sentence <- regexpr(k265WindowSentence, html, fixed = TRUE)[[1L]]
  radio <- regexpr("id=\"import-selected_organism\"", html, fixed = TRUE)[[1L]]
  filter <- regexpr("id=\"import-filter_to_organism\"", html, fixed = TRUE)[[1L]]
  expect_gt(sentence, 0L)
  expect_true(radio > 0L && filter > 0L && sentence > radio && sentence < filter)
})

# --- The observer: the module's own `process_data` observer in a Shiny test session -----------------------------
# Pattern of test-kit250-organism-dialog.R E8. The owner session is typed "DIA" (a DIA-NN import); the non-owner
# session (G6) has a plain WorkflowState and a Spectronaut peptide import, which no owner serves.

k265Diann <- function() {
  list(data = data.frame(Run = c("A", "B", "A"), Protein.Group = c("P1", "P2", "P3"),
      Stripped.Sequence = c("pep1", "pep2", "pep3"), Precursor.Normalised = c(1, 2, 3)),
    data_type = "peptide", column_mapping = list(protein_col = "Protein.Group", peptide_col = "Stripped.Sequence",
      run_col = "Run", quantity_col = "Precursor.Normalised", qvalue_col = "Q.Value", pg_qvalue_col = "PG.Q.Value"))
}

k265SpectronautPeptide <- function() {
  list(data = data.frame(Run = c("A", "B", "A"), PG.ProteinGroups = c("P1", "P2", "P3"),
      Stripped.Sequence = c("pep1", "pep2", "pep3"), Precursor.Normalised = c(1, 2, 3)),
    data_type = "peptide", column_mapping = list(protein_col = "PG.ProteinGroups", peptide_col = "Stripped.Sequence",
      run_col = "Run", quantity_col = "Precursor.Normalised"))
}

k265Events <- function() {
  ev <- new.env(parent = emptyenv())
  ev$log <- character()
  ev$notes <- list()
  ev$errors <- character()
  ev$warnings <- character()
  ev
}

k265Server <- function(root, ev, owner = TRUE, format = "diann", seams = list()) {
  imported <- if (identical(format, "diann")) k265Diann() else k265SpectronautPeptide()
  note <- function(message, type = "default", duration = 5, id = NULL, ...) {
    ev$notes[[length(ev$notes) + 1L]] <- list(message = message, type = type, duration = duration, id = id)
    ev$log <- c(ev$log, paste0("note:", type, ":", message))
  }
  function(input, output, session) {
    workflow <- shiny::reactiveValues(
      state_manager = if (owner) {
        newDIAStateAdapter(newWorkflowEntityStore(root, "kit265", "DIA"), WorkflowState$new(audit_enabled = FALSE))
      } else {
        WorkflowState$new(audit_enabled = FALSE)
      },
      tab_status = list(setup_import = "pending", design_matrix = "disabled"), processing_log = list())
    if (owner) attachDIAStateSession(workflow$state_manager, session, workflow)
    local_data <- createProtImportLocalData()
    defaults <- list(
      startProcessingFn = function(...) k265Call(startProtImportProcessing, utils::modifyList(list(...), list(
        announceStart = function(...) NULL,
        showModal = function(modal) ev$log <- c(ev$log, "showModal:processing"),
        removeNotification = function(id, ...) ev$log <- c(ev$log, paste0("removeNotification:", id))))),
      runProcessingSafelyFn = function(runProcessing, workflowData, localData, logError) {
        runProtImportProcessingSafely(runProcessing, workflowData, localData, logError = function(...) NULL,
          resetWorkflowState = function(workflowData, localData, errorMessage, logError) {
            ev$errors <- c(ev$errors, errorMessage)
            resetProtImportWorkflowStateOnError(workflowData, localData, errorMessage, logError = function(...) NULL,
              removeModal = function() ev$log <- c(ev$log, "removeModal:error"), showNotification = note,
              messageFn = function(...) NULL)
          })
      },
      runFormatDetectionFn = function(...) NULL,
      toggleMixedSpeciesInputsFn = function(...) NULL,
      readDataWithStatusFn = function(...) imported,
      loadConfigurationResourcesFn = function(workflowData, ...) workflowData$config_list <- list(fixture = TRUE),
      updateStatusFn = function(...) NULL,
      showNotificationFn = note,
      showModalFn = function(modal) ev$log <- c(ev$log, "showModal:organism"),
      completeSuccessStateFn = function(...) do.call(completeProtImportSuccessState, utils::modifyList(list(...),
        list(removeModal = function() ev$log <- c(ev$log, "removeModal:import"), showNotification = note,
          messageFn = function(...) NULL))),
      logInfo = function(...) NULL,
      logWarn = function(text) ev$warnings <- c(ev$warnings, text),
      messageFn = function(...) NULL
    )
    do.call(registerProtImportObservers, c(
      list(input = input, output = output, session = session, workflowData = workflow, localData = local_data,
        experimentPaths = list(source_dir = root), useShinyFiles = FALSE,
        searchResultsPathReactive = function() "fixture.tsv", fastaPathReactive = function() "fixture.fasta",
        activeFormatReactive = function() format),
      utils::modifyList(defaults, seams)
    ))
  }
}

# The stated clearing of a refused import (E5, E6, G6): the error shown, the window closed, nothing left on the bus.
k265ExpectRefused <- function(ev, workflow, local_data, text) {
  expect_identical(ev$errors, text)
  expect_true(paste0("note:error:Error: ", text) %in% ev$log)
  expect_true("removeModal:error" %in% ev$log)
  expect_false(any(grepl("^note:message:Data import successful", ev$log)))
  expect_null(workflow$data_tbl)
  expect_identical(workflow$tab_status$setup_import, "incomplete")
  expect_identical(workflow$tab_status$design_matrix, "disabled")
  expect_false(isTRUE(local_data$processing))
  expect_false(isTRUE(local_data$waiting_for_organism_selection))
  expect_false(isTRUE(local_data$organism_unconfirmed))
}

k265ExpectCompleted <- function(ev, workflow) {
  expect_identical(workflow$tab_status$setup_import, "complete")
  expect_identical(workflow$tab_status$design_matrix, "pending")
  expect_false(is.null(workflow$data_tbl))
  expect_true("note:message:Data import successful!" %in% ev$log)
}

k265FailingFasta <- function(...) stop("boom")

test_that("E5: the observer, ticked, FASTA processing failed: the import is refused; unticking and pressing again completes it", {
  root <- withr::local_tempdir()
  ev <- k265Events()
  shiny::testServer(k265Server(root, ev, owner = TRUE, seams = list(processFastaFn = k265FailingFasta)), {
    session$setInputs(sanitize_names = TRUE, organism_name = "Mus musculus", taxon_id = 10090,
      mixed_species_fasta = TRUE, process_data = 1L)
    k265ExpectRefused(ev, workflow, local_data, k265FastaFailed)

    # The refusal's own advice: untick and press again, in the same session.
    ev$log <- character()
    ev$errors <- character()
    session$setInputs(mixed_species_fasta = FALSE, process_data = 2L)
    expect_identical(ev$errors, character())
    k265ExpectCompleted(ev, workflow)
  })
})

test_that("E6: the observer, ticked, an error counting organisms: the same clearing, the error names the cause", {
  root <- withr::local_tempdir()
  ev <- k265Events()
  seams <- list(
    processFastaDataFn = function(workflowData, ...) workflowData$aa_seq_tbl_final <- data.frame(x = 1),
    extractOrganismsFn = function(...) stop("boom")
  )
  shiny::testServer(k265Server(root, ev, owner = TRUE, seams = seams), {
    session$setInputs(sanitize_names = TRUE, organism_name = "Mus musculus", taxon_id = 10090,
      mixed_species_fasta = TRUE, process_data = 1L)
    k265ExpectRefused(ev, workflow, local_data, k265Counting("boom"))
  })
})

test_that("E1d: the observer passes the form's organism and taxon to the route (i) notice", {
  root <- withr::local_tempdir()
  ev <- k265Events()
  seams <- list(
    processFastaDataFn = function(workflowData, ...) workflowData$aa_seq_tbl_final <- data.frame(x = 1),
    extractOrganismsFn = function(...) k265EmptyFasta()
  )
  shiny::testServer(k265Server(root, ev, owner = TRUE, seams = seams), {
    session$setInputs(sanitize_names = TRUE, organism_name = "Mus musculus", taxon_id = 10090,
      mixed_species_fasta = TRUE, process_data = 1L)
    expect_identical(ev$errors, character())
    notices <- Filter(function(n) identical(n$id, "prot_import_organism_cancelled"), ev$notes)
    expect_length(notices, 1L)
    expect_identical(
      if (length(notices)) notices[[1L]]$message else NA_character_,
      paste0(k265Head, k265Named("Mus musculus", "10090"), k265Tail)
    )
    k265ExpectCompleted(ev, workflow)
  })
})

test_that("G3: unticked with a failing FASTA step the import completes and the FASTA warning is logged as today", {
  root <- withr::local_tempdir()
  ev <- k265Events()
  shiny::testServer(k265Server(root, ev, owner = TRUE, seams = list(processFastaFn = k265FailingFasta)), {
    session$setInputs(sanitize_names = TRUE, organism_name = "Mus musculus", taxon_id = 10090,
      mixed_species_fasta = FALSE, process_data = 1L)
    expect_identical(ev$errors, character())
    expect_true("Error processing FASTA file: boom" %in% ev$warnings)
    k265ExpectCompleted(ev, workflow)
  })
})

test_that("G6: E5 on the non-owner path (no DIA owner, a Spectronaut peptide import): the same refusal and clearing", {
  root <- withr::local_tempdir()
  ev <- k265Events()
  shiny::testServer(k265Server(root, ev, owner = FALSE, format = "spectronaut",
    seams = list(processFastaFn = k265FailingFasta)), {
    session$setInputs(sanitize_names = TRUE, organism_name = "Mus musculus", taxon_id = 10090,
      spectronaut_quantity = "pep", mixed_species_fasta = TRUE, process_data = 1L)
    k265ExpectRefused(ev, workflow, local_data, k265FastaFailed)
  })
})
