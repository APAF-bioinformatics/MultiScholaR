# KIT-256: the design builder shared by every omic kept the technical-replicate selection after a successful Assign
# (the selection observer re-sends it after every design change). Adding the next sample's runs to the box and pressing
# Assign then made two biological samples one, with no refusal. The owner's decision (methods ledger VAL-249): the
# selection is emptied after each successful Assign; a refused Assign keeps it.
# What the server sends to the box is read from the GUI-entry browser double (helper-gui-entry.R), which queues every
# update*Input message: in testServer, updateSelectizeInput does not change `input`. E rows fail on the base commit at
# the ticket's base value (section 10); G rows hold on both. Each test prints what it observed ("KIT256-OBS" lines).
library(testthat)

kit256Obs <- function(id, value) {
  cat(sprintf("KIT256-OBS %s: %s\n", id, if (is.null(value)) "<no message>" else
    if (length(value) == 0L) "character(0)" else paste(value, collapse = " ")))
  invisible(value)
}

# The value of the last message the server sent for `inputId` (namespace removed) from message `from` on; NULL when
# none was sent.
kit256Last <- function(log, inputId, from = 1L) {
  updates <- log$updates[seq_along(log$updates) >= from]
  hits <- Filter(function(update) identical(sub("^.*-", "", update$id), inputId), updates)
  if (length(hits) == 0L) return(NULL)
  hits[[length(hits)]]$value
}

kit256NewLog <- function() {
  log <- new.env(parent = emptyenv())
  log$updates <- list(); log$notes <- list(); log$modals <- integer(); log$choices <- list()
  log
}

# KIT-222's clicks (test-kit222-replicate-keys.R:228-244, kit222BuilderClicks), copied and split at its technical-
# replicate clicks: the two groups (Ctrl S1 to S5, Trt S6 to S8) before, the contrast, formula and save after.
kit256SetupClicks <- function(set) {
  set(new_factor = "Ctrl"); set(add_factor = 1L)
  set(new_factor = "Trt"); set(add_factor = 2L)
  set(selected_runs = paste0("S", 1:5), factor1_select = "Ctrl", factor2_select = "", factor3_select = "",
    replicate_start = 1L)
  set(assign_metadata = 1L)
  set(selected_runs = paste0("S", 6:8), factor1_select = "Trt", replicate_start = 1L)
  set(assign_metadata = 2L)
}

kit256SaveClicks <- function(set) {
  set(contrast_group1 = "Trt", contrast_group2 = "Ctrl")
  set(add_contrast = 1L)
  set(formula_string = "~ 0 + group")
  set(save_results = 1L)
}

# The inputs of KIT-222's shared-builder save (test-kit222-replicate-keys.R:246-266), copied.
kit256SharedInputs <- function(omic, runs) {
  switch(omic,
    lipidomics = list(assays = module_ci_lipid_design_assays(layout = "single", samples = runs),
      mapping = module_ci_lipid_design_column_mapping(runs)),
    transcriptomics = list(assays = list(counts = cbind(data.frame(GeneID = c("G1", "G2"), GeneName = c("N1", "N2")),
      as.data.frame(matrix(1:16, nrow = 2L, dimnames = list(NULL, runs))))),
      mapping = list(lipid_id_col = "GeneID", annotation_col = "GeneName", sample_columns = runs)))
}

# The shared builder (lipidomics or transcriptomics configuration) under the browser double. `clicks(api)` drives it:
# `api$set(...)` sets inputs; `api$mark()` is the index of the next message; `api$apply()` applies the queued messages
# to the inputs as a browser does (the last value per input, helper-gui-entry.R:301-308); `api$design()` reads the
# builder's working design. Returns the log and the saved result.
kit256SharedBuilder <- function(omic, clicks, runs = paste0("S", 1:8)) {
  input_data <- kit256SharedInputs(omic, runs)
  log <- kit256NewLog()
  guiEntryInstallBrowserDouble(log)
  log$applied <- 0L
  seen <- new.env(parent = emptyenv())
  suppressMessages(shiny::testServer(mod_omic_design_builder_server, args = list(data_tbl = function() input_data$assays,
    config_list = shiny::reactive(list(deAnalysisParameters = list(formula_string = "~ 0 + group"))),
    column_mapping = shiny::reactive(input_data$mapping), config = omicGuiConfig(omic)), {
    session$flushReact()
    apply_updates <- function() {
      for (i in seq_len(50L)) {
        if (length(log$updates) <= log$applied) return(invisible())
        pending <- log$updates[seq.int(log$applied + 1L, length(log$updates))]
        log$applied <- length(log$updates)
        values <- stats::setNames(lapply(pending, `[[`, "value"), sub("^.*-", "", vapply(pending, `[[`, "", "id")))
        values <- values[!duplicated(names(values), fromLast = TRUE)]
        do.call(session$setInputs, values)
      }
    }
    api <- list(log = log, set = function(...) session$setInputs(...), mark = function() length(log$updates) + 1L,
      apply = apply_updates, design = function() shiny::isolate(design_matrix()),
      input = function(id) shiny::isolate(input[[id]]))
    clicks(api)
    seen$saved <- shiny::isolate(session$returned())
  }))
  list(log = log, saved = seen$saved)
}

# A DIA owner holding a freshly imported peptide table (test-kit222-replicate-keys.R:282-296, kit222ProtOwner), copied.
kit256ProtOwner <- function(label, runs, envir = parent.frame()) {
  owner <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = envir), label, "DIA"),
    WorkflowState$new(audit_enabled = FALSE))
  owner$.context$reference_states <- TRUE
  withr::defer(closeDIAStateAdapter(owner), envir = envir)
  rows <- sprintf("%s\tP1\t%d", runs, seq_along(runs))
  table <- prepareDIAParserTable(readr::read_tsv(I(paste0("Run\tProtein.Group\tPrecursor.Normalised\n",
    paste(rows, collapse = "\n"), "\n")), show_col_types = FALSE))
  mapping <- list(run_col = "Run", protein_col = "Protein.Group", quantity_col = "Precursor.Normalised")
  settings <- list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  references <- prepareDIAGUIReferences(owner, list(data_tbl = table, data_cln = table, data_type = "peptide",
    column_mapping = mapping, config_list = settings))
  owner$.context$design_adoption_pending <- NULL
  list(manager = owner, reference = references$data_tbl, mapping = mapping, settings = settings)
}

# The proteomics builder over a DIA reference (test-kit222-replicate-keys.R:299-318, kit222ProtBuilder), copied:
# `clicks(set)` drives it; returns the saved result and the notifications.
kit256ProtBuilder <- function(runs, clicks) {
  owner <- kit256ProtOwner("kit256-prot", runs)
  reference <- shiny::reactiveVal(owner$reference)
  provider <- function() resolveDIAOwnedValue(reference())
  attr(provider, "dia_reference") <- reference
  seen <- new.env(parent = emptyenv())
  seen$notes <- character()
  local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
    seen$notes <- c(seen$notes, sprintf("%s: %s", type, paste(as.character(ui), collapse = " ")))
    invisible(NULL)
  }, .package = "shiny")
  suppressMessages(shiny::testServer(mod_prot_design_builder_server, args = list(data_tbl = provider,
    config_list = shiny::reactiveVal(owner$settings), column_mapping = shiny::reactiveVal(owner$mapping)), {
    session$userData$dia_capacity_manager <- owner$manager
    session$flushReact()
    clicks(function(...) session$setInputs(...))
    seen$saved <- shiny::isolate(session$returned())
  }))
  list(saved = seen$saved, notes = seen$notes)
}

# The runs of `runs` grouped by their saved key, independent of the key labels.
kit256Partition <- function(saved, runs) {
  design <- saved$design_matrix
  keys <- split(design$Run, design$tech_rep_group)
  parts <- lapply(unname(keys), function(key) sort(intersect(key, runs)))
  parts <- parts[lengths(parts) > 0L]
  parts[order(vapply(parts, `[[`, character(1), 1L))]
}

kit256AssignFirstPair <- function(api) {
  kit256SetupClicks(api$set)
  api$set(tech_rep_samples = c("S1", "S2"), tech_rep_assignment_mode = "lowest")
  api$set(assign_tech_reps = 1L)
}

test_that("KIT-256 E1: after a successful Assign the box is emptied (lipidomics, transcriptomics)", {
  for (omic in c("lipidomics", "transcriptomics")) {
    run <- kit256SharedBuilder(omic, kit256AssignFirstPair)
    last <- kit256Obs(paste("E1", omic), kit256Last(run$log, "tech_rep_samples"))
    expect_identical(last, character(0), label = paste(omic, "last message for tech_rep_samples"))
  }
})

test_that("KIT-256 E2: after a successful Assign the box is emptied (proteomics, dia_long)", {
  log <- kit256NewLog()
  guiEntryInstallBrowserDouble(log)
  result <- kit256ProtBuilder(paste0("S", 1:8), function(set) {
    kit256SetupClicks(set)
    set(tech_rep_samples = c("S1", "S2"), tech_rep_assignment_mode = "lowest")
    set(assign_tech_reps = 1L)
  })
  expect_match(utils::tail(result$notes, 1L), "^message: ")
  last <- kit256Obs("E2 proteomics", kit256Last(log, "tech_rep_samples"))
  expect_identical(last, character(0))
})

test_that("KIT-256 E3: the next sample's runs added to what the box shows give two samples, not one", {
  for (omic in c("lipidomics", "transcriptomics")) {
    run <- kit256SharedBuilder(omic, function(api) {
      kit256AssignFirstPair(api)
      api$apply()
      shown <- as.character(api$input("tech_rep_samples"))
      kit256Obs(paste("E3", omic, "box after the first Assign"), shown)
      api$set(tech_rep_samples = union(shown, c("S3", "S4")))
      api$set(assign_tech_reps = 2L)
      kit256SaveClicks(api$set)
    })
    expect_type(run$saved, "list")
    parts <- kit256Partition(run$saved, paste0("S", 1:4))
    kit256Obs(paste("E3", omic, "keys of S1 to S4"), vapply(parts, paste, "", collapse = ","))
    expect_identical(parts, list(c("S1", "S2"), c("S3", "S4")), label = paste(omic, "samples of S1 to S4"))
  }
})

test_that("KIT-256 G1: a refused Assign sends no clear; the next design change re-sends the user's selection", {
  run <- kit256SharedBuilder("lipidomics", function(api) {
    kit256SetupClicks(api$set)
    # Too few runs.
    api$set(tech_rep_samples = "S1", tech_rep_assignment_mode = "lowest")
    mark <- api$mark()
    api$set(assign_tech_reps = 1L)
    kit256Obs("G1 too_few after the click", kit256Last(api$log, "tech_rep_samples", mark))
    expect_false(identical(kit256Last(api$log, "tech_rep_samples", mark), character(0)))
    expect_identical(tail(api$log$notes, 1L)[[1L]]$type, "warning")
    mark <- api$mark()
    api$set(sample_to_rename = "S8", new_sample_name = "S8b")
    api$set(rename_sample = 1L)
    expect_identical(kit256Obs("G1 too_few after a rename", kit256Last(api$log, "tech_rep_samples", mark)), "S1")
    # Two groups.
    api$set(tech_rep_samples = c("S5", "S6"))
    mark <- api$mark()
    api$set(assign_tech_reps = 2L)
    kit256Obs("G1 groups after the click", kit256Last(api$log, "tech_rep_samples", mark))
    expect_false(identical(kit256Last(api$log, "tech_rep_samples", mark), character(0)))
    expect_identical(tail(api$log$notes, 1L)[[1L]]$type, "warning")
    mark <- api$mark()
    api$set(sample_to_rename = "S7", new_sample_name = "S7b")
    api$set(rename_sample = 2L)
    expect_identical(kit256Obs("G1 groups after a rename", kit256Last(api$log, "tech_rep_samples", mark)),
      c("S5", "S6"))
  })
})

test_that("KIT-256 G2: after an Assign, an unrelated design change keeps a new selection (the flag is not left set)", {
  run <- kit256SharedBuilder("lipidomics", function(api) {
    kit256AssignFirstPair(api)
    api$set(tech_rep_samples = c("S5", "S6"))
    mark <- api$mark()
    api$set(sample_to_rename = "S8", new_sample_name = "S8b")
    api$set(rename_sample = 1L)
    expect_identical(kit256Obs("G2", kit256Last(api$log, "tech_rep_samples", mark)), c("S5", "S6"))
  })
})

test_that("KIT-256 G3: an Assign that leaves the design identical also empties the box", {
  run <- kit256SharedBuilder("lipidomics", function(api) {
    kit256AssignFirstPair(api)
    before <- api$design()
    api$set(tech_rep_samples = c("S1", "S2"), tech_rep_assignment_mode = "lowest")
    api$set(assign_tech_reps = 2L)
    # Precondition: the design is identical, so the row exercises the Assign observer's own clear, not the flag.
    expect_identical(api$design(), before)
    expect_identical(kit256Obs("G3", kit256Last(api$log, "tech_rep_samples")), character(0))
    # G2's check then holds.
    api$set(tech_rep_samples = c("S5", "S6"))
    mark <- api$mark()
    api$set(sample_to_rename = "S8", new_sample_name = "S8b")
    api$set(rename_sample = 1L)
    expect_identical(kit256Obs("G3 then G2", kit256Last(api$log, "tech_rep_samples", mark)), c("S5", "S6"))
  })
})

test_that("KIT-256 G4: Remove Samples still empties its own box", {
  run <- kit256SharedBuilder("lipidomics", function(api) {
    kit256SetupClicks(api$set)
    api$set(samples_to_remove = "S7")
    mark <- api$mark()
    api$set(remove_samples = 1L)
    last <- kit256Obs("G4", kit256Last(api$log, "samples_to_remove", mark))
    expect_false(is.null(last))
    expect_identical(setdiff(as.character(last), ""), character(0))
  })
})
