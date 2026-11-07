# fidelity-coverage-compare: shared
library(testthat)

getSharedProtDesignServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_design_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_design_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_design_server
}

test_that("proteomics design module server preserves public shell wiring", {
  server_fn <- getSharedProtDesignServer()
  server_env <- environment(server_fn)
  captured <- new.env(parent = emptyenv())
  captured$events <- character()

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids")
  workflow_data$design_matrix <- data.frame(Run = "S1", group = "A", stringsAsFactors = FALSE)
  workflow_data$contrasts_tbl <- data.frame(contrasts = "groupA-groupB", stringsAsFactors = FALSE)
  workflow_data$tab_status <- list(design_matrix = "pending")

  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir(),
    results_dir = tempdir()
  )
  volumes <- c(Home = tempdir())
  qc_value <- FALSE
  qc_trigger <- function(value) {
    if (!missing(value)) {
      qc_value <<- value
    }
    qc_value
  }

  input <- list()
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste("design", id, sep = "-"))

  has_target_shell <- exists(
    "registerProtDesignServerShells",
    envir = server_env,
    inherits = FALSE
  )
  if (has_target_shell) {
    local_mocked_bindings(
      registerProtDesignServerShells = function(input,
                                                output,
                                                session,
                                                workflowData,
                                                experimentPaths,
                                                volumes = NULL,
                                                qcTrigger = NULL,
                                                ...) {
        captured$shell_args <- list(
          input = input,
          output = output,
          session = session,
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          volumes = volumes,
          qcTrigger = qcTrigger
        )
        "builder-results-rv"
      },
      .env = environment()
    )
  } else if (exists("mod_prot_design_builder_server", envir = server_env, inherits = FALSE)) {
    local_mocked_bindings(
      mod_prot_design_builder_server = function(...) {
        captured$builder_registered <- TRUE
        function() NULL
      },
      .env = environment()
    )
  }

  local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
      invisible(NULL)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    reactiveVal = function(value = NULL) {
      stored <- value
      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      captured$events <- c(
        captured$events,
        paste(deparse(substitute(eventExpr)), collapse = "")
      )
      invisible(NULL)
    },
    reactive = function(expr) {
      quoted <- substitute(expr)
      env <- parent.frame()
      function() eval(quoted, env)
    },
    renderText = function(expr) {
      quoted <- substitute(expr)
      env <- parent.frame()
      function() eval(quoted, env)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      values[[1]]
    },
    outputOptions = function(...) invisible(NULL),
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL),
    showModal = function(...) invisible(NULL),
    removeModal = function(...) invisible(NULL),
    .package = "shiny"
  )
  local_mocked_bindings(
    shinyDirChoose = function(...) invisible(NULL),
    shinyFileChoose = function(...) invisible(NULL),
    .package = "shinyFiles"
  )
  local_mocked_bindings(
    renderDT = function(expr, ...) {
      quoted <- substitute(expr)
      env <- parent.frame()
      function() eval(quoted, env)
    },
    .package = "DT"
  )

  server_fn(
    id = "design",
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    volumes = volumes,
    qc_trigger = qc_trigger
  )

  expect_identical(captured$module_id, "design")

  if (has_target_shell) {
    expect_identical(captured$shell_args$input, input)
    expect_identical(captured$shell_args$output, output)
    expect_identical(captured$shell_args$session, session)
    expect_identical(captured$shell_args$workflowData, workflow_data)
    expect_identical(captured$shell_args$experimentPaths, experiment_paths)
    expect_identical(captured$shell_args$volumes, volumes)
    expect_identical(captured$shell_args$qcTrigger, qc_trigger)
  } else {
    expect_true("input$confirm_import" %in% captured$events)
    expect_true(is.function(output$data_available))
    expect_true(is.function(output$design_matrix_exists))
    expect_true(is.function(output$design_matrix_preview))
    expect_true(is.function(output$contrasts_preview))
  }
})

# KIT-103 W6: proteomics on the shared design tab (R/mod_omic_design.R) and builder server
# (R/mod_omic_design_builder.R) through the `dia_long` adapter, with its own builder UI as a plugin.
# The frozen values (UI digests, output names, texts and their order, the saved result, the data's
# digest and its cold reads) were recorded from the proteomics design files before W6 (de2d891c),
# and these blocks pass there unchanged (docs/qa/workflow-execution/kit/kit-103/w6/test-first.txt).

# The rendered HTML's sha256, with shiny's random tabset ids masked (as the entry suite does).
w6ProtDesignUiDigest <- function(tag) {
  html <- as.character(htmltools::renderTags(tag)$html)
  html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"<ID>\"", html)
  html <- gsub("tab-[0-9]+-([0-9]+)", "tab-<ID>-\\1", html)
  digest::digest(html, algo = "sha256", serialize = FALSE)
}

# The modal (by digest), notification and removal calls of a run, in order.
w6ProtDesignCaptureEvents <- function(env = parent.frame()) {
  seen <- new.env(parent = emptyenv())
  seen$events <- character()
  add <- function(event) seen$events <- c(seen$events, event)
  testthat::local_mocked_bindings(
    showNotification = function(ui, ..., type = "default", duration = 5, id = NULL) {
      add(sprintf("note[%s|%s|%s]: %s", type, if (is.null(duration)) "NULL" else duration, id %||% "",
        paste(as.character(ui), collapse = " ")))
      invisible(NULL)
    },
    removeNotification = function(id, ...) add(paste("remove_note:", id)),
    showModal = function(ui, ...) add(paste("modal", w6ProtDesignUiDigest(ui))),
    removeModal = function(...) add("remove_modal"),
    .package = "shiny", .env = env)
  seen
}

# A DIA owner holding a freshly imported peptide table (the import's references, no adoption pending).
w6ProtDesignOwner <- function(label, rows, envir = parent.frame()) {
  owner <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(.local_envir = envir), label, "DIA"),
    WorkflowState$new(audit_enabled = FALSE))
  owner$.context$reference_states <- TRUE
  withr::defer(closeDIAStateAdapter(owner), envir = envir)
  table <- prepareDIAParserTable(readr::read_tsv(I(paste0("Run\tProtein.Group\tPrecursor.Normalised\n",
    paste(rows, collapse = "\n"), "\n")), show_col_types = FALSE))
  mapping <- list(run_col = "Run", protein_col = "Protein.Group", quantity_col = "Precursor.Normalised")
  settings <- list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  references <- prepareDIAGUIReferences(owner, list(data_tbl = table, data_cln = table, data_type = "peptide",
    column_mapping = mapping, config_list = settings))
  owner$.context$design_adoption_pending <- NULL
  list(manager = owner, reference = references$data_tbl, mapping = mapping, settings = settings)
}

w6ProtDesignFrozen <- list(
  frame = "b44f4421225a1eb23923d82dda123acf353805dceeaae7ff1530e54286ce4401",
  builder = "73695d2d55355c17e51d08d7f4bc8622791d39cee144bcfce1380869f0de2da0",
  reset_modal = "569915cfbc737a3809925106609beeb062c0915cfc8bed7714d405199dcc2b7c",
  import_modal = "7b79b991cc614fcaed0e891ef3ad51deaf48b34ac022c9e02522a384c62ba5ee",
  processing_modal = "9ad0b3a27c566f3d714fe0c2f02153b20664ec93a92c7f64e1c0228502167270",
  saved_design = "ca8372feb95399452c3bdf39ee3d2d4e",
  saved_data = "d19c795299505d59ee598a34b802ca29c02297f6dcac657556ce82a59eda4319",
  saved_contrasts = "f3976dbb84b06f69bdc61b70e5189ef6",
  saved_config = "c36b9360aa43c39a08e22b919202b1a8",
  frame_outputs = c("builder-available_factors_display", "builder-contrast_factors_info", "builder-data_table",
    "builder-defined_contrasts_display", "builder-range_preview", "builder-removed_samples_display",
    "builder-replicate_inputs", "builder-tech_rep_summary", "contrasts_preview", "data_available",
    "design_matrix_exists", "design_matrix_preview", "fasta_detection_status", "import_dir_path"))

w6ProtDesignTableDigest <- function(x) {
  digest::digest(list(names = names(x), class = class(x), cols = lapply(as.list(as.data.frame(x)), identity)),
    algo = "sha256")
}

test_that("W6 the proteomics design tab and its builder plugin render today's bytes", {
  expect_identical(w6ProtDesignUiDigest(mod_prot_design_ui("design_matrix")), w6ProtDesignFrozen$frame)
  expect_identical(w6ProtDesignUiDigest(mod_prot_design_builder_ui("builder")), w6ProtDesignFrozen$builder)
})

test_that("W6 the proteomics builder over a DIA reference: texts, edits, save and reset are today's, the data cold", {
  w6_owner <- w6ProtDesignOwner("w6-builder", c("WT_2\tP1\t1", "KO_1\tP1\t2", "WT_1\tP1\t3", "KO_2\tP1\t4",
    "WT_2\tP2\t5", "KO_1\tP2\t6", "WT_1\tP2\t7", "KO_2\tP2\t8"))
  w6_reference <- shiny::reactiveVal(w6_owner$reference)
  w6_provider <- function() resolveDIAOwnedValue(w6_reference())
  attr(w6_provider, "dia_reference") <- w6_reference
  w6_seen <- w6ProtDesignCaptureEvents()
  shiny::testServer(mod_prot_design_builder_server, args = list(data_tbl = w6_provider,
    config_list = shiny::reactiveVal(w6_owner$settings), column_mapping = shiny::reactiveVal(w6_owner$mapping)), {
    session$userData$dia_capacity_manager <- w6_owner$manager
    w6_reads <- w6_owner$manager$.context$hydrates
    session$flushReact()
    expect_true(is.function(session$userData$peptide_checkpoint_design))
    session$setInputs(save_results = 1L)
    expect_null(session$returned())
    session$setInputs(formula_string = "~ 0 + group", new_factor = " WT ")
    session$setInputs(add_factor = 1L)
    session$setInputs(new_factor = "KO")
    session$setInputs(add_factor = 2L)
    expect_identical(as.character(output$available_factors_display$html), "<p>WT, KO</p>")
    session$setInputs(selected_runs = c("WT_1", "WT_2"), factor1_select = "WT", factor2_select = "",
      factor3_select = "", replicate_start = 1)
    expect_match(output$replicate_inputs$html, "Starting replicate number for 2 selected runs:", fixed = TRUE)
    expect_match(output$replicate_inputs$html, "data-testid=\"prot-design-replicate-start\"", fixed = TRUE)
    session$setInputs(assign_metadata = 1L)
    session$setInputs(selected_runs = c("KO_1", "KO_2"), factor1_select = "KO")
    session$setInputs(assign_metadata = 2L)
    session$setInputs(tech_rep_samples = "WT_1", tech_rep_assignment_mode = "lowest")
    session$setInputs(assign_tech_reps = 1L)
    session$setInputs(tech_rep_samples = c("WT_1", "KO_1"))
    session$setInputs(assign_tech_reps = 2L)
    session$setInputs(tech_rep_samples = c("KO_1", "KO_2"))
    session$setInputs(assign_tech_reps = 3L)
    expect_identical(output$tech_rep_summary,
      "Group: KO, Biological Replicate: 1\n  Samples: KO_1, KO_2\n  Technical Replicates: 1, 2")
    session$setInputs(contrast_group1 = "KO", contrast_group2 = "WT")
    session$setInputs(add_contrast = 1L)
    expect_identical(output$contrast_factors_info, paste0("Note: Contrasts will use 'group' prefix (e.g., ",
      "groupGA_Control-groupGA_Elevated)\nbased on current formula: ~ 0 + group"))
    expect_identical(as.character(output$defined_contrasts_display$html),
      "<p>\n  <code>KO_vs_WT=groupKO-groupWT</code>\n</p>")
    session$setInputs(formula_string = "~ group")
    expect_identical(output$contrast_factors_info,
      "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)")
    session$setInputs(formula_string = "~ 0 + group")
    # The design edits never read the table.
    expect_identical(w6_owner$manager$.context$hydrates, w6_reads)
    session$setInputs(samples_to_remove = "KO_2")
    session$setInputs(remove_samples = 1L)
    expect_identical(output$removed_samples_display, "Removed 1 sample(s):\nKO_2")
    session$setInputs(sample_to_rename = "WT_1", new_sample_name = "WT_A")
    session$setInputs(rename_sample = 1L)
    session$setInputs(samples_to_transform = "WT_2", transform_mode = "after_underscore", range_start = 1,
      range_end = 1)
    expect_identical(output$range_preview, "\"WT_2\" -> \"WT\"")
    session$setInputs(bulk_rename = 1L)
    session$setInputs(save_results = 2L)
    # The renames and the save rewrite the table by reference, without reading it.
    expect_identical(w6_owner$manager$.context$hydrates, w6_reads)
    w6_result <- shiny::isolate(session$returned())
    expect_identical(names(w6_result), c("design_matrix", "contrasts_tbl", "config_list", "data_reference"))
    expect_true(isSavedWorkflowReference(w6_result$data_reference))
    expect_identical(digest::digest(w6_result$design_matrix), w6ProtDesignFrozen$saved_design)
    expect_identical(typeof(w6_result$design_matrix$replicates), "integer")
    expect_identical(w6ProtDesignTableDigest(resolveDIAOwnedValue(w6_result$data_reference, w6_owner$manager)),
      w6ProtDesignFrozen$saved_data)
    expect_identical(digest::digest(w6_result$contrasts_tbl), w6ProtDesignFrozen$saved_contrasts)
    expect_identical(digest::digest(w6_result$config_list), w6ProtDesignFrozen$saved_config)
    expect_identical(diaStateOwnership(w6_owner$manager)$working_records, 0L)
    session$setInputs(reset_scope = "removed_samples")
    session$setInputs(reset_changes = 1L)
    session$setInputs(confirm_reset = 1L)
    expect_identical(output$removed_samples_display, "No samples have been removed.")
    expect_identical(as.character(output$available_factors_display$html),
      "<p>No factors defined yet (use the 'Factors' tab).</p>")
  })
  expect_identical(w6_seen$events, c(
    "note[warning|5|]: No samples have been assigned to groups. Please assign metadata before saving.",
    "note[warning|5|]: Please select at least two samples to assign as technical replicates.",
    "note[warning|5|]: All selected samples must belong to the same group. Please assign metadata first.",
    paste("note[message|5|]: Assigned 2 samples as technical replicates with biological replicate number 1",
      "and adjusted subsequent replicate numbers"),
    "note[message|5|]: Removed 1 sample(s) from analysis.",
    "note[message|5|]: Design saved successfully. You can close this builder.",
    paste("modal", w6ProtDesignFrozen$reset_modal), "remove_modal",
    "note[message|5|]: Reset of removed_samples completed."))
})

test_that("W6 the proteomics design tab's server: today's outputs and the commit hook's source-directory refusal", {
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
  w6_owner <- w6ProtDesignOwner("w6-frame", c("WT_2\tP1\t1", "KO_1\tP1\t2", "WT_1\tP1\t3", "KO_2\tP1\t4"))
  w6_workflow <- shiny::reactiveValues(data_tbl = w6_owner$reference, data_cln = w6_owner$reference,
    config_list = w6_owner$settings, column_mapping = w6_owner$mapping, design_matrix = NULL,
    contrasts_tbl = NULL, state_manager = w6_owner$manager, tab_status = list(design_matrix = "pending"))
  w6_qc <- NULL
  w6_seen <- w6ProtDesignCaptureEvents()
  suppressMessages(shiny::testServer(mod_prot_design_server, args = list(workflow_data = w6_workflow,
    experiment_paths = list(source_dir = file.path(tempdir(), "kit103-w6-absent")), volumes = c(Home = tempdir()),
    qc_trigger = function(value) w6_qc <<- value), {
    session$userData$dia_capacity_manager <- w6_owner$manager
    session$flushReact()
    expect_identical(sort(sub("^proxy[0-9]+-", "", names(session$.__enclos_env__$private$outs))),
      w6ProtDesignFrozen$frame_outputs)
    for (w6_name in c("data_available", "design_matrix_exists")) {
      w6_error <- tryCatch(output[[w6_name]], error = identity)
      expect_s3_class(w6_error, "error")
      expect_identical(conditionMessage(w6_error), "unused arguments (self, name)")
    }
    expect_identical(output$import_dir_path, "")
    session$setInputs(show_import_modal = 1L)
    session$setInputs(`builder-formula_string` = "~ 0 + group", `builder-new_factor` = "WT")
    session$setInputs(`builder-add_factor` = 1L)
    session$setInputs(`builder-new_factor` = "KO")
    session$setInputs(`builder-add_factor` = 2L)
    session$setInputs(`builder-selected_runs` = c("WT_1", "WT_2"), `builder-factor1_select` = "WT",
      `builder-factor2_select` = "", `builder-factor3_select` = "", `builder-replicate_start` = 1L)
    session$setInputs(`builder-assign_metadata` = 1L)
    session$setInputs(`builder-selected_runs` = c("KO_1", "KO_2"), `builder-factor1_select` = "KO")
    session$setInputs(`builder-assign_metadata` = 2L)
    session$setInputs(`builder-contrast_group1` = "KO", `builder-contrast_group2` = "WT")
    session$setInputs(`builder-add_contrast` = 1L)
    session$setInputs(`builder-save_results` = 1L)
  }))
  # Proteomics' own commit: the bus is written, then the missing source directory is refused; the
  # processing modal stays, the QC trigger and the tab status are untouched (today's).
  expect_identical(shiny::isolate(w6_workflow$design_matrix$Run), c("KO_1", "KO_2", "WT_1", "WT_2"))
  expect_identical(shiny::isolate(w6_workflow$design_matrix$replicates), c(1L, 2L, 1L, 2L))
  expect_true(isSavedWorkflowReference(shiny::isolate(w6_workflow$data_cln)))
  expect_identical(shiny::isolate(w6_workflow$contrasts_tbl$full_format), "KO_vs_WT=groupKO-groupWT")
  expect_identical(shiny::isolate(w6_workflow$tab_status), list(design_matrix = "pending"))
  expect_null(w6_qc)
  expect_identical(w6_seen$events, c(paste("modal", w6ProtDesignFrozen$import_modal),
    "note[message|5|]: Design saved successfully. You can close this builder.",
    paste("modal", w6ProtDesignFrozen$processing_modal),
    "note[error|15|]: Could not find the source directory to save files. Check experiment paths."))
})
