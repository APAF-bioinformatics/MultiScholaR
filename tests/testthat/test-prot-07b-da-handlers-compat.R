# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

make_mock_protein_da_object <- function(args = list()) {
  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      S1 = c(10, 11),
      S2 = c(20, 21),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "uniprot_acc",
    protein_id_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      gene_names = c("GENE1", "GENE2"),
      stringsAsFactors = FALSE
    ),
    args = args
  )
}

with_global_contrasts_tbl <- function(value, code) {
  had_existing <- exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
  if (had_existing) {
    old_value <- get("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
  }

  if (is.null(value)) {
    if (had_existing) {
      rm("contrasts_tbl", envir = .GlobalEnv)
    }
  } else {
    assign("contrasts_tbl", value, envir = .GlobalEnv)
  }

  on.exit({
    if (had_existing) {
      assign("contrasts_tbl", old_value, envir = .GlobalEnv)
    } else if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
      rm("contrasts_tbl", envir = .GlobalEnv)
    }
  }, add = TRUE)

  force(code)
}

with_function_overrides <- function(env, replacements, code) {
  names_vec <- names(replacements)
  existed <- vapply(names_vec, exists, logical(1), envir = env, inherits = FALSE)
  old_values <- lapply(names_vec, function(name) {
    if (existed[[name]]) {
      get(name, envir = env, inherits = FALSE)
    } else {
      NULL
    }
  })
  names(old_values) <- names_vec
  locked <- vapply(names_vec, function(name) {
    existed[[name]] && bindingIsLocked(name, env)
  }, logical(1))

  for (name in names_vec) {
    if (existed[[name]] && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    assign(name, replacements[[name]], envir = env)
  }

  on.exit({
    for (name in names_vec) {
      if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
        unlockBinding(name, env)
      }
      if (existed[[name]]) {
        assign(name, old_values[[name]], envir = env)
      } else if (exists(name, envir = env, inherits = FALSE)) {
        rm(list = name, envir = env)
      }
      if (locked[[name]]) {
        lockBinding(name, env)
      }
    }
  }, add = TRUE)

  force(code)
}

volcanoMainFormals <- function() {
  has_volcano_main <- exists("writeInteractiveVolcanoPlotProteomicsMain", mode = "function", inherits = TRUE)
  if (has_volcano_main) {
    names(formals(writeInteractiveVolcanoPlotProteomicsMain))
  } else {
    character()
  }
}

skipIfMissingVolcanoMainInputSupport <- function() {
  volcano_main_formals <- volcanoMainFormals()
  skip_if_not(
    any(c("da_analysis_results_list", "de_analysis_results_list") %in% volcano_main_formals),
    "writeInteractiveVolcanoPlotProteomicsMain input support unavailable in this ref"
  )
}

make_mock_prot_da_state_manager <- function(obj, fail_save = FALSE) {
  manager <- new.env(parent = emptyenv())
  manager$states <- list(
    normalized = list(config = list(source = "test"), description = "normalized state")
  )
  manager$saved <- list()
  manager$getHistory <- function() {
    c("normalized")
  }
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    if (isTRUE(fail_save)) {
      stop("state save failed", call. = FALSE)
    }
    manager$saved[[length(manager$saved) + 1L]] <- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    manager$states[[state_name]] <- list(
      object = s4_data_object,
      config = config_object,
      description = description
    )
    invisible(TRUE)
  }
  manager
}

make_mock_da_handler_result <- function(theObject, contrast_string, include_qvalue_warning = FALSE) {
  comparison <- sub("=.*$", "", contrast_string)
  warning_payload <- NULL
  if (isTRUE(include_qvalue_warning)) {
    warning_payload <- stats::setNames(
      list("qvalue failed; p.adjust fallback used"),
      contrast_string
    )
  }

  list(
    theObject = theObject,
    da_proteins_long = data.frame(
      uniprot_acc = c("P1", "P2"),
      comparison = rep(comparison, 2),
      log2FC = c(1.25, -1.1),
      fdr_qvalue = c(0.01, 0.2),
      raw_pvalue = c(0.001, 0.05),
      stringsAsFactors = FALSE
    ),
    qvalue_warnings = warning_payload
  )
}

make_prot_da_run_handler_module <- function(da_data, workflow_data, experiment_paths) {
  force(da_data)
  force(workflow_data)
  force(experiment_paths)

  function(id) {
    shiny::moduleServer(id, function(input, output, session) {
      da_server_run_analysis_handler(
        input = input,
        output = output,
        session = session,
        ns = session$ns,
        da_data = da_data,
        workflow_data = workflow_data,
        experiment_paths = experiment_paths
      )
      da_data
    })
  }
}

make_prot_da_load_session_module <- function(da_data, workflow_data, experiment_paths) {
  force(da_data)
  force(workflow_data)
  force(experiment_paths)

  function(id) {
    shiny::moduleServer(id, function(input, output, session) {
      da_server_load_session_handler(
        input = input,
        output = output,
        session = session,
        da_data = da_data,
        workflow_data = workflow_data,
        experiment_paths = experiment_paths
      )
      da_data
    })
  }
}

test_that("da_server_load_session_handler restores exported normalization session state", {
  global_names <- c("contrasts_tbl", "config_list", "uniprot_dat_cln")
  old_globals <- lapply(global_names, function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      list(exists = TRUE, value = get(name, envir = .GlobalEnv, inherits = FALSE))
    } else {
      list(exists = FALSE, value = NULL)
    }
  })
  names(old_globals) <- global_names
  withr::defer({
    for (name in global_names) {
      if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
      if (isTRUE(old_globals[[name]]$exists)) {
        assign(name, old_globals[[name]]$value, envir = .GlobalEnv)
      }
    }
  }, envir = parent.frame())

  source_dir <- tempfile("da-load-session-")
  dir.create(source_dir, recursive = TRUE)
  obj <- make_mock_protein_da_object(
    args = list(deAnalysisParameters = list(formula_string = "~ 0 + group"))
  )
  session_data <- list(
    export_timestamp = as.POSIXct("2026-04-23 10:00:00", tz = "UTC"),
    r6_current_state_name = "correlation_filtered",
    final_protein_count = 2L,
    final_sample_count = 2L,
    r6_complete_states = list(raw_data_s4 = obj, correlation_filtered = obj),
    r6_state_history = c("raw_data_s4", "correlation_filtered"),
    contrasts_tbl = data.frame(
      contrasts = "groupA-groupB",
      full_format = "A_vs_B=groupA-groupB",
      stringsAsFactors = FALSE
    ),
    design_matrix = obj@design_matrix,
    config_list = list(globalParameters = list(workflow_type = "DIA")),
    current_s4_object = obj,
    fasta_metadata = list(fasta_format = "UniProt", num_sequences = 2L),
    accession_cleanup_results = list(cleanup_applied = TRUE, aggregation_method = "sum"),
    ruv_optimization_result = list(best_k = 2L, best_percentage = 10),
    qc_params = list(protein_qc = list(minimum = 1)),
    mixed_species_analysis = list(
      enabled = TRUE,
      selected_organism = "Human",
      selected_taxon_id = "9606"
    )
  )
  saveRDS(session_data, file.path(source_dir, "filtered_session_data_latest.rds"))
  saveRDS(
    data.frame(Entry = c("P1", "P2"), gene_names = c("GENE1", "GENE2"), stringsAsFactors = FALSE),
    file.path(source_dir, "uniprot_dat_cln.RDS")
  )

  da_data <- shiny::reactiveValues(
    current_s4_object = NULL,
    contrasts_available = NULL,
    formula_from_s4 = NULL
  )
  state_manager <- new.env(parent = emptyenv())
  state_manager$states <- list()
  state_manager$state_history <- character()
  state_manager$current_state <- NULL
  workflow_data <- shiny::reactiveValues(
    state_manager = state_manager,
    tab_status = list(normalization = "pending", differential_expression = "locked"),
    state_update_trigger = NULL
  )
  notifications <- list()
  text_updates <- list()

  testthat::local_mocked_bindings(
    withProgress = function(message, value, expr) {
      force(expr)
    },
    incProgress = function(amount, detail = NULL) {
      invisible(NULL)
    },
    showNotification = function(message, type = NULL, duration = NULL, ...) {
      notifications[[length(notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    updateTextAreaInput = function(session, inputId, value = NULL, ...) {
      text_updates[[inputId]] <<- value
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    make_prot_da_load_session_module(
      da_data,
      workflow_data,
      experiment_paths = list(source_dir = source_dir)
    ),
    {
      session$setInputs(load_filtered_session = 1)
      session$flushReact()

      expect_s4_class(da_data$current_s4_object, "ProteinQuantitativeData")
      expect_equal(da_data$contrasts_available, "groupA-groupB")
      expect_equal(da_data$formula_from_s4, "~ 0 + group")
      expect_equal(text_updates$formula_string, "~ 0 + group")
      expect_equal(workflow_data$state_manager$current_state, "correlation_filtered")
      expect_equal(workflow_data$tab_status$normalization, "complete")
      expect_equal(workflow_data$tab_status$differential_expression, "pending")
      expect_s3_class(workflow_data$state_update_trigger, "POSIXct")
      expect_equal(workflow_data$fasta_metadata$fasta_format, "UniProt")
      expect_true(workflow_data$accession_cleanup_results$cleanup_applied)
      expect_equal(workflow_data$ruv_optimization_result$best_k, 2L)
      expect_true(workflow_data$mixed_species_analysis$enabled)
      expect_equal(nrow(workflow_data$uniprot_dat_cln), 2L)
      expect_true(any(vapply(
        notifications,
        function(notification) identical(notification$type, "message"),
        logical(1)
      )))
    }
  )
})

test_that("writeInteractiveVolcanoPlotProteomicsMain accepts available DA/DE inputs", {
  skipIfMissingVolcanoMainInputSupport()

  obj <- make_mock_protein_da_object()
  results_list <- list(
    theObject = obj,
    da_proteins_long = data.frame(
      uniprot_acc = c("P1", "P2"),
      comparison = c("A_vs_B", "A_vs_B"),
      log2FC = c(1.5, -1.2),
      fdr_qvalue = c(0.01, 0.02),
      raw_pvalue = c(0.001, 0.002),
      stringsAsFactors = FALSE
    ),
    de_proteins_long = data.frame(
      uniprot_acc = c("P1", "P2"),
      comparison = c("A_vs_B", "A_vs_B"),
      log2FC = c(1.5, -1.2),
      fdr_qvalue = c(0.01, 0.02),
      raw_pvalue = c(0.001, 0.002),
      stringsAsFactors = FALSE
    ),
    contrasts_results = list(fit.eb = "fit-placeholder")
  )

  captured <- NULL
  wrapper_env <- environment(writeInteractiveVolcanoPlotProteomicsMain)
  volcano_main_formals <- volcanoMainFormals()
  uses_da_input <- "da_analysis_results_list" %in% volcano_main_formals
  uses_de_input <- "de_analysis_results_list" %in% volcano_main_formals
  main_args <- list(
    theObject = obj,
    uniprot_tbl = NULL,
    args_row_id = "uniprot_acc"
  )
  if (uses_da_input && uses_de_input) {
    main_args["da_analysis_results_list"] <- list(NULL)
    main_args$de_analysis_results_list <- results_list
    main_args$de_q_val_thresh <- 0.01
  } else if (uses_da_input) {
    main_args$da_analysis_results_list <- results_list
    main_args$da_q_val_thresh <- 0.01
  } else {
    main_args$de_analysis_results_list <- results_list
    main_args$de_q_val_thresh <- 0.01
  }

  with_function_overrides(
    wrapper_env,
    list(
      writeInteractiveVolcanoPlotProteomics = function(...) {
        captured <<- list(...)
        captured
      },
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) {
        default_value
      },
      updateParamInObject = function(theObject, param_name_string) {
        theObject
      }
    ),
    {
      result <- do.call(writeInteractiveVolcanoPlotProteomicsMain, main_args)

      expect_type(result, "list")
      expect_identical(result[[1]], results_list$de_proteins_long)
      expect_identical(result$fit.eb, "fit-placeholder")
      if (uses_da_input && uses_de_input && "da_q_val_thresh" %in% names(result)) {
        expect_equal(result$da_q_val_thresh, 0.01)
      } else if ("da_q_val_thresh" %in% names(result)) {
        expect_equal(result$da_q_val_thresh, 0.05)
      } else {
        expect_equal(result$de_q_val_thresh, 0.05)
      }
      expect_equal(rownames(result$counts_tbl), c("P1", "P2"))
      expect_equal(unname(result$groups), c("A", "B"))

      if (uses_da_input && uses_de_input) {
        fallback_results <- results_list
        fallback_results$theObject <- NULL
        fallback_args <- main_args
        fallback_args$de_analysis_results_list <- fallback_results

        fallback_result <- do.call(writeInteractiveVolcanoPlotProteomicsMain, fallback_args)
        expect_type(fallback_result, "list")
        expect_equal(rownames(fallback_result$counts_tbl), c("P1", "P2"))
        expect_equal(unname(fallback_result$groups), c("A", "B"))
      }
    }
  )
})

test_that("writeInteractiveVolcanoPlotProteomicsMain errors clearly for missing inputs", {
  volcano_main_formals <- volcanoMainFormals()
  skip_if_not(
    all(c("da_analysis_results_list", "de_analysis_results_list") %in% volcano_main_formals),
    "writeInteractiveVolcanoPlotProteomicsMain legacy alias validation unavailable in this ref"
  )

  obj <- make_mock_protein_da_object()
  wrapper_env <- environment(writeInteractiveVolcanoPlotProteomicsMain)
  with_function_overrides(
    wrapper_env,
    list(
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) {
        default_value
      },
      updateParamInObject = function(theObject, param_name_string) {
        theObject
      }
    ),
    {
      expect_error(
        writeInteractiveVolcanoPlotProteomicsMain(
          da_analysis_results_list = NULL,
          de_analysis_results_list = NULL,
          theObject = obj,
          uniprot_tbl = NULL
        ),
        "A DA/DE analysis results list must be supplied."
      )

      expect_error(
        writeInteractiveVolcanoPlotProteomicsMain(
          da_analysis_results_list = list(theObject = obj),
          theObject = obj,
          uniprot_tbl = NULL
        ),
        "Results list does not contain the expected DA/DE volcano inputs."
      )
    }
  )
})

test_that("da_server_run_analysis_handler processes full-format contrasts and qvalue warnings", {
  obj <- make_mock_protein_da_object()
  da_data <- shiny::reactiveValues(
    da_results_list = NULL,
    contrasts_available = c("groupA-groupB", "groupB-groupA"),
    analysis_complete = FALSE,
    current_s4_object = obj
  )
  state_manager <- make_mock_prot_da_state_manager(obj)
  workflow_data <- shiny::reactiveValues(
    state_manager = state_manager,
    tab_status = list(differential_abundance = "pending", enrichment_analysis = "locked"),
    da_analysis_results_list = NULL,
    da_ui_params = NULL
  )
  experiment_paths <- list(
    da_output_dir = tempfile("da-output-"),
    publication_graphs_dir = tempfile("publication-graphs-")
  )
  dir.create(experiment_paths$da_output_dir, recursive = TRUE)
  dir.create(experiment_paths$publication_graphs_dir, recursive = TRUE)

  calls <- new.env(parent = emptyenv())
  calls$analysis <- list()
  calls$output <- list()
  calls$checkpoints <- list()

  analysis_stub <- function(theObject, contrasts_tbl, formula_string, da_q_val_thresh,
                            treat_lfc_cutoff, qvalue_column, raw_pvalue_column) {
    contrast_string <- contrasts_tbl$contrasts[[1]]
    calls$analysis[[length(calls$analysis) + 1L]] <- list(
      contrast_string = contrast_string,
      formula_string = formula_string,
      da_q_val_thresh = da_q_val_thresh,
      treat_lfc_cutoff = treat_lfc_cutoff,
      qvalue_column = qvalue_column,
      raw_pvalue_column = raw_pvalue_column
    )
    make_mock_da_handler_result(
      theObject = theObject,
      contrast_string = contrast_string,
      include_qvalue_warning = identical(contrast_string, "A_vs_B=groupA-groupB")
    )
  }

  output_stub <- function(theObject, da_results_list_all_contrasts, uniprot_tbl,
                          da_output_dir, publication_graphs_dir, file_prefix,
                          args_row_id, gene_names_column, uniprot_id_column) {
    calls$output[[length(calls$output) + 1L]] <- list(
      theObject = theObject,
      da_results_list_all_contrasts = da_results_list_all_contrasts,
      uniprot_tbl = uniprot_tbl,
      da_output_dir = da_output_dir,
      publication_graphs_dir = publication_graphs_dir,
      file_prefix = file_prefix,
      args_row_id = args_row_id,
      gene_names_column = gene_names_column,
      uniprot_id_column = uniprot_id_column
    )
    TRUE
  }

  checkpoint_stub <- function(value, checkpoint_id, label) {
    calls$checkpoints[[length(calls$checkpoints) + 1L]] <- list(
      value = value,
      checkpoint_id = checkpoint_id,
      label = label
    )
    invisible(value)
  }

  contrasts_tbl <- data.frame(
    contrasts = c("groupA-groupB", "groupB-groupA"),
    full_format = c("A_vs_B=groupA-groupB", "B_vs_A=groupB-groupA"),
    friendly_names = c("A_vs_B", "B_vs_A"),
    stringsAsFactors = FALSE
  )

  with_global_contrasts_tbl(contrasts_tbl, {
    with_function_overrides(
      asNamespace("MultiScholaR"),
      list(
        differentialAbundanceAnalysis = analysis_stub,
        outputDaResultsAllContrasts = output_stub,
        .capture_checkpoint = checkpoint_stub
      ),
      {
        testServer(
          make_prot_da_run_handler_module(da_data, workflow_data, experiment_paths),
          {
            session$setInputs(
              formula_string = "~ 0 + group",
              da_q_val_thresh = 0.05,
              treat_lfc_cutoff = 0.25
            )
            session$setInputs(run_da_analysis = 1)
            session$flushReact()

            expect_true(da_data$analysis_complete)
            expect_equal(names(da_data$da_results_list$individual_contrasts), contrasts_tbl$contrasts)
            expect_equal(
              unique(da_data$da_results_list$da_proteins_long$comparison),
              c("A_vs_B", "B_vs_A")
            )
            expect_length(calls$analysis, 2L)
            expect_equal(
              vapply(calls$analysis, `[[`, character(1), "contrast_string"),
              contrasts_tbl$full_format
            )
            expect_equal(workflow_data$tab_status$differential_abundance, "complete")
            expect_equal(workflow_data$tab_status$enrichment_analysis, "pending")
            expect_equal(workflow_data$da_ui_params$q_value_threshold, 0.05)
            expect_equal(workflow_data$da_ui_params$log_fold_change_cutoff, 0.25)
            expect_true(workflow_data$da_ui_params$treat_enabled)
            expect_length(state_manager$saved, 2L)
            expect_length(calls$output, 1L)
            expect_equal(names(calls$output[[1]]$da_results_list_all_contrasts), contrasts_tbl$contrasts)
            expect_length(calls$checkpoints, 1L)
            expect_equal(calls$checkpoints[[1]]$checkpoint_id, "cp07")
          }
        )
      }
    )
  })
})

test_that("da_server_run_analysis_handler auto-generates contrast format and tolerates output write failure", {
  obj <- make_mock_protein_da_object()
  da_data <- shiny::reactiveValues(
    da_results_list = NULL,
    contrasts_available = "groupA-groupB",
    analysis_complete = FALSE,
    current_s4_object = obj
  )
  state_manager <- make_mock_prot_da_state_manager(obj, fail_save = TRUE)
  workflow_data <- shiny::reactiveValues(
    state_manager = state_manager,
    tab_status = list(differential_abundance = "pending", enrichment_analysis = "locked"),
    da_analysis_results_list = NULL,
    da_ui_params = NULL
  )
  experiment_paths <- list(
    da_output_dir = tempfile("da-output-"),
    publication_graphs_dir = tempfile("publication-graphs-")
  )

  calls <- new.env(parent = emptyenv())
  calls$analysis <- list()
  calls$output_attempted <- FALSE

  analysis_stub <- function(theObject, contrasts_tbl, formula_string, da_q_val_thresh,
                            treat_lfc_cutoff, qvalue_column, raw_pvalue_column) {
    calls$analysis[[length(calls$analysis) + 1L]] <- list(
      contrast_string = contrasts_tbl$contrasts[[1]],
      formula_string = formula_string,
      da_q_val_thresh = da_q_val_thresh,
      treat_lfc_cutoff = treat_lfc_cutoff
    )
    make_mock_da_handler_result(
      theObject = theObject,
      contrast_string = contrasts_tbl$contrasts[[1]]
    )
  }

  output_stub <- function(...) {
    calls$output_attempted <- TRUE
    stop("disk unavailable", call. = FALSE)
  }

  with_global_contrasts_tbl(NULL, {
    with_function_overrides(
      asNamespace("MultiScholaR"),
      list(
        differentialAbundanceAnalysis = analysis_stub,
        outputDaResultsAllContrasts = output_stub,
        .capture_checkpoint = function(value, checkpoint_id, label) invisible(value)
      ),
      {
        testServer(
          make_prot_da_run_handler_module(da_data, workflow_data, experiment_paths),
          {
            session$setInputs(
              formula_string = "~ 0 + group",
              da_q_val_thresh = 0.1,
              treat_lfc_cutoff = 0
            )
            session$setInputs(run_da_analysis = 1)
            session$flushReact()

            expect_true(da_data$analysis_complete)
            expect_equal(calls$analysis[[1]]$contrast_string, "A_vs_B=groupA-groupB")
            expect_equal(unique(da_data$da_results_list$da_proteins_long$comparison), "A_vs_B")
            expect_true(calls$output_attempted)
            expect_length(state_manager$saved, 0L)
            expect_false(workflow_data$da_ui_params$treat_enabled)
            expect_equal(workflow_data$tab_status$differential_abundance, "complete")
          }
        )
      }
    )
  })
})

test_that("da_server_run_analysis_handler reports analysis errors without completing state", {
  obj <- make_mock_protein_da_object()
  da_data <- shiny::reactiveValues(
    da_results_list = NULL,
    contrasts_available = "groupA-groupB",
    analysis_complete = FALSE,
    current_s4_object = obj
  )
  workflow_data <- shiny::reactiveValues(
    state_manager = make_mock_prot_da_state_manager(obj),
    tab_status = list(differential_abundance = "pending", enrichment_analysis = "locked"),
    da_analysis_results_list = NULL,
    da_ui_params = NULL
  )
  experiment_paths <- list(
    da_output_dir = tempfile("da-output-"),
    publication_graphs_dir = tempfile("publication-graphs-")
  )

  with_global_contrasts_tbl(
    data.frame(
      contrasts = "groupA-groupB",
      full_format = "A_vs_B=groupA-groupB",
      friendly_names = "A_vs_B",
      stringsAsFactors = FALSE
    ),
    {
      with_function_overrides(
        asNamespace("MultiScholaR"),
        list(
          differentialAbundanceAnalysis = function(...) {
            stop("limma failed", call. = FALSE)
          },
          outputDaResultsAllContrasts = function(...) {
            stop("output should not be called", call. = FALSE)
          },
          .capture_checkpoint = function(value, checkpoint_id, label) invisible(value)
        ),
        {
          testServer(
            make_prot_da_run_handler_module(da_data, workflow_data, experiment_paths),
            {
              session$setInputs(
                formula_string = "~ 0 + group",
                da_q_val_thresh = 0.05,
                treat_lfc_cutoff = 0
              )
              session$setInputs(run_da_analysis = 1)
              session$flushReact()

              expect_false(da_data$analysis_complete)
              expect_null(da_data$da_results_list)
              expect_equal(workflow_data$tab_status$differential_abundance, "pending")
            }
          )
        }
      )
    }
  )
})

test_that("da_server_init_handlers auto-generates prefixed contrasts from workflow state", {
  with_global_contrasts_tbl(NULL, {
    selected_tab <- shiny::reactiveVal(NULL)
    obj <- make_mock_protein_da_object(
      args = list(
        deAnalysisParameters = list(
          formula_string = "~ 0 + group"
        )
      )
    )
    workflow_data <- shiny::reactiveValues(
      state_manager = list(
        current_state = "normalized",
        getState = function(state) obj
      ),
      state_update_trigger = NULL
    )

    init_test_module <- function(id, workflow_data, selected_tab) {
      shiny::moduleServer(id, function(input, output, session) {
        da_data <- shiny::reactiveValues(
          da_results_list = NULL,
          contrasts_available = NULL,
          analysis_complete = FALSE,
          current_s4_object = NULL,
          formula_from_s4 = NULL,
          current_row_clusters = NULL,
          current_col_clusters = NULL
        )

        da_server_init_handlers(input, output, session, da_data, workflow_data, selected_tab)
        da_data
      })
    }

    testServer(
      init_test_module,
      args = list(workflow_data = workflow_data, selected_tab = selected_tab),
      {
        session$setInputs(formula_string = "~ 0 + group")
        selected_tab("da")
        session$flushReact()

        expect_s4_class(da_data$current_s4_object, "ProteinQuantitativeData")
        expect_equal(da_data$formula_from_s4, "~ 0 + group")
        expect_equal(as.vector(da_data$contrasts_available), "groupA-groupB")
      }
    )
  })
})

test_that("da_server_init_handlers normalizes user contrasts to group-prefixed format", {
  with_global_contrasts_tbl(
    data.frame(comparison = "A-B", stringsAsFactors = FALSE),
    {
      selected_tab <- shiny::reactiveVal(NULL)
      obj <- make_mock_protein_da_object(
        args = list(
          deAnalysisParameters = list(
            formula_string = "~ 0 + group"
          )
        )
      )
      workflow_data <- shiny::reactiveValues(
        state_manager = list(
          current_state = "normalized",
          getState = function(state) obj
        ),
        state_update_trigger = NULL
      )

      init_test_module <- function(id, workflow_data, selected_tab) {
        shiny::moduleServer(id, function(input, output, session) {
          da_data <- shiny::reactiveValues(
            da_results_list = NULL,
            contrasts_available = NULL,
            analysis_complete = FALSE,
            current_s4_object = NULL,
            formula_from_s4 = NULL,
            current_row_clusters = NULL,
            current_col_clusters = NULL
          )

          da_server_init_handlers(input, output, session, da_data, workflow_data, selected_tab)
          da_data
        })
      }

      testServer(
        init_test_module,
        args = list(workflow_data = workflow_data, selected_tab = selected_tab),
        {
          session$setInputs(formula_string = "~ 0 + group")
          selected_tab("da")
          session$flushReact()

          expect_equal(as.vector(da_data$contrasts_available), "groupA-groupB")
        }
      )
    }
  )
})
