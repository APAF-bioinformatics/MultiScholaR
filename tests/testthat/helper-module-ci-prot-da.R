module_ci_prot_da_design <- function(samples = c("S1", "S2", "S3", "S4"),
                                     groups = c("A", "A", "B", "B"),
                                     batches = c("Batch1", "Batch2", "Batch1", "Batch2")) {
  design <- data.frame(
    Run = samples,
    group = rep(groups, length.out = length(samples)),
    batch = rep(batches, length.out = length(samples)),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$Run
  design
}

module_ci_prot_da_table <- function(samples = c("S1", "S2", "S3", "S4"),
                                    features = paste0("P", 1:6)) {
  if (length(samples) == 4L && length(features) == 6L) {
    values <- matrix(
      c(
        10, 11, 22, 23,
        12, 12, 21, 22,
        24, 23, 11, 10,
        15, 15, 15, 15,
        8, 9, 10, 11,
        30, 31, 34, 35
      ),
      nrow = length(features),
      byrow = TRUE,
      dimnames = list(features, samples)
    )
  } else {
    values <- outer(seq_along(features), seq_along(samples), function(feature_idx, sample_idx) {
      10 * feature_idx + sample_idx
    })
    dimnames(values) <- list(features, samples)
  }
  data.frame(
    uniprot_acc = features,
    as.data.frame(values, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

module_ci_prot_da_object <- function(workflow_type = "DIA",
                                     use_limpa = FALSE,
                                     report_template = paste0(workflow_type, "_report"),
                                     design = module_ci_prot_da_design(),
                                     protein_quant_table = module_ci_prot_da_table(samples = design$Run),
                                     args = NULL) {
  if (is.null(args)) {
    args <- list(
      globalParameters = list(
        workflow_type = workflow_type,
        use_limpa = use_limpa,
        report_template = report_template
      ),
      deAnalysisParameters = list(formula_string = "~ 0 + group")
    )
  }

  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = protein_quant_table,
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "uniprot_acc",
    protein_id_table = data.frame(
      uniprot_acc = protein_quant_table$uniprot_acc,
      gene_names = paste0("GENE", seq_len(nrow(protein_quant_table))),
      stringsAsFactors = FALSE
    ),
    args = args
  )
}

module_ci_prot_da_contrasts <- function(kind = "two_group") {
  switch(kind,
    two_group = data.frame(
      contrasts = "groupB-groupA",
      full_format = "B_vs_A=groupB-groupA",
      friendly_names = "B_vs_A",
      stringsAsFactors = FALSE
    ),
    multi_group = data.frame(
      contrasts = c("groupB-groupA", "groupC-groupA"),
      full_format = c("B_vs_A=groupB-groupA", "C_vs_A=groupC-groupA"),
      friendly_names = c("B_vs_A", "C_vs_A"),
      stringsAsFactors = FALSE
    ),
    batch_adjusted = data.frame(
      contrasts = "groupB-groupA",
      full_format = "B_vs_A=groupB-groupA",
      friendly_names = "B_vs_A",
      stringsAsFactors = FALSE
    ),
    reversed = data.frame(
      contrasts = "groupA-groupB",
      full_format = "A_vs_B=groupA-groupB",
      friendly_names = "A_vs_B",
      stringsAsFactors = FALSE
    ),
    duplicate = data.frame(
      contrasts = c("groupB-groupA", "groupB-groupA"),
      full_format = c("B_vs_A=groupB-groupA", "B_vs_A_copy=groupB-groupA"),
      friendly_names = c("B_vs_A", "B_vs_A_copy"),
      stringsAsFactors = FALSE
    ),
    empty = data.frame(
      contrasts = "",
      full_format = "",
      friendly_names = "",
      stringsAsFactors = FALSE
    ),
    invalid_term = data.frame(
      contrasts = "groupZ-groupA",
      full_format = "Z_vs_A=groupZ-groupA",
      friendly_names = "Z_vs_A",
      stringsAsFactors = FALSE
    ),
    raw_only = data.frame(
      contrasts = "groupB-groupA",
      stringsAsFactors = FALSE
    ),
    stop(sprintf("Unknown proteomics DA contrast fixture kind: %s", kind), call. = FALSE)
  )
}

module_ci_prot_da_session_payload <- function(workflow_type = "DIA",
                                              use_limpa = FALSE,
                                              report_template = paste0(workflow_type, "_report")) {
  object <- module_ci_prot_da_object(
    workflow_type = workflow_type,
    use_limpa = use_limpa,
    report_template = report_template
  )
  list(
    export_timestamp = as.POSIXct("2026-05-05 08:00:00", tz = "UTC"),
    r6_current_state_name = "correlation_filtered",
    final_protein_count = nrow(object@protein_quant_table),
    final_sample_count = length(setdiff(names(object@protein_quant_table), object@protein_id_column)),
    r6_complete_states = list(raw_data_s4 = object, correlation_filtered = object),
    r6_state_history = c("raw_data_s4", "correlation_filtered"),
    contrasts_tbl = module_ci_prot_da_contrasts("two_group"),
    design_matrix = object@design_matrix,
    config_list = list(globalParameters = list(workflow_type = workflow_type, use_limpa = use_limpa)),
    current_s4_object = object,
    fasta_metadata = list(fasta_format = "UniProt", num_sequences = nrow(object@protein_id_table)),
    accession_cleanup_results = list(cleanup_applied = TRUE, aggregation_method = "sum"),
    ruv_optimization_result = list(best_k = if (isTRUE(use_limpa)) 3L else 2L, best_percentage = 10),
    qc_params = list(protein_qc = list(minimum = 1)),
    mixed_species_analysis = list(enabled = FALSE)
  )
}

module_ci_prot_da_write_session <- function(source_dir, session_data = module_ci_prot_da_session_payload()) {
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(session_data, file.path(source_dir, "filtered_session_data_latest.rds"))
  saveRDS(
    data.frame(
      Entry = session_data$current_s4_object@protein_id_table$uniprot_acc,
      gene_names = session_data$current_s4_object@protein_id_table$gene_names,
      stringsAsFactors = FALSE
    ),
    file.path(source_dir, "uniprot_dat_cln.RDS")
  )
  invisible(source_dir)
}

module_ci_prot_da_state_manager <- function(object = module_ci_prot_da_object(), fail_save = FALSE) {
  manager <- new.env(parent = emptyenv())
  manager$states <- list(correlation_filtered = list(config = list(source = "module-ci"), description = "filtered"))
  manager$state_history <- "correlation_filtered"
  manager$current_state <- "correlation_filtered"
  manager$getHistory <- function() manager$state_history
  manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    if (isTRUE(fail_save)) {
      stop("state save failed", call. = FALSE)
    }
    manager$states[[state_name]] <- list(
      object = s4_data_object,
      config = config_object,
      description = description
    )
    manager$state_history <- unique(c(manager$state_history, state_name))
    manager$current_state <- state_name
    invisible(TRUE)
  }
  manager
}

module_ci_prot_da_workflow <- function(object = module_ci_prot_da_object()) {
  shiny::reactiveValues(
    state_manager = module_ci_prot_da_state_manager(object),
    tab_status = list(
      normalization = "pending",
      differential_expression = "locked",
      differential_abundance = "pending",
      enrichment_analysis = "locked"
    ),
    da_analysis_results_list = NULL,
    da_ui_params = NULL
  )
}

module_ci_prot_da_paths <- function() {
  root <- tempfile("module-ci-prot-da-")
  list(
    source_dir = file.path(root, "source"),
    da_output_dir = file.path(root, "da_output"),
    publication_graphs_dir = file.path(root, "publication_graphs")
  )
}

module_ci_prot_da_result_table <- function(variant = "standard", comparison = "B_vs_A") {
  base <- data.frame(
    uniprot_acc = paste0("P", 1:6),
    gene_names = paste0("GENE", 1:6),
    comparison = comparison,
    numerator = "B",
    denominator = "A",
    log2FC = c(1.4, -1.2, 0.2, 2.1, -0.3, 0),
    raw_pvalue = c(0.001, 0.002, 0.5, 0.01, 0.8, 0.04),
    fdr_qvalue = c(0.01, 0.02, 0.6, 0.04, 0.8, 0.049),
    fdr_value_bh_adjustment = c(0.012, 0.024, 0.6, 0.05, 0.8, 0.06),
    log2norm.S1.A = c(10, 12, 24, 15, 8, 30),
    log2norm.S2.A = c(11, 12, 23, 15, 9, 31),
    log2norm.S3.B = c(22, 21, 11, 15, 10, 34),
    log2norm.S4.B = c(23, 22, 10, 15, 11, 35),
    raw.S1.A = c(100, 120, 240, 150, 80, 300),
    raw.S2.A = c(110, 120, 230, 150, 90, 310),
    raw.S3.B = c(220, 210, 110, 150, 100, 340),
    raw.S4.B = c(230, 220, 100, 150, 110, 350),
    stringsAsFactors = FALSE
  )

  switch(variant,
    standard = base,
    no_significant = transform(base, fdr_qvalue = c(0.2, 0.3, 0.6, 0.8, 0.9, 1)),
    all_significant = transform(base, fdr_qvalue = c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006)),
    tied_pvalues = transform(base, raw_pvalue = c(0.01, 0.01, 0.01, 0.2, 0.2, 0.2), fdr_qvalue = c(0.02, 0.02, 0.02, 0.25, 0.25, 0.25)),
    missing_values = transform(base, raw_pvalue = c(0.01, NA, Inf, NaN, 0.2, 0.3), fdr_qvalue = c(0.02, NA, NA, NA, 0.3, 0.4), log2FC = c(1.1, NA, -1.1, 0, Inf, NaN)),
    small_n = base[1:2, , drop = FALSE],
    stop(sprintf("Unknown proteomics DA result fixture variant: %s", variant), call. = FALSE)
  )
}

module_ci_prot_da_assert_result_schema <- function(table) {
  testthat::expect_s3_class(table, "data.frame")
  testthat::expect_true(all(c(
    "uniprot_acc", "comparison", "log2FC", "raw_pvalue", "fdr_qvalue",
    "fdr_value_bh_adjustment", "numerator", "denominator"
  ) %in% names(table)))
  testthat::expect_true(any(grepl("^log2norm\\.", names(table))))
  testthat::expect_true(any(grepl("^raw\\.", names(table))))
  result_keys <- paste(table$uniprot_acc, table$comparison, sep = "\r")
  testthat::expect_false(anyDuplicated(result_keys) > 0L)
  invisible(TRUE)
}

module_ci_prot_da_result_list <- function(variant = "standard", comparison = "B_vs_A") {
  object <- module_ci_prot_da_object()
  list(
    theObject = object,
    da_proteins_long = module_ci_prot_da_result_table(variant, comparison),
    num_sig_da_molecules_first_go = list(
      table = countStatDaGenes(
        module_ci_prot_da_result_table(variant, comparison),
        lfc_thresh = 0,
        q_val_thresh = 0.05,
        log_fc_column = log2FC,
        q_value_column = fdr_qvalue
      ) |>
        dplyr::mutate(comparison = comparison, analysis_type = "module-ci")
    )
  )
}

module_ci_prot_da_load_module <- function(da_data, workflow_data, experiment_paths) {
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

module_ci_prot_da_run_module <- function(da_data, workflow_data, experiment_paths) {
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

module_ci_prot_da_with_globals <- function(values = list(), code) {
  global_names <- unique(c("contrasts_tbl", "config_list", "uniprot_dat_cln", names(values)))
  old <- lapply(global_names, function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      list(exists = TRUE, value = get(name, envir = .GlobalEnv, inherits = FALSE))
    } else {
      list(exists = FALSE, value = NULL)
    }
  })
  names(old) <- global_names

  for (name in global_names) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }
  for (name in names(values)) {
    assign(name, values[[name]], envir = .GlobalEnv)
  }

  on.exit({
    for (name in global_names) {
      if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
      if (isTRUE(old[[name]]$exists)) {
        assign(name, old[[name]]$value, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  force(code)
}

module_ci_prot_da_with_namespace_bindings <- function(replacements, code) {
  package_env <- asNamespace("MultiScholaR")
  imports_env <- parent.env(package_env)
  names_vec <- names(replacements)
  target_envs <- lapply(names_vec, function(name) {
    if (exists(name, envir = package_env, inherits = FALSE)) {
      package_env
    } else {
      imports_env
    }
  })
  names(target_envs) <- names_vec
  existed <- vapply(names_vec, function(name) {
    exists(name, envir = target_envs[[name]], inherits = FALSE)
  }, logical(1))
  old_values <- lapply(names_vec, function(name) {
    if (existed[[name]]) {
      get(name, envir = target_envs[[name]], inherits = FALSE)
    } else {
      NULL
    }
  })
  names(old_values) <- names_vec
  locked <- vapply(names_vec, function(name) {
    existed[[name]] && bindingIsLocked(name, target_envs[[name]])
  }, logical(1))

  for (name in names_vec) {
    target_env <- target_envs[[name]]
    if (existed[[name]] && bindingIsLocked(name, target_env)) {
      unlockBinding(name, target_env)
    }
    assign(name, replacements[[name]], envir = target_env)
  }

  on.exit({
    for (name in names_vec) {
      target_env <- target_envs[[name]]
      if (exists(name, envir = target_env, inherits = FALSE) && bindingIsLocked(name, target_env)) {
        unlockBinding(name, target_env)
      }
      if (existed[[name]]) {
        assign(name, old_values[[name]], envir = target_env)
      } else if (exists(name, envir = target_env, inherits = FALSE)) {
        rm(list = name, envir = target_env)
      }
      if (locked[[name]]) {
        lockBinding(name, target_env)
      }
    }
  }, add = TRUE)

  force(code)
}

module_ci_prot_da_shiny_mocks <- function(captured = new.env(parent = emptyenv()),
                                          .local_envir = parent.frame()) {
  captured$notifications <- list()
  captured$text_updates <- list()
  captured$select_updates <- list()

  testthat::local_mocked_bindings(
    withProgress = function(message, value, expr) force(expr),
    incProgress = function(amount, detail = NULL) invisible(NULL),
    showNotification = function(message, type = NULL, duration = NULL, id = NULL, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration,
        id = id
      )
      invisible(NULL)
    },
    removeNotification = function(id) invisible(NULL),
    updateTextAreaInput = function(session, inputId, value = NULL, ...) {
      captured$text_updates[[inputId]] <<- value
      invisible(NULL)
    },
    updateSelectInput = function(session, inputId, choices = NULL, ...) {
      captured$select_updates[[inputId]] <<- choices
      invisible(NULL)
    },
    showModal = function(ui, ...) {
      captured$modal <<- ui
      invisible(NULL)
    },
    removeModal = function(...) invisible(NULL),
    .package = "shiny",
    .env = .local_envir
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger",
    .env = .local_envir
  )
  captured
}

module_ci_prot_da_render_mocks <- function(captured = new.env(parent = emptyenv()),
                                           .local_envir = parent.frame()) {
  captured$notifications <- list()
  captured$datatable_calls <- list()

  testthat::local_mocked_bindings(
    observeEvent = function(eventExpr, handlerExpr, ..., ignoreInit = FALSE) {
      value <- eval(substitute(eventExpr), parent.frame())
      triggered <- if (is.list(value)) {
        all(vapply(value, function(item) !is.null(item), logical(1)))
      } else {
        !is.null(value) && !identical(value, FALSE)
      }
      if (triggered) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value) || identical(value, FALSE)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1L]])
    },
    renderPlot = function(expr) eval(substitute(expr), parent.frame()),
    renderPrint = function(expr) paste(capture.output(eval(substitute(expr), parent.frame())), collapse = "\n"),
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    renderUI = function(expr) eval(substitute(expr), parent.frame()),
    showNotification = function(message, type = NULL, duration = NULL, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    downloadHandler = function(filename, content, ...) {
      list(filename = filename, content = content)
    },
    .package = "shiny",
    .env = .local_envir
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr) eval(substitute(expr), parent.frame()),
    datatable = function(data, options = NULL, extensions = NULL, ...) {
      captured$datatable_calls[[length(captured$datatable_calls) + 1L]] <<- list(
        data = data,
        options = options,
        extensions = extensions
      )
      structure(list(data = data, options = options, extensions = extensions), class = "module_ci_datatable")
    },
    formatRound = function(table, columns, digits, ...) {
      table$format_round <- list(columns = columns, digits = digits)
      table
    },
    .package = "DT",
    .env = .local_envir
  )
  captured
}

module_ci_prot_da_render_state <- function(table = module_ci_prot_da_result_table("standard")) {
  da_data <- new.env(parent = emptyenv())
  da_data$current_s4_object <- module_ci_prot_da_object()
  da_data$da_results_list <- list(
    da_proteins_long = table,
    theObject = da_data$current_s4_object
  )
  da_data$current_row_clusters <- NULL
  da_data$current_col_clusters <- NULL
  da_data$current_heatmap_plot <- NULL
  da_data
}
