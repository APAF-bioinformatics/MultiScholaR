library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "mod_metab_design_import_helpers.R"),
    file.path(repo_root, "R", "mod_metab_design.R")
  ),
  symbols = c(
    "registerMetabDesignPreviewOutputs",
    "registerMetabDesignStateOutputs",
    "registerMetabDesignBuilderModule",
    "registerMetabDesignBuilderResultsObserver",
    "initializeMetabDesignImportBootstrap",
    "registerMetabDesignImportModalShell",
    "resolveMetabDesignImportPreflight",
    "hydrateMetabDesignImportArtifacts",
    "hydrateMetabDesignImportMetadata",
    "hydrateMetabDesignImportWorkflowState",
    "createMetabDesignImportedS4Object",
    "saveMetabDesignImportedS4State",
    "initializeMetabDesignImportedQcBaseline",
    "completeMetabDesignImportedPostCheckpoint",
    "registerMetabDesignImportObserverShell",
    "mod_metab_design_server"
  ),
  env = environment()
)

test_that("metabolomics design preview seam preserves saved preview registration", {
  workflow_state <- list(
    design_matrix = data.frame(
      Run = c("Sample1", "Sample2"),
      group = c("Control", "Treatment"),
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    data_tbl = list(
      LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE),
      GCMS = data.frame(feature = "M2", stringsAsFactors = FALSE)
    )
  )

  output <- new.env(parent = emptyenv())

  fakeReq <- function(...) {
    values <- list(...)
    for (value in values) {
      if (is.null(value)) {
        stop("required value missing", call. = FALSE)
      }
    }
    invisible(values[[1]])
  }

  fakeRenderDt <- function(expr, options) {
    list(
      kind = "dt",
      value = eval(substitute(expr), parent.frame()),
      options = options
    )
  }

  fakeRenderText <- function(expr) {
    list(
      kind = "text",
      value = eval(substitute(expr), parent.frame())
    )
  }

  registerMetabDesignPreviewOutputs(
    output = output,
    workflowData = workflow_state,
    renderDt = fakeRenderDt,
    renderText = fakeRenderText,
    req = fakeReq
  )

  expect_identical(output$design_matrix_preview$kind, "dt")
  expect_identical(output$design_matrix_preview$value, workflow_state$design_matrix)
  expect_identical(
    output$design_matrix_preview$options,
    list(pageLength = 5, scrollX = TRUE)
  )

  expect_identical(output$contrasts_preview$kind, "dt")
  expect_identical(output$contrasts_preview$value, workflow_state$contrasts_tbl)
  expect_identical(
    output$contrasts_preview$options,
    list(pageLength = 5, scrollX = TRUE)
  )

  expect_identical(output$assays_preview$kind, "text")
  expect_identical(
    output$assays_preview$value,
    "Included assays: LCMS_Pos, GCMS"
  )
})

test_that("metabolomics design state seam preserves reactive output registration", {
  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    )
  )

  output <- new.env(parent = emptyenv())

  fakeReactive <- function(expr) {
    list(
      kind = "reactive",
      value = eval(substitute(expr), parent.frame())
    )
  }

  fakeOutputOptions <- function(output, name, suspendWhenHidden) {
    output[[paste0(name, "_options")]] <- list(
      suspendWhenHidden = suspendWhenHidden
    )
    invisible(output)
  }

  registerMetabDesignStateOutputs(
    output = output,
    workflowData = workflow_state,
    reactive = fakeReactive,
    outputOptions = fakeOutputOptions
  )

  expect_identical(output$data_available$kind, "reactive")
  expect_true(output$data_available$value)
  expect_identical(
    output$data_available_options,
    list(suspendWhenHidden = FALSE)
  )

  expect_identical(output$design_matrix_exists$kind, "reactive")
  expect_true(output$design_matrix_exists$value)
  expect_identical(
    output$design_matrix_exists_options,
    list(suspendWhenHidden = FALSE)
  )
})

test_that("metabolomics design builder seam preserves builder module registration", {
  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = c("Sample1", "Sample2")),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    )
  )

  fakeReactive <- function(value) {
    list(kind = "reactive", value = value)
  }

  captured_call <- NULL
  fakeBuilderServer <- function(
      id,
      data_tbl,
      config_list,
      column_mapping,
      existing_design_matrix,
      existing_contrasts
  ) {
    captured_call <<- list(
      id = id,
      data_tbl = data_tbl,
      config_list = config_list,
      column_mapping = column_mapping,
      existing_design_matrix = existing_design_matrix,
      existing_contrasts = existing_contrasts
    )

    "builder-results"
  }

  result <- registerMetabDesignBuilderModule(
    workflowData = workflow_state,
    builderServerExists = TRUE,
    builderServerFn = fakeBuilderServer,
    reactiveFn = fakeReactive
  )

  expect_identical(result, "builder-results")
  expect_identical(captured_call$id, "builder")
  expect_identical(captured_call$data_tbl, fakeReactive(workflow_state$data_tbl))
  expect_identical(captured_call$config_list, fakeReactive(workflow_state$config_list))
  expect_identical(captured_call$column_mapping, fakeReactive(workflow_state$column_mapping))
  expect_identical(
    captured_call$existing_design_matrix,
    fakeReactive(workflow_state$design_matrix)
  )
  expect_identical(
    captured_call$existing_contrasts,
    fakeReactive(workflow_state$contrasts_tbl)
  )
})

test_that("metabolomics design builder seam preserves missing-builder fallback", {
  fallback_values <- list()
  fakeReactiveVal <- function(value) {
    fallback_values <<- c(fallback_values, list(value))
    list(kind = "reactiveVal", value = value)
  }

  result <- registerMetabDesignBuilderModule(
    workflowData = list(),
    builderServerExists = FALSE,
    reactiveValFn = fakeReactiveVal
  )

  expect_identical(result, list(kind = "reactiveVal", value = NULL))
  expect_length(fallback_values, 1)
  expect_null(fallback_values[[1]])
})

test_that("metabolomics design builder-results seam preserves observer delegation", {
  builder_results <- list(
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    data_cln = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    config_list = list(normalization = "none")
  )

  builderResultsRv <- function() {
    builder_results
  }

  observed <- list()
  fakeObserveEvent <- function(eventExpr, handlerExpr, ignoreNULL) {
    observed$event <<- eval(substitute(eventExpr), parent.frame())
    observed$ignoreNULL <<- ignoreNULL
    eval(substitute(handlerExpr), parent.frame())
    "observer-registered"
  }

  req_calls <- list()
  fakeReq <- function(value) {
    req_calls <<- c(req_calls, list(value))
    if (is.null(value)) {
      stop("required value missing", call. = FALSE)
    }
    invisible(value)
  }

  delegated_call <- NULL
  fakeRunFlow <- function(results, workflowData, experimentPaths, qcTrigger) {
    delegated_call <<- list(
      results = results,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      qcTrigger = qcTrigger
    )

    "delegated"
  }

  workflow_state <- list(tab_status = list(design_matrix = "pending"))
  qc_trigger <- function(value) value
  experiment_paths <- list(source_dir = tempfile("metab-design-source-"))

  result <- registerMetabDesignBuilderResultsObserver(
    builderResultsRv = builderResultsRv,
    workflowData = workflow_state,
    experimentPaths = experiment_paths,
    qcTrigger = qc_trigger,
    observeEventFn = fakeObserveEvent,
    reqFn = fakeReq,
    runBuilderResultsFlow = fakeRunFlow
  )

  expect_identical(result, "observer-registered")
  expect_identical(observed$event, builder_results)
  expect_true(observed$ignoreNULL)
  expect_length(req_calls, 1)
  expect_identical(req_calls[[1]], builder_results)
  expect_identical(delegated_call$results, builder_results)
  expect_identical(delegated_call$workflowData, workflow_state)
  expect_identical(delegated_call$experimentPaths, experiment_paths)
  expect_identical(delegated_call$qcTrigger, qc_trigger)
})

test_that("metabolomics design import bootstrap seam preserves volume resolution and chooser registration", {
  captured <- new.env(parent = emptyenv())

  input <- list(import_dir = "dir-token")
  session <- list(ns = function(id) paste0("metab-", id))
  experiment_paths <- list(base_dir = tempdir())

  result <- initializeMetabDesignImportBootstrap(
    input = input,
    session = session,
    experimentPaths = experiment_paths,
    volumes = c(Home = "/project/home"),
    isolateFn = function(expr) eval(substitute(expr), parent.frame()),
    dirExistsFn = function(path) identical(path, experiment_paths$base_dir),
    dirChooseFn = function(input, id, roots, session) {
      captured$dir_choose <- list(
        input = input,
        id = id,
        roots = roots,
        session = session
      )
      invisible(NULL)
    },
    logInfo = function(message) {
      captured$log_message <- message
      invisible(NULL)
    }
  )

  expect_identical(
    unname(result),
    c(experiment_paths$base_dir, "/project/home")
  )
  expect_identical(names(result), c("Project Base Dir", "Home"))
  expect_identical(captured$dir_choose$id, "import_dir")
  expect_identical(captured$dir_choose$input, input)
  expect_identical(captured$dir_choose$session, session)
  expect_identical(captured$dir_choose$roots, result)
  expect_identical(
    captured$log_message,
    paste("Added base_dir to volumes:", experiment_paths$base_dir)
  )
})

test_that("metabolomics design import modal-shell seam preserves modal registration and path rendering", {
  captured <- new.env(parent = emptyenv())
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token"
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  resolved_volumes <- c("Project Base Dir" = tempdir(), Home = "/project/home")

  result <- registerMetabDesignImportModalShell(
    input = input,
    output = output,
    session = session,
    resolvedVolumes = resolved_volumes,
    observeEventFn = function(eventExpr, handlerExpr, ...) {
      captured$event <- eval(substitute(eventExpr), parent.frame())
      eval(substitute(handlerExpr), parent.frame())
      "observer-registered"
    },
    showModalFn = function(x) {
      captured$modal <- x
      invisible(NULL)
    },
    modalDialogFn = function(...) list(...),
    paragraphFn = function(...) paste(...),
    helpTextFn = function(...) paste(...),
    dirButtonFn = function(id, label, title) {
      list(kind = "dirButton", id = id, label = label, title = title)
    },
    verbatimTextOutputFn = function(id, placeholder = FALSE) {
      list(kind = "verbatim", id = id, placeholder = placeholder)
    },
    tagListFn = function(...) list(...),
    modalButtonFn = function(label) list(kind = "modalButton", label = label),
    actionButtonFn = function(inputId, label, class = NULL, icon = NULL) {
      list(kind = "actionButton", inputId = inputId, label = label, class = class)
    },
    renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
    reqFn = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    parseDirPathFn = function(roots, selection) {
      captured$parse <- list(roots = roots, selection = selection)
      "/tmp/imported-design"
    }
  )

  expect_identical(result, output)
  expect_identical(captured$event, 1)
  expect_identical(captured$parse$roots, resolved_volumes)
  expect_identical(captured$parse$selection, "dir-token")
  expect_identical(output$import_dir_path, "/tmp/imported-design")
  expect_identical(captured$modal$title, "Import Existing Design Matrix")
  expect_identical(captured$modal[[4]]$id, "metab-import_dir")
  expect_identical(captured$modal[[5]]$id, "metab-import_dir_path")
  expect_identical(captured$modal$footer[[2]]$inputId, "metab-confirm_import")
})

test_that("metabolomics design import preflight seam preserves path resolution and required-file validation", {
  captured <- new.env(parent = emptyenv())
  input <- list(import_dir = "dir-token")
  resolved_volumes <- c("Project Base Dir" = tempdir(), Home = "/project/home")
  req_calls <- list()

  result <- resolveMetabDesignImportPreflight(
    input = input,
    resolvedVolumes = resolved_volumes,
    reqFn = function(value) {
      req_calls <<- c(req_calls, list(value))
      if (is.null(value)) {
        stop("required value missing", call. = FALSE)
      }
      invisible(value)
    },
    parseDirPathFn = function(roots, selection) {
      captured$parse <- list(roots = roots, selection = selection)
      "/tmp/imported-design"
    },
    fileExistsFn = function(path) {
      captured$file_checks <- c(captured$file_checks, path)
      TRUE
    }
  )

  expect_true(result$ok)
  expect_identical(req_calls, list("dir-token", "/tmp/imported-design"))
  expect_identical(captured$parse$roots, resolved_volumes)
  expect_identical(captured$parse$selection, "dir-token")
  expect_identical(result$importPath, "/tmp/imported-design")
  expect_identical(result$designFile, "/tmp/imported-design/design_matrix.tab")
  expect_identical(result$contrastFile, "/tmp/imported-design/contrast_strings.tab")
  expect_identical(result$manifestFile, "/tmp/imported-design/assay_manifest.txt")
  expect_identical(result$colMapFile, "/tmp/imported-design/column_mapping.json")
  expect_identical(result$configFile, "/tmp/imported-design/config.ini")
  expect_identical(
    captured$file_checks,
    c(
      "/tmp/imported-design/design_matrix.tab",
      "/tmp/imported-design/assay_manifest.txt"
    )
  )
})

test_that("metabolomics design import preflight seam preserves missing-file failure contract", {
  input <- list(import_dir = "dir-token")
  resolved_volumes <- c("Project Base Dir" = tempdir(), Home = "/project/home")

  missing_design <- resolveMetabDesignImportPreflight(
    input = input,
    resolvedVolumes = resolved_volumes,
    reqFn = function(value) invisible(value),
    parseDirPathFn = function(...) "/tmp/imported-design",
    fileExistsFn = function(path) !grepl("design_matrix\\.tab$", path)
  )

  expect_false(missing_design$ok)
  expect_identical(
    missing_design$errorMessage,
    "Import failed: 'design_matrix.tab' not found in the selected directory."
  )

  missing_manifest <- resolveMetabDesignImportPreflight(
    input = input,
    resolvedVolumes = resolved_volumes,
    reqFn = function(value) invisible(value),
    parseDirPathFn = function(...) "/tmp/imported-design",
    fileExistsFn = function(path) !grepl("assay_manifest\\.txt$", path)
  )

  expect_false(missing_manifest$ok)
  expect_identical(
    missing_manifest$errorMessage,
    "Import failed: 'assay_manifest.txt' not found. Cannot determine which assay files to load."
  )
})

test_that("metabolomics design import hydration seam preserves config and artifact loading", {
  captured <- new.env(parent = emptyenv())
  workflow_state <- new.env(parent = emptyenv())
  workflow_state$config_list <- NULL
  experiment_paths <- list(source_dir = "/tmp/source-dir")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    stringsAsFactors = FALSE
  )
  assay_tables <- list(
    LCMS_Pos = data.frame(feature = "M1", Sample1 = 1, stringsAsFactors = FALSE),
    GCMS = data.frame(feature = "M2", Sample1 = 2, stringsAsFactors = FALSE)
  )

  result <- hydrateMetabDesignImportArtifacts(
    workflowData = workflow_state,
    experimentPaths = experiment_paths,
    importPath = "/tmp/imported-design",
    designFile = "/tmp/imported-design/design_matrix.tab",
    manifestFile = "/tmp/imported-design/assay_manifest.txt",
    configFile = "/tmp/imported-design/config.ini",
    readConfigFn = function(file) {
      captured$config_read <- file
      list(globalParameters = list(workflow_type = "metabolomics"))
    },
    readTabularFn = function(file, show_col_types = FALSE) {
      captured$tabular_reads <- c(captured$tabular_reads, file)
      expect_false(show_col_types)

      if (grepl("design_matrix\\.tab$", file)) {
        imported_design
      } else if (grepl("data_cln_LCMS_Pos\\.tab$", file)) {
        assay_tables$LCMS_Pos
      } else if (grepl("data_cln_GCMS\\.tab$", file)) {
        assay_tables$GCMS
      } else {
        stop(sprintf("unexpected read path: %s", file))
      }
    },
    readLinesFn = function(file) {
      captured$manifest_read <- file
      c("LCMS_Pos", "", "GCMS")
    },
    fileExistsFn = function(path) {
      captured$file_checks <- c(captured$file_checks, path)
      TRUE
    },
    assignFn = function(x, value, envir) {
      captured$assign <- list(name = x, value = value, envir = envir)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured$log_warn <- c(captured$log_warn, message)
      invisible(NULL)
    },
    messageFn = function(...) {
      captured$messages <- c(captured$messages, paste(..., collapse = " "))
      invisible(NULL)
    }
  )

  expect_identical(
    workflow_state$config_list,
    list(globalParameters = list(workflow_type = "metabolomics"))
  )
  expect_identical(captured$config_read, "/tmp/imported-design/config.ini")
  expect_identical(captured$manifest_read, "/tmp/imported-design/assay_manifest.txt")
  expect_identical(captured$assign$name, "config_list")
  expect_identical(captured$assign$value, workflow_state$config_list)
  expect_identical(captured$assign$envir, .GlobalEnv)
  expect_identical(result$importedDesign, imported_design)
  expect_identical(result$assayNames, c("LCMS_Pos", "GCMS"))
  expect_identical(names(result$assayList), c("LCMS_Pos", "GCMS"))
  expect_identical(result$assayList$LCMS_Pos, assay_tables$LCMS_Pos)
  expect_identical(result$assayList$GCMS, assay_tables$GCMS)
  expect_identical(
    captured$tabular_reads,
    c(
      "/tmp/imported-design/design_matrix.tab",
      "/tmp/imported-design/data_cln_LCMS_Pos.tab",
      "/tmp/imported-design/data_cln_GCMS.tab"
    )
  )
})

test_that("metabolomics design import hydration seam preserves missing-assay failure contract", {
  workflow_state <- new.env(parent = emptyenv())
  workflow_state$config_list <- NULL
  experiment_paths <- list(source_dir = "/tmp/source-dir")

  expect_error(
    hydrateMetabDesignImportArtifacts(
      workflowData = workflow_state,
      experimentPaths = experiment_paths,
      importPath = "/tmp/imported-design",
      designFile = "/tmp/imported-design/design_matrix.tab",
      manifestFile = "/tmp/imported-design/assay_manifest.txt",
      configFile = "/tmp/imported-design/config.ini",
      readConfigFn = function(file) list(file = file),
      readTabularFn = function(file, show_col_types = FALSE) {
        data.frame(Run = "Sample1", stringsAsFactors = FALSE)
      },
      readLinesFn = function(file) "LCMS_Pos",
      fileExistsFn = function(path) !grepl("data_cln_LCMS_Pos\\.tab$", path),
      assignFn = function(...) invisible(NULL),
      logInfo = function(...) invisible(NULL),
      logWarn = function(...) invisible(NULL),
      messageFn = function(...) invisible(NULL)
    ),
    "Assay data file not found: /tmp/imported-design/data_cln_LCMS_Pos.tab"
  )
})

test_that("metabolomics design import metadata seam infers column mapping and loads contrasts", {
  captured <- new.env(parent = emptyenv())
  captured$file_checks <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$messages <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$column_mapping <- NULL
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      "Peak ID" = "M1",
      "Metabolite name" = "Met1",
      Sample1 = 10,
      Sample2 = 20,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )

  result <- hydrateMetabDesignImportMetadata(
    workflowData = workflow_state,
    importedDesign = imported_design,
    assayList = assay_list,
    colMapFile = "/tmp/imported-design/column_mapping.json",
    contrastFile = "/tmp/imported-design/contrast_strings.tab",
    fileExistsFn = function(path) {
      captured$file_checks <- c(captured$file_checks, path)
      grepl("contrast_strings\\.tab$", path)
    },
    readJsonFn = function(path) stop(sprintf("unexpected json read: %s", path)),
    readLinesFn = function(path) {
      captured$contrast_read <- path
      c("groupTreatment-groupControl", "groupDose-groupControl")
    },
    logInfo = function(message) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured$log_warn <- c(captured$log_warn, message)
      invisible(NULL)
    },
    messageFn = function(...) {
      captured$messages <- c(captured$messages, paste(..., collapse = " "))
      invisible(NULL)
    }
  )

  expect_identical(
    workflow_state$column_mapping,
    list(
      metabolite_id_col = "Peak ID",
      annotation_col = "Metabolite name",
      sample_columns = c("Sample1", "Sample2"),
      is_pattern = NA_character_
    )
  )
  expect_identical(result$columnMapping, workflow_state$column_mapping)
  expect_identical(captured$contrast_read, "/tmp/imported-design/contrast_strings.tab")
  expect_identical(
    result$importedContrasts,
    data.frame(
      contrasts = c("groupTreatment-groupControl", "groupDose-groupControl"),
      friendly_names = c("Treatment_vs_Control", "Dose_vs_Control"),
      full_format = c(
        "Treatment_vs_Control=groupTreatment-groupControl",
        "Dose_vs_Control=groupDose-groupControl"
      ),
      stringsAsFactors = FALSE
    )
  )
  expect_identical(
    captured$log_warn,
    "column_mapping.json not found. Inferring from data structure."
  )
  expect_identical(
    captured$file_checks,
    c(
      "/tmp/imported-design/column_mapping.json",
      "/tmp/imported-design/contrast_strings.tab"
    )
  )
})

test_that("metabolomics design import metadata seam loads column mapping json and preserves missing-contrast contract", {
  captured <- new.env(parent = emptyenv())
  captured$file_checks <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$messages <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$column_mapping <- NULL
  imported_design <- data.frame(Run = "Sample1", stringsAsFactors = FALSE)
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = "M1",
      Sample1 = 10,
      stringsAsFactors = FALSE
    )
  )
  expected_col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = "Sample1",
    is_pattern = "IS_"
  )

  result <- hydrateMetabDesignImportMetadata(
    workflowData = workflow_state,
    importedDesign = imported_design,
    assayList = assay_list,
    colMapFile = "/tmp/imported-design/column_mapping.json",
    contrastFile = "/tmp/imported-design/contrast_strings.tab",
    fileExistsFn = function(path) {
      captured$file_checks <- c(captured$file_checks, path)
      grepl("column_mapping\\.json$", path)
    },
    readJsonFn = function(path) {
      captured$json_read <- path
      expected_col_map
    },
    readLinesFn = function(path) stop(sprintf("unexpected contrast read: %s", path)),
    logInfo = function(message) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      captured$log_warn <- c(captured$log_warn, message)
      invisible(NULL)
    },
    messageFn = function(...) {
      captured$messages <- c(captured$messages, paste(..., collapse = " "))
      invisible(NULL)
    }
  )

  expect_identical(captured$json_read, "/tmp/imported-design/column_mapping.json")
  expect_identical(workflow_state$column_mapping, expected_col_map)
  expect_identical(result$columnMapping, expected_col_map)
  expect_null(result$importedContrasts)
  expect_identical(captured$log_warn, character())
  expect_identical(captured$log_info, c("Loading column_mapping.json", "No contrast_strings.tab found."))
  expect_true(any(grepl("Loaded column_mapping from JSON", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("No contrasts loaded", captured$messages, fixed = TRUE)))
})

test_that("metabolomics design import workflow-state seam hydrates workflow fields and tech-rep groups", {
  captured <- new.env(parent = emptyenv())
  captured$log_info <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$design_matrix <- NULL
  workflow_state$contrasts_tbl <- NULL
  workflow_state$data_tbl <- NULL
  workflow_state$data_cln <- NULL
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )

  result <- hydrateMetabDesignImportWorkflowState(
    workflowData = workflow_state,
    importedDesign = imported_design,
    assayList = assay_list,
    importedContrasts = imported_contrasts,
    assignFn = function(name, value, envir) {
      captured$assigned <- list(name = name, value = value, envir = envir)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    }
  )

  expect_identical(workflow_state$data_tbl, assay_list)
  expect_identical(workflow_state$data_cln, assay_list)
  expect_identical(workflow_state$contrasts_tbl, imported_contrasts)
  expect_identical(workflow_state$design_matrix$tech_rep_group, c("Control_R1", "Treatment_R2"))
  expect_identical(result$designMatrix, workflow_state$design_matrix)
  expect_identical(result$assayList, assay_list)
  expect_identical(result$importedContrasts, imported_contrasts)
  expect_identical(
    captured$assigned,
    list(name = "contrasts_tbl", value = imported_contrasts, envir = .GlobalEnv)
  )
  expect_identical(captured$log_info, "Saved contrasts_tbl to global environment.")
})

test_that("metabolomics design import S4 seam preserves argument assembly and invocation", {
  captured <- new.env(parent = emptyenv())
  captured$messages <- character()
  captured$log_info <- character()

  workflow_state <- list(
    design_matrix = data.frame(
      Run = c("Sample1", "Sample2"),
      group = c("Control", "Treatment"),
      tech_rep_group = c("Control_R1", "Treatment_R2"),
      stringsAsFactors = FALSE
    ),
    config_list = list(normalization = "none")
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )

  result <- createMetabDesignImportedS4Object(
    workflowData = workflow_state,
    assayList = assay_list,
    colMap = col_map,
    createS4Fn = function(...) {
      captured$create_args <- list(...)
      structure(list(ok = TRUE), class = "MetaboliteAssayData")
    },
    logInfo = function(message) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    messageFn = function(...) {
      captured$messages <- c(captured$messages, paste(..., collapse = " "))
      invisible(NULL)
    }
  )

  expect_s3_class(result, "MetaboliteAssayData")
  expect_identical(captured$create_args$metabolite_data, assay_list)
  expect_identical(captured$create_args$design_matrix, workflow_state$design_matrix)
  expect_identical(captured$create_args$metabolite_id_column, "Feature")
  expect_identical(captured$create_args$annotation_id_column, NA_character_)
  expect_identical(captured$create_args$sample_id, "Run")
  expect_identical(captured$create_args$group_id, "group")
  expect_identical(captured$create_args$technical_replicate_id, "tech_rep_group")
  expect_identical(captured$create_args$database_identifier_type, "Unknown")
  expect_identical(captured$create_args$internal_standard_regex, NA_character_)
  expect_identical(captured$create_args$args, workflow_state$config_list)
  expect_identical(captured$log_info, "Creating MetaboliteAssayData S4 object from imported data")
  expect_true(any(grepl("DEBUG66: assay_list names: LCMS_Pos", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("DEBUG66: Calling createMetaboliteAssayData()", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("DEBUG66: S4 object class: MetaboliteAssayData", captured$messages, fixed = TRUE)))
})

test_that("metabolomics design imported state-manager seam preserves initialization and saveState contract", {
  captured <- new.env(parent = emptyenv())
  captured$messages <- character()

  workflow_state <- new.env(parent = emptyenv())
  workflow_state$config_list <- list(normalization = "none")
  workflow_state$state_manager <- NULL

  fake_manager <- new.env(parent = emptyenv())
  fake_manager$saveState <- function(...) {
    captured$save_args <- list(...)
    invisible(NULL)
  }

  result <- saveMetabDesignImportedS4State(
    workflowData = workflow_state,
    s4Object = structure(list(ok = TRUE), class = "MetaboliteAssayData"),
    stateManagerFactory = function(omicsType) {
      captured$factory_arg <- omicsType
      fake_manager
    },
    messageFn = function(...) {
      captured$messages <- c(captured$messages, paste(..., collapse = " "))
      invisible(NULL)
    }
  )

  expect_identical(captured$factory_arg, "metabolomics")
  expect_identical(workflow_state$state_manager, fake_manager)
  expect_identical(result, fake_manager)
  expect_identical(captured$save_args$state_name, "metab_raw_data_s4")
  expect_s3_class(captured$save_args$s4_data_object, "MetaboliteAssayData")
  expect_identical(captured$save_args$config_object, workflow_state$config_list)
  expect_identical(
    captured$save_args$description,
    "MetaboliteAssayData S4 object created from imported design"
  )
  expect_identical(
    captured$messages,
    c(
      "DEBUG66: Initializing state manager...",
      "DEBUG66: Created new WorkflowState",
      "DEBUG66: Saving state to state_manager...",
      "DEBUG66: State saved successfully"
    )
  )
})

test_that("metabolomics design imported QC-baseline seam preserves raw-data updateFiltering contract", {
  success_capture <- new.env(parent = emptyenv())
  success_capture$log_info <- character()
  success_capture$log_warn <- character()
  fake_s4 <- structure(list(ok = TRUE), class = "MetaboliteAssayData")

  success_result <- initializeMetabDesignImportedQcBaseline(
    s4Object = fake_s4,
    updateFilteringFn = function(...) {
      success_capture$update_args <- list(...)
      invisible(NULL)
    },
    logInfo = function(message) {
      success_capture$log_info <- c(success_capture$log_info, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      success_capture$log_warn <- c(success_capture$log_warn, message)
      invisible(NULL)
    }
  )

  expect_true(success_result)
  expect_identical(success_capture$update_args$theObject, fake_s4)
  expect_identical(success_capture$update_args$step_name, "1_Raw_Data")
  expect_identical(success_capture$update_args$omics_type, "metabolomics")
  expect_false(success_capture$update_args$return_grid)
  expect_true(success_capture$update_args$overwrite)
  expect_identical(
    success_capture$log_info,
    "Initialized QC progress tracking with raw data baseline"
  )
  expect_identical(success_capture$log_warn, character())

  failure_capture <- new.env(parent = emptyenv())
  failure_capture$log_info <- character()
  failure_capture$log_warn <- character()

  failure_result <- initializeMetabDesignImportedQcBaseline(
    s4Object = fake_s4,
    updateFilteringFn = function(...) stop("qc update failed"),
    logInfo = function(message) {
      failure_capture$log_info <- c(failure_capture$log_info, message)
      invisible(NULL)
    },
    logWarn = function(message) {
      failure_capture$log_warn <- c(failure_capture$log_warn, message)
      invisible(NULL)
    }
  )

  expect_false(failure_result)
  expect_identical(failure_capture$log_info, character())
  expect_identical(
    failure_capture$log_warn,
    "Could not initialize QC progress: qc update failed"
  )
})

test_that("metabolomics design imported completion seam preserves qc-trigger, tab-status, and notification handoff", {
  capture <- new.env(parent = emptyenv())
  capture$log_info <- character()
  capture$notifications <- list()

  workflow_state <- list(
    tab_status = list(design_matrix = "pending", qc = "pending")
  )
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    stringsAsFactors = FALSE
  )

  fake_qc_trigger <- function(value) {
    capture$qc_trigger_value <- value
    invisible(NULL)
  }

  result <- completeMetabDesignImportedPostCheckpoint(
    workflowData = workflow_state,
    assayNames = c("LCMS_Pos", "GCMS"),
    importedDesign = imported_design,
    qcTrigger = fake_qc_trigger,
    logInfo = function(message) {
      capture$log_info <- c(capture$log_info, message)
      invisible(NULL)
    },
    removeNotification = function(id) {
      capture$removed_notification <- id
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      capture$notifications[[length(capture$notifications) + 1]] <- list(
        ui = ui,
        type = type,
        duration = duration,
        id = id
      )
      invisible(NULL)
    }
  )

  expect_identical(
    capture$log_info,
    c(
      "S4 object saved to state manager as 'metab_raw_data_s4'",
      "Import complete - user can proceed to QC"
    )
  )
  expect_true(isTRUE(capture$qc_trigger_value))
  expect_identical(workflow_state$tab_status$design_matrix, "pending")
  expect_identical(workflow_state$tab_status$qc, "pending")
  expect_identical(result$tabStatus$design_matrix, "complete")
  expect_identical(result$tabStatus$qc, "pending")
  expect_identical(result$successMessage, "Design imported successfully! Loaded 2 assays with 2 samples.")
  expect_identical(capture$removed_notification, "importing_design")
  expect_identical(
    capture$notifications,
    list(list(
      ui = "Design imported successfully! Loaded 2 assays with 2 samples.",
      type = "message",
      duration = NULL,
      id = NULL
    ))
  )

  null_capture <- new.env(parent = emptyenv())
  null_capture$notifications <- list()
  workflow_state_null <- list(
    tab_status = list(design_matrix = "pending")
  )

  completeMetabDesignImportedPostCheckpoint(
    workflowData = workflow_state_null,
    assayNames = "LCMS_Pos",
    importedDesign = imported_design[1, , drop = FALSE],
    qcTrigger = NULL,
    logInfo = function(message) {
      null_capture$log_info <- c(null_capture$log_info, message)
      invisible(NULL)
    },
    removeNotification = function(id) {
      null_capture$removed_notification <- id
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      null_capture$notifications[[length(null_capture$notifications) + 1]] <- list(
        ui = ui,
        type = type,
        duration = duration,
        id = id
      )
      invisible(NULL)
    }
  )

  expect_false(exists("qc_trigger_value", envir = null_capture, inherits = FALSE))
  expect_identical(workflow_state_null$tab_status$design_matrix, "pending")
  expect_identical(
    null_capture$notifications[[1]]$ui,
    "Design imported successfully! Loaded 1 assays with 1 samples."
  )
})

test_that("metabolomics design import observer-shell seam preserves confirm_import orchestration", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()

  input <- list(confirm_import = 1)
  resolved_volumes <- c("Project Base Dir" = tempdir(), Home = "/project/home")
  workflow_state <- list(tab_status = list(design_matrix = "pending"))
  experiment_paths <- list(source_dir = "/tmp/source-dir")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    stringsAsFactors = FALSE
  )
  assay_names <- c("LCMS_Pos", "GCMS")
  assay_list <- list(
    LCMS_Pos = data.frame(Feature = "M1", stringsAsFactors = FALSE),
    GCMS = data.frame(Feature = "M2", stringsAsFactors = FALSE)
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )
  fake_s4 <- structure(list(ok = TRUE), class = "MetaboliteAssayData")
  fake_qc_trigger <- function(value) {
    captured$qc_trigger_value <- value
    invisible(NULL)
  }

  result <- registerMetabDesignImportObserverShell(
    input = input,
    resolvedVolumes = resolved_volumes,
    workflowData = workflow_state,
    experimentPaths = experiment_paths,
    qcTrigger = fake_qc_trigger,
    observeEventFn = function(eventExpr, handlerExpr, ...) {
      captured$event <- eval(substitute(eventExpr), parent.frame())
      eval(substitute(handlerExpr), parent.frame())
      "observer-registered"
    },
    resolveImportPreflight = function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    hydrateImportArtifacts = function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$artifacts <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = assay_names,
        assayList = assay_list
      )
    },
    hydrateImportMetadata = function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      list(
        columnMapping = col_map,
        importedContrasts = imported_contrasts
      )
    },
    hydrateImportWorkflowState = function(workflowData, importedDesign, assayList, importedContrasts, ...) {
      captured$workflow_state <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        importedContrasts = importedContrasts
      )
      invisible(NULL)
    },
    createImportedS4Object = function(workflowData, assayList, colMap, ...) {
      captured$create_s4 <- list(
        workflowData = workflowData,
        assayList = assayList,
        colMap = colMap
      )
      fake_s4
    },
    saveImportedS4State = function(workflowData, s4Object, ...) {
      captured$save_state <- list(
        workflowData = workflowData,
        s4Object = s4Object
      )
      invisible("saved")
    },
    initializeImportedQcBaseline = function(s4Object, ...) {
      captured$qc_baseline <- list(s4Object = s4Object)
      invisible(TRUE)
    },
    completeImportedPostCheckpoint = function(workflowData, assayNames, importedDesign, qcTrigger, ...) {
      captured$completion <- list(
        workflowData = workflowData,
        assayNames = assayNames,
        importedDesign = importedDesign,
        qcTrigger = qcTrigger
      )
      invisible(list(done = TRUE))
    },
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notifications[[length(captured$notifications) + 1]] <- list(
        ui = ui,
        type = type,
        duration = duration,
        id = id
      )
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    logError = function(message) {
      captured$log_error <- message
      invisible(NULL)
    }
  )

  expect_identical(result, "observer-registered")
  expect_identical(captured$event, 1)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, resolved_volumes)
  expect_identical(captured$artifacts$workflowData, workflow_state)
  expect_identical(captured$artifacts$experimentPaths, experiment_paths)
  expect_identical(captured$artifacts$importPath, "/tmp/imported-design")
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$metadata$importedDesign, imported_design)
  expect_identical(captured$metadata$assayList, assay_list)
  expect_identical(captured$metadata$colMapFile, "/tmp/imported-design/column_mapping.json")
  expect_identical(captured$metadata$contrastFile, "/tmp/imported-design/contrast_strings.tab")
  expect_identical(captured$workflow_state$workflowData, workflow_state)
  expect_identical(captured$workflow_state$importedDesign, imported_design)
  expect_identical(captured$workflow_state$assayList, assay_list)
  expect_identical(captured$workflow_state$importedContrasts, imported_contrasts)
  expect_identical(captured$create_s4$workflowData, workflow_state)
  expect_identical(captured$create_s4$assayList, assay_list)
  expect_identical(captured$create_s4$colMap, col_map)
  expect_identical(captured$save_state$workflowData, workflow_state)
  expect_identical(captured$save_state$s4Object, fake_s4)
  expect_identical(captured$qc_baseline$s4Object, fake_s4)
  expect_identical(captured$completion$workflowData, workflow_state)
  expect_identical(captured$completion$assayNames, assay_names)
  expect_identical(captured$completion$importedDesign, imported_design)
  expect_identical(captured$completion$qcTrigger, fake_qc_trigger)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(
    captured$notifications,
    list(list(
      ui = "Importing design files...",
      type = NULL,
      duration = NULL,
      id = "importing_design"
    ))
  )
  expect_false(exists("remove_notification", envir = captured, inherits = FALSE))
  expect_false(exists("log_error", envir = captured, inherits = FALSE))
})

test_that("metabolomics design import observer-shell seam preserves shared error notification handoff", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()

  result <- registerMetabDesignImportObserverShell(
    input = list(confirm_import = 1),
    resolvedVolumes = c("Project Base Dir" = tempdir(), Home = "/project/home"),
    workflowData = list(tab_status = list(design_matrix = "pending")),
    experimentPaths = list(source_dir = "/tmp/source-dir"),
    observeEventFn = function(eventExpr, handlerExpr, ...) {
      captured$event <- eval(substitute(eventExpr), parent.frame())
      eval(substitute(handlerExpr), parent.frame())
      "observer-registered"
    },
    resolveImportPreflight = function(...) {
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    hydrateImportArtifacts = function(...) stop("artifact-stop"),
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notifications[[length(captured$notifications) + 1]] <- list(
        ui = ui,
        type = type,
        duration = duration,
        id = id
      )
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    logError = function(message) {
      captured$log_error <- message
      invisible(NULL)
    }
  )

  expect_identical(result, "observer-registered")
  expect_identical(captured$event, 1)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: artifact-stop")
  expect_identical(
    captured$notifications,
    list(
      list(
        ui = "Importing design files...",
        type = NULL,
        duration = NULL,
        id = "importing_design"
      ),
      list(
        ui = "Error during import: artifact-stop",
        type = "error",
        duration = 15,
        id = NULL
      )
    )
  )
})

test_that("metabolomics design server delegates import observer-shell seam after modal shell", {
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  fake_qc_trigger <- function(value) value

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_import_observer <- exists("registerMetabDesignImportObserverShell", envir = server_env, inherits = FALSE)
  if (had_import_observer) {
    original_import_observer <- get("registerMetabDesignImportObserverShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportObserverShell",
    function(input, resolvedVolumes, workflowData, experimentPaths, qcTrigger, ...) {
      captured$import_observer <- list(
        input = input,
        resolvedVolumes = resolvedVolumes,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        qcTrigger = qcTrigger
      )
      invisible("observer-registered")
    },
    envir = server_env
  )
  on.exit({
    if (had_import_observer) {
      assign("registerMetabDesignImportObserverShell", original_import_observer, envir = server_env)
    } else if (exists("registerMetabDesignImportObserverShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportObserverShell", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      invisible(NULL)
    },
    req = function(...) invisible(list(...)[[1]]),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home"),
    qc_trigger = fake_qc_trigger
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$import_observer$input, input)
  expect_identical(captured$import_observer$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$import_observer$workflowData, workflow_state)
  expect_identical(captured$import_observer$experimentPaths, experiment_paths)
  expect_identical(captured$import_observer$qcTrigger, fake_qc_trigger)
})

test_that("metabolomics design server delegates import bootstrap, modal-shell, and preflight seams", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = FALSE,
        errorMessage = "preflight-failed"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )
  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$bootstrap$session, session)
  expect_identical(captured$bootstrap$experimentPaths, experiment_paths)
  expect_identical(captured$bootstrap$volumes, c(Home = "/project/home"))
  expect_identical(captured$modal_shell$input, input)
  expect_identical(captured$modal_shell$output, output)
  expect_identical(captured$modal_shell$session, session)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$log_error, "preflight-failed")
  expect_identical(
    captured$notification,
    list(ui = "preflight-failed", type = "error", duration = 10, id = NULL)
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates imported QC-baseline seam after state-manager save", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    friendly_names = "Treatment_vs_Control",
    full_format = "Treatment_vs_Control=groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )
  fake_s4 <- structure(list(ok = TRUE), class = "MetaboliteAssayData")

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = "LCMS_Pos",
        assayList = assay_list
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)
  had_metadata <- exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  if (had_metadata) {
    original_metadata <- get("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportMetadata",
    function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      list(
        columnMapping = col_map,
        importedContrasts = imported_contrasts
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_metadata) {
      assign("hydrateMetabDesignImportMetadata", original_metadata, envir = server_env)
    } else if (exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportMetadata", envir = server_env)
    }
  }, add = TRUE)
  had_workflow_state <- exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  if (had_workflow_state) {
    original_workflow_state <- get("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportWorkflowState",
    function(workflowData, importedDesign, assayList, importedContrasts, ...) {
      captured$workflow_state <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        importedContrasts = importedContrasts
      )
      invisible(NULL)
    },
    envir = server_env
  )
  on.exit({
    if (had_workflow_state) {
      assign("hydrateMetabDesignImportWorkflowState", original_workflow_state, envir = server_env)
    } else if (exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportWorkflowState", envir = server_env)
    }
  }, add = TRUE)
  had_create_s4 <- exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  if (had_create_s4) {
    original_create_s4 <- get("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  }
  assign(
    "createMetabDesignImportedS4Object",
    function(workflowData, assayList, colMap, ...) {
      captured$create_s4 <- list(
        workflowData = workflowData,
        assayList = assayList,
        colMap = colMap
      )
      fake_s4
    },
    envir = server_env
  )
  on.exit({
    if (had_create_s4) {
      assign("createMetabDesignImportedS4Object", original_create_s4, envir = server_env)
    } else if (exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)) {
      rm("createMetabDesignImportedS4Object", envir = server_env)
    }
  }, add = TRUE)
  had_save_state <- exists("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)
  if (had_save_state) {
    original_save_state <- get("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)
  }
  assign(
    "saveMetabDesignImportedS4State",
    function(workflowData, s4Object, ...) {
      captured$save_state <- list(
        workflowData = workflowData,
        s4Object = s4Object
      )
      invisible("saved")
    },
    envir = server_env
  )
  on.exit({
    if (had_save_state) {
      assign("saveMetabDesignImportedS4State", original_save_state, envir = server_env)
    } else if (exists("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)) {
      rm("saveMetabDesignImportedS4State", envir = server_env)
    }
  }, add = TRUE)
  had_qc_baseline <- exists("initializeMetabDesignImportedQcBaseline", envir = server_env, inherits = FALSE)
  if (had_qc_baseline) {
    original_qc_baseline <- get("initializeMetabDesignImportedQcBaseline", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportedQcBaseline",
    function(s4Object, ...) {
      captured$qc_baseline <- list(
        s4Object = s4Object
      )
      stop("qc-baseline-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_qc_baseline) {
      assign("initializeMetabDesignImportedQcBaseline", original_qc_baseline, envir = server_env)
    } else if (exists("initializeMetabDesignImportedQcBaseline", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportedQcBaseline", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$workflow_state$workflowData, workflow_state)
  expect_identical(captured$create_s4$workflowData, workflow_state)
  expect_identical(captured$save_state$workflowData, workflow_state)
  expect_identical(captured$save_state$s4Object, fake_s4)
  expect_identical(captured$qc_baseline$s4Object, fake_s4)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: qc-baseline-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: qc-baseline-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates imported completion seam after QC baseline setup", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    friendly_names = "Treatment_vs_Control",
    full_format = "Treatment_vs_Control=groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  assay_names <- "LCMS_Pos"
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )
  fake_s4 <- structure(list(ok = TRUE), class = "MetaboliteAssayData")
  fake_qc_trigger <- function(value) {
    captured$qc_trigger_value <- value
    invisible(NULL)
  }

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = assay_names,
        assayList = assay_list
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)
  had_metadata <- exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  if (had_metadata) {
    original_metadata <- get("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportMetadata",
    function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      list(
        columnMapping = col_map,
        importedContrasts = imported_contrasts
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_metadata) {
      assign("hydrateMetabDesignImportMetadata", original_metadata, envir = server_env)
    } else if (exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportMetadata", envir = server_env)
    }
  }, add = TRUE)
  had_workflow_state <- exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  if (had_workflow_state) {
    original_workflow_state <- get("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportWorkflowState",
    function(workflowData, importedDesign, assayList, importedContrasts, ...) {
      captured$workflow_state <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        importedContrasts = importedContrasts
      )
      invisible(NULL)
    },
    envir = server_env
  )
  on.exit({
    if (had_workflow_state) {
      assign("hydrateMetabDesignImportWorkflowState", original_workflow_state, envir = server_env)
    } else if (exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportWorkflowState", envir = server_env)
    }
  }, add = TRUE)
  had_create_s4 <- exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  if (had_create_s4) {
    original_create_s4 <- get("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  }
  assign(
    "createMetabDesignImportedS4Object",
    function(workflowData, assayList, colMap, ...) {
      captured$create_s4 <- list(
        workflowData = workflowData,
        assayList = assayList,
        colMap = colMap
      )
      fake_s4
    },
    envir = server_env
  )
  on.exit({
    if (had_create_s4) {
      assign("createMetabDesignImportedS4Object", original_create_s4, envir = server_env)
    } else if (exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)) {
      rm("createMetabDesignImportedS4Object", envir = server_env)
    }
  }, add = TRUE)
  had_save_state <- exists("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)
  if (had_save_state) {
    original_save_state <- get("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)
  }
  assign(
    "saveMetabDesignImportedS4State",
    function(workflowData, s4Object, ...) {
      captured$save_state <- list(
        workflowData = workflowData,
        s4Object = s4Object
      )
      invisible("saved")
    },
    envir = server_env
  )
  on.exit({
    if (had_save_state) {
      assign("saveMetabDesignImportedS4State", original_save_state, envir = server_env)
    } else if (exists("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)) {
      rm("saveMetabDesignImportedS4State", envir = server_env)
    }
  }, add = TRUE)
  had_qc_baseline <- exists("initializeMetabDesignImportedQcBaseline", envir = server_env, inherits = FALSE)
  if (had_qc_baseline) {
    original_qc_baseline <- get("initializeMetabDesignImportedQcBaseline", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportedQcBaseline",
    function(s4Object, ...) {
      captured$qc_baseline <- list(
        s4Object = s4Object
      )
      invisible(TRUE)
    },
    envir = server_env
  )
  on.exit({
    if (had_qc_baseline) {
      assign("initializeMetabDesignImportedQcBaseline", original_qc_baseline, envir = server_env)
    } else if (exists("initializeMetabDesignImportedQcBaseline", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportedQcBaseline", envir = server_env)
    }
  }, add = TRUE)
  had_completion <- exists("completeMetabDesignImportedPostCheckpoint", envir = server_env, inherits = FALSE)
  if (had_completion) {
    original_completion <- get("completeMetabDesignImportedPostCheckpoint", envir = server_env, inherits = FALSE)
  }
  assign(
    "completeMetabDesignImportedPostCheckpoint",
    function(workflowData, assayNames, importedDesign, qcTrigger, ...) {
      captured$completion <- list(
        workflowData = workflowData,
        assayNames = assayNames,
        importedDesign = importedDesign,
        qcTrigger = qcTrigger
      )
      stop("import-completion-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_completion) {
      assign("completeMetabDesignImportedPostCheckpoint", original_completion, envir = server_env)
    } else if (exists("completeMetabDesignImportedPostCheckpoint", envir = server_env, inherits = FALSE)) {
      rm("completeMetabDesignImportedPostCheckpoint", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home"),
    qc_trigger = fake_qc_trigger
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$workflow_state$workflowData, workflow_state)
  expect_identical(captured$create_s4$workflowData, workflow_state)
  expect_identical(captured$save_state$workflowData, workflow_state)
  expect_identical(captured$save_state$s4Object, fake_s4)
  expect_identical(captured$qc_baseline$s4Object, fake_s4)
  expect_identical(captured$completion$workflowData, workflow_state)
  expect_identical(captured$completion$assayNames, assay_names)
  expect_identical(captured$completion$importedDesign, imported_design)
  expect_identical(captured$completion$qcTrigger, fake_qc_trigger)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: import-completion-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: import-completion-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates import S4 creation seam after workflow-state hydration", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    friendly_names = "Treatment_vs_Control",
    full_format = "Treatment_vs_Control=groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = "LCMS_Pos",
        assayList = assay_list
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)
  had_metadata <- exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  if (had_metadata) {
    original_metadata <- get("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportMetadata",
    function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      list(
        columnMapping = col_map,
        importedContrasts = imported_contrasts
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_metadata) {
      assign("hydrateMetabDesignImportMetadata", original_metadata, envir = server_env)
    } else if (exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportMetadata", envir = server_env)
    }
  }, add = TRUE)
  had_workflow_state <- exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  if (had_workflow_state) {
    original_workflow_state <- get("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportWorkflowState",
    function(workflowData, importedDesign, assayList, importedContrasts, ...) {
      captured$workflow_state <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        importedContrasts = importedContrasts
      )
      invisible(NULL)
    },
    envir = server_env
  )
  on.exit({
    if (had_workflow_state) {
      assign("hydrateMetabDesignImportWorkflowState", original_workflow_state, envir = server_env)
    } else if (exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportWorkflowState", envir = server_env)
    }
  }, add = TRUE)
  had_create_s4 <- exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  if (had_create_s4) {
    original_create_s4 <- get("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  }
  assign(
    "createMetabDesignImportedS4Object",
    function(workflowData, assayList, colMap, ...) {
      captured$create_s4 <- list(
        workflowData = workflowData,
        assayList = assayList,
        colMap = colMap
      )
      stop("s4-create-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_create_s4) {
      assign("createMetabDesignImportedS4Object", original_create_s4, envir = server_env)
    } else if (exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)) {
      rm("createMetabDesignImportedS4Object", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$workflow_state$workflowData, workflow_state)
  expect_identical(captured$workflow_state$importedDesign, imported_design)
  expect_identical(captured$workflow_state$assayList, assay_list)
  expect_identical(captured$workflow_state$importedContrasts, imported_contrasts)
  expect_identical(captured$create_s4$workflowData, workflow_state)
  expect_identical(captured$create_s4$assayList, assay_list)
  expect_identical(captured$create_s4$colMap, col_map)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: s4-create-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: s4-create-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates imported state-manager save seam after S4 creation", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    friendly_names = "Treatment_vs_Control",
    full_format = "Treatment_vs_Control=groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )
  fake_s4 <- structure(list(ok = TRUE), class = "MetaboliteAssayData")

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = "LCMS_Pos",
        assayList = assay_list
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)
  had_metadata <- exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  if (had_metadata) {
    original_metadata <- get("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportMetadata",
    function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      list(
        columnMapping = col_map,
        importedContrasts = imported_contrasts
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_metadata) {
      assign("hydrateMetabDesignImportMetadata", original_metadata, envir = server_env)
    } else if (exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportMetadata", envir = server_env)
    }
  }, add = TRUE)
  had_workflow_state <- exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  if (had_workflow_state) {
    original_workflow_state <- get("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportWorkflowState",
    function(workflowData, importedDesign, assayList, importedContrasts, ...) {
      captured$workflow_state <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        importedContrasts = importedContrasts
      )
      invisible(NULL)
    },
    envir = server_env
  )
  on.exit({
    if (had_workflow_state) {
      assign("hydrateMetabDesignImportWorkflowState", original_workflow_state, envir = server_env)
    } else if (exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportWorkflowState", envir = server_env)
    }
  }, add = TRUE)
  had_create_s4 <- exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  if (had_create_s4) {
    original_create_s4 <- get("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)
  }
  assign(
    "createMetabDesignImportedS4Object",
    function(workflowData, assayList, colMap, ...) {
      captured$create_s4 <- list(
        workflowData = workflowData,
        assayList = assayList,
        colMap = colMap
      )
      fake_s4
    },
    envir = server_env
  )
  on.exit({
    if (had_create_s4) {
      assign("createMetabDesignImportedS4Object", original_create_s4, envir = server_env)
    } else if (exists("createMetabDesignImportedS4Object", envir = server_env, inherits = FALSE)) {
      rm("createMetabDesignImportedS4Object", envir = server_env)
    }
  }, add = TRUE)
  had_save_state <- exists("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)
  if (had_save_state) {
    original_save_state <- get("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)
  }
  assign(
    "saveMetabDesignImportedS4State",
    function(workflowData, s4Object, ...) {
      captured$save_state <- list(
        workflowData = workflowData,
        s4Object = s4Object
      )
      stop("state-save-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_save_state) {
      assign("saveMetabDesignImportedS4State", original_save_state, envir = server_env)
    } else if (exists("saveMetabDesignImportedS4State", envir = server_env, inherits = FALSE)) {
      rm("saveMetabDesignImportedS4State", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$workflow_state$workflowData, workflow_state)
  expect_identical(captured$create_s4$workflowData, workflow_state)
  expect_identical(captured$create_s4$assayList, assay_list)
  expect_identical(captured$create_s4$colMap, col_map)
  expect_identical(captured$save_state$workflowData, workflow_state)
  expect_identical(captured$save_state$s4Object, fake_s4)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: state-save-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: state-save-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates import workflow-state seam after metadata hydration", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  imported_contrasts <- data.frame(
    contrasts = "groupTreatment-groupControl",
    friendly_names = "Treatment_vs_Control",
    full_format = "Treatment_vs_Control=groupTreatment-groupControl",
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )
  col_map <- list(
    metabolite_id_col = "Feature",
    annotation_col = "Name",
    sample_columns = c("Sample1", "Sample2"),
    is_pattern = NA_character_
  )

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = "LCMS_Pos",
        assayList = assay_list
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)
  had_metadata <- exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  if (had_metadata) {
    original_metadata <- get("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportMetadata",
    function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      list(
        columnMapping = col_map,
        importedContrasts = imported_contrasts
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_metadata) {
      assign("hydrateMetabDesignImportMetadata", original_metadata, envir = server_env)
    } else if (exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportMetadata", envir = server_env)
    }
  }, add = TRUE)
  had_workflow_state <- exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  if (had_workflow_state) {
    original_workflow_state <- get("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportWorkflowState",
    function(workflowData, importedDesign, assayList, importedContrasts, ...) {
      captured$workflow_state <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        importedContrasts = importedContrasts
      )
      stop("workflow-state-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_workflow_state) {
      assign("hydrateMetabDesignImportWorkflowState", original_workflow_state, envir = server_env)
    } else if (exists("hydrateMetabDesignImportWorkflowState", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportWorkflowState", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$metadata$importedDesign, imported_design)
  expect_identical(captured$metadata$assayList, assay_list)
  expect_identical(captured$workflow_state$workflowData, workflow_state)
  expect_identical(captured$workflow_state$importedDesign, imported_design)
  expect_identical(captured$workflow_state$assayList, assay_list)
  expect_identical(captured$workflow_state$importedContrasts, imported_contrasts)
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: workflow-state-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: workflow-state-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates import metadata hydration seam after artifact loading", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")
  imported_design <- data.frame(
    Run = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      Feature = c("M1", "M2"),
      Sample1 = c(10, 20),
      Sample2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      list(
        importedDesign = imported_design,
        assayNames = "LCMS_Pos",
        assayList = assay_list
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)
  had_metadata <- exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  if (had_metadata) {
    original_metadata <- get("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportMetadata",
    function(workflowData, importedDesign, assayList, colMapFile, contrastFile, ...) {
      captured$metadata <- list(
        workflowData = workflowData,
        importedDesign = importedDesign,
        assayList = assayList,
        colMapFile = colMapFile,
        contrastFile = contrastFile
      )
      stop("metadata-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_metadata) {
      assign("hydrateMetabDesignImportMetadata", original_metadata, envir = server_env)
    } else if (exists("hydrateMetabDesignImportMetadata", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportMetadata", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$hydration$experimentPaths, experiment_paths)
  expect_identical(captured$hydration$importPath, "/tmp/imported-design")
  expect_identical(captured$metadata$workflowData, workflow_state)
  expect_identical(captured$metadata$importedDesign, imported_design)
  expect_identical(captured$metadata$assayList, assay_list)
  expect_identical(captured$metadata$colMapFile, "/tmp/imported-design/column_mapping.json")
  expect_identical(captured$metadata$contrastFile, "/tmp/imported-design/contrast_strings.tab")
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: metadata-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: metadata-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})

test_that("metabolomics design server delegates import artifact hydration seam after preflight", {
  captured <- new.env(parent = emptyenv())
  captured$events <- character()
  server_env <- environment(mod_metab_design_server)

  workflow_state <- list(
    data_tbl = list(LCMS_Pos = data.frame(feature = "M1", stringsAsFactors = FALSE)),
    config_list = list(normalization = "none"),
    column_mapping = list(sample_columns = "Sample1"),
    design_matrix = data.frame(
      Run = "Sample1",
      group = "Control",
      stringsAsFactors = FALSE
    ),
    contrasts_tbl = data.frame(
      contrasts = "groupTreatment-groupControl",
      stringsAsFactors = FALSE
    ),
    tab_status = list(design_matrix = "pending")
  )
  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir()
  )
  input <- list(
    show_import_modal = 1,
    import_dir = "dir-token",
    confirm_import = 0
  )
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste0("metab-", id))
  bootstrap_roots <- c("Project Base Dir" = experiment_paths$base_dir, Home = "/project/home")

  had_builder <- exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  if (had_builder) {
    original_builder <- get("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)
  }
  assign(
    "mod_metab_design_builder_server",
    function(...) function() NULL,
    envir = server_env
  )
  on.exit({
    if (had_builder) {
      assign("mod_metab_design_builder_server", original_builder, envir = server_env)
    } else if (exists("mod_metab_design_builder_server", envir = server_env, inherits = FALSE)) {
      rm("mod_metab_design_builder_server", envir = server_env)
    }
  }, add = TRUE)
  had_bootstrap <- exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  if (had_bootstrap) {
    original_bootstrap <- get("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)
  }
  assign(
    "initializeMetabDesignImportBootstrap",
    function(input, session, experimentPaths, volumes) {
      captured$bootstrap <- list(
        input = input,
        session = session,
        experimentPaths = experimentPaths,
        volumes = volumes
      )
      bootstrap_roots
    },
    envir = server_env
  )
  on.exit({
    if (had_bootstrap) {
      assign("initializeMetabDesignImportBootstrap", original_bootstrap, envir = server_env)
    } else if (exists("initializeMetabDesignImportBootstrap", envir = server_env, inherits = FALSE)) {
      rm("initializeMetabDesignImportBootstrap", envir = server_env)
    }
  }, add = TRUE)
  had_modal_shell <- exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  if (had_modal_shell) {
    original_modal_shell <- get("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)
  }
  assign(
    "registerMetabDesignImportModalShell",
    function(input, output, session, resolvedVolumes) {
      captured$modal_shell <- list(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolvedVolumes
      )
      invisible(output)
    },
    envir = server_env
  )
  on.exit({
    if (had_modal_shell) {
      assign("registerMetabDesignImportModalShell", original_modal_shell, envir = server_env)
    } else if (exists("registerMetabDesignImportModalShell", envir = server_env, inherits = FALSE)) {
      rm("registerMetabDesignImportModalShell", envir = server_env)
    }
  }, add = TRUE)
  had_preflight <- exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  if (had_preflight) {
    original_preflight <- get("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)
  }
  assign(
    "resolveMetabDesignImportPreflight",
    function(input, resolvedVolumes, ...) {
      captured$preflight <- list(
        input = input,
        resolvedVolumes = resolvedVolumes
      )
      list(
        ok = TRUE,
        importPath = "/tmp/imported-design",
        designFile = "/tmp/imported-design/design_matrix.tab",
        contrastFile = "/tmp/imported-design/contrast_strings.tab",
        manifestFile = "/tmp/imported-design/assay_manifest.txt",
        colMapFile = "/tmp/imported-design/column_mapping.json",
        configFile = "/tmp/imported-design/config.ini"
      )
    },
    envir = server_env
  )
  on.exit({
    if (had_preflight) {
      assign("resolveMetabDesignImportPreflight", original_preflight, envir = server_env)
    } else if (exists("resolveMetabDesignImportPreflight", envir = server_env, inherits = FALSE)) {
      rm("resolveMetabDesignImportPreflight", envir = server_env)
    }
  }, add = TRUE)
  had_hydration <- exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  if (had_hydration) {
    original_hydration <- get("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)
  }
  assign(
    "hydrateMetabDesignImportArtifacts",
    function(workflowData, experimentPaths, importPath, designFile, manifestFile, configFile, ...) {
      captured$hydration <- list(
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        importPath = importPath,
        designFile = designFile,
        manifestFile = manifestFile,
        configFile = configFile
      )
      stop("hydration-stop")
    },
    envir = server_env
  )
  on.exit({
    if (had_hydration) {
      assign("hydrateMetabDesignImportArtifacts", original_hydration, envir = server_env)
    } else if (exists("hydrateMetabDesignImportArtifacts", envir = server_env, inherits = FALSE)) {
      rm("hydrateMetabDesignImportArtifacts", envir = server_env)
    }
  }, add = TRUE)

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
      captured$events <- c(captured$events, event_label)

      if (event_label %in% c("input$show_import_modal", "input$confirm_import")) {
        eval(substitute(handlerExpr), parent.frame())
      }

      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(...) invisible(NULL),
    reactiveVal = function(value = NULL) function() value,
    removeModal = function() {
      captured$remove_modal <- TRUE
      invisible(NULL)
    },
    showNotification = function(ui, type = NULL, duration = NULL, id = NULL) {
      captured$notification <- list(ui = ui, type = type, duration = duration, id = id)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$remove_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr, options = NULL) {
      list(value = substitute(expr), options = options)
    },
    .package = "DT"
  )
  testthat::local_mocked_bindings(
    log_error = function(message, ...) {
      captured$log_error <- message
      invisible(NULL)
    },
    .package = "logger"
  )

  mod_metab_design_server(
    id = "design",
    workflow_data = workflow_state,
    experiment_paths = experiment_paths,
    volumes = c(Home = "/project/home")
  )

  expect_identical(captured$module_id, "design")
  expect_identical(captured$bootstrap$input, input)
  expect_identical(captured$modal_shell$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$preflight$input, input)
  expect_identical(captured$preflight$resolvedVolumes, bootstrap_roots)
  expect_identical(captured$hydration$workflowData, workflow_state)
  expect_identical(captured$hydration$experimentPaths, experiment_paths)
  expect_identical(captured$hydration$importPath, "/tmp/imported-design")
  expect_identical(captured$hydration$designFile, "/tmp/imported-design/design_matrix.tab")
  expect_identical(captured$hydration$manifestFile, "/tmp/imported-design/assay_manifest.txt")
  expect_identical(captured$hydration$configFile, "/tmp/imported-design/config.ini")
  expect_true(isTRUE(captured$remove_modal))
  expect_identical(captured$remove_notification, "importing_design")
  expect_identical(captured$log_error, "Error during import: hydration-stop")
  expect_identical(
    captured$notification,
    list(
      ui = "Error during import: hydration-stop",
      type = "error",
      duration = 15,
      id = NULL
    )
  )
  expect_identical(
    captured$events,
    c("input$confirm_import", "builderResultsRv()")
  )
})
