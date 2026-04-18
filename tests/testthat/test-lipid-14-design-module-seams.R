library(testthat)

source(test_path("..", "..", "R", "mod_lipid_design_preview_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_builder_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_import_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_ui.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design_server.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_design.R"), local = environment())

test_that("formatLipidDesignAssaysPreview keeps assay list text stable", {
    expect_identical(
        formatLipidDesignAssaysPreview(list(LCMS_Pos = data.frame(), LCMS_Neg = data.frame())),
        "Included assays: LCMS_Pos, LCMS_Neg"
    )

    expect_identical(
        formatLipidDesignAssaysPreview(list()),
        "Included assays: "
    )
})

test_that("registerLipidDesignPreviewOutputs keeps preview output handoff stable", {
    output <- new.env(parent = emptyenv())
    req_calls <- character()

    render_dt <- function(expr, options = NULL) {
        structure(
            list(
                kind = "dt",
                value = eval.parent(substitute(expr)),
                options = options
            ),
            class = "render_stub"
        )
    }

    render_text <- function(expr) {
        structure(
            list(
                kind = "text",
                value = eval.parent(substitute(expr))
            ),
            class = "render_stub"
        )
    }

    req <- function(value) {
        req_calls <<- c(req_calls, deparse(substitute(value), nlines = 1L))
        value
    }

    workflow_data <- list(
        design_matrix = data.frame(Run = "Sample_1", group = "Control", stringsAsFactors = FALSE),
        contrasts_tbl = data.frame(contrasts = "Treatment-Control", stringsAsFactors = FALSE),
        data_tbl = list(LCMS_Pos = data.frame(), LCMS_Neg = data.frame())
    )

    result <- registerLipidDesignPreviewOutputs(
        output = output,
        workflowData = workflow_data,
        renderDt = render_dt,
        renderText = render_text,
        req = req
    )

    expect_identical(result, output)
    expect_identical(output$design_matrix_preview$value, workflow_data$design_matrix)
    expect_identical(output$contrasts_preview$value, workflow_data$contrasts_tbl)
    expect_identical(output$assays_preview$value, "Included assays: LCMS_Pos, LCMS_Neg")
    expect_identical(output$design_matrix_preview$options, list(pageLength = 5, scrollX = TRUE))
    expect_identical(output$contrasts_preview$options, list(pageLength = 5, scrollX = TRUE))
    expect_identical(
        req_calls,
        c("workflowData$design_matrix", "workflowData$contrasts_tbl", "workflowData$data_tbl")
    )
})

test_that("registerLipidDesignBuilderModule keeps builder registration handoff stable", {
    reactive_calls <- character()

    reactive_fn <- function(expr) {
        reactive_calls <<- c(reactive_calls, deparse(substitute(expr), nlines = 1L))
        eval.parent(substitute(expr))
    }

    workflow_data <- list(
        data_tbl = list(LCMS_Pos = data.frame()),
        config_list = list(mode = "positive"),
        column_mapping = list(sample_id = "Run"),
        design_matrix = data.frame(Run = "Sample_1", group = "Control", stringsAsFactors = FALSE),
        contrasts_tbl = data.frame(contrasts = "Treatment-Control", stringsAsFactors = FALSE)
    )

    builder_result <- list(saved = TRUE)
    builder_calls <- list()

    builder_server_fn <- function(
        moduleId,
        data_tbl,
        config_list,
        column_mapping,
        existing_design_matrix,
        existing_contrasts
    ) {
        builder_calls <<- list(
            moduleId = moduleId,
            data_tbl = data_tbl,
            config_list = config_list,
            column_mapping = column_mapping,
            existing_design_matrix = existing_design_matrix,
            existing_contrasts = existing_contrasts
        )
        builder_result
    }

    result <- registerLipidDesignBuilderModule(
        workflowData = workflow_data,
        builderServerExists = TRUE,
        builderServerFn = builder_server_fn,
        reactiveFn = reactive_fn,
        reactiveValFn = function(value) stop("reactiveVal should not be used when builder exists")
    )

    expect_identical(result, builder_result)
    expect_identical(builder_calls$moduleId, "builder")
    expect_identical(builder_calls$data_tbl, workflow_data$data_tbl)
    expect_identical(builder_calls$config_list, workflow_data$config_list)
    expect_identical(builder_calls$column_mapping, workflow_data$column_mapping)
    expect_identical(builder_calls$existing_design_matrix, workflow_data$design_matrix)
    expect_identical(builder_calls$existing_contrasts, workflow_data$contrasts_tbl)
    expect_identical(
        reactive_calls,
        c(
            "workflowData$data_tbl",
            "workflowData$config_list",
            "workflowData$column_mapping",
            "workflowData$design_matrix",
            "workflowData$contrasts_tbl"
        )
    )
})

test_that("registerLipidDesignBuilderModule falls back to a NULL reactive when builder is unavailable", {
    reactive_val_calls <- list()

    reactive_val_fn <- function(value) {
        reactive_val_calls <<- c(reactive_val_calls, list(value))
        function() value
    }

    result <- registerLipidDesignBuilderModule(
        workflowData = list(),
        builderServerExists = FALSE,
        builderServerFn = function(...) stop("builder server should not be called"),
        reactiveFn = function(expr) stop("reactive should not be used when builder is unavailable"),
        reactiveValFn = reactive_val_fn
    )

    expect_length(reactive_val_calls, 1L)
    expect_null(reactive_val_calls[[1]])
    expect_identical(result(), NULL)
})

test_that("registerLipidDesignBuilderResultsObserver keeps observer handoff stable", {
    builder_result <- list(
        design_matrix = data.frame(Run = "Sample_1", group = "Control", stringsAsFactors = FALSE),
        data_cln = list(LCMS_Pos = data.frame(Run = "Sample_1", stringsAsFactors = FALSE)),
        contrasts_tbl = data.frame(contrasts = "Treatment-Control", stringsAsFactors = FALSE),
        config_list = list(mode = "positive")
    )
    observed_event <- NULL
    observed_ignore_null <- NULL
    req_calls <- list()
    shell_calls <- list()

    observe_event_fn <- function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
        observed_event <<- eval.parent(substitute(eventExpr))
        observed_ignore_null <<- ignoreNULL
        eval.parent(substitute(handlerExpr))
        invisible(NULL)
    }

    req_fn <- function(value) {
        req_calls <<- c(req_calls, list(value))
        value
    }

    run_builder_observer_shell <- function(results, workflowData, experimentPaths, qcTrigger = NULL) {
        shell_calls <<- list(
            results = results,
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            qcTrigger = qcTrigger
        )
    }

    builder_results_rv <- function() builder_result
    workflow_data <- list(design_matrix = NULL)
    experiment_paths <- list(source_dir = tempdir())
    qc_trigger <- function(value) value

    registerLipidDesignBuilderResultsObserver(
        builderResultsRv = builder_results_rv,
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        qcTrigger = qc_trigger,
        observeEventFn = observe_event_fn,
        reqFn = req_fn,
        runBuilderObserverShell = run_builder_observer_shell
    )

    expect_identical(observed_event, builder_result)
    expect_true(observed_ignore_null)
    expect_length(req_calls, 1L)
    expect_identical(req_calls[[1]], builder_result)
    expect_identical(shell_calls$results, builder_result)
    expect_identical(shell_calls$workflowData, workflow_data)
    expect_identical(shell_calls$experimentPaths, experiment_paths)
    expect_identical(shell_calls$qcTrigger, qc_trigger)
})

test_that("initializeLipidDesignImportBootstrap keeps import bootstrap handoff stable", {
    dir_choose_calls <- list()
    log_messages <- character()
    isolate_calls <- 0L

    input <- list(import_dir = NULL)
    session <- list(id = "lipid-design")
    experiment_paths <- list(base_dir = "/tmp/lipid-design-project")

    result <- initializeLipidDesignImportBootstrap(
        input = input,
        session = session,
        experimentPaths = experiment_paths,
        volumes = function() c(Home = "/home/tester"),
        dirChooseFn = function(input, id, roots, session) {
            dir_choose_calls <<- list(
                input = input,
                id = id,
                roots = roots,
                session = session
            )
        },
        isolateFn = function(expr) {
            isolate_calls <<- isolate_calls + 1L
            eval.parent(substitute(expr))
        },
        getVolumesFn = function() stop("getVolumes should not be used when volumes is supplied"),
        dirExistsFn = function(path) identical(path, experiment_paths$base_dir),
        logInfo = function(message) {
            log_messages <<- c(log_messages, message)
        }
    )

    expect_identical(
        result$resolvedVolumes,
        c("Project Base Dir" = experiment_paths$base_dir, Home = "/home/tester")
    )
    expect_identical(dir_choose_calls$input, input)
    expect_identical(dir_choose_calls$id, "import_dir")
    expect_identical(dir_choose_calls$roots, result$resolvedVolumes)
    expect_identical(dir_choose_calls$session, session)
    expect_identical(isolate_calls, 1L)
    expect_identical(log_messages, paste("Added base_dir to volumes:", experiment_paths$base_dir))
})

test_that("buildLipidDesignImportModal keeps import modal contract stable", {
    modal_call <- NULL

    result <- buildLipidDesignImportModal(
        ns = function(id) paste0("lipid-", id),
        modalDialogFn = function(...) {
            modal_call <<- list(...)
            structure(list(kind = "modal", args = list(...)), class = "modal_stub")
        },
        paragraphFn = function(text) list(kind = "p", text = text),
        helpTextFn = function(text) list(kind = "helpText", text = text),
        dirButtonFn = function(id, label, title) list(kind = "dirButton", id = id, label = label, title = title),
        verbatimTextOutputFn = function(id, placeholder = FALSE) {
            list(kind = "verbatimTextOutput", id = id, placeholder = placeholder)
        },
        tagListFn = function(...) list(...),
        modalButtonFn = function(label) list(kind = "modalButton", label = label),
        actionButtonFn = function(id, label, class = NULL) {
            list(kind = "actionButton", id = id, label = label, class = class)
        }
    )

    expect_s3_class(result, "modal_stub")
    expect_identical(modal_call$title, "Import Existing Design Matrix")
    expect_identical(
        modal_call[[2]],
        list(
            kind = "p",
            text = "Select the folder containing 'design_matrix.tab', assay data files, and optionally 'contrast_strings.tab'."
        )
    )
    expect_identical(
        modal_call[[3]],
        list(
            kind = "helpText",
            text = "Required files: design_matrix.tab, assay_manifest.txt, data_cln_*.tab files"
        )
    )
    expect_identical(
        modal_call[[4]],
        list(kind = "dirButton", id = "lipid-import_dir", label = "Select Folder", title = "Choose a directory")
    )
    expect_identical(
        modal_call[[5]],
        list(kind = "verbatimTextOutput", id = "lipid-import_dir_path", placeholder = TRUE)
    )
    expect_identical(
        modal_call$footer,
        list(
            list(kind = "modalButton", label = "Cancel"),
            list(kind = "actionButton", id = "lipid-confirm_import", label = "Import", class = "btn-primary")
        )
    )
})

test_that("registerLipidDesignImportModalShell keeps import modal handoff stable", {
    output <- new.env(parent = emptyenv())
    observed_event <- NULL
    shown_modal <- NULL
    req_calls <- character()
    parse_calls <- list()
    modal_builds <- 0L

    input <- list(
        show_import_modal = TRUE,
        import_dir = "token"
    )
    session <- list(
        ns = function(id) paste0("lipid-", id)
    )
    resolved_volumes <- c("Project Base Dir" = "/tmp/lipid-design-project")

    observe_event_fn <- function(eventExpr, handlerExpr, ignoreNULL = FALSE) {
        observed_event <<- eval.parent(substitute(eventExpr))
        eval.parent(substitute(handlerExpr))
        invisible(NULL)
    }

    render_text_fn <- function(expr) {
        structure(
            list(value = eval.parent(substitute(expr))),
            class = "render_stub"
        )
    }

    req_fn <- function(value) {
        req_calls <<- c(req_calls, deparse(substitute(value), nlines = 1L))
        value
    }

    parse_dir_path_fn <- function(roots, selection) {
        parse_calls <<- list(
            roots = roots,
            selection = selection
        )
        "/tmp/lipid-design-project/imported-design"
    }

    build_import_modal <- function(ns) {
        modal_builds <<- modal_builds + 1L
        list(
            kind = "modal",
            importDirId = ns("import_dir"),
            importPathId = ns("import_dir_path"),
            confirmId = ns("confirm_import")
        )
    }

    result <- registerLipidDesignImportModalShell(
        input = input,
        output = output,
        session = session,
        resolvedVolumes = resolved_volumes,
        observeEventFn = observe_event_fn,
        showModalFn = function(modal) {
            shown_modal <<- modal
        },
        renderTextFn = render_text_fn,
        reqFn = req_fn,
        parseDirPathFn = parse_dir_path_fn,
        buildImportModal = build_import_modal
    )

    expect_identical(result, output)
    expect_true(observed_event)
    expect_identical(modal_builds, 1L)
    expect_identical(
        shown_modal,
        list(
            kind = "modal",
            importDirId = "lipid-import_dir",
            importPathId = "lipid-import_dir_path",
            confirmId = "lipid-confirm_import"
        )
    )
    expect_identical(req_calls, "input$import_dir")
    expect_identical(parse_calls$roots, resolved_volumes)
    expect_identical(parse_calls$selection, input$import_dir)
    expect_identical(output$import_dir_path$value, "/tmp/lipid-design-project/imported-design")
})

test_that("registerLipidDesignImportConfirmationObserver keeps import confirmation handoff stable", {
    observed_event <- NULL
    req_calls <- character()
    parse_calls <- list()
    shell_calls <- list()

    input <- list(
        confirm_import = 1L,
        import_dir = "token"
    )
    resolved_volumes <- c("Project Base Dir" = "/tmp/lipid-design-project")
    workflow_data <- list(tab_status = list(design_matrix = "pending"))
    experiment_paths <- list(source_dir = "/tmp/lipid-design-project")
    qc_trigger <- function(value) value

    observe_event_fn <- function(eventExpr, handlerExpr, ...) {
        observed_event <<- eval.parent(substitute(eventExpr))
        eval.parent(substitute(handlerExpr))
        invisible(NULL)
    }

    req_fn <- function(value) {
        req_calls <<- c(req_calls, deparse(substitute(value), nlines = 1L))
        value
    }

    parse_dir_path_fn <- function(roots, selection) {
        parse_calls <<- list(
            roots = roots,
            selection = selection
        )
        "/tmp/lipid-design-project/imported-design"
    }

    run_import_confirmation_shell <- function(workflowData, experimentPaths, importPath, qcTrigger = NULL, ...) {
        shell_calls <<- list(
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            importPath = importPath,
            qcTrigger = qcTrigger
        )
        invisible(NULL)
    }

    registerLipidDesignImportConfirmationObserver(
        input = input,
        resolvedVolumes = resolved_volumes,
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        qcTrigger = qc_trigger,
        observeEventFn = observe_event_fn,
        reqFn = req_fn,
        parseDirPathFn = parse_dir_path_fn,
        runImportConfirmationShell = run_import_confirmation_shell
    )

    expect_identical(observed_event, 1L)
    expect_identical(req_calls, c("input$import_dir", "importPath"))
    expect_identical(parse_calls$roots, resolved_volumes)
    expect_identical(parse_calls$selection, input$import_dir)
    expect_identical(shell_calls$workflowData, workflow_data)
    expect_identical(shell_calls$experimentPaths, experiment_paths)
    expect_identical(shell_calls$importPath, "/tmp/lipid-design-project/imported-design")
    expect_identical(shell_calls$qcTrigger, qc_trigger)
})
