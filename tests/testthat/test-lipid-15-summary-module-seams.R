library(testthat)

source(test_path("..", "..", "R", "mod_lipid_summary_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_summary.R"), local = environment())

test_that("buildLipidSummaryTemplateStatus keeps missing-project warnings stable", {
    expect_identical(
        buildLipidSummaryTemplateStatus(
            projectDirs = list(other = list(base_dir = tempdir())),
            omicType = "lipidomics"
        ),
        "[WARNING] Project directories not available"
    )
})

test_that("buildLipidSummaryTemplateStatus keeps template status text stable", {
    base_dir <- file.path(tempdir(), "lipid-summary-template-status")
    template_dir <- file.path(base_dir, "scripts", "lipidomics")
    dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)

    project_dirs <- list(
        lipidomics = list(base_dir = base_dir)
    )

    expect_identical(
        buildLipidSummaryTemplateStatus(
            projectDirs = project_dirs,
            omicType = "lipidomics"
        ),
        "[WARNING] Report template will be downloaded when generating report"
    )

    template_path <- file.path(template_dir, "lipidomics_report.rmd")
    writeLines("---", template_path)

    expect_identical(
        buildLipidSummaryTemplateStatus(
            projectDirs = project_dirs,
            omicType = "lipidomics"
        ),
        "Template: Lipidomics Report [OK]"
    )
})

test_that("registerLipidSummaryTemplateStatusOutput keeps render handoff stable", {
    output <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- list()
    project_dirs <- list(
        lipidomics = list(base_dir = "/tmp/project")
    )

    result <- registerLipidSummaryTemplateStatusOutput(
        output = output,
        projectDirs = project_dirs,
        omicType = "lipidomics",
        renderTextFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render_text"), class = "mock_render_text")
        },
        reqFn = function(value) {
            captured$reqValue <<- value
            value
        },
        buildStatusFn = function(projectDirsArg, omicTypeArg) {
            captured$projectDirs <<- projectDirsArg
            captured$omicType <<- omicTypeArg
            "template-status"
        }
    )

    expect_identical(result, output)
    expect_identical(render_calls, "template-status")
    expect_identical(captured$reqValue, project_dirs)
    expect_identical(captured$projectDirs, project_dirs)
    expect_identical(captured$omicType, "lipidomics")
    expect_s3_class(output$template_status, "mock_render_text")
})

test_that("registerLipidSummaryInitialOutputs keeps wrapper output setup stable", {
    output <- new.env(parent = emptyenv())
    render_calls <- character()
    reactive_calls <- logical()
    output_options <- list()

    result <- registerLipidSummaryInitialOutputs(
        output = output,
        renderTextFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render_text"), class = "mock_render_text")
        },
        reactiveFn = function(expr) {
            reactive_calls <<- c(reactive_calls, force(expr))
            structure(list(kind = "mock_reactive"), class = "mock_reactive")
        },
        outputOptionsFn = function(outputArg, name, suspendWhenHidden) {
            output_options <<- c(output_options, list(list(output = outputArg, name = name, suspendWhenHidden = suspendWhenHidden)))
            invisible(NULL)
        }
    )

    expect_identical(result, output)
    expect_identical(render_calls, "Ready to save workflow parameters and generate report")
    expect_identical(reactive_calls, FALSE)
    expect_s3_class(output$session_summary, "mock_render_text")
    expect_s3_class(output$report_ready, "mock_reactive")
    expect_identical(
        output_options,
        list(list(output = output, name = "report_ready", suspendWhenHidden = FALSE))
    )
})

test_that("initializeLipidSummarySessionBootstrap keeps bootstrap handoff stable", {
    prefills <- list()
    reactive_seed <- NULL

    result <- initializeLipidSummarySessionBootstrap(
        session = "mock-session",
        experimentLabel = "lipid-study",
        updateTextInputFn = function(session, inputId, value) {
            prefills <<- c(prefills, list(list(session = session, inputId = inputId, value = value)))
            invisible(NULL)
        },
        reactiveValuesFn = function(...) {
            reactive_seed <<- list(...)
            structure(list(kind = "mock_reactive_values"), class = "mock_reactive_values")
        }
    )

    expect_s3_class(result, "mock_reactive_values")
    expect_identical(
        prefills,
        list(list(session = "mock-session", inputId = "experiment_label", value = "lipid-study"))
    )
    expect_identical(
        reactive_seed,
        list(
            workflow_args_saved = FALSE,
            files_copied = FALSE,
            report_generated = FALSE,
            report_path = NULL
        )
    )
})

test_that("initializeLipidSummarySessionBootstrap skips label prefill when absent", {
    prefill_calls <- 0L

    result <- initializeLipidSummarySessionBootstrap(
        session = "mock-session",
        experimentLabel = NULL,
        updateTextInputFn = function(session, inputId, value) {
            prefill_calls <<- prefill_calls + 1L
            invisible(NULL)
        },
        reactiveValuesFn = function(...) list(...)
    )

    expect_identical(prefill_calls, 0L)
    expect_identical(
        result,
        list(
            workflow_args_saved = FALSE,
            files_copied = FALSE,
            report_generated = FALSE,
            report_path = NULL
        )
    )
})

test_that("collectLipidSummaryWorkflowArgsContext keeps state resolution stable", {
    if (!methods::isClass("MockLipidSummaryS4")) {
        methods::setClass("MockLipidSummaryS4", slots = c(args = "list"))
    }

    final_s4_object <- methods::new("MockLipidSummaryS4", args = list(normalization = list(method = "median")))
    assigned_env <- new.env(parent = emptyenv())
    log_messages <- character()

    workflow_data <- list(
        state_manager = list(
            getHistory = function() c("lipid_qc_complete", "lipid_normalized"),
            getState = function(state) {
                expect_identical(state, "lipid_normalized")
                final_s4_object
            }
        ),
        contrasts_tbl = data.frame(lhs = "A", rhs = "B", stringsAsFactors = FALSE),
        config_list = list(alpha = 1, beta = 2)
    )

    result <- collectLipidSummaryWorkflowArgsContext(
        workflowData = workflow_data,
        assignFn = function(x, value, envir) {
            assign(x, value, envir = envir)
            invisible(NULL)
        },
        globalEnv = assigned_env,
        catFn = function(...) {
            log_messages <<- c(log_messages, paste0(..., collapse = ""))
            invisible(NULL)
        }
    )

    expect_identical(result$finalS4Object, final_s4_object)
    expect_identical(result$contrastsTbl, workflow_data$contrasts_tbl)
    expect_identical(get("config_list", envir = assigned_env), workflow_data$config_list)
    expect_true(any(grepl("Retrieved DATA S4 object from state 'lipid_normalized'", log_messages, fixed = TRUE)))
    expect_true(any(grepl("Using contrasts_tbl from workflow_data", log_messages, fixed = TRUE)))
})

test_that("handleLipidSummarySaveWorkflowArgs keeps success side effects stable", {
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- FALSE
    output <- new.env(parent = emptyenv())
    notifications <- list()
    created_dirs <- character()
    create_args <- NULL
    saved <- list()

    result <- handleLipidSummarySaveWorkflowArgs(
        input = list(
            experiment_label = "lipid-study",
            description = "workflow summary"
        ),
        values = values,
        output = output,
        projectDirs = list(
            lipidomics = list(
                source_dir = "/tmp/lipid-source",
                base_dir = "/tmp/lipid-base"
            )
        ),
        omicType = "lipidomics",
        workflowData = list(marker = "workflow-data"),
        collectContextFn = function(workflowData, catFn) {
            expect_identical(workflowData$marker, "workflow-data")
            list(
                finalS4Object = list(kind = "final-s4"),
                contrastsTbl = list(kind = "contrasts")
            )
        },
        createWorkflowArgsFn = function(...) {
            create_args <<- list(...)
            "/tmp/lipid-source/study_parameters.txt"
        },
        saveRDSFn = function(object, file) {
            saved <<- list(object = object, file = file)
            invisible(file)
        },
        showNotificationFn = function(message, type, duration = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type, duration = duration)))
            invisible(NULL)
        },
        renderTextFn = function(expr) force(expr),
        dirExistsFn = function(path) FALSE,
        dirCreateFn = function(path, recursive, showWarnings) {
            created_dirs <<- c(created_dirs, path)
            invisible(TRUE)
        },
        timeFn = function() as.POSIXct("2026-04-16 13:00:00", tz = "UTC"),
        catFn = function(...) invisible(NULL)
    )

    expect_true(values$workflow_args_saved)
    expect_identical(result$studyParamsFile, "/tmp/lipid-source/study_parameters.txt")
    expect_identical(result$s4Filename, "lipidomics_lipid-study_final_s4.RDS")
    expect_identical(result$s4Filepath, "/tmp/lipid-base/integration/lipidomics_lipid-study_final_s4.RDS")
    expect_identical(created_dirs, "/tmp/lipid-base/integration")
    expect_identical(create_args$workflow_name, "lipid-study")
    expect_identical(create_args$description, "workflow summary")
    expect_identical(create_args$source_dir_path, "/tmp/lipid-source")
    expect_identical(create_args$final_s4_object, list(kind = "final-s4"))
    expect_identical(create_args$contrasts_tbl, list(kind = "contrasts"))
    expect_identical(saved$object, list(kind = "final-s4"))
    expect_identical(saved$file, result$s4Filepath)
    expect_match(output$session_summary, "Study parameters created for: lipid-study", fixed = TRUE)
    expect_identical(
        notifications,
        list(
            list(message = "Saved Integration S4 Object", type = "message", duration = NULL),
            list(message = "Study parameters saved successfully", type = "message", duration = NULL)
        )
    )
})

test_that("registerLipidSummarySaveWorkflowArgsObserver keeps the observer shell handoff stable", {
    input <- list(
        save_workflow_args = 12L,
        experiment_label = "lipid-study",
        description = "workflow summary"
    )
    observer_event <- NULL
    req_value <- NULL
    handler_args <- NULL
    values <- list(workflow_args_saved = FALSE)
    output <- new.env(parent = emptyenv())
    project_dirs <- list(lipidomics = list(source_dir = tempdir()))
    workflow_data <- list(marker = "workflow-data")

    result <- registerLipidSummarySaveWorkflowArgsObserver(
        input = input,
        values = values,
        output = output,
        projectDirs = project_dirs,
        omicType = "lipidomics",
        workflowData = workflow_data,
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        reqFn = function(value) {
            req_value <<- value
            value
        },
        handleSaveFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 12L)
    expect_identical(req_value, "lipid-study")
    expect_setequal(names(handler_args), c("input", "values", "output", "projectDirs", "omicType", "workflowData"))
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$values, values)
    expect_identical(handler_args$output, output)
    expect_identical(handler_args$projectDirs, project_dirs)
    expect_identical(handler_args$omicType, "lipidomics")
    expect_identical(handler_args$workflowData, workflow_data)
})

test_that("handleLipidSummaryCopyToPublication keeps success side effects stable", {
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- FALSE
    values$files_copied <- FALSE
    output <- new.env(parent = emptyenv())
    progress_messages <- character()
    written_file <- NULL
    read_paths <- character()
    assigned_env <- new.env(parent = emptyenv())
    copy_args <- NULL
    notifications <- list()

    source_dir <- file.path(tempdir(), "lipid-summary-copy-source")
    design_path <- file.path(source_dir, "design_matrix.tab")
    contrasts_path <- file.path(source_dir, "contrasts_tbl.tab")
    basic_params_path <- file.path(source_dir, "study_parameters.txt")

    result <- handleLipidSummaryCopyToPublication(
        input = list(
            experiment_label = "lipid-study",
            description = "copy summary"
        ),
        values = values,
        output = output,
        projectDirs = list(
            lipidomics = list(source_dir = source_dir)
        ),
        omicType = "lipidomics",
        workflowData = list(),
        withProgressFn = function(message, expr) {
            progress_messages <<- c(progress_messages, message)
            force(expr)
        },
        fileExistsFn = function(path) {
            path %in% c(design_path, contrasts_path)
        },
        writeLinesFn = function(text, con) {
            written_file <<- list(text = text, con = con)
            invisible(NULL)
        },
        timeFn = function() as.POSIXct("2026-04-16 16:00:00", tz = "UTC"),
        readTsvFn = function(path, show_col_types = FALSE) {
            read_paths <<- c(read_paths, path)
            if (identical(path, design_path)) {
                return(data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE))
            }

            if (identical(path, contrasts_path)) {
                return(data.frame(lhs = "A", rhs = "B", stringsAsFactors = FALSE))
            }

            stop("unexpected path")
        },
        existsFn = function(name, envir) FALSE,
        assignFn = function(x, value, envir) {
            assign(x, value, envir = envir)
            invisible(NULL)
        },
        globalEnv = assigned_env,
        copyResultsSummaryFn = function(...) {
            copy_args <<- list(...)
            invisible(TRUE)
        },
        showNotificationFn = function(message, type, duration = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type, duration = duration)))
            invisible(NULL)
        },
        renderTextFn = function(expr) force(expr),
        logErrorFn = function(message) stop(message),
        skipFormatterFn = identity,
        catFn = function(...) invisible(NULL),
        tracebackFn = function() invisible(NULL)
    )

    expect_identical(progress_messages, "Copying files to publication directory...")
    expect_true(values$workflow_args_saved)
    expect_true(values$files_copied)
    expect_identical(written_file$con, basic_params_path)
    expect_match(written_file$text, "Workflow Name: lipid-study", fixed = TRUE)
    expect_identical(read_paths, c(design_path, contrasts_path))
    expect_identical(get("project_dirs", envir = assigned_env), list(lipidomics = list(source_dir = source_dir)))
    expect_identical(copy_args$omic_type, "lipidomics")
    expect_identical(copy_args$experiment_label, "lipid-study")
    expect_true(copy_args$force)
    expect_identical(
        copy_args$design_matrix,
        data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE)
    )
    expect_identical(
        copy_args$contrasts_tbl,
        data.frame(lhs = "A", rhs = "B", stringsAsFactors = FALSE)
    )
    expect_identical(
        result,
        list(
            contrastsTbl = data.frame(lhs = "A", rhs = "B", stringsAsFactors = FALSE),
            designMatrix = data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE)
        )
    )
    expect_identical(output$copy_status, "Files copied to publication directory successfully [OK]")
    expect_match(output$session_summary, "Files copied \\[OK\\]", perl = TRUE)
    expect_identical(
        notifications,
        list(list(
            message = "Publication files copied",
            type = "message",
            duration = NULL
        ))
    )
})

test_that("registerLipidSummaryCopyToPublicationObserver keeps the observer shell handoff stable", {
    input <- list(
        copy_to_publication = 29L,
        experiment_label = "lipid-study",
        description = "copy summary"
    )
    values <- list(workflow_args_saved = TRUE, files_copied = FALSE)
    output <- new.env(parent = emptyenv())
    observer_event <- NULL
    observer_ignore_init <- NULL
    req_value <- NULL
    handler_args <- NULL
    project_dirs <- list(lipidomics = list(source_dir = tempdir()))
    workflow_data <- list(marker = "workflow-data")

    result <- registerLipidSummaryCopyToPublicationObserver(
        input = input,
        values = values,
        output = output,
        projectDirs = project_dirs,
        omicType = "lipidomics",
        workflowData = workflow_data,
        observeEventFn = function(event, handler, ignoreInit = FALSE) {
            observer_event <<- event
            observer_ignore_init <<- ignoreInit
            force(handler)
            "registered"
        },
        reqFn = function(value) {
            req_value <<- value
            value
        },
        handleCopyFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 29L)
    expect_true(observer_ignore_init)
    expect_identical(req_value, "lipid-study")
    expect_setequal(names(handler_args), c("input", "values", "output", "projectDirs", "omicType", "workflowData"))
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$values, values)
    expect_identical(handler_args$output, output)
    expect_identical(handler_args$projectDirs, project_dirs)
    expect_identical(handler_args$omicType, "lipidomics")
    expect_identical(handler_args$workflowData, workflow_data)
})

test_that("handleLipidSummaryGenerateReport keeps success side effects stable", {
    values <- new.env(parent = emptyenv())
    values$files_copied <- TRUE
    values$report_generated <- FALSE
    values$report_path <- NULL

    output <- new.env(parent = emptyenv())
    notifications <- list()
    info_messages <- character()
    progress_messages <- character()
    progress_updates <- list()
    output_options <- list()
    render_args <- NULL
    download_handler <- NULL

    base_dir <- file.path(tempdir(), "lipid-summary-report-success")
    template_dir <- file.path(base_dir, "scripts", "lipidomics")
    dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)

    template_path <- file.path(template_dir, "lipidomics_report.rmd")
    writeLines("---", template_path)

    rendered_path <- file.path(base_dir, "reports", "lipid-report.html")
    dir.create(dirname(rendered_path), recursive = TRUE, showWarnings = FALSE)
    writeLines("<html></html>", rendered_path)

    handleLipidSummaryGenerateReport(
        input = list(
            experiment_label = "lipid-study",
            description = "report summary"
        ),
        values = values,
        output = output,
        projectDirs = list(
            lipidomics = list(base_dir = base_dir)
        ),
        omicType = "lipidomics",
        withProgressFn = function(message, expr) {
            progress_messages <<- c(progress_messages, message)
            force(expr)
        },
        showNotificationFn = function(message, type, duration = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type, duration = duration)))
            invisible(NULL)
        },
        incProgressFn = function(value, detail = NULL) {
            progress_updates <<- c(progress_updates, list(list(value = value, detail = detail)))
            invisible(NULL)
        },
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
            invisible(message)
        },
        logErrorFn = function(message) stop(message),
        renderReportExistsFn = function() TRUE,
        renderReportFn = function(...) {
            render_args <<- list(...)
            rendered_path
        },
        reactiveFn = function(expr) force(expr),
        outputOptionsFn = function(outputArg, name, suspendWhenHidden) {
            output_options <<- c(output_options, list(list(name = name, suspendWhenHidden = suspendWhenHidden)))
            invisible(NULL)
        },
        downloadHandlerFn = function(filename, content) {
            download_handler <<- list(filename = filename, content = content)
            "download-handler"
        },
        renderTextFn = function(expr) force(expr),
        timeFn = function() as.POSIXct("2026-04-16 14:30:00", tz = "UTC"),
        catFn = function(...) invisible(NULL),
        printFn = function(...) invisible(NULL)
    )

    expect_true(values$report_generated)
    expect_identical(values$report_path, rendered_path)
    expect_identical(progress_messages, "Generating report...")
    expect_length(progress_updates, 1)
    expect_identical(progress_updates[[1]]$value, 0.5)
    expect_identical(progress_updates[[1]]$detail, "Rendering report...")
    expect_identical(
        render_args,
        list(
            omic_type = "lipidomics",
            experiment_label = "lipid-study",
            rmd_filename = "lipidomics_report.rmd"
        )
    )
    expect_identical(info_messages[1], "Calling RenderReport with omic_type: {omicType}, experiment_label: {input$experiment_label}")
    expect_identical(info_messages[2], "RenderReport returned path: {renderedPath}")
    expect_identical(output$report_ready, TRUE)
    expect_identical(output$download_report, "download-handler")
    expect_identical(download_handler$filename(), "lipid-report.html")
    expect_match(output$session_summary, "Report generated \\[OK\\]", perl = TRUE)
    expect_match(output$session_summary, rendered_path, fixed = TRUE)
    expect_identical(
        output_options,
        list(list(name = "report_ready", suspendWhenHidden = FALSE))
    )
    expect_identical(
        notifications,
        list(list(
            message = "Report generated successfully!",
            type = "message",
            duration = NULL
        ))
    )
})

test_that("registerLipidSummaryGenerateReportObserver keeps the observer shell handoff stable", {
    input <- list(
        generate_report = 33L,
        experiment_label = "lipid-study"
    )
    values <- list(files_copied = TRUE, report_generated = FALSE)
    output <- new.env(parent = emptyenv())
    observer_event <- NULL
    req_values <- list()
    handler_args <- NULL

    result <- registerLipidSummaryGenerateReportObserver(
        input = input,
        values = values,
        output = output,
        projectDirs = list(lipidomics = list(base_dir = tempdir())),
        omicType = "lipidomics",
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        reqFn = function(value) {
            req_values <<- c(req_values, list(value))
            value
        },
        handleGenerateReportFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 33L)
    expect_identical(req_values, list("lipid-study", TRUE))
    expect_setequal(names(handler_args), c("input", "values", "output", "projectDirs", "omicType"))
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$values, values)
    expect_identical(handler_args$output, output)
    expect_identical(handler_args$projectDirs$lipidomics$base_dir, tempdir())
    expect_identical(handler_args$omicType, "lipidomics")
})

test_that("handleLipidSummaryPushToGithub keeps success side effects stable", {
    output <- new.env(parent = emptyenv())
    notifications <- list()
    progress_messages <- character()
    options_args <- NULL
    push_args <- NULL

    handleLipidSummaryPushToGithub(
        input = list(
            enable_github = TRUE,
            github_org = "apaf-bioinformatics",
            github_email = "team@example.org",
            github_username = "lipid-user",
            project_id = "proj-001",
            experiment_label = "lipid-study",
            description = "github push"
        ),
        values = list(report_generated = TRUE),
        output = output,
        projectDirs = list(
            lipidomics = list(base_dir = tempdir())
        ),
        omicType = "lipidomics",
        withProgressFn = function(message, expr) {
            progress_messages <<- c(progress_messages, message)
            force(expr)
        },
        optionsFn = function(...) {
            options_args <<- list(...)
            invisible(NULL)
        },
        pushProjectFn = function(...) {
            push_args <<- list(...)
            invisible(TRUE)
        },
        showNotificationFn = function(message, type, duration = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type, duration = duration)))
            invisible(NULL)
        },
        logErrorFn = function(message) stop(message),
        skipFormatterFn = identity,
        renderTextFn = function(expr) force(expr),
        timeFn = function() as.POSIXct("2026-04-16 15:00:00", tz = "UTC")
    )

    expect_identical(progress_messages, "Pushing to GitHub...")
    expect_identical(
        options_args,
        list(
            github_org = "apaf-bioinformatics",
            github_user_email = "team@example.org",
            github_user_name = "lipid-user"
        )
    )
    expect_identical(
        push_args,
        list(
            project_dirs = list(lipidomics = list(base_dir = tempdir())),
            omic_type = "lipidomics",
            experiment_label = "lipid-study",
            project_id = "proj-001"
        )
    )
    expect_match(output$session_summary, "GitHub pushed \\[OK\\]", perl = TRUE)
    expect_match(output$session_summary, "2026-04-16 15:00:00", fixed = TRUE)
    expect_identical(
        notifications,
        list(list(
            message = "Successfully pushed to GitHub",
            type = "message",
            duration = NULL
        ))
    )
})

test_that("registerLipidSummaryPushToGithubObserver keeps the observer shell handoff stable", {
    input <- list(
        push_to_github = 41L,
        enable_github = TRUE,
        github_org = "apaf-bioinformatics",
        github_email = "team@example.org",
        github_username = "lipid-user",
        project_id = "proj-001",
        experiment_label = "lipid-study"
    )
    values <- list(report_generated = TRUE)
    output <- new.env(parent = emptyenv())
    observer_event <- NULL
    req_values <- list()
    handler_args <- NULL

    result <- registerLipidSummaryPushToGithubObserver(
        input = input,
        values = values,
        output = output,
        projectDirs = list(lipidomics = list(base_dir = tempdir())),
        omicType = "lipidomics",
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        reqFn = function(...) {
            args <- list(...)
            req_values <<- c(req_values, list(args))
            if (length(args) == 1L) {
                args[[1L]]
            } else {
                args
            }
        },
        handlePushFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 41L)
    expect_identical(
        req_values,
        list(
            list(TRUE, "apaf-bioinformatics", "team@example.org", "lipid-user", "proj-001"),
            list(TRUE)
        )
    )
    expect_setequal(names(handler_args), c("input", "values", "output", "projectDirs", "omicType"))
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$values, values)
    expect_identical(handler_args$output, output)
    expect_identical(handler_args$projectDirs$lipidomics$base_dir, tempdir())
    expect_identical(handler_args$omicType, "lipidomics")
})

test_that("buildLipidSummarySessionState keeps export payload stable", {
    fixed_time <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")
    input <- list(
        experiment_label = "lipid-study",
        description = "session export"
    )
    values <- list(
        workflow_args_saved = TRUE,
        files_copied = FALSE,
        report_generated = TRUE,
        report_path = "/tmp/report.html"
    )
    project_dirs <- list(
        lipidomics = list(source_dir = "/tmp/lipid-source")
    )

    result <- buildLipidSummarySessionState(
        input = input,
        values = values,
        projectDirs = project_dirs,
        omicType = "lipidomics",
        timeFn = function() fixed_time
    )

    expect_identical(
        result,
        list(
            experiment_label = "lipid-study",
            description = "session export",
            timestamp = fixed_time,
            omic_type = "lipidomics",
            workflow_args_saved = TRUE,
            files_copied = FALSE,
            report_generated = TRUE,
            report_path = "/tmp/report.html",
            project_dirs = project_dirs
        )
    )
})

test_that("handleLipidSummaryExportSessionState keeps export side effects stable", {
    saved <- list()
    notifications <- list()
    info_message <- NULL
    error_message <- NULL
    project_dirs <- list(
        lipidomics = list(source_dir = tempdir())
    )

    result <- handleLipidSummaryExportSessionState(
        input = list(experiment_label = "lipid-study", description = "session export"),
        values = list(
            workflow_args_saved = TRUE,
            files_copied = TRUE,
            report_generated = FALSE,
            report_path = NULL
        ),
        projectDirs = project_dirs,
        omicType = "lipidomics",
        saveRDSFn = function(object, file) {
            saved <<- list(object = object, file = file)
            invisible(file)
        },
        showNotificationFn = function(message, type, duration = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type, duration = duration)))
            invisible(NULL)
        },
        logInfoFn = function(message) {
            info_message <<- message
            invisible(message)
        },
        logErrorFn = function(message) {
            error_message <<- message
            invisible(message)
        },
        dateFn = function() as.Date("2026-04-16"),
        buildSessionStateFn = function(...) list(kind = "session-state")
    )

    expect_identical(
        result,
        file.path(tempdir(), "session_state_2026-04-16.RDS")
    )
    expect_identical(saved$object, list(kind = "session-state"))
    expect_identical(saved$file, result)
    expect_identical(
        notifications,
        list(list(
            message = paste("Session state exported to:", result),
            type = "message",
            duration = NULL
        ))
    )
    expect_identical(info_message, paste("Session state exported to:", result))
    expect_null(error_message)
})

test_that("registerLipidSummaryExportSessionObserver keeps the observer shell handoff stable", {
    input <- list(
        export_session_state = 27L,
        experiment_label = "lipid-study"
    )
    observer_event <- NULL
    req_value <- NULL
    handler_args <- NULL

    result <- registerLipidSummaryExportSessionObserver(
        input = input,
        values = list(report_generated = FALSE),
        projectDirs = list(lipidomics = list(source_dir = tempdir())),
        omicType = "lipidomics",
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        reqFn = function(value) {
            req_value <<- value
            value
        },
        handleExportFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 27L)
    expect_identical(req_value, "lipid-study")
    expect_setequal(names(handler_args), c("input", "values", "projectDirs", "omicType"))
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$values$report_generated, FALSE)
    expect_identical(handler_args$projectDirs$lipidomics$source_dir, tempdir())
    expect_identical(handler_args$omicType, "lipidomics")
})
