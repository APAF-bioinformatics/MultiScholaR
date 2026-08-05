bootstrapProtSummaryCopyFallbackStudyParams <- function(values,
                                                        projectDirs,
                                                        omicType = "proteomics",
                                                        experimentLabel,
                                                        description,
                                                        fileExistsFn = file.exists,
                                                        writeLinesFn = writeLines,
                                                        timestampFn = Sys.time,
                                                        logErrorFn = function(message) {
                                                          logger::log_error(
                                                            logger::skip_formatter(message)
                                                          )
                                                        },
                                                        catFn = cat) {
  if (isTRUE(values$workflow_args_saved)) {
    return(FALSE)
  }

  basicParamsFile <- file.path(
    projectDirs[[omicType]]$source_dir,
    "study_parameters.txt"
  )

  if (fileExistsFn(basicParamsFile)) {
    return(FALSE)
  }

  catFn("   SESSION SUMMARY: Creating basic study_parameters.txt as fallback\n")

  basicContent <- paste(
    "Study Parameters",
    "================",
    "",
    paste("Workflow Name:", experimentLabel),
    paste("Description:", description),
    paste("Timestamp:", timestampFn()),
    paste("Note: Some parameters could not be saved due to serialization issues"),
    sep = "\n"
  )

  tryCatch({
    writeLinesFn(basicContent, basicParamsFile)
    values$workflow_args_saved <- TRUE
    TRUE
  }, error = function(e) {
    logErrorFn(paste("Failed to create basic study_parameters.txt:", e$message))
    FALSE
  })
}

prepareProtSummaryCopyInputs <- function(workflowData,
                                         projectDirs,
                                         omicType = "proteomics",
                                         readTsvFn = readr::read_tsv,
                                         catFn = cat) {
  contrastsTbl <- NULL
  designMatrix <- NULL

  if (!is.null(workflowData)) {
    if (!is.null(workflowData$contrasts_tbl)) {
      contrastsTbl <- workflowData$contrasts_tbl
      catFn("   SESSION SUMMARY Step: Got contrasts_tbl from workflow_data\n")
    }
    if (!is.null(workflowData$design_matrix)) {
      designMatrix <- workflowData$design_matrix
      catFn("   SESSION SUMMARY Step: Got design_matrix from workflow_data\n")
    }
  }

  if (is.null(designMatrix)) {
    designMatrixFile <- file.path(projectDirs[[omicType]]$source_dir, "design_matrix.tab")
    if (file.exists(designMatrixFile)) {
      designMatrix <- readTsvFn(designMatrixFile, show_col_types = FALSE)
      catFn("   SESSION SUMMARY Step: Loaded design_matrix from file\n")
    }
  }

  if (is.null(contrastsTbl)) {
    contrastsFile <- file.path(projectDirs[[omicType]]$source_dir, "contrasts_tbl.tab")
    if (file.exists(contrastsFile)) {
      contrastsTbl <- readTsvFn(contrastsFile, show_col_types = FALSE)
      catFn("   SESSION SUMMARY Step: Loaded contrasts_tbl from file\n")
    }
  }

  list(contrastsTbl = contrastsTbl, designMatrix = designMatrix)
}

runProtSummaryPublicationCopy <- function(output,
                                          values,
                                          projectDirs,
                                          omicType = "proteomics",
                                          experimentLabel,
                                          description,
                                          contrastsTbl = NULL,
                                          designMatrix = NULL,
                                          existsFn = exists,
                                          assignFn = assign,
                                          copyFn = copyToResultsSummary,
                                          renderTextFn = shiny::renderText,
                                          showNotificationFn = shiny::showNotification,
                                          timestampFn = Sys.time,
                                          catFn = cat) {
  if (!existsFn("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
    catFn("   SESSION SUMMARY: Setting project_dirs in global environment\n")
    assignFn("project_dirs", projectDirs, envir = .GlobalEnv)
  }

  catFn(
    "   SESSION SUMMARY: project_dirs keys:",
    paste(names(projectDirs), collapse = ", "),
    "\n"
  )
  catFn(
    "   SESSION SUMMARY: Using omic_type =",
    omicType,
    "experiment_label =",
    experimentLabel,
    "\n"
  )

  copyArgs <- list(
    omic_type = omicType,
    experiment_label = experimentLabel,
    contrasts_tbl = contrastsTbl,
    design_matrix = designMatrix,
    force = TRUE
  )
  do.call(copyFn, copyArgs)

  values$files_copied <- TRUE

  output$copy_status <- renderTextFn(
    "Files copied to publication directory successfully [OK]"
  )
  showNotificationFn("Publication files copied", type = "message")

  output$session_summary <- renderTextFn({
    paste(
      "Workflow args created for:",
      experimentLabel,
      "\nDescription:",
      description,
      "\nTimestamp:",
      timestampFn(),
      "\nStatus: Arguments saved [OK], Files copied [OK]"
    )
  })

  copyArgs
}

handleProtSummaryPublicationCopyError <- function(output,
                                                  error,
                                                  renderTextFn = shiny::renderText,
                                                  showNotificationFn = shiny::showNotification,
                                                  logErrorFn = function(message) {
                                                    logger::log_error(
                                                      logger::skip_formatter(message)
                                                    )
                                                  },
                                                  catFn = cat,
                                                  tracebackFn = traceback) {
  output$copy_status <- renderTextFn(paste("Error:", error$message))
  showNotificationFn(
    paste("Copy error:", error$message),
    type = "error",
    duration = 10
  )
  logErrorFn(paste("Failed to copy files:", error$message))
  catFn("   SESSION SUMMARY ERROR:", error$message, "\n")
  catFn("   SESSION SUMMARY Traceback:\n")
  tracebackFn()

  FALSE
}

runProtSummaryGithubPush <- function(projectDirs,
                                     omicType = "proteomics",
                                     experimentLabel,
                                     githubOrg,
                                     githubEmail,
                                     githubUsername,
                                     projectId,
                                     optionsFn = options,
                                     pushFn = pushProjectToGithubFromDirs) {
  githubOptions <- list(
    github_org = githubOrg,
    github_user_email = githubEmail,
    github_user_name = githubUsername
  )
  do.call(optionsFn, githubOptions)

  pushArgs <- list(
    project_dirs = projectDirs,
    omic_type = omicType,
    experiment_label = experimentLabel,
    project_id = projectId
  )
  do.call(pushFn, pushArgs)

  list(githubOptions = githubOptions, pushArgs = pushArgs)
}

completeProtSummaryGithubPush <- function(output,
                                          projectDirs,
                                          omicType = "proteomics",
                                          experimentLabel,
                                          description,
                                          githubOrg,
                                          githubEmail,
                                          githubUsername,
                                          projectId,
                                          pushGithubFn = runProtSummaryGithubPush,
                                          renderTextFn = shiny::renderText,
                                          showNotificationFn = shiny::showNotification,
                                          timestampFn = Sys.time,
                                          logErrorFn = function(message) {
                                            logger::log_error(
                                              logger::skip_formatter(message)
                                            )
                                          }) {
  tryCatch({
    pushGithubFn(
      projectDirs = projectDirs,
      omicType = omicType,
      experimentLabel = experimentLabel,
      githubOrg = githubOrg,
      githubEmail = githubEmail,
      githubUsername = githubUsername,
      projectId = projectId
    )

    showNotificationFn("Successfully pushed to GitHub", type = "message")

    output$session_summary <- renderTextFn({
      paste("Workflow args created for:", experimentLabel,
            "\nDescription:", description,
            "\nTimestamp:", timestampFn(),
            "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]")
    })

    TRUE
  }, error = function(e) {
    showNotificationFn(
      paste("GitHub push failed:", e$message),
      type = "error",
      duration = 10
    )
    logErrorFn(paste("Failed to push to GitHub:", e$message))
    FALSE
  })
}

