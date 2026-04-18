# Keep the observer shell top-level so later waves can move it without
# reopening the run-normalization body itself.
runMetabNormNormalizationObserverShell <- function(
    workflowData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    runPipelineFn,
    logErrorFn = logger::log_error
) {
    reqFn(workflowData$state_manager)
    addLogFn("Starting normalization pipeline...")

    shellState <- withProgressFn(
        message = "Running normalization pipeline..."
        , value = 0
        , {
            tryCatch({
                pipelineState <- runPipelineFn()

                addLogFn("Normalization pipeline complete!")
                showNotificationFn(
                    "Normalization pipeline complete!"
                    , type = "message"
                )

                invisible(list(
                    status = "success"
                    , pipelineState = pipelineState
                ))
            }, error = function(e) {
                errorMessage <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

                addLogFn(paste("ERROR:", errorMessage))
                logErrorFn(paste("Normalization pipeline error:", errorMessage))
                showNotificationFn(paste("Error:", errorMessage), type = "error")

                invisible(list(
                    status = "error"
                    , errorMessage = errorMessage
                ))
            })
        }
    )

    invisible(shellState)
}

