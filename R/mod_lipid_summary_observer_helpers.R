initializeLipidSummarySessionBootstrap <- function(
    session,
    experimentLabel = NULL,
    updateTextInputFn = shiny::updateTextInput,
    reactiveValuesFn = shiny::reactiveValues
) {
    if (!is.null(experimentLabel)) {
        updateTextInputFn(session, "experiment_label", value = experimentLabel)
    }

    reactiveValuesFn(
        workflow_args_saved = FALSE,
        files_copied = FALSE,
        report_generated = FALSE,
        report_path = NULL
    )
}

registerLipidSummaryInitialOutputs <- function(
    output,
    renderTextFn = shiny::renderText,
    reactiveFn = shiny::reactive,
    outputOptionsFn = shiny::outputOptions
) {
    output$session_summary <- renderTextFn({
        "Ready to save workflow parameters and generate report"
    })

    output$report_ready <- reactiveFn({ FALSE })
    outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)

    output
}
