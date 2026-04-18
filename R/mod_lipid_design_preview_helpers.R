formatLipidDesignAssaysPreview <- function(dataTbl) {
    assayNames <- names(dataTbl)
    paste("Included assays:", paste(assayNames, collapse = ", "))
}

registerLipidDesignPreviewOutputs <- function(
    output,
    workflowData,
    renderDt = DT::renderDT,
    renderText = shiny::renderText,
    req = shiny::req,
    assayPreviewFormatter = formatLipidDesignAssaysPreview
) {
    output$design_matrix_preview <- renderDt({
        req(workflowData$design_matrix)
        workflowData$design_matrix
    }, options = list(pageLength = 5, scrollX = TRUE))

    output$contrasts_preview <- renderDt({
        req(workflowData$contrasts_tbl)
        workflowData$contrasts_tbl
    }, options = list(pageLength = 5, scrollX = TRUE))

    output$assays_preview <- renderText({
        req(workflowData$data_tbl)
        assayPreviewFormatter(workflowData$data_tbl)
    })

    invisible(output)
}

