# KIT-147: since contract v2 the integration input's reader returns no matrix; a data set's matrix is read through the
# contract's view reader (`openIntegrationInputView()`). These give the tests that compare written values the matrix
# whole, as the v1 reader did.

# One data set's whole matrix (KIT-153: `matrix = "ruv_corrected"` reads its RUV-corrected matrix).
integrationInputMatrix <- function(dir, view, input = NULL, matrix = "uncorrected") {
    openIntegrationInputView(dir, view, input, matrix)$read()
}

# The reader's value with each view's `matrix` added (KIT-153: and `corrected_matrix` where the view has one).
readIntegrationInputWithMatrices <- function(dir) {
    input <- readIntegrationInput(dir)
    for (view in names(input$views)) {
        input$views[[view]]$matrix <- integrationInputMatrix(dir, view, input)
        if (isTRUE(input$views[[view]]$corrected)) {
            input$views[[view]]$corrected_matrix <- integrationInputMatrix(dir, view, input, "ruv_corrected")
        }
    }
    input
}
