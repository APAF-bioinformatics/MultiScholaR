# KIT-147: since contract v2 the integration input's reader returns no matrix; a data set's matrix is read through the
# contract's view reader (`openIntegrationInputView()`). These give the tests that compare written values the matrix
# whole, as the v1 reader did.

# One data set's whole matrix.
integrationInputMatrix <- function(dir, view, input = NULL) openIntegrationInputView(dir, view, input)$read()

# The reader's value with each view's `matrix` added.
readIntegrationInputWithMatrices <- function(dir) {
    input <- readIntegrationInput(dir)
    for (view in names(input$views)) input$views[[view]]$matrix <- integrationInputMatrix(dir, view, input)
    input
}
