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

# KIT-184 K4 (gate G1): the removed factors of a written data set, checked against its two Parquet matrices read back
# and not against the record alone: the record's status and number of factors; the number of features with a finite
# value in every run of both matrices; and, for each of those features, the row-centred difference (normalised minus
# corrected) projected off the written factors, whose residual must be below 1e-8 of the feature's difference (of the
# largest difference times 1e-6 for a feature from which almost nothing was removed) when the factors are exact. The
# measured s_(k+1)/s_1 of that difference, the share and the largest residual are printed for the ticket's record.
# Returns the written factors (runs x k), invisibly.
expectIntegrationRemovedFactors <- function(dir, view, k, status = "exact", input = readIntegrationInput(dir)) {
    k <- as.integer(k)
    entry <- Filter(function(item) identical(item$name, view), input$manifest$views)[[1L]]
    record <- entry$removed_factors
    testthat::expect_identical(record$status, status, info = view)
    testthat::expect_identical(as.integer(record$k_declared), k, info = view)
    testthat::expect_identical(as.integer(record$k), k, info = view)
    testthat::expect_identical(entry$removed_factors_file, sprintf("removed_%s.parquet", view), info = view)
    normalised <- integrationInputMatrix(dir, view, input)
    corrected <- integrationInputMatrix(dir, view, input, "ruv_corrected")
    runs <- names(normalised)[-1L]
    u <- as.matrix(normalised[runs])
    v <- as.matrix(corrected[runs])
    complete <- rowSums(!is.finite(u)) == 0L & rowSums(!is.finite(v)) == 0L
    testthat::expect_identical(as.integer(record$features), sum(complete), info = view)
    d <- (u - v)[complete, , drop = FALSE]
    d <- d - rowMeans(d)
    s <- svd(d, nu = 0L, nv = 0L)$d
    ratio <- if (length(s) > k) s[[k + 1L]] / s[[1L]] else 0
    factors <- readIntegrationInputRemovedFactors(dir, view, input)
    testthat::expect_identical(dim(factors), c(length(runs), k), info = view)
    testthat::expect_identical(rownames(factors), runs, info = view)
    norms <- sqrt(rowSums(d^2))
    residual <- sqrt(rowSums((d - d %*% factors %*% t(factors))^2))
    worst <- max(residual / pmax(norms, 1e-6 * max(norms)))
    cat(sprintf(paste0("\nKIT-184 G1 %s: status %s, k declared %d, written %d, rank %s, features %d of %d, s_(k+1)/s_1 %.3g, ",
        "share %.9g, largest relative residual %.3g\n"), view, record$status, k, as.integer(record$k),
        format(record$rank), sum(complete), nrow(u), ratio, as.numeric(record$share), worst))
    if (identical(status, "exact")) {
        testthat::expect_lt(worst, 1e-8, label = sprintf("largest relative residual off the written factors of %s (%g)", view, worst))
        testthat::expect_lte(ratio, 1e-8, label = sprintf("s_(k+1)/s_1 of %s (%g)", view, ratio))
    }
    invisible(factors)
}
