# KIT-149: MOFA2 trained through the partitioned alignment against the whole-matrix alignment on the KIT-137 tier-1
# corpus, each variant on its own store. Skipped where MOFA2 or a Python with mofapy2 is absent; run on the MOFA2 lane
# (its own file, run first, so that it runs in a process where nothing else has trained).

intPartMofaCorpus <- function() file.path(normalizePath(testthat::test_path("..", ".."), mustWork = FALSE), "tests", "testdata",
    "integration", "corpus")
intPartMofaOmics <- c("proteomics", "transcriptomics", "lipidomics", "metabolomics")
intPartMofaWarnings <- function(messages) messages[startsWith(messages, "WARNING:")]

# The same data handed over, the same model, value for value (identical, one platform); the stored-result test of
# test-integration-mofa-model.R runs the default variant.
test_that("KIT-149: MOFA2 trained through the partitioned alignment gives the whole-matrix alignment's model, value for value", {
    testthat::skip_if_not(dir.exists(intPartMofaCorpus()), "the tier-1 corpus is not present")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
    trained <- lapply(c(whole = "whole", partitioned = "partitioned"), function(variant) {
        manager <- newIntegrationStateOwner(file.path(tempfile(paste0("kit149-", variant, "-")), "project"))
        on.exit(closeWorkflowStateOwner(manager), add = TRUE)
        importIntegrationInputs(manager, file.path(intPartMofaCorpus(), intPartMofaOmics))
        run <- runIntegrationModel(manager, list(variant = variant))
        value <- run$value
        aligned <- lapply(stats::setNames(names(value$aligned_values), names(value$aligned_values)), function(view) {
            integrationModelViewTable(manager, value, "aligned_values", view)
        })
        comparable <- integrationModelComparable(manager, value)
        c(list(status = run$status, aligned_values = aligned, warnings = intPartMofaWarnings(value$messages)),
            comparable[c("factors", "weights", "variance_explained", "variance_total", "samples", "view_summary", "elbo",
                "versions")])
    })
    expect_identical(trained$whole$status, "computed")
    differences <- vapply(c("factors", "weights", "variance_explained"), function(name) {
        a <- trained$whole[[name]]
        b <- trained$partitioned[[name]]
        numeric <- names(a)[vapply(a, is.numeric, logical(1))]
        if (!identical(dim(a), dim(b))) return(Inf)
        max(abs(as.matrix(a[numeric]) - as.matrix(b[numeric])))
    }, numeric(1))
    message(sprintf("KIT-149 partitioned against whole: largest absolute difference %s; identical %s",
        paste(sprintf("%s %.3g", names(differences), differences), collapse = ", "),
        paste(vapply(names(trained$whole), function(name) identical(trained$partitioned[[name]], trained$whole[[name]]),
            logical(1)), collapse = ", ")))
    for (name in names(trained$whole)) expect_identical(trained$partitioned[[name]], trained$whole[[name]], info = name)
})
