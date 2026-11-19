# KIT-129 item 5 (the KIT-145 observation): with no session loaded, the lipidomics and metabolomics DA tab falls back to
# the owner's head only when the normalisation tab is complete; otherwise it refuses with the group's `no_data` text and
# runs nothing. Transcriptomics and the other omics are as before.

omicDaFallbackRun <- function(omic_type, normalization, class = "LipidomicsAssayData") {
    calls <- new.env()
    calls$notifications <- character()
    calls$runs <- 0L
    testthat::local_mocked_bindings(
        omicDaLoadedMetadata = function(...) NULL,
        omicDaReadSource = function(...) list(class = class),
        executeOmicDaRun = function(...) {
            calls$runs <- calls$runs + 1L
            list(status = "success")
        },
        .env = environment(runOmicDaAnalysis))
    testthat::local_mocked_bindings(showNotification = function(ui, ...) {
        calls$notifications <- c(calls$notifications, as.character(ui))
        invisible(NULL)
    }, removeNotification = function(...) invisible(NULL), .package = "shiny")
    workflowData <- list(contrasts_tbl = data.frame(contrasts = "groupB-groupA"),
        tab_status = list(quality_control = "complete", normalization = normalization))
    da <- list(source = list(class = class, no_data = "No data loaded."), model = list(editable = FALSE), preflight_first = FALSE,
        preflight = NULL)
    value <- runOmicDaAnalysis(list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), NULL, workflowData,
        list(current_state_name = NULL), list(), da, list(omic_type = omic_type))
    list(value = value, runs = calls$runs, notifications = calls$notifications)
}

test_that("lipidomics and metabolomics refuse the head fallback until the normalisation tab is complete", {
    for (omic in list(c("lipidomics", "LipidomicsAssayData"), c("metabolomics", "MetaboliteAssayData"))) {
        refused <- omicDaFallbackRun(omic[[1L]], "pending", omic[[2L]])
        expect_null(refused$value, info = omic[[1L]])
        expect_identical(refused$runs, 0L, info = omic[[1L]])
        expect_identical(refused$notifications, "No data loaded.", info = omic[[1L]])
        allowed <- omicDaFallbackRun(omic[[1L]], "complete", omic[[2L]])
        expect_identical(allowed$runs, 1L, info = omic[[1L]])
    }
})

test_that("the other omics keep the head fallback whatever the normalisation tab says", {
    kept <- omicDaFallbackRun("transcriptomics", "pending", "TranscriptomicsAssayData")
    expect_identical(kept$runs, 1L)
})
