# The parity stage's comparisons (KIT-159 7.9, G13; VAL-14): the browser's own record of the
# values it sent against the driver's ledger, the browser replaying the driver's actions, and the
# bundle driver's DA table against L1's (digest, then fingerprints at 10 significant digits; A32).

test_that("the browser's inputValues lose their type suffixes", {
    js <- list(`proteomics_workflow-norm-ruv_k:shiny.number` = 3, `proteomics_workflow-x` = "a",
               `.clientdata_output_x_width` = 100)
    vals <- valBrowserValues(js, prefix = "proteomics_workflow-")
    expect_identical(vals, list(`norm-ruv_k` = 3, x = "a"))
})

test_that("the ledgers are compared per action and input, with excluded client-state inputs", {
    driver <- data.frame(step = "norm", action = "run_normalization",
                         input_id = c("norm-ruv_k", "norm-method", "norm-plot_zoom"),
                         value = c("3", "median", "1"), source = "ui_default", location = "ui",
                         stringsAsFactors = FALSE)
    browser <- data.frame(action = "run_normalization",
                          input_id = c("norm-ruv_k", "norm-method", "norm-plot_zoom"),
                          value = c("3", "median", "2"), stringsAsFactors = FALSE)
    exclusions <- data.frame(input_id = "norm-plot_zoom", class = "client_state")
    cmp <- valCompareLedgers(driver, browser, exclusions)
    expect_named(cmp, c("action", "input_id", "driver_value", "browser_value", "equal",
                        "excluded_class"))
    expect_true(all(cmp$equal[is.na(cmp$excluded_class)]))
    expect_identical(cmp$excluded_class[cmp$input_id == "norm-plot_zoom"], "client_state")
    browser$value[[1L]] <- "4"
    expect_error(valCompareLedgers(driver, browser, exclusions), class = "val_parity_differs")
    only_driver <- browser[-2L, ]
    only_driver$value[[1L]] <- "3"
    expect_error(valCompareLedgers(driver, only_driver, exclusions), class = "val_parity_differs")
})

test_that("the browser replays the driver's actions in order and records its values before each", {
    calls <- character()
    sent <- list(`app-norm-ruv_k:shiny.number` = 3)
    app <- list(
        set_inputs = function(...) {
            calls <<- c(calls, paste("set", names(list(...))))
            invisible(NULL)
        },
        upload_file = function(...) {
            calls <<- c(calls, paste("upload", names(list(...))))
            invisible(NULL)
        },
        click = function(input) {
            calls <<- c(calls, paste("click", input))
            invisible(NULL)
        },
        get_js = function(script) sent,
        wait_for_idle = function(...) invisible(NULL)
    )
    actions <- data.frame(step = c("import", "import", "norm"), kind = c("upload", "set", "press"),
                          input_id = c("import-file", "import-format", "norm-run"),
                          value = c("/x/report.tsv", "diann", "1"), stringsAsFactors = FALSE)
    rec <- valParityBrowser(app, actions, prefix = "app-")
    expect_identical(calls, c("upload app-import-file", "set app-import-format", "click app-norm-run"))
    expect_named(rec, c("action", "input_id", "value"))
    expect_identical(unique(rec$action), c("import-file", "import-format", "norm-run"))
    expect_identical(rec$value[rec$input_id == "norm-ruv_k"][[1L]], "3")
})

test_that("the bundle DA table against L1's: digest, then the values under SWD-083 (VAL-133 (1))", {
    da <- data.frame(feature_id = c("a", "b"), log2fc = c(1.123456789012, 2), p = c(0.1, 0.2))
    l1 <- list(digest = valMaskedDigest(da), fingerprint = valFingerprint(da), table = da,
               keys = "feature_id", q_col = "p", fc_col = "log2fc")
    expect_identical(valCompareParityDa(da, l1)$status, "identical")
    near <- da
    near$log2fc[[1L]] <- near$log2fc[[1L]] * (1 + 1e-14)
    expect_identical(valCompareParityDa(near, l1)$status, "within_tolerance")
    far <- da
    far$log2fc[[1L]] <- 1.2
    expect_identical(valCompareParityDa(far, l1)$status, "differs")
    no_table <- l1[c("digest", "fingerprint")]
    expect_error(valCompareParityDa(near, no_table), class = "val_spec_invalid")
})
