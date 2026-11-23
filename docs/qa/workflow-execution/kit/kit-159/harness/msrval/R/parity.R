# The parity stage's comparisons (KIT-159 7.9, G13; VAL-14): the browser's record of the values it
# sent against the driver's ledger, the browser replaying the driver's recorded actions on the
# installed application, and the bundle driver's DA table against L1's.

#' The browser's own record of the values it has sent (`Shiny.shinyapp.$inputValues`), with type
#' suffixes after the colon removed, client-data entries dropped and the module prefix stripped
#'
#' @param js_values Named list from `AppDriver$get_js("Shiny.shinyapp.$inputValues")`.
#' @param prefix The omic module's input prefix.
#' @return Named list by input id relative to the module.
#' @export
valBrowserValues <- function(js_values, prefix) {
    ids <- sub(":.*$", "", names(js_values))
    keep <- startsWith(ids, prefix)
    out <- js_values[keep]
    names(out) <- substring(ids[keep], nchar(prefix) + 1L)
    out
}

#' Drive the installed application with the driver's recorded actions, recording before each
#' action the values the browser has sent (7.9 (c))
#'
#' @param app An `shinytest2::AppDriver` (or an object with its `set_inputs`, `upload_file`,
#'   `click`, `get_js` and `wait_for_idle` methods).
#' @param actions The driver's action log (`step, kind, input_id, value`).
#' @param prefix The module prefix in the application's input ids.
#' @return A data frame: `action`, `input_id`, `value` (the browser's ledger).
#' @export
valParityBrowser <- function(app, actions, prefix) {
    rows <- list()
    for (i in seq_len(nrow(actions))) {
        a <- actions[i, ]
        sent <- valBrowserValues(app$get_js("Shiny.shinyapp.$inputValues"), prefix)
        if (length(sent)) {
            rows[[length(rows) + 1L]] <- data.frame(
                action = a$input_id, input_id = names(sent),
                value = vapply(sent, function(v) {
                    if (is.null(v)) "NULL" else paste(format(unlist(v), digits = 15), collapse = "|")
                }, character(1), USE.NAMES = FALSE),
                stringsAsFactors = FALSE)
        }
        id <- paste0(prefix, a$input_id)
        switch(a$kind,
               upload = do.call(app$upload_file, stats::setNames(list(a$value), id)),
               set = do.call(app$set_inputs, stats::setNames(list(a$value), id)),
               press = app$click(input = id),
               valAbort("unknown action kind {a$kind}", "val_spec_invalid"))
        app$wait_for_idle()
    }
    if (length(rows)) do.call(rbind, rows) else
        data.frame(action = character(), input_id = character(), value = character(),
                   stringsAsFactors = FALSE)
}

#' Compare the driver's ledger with the browser's, per action and input (7.9 (d))
#'
#' Inputs in the committed exclusion list (client state) are reported with their class and not
#' compared. Any other difference, or an input on one side only, is `val_parity_differs`.
#'
#' @param driver The driver's ledger rows (`action, input_id, value`, ...).
#' @param browser The browser's rows (`action, input_id, value`).
#' @param exclusions Data frame `input_id, class`.
#' @return The comparison (`action, input_id, driver_value, browser_value, equal,
#'   excluded_class`), the content of `parity.tsv`.
#' @export
valCompareLedgers <- function(driver, browser, exclusions) {
    d <- unique(driver[, c("action", "input_id", "value")])
    b <- unique(browser[, c("action", "input_id", "value")])
    m <- merge(d, b, by = c("action", "input_id"), all = TRUE, suffixes = c("_driver", "_browser"))
    cmp <- data.frame(action = m$action, input_id = m$input_id, driver_value = m$value_driver,
                      browser_value = m$value_browser,
                      equal = !is.na(m$value_driver) & !is.na(m$value_browser) &
                          m$value_driver == m$value_browser,
                      excluded_class = exclusions$class[match(m$input_id, exclusions$input_id)],
                      stringsAsFactors = FALSE)
    bad <- cmp[!cmp$equal & is.na(cmp$excluded_class), , drop = FALSE]
    if (nrow(bad)) {
        valAbort("{nrow(bad)} input value{?s} differ between the driver and the browser",
                 "val_parity_differs")
    }
    cmp
}

#' The bundle driver's DA table against L1's for the same input (7.9 (e); A32): masked digest,
#' then the fingerprint's hash at 10 significant digits
#'
#' @param da The bundle driver's DA table.
#' @param l1 L1's committed `list(digest, fingerprint)` for that input.
#' @return A list: `status` (`identical`, `within_rounding`, `differs`), `digest`, `sha256_10`.
#' @export
valCompareParityDa <- function(da, l1) {
    d <- valMaskedDigest(da)
    f <- valFingerprint(da)
    status <- if (identical(d, l1$digest)) {
        "identical"
    } else if (identical(f$sha256_10, l1$fingerprint$sha256_10)) {
        "within_rounding"
    } else {
        "differs"
    }
    list(status = status, digest = d, sha256_10 = f$sha256_10)
}
