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
                action = a$input_id, input_id = names(sent),  # nolint
                value = vapply(sent, function(v) {
                    if (is.null(v)) {
                        "NULL"
                    } else {
                        paste(format(unlist(v), digits = 15), collapse = "|")
                    }
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

#' The bundle driver's DA table against L1's for the same input (7.9 (e); A32; VAL-133 (1))
#'
#' Equal masked digests: `identical`. Otherwise the values against L1's table under SWD-083
#' ([valCompareDaValues()]: keys, order, missing values, calls and text columns exactly; numbers
#' within 1e-6 relative or 1e-9 absolute): `within_tolerance`, else `differs`. Rounding to 10
#' significant digits no longer decides (VAL-133: it does not absorb the measured DIA-with-limpa
#' difference between the two maths libraries); `sha256_10` is still reported.
#'
#' @param da The bundle driver's DA table.
#' @param l1 L1's committed `list(digest, fingerprint, table, keys, q_col, fc_col)` for that input.
#' @return A list: `status` (`identical`, `within_tolerance`, `differs`), `digest`, `sha256_10`,
#'   `comparison`.
#' @export
valCompareParityDa <- function(da, l1) {
    d <- valMaskedDigest(da)
    f <- valFingerprint(da)
    if (identical(d, l1$digest)) {
        return(list(status = "identical", digest = d, sha256_10 = f$sha256_10, comparison = NULL))
    }
    if (is.null(l1$table) || is.null(l1$keys) || is.null(l1$q_col) || is.null(l1$fc_col)) {
        valAbort("L1's table, keys, q and fold-change columns are needed to compare values",
                 "val_spec_invalid")
    }
    cmp <- valCompareDaValues(da, l1$table, keys = l1$keys, q_col = l1$q_col,
                              fc_col = l1$fc_col, mode = "swd083")
    list(status = if (identical(cmp$verdict, "outside")) "differs" else "within_tolerance",
         digest = d, sha256_10 = f$sha256_10, comparison = cmp)
}

#' The Linux bundle's GUI-entry lane steps excluded by name from any installer acceptance
#' (VAL-133; KIT-199 implementation audit O4): the DIA-with-limpa DA table is written with one
#' value in another notation under OpenBLAS, so these six steps stay red on the bundle permanently
#' @noRd
PARITY_EXCLUDED_STEPS <- c("da", "summary_refusal", "summary_save", "summary_copy",
                           "summary_report", "summary_export")
PARITY_EXCLUDED_REASON <- paste("VAL-133: the bundle's DIA-with-limpa DA table differs in one",
                                "value's notation (OpenBLAS); excluded by name from installer",
                                "acceptance")
PARITY_INSTALLER_EXCLUDED <- data.frame(lane = "prot_dia_limpa", step = PARITY_EXCLUDED_STEPS,
                                        reason = PARITY_EXCLUDED_REASON,
                                        stringsAsFactors = FALSE)

#' Rows of a bundle GUI-entry lane summary that an installer acceptance may read (VAL-133)
#'
#' @param rows The lane summary (`lane`, `step`, ...).
#' @return A list: `kept` (rows outside the exclusion), `excluded` (the excluded rows with their
#'   `reason`).
#' @export
valInstallerAcceptanceRows <- function(rows) {
    k <- match(paste(rows$lane, rows$step, sep = "\r"),
               paste(PARITY_INSTALLER_EXCLUDED$lane, PARITY_INSTALLER_EXCLUDED$step, sep = "\r"))
    excluded <- rows[!is.na(k), , drop = FALSE]
    excluded$reason <- PARITY_INSTALLER_EXCLUDED$reason[k[!is.na(k)]]
    kept <- rows[is.na(k), , drop = FALSE]
    rownames(kept) <- NULL
    rownames(excluded) <- NULL
    list(kept = kept, excluded = excluded)
}

PARITY_TWO_SPECIES <- "prot_dia_two_species"
PARITY_WINDOWS_EXCLUDED <- c("prot_dia", "prot_dia_limpa", PARITY_TWO_SPECIES)

#' Read and check a parity request (7.9 item 1; `msrval.parity_request.v1`)
#'
#' Routes are the ten of section 7.4 and the two-species route; Windows excludes the DIA routes
#' (KIT-134 D4). Commits are 40 hex digits; each declared generated input carries a SHA-256.
#'
#' @param path The request file.
#' @return The request (a list).
#' @export
valParityRequest <- function(path) {
    req <- jsonlite::read_json(path, simplifyVector = FALSE)
    bad <- function(what) valAbort("parity request: {what}", "val_spec_invalid")
    if (!identical(req[["schema"]], "msrval.parity_request.v1")) bad("schema")
    if (!req[["platform"]] %in% c("linux", "windows")) bad("platform")
    routes <- unlist(req[["routes"]])
    if (length(routes) == 0L || !all(routes %in% c(valRouteTable()$route, PARITY_TWO_SPECIES))) {
        bad("an unknown route")
    }
    if (identical(req[["platform"]], "windows") && any(routes %in% PARITY_WINDOWS_EXCLUDED)) {
        bad("a DIA route on Windows")
    }
    for (f in c("release_commit", "harness_commit", "harness_tree")) {
        if (!grepl("^[0-9a-f]{40}$", req[[f]] %||% "")) bad(f)
    }
    for (x in req[["inputs"]]) {
        if (!grepl("^[0-9a-f]{64}$", x[["sha256"]] %||% "") || !nzchar(x[["path"]] %||% "")) {
            bad("an input without a path or SHA-256")
        }
    }
    req
}

#' The installed bundle's commit against the release commit (7.9 item 3; KIT-176 A02)
#'
#' Accepted when `VERSION.json`'s `package_commit` equals the release commit, or when
#' `git diff --quiet` finds no difference between them under `pkg_paths` (KIT-176's `PKG_PATHS`,
#' read from `packaging/package-commit.sh` by the caller, never a copy).
#'
#' @param version The bundle's `VERSION.json` (a list).
#' @param release_commit The request's release commit.
#' @param pkg_paths KIT-176's `PKG_PATHS`, split into paths.
#' @param repo The checkout (full history).
#' @return `equal` or `same_package_paths`; otherwise `val_product_tree_mismatch`.
#' @export
valCheckBundleCommit <- function(version, release_commit, pkg_paths, repo = ".") {
    pc <- version[["package_commit"]] %||% ""
    if (identical(pc, release_commit)) return("equal")
    st <- valGit(repo, "diff", "--quiet", pc, release_commit, "--", pkg_paths)$status
    if (identical(as.integer(st), 0L)) return("same_package_paths")
    valAbort("the bundle's package commit {pc} differs from the release commit under PKG_PATHS",
             "val_product_tree_mismatch")
}

#' Check the workflow's guard for the parity stage (7.9 items 1, 5)
#'
#' The stage `parity` is offered; outside the snapshot job exactly one upload exists in each of the
#' `linux` and `windows` jobs, each a step named `parity-results` whose `if:` requires stage
#' `parity` and `steps.parity-check.outcome == 'success'` and whose `path:` is `val-out/parity/`;
#' each of those jobs has a step `parity-check`.
#'
#' @param path `.github/workflows/distributable.yml`.
#' @return `TRUE`, or `val_spec_invalid`.
#' @export
valCheckParityWorkflow <- function(path) {
    lines <- readLines(path, warn = FALSE)
    bad <- function(what) valAbort("distributable.yml: {what}", "val_spec_invalid")
    if (!any(grepl("options: \\[.*\\bparity\\b.*\\]", lines, perl = TRUE))) {
        bad("no stage parity")
    }
    job <- character(length(lines))
    cur <- ""
    for (i in seq_along(lines)) {
        if (grepl("^  [a-z][a-z0-9-]*:\\s*$", lines[[i]])) {
            cur <- sub(":\\s*$", "", trimws(lines[[i]]))
        }
        job[[i]] <- cur
    }
    starts <- grep("^\\s*- (name|uses|id):", lines)
    step_of <- function(i) {
        s <- max(starts[starts <= i])
        e <- min(c(starts[starts > s], length(lines) + 1L)) - 1L
        lines[s:e]
    }
    ups <- grep("uses: *actions/upload-artifact", lines)
    outside <- ups[job[ups] != "snapshot"]
    for (j in c("linux", "windows")) {
        u <- outside[job[outside] == j]
        if (length(u) != 1L) bad(sprintf("%d uploads in job %s (one admitted)", length(u), j))
        step <- step_of(u)
        if (!any(grepl("name: *parity-results\\s*$", step))) {
            bad(paste("the", j, "upload is not parity-results"))
        }
        cond <- step[grepl("^\\s*if:", step)]
        guard <- "steps.parity-check.outcome == 'success'"
        if (length(cond) != 1L || !grepl(guard, cond, fixed = TRUE) ||
                !grepl("parity", sub(guard, "", cond, fixed = TRUE))) {
            bad(paste("the", j, "upload's if: lacks stage parity or the parity-check condition"))
        }
        if (!any(grepl("path: *val-out/parity/\\s*$", step))) bad(paste("the", j, "upload's path"))
        if (!any(job == j & grepl("id: *parity-check\\s*$", lines))) {
            bad(paste("no parity-check in", j))
        }
    }
    if (length(outside) != 2L) bad("an upload outside the snapshot job and the two parity-results")
    TRUE
}

#' Check a parity request's copied L1 values against L1's committed evidence (7.9 item 1, r1.12)
#'
#' The request names `l1_path`, L1's committed digests and fingerprints per route; the values the
#' request copies (`l1`) must equal that file at the given commit.
#'
#' @param req The request (from [valParityRequest()]).
#' @param repo The checkout (full history).
#' @param at The commit to read the evidence at (default the request's release commit).
#' @return `TRUE`, or `val_spec_invalid`.
#' @export
valCheckParityL1 <- function(req, repo = ".", at = req[["release_commit"]]) {
    bad <- function(what) valAbort("parity request: {what}", "val_spec_invalid")
    path <- req[["l1_path"]] %||% ""
    if (!nzchar(path)) bad("no l1_path")
    text <- valGitShow(repo, at, path)
    if (length(text) == 0L || !nzchar(paste(text, collapse = ""))) bad("{path} is not committed")
    committed <- jsonlite::fromJSON(paste(text, collapse = "\n"), simplifyVector = FALSE)
    for (route in names(req[["l1"]])) {
        got <- req[["l1"]][[route]]
        want <- committed[[route]]
        same_fp <- identical(got[["fingerprint"]][["sha256_10"]],
                             want[["fingerprint"]][["sha256_10"]])
        if (is.null(want) || !identical(got[["digest"]], want[["digest"]]) || !same_fp) {
            bad(paste("the copied L1 values of", route, "differ from", path))
        }
    }
    TRUE
}
