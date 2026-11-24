# G1(a)'s control-record loader, the release-commit value comparison and the part (1) and (2)
# runners (KIT-159 section 2; KIT-199 12.3; VAL-120 (2), VAL-133; I15 of the pre-release audit).

#' G1(a)'s reference point: run 44953187101's base
#' @noRd
G1A_SINCE <- "12d68ca78442a50b46f52e33525dd85cf86df5e7"

#' The DIA-with-limpa DA artifact the release commit compares by value (KIT-199 12.3, D7)
#' @noRd
G1A_DA_KEY <- "results/proteomics/da_proteins/da_proteins_groupKO-groupWT_long_annot.tsv"

#' The environment keys D7's exact mode requires equal to the link's control jobs (KIT-199 I04)
#' @noRd
G1A_D7_ENV_KEYS <- c("image_version", "c_library_version", "libblas3_package",
                     "liblapack3_package", "blas_library", "lapack_library", "r_version",
                     "rust_toolchain_version", "pin_fastduckdb", "pin_fastlimpadpc")

#' SWD-083: doubles agree within 1e-6 relative or 1e-9 absolute (`compare-results.R`)
#' @noRd
SWD083_REL <- 1e-6
SWD083_ABS <- 1e-9

#' The app's q threshold for the calls (`R/func_prot_da_model_methods.R:597-608`)
#' @noRd
DA_Q_THRESHOLD <- 0.05

#' Lines of a file, or of text already read
#' @noRd
valLinesOf <- function(x) {
    if (length(x) == 1L && !grepl("\n", x, fixed = TRUE) && file.exists(x)) {
        return(readLines(x, warn = FALSE))
    }
    unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE))
}

#' Read a tab-separated table as text, cell for cell (no type conversion, no NA strings), skipping
#' comment lines that start with `#`
#'
#' @param x A path, or the file's text.
#' @return A data frame of character columns; attribute `comments` holds the comment lines.
#' @export
valReadTextTable <- function(x) {
    lines <- valLinesOf(x)
    comments <- grep("^#", lines, value = TRUE)
    lines <- lines[!grepl("^#", lines)]
    lines <- lines[cumsum(nzchar(lines)) > 0L]
    if (length(lines) == 0L) valAbort("an empty table", "val_spec_invalid")
    split_row <- function(l) {
        f <- strsplit(paste0(l, "\tEND"), "\t", fixed = TRUE)[[1L]]
        f[-length(f)]
    }
    header <- split_row(lines[[1L]])
    # A cell may hold line breaks written raw (KIT-172's `declared-keys.tsv`, DA summary texts):
    # a line with fewer fields than the header continues on the next line.
    rows <- character()
    buf <- NULL
    for (l in lines[-1L]) {
        buf <- if (is.null(buf)) l else paste(buf, l, sep = "\n")
        if (length(split_row(buf)) >= length(header)) {
            rows <- c(rows, buf)
            buf <- NULL
        }
    }
    if (!is.null(buf) && nzchar(buf)) rows <- c(rows, buf)
    fields <- c(list(header), lapply(rows, split_row))
    if (any(lengths(fields) != length(header))) {
        valAbort("a table row with another field count than its header", "val_spec_invalid")
    }
    body <- fields[-1L]
    out <- as.data.frame(stats::setNames(lapply(seq_along(header), function(j) {
        vapply(body, function(row) row[[j]], character(1))
    }), header), stringsAsFactors = FALSE, check.names = FALSE)
    attr(out, "comments") <- comments
    out
}

#' The classes of a ticket's own change in a record: every `; `-separated part tagged `(<ticket>)`
#' @noRd
valTicketClasses <- function(classes, ticket) {
    u <- unique(classes[!is.na(classes)])
    tag <- sprintf("(%s)", ticket)
    u[vapply(strsplit(u, "; ", fixed = TRUE), function(p) all(endsWith(p, tag)), logical(1))]
}

#' Build one lane's chain link from a ticket's control record (G1(a); KIT-199 12.3)
#'
#' Reads `gui-entry-control.tsv` in either committed column form: `a`/`b` (KIT-177, 180, 199: `a`
#' the earlier side) or `runner`/`oracle` (KIT-172: `oracle` the earlier side). The comparison of
#' the stored file against the base (`frozen_vs_base`, KIT-172 `frozen-base`) gives the base
#' differences (`base`, `live` values kept for link (i)); the base against the head
#' (`base_vs_head`; KIT-172 `base-head0`, the head at the shipped defaults, its declared-changes
#' file) gives the differences before and from the changed step. Other comparisons (KIT-172's
#' `base-head1`, `head1-head0`) are not part of the link. The changed step is the record's
#' `changed_step` column, or its comment line; a row is before it by the record's
#' `before_changed_step` or `phase` column, else by `lane_steps`. `declared-keys.tsv` gives
#' `declared` (keys with a leading `outputs/`, `states/` or `artifacts/` field are read without
#' it, as the record writes them).
#'
#' @param control The record (path or text).
#' @param declared `declared-keys.tsv` (path or text), or `NULL`.
#' @param ticket `KIT-<n>`.
#' @param lane The lane.
#' @param lane_steps The lane's steps in order (needed when the record has no position column).
#' @param change_classes The classes of the ticket's change; by default every class of the record
#'   whose parts are all tagged with the ticket.
#' @return A link for [valReadChain()].
#' @export
valReadControlLink <- function(control, declared = NULL, ticket, lane, lane_steps = NULL,
                               change_classes = NULL) {
    tab <- valReadTextTable(control)
    if (all(c("a", "b") %in% names(tab))) {
        earlier <- tab$a
        later <- tab$b
    } else if (all(c("runner", "oracle") %in% names(tab))) {
        earlier <- tab$oracle
        later <- tab$runner
    } else {
        valAbort("{ticket}'s control record has neither a/b nor runner/oracle", "val_spec_invalid")
    }
    tab$earlier <- ifelse(earlier == "<absent>", NA_character_, earlier)
    tab$later <- ifelse(later == "<absent>", NA_character_, later)
    tab <- tab[tab$lane == lane, , drop = FALSE]
    cmps <- unique(tab$comparison)
    base_cmp <- intersect(c("frozen_vs_base", "frozen-base"), cmps)
    head_cmp <- intersect(c("base_vs_head", "base-head0"), cmps)[1L]
    changed <- unique(tab$changed_step[!is.na(tab$changed_step %||% NA)])
    if (length(changed) == 0L) {
        com <- attr(tab, "comments") %||% attr(valReadTextTable(control), "comments")
        hit <- regmatches(com, regexpr("changed step per lane: [A-Za-z_]+", com))
        changed <- unique(sub("changed step per lane: ", "", hit, fixed = TRUE))
    }
    if (length(changed) != 1L) {
        valAbort("{ticket}'s record names no single changed step for {lane}", "val_spec_invalid")
    }
    bd <- tab[tab$comparison %in% base_cmp, , drop = FALSE]
    hd <- tab[tab$comparison %in% head_cmp, , drop = FALSE]
    before <- if (!is.null(hd$before_changed_step)) {
        toupper(hd$before_changed_step) == "TRUE"
    } else if (!is.null(hd$phase)) {
        startsWith(hd$phase, "before_")
    } else {
        if (is.null(lane_steps)) valAbort("no step order for {lane}", "val_spec_invalid")
        match(hd$step, lane_steps) < match(changed, lane_steps)
    }
    before <- before %in% TRUE
    dk <- data.frame(lane = character(), step = character(), key = character(),
                     expected = character(), stringsAsFactors = FALSE)
    if (!is.null(declared)) {
        dk <- valReadTextTable(declared)
        mine <- dk$lane == lane
        if (!is.null(dk$ticket)) mine <- mine & dk$ticket %in% ticket
        dk <- dk[mine, , drop = FALSE]
        dk$key <- sub("^(outputs|states|artifacts)/", "", dk$key)
    }
    cls <- c(bd$class, hd$class)
    list(ticket = ticket, lane = lane, changed_step = changed,
         base_diffs = data.frame(step = bd$step, field = bd$field, key = bd$key, class = bd$class,
                                 base = bd$earlier, live = bd$later, stringsAsFactors = FALSE),
         before_diffs = data.frame(step = hd$step[before], field = hd$field[before],
                                   key = hd$key[before], class = hd$class[before],
                                   head_value = hd$later[before], stringsAsFactors = FALSE),
         after_diffs = data.frame(step = hd$step[!before], field = hd$field[!before],
                                  key = hd$key[!before], class = hd$class[!before],
                                  stringsAsFactors = FALSE),
         declared = dk,
         change_classes = change_classes %||% valTicketClasses(cls, ticket),
         comparisons = list(base = base_cmp, head = head_cmp,
                            ignored = setdiff(cmps, c(base_cmp, head_cmp))))
}

#' The key set of a reference row, as `field/key` strings
#'
#' Reads both committed forms: `fields`/`keys` (KIT-172, G1(r)'s own: keys `;`-separated, each
#' `field/key`, or bare when one field is named) and `differing_fields_and_keys` (KIT-180, 199:
#' `field: head only k,k; frozen only k,k` or `field: <key> frozen <d> live <d>`).
#' @noRd
valReferenceKeys <- function(row) {
    if (!is.null(row$keys)) {
        fields <- strsplit(row$fields %||% "", ",", fixed = TRUE)[[1L]]
        k <- strsplit(row$keys, "[;,]")[[1L]]
        k <- trimws(k[nzchar(trimws(k))])
        known <- paste0(unique(c(fields, "artifacts", "outputs", "states")), "/")
        prefixed <- vapply(k, function(x) any(startsWith(x, known)), logical(1))
        if (length(fields) == 1L) k[!prefixed] <- paste0(fields, "/", k[!prefixed])
        return(sort(unique(k)))
    }
    out <- character()
    for (seg in strsplit(row$differing_fields_and_keys, "; (?=[a-z_]+: )", perl = TRUE)[[1L]]) {
        field <- sub(":.*$", "", seg)
        rest <- sub("^[a-z_]+: ", "", seg)
        if (grepl("^(head|frozen) only ", rest)) {
            parts <- strsplit(rest, "; ", fixed = TRUE)[[1L]]
            k <- unlist(strsplit(sub("^(head|frozen) only ", "", parts), ",", fixed = TRUE))
        } else {
            k <- sub(" frozen .*$", "", rest)
        }
        out <- c(out, paste0(field, "/", trimws(k)))
    }
    sort(unique(out))
}

#' Read committed reference rows (`g1a-reference*.tsv`, `k156-gui-entry-failures.tsv`)
#'
#' @param x A path, or the file's text.
#' @return `lane, step, class`, `keyset` (a list of `field/key` sets) and `truncated` (the source
#'   listed only part of the names).
#' @export
valReadReferenceRows <- function(x) {
    t <- valReadTextTable(x)
    out <- data.frame(lane = t$lane, step = t$step, class = t$class, stringsAsFactors = FALSE)
    out$keyset <- lapply(seq_len(nrow(t)), function(i) valReferenceKeys(t[i, , drop = FALSE]))
    out$truncated <- if (is.null(t$names_list_truncated)) rep(FALSE, nrow(t)) else
        t$names_list_truncated %in% "yes"
    out
}

#' Compare a DA table with a reference table by value (KIT-199 12.3 D7; VAL-133; G13)
#'
#' Exactly: the header, the keys and their order, the missing values, the calls (`q < 0.05`, up
#' when the fold change is at least 0, else down) and every non-numeric column. Numbers: under
#' SWD-083 (`mode = "swd083"`) or equal (`mode = "exact"`). A column is numeric when every
#' non-missing cell of both tables reads as a number.
#'
#' @param a,b The two tables (text as written, or numeric columns).
#' @param keys,q_col,fc_col Key, q-value and fold-change columns.
#' @param mode `swd083` or `exact`.
#' @return A list: `verdict` (`identical`, `agree`, `outside`), the exact checks,
#'   `calls_changed`, and `columns` (per numeric column: compared, differ, max_abs, max_rel,
#'   outside).
#' @export
valCompareDaValues <- function(a, b, keys = c("Protein.Ids", "comparison"), q_col = "fdr_qvalue",
                               fc_col = "log2FC", mode = c("swd083", "exact")) {
    mode <- match.arg(mode)
    a <- as.data.frame(a, stringsAsFactors = FALSE)
    b <- as.data.frame(b, stringsAsFactors = FALSE)
    need <- c(keys, q_col, fc_col)
    if (!all(need %in% names(a)) || !all(need %in% names(b))) {
        valAbort("the tables lack a key, q or fold-change column", "val_spec_invalid")
    }
    num <- function(x) {
        if (is.numeric(x)) return(as.numeric(x))
        x <- as.character(x)
        x[x %in% c("", "NA", "NaN")] <- NA_character_
        suppressWarnings(as.numeric(x))
    }
    isMissing <- function(x) {
        if (is.numeric(x)) is.na(x) else as.character(x) %in% c("", "NA", "NaN")
    }
    keyOf <- function(t) do.call(paste, c(lapply(t[keys], as.character), sep = "\r"))
    ka <- keyOf(a)
    kb <- keyOf(b)
    header <- identical(names(a), names(b))
    order <- identical(ka, kb)
    ia <- match(intersect(ka, kb), ka)
    ib <- match(intersect(ka, kb), kb)
    callsOf <- function(t) {
        q <- num(t[[q_col]])
        l <- num(t[[fc_col]])
        hit <- !is.na(q) & q < DA_Q_THRESHOLD & !is.na(l)
        ifelse(hit & l >= 0, "up", ifelse(hit, "down", "none"))
    }
    calls_changed <- sum(callsOf(a)[ia] != callsOf(b)[ib])
    cols <- setdiff(intersect(names(a), names(b)), keys)
    numeric_col <- vapply(cols, function(cl) {
        xa <- a[[cl]]
        xb <- b[[cl]]
        all(isMissing(xa) | !is.na(num(xa))) && all(isMissing(xb) | !is.na(num(xb)))
    }, logical(1))
    disc <- sum(vapply(cols[!numeric_col], function(cl) {
        sum(as.character(a[[cl]][ia]) != as.character(b[[cl]][ib]))
    }, numeric(1)))
    missing_same <- TRUE
    rows <- lapply(cols[numeric_col], function(cl) {
        x <- num(a[[cl]][ia])
        y <- num(b[[cl]][ib])
        if (!identical(is.na(x), is.na(y))) missing_same <<- FALSE
        both <- !is.na(x) & !is.na(y)
        d <- abs(x[both] - y[both])
        m <- pmax(abs(x[both]), abs(y[both]))
        rel <- ifelse(m > 0, d / m, 0)
        bad <- if (mode == "exact") d != 0 else !(d <= SWD083_REL * m | d <= SWD083_ABS)
        data.frame(column = cl, compared = sum(both), differ = sum(d != 0),
                   max_abs = if (length(d)) max(d) else 0,
                   max_rel = if (length(rel)) max(rel) else 0,
                   outside = sum(bad) + sum(xor(is.na(x), is.na(y))), stringsAsFactors = FALSE)
    })
    columns <- do.call(rbind, rows) %||%
        data.frame(column = character(), compared = integer(), differ = integer(),
                   max_abs = numeric(), max_rel = numeric(), outside = integer())
    exact_parts <- header && order && nrow(a) == nrow(b) && length(ia) == nrow(a) &&
        calls_changed == 0L && missing_same && disc == 0
    verdict <- if (!exact_parts || any(columns$outside > 0)) "outside" else
        if (all(columns$differ == 0)) "identical" else "agree"
    list(verdict = verdict, mode = mode, header_identical = header, key_order_identical = order,
         keys_only_one_side = length(setdiff(union(ka, kb), intersect(ka, kb))),
         missing_values_same = missing_same, discrete_mismatches = as.integer(disc),
         calls_changed = as.integer(calls_changed), columns = columns)
}

#' Read the environment record of a control job (`environment.tsv`: key, value)
#' @noRd
valReadEnvironment <- function(path) {
    t <- valReadTextTable(path)
    stats::setNames(t$value, t$key)
}

#' The DA reference a re-freezing link supplies (KIT-199: `da-runner.tsv`, its two control jobs'
#' `environment.tsv` and `da-runner.tsv`)
#'
#' @param table The link's `da-runner.tsv`.
#' @param envs The control jobs' `environment.tsv` files.
#' @param controls The control jobs' own `da-runner.tsv` files (raw agreement, agreement 4).
#' @return A list: `table`, `cpus`, `env` (the D7 keys where the control jobs agree, else `NA`),
#'   `raw_agreed`.
#' @export
valReadDaLink <- function(table, envs, controls) {
    e <- lapply(envs, valReadEnvironment)
    env <- vapply(G1A_D7_ENV_KEYS, function(k) {
        v <- unique(vapply(e, function(x) unname(x[k] %||% NA_character_), character(1)))
        if (length(v) == 1L) v else NA_character_
    }, character(1))
    sums <- vapply(controls, valSha256, character(1))
    list(table = valReadTextTable(table),
         cpus = unique(vapply(e, function(x) unname(x[["cpu_model"]]), character(1))),
         env = env, raw_agreed = length(unique(sums)) == 1L)
}

#' D7's mode for the G1(a) job (KIT-199 12.3 as corrected by its implementation audit I04)
#'
#' `exact` only when the job's CPU model is one of the link's control jobs', those jobs agreed byte
#' for byte, every key of `G1A_D7_ENV_KEYS` equals the link's, and the lockfile is unchanged since
#' the link's commit; otherwise `swd083`.
#'
#' @param job Named values of the G1(a) job (`cpu_model` and the D7 keys).
#' @param link From [valReadDaLink()].
#' @param lock_unchanged `packaging/multischolar.lock.json` unchanged since the link's commit.
#' @return A list: `mode`, `reasons`.
#' @export
valG1aDaMode <- function(job, link, lock_unchanged) {
    job <- unlist(job)
    why <- character()
    if (!isTRUE(link$raw_agreed)) {
        why <- c(why, "the link's control jobs did not agree byte for byte")
    }
    if (!(job[["cpu_model"]] %||% "") %in% link$cpus) why <- c(why, "another CPU model")
    for (k in G1A_D7_ENV_KEYS) {
        if (is.na(link$env[[k]]) || !identical(unname(job[k]), unname(link$env[[k]]))) {
            why <- c(why, paste("environment", k))
        }
    }
    if (!isTRUE(lock_unchanged)) why <- c(why, "the lockfile changed")
    list(mode = if (length(why)) "swd083" else "exact", reasons = why)
}

#' The release-commit verdict of the DIA-with-limpa lane (G1(a) part (2), VAL-120 (2), VAL-133 (2))
#'
#' First the masked captures against the re-frozen file; a difference in any capture other than
#' the DA artifact is outside the classes (`not_certified`). Then, always, the live DA table's
#' values against the link's table in `mode`: keys, order, missing values and calls exactly,
#' numbers exactly or under SWD-083. Values that agree certify the lane; when the captures
#' differed, with the difference printed (`difference`).
#'
#' @param live_diffs The lane's differences from the re-frozen file (`field`, `key`, ...).
#' @param live_table The live DA table (text).
#' @param link From [valReadDaLink()] (its `table`).
#' @param mode From [valG1aDaMode()].
#' @return A list: `status`, `masked_identical`, `mode`, `reason`, `comparison`, `difference`.
#' @export
valG1aReleaseValues <- function(live_diffs, live_table, link, mode) {
    other <- live_diffs[!(live_diffs$field == "artifacts" & live_diffs$key == G1A_DA_KEY), ,
                        drop = FALSE]
    cmp <- valCompareDaValues(live_table, link$table, mode = mode)
    status <- "certified"
    reason <- ""
    if (nrow(other)) {
        status <- "not_certified"
        reason <- "captures other than the DA table differ from the re-frozen file"
    } else if (identical(cmp$verdict, "outside")) {
        status <- "not_certified"
        reason <- "the DA table's values disagree"
    }
    list(status = status, masked_identical = nrow(live_diffs) == 0L, mode = mode, reason = reason,
         comparison = cmp, difference = cmp$columns[cmp$columns$differ > 0, , drop = FALSE])
}

#' `compareLane`, extracted from the checkout's `packaging/verify/run-lanes.R` (never copied)
#'
#' @param pkg_root The package checkout.
#' @return The function `compareLane(lane, recorded, frozen)`.
#' @export
valLaneComparator <- function(pkg_root) {
    env <- valExtractHelpers(file.path(pkg_root, "packaging", "verify", "run-lanes.R"),
                             "compareLane", parent = globalenv())
    env$compareLane
}

#' A lane's live differences from an oracle, in the classifier's form
#' @noRd
valG1aLaneDiffs <- function(lane, live, frozen, compare_lane) {
    rows <- compare_lane(lane, live, frozen)
    if (is.null(rows) || nrow(rows) == 0L) {
        return(data.frame(lane = character(), step = character(), field = character(),
                          key = character(), base = character(), live = character(),
                          stringsAsFactors = FALSE))
    }
    na <- function(x) ifelse(x == "<absent>", NA_character_, x)
    data.frame(lane = rows$lane, step = rows$step, field = rows$field, key = rows$key,
               base = na(rows$oracle), live = na(rows$runner), stringsAsFactors = FALSE)
}

#' Part (1): one lane's captures, installed against `load_all`; on `prot_dia` the steps after `da`
#' are not compared (enrichment)
#'
#' @param lane The lane.
#' @param installed,load_all The two runs' steps.
#' @return A one-row data frame: `lane`, `status` (`identical`, `differs`), `first_step`.
#' @export
valG1aCompareModes <- function(lane, installed, load_all) {
    cut <- function(s) {
        if (!identical(lane, "prot_dia") || !"da" %in% names(s)) return(s)
        s[seq_len(match("da", names(s)))]
    }
    a <- cut(installed)
    b <- cut(load_all)
    first <- NA_character_
    for (st in union(names(b), names(a))) {
        if (!identical(a[[st]], b[[st]])) {
            first <- st
            break
        }
    }
    if (is.na(first) && !identical(names(a), names(b))) first <- "<step order>"
    data.frame(lane = lane, status = if (is.na(first)) "identical" else "differs",
               first_step = first, stringsAsFactors = FALSE)
}

#' The child process of a part (1) lane run (inside `testthat::test_that()`, as the controls of
#' KIT-180 R4b and KIT-199 ran; the DA copy hook is G1D's, `g1d.R:119-142`)
#' @noRd
G1A_CHILD <- r"---(
args <- commandArgs(TRUE)
pkg_root <- args[[1L]]; lane <- args[[2L]]; mode <- args[[3L]]; copy <- args[[4L]]
suppressPackageStartupMessages(library(testthat))
if (identical(mode, "load_all")) {
    suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE, export_all = TRUE, helpers = TRUE))
    env <- new.env(parent = globalenv())
    td <- file.path(pkg_root, "tests", "testthat")
    withr::with_dir(td, withr::with_envvar(c(TESTTHAT = "true"),
        suppressMessages(testthat::source_test_helpers(td, env = env))))
} else {
    suppressPackageStartupMessages(library(MultiScholaR))
    env <- msrval::valGuiHelpers(pkg_root)
}
run_lane <- get("guiEntryRunLane", envir = env)
fault <- NULL
if (nzchar(copy)) {
    fault <- list(step = "summary_refusal", fn = function() {
        caller <- parent.frame()
        step <- get("step", envir = caller)
        api <- get("api", envir = caller)
        base_dir <- get("project_dirs", envir = globalenv())[["proteomics"]][["base_dir"]]
        src <- file.path(base_dir, args[[5L]])
        if (!file.exists(src) || !file.copy(src, copy, overwrite = TRUE)) {
            stop("the DA table could not be copied after step da")
        }
        step$run(api)
    })
}
ok <- testthat::test_that(paste("G1(a) part 1", mode, lane), {
    run_lane(lane, fault = fault)
    testthat::succeed()
})
quit(save = "no", status = if (isTRUE(ok)) 0L else 1L)
)---"

#' Run one GUI-entry lane in a child process, on the installed package through the extracted
#' helpers (`installed`) or under `pkgload::load_all` (`load_all`, the oracles' load mode)
#'
#' @param pkg_root The package checkout (the child's working directory).
#' @param lane The lane.
#' @param mode `installed` or `load_all`.
#' @param out_dir Where the record is written (`<out_dir>/<mode>/<lane>.json`).
#' @param copy_da A path to copy the DIA-with-limpa DA table to after step `da`, or `NULL`.
#' @return The record's path; `val_g1a_lane_failed` when the child fails.
#' @export
valG1aRunLane <- function(pkg_root, lane, mode, out_dir, copy_da = NULL) {
    rec <- file.path(out_dir, mode)
    dir.create(rec, recursive = TRUE, showWarnings = FALSE)
    script <- tempfile(fileext = ".R")
    writeLines(G1A_CHILD, script)
    res <- processx::run(file.path(R.home("bin"), "Rscript"),
                         c("--vanilla", script, pkg_root, lane, mode, copy_da %||% "",
                           G1A_DA_KEY),
                         wd = pkg_root, error_on_status = FALSE, timeout = 3600,
                         env = c("current", MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = rec))
    path <- file.path(rec, paste0(lane, ".json"))
    writeLines(c(res$stdout, res$stderr), file.path(out_dir, sprintf("%s-%s.log", lane, mode)))
    if (res$status != 0L || !file.exists(path)) {
        valAbort("lane {lane} ({mode}) did not complete (exit {res$status})",
                 "val_g1a_lane_failed")
    }
    path
}

#' G1(a) part (1): every lane twice, installed against `load_all`, every compared step identical
#'
#' @param pkg_root The package checkout.
#' @param lanes The lanes.
#' @param out_dir The output folder (`g1a-part1.tsv`; records under `installed/`, `load_all/`).
#' @param run_lane The lane runner ([valG1aRunLane()]).
#' @return The table (`lane`, `status`: `identical`, `differs` or `error`, `first_step`); a
#'   difference stops G1(a) (the caller reports the first differing step, D3).
#' @export
valG1aPart1 <- function(pkg_root, lanes, out_dir, run_lane = valG1aRunLane) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    rows <- lapply(lanes, function(lane) {
        copy <- if (identical(lane, "prot_dia_limpa")) file.path(out_dir, "da-live.tsv") else NULL
        got <- tryCatch(list(
            installed = run_lane(pkg_root, lane, "installed", out_dir, copy_da = copy),
            load_all = run_lane(pkg_root, lane, "load_all", out_dir)
        ), error = function(e) e)
        if (inherits(got, "error")) {
            return(data.frame(lane = lane, status = "error", first_step = conditionMessage(got),
                              stringsAsFactors = FALSE))
        }
        steps <- lapply(got, function(p) jsonlite::read_json(p, simplifyVector = FALSE)$steps)
        valG1aCompareModes(lane, steps$installed, steps$load_all)
    })
    out <- do.call(rbind, rows)
    utils::write.table(out, file.path(out_dir, "g1a-part1.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "")
    out
}

#' The commits that changed a lane's oracle since the reference point, oldest first, on the first
#' parent line (merges of re-freezing tickets), each with its ticket's control record
#' @noRd
valG1aChainCommits <- function(repo, lane, at, since = G1A_SINCE) {
    path <- sprintf("tests/testdata/gui-entry/%s.json", lane)
    commits <- valGit(repo, "log", "--first-parent", "--reverse", "--format=%H",
                      paste0(since, "..", at), "--", path)$out
    commits <- commits[nzchar(commits)]
    lapply(commits, function(c) {
        changed <- valGit(repo, "diff", "--name-only", paste0(c, "^1"), c)$out
        rec <- grep("^docs/qa/workflow-execution/kit/kit-[0-9]+/gui-entry-control[.]tsv$", changed,
                    value = TRUE)
        list(commit = c, records = rec)
    })
}

#' G1(a) part (2): each lane against the oracle of the commit under test, through its re-freeze
#' chain (KIT-159 section 2)
#'
#' A lane whose oracle no re-freeze changed since `since`: its live differences classed
#' ([valClassifyLiveDiffs()]) and read against the reference counts ([valReadChain()]; `stop`
#' outside them, or `not_certified (changed, not re-frozen: KIT-<n>)` for a lane `not_refrozen`
#' names, C51). A re-frozen lane: every chain link built from its ticket's record at the merge
#' ([valReadControlLink()]; a merge without one is a link without a record) and read in order; the
#' live run is then expected identical to the re-frozen file, except on `prot_dia_limpa` when
#' `da` is given, where [valG1aReleaseValues()] decides (VAL-133 (2)). On `prot_dia` the steps
#' after `da` are not compared.
#'
#' @param repo The checkout (full history).
#' @param at The commit under test.
#' @param lanes The lanes.
#' @param live Named list: each lane's live steps.
#' @param reference The reference rows (G1(r)'s, with the adopted rows).
#' @param compare_lane From [valLaneComparator()].
#' @param since The reference point.
#' @param fingerprints Committed reference fingerprints, or `NULL`.
#' @param da For `prot_dia_limpa`: `list(table, link, mode)`, or `NULL`.
#' @param not_refrozen Named character: lane -> the ticket that changed it without a re-freeze.
#' @return A data frame: `lane`, `status`, `reason`, `chain`, `differences`, `classes`.
#' @export
valG1aPart2 <- function(repo, at, lanes, live, reference, compare_lane, since = G1A_SINCE,
                        fingerprints = NULL, da = NULL, not_refrozen = NULL) {
    rows <- lapply(lanes, function(lane) {
        text <- valGitShow(repo, at, sprintf("tests/testdata/gui-entry/%s.json", lane))
        if (is.null(text)) valAbort("no oracle of {lane} at {at}", "val_spec_invalid")
        frozen <- jsonlite::fromJSON(text, simplifyVector = FALSE)$steps
        lv <- live[[lane]]
        if (identical(lane, "prot_dia") && "da" %in% names(frozen)) {
            keep <- names(frozen)[seq_len(match("da", names(frozen)))]
            frozen <- frozen[keep]
            lv <- lv[intersect(names(lv), keep)]
        }
        diffs <- valG1aLaneDiffs(lane, lv, frozen, compare_lane)
        chain <- valG1aChainCommits(repo, lane, at, since)
        tickets <- character()
        if (length(chain) == 0L) {
            cl <- valClassifyLiveDiffs(diffs, stats::setNames(list(names(frozen)), lane),
                                       fingerprints = fingerprints)
            res <- valReadChain(lane, list(), reference, live_diffs = cl)
            status <- res$status
            if (status == "stop" && lane %in% names(not_refrozen)) {
                status <- sprintf("not_certified (changed, not re-frozen: %s)",
                                  not_refrozen[[lane]])
            }
            reason <- res$reason
            classes <- paste(sprintf("%s=%d", names(valStepClassCounts(cl)),
                                     valStepClassCounts(cl)), collapse = ";")
        } else {
            links <- lapply(chain, function(ln) {
                if (length(ln$records) != 1L) return(NULL)
                dir <- dirname(ln$records)
                ticket <- sub("^kit-", "KIT-", basename(dir))
                ctl <- valGitShow(repo, ln$commit, ln$records)
                if (!any(grepl(paste0("\t", lane, "\t"), valLinesOf(ctl), fixed = TRUE))) {
                    return(NULL)
                }
                tickets <<- c(tickets, ticket)
                keys <- valGitShow(repo, ln$commit, file.path(dir, "declared-keys.tsv"))
                valReadControlLink(ctl, declared = keys, ticket = ticket, lane = lane,
                                   lane_steps = names(frozen))
            })
            res <- valReadChain(lane, links, reference)
            status <- res$status
            reason <- res$reason
            classes <- ""
            if (identical(status, "certified")) {
                if (identical(lane, "prot_dia_limpa") && !is.null(da)) {
                    v <- valG1aReleaseValues(diffs, da$table, da$link, da$mode)
                    status <- v$status
                    reason <- v$reason
                } else if (nrow(diffs)) {
                    status <- "not_certified"
                    reason <- "the live run differs from the re-frozen file"
                }
            }
        }
        data.frame(lane = lane, status = status, reason = reason,
                   chain = paste(tickets, collapse = ","), differences = nrow(diffs),
                   classes = classes, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
