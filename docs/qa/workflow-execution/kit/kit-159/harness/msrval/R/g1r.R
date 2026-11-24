# G1(r)'s own reference rows (KIT-159 section 2, G1(r); re-audit 2, R2.6 item 1): a reader of the
# committed source blocks (`g1a-reference-source.txt`, the 69 GUI-entry blocks of run
# 44953187101's failures.txt) into valClassifyLiveDiffs()'s input, and the rows in the column form
# of KIT-172's `k156-gui-entry-failures.tsv`.

G1R_INFO <- "^[A-Za-z0-9_]+ at step '[^']+' \\(oracle [0-9a-f]+\\)$"

#' Read the GUI-entry failure blocks of a testthat failures file (G1(r))
#'
#' A block runs from a line `file: <name>` to the next. Its message is the lines after `message:`
#' up to the line starting `trace`, without their two-space indent. Its lane, step and oracle come
#' from the info line `<lane> at step '<step>' (oracle <sha>)` (`helper-gui-entry.R:414`).
#'
#' @param x A path, or the file's lines.
#' @return A list of blocks, each a list: `file`, `lane`, `step`, `oracle`, `message`.
#' @export
valReadG1rBlocks <- function(x) {
    lines <- valLinesOf(x)
    starts <- which(startsWith(lines, "file: "))
    if (length(starts) == 0L) valAbort("the source holds no block", "val_g1r_unreadable")
    ends <- c(starts[-1L] - 1L, length(lines))
    lapply(seq_along(starts), function(i) {
        b <- lines[starts[[i]]:ends[[i]]]
        m_at <- which(b == "message:")
        if (length(m_at) != 1L) valAbort("a block without one message", "val_g1r_unreadable")
        t_at <- which(startsWith(b, "trace"))
        t_at <- t_at[t_at > m_at]
        stop_at <- if (length(t_at)) t_at[[1L]] - 1L else length(b)
        msg <- if (stop_at > m_at) sub("^  ", "", b[(m_at + 1L):stop_at]) else character()
        info <- msg[grepl(G1R_INFO, msg)]
        if (length(info) != 1L) {
            valAbort("a block names no lane and step", "val_g1r_unreadable")
        }
        list(file = sub("^file: ", "", b[[1L]]),
             lane = sub(" at step .*$", "", info),
             step = sub("^.* at step '([^']+)'.*$", "\\1", info),
             oracle = sub("^.*\\(oracle ([0-9a-f]+)\\)$", "\\1", info),
             message = msg)
    })
}

#' One side of a printed comparison: `actual$<field>$<key>` (or `[[i]]...` after the field)
#' @noRd
valG1rRef <- function(ref) {
    m <- regmatches(ref, regexec("^(actual|expected)\\$([A-Za-z_]+)(.*)$", ref))[[1L]]
    if (length(m) != 4L) valAbort("an unreadable reference in a comparison", "val_g1r_unreadable")
    key <- sub("^\\$", "", m[[4L]])
    field <- m[[3L]]
    if (identical(field, "artifacts") && identical(key, "<store entries>")) field <- "store_entries"
    list(side = m[[2L]], field = field, key = key)
}

#' A printed value without its enclosing double quotes
#' @noRd
valG1rUnquote <- function(v) {
    if (grepl("^\".*\"$", v) && nchar(v) >= 2L) substr(v, 2L, nchar(v) - 1L) else v
}

#' The differences of one block, in the four waldo forms (re-audit 2, R2.6)
#' @noRd
valG1rBlockDiffs <- function(block) {
    msg <- block$message
    where <- paste(block$lane, block$step)
    rows <- list()
    put <- function(field, key, side, value) {
        id <- paste(field, key, sep = "\r")
        r <- rows[[id]]
        if (is.null(r)) {
            r <- list(field = field, key = key, base = NA_character_, live = NA_character_)
        }
        if (!is.na(value)) {
            col <- if (identical(side, "actual")) "live" else "base"
            r[[col]] <- if (is.na(r[[col]])) value else paste(r[[col]], value, sep = "\n")
        }
        rows[[id]] <<- r
    }
    unreadable <- function() {
        valAbort("a line outside the four forms in {where}", "val_g1r_unreadable")
    }
    blank <- function(l) !nzchar(trimws(l))
    match1 <- function(re, l) regmatches(l, regexec(re, l, perl = TRUE))[[1L]]
    re_is <- "^`((actual|expected)\\$[^`]*)` is (.*)$"
    re_vs <- "^((actual|expected)\\$\\S+) vs ((actual|expected)\\$\\S+)$"
    re_names <- paste0("^names\\(actual\\$([A-Za-z_]+)\\)(\\[[0-9]+:[0-9]+\\])? vs ",
                       "names\\(expected\\$([A-Za-z_]+)\\)(\\[[0-9]+:[0-9]+\\])?$")
    re_lines <- paste0("^lines\\(actual\\$([A-Za-z_]+)\\$(.+?)\\)(\\[[0-9]+:[0-9]+\\])? vs ",
                       "lines\\(expected\\$([A-Za-z_]+)\\$(.+?)\\)(\\[[0-9]+:[0-9]+\\])?$")
    # the value lines of a `vs` or `lines` comparison, up to the next blank line
    values <- function(i, field, key, context) {
        while (i <= length(msg) && !blank(msg[[i]])) {
            v <- msg[[i]]
            if (startsWith(v, "- ")) {
                put(field, key, "actual", valG1rUnquote(substring(v, 3L)))
            } else if (startsWith(v, "+ ")) {
                put(field, key, "expected", valG1rUnquote(substring(v, 3L)))
            } else if (!(context && startsWith(v, "  "))) {
                unreadable()
            }
            i <- i + 1L
        }
        i
    }
    i <- 1L
    while (i <= length(msg)) {
        l <- msg[[i]]
        if (blank(l) || identical(l, "Differences:") ||
                startsWith(l, "Expected `steps[[id]]`") || grepl(G1R_INFO, l)) {
            i <- i + 1L
            next
        }
        m <- match1(re_is, l)
        if (length(m)) {
            r <- valG1rRef(m[[2L]])
            desc <- m[[4L]]
            if (!nzchar(r$key)) {
                if (!grepl("^length [0-9]+$", desc)) unreadable()
            } else {
                cv <- match1("^a character vector \\('(.*)'\\)$", desc)
                value <- if (identical(desc, "absent")) {
                    NA_character_
                } else if (length(cv)) {
                    cv[[2L]]
                } else {
                    paste0("<", desc, ">")
                }
                put(r$field, r$key, r$side, value)
            }
            i <- i + 1L
            next
        }
        m <- match1(re_vs, l)
        if (length(m)) {
            a <- valG1rRef(m[[2L]])
            e <- valG1rRef(m[[4L]])
            if (!identical(a$side, "actual") || !identical(e$side, "expected") ||
                    !identical(a$field, e$field) || !identical(a$key, e$key) || !nzchar(a$key)) {
                unreadable()
            }
            put(a$field, a$key, "actual", NA_character_)
            i <- values(i + 1L, a$field, a$key, context = FALSE)
            next
        }
        m <- match1(re_names, l)
        if (length(m)) {
            if (!identical(m[[2L]], m[[4L]])) unreadable()
            i <- i + 1L
            while (i <= length(msg) && !blank(msg[[i]])) {
                q <- match1("^([-+ ]) \"(.*)\"$", msg[[i]])
                if (!length(q)) unreadable()
                if (!identical(q[[2L]], " ")) {
                    f <- m[[2L]]
                    if (identical(f, "artifacts") && identical(q[[3L]], "<store entries>")) {
                        f <- "store_entries"
                    }
                    put(f, q[[3L]], if (identical(q[[2L]], "-")) "actual" else "expected",
                        NA_character_)
                }
                i <- i + 1L
            }
            next
        }
        m <- match1(re_lines, l)
        if (length(m)) {
            if (!identical(m[[2L]], m[[5L]]) || !identical(m[[3L]], m[[6L]])) unreadable()
            put(m[[2L]], m[[3L]], "actual", NA_character_)
            i <- values(i + 1L, m[[2L]], m[[3L]], context = TRUE)
            next
        }
        unreadable()
    }
    if (length(rows) == 0L) valAbort("a block gives no difference: {where}", "val_g1r_unreadable")
    col <- function(name) unname(vapply(rows, `[[`, "", name))
    data.frame(lane = block$lane, step = block$step, field = col("field"), key = col("key"),
               base = col("base"), live = col("live"), oracle = block$oracle,
               stringsAsFactors = FALSE)
}

#' The differences of G1(r)'s source blocks, in valClassifyLiveDiffs()'s input form
#'
#' One row per lane, step, field and key: `base` is the frozen oracle's value (waldo's
#' `expected`), `live` the run's (`actual`), `NA` where the key is on one side only; a `lines`
#' comparison's differing lines joined by line breaks. The four forms (re-audit 2, R2.6):
#' `` `actual$<field>$<key>` is ... `` / `is absent`; `actual$<field>$<key> vs expected$...` with
#' `-` and `+` values; `names(actual$<field>)[a:b] vs ...` windows; and
#' `lines(actual$outputs$<key>)[a:b] vs ...`. An artifacts entry `<store entries>` is field
#' `store_entries`. A line in no form, or a block with no difference, is `val_g1r_unreadable`
#' (the gate's "fields not readable").
#'
#' @param blocks From [valReadG1rBlocks()].
#' @return A data frame: `lane`, `step`, `field`, `key`, `base`, `live`, `oracle`.
#' @export
valG1rDiffs <- function(blocks) {
    d <- do.call(rbind, lapply(blocks, valG1rBlockDiffs))
    rownames(d) <- NULL
    d
}

#' G1(r)'s reference rows from classed differences, in KIT-172's column form
#'
#' One row per lane, step and class, in source order: `fields` (comma-separated) and `keys`
#' (`field/key`, semicolon-separated), each sorted in C order; `names_list_truncated` is `no` (the
#' source lists every key).
#'
#' @param diffs [valG1rDiffs()] rows classed by [valClassifyLiveDiffs()].
#' @param source,source_sha256 The source's path text and SHA-256.
#' @return A data frame: `lane`, `step`, `class`, `fields`, `keys`, `names_list_truncated`,
#'   `oracle_sha`, `source`, `source_sha256`.
#' @export
valG1rRows <- function(diffs, source, source_sha256) {
    if (anyNA(diffs$class)) valAbort("a difference without a class of G1(a)", "val_g1r_no_class")
    g <- paste(diffs$lane, diffs$step, diffs$class, sep = "\r")
    rows <- lapply(unique(g), function(k) {
        d <- diffs[g == k, , drop = FALSE]
        data.frame(lane = d$lane[[1L]], step = d$step[[1L]], class = d$class[[1L]],
                   fields = paste(sort(unique(d$field), method = "radix"), collapse = ","),
                   keys = paste(sort(unique(paste0(d$field, "/", d$key)), method = "radix"),
                                collapse = ";"),
                   names_list_truncated = "no", oracle_sha = d$oracle[[1L]], source = source,
                   source_sha256 = source_sha256, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}
