# run-tests: the shards' tables, failures and logs assembled into one results directory, the base-versus-head
# comparison when the request has a base commit, and the summary.
#
# Usage: Rscript --vanilla assemble.R <downloaded dir> <results dir> <requested files, one per line>
# Reads <downloaded dir>/shard-<side>-<k>/ (tests.tsv, failures/, logs/, results/*.rds, withheld.txt) for the sides
# head and base, check-<side>/ and scripts/ when present. Writes to <results dir>:
#   tests.tsv, failures.txt, logs/                   head (the commit under test), in the requested order; a requested
#                                                    file without a row (its shard job lost) is recorded as crashed
#   tests-base.tsv, failures-base.txt, logs-base/     the same for the base commit (a file the base tree lacks: absent)
#   compare.tsv, compare-failures.txt                 the comparison (below)
#   check-summary.txt, check/, check-base/, check-compare.txt, scripts/, withheld.txt, summary.md
#
# The comparison (the rule for runners: stored hashes recorded elsewhere cannot be judged by pass or fail; a change is
# clean when no test that passes on base fails on head and every expectation that fails on both fails with the same
# message). Per file, tests are matched by name (and occurrence); messages are compared after the job-specific parts
# (temporary paths, the checkout path) are masked. Verdict per file, the first that applies:
#   only-on-head      the base tree lacks the file (its head result is in the row)
#   new-failure       a test fails on head and passes on base, or is new on head and fails, or head did not complete
#                     (crashed, timeout) where base did; also a failing head file where base did not complete
#   message-changed   a test fails on both with a different set of messages
#   fixed             a test fails on base and passes on head (and nothing above)
#   count-changed     the counts differ (passed, failed, errors, skipped) with none of the above
#   same              identical counts and identical failure messages
# A failure identical on both sides whose message names a digest or hash but shows no hexadecimal value of 32 or more
# digits is noted (`note`): the identity of its actual value cannot be read from the message.
# Exit status 1 when a head file crashed or timed out or a script failed; 0 otherwise.

arguments <- commandArgs(trailingOnly = TRUE)
downloaded <- arguments[[1L]]; results <- arguments[[2L]]
requested <- readLines(arguments[[3L]], warn = FALSE)
requested <- requested[nzchar(requested)]
dir.create(results, recursive = TRUE, showWarnings = FALSE)
columns <- c("file", "passed", "failed", "errors", "skipped", "warnings", "seconds", "status")
COUNTS <- c("passed", "failed", "errors", "skipped")

shardsOf <- function(side) sort(list.files(downloaded, pattern = sprintf("^shard-%s-[0-9]+$", side), full.names = TRUE))
read <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, quote = "", na.strings = character(),
    colClasses = c(file = "character", status = "character"))

assembleSide <- function(side) {
    shards <- shardsOf(side)
    if (!length(shards)) return(NULL)
    suffix <- if (side == "head") "" else "-base"
    rows <- do.call(rbind, lapply(shards, function(shard) {
        path <- file.path(shard, "tests.tsv")
        if (file.exists(path)) read(path)[columns]
    }))
    if (is.null(rows)) rows <- stats::setNames(data.frame(matrix(nrow = 0L, ncol = 8L)), columns)
    lost <- setdiff(requested, rows$file)
    if (length(lost)) {
        rows <- rbind(rows, data.frame(file = lost, passed = 0L, failed = 0L, errors = 0L, skipped = 0L, warnings = 0L,
            seconds = NA, status = "crashed", stringsAsFactors = FALSE))
    }
    rows <- rows[match(requested, rows$file), columns]
    utils::write.table(rows, file.path(results, paste0("tests", suffix, ".tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "")
    failures <- character()
    for (file in requested) {
        for (shard in shards) {
            path <- file.path(shard, "failures", sub("[.]R$", ".txt", file))
            if (file.exists(path)) failures <- c(failures, readLines(path, warn = FALSE))
        }
        if (file %in% lost) failures <- c(failures, sprintf("file: %s", file), "test: <shard job>", "kind: crashed",
            "message:", "  no result: the shard job that held this file failed or was cancelled (see the run's job log)", "")
    }
    writeLines(if (length(failures)) failures else "no failing expectations", file.path(results, paste0("failures", suffix, ".txt")))
    logs_dir <- file.path(results, paste0("logs", suffix))
    dir.create(logs_dir, showWarnings = FALSE)
    detail <- list()
    withheld <- character()
    for (shard in shards) {
        file.copy(list.files(file.path(shard, "logs"), full.names = TRUE), logs_dir, overwrite = TRUE)
        for (path in list.files(file.path(shard, "results"), pattern = "[.]rds$", full.names = TRUE)) {
            detail[[paste0(sub("[.]rds$", "", basename(path)), ".R")]] <- readRDS(path)
        }
        if (file.exists(file.path(shard, "withheld.txt"))) {
            withheld <- c(withheld, paste0(basename(shard), "/", readLines(file.path(shard, "withheld.txt"))))
        }
    }
    list(rows = rows, detail = detail, withheld = withheld)
}

mask <- function(text) {
    text <- gsub("/tmp/Rtmp[A-Za-z0-9]+", "<tmp>", text)
    text <- gsub("/home/runner/work/[^/]+/[^/]+/pkg", "<pkg>", text)
    text <- gsub("/home/runner/work/_temp", "<runner-temp>", text)
    gsub("file[0-9a-f]{6,}", "file<id>", text)
}
keyed <- function(names) paste0(names, " #", stats::ave(seq_along(names), names, FUN = seq_along))
digestWithoutValue <- function(messages) {
    any(grepl("digest|sha256|sha1|md5|hash", messages, ignore.case = TRUE)) && !any(grepl("[0-9a-f]{32,}", messages))
}

compareFile <- function(b, h, bd, hd) {
    out <- list(verdict = "same", tests = character(), note = character(), pairs = list())
    add <- function(verdict, test, base_messages, head_messages, label = verdict) {
        out$tests <<- c(out$tests, sprintf("%s: %s", label, test))
        out$pairs[[length(out$pairs) + 1L]] <<- list(verdict = verdict, test = test, base = base_messages, head = head_messages)
    }
    incomplete <- c("crashed", "timeout")
    if (b$status == "absent") { out$verdict <- "only-on-head"; return(out) }
    if (h$status %in% incomplete || b$status %in% incomplete) {
        if (h$status == b$status) { out$note <- sprintf("%s on both", h$status); return(out) }
        if (h$status %in% incomplete) { out$verdict <- "new-failure"; out$note <- sprintf("head %s, base %s", h$status, b$status); return(out) }
        out$note <- sprintf("base %s: head's tests cannot be matched", b$status)
        out$verdict <- if (h$status == "ok") "fixed" else "new-failure"
        return(out)
    }
    counts_differ <- any(unlist(b[COUNTS]) != unlist(h[COUNTS]))
    if (is.null(bd) || is.null(hd)) {
        out$note <- "per-test results missing (withheld or not written): counts compared only"
        out$verdict <- if (counts_differ) "count-changed" else "same"
        return(out)
    }
    bt <- bd$tests; ht <- hd$tests
    bt$key <- keyed(bt$test); ht$key <- keyed(ht$test)
    messagesOf <- function(m, test) sort(mask(m$message[m$test == test]))
    fails <- function(t) t$failed + t$errors > 0L
    found <- character()
    for (i in seq_len(nrow(ht))) {
        test <- ht$test[[i]]
        j <- match(ht$key[[i]], bt$key)
        if (!fails(ht[i, ])) {
            if (!is.na(j) && fails(bt[j, ])) { found <- c(found, "fixed"); add("fixed", test, messagesOf(bd$messages, test), character()) }
            next
        }
        if (is.na(j)) {
            found <- c(found, "new-failure")
            add("new-failure", test, "<no such test on base>", messagesOf(hd$messages, test), "new-failure (test only on head)")
        } else if (!fails(bt[j, ])) {
            found <- c(found, "new-failure"); add("new-failure", test, character(), messagesOf(hd$messages, test))
        } else {
            bmsg <- messagesOf(bd$messages, test); hmsg <- messagesOf(hd$messages, test)
            if (!identical(bmsg, hmsg)) {
                found <- c(found, "message-changed"); add("message-changed", test, bmsg, hmsg)
            } else if (digestWithoutValue(hmsg)) {
                out$note <- c(out$note, sprintf("identical message without the actual digest value: %s", test))
            }
        }
    }
    for (j in seq_len(nrow(bt))) {
        if (fails(bt[j, ]) && is.na(match(bt$key[[j]], ht$key))) {
            found <- c(found, "fixed")
            add("fixed", bt$test[[j]], messagesOf(bd$messages, bt$test[[j]]), "<no such test on head>", "fixed (test gone on head)")
        }
    }
    out$verdict <- if ("new-failure" %in% found) "new-failure" else if ("message-changed" %in% found) "message-changed" else
        if ("fixed" %in% found) "fixed" else if (counts_differ) "count-changed" else "same"
    if (counts_differ) {
        moved <- merge(bt[c("key", "passed", "skipped")], ht[c("key", "passed", "skipped")], by = "key", all = TRUE,
            suffixes = c(".base", ".head"))
        moved <- moved[is.na(moved$passed.base) | is.na(moved$passed.head) | moved$passed.base != moved$passed.head |
            moved$skipped.base != moved$skipped.head, ]
        if (nrow(moved)) out$tests <- c(out$tests, sprintf("counts: %s (passed %s -> %s, skipped %s -> %s)", sub(" #1$", "", moved$key),
            moved$passed.base, moved$passed.head, moved$skipped.base, moved$skipped.head))
    }
    out
}

head_side <- assembleSide("head")
base_side <- assembleSide("base")
withheld <- c(head_side$withheld, base_side$withheld)
compare <- NULL
if (!is.null(base_side) && !is.null(head_side)) {
    rows <- list(); pairs <- character()
    for (i in seq_along(requested)) {
        file <- requested[[i]]
        b <- base_side$rows[i, ]; h <- head_side$rows[i, ]
        result <- compareFile(b, h, base_side$detail[[file]], head_side$detail[[file]])
        rows[[i]] <- data.frame(file = file, base_passed = b$passed, base_failed = b$failed, base_errors = b$errors,
            base_skipped = b$skipped, base_status = b$status, head_passed = h$passed, head_failed = h$failed,
            head_errors = h$errors, head_skipped = h$skipped, head_status = h$status, verdict = result$verdict,
            tests = gsub("[\t\r\n]+", " ", paste(result$tests, collapse = "; ")),
            note = gsub("[\t\r\n]+", " ", paste(result$note, collapse = "; ")), stringsAsFactors = FALSE)
        for (pair in result$pairs) {
            lines <- function(m, empty) if (length(m)) paste0("  ", unlist(strsplit(m, "\n", fixed = TRUE))) else empty
            pairs <- c(pairs, sprintf("== %s | %s | %s", file, pair$verdict, pair$test),
                "-- base:", lines(pair$base, "  <passes>"), "-- head:", lines(pair$head, "  <passes>"), "")
        }
    }
    compare <- do.call(rbind, rows)
    utils::write.table(compare, file.path(results, "compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(if (length(pairs)) pairs else "no non-identical failures", file.path(results, "compare-failures.txt"))
}

for (part in c("check-head", "check-base", "scripts")) {
    source_dir <- file.path(downloaded, part)
    if (!dir.exists(source_dir)) next
    to <- file.path(results, switch(part, "check-head" = "check", "check-base" = "check-base", scripts = "scripts"))
    dir.create(to, showWarnings = FALSE)
    file.copy(list.files(source_dir, full.names = TRUE), to, recursive = TRUE)
    if (file.exists(file.path(to, "withheld.txt"))) withheld <- c(withheld, paste0(part, "/", readLines(file.path(to, "withheld.txt"))))
    unlink(file.path(to, "withheld.txt"))
}
if (file.exists(file.path(results, "check", "check-summary.txt"))) {
    file.copy(file.path(results, "check", "check-summary.txt"), file.path(results, "check-summary.txt"), overwrite = TRUE)
}
check_compare <- NULL
findings <- function(dir) if (file.exists(file.path(dir, "findings.rds"))) readRDS(file.path(dir, "findings.rds"))
fh <- findings(file.path(results, "check")); fb <- findings(file.path(results, "check-base"))
if (!is.null(fh) && !is.null(fb)) {
    norm <- function(f) {
        if (!length(f$texts)) return(character())
        text <- gsub("\\[[0-9]+s(/[0-9]+s)?\\] ", "", f$texts)
        paste0("[", f$levels, "] ", mask(gsub("/home/runner/work/_temp/check-work", "<check>", text)))
    }
    nh <- norm(fh); nb <- norm(fb)
    count <- function(f, level) sum(f$levels == level)
    check_compare <- c(
        sprintf("base: ERRORs %d, WARNINGs %d, NOTEs %d (check exit %s)", count(fb, "ERROR"), count(fb, "WARNING"), count(fb, "NOTE"), fb$status),
        sprintf("head: ERRORs %d, WARNINGs %d, NOTEs %d (check exit %s)", count(fh, "ERROR"), count(fh, "WARNING"), count(fh, "NOTE"), fh$status),
        sprintf("findings identical on both: %d; only on head: %d; only on base: %d", length(intersect(nh, nb)),
            length(setdiff(nh, nb)), length(setdiff(nb, nh))), "",
        "== only on head (new)", if (length(setdiff(nh, nb))) c(setdiff(nh, nb), "") else "none", "",
        "== only on base (gone)", if (length(setdiff(nb, nh))) c(setdiff(nb, nh), "") else "none")
    writeLines(check_compare, file.path(results, "check-compare.txt"))
}
if (length(withheld)) writeLines(withheld, file.path(results, "withheld.txt"))

scripts <- if (file.exists(file.path(results, "scripts", "scripts.tsv"))) {
    utils::read.delim(file.path(results, "scripts", "scripts.tsv"), stringsAsFactors = FALSE, quote = "")
}
rows <- if (!is.null(head_side)) head_side$rows else stats::setNames(data.frame(matrix(nrow = 0L, ncol = 8L)), columns)
bad <- rows$status %in% c("crashed", "timeout")
scripts_bad <- !is.null(scripts) && any(scripts$status != "ok")

tableTitle <- function(rows, title) c(
    sprintf("### %s: %d files; passed %d, failed %d, errors %d, skipped %d; ok %d, failed %d, timeout %d, crashed %d%s", title,
        nrow(rows), sum(rows$passed), sum(rows$failed), sum(rows$errors), sum(rows$skipped), sum(rows$status == "ok"),
        sum(rows$status == "failed"), sum(rows$status == "timeout"), sum(rows$status == "crashed"),
        if (any(rows$status == "absent")) sprintf(", absent %d", sum(rows$status == "absent")) else ""), "")
md <- tableTitle(rows, "Tests (head)")
if (!is.null(compare)) {
    verdicts <- table(factor(compare$verdict, levels = c("same", "new-failure", "message-changed", "count-changed", "fixed", "only-on-head")))
    shown <- compare[compare$verdict != "same" | nzchar(compare$note), ]
    md <- c(md, tableTitle(base_side$rows, "Tests (base)"),
        sprintf("### Base versus head: %s", paste(sprintf("%s %d", names(verdicts), as.integer(verdicts)), collapse = ", ")),
        if (verdicts[["new-failure"]] + verdicts[["message-changed"]] == 0L) "**Clean under the rule** (no new failure, no changed failure message)" else
            "**Not clean under the rule**: see the rows below and compare-failures.txt", "",
        "| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |", "|---|---|---|---|---|---|",
        sprintf("| %s | %d/%d/%d/%d %s | %d/%d/%d/%d %s | %s | %s | %s |", shown$file, shown$base_passed, shown$base_failed,
            shown$base_errors, shown$base_skipped, shown$base_status, shown$head_passed, shown$head_failed, shown$head_errors,
            shown$head_skipped, shown$head_status, shown$verdict, substr(gsub("\\|", "/", shown$tests), 1L, 400L),
            gsub("\\|", "/", shown$note)), "")
} else {
    md <- c(md, "| file | passed | failed | errors | skipped | warnings | seconds | status |", "|---|---:|---:|---:|---:|---:|---:|---|",
        sprintf("| %s | %d | %d | %d | %d | %d | %s | %s |", rows$file, rows$passed, rows$failed, rows$errors, rows$skipped,
            rows$warnings, ifelse(is.na(rows$seconds), "", rows$seconds), rows$status), "")
}
if (!is.null(check_compare)) {
    md <- c(md, "### R CMD check, base versus head", "```", utils::head(check_compare, 80L), "```", "")
} else if (file.exists(file.path(results, "check-summary.txt"))) {
    md <- c(md, "### R CMD check", "```", utils::head(readLines(file.path(results, "check-summary.txt")), 60L), "```", "")
}
if (!is.null(scripts)) {
    md <- c(md, "### Scripts", "| # | path | exit | seconds | status |", "|---:|---|---:|---:|---|",
        sprintf("| %s | %s | %s | %s | %s |", scripts$index, scripts$path, scripts$exit, scripts$seconds, scripts$status), "")
}
if (length(withheld)) md <- c(md, sprintf("Withheld from the artifact (matched a forbidden pattern): %d files; see withheld.txt", length(withheld)), "")
md <- c(md, if (any(bad) || scripts_bad) "**Run outcome: failure** (a head file crashed or timed out, or a script failed)" else
    "**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)")
writeLines(md, file.path(results, "summary.md"))
cat(md, sep = "\n")
quit(save = "no", status = if (any(bad) || scripts_bad) 1L else 0L)
