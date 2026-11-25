# KIT-270: the three "no Windows path-length hint" guards of tests/testthat/test-workflow-figure-paths.R, read from git
# at a revision and evaluated on their own, in three modes.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-270/guards-check.R <revisions> <modes> <repeats>
#   <revisions>  comma-separated git revisions (full hashes)
#   <modes>      comma-separated: repeat | seeded | mutant
#   <repeats>    iterations for `repeat` (ignored by the other modes)
# Run from the package root after the package is installed, with full git history.
#   repeat  both test_that blocks that hold the guards, <repeats> times in this one process, each under a fresh
#           ListReporter: counts the iterations whose guard failed on a message with no hint (a temporary folder holding
#           259) and the iterations that failed in any other way.
#   seeded  the blocks once, every withr::local_tempdir() made with pattern "file259x" (a folder holding 259).
#   mutant  the blocks once, the guards' messages made to carry the hint (windows = TRUE; `long` for `empty`).
# Output: kit270-out/<mode>-<revision, 8 characters>.tsv and the same on the console.

file_path <- "tests/testthat/test-workflow-figure-paths.R"
test_names <- c(
    "the figure check names the missing files and, on Windows, the over-long path",
    "a failed copy says, on Windows, that the project's path is too long with long paths disabled"
)
helper_name <- "figurePathText"

state <- new.env()

# Replace each call in `x` by f(call), then walk into the result.
walkCalls <- function(x, f) {
    if (!is.call(x)) return(x)
    x <- f(x)
    if (!is.call(x)) return(x)
    for (i in seq_along(x)) {
        if (is.name(x[[i]]) && !nzchar(as.character(x[[i]]))) next
        element <- x[[i]]
        if (is.call(element)) x[[i]] <- walkCalls(element, f)
    }
    x
}

isGuard <- function(call) {
    is.call(call) && length(call) >= 3L && identical(call[[1L]], quote(expect_no_match)) &&
        is.character(call[[3L]]) && startsWith(call[[3L]], "259")
}

withMutation <- function(call, id) {
    force_windows <- function(node) {
        walkCalls(node, function(inner) {
            nms <- names(inner)
            if (!is.null(nms)) for (i in which(nms == "windows")) if (identical(inner[[i]], FALSE)) inner[[i]] <- TRUE
            inner
        })
    }
    call[[2L]] <- force_windows(call[[2L]])
    if (id == 2L) call[[2L]] <- do.call(substitute, list(call[[2L]], list(empty = quote(long))))
    call
}

isLocalTempdir <- function(call) {
    identical(call[[1L]], quote(withr::local_tempdir)) && length(call) == 1L
}

# The selected test calls and the helper, from the file at `revision`, with the guards numbered in file order.
selectCode <- function(revision, mode) {
    text <- system2("git", c("show", paste0(revision, ":", file_path)), stdout = TRUE)
    exprs <- parse(text = text, keep.source = FALSE)
    counter <- 0L
    patterns <- character()
    rewrite <- function(call) {
        if (isGuard(call)) {
            counter <<- counter + 1L
            patterns[[counter]] <<- call[[3L]]
            if (mode == "mutant") call <- withMutation(call, counter)
            call[[1L]] <- quote(.guard)
            call <- as.call(append(as.list(call), counter, after = 1L))
        } else if (mode == "seeded" && isLocalTempdir(call)) {
            call <- quote(withr::local_tempdir(pattern = "file259x"))
        }
        call
    }
    tests <- list()
    helper <- NULL
    for (expr in exprs) {
        if (!is.call(expr)) next
        if (identical(expr[[1L]], quote(test_that)) && is.character(expr[[2L]]) && expr[[2L]] %in% test_names) {
            tests[[length(tests) + 1L]] <- walkCalls(expr, rewrite)
        } else if (identical(expr[[1L]], quote(`<-`)) && identical(expr[[2L]], as.name(helper_name))) {
            helper <- expr
        }
    }
    if (length(tests) != length(test_names) || is.null(helper) || counter != 3L)
        stop("expected 2 test blocks, the helper and 3 guards at ", revision, "; found ", length(tests), ", ",
            !is.null(helper), ", ", counter)
    list(tests = tests, helper = helper, patterns = patterns)
}

# The guard: evaluates the message, records the expectation's outcome and the raw text, passes the expectation on.
.guard <- function(id, object, regexp, ...) {
    text <- object
    outcome <- NA_character_
    on.exit(state$log[[length(state$log) + 1L]] <- list(id = id, outcome = outcome, text = text, regexp = regexp),
        add = TRUE)
    withCallingHandlers(
        testthat::expect_no_match(text, regexp, ...),
        expectation = function(e) {
            outcome <<- if (inherits(e, "expectation_failure")) "fail" else if (inherits(e, "expectation_success")) "pass" else class(e)[[1L]]
        }
    )
    invisible(text)
}

# One pass over the selected blocks under a fresh ListReporter.
runOnce <- function(code) {
    state$log <- list()
    env <- new.env(parent = asNamespace("MultiScholaR"))
    env$.guard <- .guard
    eval(code$helper, env)
    reporter <- testthat::ListReporter$new()
    testthat::with_reporter(reporter, for (call in code$tests) eval(call, env))
    kinds <- unlist(lapply(reporter$get_results(), function(test) vapply(test$results, function(r) class(r)[[1L]], "")))
    list(log = state$log, broken = sum(kinds %in% c("expectation_failure", "expectation_error")),
        errors = sum(kinds == "expectation_error"), passed = sum(kinds == "expectation_success"))
}

guardRows <- function(revision, mode, code, result) {
    do.call(rbind, lapply(1:3, function(id) {
        entry <- Filter(function(e) e$id == id, result$log)
        entry <- if (length(entry)) entry[[1L]] else list(outcome = "not reached", text = "", regexp = NA_character_)
        data.frame(revision = revision, mode = mode, guard = id, pattern = code$patterns[[id]], outcome = entry$outcome,
            hint_in_message = grepl("259 characters", entry$text, fixed = TRUE),
            folder_holds_259 = grepl("/file[0-9a-fx]*259", entry$text),
            stringsAsFactors = FALSE)
    }))
}

main <- function() {
    library(testthat, warn.conflicts = FALSE)
    testthat::local_edition(3)
    arguments <- commandArgs(trailingOnly = TRUE)
    if (length(arguments) != 3L) stop("usage: Rscript guards-check.R <revisions> <modes> <repeats>")
    revisions <- strsplit(arguments[[1L]], ",", fixed = TRUE)[[1L]]
    modes <- strsplit(arguments[[2L]], ",", fixed = TRUE)[[1L]]
    repeats <- as.integer(arguments[[3L]])
    suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
    dir.create("kit270-out", showWarnings = FALSE)
    cat(sprintf("pid %d, hex %x, session folder %s\n", Sys.getpid(), Sys.getpid(), basename(tempdir())))
    for (revision in revisions) for (mode in modes) {
        code <- selectCode(revision, mode)
        started <- Sys.time()
        if (mode == "repeat") {
            f1 <- 0L
            f2 <- 0L
            guard_failures <- c(0L, 0L, 0L)
            for (i in seq_len(repeats)) {
                result <- runOnce(code)
                failed_guards <- vapply(result$log, function(e) identical(e$outcome, "fail"), NA)
                no_hint <- vapply(result$log, function(e) !grepl("259 characters", e$text, fixed = TRUE), NA)
                for (e in result$log[failed_guards]) guard_failures[[e$id]] <- guard_failures[[e$id]] + 1L
                folder_failure <- any(failed_guards & no_hint)
                other <- result$broken - sum(failed_guards) > 0L || any(failed_guards & !no_hint) || length(result$log) != 3L
                f1 <- f1 + as.integer(folder_failure)
                f2 <- f2 + as.integer(other)
            }
            table <- data.frame(revision = revision, mode = mode, pid = Sys.getpid(), pid_hex = sprintf("%x", Sys.getpid()),
                session_folder = basename(tempdir()), repeats = repeats, f1_iterations = f1, f2_iterations = f2,
                guard1_failures = guard_failures[[1L]], guard2_failures = guard_failures[[2L]],
                guard3_failures = guard_failures[[3L]], seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs"))),
                stringsAsFactors = FALSE)
        } else {
            result <- runOnce(code)
            table <- guardRows(revision, mode, code, result)
            table$other_broken <- result$broken - sum(table$outcome == "fail")
            table$other_errors <- result$errors
            table$expectations_passed <- result$passed
        }
        out <- file.path("kit270-out", sprintf("%s-%s.tsv", mode, substr(revision, 1L, 8L)))
        utils::write.table(table, out, sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("== %s %s (%s)\n", mode, substr(revision, 1L, 8L), out))
        print(table, row.names = FALSE)
    }
}

main()
