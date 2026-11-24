# KIT-184: the control before the re-freeze of tests/testdata/gui-entry/prot_dia_limpa.json (owner's decision VAL-178
# (2); implementation audit AUDIT-KIT-184-impl-20261125, section 5). Modelled on kit-199/control.R (the same lane) and
# kit-180/gui-entry-control.R; no file of the package is edited and nothing is written into tests/.
#
# Usage (run-tests scripts job, from the package root of the merge commit, full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-184/gui-entry-control.R <staging commit> <out dir>
# Internal (one child process per lane run, its working directory at that side's tree root):
#   Rscript --vanilla gui-entry-control.R child <run|freeze> <steps file> <oracle dir or ->
#
# Three lane runs, one after another, each in its own process, each inside testthat::test_that() as the test file runs
# the lane (the helpers sourced as testthat::test_file() sources them):
#   base    a detached worktree of <staging commit> (with the workflow's native step for a package tree):
#           guiEntryRunLane("prot_dia_limpa"); its steps written to <out>/steps-base.json.
#   head-1  this checkout: guiEntryExpectLane("prot_dia_limpa") with MULTISCHOLAR_GUI_ENTRY_FREEZE=1 and
#           MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=<out>/oracle, so the frozen file is written to
#           <out>/oracle/prot_dia_limpa.json, never into the checkout; its steps to <out>/steps-head-1.json.
#   head-2  this checkout again: guiEntryRunLane("prot_dia_limpa"); its steps to <out>/steps-head-2.json.
# The stored reference is the staging commit's tests/testdata/gui-entry/prot_dia_limpa.json. Beside its steps it holds
# the commit, the date and R's version of its freeze, which differ from any later dump by construction; so "equal to
# the stored file" is checked on its `steps` member, written alone by the helper's own writer to
# <out>/steps-stored.json, byte for byte (SHA-256), as kit-199/control.R did.
#
# The three verdicts (<out>/verdict.tsv):
#   base_equals_stored         steps-base.json and steps-stored.json have the same SHA-256
#   head_dumps_identical       steps-head-1.json and steps-head-2.json have the same SHA-256
#   head_differs_in_four_only  between the stored steps and head-1's, the differing leaves are exactly the four added
#                              artifacts integration/proteomics/removed_proteomics.parquet = "<present>" of the steps
#                              summary_save, summary_copy, summary_report and summary_export
# and, as checks of the instrument: the frozen file's steps equal head-1's; the checkout's stored file is untouched.
# Also written: sha256.tsv, differences.tsv (every differing leaf, stored against head-1), status.tsv (one row per
# stage). Exit status non-zero only when a stage could not run; a verdict that does not hold is a result, in
# verdict.tsv.

LANE <- "prot_dia_limpa"
ADDED <- "integration/proteomics/removed_proteomics.parquet"
ADDED_STEPS <- c("summary_save", "summary_copy", "summary_report", "summary_export")
arguments <- commandArgs(trailingOnly = TRUE)

writeSteps <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE, null = "null", na = "string", digits = NA,
    pretty = TRUE)

if (length(arguments) == 4L && identical(arguments[[1L]], "child")) {
    mode <- arguments[[2L]]
    steps_path <- arguments[[3L]]
    oracle_dir <- arguments[[4L]]
    suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
    suppressPackageStartupMessages(library(testthat))
    # testthat::test_file() sources the helpers again into the test environment, with TESTTHAT set and the working
    # directory at the test directory (kit-199/control.R).
    helper_env <- new.env(parent = globalenv())
    test_dir <- normalizePath("tests/testthat")
    withr::with_dir(test_dir, withr::with_envvar(c(TESTTHAT = "true"),
        suppressMessages(testthat::source_test_helpers(test_dir, env = helper_env))))
    run_lane <- get("guiEntryRunLane", envir = helper_env, mode = "function")
    expect_lane <- get("guiEntryExpectLane", envir = helper_env, mode = "function")
    comparable <- get("guiEntryComparable", envir = helper_env, mode = "function")
    captured <- NULL
    problem <- ""
    if (identical(mode, "freeze")) {
        Sys.setenv(MULTISCHOLAR_GUI_ENTRY_FREEZE = "1", MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR = oracle_dir)
        # guiEntryExpectLane resolves guiEntryRunLane in helper_env: the record is kept on its way through.
        assign("guiEntryRunLane", function(lane_id, fault = NULL) {
            record <- run_lane(lane_id, fault = fault)
            captured <<- record
            record
        }, envir = helper_env)
        passed <- tryCatch(testthat::test_that(paste("KIT-184 control, freeze", LANE), expect_lane(LANE)),
            error = function(e) {
                problem <<- conditionMessage(e)
                FALSE
            })
    } else {
        passed <- tryCatch(testthat::test_that(paste("KIT-184 control, run", LANE), {
            captured <<- run_lane(LANE)
            testthat::succeed()
        }), error = function(e) {
            problem <<- conditionMessage(e)
            FALSE
        })
    }
    if (!is.null(captured)) writeSteps(comparable(captured), steps_path)
    ok <- !is.null(captured) && isTRUE(passed)
    cat("KIT-184 child", mode, if (ok) "ok" else "FAILED", substr(gsub("[\t\r\n]+", " ", problem), 1L, 200L), "\n")
    quit(status = if (ok) 0L else 1L)
}

stopifnot(length(arguments) == 2L, grepl("^[0-9a-f]{40}$", arguments[[1L]]))
base_commit <- arguments[[1L]]
out <- arguments[[2L]]
dir.create(file.path(out, "oracle"), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
root <- normalizePath(".")
script <- file.path(root, "docs", "qa", "workflow-execution", "kit", "kit-184", "gui-entry-control.R")
stopifnot(file.exists(script))
rscript <- file.path(R.home("bin"), "Rscript")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
work <- file.path(tempdir(), "k184-control")
dir.create(work)
cat("KIT-184 control: head", head_commit, "base", base_commit, "\n")

status <- list()
note <- function(stage, outcome, detail = "") {
    status[[length(status) + 1L]] <<- data.frame(stage = stage, outcome = outcome, detail = detail, stringsAsFactors = FALSE)
    cat(sprintf("%s: %s %s\n", stage, outcome, detail))
}
finish <- function(code) {
    utils::write.table(do.call(rbind, status), file.path(out, "status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    quit(status = code)
}
sha256 <- function(path) if (file.exists(path)) digest::digest(file = path, algo = "sha256") else NA_character_

# The base tree, with the workflow's native step for a package tree (run-tests.yml, "Install package_commit").
base_tree <- file.path(work, "base")
made <- system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = FALSE, stderr = FALSE)
if (!identical(as.integer(made), 0L) || !dir.exists(file.path(base_tree, "tests", "testthat"))) {
    note("worktree", "could not run", "git worktree add failed")
    finish(1L)
}
note("worktree", "ok", base_commit)
manifest <- file.path(base_tree, "inst", "duckdb", "linux_amd64", "manifest.rds")
if (file.exists(file.path(base_tree, "src", "build_dia_rsum.R")) && !file.exists(manifest)) {
    withr::with_dir(base_tree, suppressWarnings(system2(rscript, c("--vanilla", "src/build_dia_rsum.R"), stdout = FALSE,
        stderr = FALSE)))
    if (!file.exists(manifest)) {
        note("native", "could not run", "the native step did not complete on the base tree")
        finish(1L)
    }
    note("native", "ok", "")
} else {
    note("native", "ok", "not needed")
}

# One lane run in its own process, its working directory at `tree`. Only the child's own line is printed.
runSide <- function(side, tree, mode, oracle_dir = "-") {
    steps_path <- file.path(out, paste0("steps-", side, ".json"))
    started <- Sys.time()
    console <- withr::with_dir(tree, suppressWarnings(system2(rscript, c("--vanilla", shQuote(script), "child", mode,
        shQuote(steps_path), shQuote(oracle_dir)), stdout = TRUE, stderr = TRUE)))
    code <- attr(console, "status")
    code <- if (is.null(code)) 0L else as.integer(code)
    line <- grep("^KIT-184 child", console, value = TRUE)
    seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")))
    ok <- identical(code, 0L) && file.exists(steps_path)
    note(paste("lane", side), if (ok) "ok" else "could not run", sprintf("%s (%d s)", if (length(line)) line[[length(line)]] else
        "the child printed no result line", seconds))
    ok
}
stored_path <- file.path("tests", "testdata", "gui-entry", paste0(LANE, ".json"))
stored_before <- sha256(file.path(root, stored_path))
if (!runSide("base", base_tree, "run")) finish(1L)
if (!runSide("head-1", root, "freeze", file.path(out, "oracle"))) finish(1L)
if (!runSide("head-2", root, "run")) finish(1L)
frozen_path <- file.path(out, "oracle", paste0(LANE, ".json"))
if (!file.exists(frozen_path)) {
    note("freeze", "could not run", "the frozen file was not written")
    finish(1L)
}

# The stored reference of the staging commit, and the frozen file: their `steps` members written alone.
readSteps <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)$steps
stored <- readSteps(file.path(base_tree, stored_path))
writeSteps(stored, file.path(out, "steps-stored.json"))
frozen <- jsonlite::read_json(frozen_path, simplifyVector = FALSE)
writeSteps(frozen$steps, file.path(out, "steps-frozen.json"))
head1 <- jsonlite::read_json(file.path(out, "steps-head-1.json"), simplifyVector = FALSE)

# Every leaf of a steps list as a named character vector, path -> value (the value as JSON text).
asText <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "string", digits = NA))
leaves <- function(x, prefix = character()) {
    if (is.list(x) && length(x)) {
        keys <- names(x)
        if (is.null(keys)) keys <- as.character(seq_along(x))
        return(unlist(lapply(seq_along(x), function(i) leaves(x[[i]], c(prefix, keys[[i]]))), use.names = TRUE))
    }
    stats::setNames(asText(x), paste(prefix, collapse = "\t"))
}
a <- leaves(stored)
b <- leaves(head1)
stopifnot(!anyDuplicated(names(a)), !anyDuplicated(names(b)))
paths <- union(names(a), names(b))
in_stored <- unname(a[match(paths, names(a))])
in_head <- unname(b[match(paths, names(b))])
differs <- is.na(in_stored) != is.na(in_head) | (!is.na(in_stored) & !is.na(in_head) & in_stored != in_head)
shown <- function(value) ifelse(is.na(value), "<absent>", substr(value, 1L, 160L))
differences <- data.frame(path = gsub("\t", " > ", paths[differs], fixed = TRUE), stored = shown(in_stored[differs]),
    head = shown(in_head[differs]), stringsAsFactors = FALSE, row.names = NULL)
utils::write.table(differences, file.path(out, "differences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
expected <- paste(ADDED_STEPS, "artifacts", ADDED, sep = "\t")
at <- match(expected, paths)
four_only <- setequal(paths[differs], expected) && !anyNA(at) && all(is.na(in_stored[at])) &&
    all(in_head[at] %in% asText("<present>"))
differs <- paths[differs]

files <- c("steps-stored.json", "steps-base.json", "steps-head-1.json", "steps-head-2.json", "steps-frozen.json")
hashes <- vapply(file.path(out, files), sha256, character(1))
names(hashes) <- files
utils::write.table(data.frame(file = c(files, paste0("oracle/", LANE, ".json"), paste0(stored_path, " (checkout, before)"),
    paste0(stored_path, " (checkout, after)"), paste0(stored_path, " (staging commit)")),
    sha256 = c(unname(hashes), sha256(frozen_path), stored_before, sha256(file.path(root, stored_path)),
        sha256(file.path(base_tree, stored_path)))), file.path(out, "sha256.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
verdict <- data.frame(check = c("base_equals_stored", "head_dumps_identical", "head_differs_in_four_only",
    "frozen_file_steps_equal_head", "frozen_file_names_this_commit", "checkout_reference_untouched", "differing_leaves"),
    value = c(as.character(c(identical(hashes[["steps-base.json"]], hashes[["steps-stored.json"]]),
        identical(hashes[["steps-head-1.json"]], hashes[["steps-head-2.json"]]), four_only,
        identical(hashes[["steps-frozen.json"]], hashes[["steps-head-1.json"]]),
        identical(as.character(frozen$source_sha), head_commit[[1L]]),
        identical(stored_before, sha256(file.path(root, stored_path))))), as.character(length(differs))),
    stringsAsFactors = FALSE)
utils::write.table(verdict, file.path(out, "verdict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(verdict, row.names = FALSE)
note("comparison", "ok", sprintf("%d differing leaves", length(differs)))
invisible(system2("git", c("worktree", "remove", "--force", shQuote(base_tree)), stdout = FALSE, stderr = FALSE))
finish(0L)
