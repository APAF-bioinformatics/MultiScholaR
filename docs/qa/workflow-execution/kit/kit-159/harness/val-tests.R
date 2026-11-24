# KIT-159 harness tests (section 7.2: `val-tests.R`). Run from the package root by a `run-tests`
# scripts request:
#
#   Rscript docs/qa/workflow-execution/kit/kit-159/harness/val-tests.R <out dir> [options]
#
#   (no option)           msrval's own tests -> <out>/harness-tests.tsv, harness-failures.txt
#   --files <names>       only these test files (comma-separated names after "test-", without ".R")
#   --covr                also covr line coverage of msrval -> <out>/coverage.txt
#   --lint                also lintr on msrval (lintr from `extra_packages`) -> <out>/lint.txt
#   --time-regime <id>    time one coverage regime in one process and in four concurrent chunks
#                         (f_sim, G14 (c)) -> <out>/fsim.tsv; with --nsim <n>; and a cost probe
#   --probe-cells <ids>   cost probe of the named cells only -> <out>/cost-probe.tsv
#                         (--probe-nsim, --probe-target-sims)
#   --coverage-plan <f>   run the plan's tasks of --request <id> in --chunks <k> concurrent
#                         processes -> coverage-detail.tsv, coverage-table.tsv, coverage-tasks.tsv
#   --ticket <dir>        run a ticket's pinned tests in <dir>/tests -> <out>/ticket-tests.tsv;
#                         with --declaration <json>, its coverage_regimes -> ticket-coverage.tsv
#
# msrval is installed from this checkout into $RUNNER_TEMP/val-lib first (7.8 step 1). Only fixed
# template lines are printed; the tests use synthetic data only.

harness_dir <- "docs/qa/workflow-execution/kit/kit-159/harness"
pkg_dir <- file.path(harness_dir, "msrval")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("usage: val-tests.R <out dir> [options]", call. = FALSE)
out_dir <- args[[1L]]
opts <- args[-1L]
optValue <- function(name, default = NULL) {
    i <- match(name, opts)
    if (is.na(i) || i == length(opts)) return(default)
    opts[[i + 1L]]
}
has <- function(name) name %in% opts

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
tmp_root <- Sys.getenv("RUNNER_TEMP", tempdir())
val_lib <- file.path(tmp_root, "val-lib")
dir.create(val_lib, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf("[val-tests] %s\n", paste0(...)))

# --- Install msrval from this checkout ---
r_bin <- file.path(R.home("bin"), "R")
install_log <- file.path(out_dir, "install.log")
status <- system2(r_bin, c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", val_lib),
                           pkg_dir), stdout = install_log, stderr = install_log)
say("install msrval: exit ", status)
if (!identical(as.integer(status), 0L)) quit(status = 2L)
.libPaths(c(val_lib, .libPaths()))
Sys.setenv(VAL_HARNESS_DIR = normalizePath(harness_dir))
src <- sort(list.files(pkg_dir, recursive = TRUE, full.names = TRUE))
utils::write.table(data.frame(path = sub(paste0("^", harness_dir, "/"), "", src),
                              sha256 = vapply(src, function(f) digest::digest(file = f,
                                                                             algo = "sha256"),
                                              character(1))),
                   file.path(out_dir, "code-hashes.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)
Sys.setenv(R_LIBS = paste(c(val_lib, .libPaths()), collapse = .Platform$path.sep))

# --- msrval's own tests ---
runTests <- function(test_path, package, tsv_name, failures_name) {
    started <- Sys.time()
    load <- if (is.null(package)) "none" else "installed"
    only <- optValue("--files")
    filter <- if (is.null(only)) NULL else paste0("^(", gsub(",", "|", only), ")$")
    res <- testthat::test_dir(test_path, package = package, load_package = load, filter = filter,
                              reporter = testthat::ListReporter$new(), stop_on_failure = FALSE,
                              stop_on_warning = FALSE)
    df <- as.data.frame(res)
    rows <- data.frame(
        file = df$file, test = df$test, expectations = df$nb, passed = df$nb - df$failed,
        failed = df$failed, skipped = as.integer(df$skipped), error = as.integer(df$error),
        warnings = df$warning, seconds = round(df$real, 3), stringsAsFactors = FALSE
    )
    utils::write.table(rows, file.path(out_dir, tsv_name), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    msgs <- character()
    for (r in res) {
        bad <- Filter(function(e) inherits(e, c("expectation_failure", "expectation_error")),
                      r$results)
        if (length(bad) == 0L) next
        first <- substr(conditionMessage(bad[[1L]]), 1L, 600L)
        msgs <- c(msgs, sprintf("== %s :: %s (%d not passed)\n%s\n", r$file, r$test,
                                length(bad), first))
    }
    writeLines(msgs, file.path(out_dir, failures_name))
    say(tsv_name, ": ", nrow(rows), " tests, ", sum(rows$failed > 0 | rows$error > 0),
        " not passing, ", sum(rows$skipped > 0), " skipped, ",
        round(as.numeric(difftime(Sys.time(), started, units = "secs"))), " s")
    rows
}

if (!is.null(optValue("--ticket"))) {
    ticket_dir <- optValue("--ticket")
    runTests(file.path(ticket_dir, "tests"), NULL, "ticket-tests.tsv", "ticket-failures.txt")
    # The declaration's coverage_regimes (4.4, I03): ticket-coverage.tsv and its cells file.
    decl_path <- optValue("--declaration")
    if (!is.null(decl_path)) {
        decl <- jsonlite::read_json(decl_path, simplifyVector = FALSE)
        tc <- msrval::valTicketCoverage(decl, dirname(decl_path))
        if (!is.null(tc)) {
            utils::write.table(tc$coverage, file.path(out_dir, "ticket-coverage.tsv"), sep = "\t",
                               quote = FALSE, row.names = FALSE, na = "NA")
            utils::write.table(tc$cells, file.path(out_dir, "ticket-coverage-cells.tsv"),
                               sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
            say("ticket coverage: ", nrow(tc$coverage), " rows")
        }
    }
} else {
    runTests(file.path(pkg_dir, "tests", "testthat"), "msrval", "harness-tests.tsv",
             "harness-failures.txt")
}

# --- Lint ---
if (has("--lint")) {
    linters <- lintr::linters_with_defaults(
        line_length_linter = lintr::line_length_linter(100L),
        indentation_linter = lintr::indentation_linter(indent = 4L),
        # constants are UPPER_SNAKE_CASE by conventions/R.md (naming table)
        object_name_linter = lintr::object_name_linter(styles = c("camelCase", "snake_case",
                                                                  "SNAKE_CASE"))
    )
    lints <- lintr::lint_dir(file.path(pkg_dir, "R"), linters = linters)
    lint_df <- as.data.frame(lints)
    if (nrow(lint_df) == 0L) {
        writeLines("no lints", file.path(out_dir, "lint.txt"))
    } else {
        utils::write.table(lint_df[, c("filename", "line_number", "linter", "message")],
                           file.path(out_dir, "lint.txt"), sep = "\t", quote = FALSE,
                           row.names = FALSE)
    }
    say("lintr: ", nrow(lint_df), " lints")
}

# --- Line coverage ---
if (has("--covr")) {
    # Lint runs first (it takes seconds), so a covr failure or a time-out cannot lose it; a covr
    # error is written to coverage.txt (its class only on the console) instead of ending the script.
    cov <- tryCatch(covr::package_coverage(pkg_dir, type = "tests", quiet = TRUE),
                    error = function(e) e)
    if (inherits(cov, "error")) {
        writeLines(c("error", substr(conditionMessage(cov), 1L, 2000L)),
                   file.path(out_dir, "coverage.txt"))
        say("covr: error (", class(cov)[[1L]], ")")
    } else {
        per_file <- covr::tally_coverage(cov, by = "line")
        by_file <- tapply(per_file$value > 0, per_file$filename, mean)
        lines <- c(sprintf("total\t%.2f", covr::percent_coverage(cov)),
                   sprintf("%s\t%.2f", names(by_file), 100 * by_file))
        writeLines(lines, file.path(out_dir, "coverage.txt"))
        say("covr: ", lines[[1L]])
    }
}

# --- Coverage regimes (T-COV) and the chunk speed-up on simulations (f_sim) ---
if (!is.null(optValue("--time-regime")) || !is.null(optValue("--probe-cells"))) {
    probe <- optValue("--probe-cells")
    msrval::valTimeRegime(optValue("--time-regime", ""),
                          nsim = as.integer(optValue("--nsim", "200")), out_dir = out_dir,
                          chunks = 4L,
                          probe_nsim = as.integer(optValue("--probe-nsim", "20")),
                          probe_target_sims = as.integer(optValue("--probe-target-sims", "200")),
                          probe_cells = if (is.null(probe)) NULL else strsplit(probe, ",")[[1L]])
    say("f_sim and cost probe written")
}
if (!is.null(optValue("--coverage-plan"))) {
    msrval::valRunCoverage(optValue("--coverage-plan"), request = optValue("--request", "all"),
                           nsim = as.integer(optValue("--nsim", "2000")),
                           chunks = as.integer(optValue("--chunks", "4")), out_dir = out_dir)
    say("coverage-table.tsv written")
}
