# KIT-134 P5 and P6: the browser e2e lanes (tests/testthat/test-e2e-<lane>.R, shinytest2 and the
# runner's Chrome) against the bundle's installed MultiScholaR and library, with the test tools in a
# separate library after the bundle's. Build-time verification, not package code.
#
# Usage: Rscript --vanilla run-e2e.R <checkout> <test-tool library> <output directory> <lane>...
# <lane> is the test file's name without "test-e2e-" and ".R" (lipidomics-canonical, ...). Run with the
# bundle's launcher, with MULTISCHOLAR_E2E_INSTALLED=1 set here: the application the browser drives is
# library(MultiScholaR) from the bundle, not the source tree (tests/testthat/helper-e2e-browser.R).
#
# Each lane is one test_file() of the installed package. Written: e2e-summary.tsv (passed, failed,
# errors, skipped, and the laptop's last recorded count, KIT-118's exit at e94b92b1) and
# e2e-failures.tsv (text only). Exits 1 when a lane fails, errors or passes nothing; a pass
# count that differs from the laptop's recorded one is printed and written, not decided here (the
# record is of an earlier tree). Exits 2 when the lanes could not be run.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 4L) stop("usage: Rscript --vanilla run-e2e.R <checkout> <test-tool library> <output directory> <lane>...")
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
tools_lib <- normalizePath(arguments[[2L]], mustWork = TRUE)
out_dir <- arguments[[3L]]
lanes <- arguments[-(1:3)]
dir.create(file.path(out_dir, "artifacts"), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
LAPTOP <- c("lipidomics-canonical" = 88L, "metabolomics-lc-gc" = 190L, "transcriptomics-screen" = 16L,
    "proteomics-dia" = 42L, "proteomics-dia-limpa" = 86L, "proteomics-enrichment" = 98L, "proteomics-tmt-lfq" = 153L)

bundle_libs <- .libPaths()
.libPaths(c(bundle_libs, tools_lib))
shadowing <- intersect(rownames(utils::installed.packages(tools_lib)), rownames(utils::installed.packages(bundle_libs)))
cat(sprintf("library paths: %s\n", paste(.libPaths(), collapse = " | ")))
if (length(shadowing)) {
    cat("the test-tool library holds packages the bundle also has:", paste(shadowing, collapse = ", "), "\n")
    quit(save = "no", status = 2L)
}
suppressPackageStartupMessages(library(MultiScholaR))
for (package in c("MultiScholaR", "testthat", "shinytest2", "chromote", "fastduckdb", "fastlimpadpc", "duckdb")) {
    path <- find.package(package, quiet = TRUE)
    cat(if (length(path)) sprintf("%s %s (%s)", package, utils::packageVersion(package), dirname(path)) else paste(package, "absent"), "\n")
}
chrome <- tryCatch(as.character(chromote::find_chrome()), error = function(e) "")
cat(sprintf("browser: %s\nlocale: %s\n", if (length(chrome) && nzchar(chrome)) chrome else "none found", Sys.getlocale()))
if (!length(chrome) || !nzchar(chrome)) quit(save = "no", status = 2L)
reached <- tryCatch({
    close(socketConnection("cran.r-project.org", 80L, open = "r+", blocking = TRUE, timeout = 10))
    TRUE
}, error = function(e) FALSE, warning = function(w) FALSE)
cat(sprintf("network: %s\n", if (reached) "OPEN (cran.r-project.org reached)" else "closed (cran.r-project.org not reached)"))
Sys.setenv(MULTISCHOLAR_E2E_INSTALLED = "1", NOT_CRAN = "true", MULTISCHOLAR_E2E_ARTIFACT_DIR = file.path(out_dir, "artifacts"))
truncate <- function(x, n = 1500L) {
    x <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(x)), collapse = " "))
    if (nchar(x) > n) paste0(substr(x, 1L, n), "...") else x
}

summaries <- list()
failures <- list()
for (lane in lanes) {
    file <- file.path(checkout, "tests", "testthat", sprintf("test-e2e-%s.R", lane))
    if (!file.exists(file)) {
        cat("no test file for lane", lane, "\n")
        quit(save = "no", status = 2L)
    }
    started <- Sys.time()
    reporter <- testthat::ListReporter$new()
    outcome <- tryCatch({
        testthat::test_file(file, reporter = reporter, stop_on_failure = FALSE, package = "MultiScholaR", load_package = "installed")
        ""
    }, error = function(e) conditionMessage(e))
    expectations <- unlist(lapply(reporter$get_results(), function(test) test$results), recursive = FALSE)
    classes <- vapply(expectations, function(x) class(x)[[1L]], character(1))
    for (x in expectations[!classes %in% "expectation_success"]) {
        failures[[length(failures) + 1L]] <- data.frame(lane = lane, class = class(x)[[1L]],
            message = truncate(conditionMessage(x)), stringsAsFactors = FALSE)
    }
    laptop <- if (lane %in% names(LAPTOP)) LAPTOP[[lane]] else NA_integer_
    summaries[[lane]] <- data.frame(lane = lane, passed = sum(classes == "expectation_success"),
        failed = sum(classes == "expectation_failure"), errors = sum(classes == "expectation_error"),
        skipped = sum(classes == "expectation_skip"), laptop_passed = laptop,
        seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
        run_error = truncate(outcome, 300L), stringsAsFactors = FALSE)
    cat(sprintf("E2E %s: pass %d fail %d error %d skip %d; the laptop's recorded pass count %s; %.0f s\n", lane,
        summaries[[lane]]$passed, summaries[[lane]]$failed, summaries[[lane]]$errors, summaries[[lane]]$skipped,
        laptop, summaries[[lane]]$seconds))
}
summary <- do.call(rbind, summaries)
utils::write.table(summary, file.path(out_dir, "e2e-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(if (length(failures)) do.call(rbind, failures) else data.frame(lane = character(), class = character(), message = character()),
    file.path(out_dir, "e2e-failures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (x in utils::head(failures, 12L)) cat("  ", x$lane, x$class, ":", substr(x$message, 1L, 600L), "\n")
passed <- all(summary$failed == 0L & summary$errors == 0L & summary$passed > 0L & !nzchar(summary$run_error))
differs <- summary$lane[!is.na(summary$laptop_passed) & summary$passed != summary$laptop_passed]
if (length(differs)) cat("pass counts that differ from the laptop's recorded ones:", paste(differs, collapse = ", "), "\n")
cat(if (passed) "KIT-134 E2E PASSED" else "KIT-134 E2E FAILED", "\n")
quit(save = "no", status = if (passed) 0L else 1L)
