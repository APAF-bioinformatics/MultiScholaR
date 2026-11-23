# KIT-134 P5/P6 (release-candidate run, branch run-k134-rc): one shard of the package's test suite ON THE BUNDLE: the
# bundle's own R, library and installed MultiScholaR (testthat::test_file with load_package = "installed"), not
# pkgload::load_all. Build-time verification, not package code.
#
# Usage: Rscript --vanilla run-suite.R <checkout> <test-tool library> <out dir> <shard> <shards> <deadline, epoch s>
#        <per-file timeout, minutes>
# Run with the bundle's launcher (Linux bin/Rscript; Windows packaging/verify/windows-clean-run.ps1 -NoBundleArgument),
# so every child R is the bundle's R in the bundle's environment.
#
# The files of tests/testthat are dealt to the shards by the longest-first rule over the seconds each file took in the
# runner reference run 44438735444 (packaging/verify/suite-seconds.tsv; 60 s for a file not listed), the same deal on
# every shard. Each file runs in its own R process (as the laptop's per-file runs and the run-tests workflow), killed
# at the per-file timeout; a file that would start after the deadline is not started (status timeout). The test tools
# (testthat and sessioninfo, as in the run-tests library) are a separate library after the bundle's; R_LIBS carries
# both to the processes the tests start.
# Prints SUITE lines (one per file: file, passed, failed, errors, skipped, warnings, seconds, status) and SUITEFAIL
# lines (file, test, kind, the message's first 600 characters on one line); writes <out dir>/tests.tsv.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) == 8L) arguments <- arguments[-1L] # a runner that passes the bundle folder first
if (length(arguments) != 7L) stop("usage: Rscript --vanilla run-suite.R <checkout> <tools> <out> <shard> <shards> <deadline> <timeout min>")
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
tools_lib <- normalizePath(arguments[[2L]], mustWork = TRUE)
out_dir <- arguments[[3L]]
shard <- as.integer(arguments[[4L]])
shards <- as.integer(arguments[[5L]])
deadline <- as.numeric(arguments[[6L]])
per_file <- as.numeric(arguments[[7L]])
for (folder in c("rows", "failures", "logs", "results")) dir.create(file.path(out_dir, folder), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)

bundle_libs <- .libPaths()
shadowing <- intersect(rownames(utils::installed.packages(tools_lib)), rownames(utils::installed.packages(bundle_libs)))
if (length(shadowing)) {
    cat("the test-tool library holds packages the bundle also has:", paste(shadowing, collapse = ", "), "\n")
    quit(save = "no", status = 2L)
}
cat(sprintf("SUITE SETUP R %s at %s; library %s; test tools %s (%s)\n", getRversion(), R.home(),
    paste(bundle_libs, collapse = " | "), tools_lib, paste(sort(rownames(utils::installed.packages(tools_lib))), collapse = ", ")))
cat(sprintf("SUITE SETUP MultiScholaR %s from %s\n", utils::packageVersion("MultiScholaR"), find.package("MultiScholaR")))

files <- sort(list.files(file.path(checkout, "tests", "testthat"), pattern = "^test-.*[.]R$"), method = "radix")
# run-win-k147: MSR_SUITE_FILES (comma-separated) limits the run to those files; a name with a "/" is a path under the
# checkout (the lifecycle probe), any other a file of tests/testthat. A name that does not exist stops the run.
selected <- strsplit(Sys.getenv("MSR_SUITE_FILES", ""), ",", fixed = TRUE)[[1L]]
selected <- trimws(selected[nzchar(trimws(selected))])
if (length(selected)) {
    missing <- selected[!file.exists(ifelse(grepl("/", selected, fixed = TRUE), file.path(checkout, selected),
        file.path(checkout, "tests", "testthat", selected)))]
    if (length(missing)) stop("MSR_SUITE_FILES names files that do not exist: ", paste(missing, collapse = ", "))
    files <- sort(selected, method = "radix")
    cat(sprintf("SUITE SETUP selected %d files: %s\n", length(files), paste(files, collapse = ", ")))
}
# run-win-k147: MSR_SUITE_TMP gives each file a fresh, empty TMPDIR/TMP/TEMP (<root>/<n>); what is left in it after the
# file's process ends is printed as SUITETMP lines (staging or session folders not removed).
tmp_root <- Sys.getenv("MSR_SUITE_TMP", "")
reference <- utils::read.delim(file.path(checkout, "packaging", "verify", "suite-seconds.tsv"), stringsAsFactors = FALSE)
seconds <- reference$seconds[match(files, reference$file)]
seconds[is.na(seconds)] <- 60
order <- order(-seconds, files, method = "radix")
load <- numeric(shards)
assigned <- integer(length(files))
for (i in order) {
    k <- which.min(load)
    assigned[[i]] <- k
    load[[k]] <- load[[k]] + seconds[[i]]
}
mine <- files[assigned == shard]
cat(sprintf("SUITE SETUP shard %d of %d: %d of %d files, %.0f reference minutes (shards: %s)\n", shard, shards,
    length(mine), length(files), load[[shard]] / 60, paste(round(load / 60), collapse = " ")))

Sys.setenv(R_LIBS = paste(c(bundle_libs, tools_lib), collapse = .Platform$path.sep), NOT_CRAN = "true")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
child <- file.path(checkout, "packaging", "verify", "run-suite-file.R")
row <- function(file, status, seconds) data.frame(file = file, passed = 0L, failed = 0L, errors = 0L, skipped = 0L,
    warnings = 0L, seconds = seconds, status = status, stringsAsFactors = FALSE)
rows <- list()
for (file in mine) {
    stem <- sub("[.]R$", "", basename(file))
    path <- if (grepl("/", file, fixed = TRUE)) file.path(checkout, file) else file.path(checkout, "tests", "testthat", file)
    if (nzchar(tmp_root)) {
        file_tmp <- file.path(tmp_root, match(file, mine))
        dir.create(file_tmp, recursive = TRUE, showWarnings = FALSE)
        file_tmp <- normalizePath(file_tmp)
        Sys.setenv(TMPDIR = file_tmp, TMP = file_tmp, TEMP = file_tmp)
    }
    log <- file.path(out_dir, "logs", paste0(stem, ".log"))
    left <- (deadline - as.numeric(Sys.time())) / 60
    if (left < 3) {
        rows[[file]] <- row(file, "timeout", 0)
        cat(sprintf("SUITE NOTE %s not started: the job's time is spent\n", file))
        next
    }
    limit <- min(per_file, left - 1)
    started <- Sys.time()
    outcome <- tryCatch(processx::run(rscript, c("--vanilla", child, checkout, tools_lib, out_dir,
        path), error_on_status = FALSE, timeout = limit * 60,
        stdout = log, stderr = "2>&1", cleanup_tree = TRUE, wd = checkout),
        system_command_timeout_error = function(e) list(status = NA_integer_, timeout = TRUE),
        error = function(e) list(status = NA_integer_, timeout = FALSE, error = conditionMessage(e)))
    took <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
    row_file <- file.path(out_dir, "rows", paste0(stem, ".tsv"))
    if (isTRUE(outcome$timeout)) {
        unlink(row_file)
        rows[[file]] <- row(file, "timeout", took)
    } else if (!file.exists(row_file)) {
        rows[[file]] <- row(file, "crashed", took)
        tail <- if (file.exists(log)) utils::tail(readLines(log, warn = FALSE), 15L) else "no log"
        cat(sprintf("SUITEFAIL\t%s\t<process>\tcrashed\texit %s%s; last lines: %s\n", file, paste(outcome$status),
            if (!is.null(outcome$error)) paste0(" (", outcome$error, ")") else "",
            substr(gsub("[\t\r\n]+", " ", paste(tail, collapse = " / ")), 1L, 1500L)))
    } else {
        rows[[file]] <- utils::read.delim(row_file, stringsAsFactors = FALSE)
    }
    r <- rows[[file]]
    cat(sprintf("SUITE\t%s\t%d\t%d\t%d\t%d\t%d\t%.1f\t%s\n", r$file, r$passed, r$failed, r$errors, r$skipped, r$warnings,
        took, r$status))
    if (nzchar(tmp_root)) {
        left <- list.files(file_tmp, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)
        cat(sprintf("SUITETMP\t%s\t%d entries left in %s\n", file, length(left), file_tmp))
        for (entry in utils::head(left, 40L)) cat(sprintf("SUITETMP\t%s\t%s\n", file, entry))
    }
    failures <- file.path(out_dir, "results", paste0(stem, ".rds"))
    if (file.exists(failures)) {
        messages <- readRDS(failures)$messages
        for (j in seq_len(nrow(messages))) cat(sprintf("SUITEFAIL\t%s\t%s\t%s\t%s\n", file,
            substr(gsub("[\t\r\n]+", " ", messages$test[[j]]), 1L, 200L), messages$kind[[j]],
            substr(gsub("[\t\r\n]+", " ", messages$message[[j]]), 1L, 600L)))
    }
}
suite_table <- do.call(rbind, rows)
utils::write.table(suite_table, file.path(out_dir, "tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("SUITE TOTAL shard %d: %d files; passed %d, failed %d, errors %d, skipped %d; status %s\n", shard,
    nrow(suite_table), sum(suite_table$passed), sum(suite_table$failed), sum(suite_table$errors), sum(suite_table$skipped),
    paste(names(table(suite_table$status)), table(suite_table$status), sep = " ", collapse = ", ")))
quit(save = "no", status = 0L)
