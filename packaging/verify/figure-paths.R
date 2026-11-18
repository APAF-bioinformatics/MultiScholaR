# KIT-134 / KIT-052 (Windows figure paths): one case of the figure-path confirmation with the installed
# application's R. Build-time verification, not package code; nothing in R/ or tests/ is changed.
#
# Usage: Rscript --vanilla figure-paths.R <checkout> <test-tool library> <output directory> <case> <item>...
# Items:
#   devices        each way a figure or file is written (grDevices::png of each type, ragg::agg_png,
#                  grDevices::pdf, ggplot2::ggsave, writeLines, file.copy) to paths of 240 to 300
#                  characters, in a short folder (MSR_PROBE_DIR) with long file names, so that only the
#                  file's path length varies; whether the file is then in the folder and not empty
#   figure_tests   tests/testthat/test-workflow-figure-paths.R (KIT-052's tests), whole file
#   trx_worker_da  test-module-ci-trx-da.R, MCI-031.2: transcriptomics DA on the managed route, the worker's
#                  entry point and the acceptance in this process; the worker draws the DA figures in its
#                  request folder (<store>/.workflow-export/wfx-<32>/da-publication|omics-publication)
#   prot_worker_da test-omic-da-module.R, KIT-112 W3, one contrast: proteomics DA on the managed route
#                  (DA-contrast-output, then da_publication in the worker's request folder)
#   prot_tmt_lane  the GUI-entry lane prot_tmt (QC with the Q-value figures, DA in process), against its oracle
# The project root's length is set by the caller through the temporary folder (TMPDIR), under which every
# fixture makes its project (withr::local_tempdir). Each item's figures are printed by figure-trace.R.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 5L) stop("usage: Rscript --vanilla figure-paths.R <checkout> <test-tool library> <output directory> <case> <item>...")
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
tools_lib <- normalizePath(arguments[[2L]], mustWork = TRUE)
out_dir <- arguments[[3L]]
case <- arguments[[4L]]
items <- arguments[-(1:4)]
dir.create(file.path(out_dir, "records"), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
.libPaths(c(.libPaths(), tools_lib))
suppressPackageStartupMessages(library(MultiScholaR))
source(file.path(checkout, "packaging", "verify", "figure-trace.R"))
cat(sprintf("CASE %s: R %s; tempdir %s (%d characters); TMPDIR %s\n", case, getRversion(), tempdir(), nchar(tempdir()),
    Sys.getenv("TMPDIR")))
if (.Platform$OS.type == "windows") {
    value <- tryCatch(utils::readRegistry("SYSTEM\\CurrentControlSet\\Control\\FileSystem", "HLM")$LongPathsEnabled,
        error = function(error) conditionMessage(error))
    cat(sprintf("CASE %s: LongPathsEnabled as this R process reads it: %s\n", case, paste(value, collapse = " ")))
}
kit134FigureTraceInstall()
Sys.setenv(NOT_CRAN = "true", MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = file.path(out_dir, "records"))

runBlock <- function(file, desc) {
    reporter <- testthat::ListReporter$new()
    outcome <- tryCatch({
        testthat::test_file(file.path(checkout, "tests", "testthat", file), desc = desc, reporter = reporter,
            stop_on_failure = FALSE, package = "MultiScholaR", load_package = "installed")
        ""
    }, error = function(error) conditionMessage(error))
    results <- unlist(lapply(reporter$get_results(), function(test) test$results), recursive = FALSE)
    classes <- vapply(results, function(x) class(x)[[1L]], character(1))
    cat(sprintf("TEST %s %s: pass %d fail %d error %d skip %d%s\n", case, .kit134FigureTrace$item, sum(classes == "expectation_success"),
        sum(classes == "expectation_failure"), sum(classes == "expectation_error"), sum(classes == "expectation_skip"),
        if (nzchar(outcome)) paste0("; run error: ", outcome) else ""))
    for (x in utils::head(results[!classes %in% "expectation_success"], 4L))
        cat("  ", class(x)[[1L]], ":", substr(gsub("[\r\n\t]+", " ", conditionMessage(x)), 1L, 700L), "\n")
}

probeDevices <- function() {
    base <- Sys.getenv("MSR_PROBE_DIR")
    folder <- file.path(base, case)
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    folder <- gsub("\\\\", "/", normalizePath(folder))
    source <- file.path(folder, "source.txt")
    writeLines("source", source)
    figure <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
    writers <- list(
        "grDevices::png type=windows" = function(path) { grDevices::png(path, type = "windows"); plot(1:3); grDevices::dev.off() },
        "grDevices::png type=cairo" = function(path) { grDevices::png(path, type = "cairo"); plot(1:3); grDevices::dev.off() },
        "ragg::agg_png" = function(path) { ragg::agg_png(path); plot(1:3); grDevices::dev.off() },
        "grDevices::pdf" = function(path) { grDevices::pdf(path); plot(1:3); grDevices::dev.off() },
        "ggplot2::ggsave png" = function(path) ggplot2::ggsave(path, figure, width = 3, height = 3, dpi = 72),
        "ggplot2::ggsave pdf" = function(path) ggplot2::ggsave(path, figure, width = 3, height = 3),
        "writeLines" = function(path) writeLines("text", path),
        "file.copy" = function(path) file.copy(source, path))
    rows <- list()
    for (target in c(240L, 255L, 259L, 260L, 261L, 270L, 300L)) {
        # No name over 150 characters (Windows allows 255 per name): a subfolder takes the rest, and it
        # stays shorter than 248 characters, so it can be made and listed without long paths.
        rest <- target - nchar(folder) - 1L
        directory <- folder
        if (rest > 150L) {
            directory <- file.path(folder, paste0("d", target, strrep("y", rest - 152L - nchar(paste0("d", target)))))
            dir.create(directory, showWarnings = FALSE)
            rest <- 150L
        }
        for (writer in names(writers)) {
            extension <- if (grepl("pdf", writer)) ".pdf" else if (writer %in% c("writeLines", "file.copy")) ".txt" else ".png"
            stem <- sprintf("%s-%d-", gsub("[^A-Za-z0-9]+", "_", writer), target)
            name <- paste0(stem, strrep("x", rest - nchar(stem) - nchar(extension)), extension)
            path <- file.path(directory, name)
            message <- ""
            withCallingHandlers(tryCatch(writers[[writer]](path), error = function(error) message <<- paste("error:", conditionMessage(error))),
                warning = function(warning) {
                    message <<- paste(message, "warning:", conditionMessage(warning))
                    invokeRestart("muffleWarning")
                })
            while (grDevices::dev.cur() > 1L) grDevices::dev.off()
            listed <- name %in% list.files(directory)
            bytes <- suppressWarnings(file.size(path))
            rows[[length(rows) + 1L]] <- data.frame(case = case, writer = writer, length = nchar(path), exists = listed, bytes = bytes,
                condition = substr(gsub("[\r\n\t]+", " ", trimws(message)), 1L, 200L), stringsAsFactors = FALSE)
        }
    }
    rows <- do.call(rbind, rows)
    utils::write.table(rows, file.path(out_dir, sprintf("devices-%s.tsv", case)), sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("DEVICES %s: folder %s (%d characters)\n", case, folder, nchar(folder)))
    cat(paste("DEVROW", apply(rows, 1L, function(r) paste(trimws(r), collapse = "\t"))), sep = "\n")
}

for (item in items) {
    .kit134FigureTrace$item <- item
    started <- Sys.time()
    switch(item,
        devices = probeDevices(),
        figure_tests = runBlock("test-workflow-figure-paths.R", NULL),
        trx_worker_da = runBlock("test-module-ci-trx-da.R",
            "MCI-031.2 several pairwise contrasts run on one fit and are offered by the tab; the managed tab equals the in-process one"),
        prot_worker_da = runBlock("test-omic-da-module.R",
            "KIT-112 W3 (D3): a managed run with one contrast dispatches DA-contrast-output then da_publication and gives the in-process run's manifest, tab and bus"),
        # The lipid and metabolomics worker routes are the "(b-managed)" cases of test-omic-da-module.R, made in a
        # top-level loop: testthat's desc filter runs every top-level call that is not a test_that block, so naming
        # the one top-level lipid managed test below runs that test plus the file's loops ((b) in process and
        # (b-managed) lipid, metab and trx), not the whole file.
        omic_worker_da = runBlock("test-omic-da-module.R",
            "(W3, D3) lipid managed: the volcano is <req> until ticked, then the worker's lipid_volcano_page iframe; the tab holds no state, widget or page"),
        prot_tmt_lane = {
            lines <- readLines(file.path(checkout, "tests", "testthat", "test-gui-entry-proteomics.R"))
            block <- gsub("^.*test_that\\(\"|\"$", "", regmatches(lines, regexpr("test_that\\(\"KIT-101 entry lane prot_tmt:[^\"]+\"", lines)))
            runBlock("test-gui-entry-proteomics.R", block)
        },
        cat("unknown item", item, "\n"))
    cat(sprintf("ITEM %s %s: %.0f s\n", case, item, as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
kit134FigureTraceReport(file.path(out_dir, sprintf("figures-%s.tsv", case)), case)
