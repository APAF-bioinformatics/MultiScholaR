# KIT-203 S5: a record made on the base commit with the real kernel replays on this commit with the q-value, its
# declared names, values, product and files unchanged. Evidence code under docs/qa: no package code. Runs on a
# GitHub-hosted runner (run-tests `scripts`) after the commit under test (head) is installed, from the package root.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-203/replay-k203.R <out dir> <base commit>
#
# 1. The base commit's tree is installed into a private library; KIT-173's capture-record.R runs on it (one proteomics
#    project on the clusterProfiler route, taxon 999999, the T7 fixture, the tab's default propagation; no service
#    call), writing the project store, the record and the evidence.
# 2. On head (this process), the record is replayed on a copy of that store (helper goPropReplay, the kernel the GUI
#    selects outside test mode). The base record is made after VAL-61 (geneID cells sorted), so every table and every
#    results TSV is compared exactly.
# Writes <out dir>/s5.tsv (check, pass, detail) and s5-base.tsv (which code made the capture). Exit status 0.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
out <- normalizePath(args[[1L]], mustWork = FALSE)
base <- args[[2L]]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k203-s5")
src <- file.path(work, "base-src")
lib <- file.path(work, "base-lib")
capture <- file.path(work, "capture")
for (d in c(src, lib, capture)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
base_check <- list(base_commit = base)
base_check$archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
base_check$install_exit <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load",
    paste0("--library=", shQuote(lib)), shQuote(src)), stdout = file.path(work, "install.log"),
    stderr = file.path(work, "install.log"))
probe <- system2(rscript, c("-e", shQuote(paste("suppressPackageStartupMessages(library(MultiScholaR));",
    "cat(find.package('MultiScholaR'), 'significance_rule' %in% names(formals(processEnrichments)), sep = '\\t')"))),
    stdout = TRUE, env = paste0("R_LIBS=", lib))
base_check$base_package_path_is_private_library <- startsWith(probe[[length(probe)]], lib)
base_check$base_has_significance_rule <- sub("^.*\t", "", probe[[length(probe)]])
base_check$capture_exit <- system2(rscript, c(file.path("docs", "qa", "workflow-execution", "kit", "kit-173",
    "capture-record.R"), shQuote(capture)), stdout = file.path(work, "capture.log"), stderr = file.path(work, "capture.log"),
    env = paste0("R_LIBS=", lib))
utils::write.table(data.frame(measure = names(base_check), value = vapply(base_check, as.character, "")),
    file.path(out, "s5-base.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(base_check)

checks <- list()
check <- function(id, pass, detail = "") checks[[length(checks) + 1L]] <<- data.frame(check = id, pass = isTRUE(pass),
    detail = substr(gsub("[\t\r\n]+", " ", detail), 1L, 300L), stringsAsFactors = FALSE)
result <- tryCatch({
    suppressPackageStartupMessages(library(MultiScholaR))
    Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
    h <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
    scratch <- file.path(work, "replay")
    project <- file.path(scratch, "origin")
    dir.create(project, recursive = TRUE)
    utils::untar(file.path(capture, "project.tar"), exdir = project, tar = "internal")
    record <- jsonlite::read_json(file.path(capture, "record.json"), simplifyVector = FALSE)
    base_evidence <- readRDS(file.path(capture, "evidence.rds"))
    base_product <- readRDS(file.path(capture, "product.rds"))
    replayed <- NULL
    utils::capture.output(replayed <- suppressMessages(suppressWarnings(h$goPropReplay(project, record,
        file.path(scratch, "run")))))
    step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
    recorded_names <- vapply(step$parameters, `[[`, "", "name")
    check("states", identical(replayed$states, c("initial", "correlation_filtered", "enrichment_completed")),
        paste(replayed$states, collapse = ","))
    check("record declares no significance_rule", !"significance_rule" %in% recorded_names,
        paste(recorded_names, collapse = ","))
    check("declared names = base's", identical(names(replayed$provenance$values), base_product$provenance_names),
        paste(names(replayed$provenance$values), collapse = ","))
    check("declared values (product parameters) = base's", identical(replayed$product$parameters,
        base_product$parameters))
    check("service = base's", identical(replayed$product$service, base_product$service),
        paste(replayed$product$service, collapse = ","))
    check("contrasts = base's", identical(replayed$product$contrasts, base_product$contrasts))
    evidence <- h$goPropRunEvidence(replayed$state, replayed$pathway_dir)
    check("result tables identical (exact)", identical(evidence$tables, base_evidence$tables))
    check("plot inputs identical", identical(evidence$plot_inputs, base_evidence$plot_inputs))
    check("pathway file names", identical(evidence$files$file, base_evidence$files$file),
        paste(evidence$files$file, collapse = ","))
    check("every file's SHA-256 (TSV) or existence", identical(evidence$files$sha256, base_evidence$files$sha256))
    tables <- unlist(lapply(evidence$tables, function(contrast) lapply(contrast, function(t) if (!is.null(t)) names(t))))
    check("no family column in the replayed tables", !any(c("p_adjust_family", "family_size") %in% tables))
    "ok"
}, error = function(e) paste("instrument:", conditionMessage(e)))
table <- if (length(checks)) do.call(rbind, checks) else data.frame(check = character(), pass = logical(),
    detail = character())
if (!identical(result, "ok")) table <- rbind(table, data.frame(check = "script", pass = FALSE, detail = result))
utils::write.table(table, file.path(out, "s5.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table)
