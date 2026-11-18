# run-tests: the shards' tables, failures and logs assembled into one results directory, and the summary.
#
# Usage: Rscript --vanilla assemble.R <downloaded dir> <results dir> <requested files, one per line>
# Reads <downloaded dir>/shard-*/ (tests.tsv, failures/, logs/, withheld.txt), check/ and scripts/ when present.
# Writes <results dir>/tests.tsv (in the requested order; a requested file without a row, its shard job lost, is
# recorded as crashed), failures.txt, logs/, check-summary.txt, scripts/, withheld.txt and summary.md.
# Exit status 1 when a file crashed or timed out or a script failed; 0 otherwise (failing expectations are results).

arguments <- commandArgs(trailingOnly = TRUE)
downloaded <- arguments[[1L]]; results <- arguments[[2L]]
requested <- readLines(arguments[[3L]], warn = FALSE)
requested <- requested[nzchar(requested)]
dir.create(file.path(results, "logs"), recursive = TRUE, showWarnings = FALSE)

shards <- sort(list.files(downloaded, pattern = "^shard-[0-9]+$", full.names = TRUE))
read <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, quote = "", na.strings = character(),
    colClasses = c(file = "character", status = "character"))
rows <- do.call(rbind, lapply(shards, function(shard) {
    path <- file.path(shard, "tests.tsv")
    if (file.exists(path)) read(path)
}))
columns <- c("file", "passed", "failed", "errors", "skipped", "warnings", "seconds", "status")
if (is.null(rows)) rows <- stats::setNames(data.frame(matrix(nrow = 0L, ncol = 8L)), columns)
lost <- setdiff(requested, rows$file)
if (length(lost)) {
    rows <- rbind(rows, data.frame(file = lost, passed = 0L, failed = 0L, errors = 0L, skipped = 0L, warnings = 0L,
        seconds = NA, status = "crashed", stringsAsFactors = FALSE))
}
rows <- rows[match(requested, rows$file), columns]
utils::write.table(rows, file.path(results, "tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

failures <- character()
for (file in requested) {
    for (shard in shards) {
        path <- file.path(shard, "failures", sub("[.]R$", ".txt", file))
        if (file.exists(path)) failures <- c(failures, readLines(path, warn = FALSE))
    }
    if (file %in% lost) failures <- c(failures, sprintf("file: %s", file), "test: <shard job>", "kind: crashed",
        "message:", "  no result: the shard job that held this file failed or was cancelled (see the run's job log)", "")
}
writeLines(if (length(failures)) failures else "no failing expectations", file.path(results, "failures.txt"))
withheld <- character()
for (shard in shards) {
    logs <- list.files(file.path(shard, "logs"), full.names = TRUE)
    file.copy(logs, file.path(results, "logs"), overwrite = TRUE)
    if (file.exists(file.path(shard, "withheld.txt"))) withheld <- c(withheld, readLines(file.path(shard, "withheld.txt")))
}
for (part in c("check", "scripts")) {
    source_dir <- file.path(downloaded, part)
    if (!dir.exists(source_dir)) next
    file.copy(source_dir, results, recursive = TRUE)
    if (file.exists(file.path(source_dir, "withheld.txt"))) {
        withheld <- c(withheld, paste0(part, "/", readLines(file.path(source_dir, "withheld.txt"))))
    }
}
if (file.exists(file.path(results, "check", "check-summary.txt"))) {
    file.copy(file.path(results, "check", "check-summary.txt"), file.path(results, "check-summary.txt"), overwrite = TRUE)
}
unlink(c(file.path(results, "check", "withheld.txt"), file.path(results, "scripts", "withheld.txt")))
if (length(withheld)) writeLines(withheld, file.path(results, "withheld.txt"))

scripts <- if (file.exists(file.path(results, "scripts", "scripts.tsv"))) {
    utils::read.delim(file.path(results, "scripts", "scripts.tsv"), stringsAsFactors = FALSE, quote = "")
}
bad <- rows$status %in% c("crashed", "timeout")
scripts_bad <- !is.null(scripts) && any(scripts$status != "ok")

md <- c(sprintf("### Tests: %d files; passed %d, failed %d, errors %d, skipped %d; ok %d, failed %d, timeout %d, crashed %d",
        nrow(rows), sum(rows$passed), sum(rows$failed), sum(rows$errors), sum(rows$skipped), sum(rows$status == "ok"),
        sum(rows$status == "failed"), sum(rows$status == "timeout"), sum(rows$status == "crashed")), "",
    "| file | passed | failed | errors | skipped | warnings | seconds | status |", "|---|---:|---:|---:|---:|---:|---:|---|",
    sprintf("| %s | %d | %d | %d | %d | %d | %s | %s |", rows$file, rows$passed, rows$failed, rows$errors, rows$skipped,
        rows$warnings, ifelse(is.na(rows$seconds), "", rows$seconds), rows$status), "")
if (file.exists(file.path(results, "check-summary.txt"))) {
    md <- c(md, "### R CMD check", "```", utils::head(readLines(file.path(results, "check-summary.txt")), 60L), "```", "")
}
if (!is.null(scripts)) {
    md <- c(md, "### Scripts", "| # | path | exit | seconds | status |", "|---:|---|---:|---:|---|",
        sprintf("| %s | %s | %s | %s | %s |", scripts$index, scripts$path, scripts$exit, scripts$seconds, scripts$status), "")
}
if (length(withheld)) md <- c(md, sprintf("Withheld from the artifact (matched a forbidden pattern): %d files; see withheld.txt", length(withheld)), "")
md <- c(md, if (any(bad) || scripts_bad) "**Run outcome: failure** (a file crashed or timed out, or a script failed)" else
    "**Run outcome: success** (every file ran to completion; failing expectations, if any, are in the table)")
writeLines(md, file.path(results, "summary.md"))
cat(md, sep = "\n")
quit(save = "no", status = if (any(bad) || scripts_bad) 1L else 0L)
