# KIT-082: how long each new worker display action took on a journey member, request written to journal last written. Times only.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-082/site-durations.R <member dir>...
sites <- "_(filtering_metrics|final_pca|volcano_page|session_export|integration_export)$"
for (member in commandArgs(TRUE)) {
    requests <- list.files(member, "^request[.]rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
    requests <- requests[grepl("/[.]workflow-export/wfx-", requests)]
    for (path in requests) {
        request <- tryCatch(readRDS(path), error = function(error) NULL)
        if (is.null(request$operation) || !grepl(sites, request$operation)) next
        seconds <- as.numeric(difftime(file.mtime(file.path(dirname(path), "journal.json")), file.mtime(path), units = "secs"))
        cat(sprintf("%s\t%s\t%.1f\n", basename(dirname(member)), request$operation, seconds))
    }
}
