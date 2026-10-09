# KIT journey member driver entry (KIT-005 correction 3). One fresh process per journey member, so the
# pinned browser helpers' wait budget (which counts this process's elapsed time against the diagnostic
# ceiling) starts at zero for every member and covers the whole journey ceiling. The browser and the
# in-run fingerprints run here, outside the measured slice; the coordinator waits for this process.
# Usage: Rscript workflow_kit_member.R <library> <snapshot> <member_dir>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
loadNamespace("MultiScholaR")
api <- function(name) get(name, asNamespace("MultiScholaR"))
snapshot <- args[[2L]]
member_dir <- args[[3L]]
helpers <- new.env(parent = globalenv())
sys.source(file.path(snapshot, "tools/profiling/workflow_checkpoint_dia_driver.R"), helpers)
sys.source(file.path(snapshot, "tools/profiling/workflow_kit.R"), helpers)
sys.source(file.path(snapshot, "tools/profiling/workflow_kit_driver.R"), helpers)
request <- readRDS(file.path(member_dir, "driver-request.rds"))
stopifnot(identical(request$schema, 1L), is.numeric(request$deadline_uptime))
# The pinned helpers' wait budget is a diagnostic ceiling (at most 2 h) minus this process's elapsed time;
# a 300-sample journey runs longer than that. Every helper resolves `checkpointWaitBudget` in this
# environment, so the member's own ceiling (the sampler's time guard) is the budget instead
# (KIT-005 correction 5); control waits (element present, tab switched) get up to 120 s.
helpers$checkpointWaitBudget <- function(control = FALSE, elapsed_seconds = NULL) {
    remaining <- request$deadline_uptime - helpers$kitClock()
    if (remaining <= 0) stop("KIT member ceiling exhausted")
    as.integer(max(1, floor(min(remaining, if (control) 120 else remaining) * 1000)))
}
result <- tryCatch(helpers$runKitMember(helpers, member_dir, request$application_root, args[[1L]], request$corpus,
    request$journey, request$fault, api), error = function(error) {
        record <- list(status = "driver_failed", error = conditionMessage(error), corpus = request$journey$corpus_key)
        saveRDS(record, file.path(member_dir, "journey-result.rds"), compress = FALSE, version = 3)
        record
    })
cat("MEMBER", result$status, "\n")
quit(save = "no", status = 0L)
