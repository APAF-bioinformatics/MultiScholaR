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
stopifnot(identical(request$schema, 1L))
result <- tryCatch(helpers$runKitMember(helpers, member_dir, request$application_root, args[[1L]], request$corpus,
    request$journey, request$fault, api), error = function(error) {
        record <- list(status = "driver_failed", error = conditionMessage(error), corpus = request$journey$corpus_key)
        saveRDS(record, file.path(member_dir, "journey-result.rds"), compress = FALSE, version = 3)
        record
    })
cat("MEMBER", result$status, "\n")
quit(save = "no", status = 0L)
