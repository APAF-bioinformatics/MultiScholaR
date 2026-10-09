# KIT-001 attempt 4 correction (2026-10-09): the driver's fingerprinting overran the C stack on the
# RUV diagnostic product (kitMaskPaths recursion; fixed in 99487853). The in-run per-step snapshots
# (head sha256, state, entity id, rows) are the driver's; the fingerprints are recomputed here from
# the immutable store with the corrected harness. The driver's originals are kept beside the files.
# Usage: Rscript attempt4-post-hoc-fingerprints.R <kit_001 library> <member dir>
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
.libPaths(c(args[[1L]], .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
env <- new.env(parent = globalenv())
sys.source("tools/profiling/workflow_kit.R", env)
member <- args[[2L]]
journey <- readRDS(file.path(member, "journey-result.rds"))
stopifnot(!file.exists(file.path(member, "journey-result.driver.rds")), identical(journey$status, "complete"))
file.copy(file.path(member, "journey-result.rds"), file.path(member, "journey-result.driver.rds"))
file.copy(file.path(member, "operation-receipt.rds"), file.path(member, "operation-receipt.driver.rds"))
fingerprints <- env$kitFingerprintStore(journey$store, api)
head <- env$kitReadHead(journey$store, api)
journey$final_head_sha256 <- fingerprints$head_sha256
journey$product_fingerprints <- fingerprints$products
`%||%` <- function(x, y) if (is.null(x)) y else x
for (id in names(journey$populations)) {
    state <- journey$populations[[id]]$state %||% journey$populations[[id]]$expected_state
    if (is.null(state) || is.null(fingerprints$states[[state]])) next
    journey$populations[[id]]$state_as_expected <- identical(state, journey$populations[[id]]$expected_state)
    if (is.null(journey$populations[[id]]$state)) {
        journey$populations[[id]]$state <- state
        journey$populations[[id]]$store <- journey$store
        journey$populations[[id]]$entity_id <- head$states[[state]]
        journey$populations[[id]]$rows <- head$rows[[state]]
        journey$populations[[id]]$snapshot_time <- "after_adoption"
    } else journey$populations[[id]]$snapshot_time <- "after_step"
    journey$populations[[id]]$fingerprint <- fingerprints$states[[state]]
}
steps <- env$kitJourneySteps()
journey$populations_complete <- all(vapply(steps, function(step) is.null(step$state) ||
    !is.null(journey$populations[[step$id]]$fingerprint), logical(1)))
journey$fingerprint_error <- NULL
journey$correction <- list(kind = "post_hoc_fingerprints",
    reason = "driver fingerprinting overran the C stack on the RUV diagnostic product; recomputed from the immutable store with the corrected harness (99487853)",
    harness_commit = system("git rev-parse HEAD", intern = TRUE), when = format(Sys.time()))
saveRDS(journey, file.path(member, "journey-result.rds"), compress = FALSE, version = 3)
operation <- readRDS(file.path(member, "operation-receipt.rds"))
operation$populations_complete <- journey$populations_complete
operation$corrections <- list(journey$correction)
saveRDS(operation, file.path(member, "operation-receipt.rds"))
in_run <- sum(vapply(journey$populations, function(p) identical(p$snapshot_time, "after_step"), logical(1)))
cat("populations_complete:", journey$populations_complete, "; in-run snapshots:", in_run,
    "; final head equals the DA-step snapshot:", identical(journey$populations$da$head_sha256, fingerprints$head_sha256), "\n")
