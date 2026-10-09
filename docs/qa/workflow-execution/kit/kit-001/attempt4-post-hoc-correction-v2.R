# KIT-001 attempt 4 correction, successor record (audit PA-01/PA-02). This is the complete correction
# as it ran on 2026-10-09 with the harness at the content later committed as dfd6aacd (the
# `journey$correction$harness_commit` field in the member says 3512f84d because the attach-loop fix
# was applied before the run and committed after). The original fingerprint script is kept unchanged.
# Not re-runnable against the recorded member: both stages refuse when their outputs already exist.
# Usage: Rscript attempt4-post-hoc-correction-v2.R <kit_001 library> <member dir> <repo> <stage dir>
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
library <- args[[1L]]; member <- args[[2L]]; repo <- args[[3L]]; stage <- args[[4L]]
.libPaths(c(library, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
env <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), env)
`%||%` <- function(x, y) if (is.null(x)) y else x
# Stage 1: fingerprints from the immutable store, attached to the driver's in-run snapshots.
journey <- readRDS(file.path(member, "journey-result.rds"))
stopifnot(!file.exists(file.path(member, "journey-result.driver.rds")), identical(journey$status, "complete"))
file.copy(file.path(member, "journey-result.rds"), file.path(member, "journey-result.driver.rds"))
file.copy(file.path(member, "operation-receipt.rds"), file.path(member, "operation-receipt.driver.rds"))
fingerprints <- env$kitFingerprintStore(journey$store, api)
head <- env$kitReadHead(journey$store, api)
journey$final_head_sha256 <- fingerprints$head_sha256
journey$product_fingerprints <- fingerprints$products
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
    reason = "driver fingerprinting overran the C stack on the RUV diagnostic product; recomputed from the immutable store with the corrected harness",
    harness_commit = "dfd6aacd (content as run; see outcome.md)", when = format(Sys.time()))
saveRDS(journey, file.path(member, "journey-result.rds"), compress = FALSE, version = 3)
operation <- readRDS(file.path(member, "operation-receipt.rds"))
operation$populations_complete <- journey$populations_complete
operation$corrections <- list(journey$correction)
saveRDS(operation, file.path(member, "operation-receipt.rds"))
# Stage 2: independent verifier from the corrected harness (outputs beside the member).
verify_dir <- file.path(member, "verify2")
dir.create(verify_dir)
stopifnot(!file.exists(file.path(member, "verification.rds")))
verification <- env$verifyKitMember(verify_dir, repo, "kit:lifecycle:complete", api)
stopifnot(identical(verification$status, "complete"))
# Stage 3: decision on the amended receipt.
receipt <- readRDS(file.path(stage, "stage.rds"))
receipt$cases$member$operation <- operation
receipt$cases$member$verifier_completion <- list(accepted = TRUE, reason = NULL,
    corrected = "verifier rerun from the corrected harness after the post-hoc fingerprints")
bindings <- env$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
cases <- env$kitCases(bindings, "kit-lifecycle-complete")
file.copy(file.path(stage, "decision.rds"), file.path(stage, "decision.driver.rds"))
decision <- env$kitLifecycleDecision(receipt, cases, bindings)
receipt$decision <- decision
receipt$corrections <- operation$corrections
saveRDS(decision, file.path(stage, "decision.rds"))
saveRDS(receipt, file.path(stage, "stage.rds"))
cat("DECISION", decision$status, decision$failing_rule %||% "", "\n")
