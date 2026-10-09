# KIT-005 correction 1 (2026-10-10): the kit-lifecycle-operation-complete member measured validly
# (status complete, 38 s, peak 1.53 GiB, zero swap, every P14 phase complete, lifecycle tests passed),
# but the verifier compared the output pages against the oracle's page *count* key (`screen_pages`)
# instead of the page hashes (`screen_page_hashes`). The member is not rerun; the verifier is rerun
# from the corrected harness (`657e6ea7`) against the stage's frozen bindings, the failed verification is
# kept beside it, and the decision is recomputed with the stage's own decision rule.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-005/operation-complete-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage <- "kit-lifecycle-operation-complete"
stage_root <- file.path(bindings$root, stage)
snapshot <- file.path(stage_root, "source")
owned <- file.path(stage_root, "member")
failed <- file.path(owned, "verification.rds")
stopifnot(file.exists(failed), identical(readRDS(failed)$status, "failed"))
stopifnot(file.rename(failed, file.path(owned, "verification.harness-failed.rds")))
verify_dir <- file.path(owned, "verify-correction-1")
dir.create(verify_dir)
# Bindings from the stage's frozen snapshot (the member's binding digest); harness from the corrected repo.
verification <- helpers$verifyKitMember(verify_dir, snapshot, "kit:lifecycle:operation:complete", api)
stopifnot(identical(verification$status, "complete"))
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
receipt$cases$member$verifier_completion <- list(accepted = TRUE, reproduced = NULL, reason = NULL,
    correction = list(id = 1L, harness_commit = system("git rev-parse --short HEAD", intern = TRUE),
        harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
        reason = "verifier compared page hashes against the oracle page count key; member measurement untouched"))
receipt$decision_verifier_failed <- receipt$decision
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), stage)
receipt$decision <- helpers$kitOperationLifecycleDecision(receipt, cases, bindings)
receipt$corrections <- c(receipt$corrections, list(receipt$cases$member$verifier_completion$correction))
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
cat("DECISION", receipt$decision$status, receipt$decision$failing_rule %||% "", "\n")
