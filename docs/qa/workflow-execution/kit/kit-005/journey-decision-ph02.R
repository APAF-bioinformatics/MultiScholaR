# KIT-005 audit remediation PH-02 (2026-10-11): the production-default rule of the stage decision follows the
# acceptance as amended by the user ("completes, or ends at a guard, with its per-step peaks recorded").
# default-1 ended at the memory guard in imputation; its verification (which stopped on `status == complete`)
# is kept as verification.harness-failed.rds and rerun from the corrected verifier; the stage decision is
# recomputed from the existing receipts (the v1.3 decision kept as decision_before_ph02). No member is rerun.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-005/journey-decision-ph02.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage_root <- file.path(bindings$root, "kit-journey")
# The members' binding digests were taken from the bindings frozen at launch (the stage's source snapshot);
# the verifier and the cases read those, the corrected harness comes from the repository.
snapshot <- file.path(stage_root, "source")
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"), identical(receipt$decision$status, "baseline_recorded"),
    isFALSE(receipt$decision$production_default$recorded))
owned <- file.path(stage_root, "default-1")
failed <- file.path(owned, "verification.rds")
kept <- file.path(owned, "verification.harness-failed.rds")
run_dir <- file.path(owned, "verify-ph02")
if (file.exists(failed)) {
    stopifnot(identical(readRDS(failed)$status, "failed"))
    if (identical(readRDS(failed)$run_dir, run_dir)) unlink(failed) else stopifnot(!file.exists(kept), file.rename(failed, kept))
}
dir.create(run_dir, showWarnings = FALSE)
verification <- helpers$verifyKitMember(run_dir, snapshot, "kit:journey:default:1", api)
stopifnot(identical(verification$status, "complete"), identical(verification$guard$outcome, "memory_guard"))
receipt$cases[["default-1"]]$verifier_completion <- list(accepted = TRUE, reproduced = NULL, reason = NULL,
    correction = list(id = 8L, reason = "KIT-005 audit PH-02: the production-default rule and verifier accept a member that ends at a protocol guard with its per-step peaks (acceptance as amended by the user 2026-10-11); member measurement untouched"))
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), "kit-journey")
receipt$decision_before_ph02 <- receipt$decision
receipt$decision <- helpers$kitJourneyDecision(receipt, cases, bindings)
receipt$decision$protocol <- list(version = "1.3", sha256 = bindings$protocol_sha256,
    note = "recomputed from the attempt 6 receipts after audit remediation PH-02 (production-default guard rule); the earlier v1.3 decision kept as decision_before_ph02")
receipt$corrections <- c(receipt$corrections, list(receipt$cases[["default-1"]]$verifier_completion$correction))
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
cat("DECISION", receipt$decision$status, "\n")
str(receipt$decision$production_default)
str(verification$guard)
