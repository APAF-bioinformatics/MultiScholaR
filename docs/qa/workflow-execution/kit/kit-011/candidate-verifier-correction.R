# KIT-011 correction 1 (2026-10-11, verifier only): candidate-1 completed the whole 300-sample journey on the
# kit_011 build and every state matched the KIT-005 oracles, but the `normalization_science` product differed
# because its nested checkpoint records carry the KIT-002 `declaration` field (additive storage format) that the
# oracles, made on kit_001 before KIT-002, do not; every scientific field (best_k, controls, percentage,
# parameters, ruv_result) is identical. The canonical form now drops that field. The member is not rerun: its
# failed verification is kept beside the rerun, the verifier runs again from the corrected harness against the
# stage's frozen bindings, and the decision is recomputed over both members.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-011/candidate-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage_root <- file.path(bindings$root, "kit-journey-011")
snapshot <- file.path(stage_root, "source")
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"), identical(receipt$decision$status, "invalid"),
    identical(unname(receipt$decision$failing_rule), "verifier_accepted"))
owned <- file.path(stage_root, "candidate-1")
failed <- file.path(owned, "verification.rds")
kept <- file.path(owned, "verification.harness-failed.rds")
run_dir <- file.path(owned, "verify-correction-1")
if (file.exists(failed)) {
    first <- readRDS(failed)
    stopifnot(identical(first$status, "failed"))
    # The launch verification is kept; a failed attempt of this correction is discarded.
    if (identical(first$run_dir, run_dir)) unlink(failed) else stopifnot(
        identical(first$oracle_mismatches, "normalisation_ruv/products/normalization_science"), !file.exists(kept), file.rename(failed, kept))
}
stopifnot(file.exists(kept))
dir.create(run_dir, showWarnings = FALSE)
verification <- helpers$verifyKitMember(run_dir, snapshot, "kit:journey-011:candidate:1", api)
stopifnot(identical(verification$status, "complete"), !length(verification$oracle_mismatches))
cat("in-run fingerprints differing from the recomputed ones only by the canonical-form correction:",
    paste(verification$canonical_form_differences, collapse = ", "), "\n")
correction <- list(id = 1L, harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
    reason = "canonical fingerprint dropped the KIT-002 record field `declaration`, absent from the kit_001 oracles and present in every later build's nested records; verifier rerun only, member measurement untouched")
receipt$cases[["candidate-1"]]$verifier_completion <- list(accepted = TRUE, reproduced = NULL, reason = NULL, correction = correction)
receipt$cases[["candidate-1"]]$verifier_exit <- 0L
receipt$corrections <- c(receipt$corrections, list(correction))
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), "kit-journey-011")
receipt$decisions_before_correction <- c(receipt$decisions_before_correction, list(receipt$decision))
receipt$decision <- helpers$kitJourneyTicketDecision(receipt, cases, bindings)
receipt$decision$protocol <- list(version = "1.3", sha256 = bindings$protocol_sha256,
    note = "recomputed after correction 1 (verifier only); the launch decision kept in decisions_before_correction")
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
d <- receipt$decision
cat("DECISION", d$status, paste(d$failing_rule, collapse = ","), "| control", d$control_outcome, "lower bound", d$control_peak_is_lower_bound,
    "| whole-slice anon ratio", round(d$whole_slice_anon_ratio, 3), "| guarded steps", paste(d$steps_guarded, collapse = ","), "\n")
print(round(unlist(d$step_ratios), 3))
cat("identity products reported:", length(verification$operation_identity_products), "\n")
