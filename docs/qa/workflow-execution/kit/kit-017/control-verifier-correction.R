# KIT-017 correction 1 (2026-10-15, verifier only): control-1 (`kit_018`) completed the whole 300-sample journey and
# every state matched the KIT-005 oracles, but the verifier refused on `da/products/names`: the DA state carries the
# interactive volcano page as a rendered asset (`rendered_asset:<id>`, KIT-016), a product the kit_001 oracle does
# not have, and the name comparison required the two sets to be equal. The comparison now reports a product the
# oracle lacks as an additional product and still judges every oracle product (a missing one mismatches). The member
# is not rerun: its failed verification is kept beside the rerun, the verifier runs again from the corrected harness
# against the stage's frozen bindings, and the decision is recomputed over both members.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-017/control-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage_root <- file.path(bindings$root, "kit-journey-017")
snapshot <- file.path(stage_root, "source")
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"), identical(receipt$decision$status, "invalid"))
owned <- file.path(stage_root, "control-1")
failed <- file.path(owned, "verification.rds")
kept <- file.path(owned, "verification.harness-failed.rds")
run_dir <- file.path(owned, "verify-correction-1")
if (file.exists(failed)) {
    first <- readRDS(failed)
    stopifnot(identical(first$status, "failed"))
    if (identical(first$run_dir, run_dir)) unlink(failed) else stopifnot(
        identical(first$oracle_mismatches, "da/products/names"), !file.exists(kept), file.rename(failed, kept))
}
stopifnot(file.exists(kept))
dir.create(run_dir, showWarnings = FALSE)
verification <- helpers$verifyKitMember(run_dir, snapshot, "kit:journey-017:control:1", api)
stopifnot(identical(verification$status, "complete"), !length(verification$oracle_mismatches))
cat("additional products reported:", paste(verification$additional_products, collapse = ", "), "\n")
cat("canonical-form differences:", paste(verification$canonical_form_differences, collapse = ", "), "\n")
correction <- list(id = 1L, member = "control-1",
    harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
    reason = "the DA state's rendered volcano page (KIT-016) is a product the kit_001 oracle does not carry; the comparison reports products the oracle lacks and judges every oracle product; verifier rerun only, member measurement kept",
    additional_products = verification$additional_products)
receipt$cases[["control-1"]]$verifier_completion <- list(accepted = TRUE, reproduced = NULL, reason = NULL, correction = correction)
receipt$cases[["control-1"]]$verifier_exit <- 0L
receipt$corrections <- c(receipt$corrections, list(correction))
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), "kit-journey-017")
receipt$decisions_before_correction <- c(receipt$decisions_before_correction, list(receipt$decision))
receipt$decision <- helpers$kitJourneyTicketDecision(receipt, cases, bindings)
receipt$decision$protocol <- list(version = "1.4", sha256 = bindings$protocol_sha256,
    note = "recomputed after correction 1 (verifier only, control-1); the launch decision kept in decisions_before_correction")
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
d <- receipt$decision
cat("DECISION", d$status, paste(names(d$failing_rule), d$failing_rule, collapse = ","), "| control", d$control_outcome %||% "", "\n")
