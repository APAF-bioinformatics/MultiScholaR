# KIT-006 correction 1 (2026-10-11): the sample-filter verifier required `execution$admission$batches`, a field
# only the adaptive admission writes; the control build's fixed-three policy records `count_batch_sizes` alone,
# so control-2 (complete, 30 s, exact pages) was refused by the verifier, not by any measurement. The member is
# not rerun: its failed verification is kept beside the rerun, the verifier runs again from the corrected harness
# against the stage's rerun bindings, and the decision is recomputed over all four members.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-006/control-2-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage_root <- file.path(bindings$root, "kit-006")
snapshot <- file.path(stage_root, "source.rerun1")
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"), identical(receipt$decision$status, "invalid"),
    identical(unname(receipt$decision$failing_rule), "verifier_accepted"))
owned <- file.path(stage_root, "control-2")
failed <- file.path(owned, "verification.rds")
kept <- file.path(owned, "verification.harness-failed.rds")
stopifnot(file.exists(failed), identical(readRDS(failed)$status, "failed"), !file.exists(kept), file.rename(failed, kept))
run_dir <- file.path(owned, "verify-correction-1")
dir.create(run_dir)
verification <- helpers$verifyKitMember(run_dir, snapshot, "kit:006:control:2", api)
stopifnot(identical(verification$status, "complete"), identical(verification$batches, 3L))
correction <- list(id = 1L, harness_commit = trimws(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE)),
    harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
    reason = "sample-filter verifier required admission$batches, absent from the control's fixed-three policy (count_batch_sizes only); verifier rerun only, member measurement untouched")
receipt$cases[["control-2"]]$verifier_completion <- list(accepted = TRUE, reproduced = NULL, reason = NULL, correction = correction)
receipt$cases[["control-2"]]$verifier_exit <- 0L
receipt$corrections <- c(receipt$corrections, list(correction))
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), "kit-006")
receipt$decisions_before_rerun <- c(receipt$decisions_before_rerun, list(receipt$decision))
receipt$decision <- helpers$kitOperationDecision(receipt, cases, bindings)
receipt$decision$protocol <- list(version = "1.3", sha256 = bindings$protocol_sha256,
    note = "recomputed after correction 1 (verifier only); attempt 1 and rerun 1 decisions kept in decisions_before_rerun")
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
cat("DECISION", receipt$decision$status, paste(receipt$decision$failing_rule, collapse = ","), "\n")
str(verification$filter); cat("batches", verification$batches, "\n")
str(lapply(receipt$decision$pairs, `[`, c("peak_ratio", "peak_anon_ratio", "not_worse", "batches")))
