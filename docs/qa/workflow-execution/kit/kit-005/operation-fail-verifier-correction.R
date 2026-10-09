# KIT-005 correction 2 (2026-10-10): the kit-lifecycle-operation-fail member measured validly (the worker
# failed inside the aggregate stage as injected; 28 s, zero swap, slice stopped) and its verifier passed
# every check, but (a) the verifier saved its record from the pre-return variable, so the completed
# verification was recorded as failed, and (b) the coordinator compared the state checkpoint id with the
# bound input component id, so `head_unchanged` was FALSE although the component was unchanged. The member
# is not rerun: the verifier is rerun from the corrected harness (`61e8e3a6`), the head check is recomputed
# from the member's immutable project store, the originals are kept, and the decision is recomputed.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-005/operation-fail-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage <- "kit-lifecycle-operation-fail"
stage_root <- file.path(bindings$root, stage)
snapshot <- file.path(stage_root, "source")
owned <- file.path(stage_root, "member")
failed <- file.path(owned, "verification.rds")
stopifnot(file.exists(failed), identical(readRDS(failed)$status, "failed"))
stopifnot(file.rename(failed, file.path(owned, "verification.harness-failed.rds")))
verify_dir <- file.path(owned, "verify-correction-2")
dir.create(verify_dir)
verification <- helpers$verifyKitMember(verify_dir, snapshot, "kit:lifecycle:operation:fail", api)
stopifnot(identical(verification$status, "complete"))
operation_path <- file.path(owned, "operation-receipt.rds")
stopifnot(file.copy(operation_path, file.path(owned, "operation-receipt.pre-correction-2.rds")))
operation <- readRDS(operation_path)
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), stage)
current <- helpers$kitCurrentInputComponentId(file.path(owned, "project/.workflow-entities"), api)
operation$head_unchanged_original <- operation$head_unchanged
operation$head_unchanged <- identical(current, cases[[1L]]$input$entity_id)
operation$correction <- list(id = 2L, field = "head_unchanged", recomputed_from = "member project store (immutable)",
    harness_commit = system("git rev-parse --short HEAD", intern = TRUE))
saveRDS(operation, operation_path)
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
receipt$cases$member$operation <- operation
receipt$cases$member$verifier_completion <- list(accepted = TRUE, reproduced = NULL, reason = NULL,
    correction = list(id = 2L, harness_commit = operation$correction$harness_commit,
        harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
        reason = "verification record saved from the pre-return variable; head check compared the wrong identity; member measurement untouched"))
receipt$decision_verifier_failed <- receipt$decision
receipt$decision <- helpers$kitOperationLifecycleDecision(receipt, cases, bindings)
receipt$corrections <- c(receipt$corrections, list(receipt$cases$member$verifier_completion$correction))
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
cat("head_unchanged", operation$head_unchanged, "DECISION", receipt$decision$status, receipt$decision$failing_rule %||% "", "\n")
