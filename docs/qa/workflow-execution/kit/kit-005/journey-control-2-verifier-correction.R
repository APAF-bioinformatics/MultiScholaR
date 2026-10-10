# KIT-005 correction 7 (2026-10-10): the kit-journey attempt 6 member control-2 completed the whole
# 300-sample journey (3941 s; no sampling gap above 250 ms) and reproduced every scientific fingerprint of
# control-1 (0 mismatches over 18 steps), but the verifier's between-member file check compared raw
# paths, and 956 paths differ only by a timestamped directory, request-id staging folders or the
# session export's run id and export time. The member is not rerun: the verifier is rerun from the
# corrected harness (`6b598edf`) against the stage's frozen bindings, the failed verification is kept
# beside it, and the stage decision is recomputed with the stage's own decision rule.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-005/journey-control-2-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage <- "kit-journey"
stage_root <- file.path(bindings$root, stage)
snapshot <- file.path(stage_root, "source")
owned <- file.path(stage_root, "control-2")
failed <- file.path(owned, "verification.rds")
harness_failed <- file.path(owned, "verification.harness-failed.rds")
stopifnot(file.exists(harness_failed), identical(readRDS(harness_failed)$status, "failed"),
    length(readRDS(harness_failed)$reproduce_mismatches) == 0L, length(readRDS(harness_failed)$file_mismatches) > 0L)
# The failed verification already proves the within-member integrity: every in-run fingerprint was
# recomputed from the immutable store and matched, and every file hash matched the in-run record, before
# the between-member check compared raw paths. The corrected between-member check is applied here to
# both members' files re-hashed with the corrected harness (the files are immutable).
first <- readRDS(file.path(stage_root, "control-1", "journey-result.rds"))
second <- readRDS(file.path(owned, "journey-result.rds"))
reproduce_mismatches <- helpers$kitReproduceMismatches(first$populations, second$populations)
files_first <- helpers$kitFingerprintFiles(first$project)
files_second <- helpers$kitFingerprintFiles(second$project)
file_mismatches <- helpers$kitFileMismatches(files_first, files_second)
verification <- list(status = "complete", member_id = "control-2", selector = "kit:journey:control:2", kind = "journey",
    within_member = "verification.harness-failed.rds: recomputed fingerprints and file hashes matched the in-run record",
    reproduce_mismatches = reproduce_mismatches, file_mismatches = file_mismatches,
    reproduced = !length(reproduce_mismatches) && !length(file_mismatches),
    files = nrow(files_second), populations = sum(vapply(second$populations, function(p) !is.null(p$fingerprint), logical(1))),
    correction = 7L)
cat("reproduce mismatches:", length(reproduce_mismatches), "| file mismatches:", length(file_mismatches), "\n")
if (length(file_mismatches)) print(file_mismatches)
stopifnot(isTRUE(verification$reproduced))
saveRDS(verification, failed)
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"))
receipt$cases[["control-2"]]$verifier_completion <- list(accepted = TRUE, reproduced = TRUE, reason = NULL,
    correction = list(id = 7L, harness_commit = system("git rev-parse --short HEAD", intern = TRUE),
        harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
        reason = "between-member file check compared raw paths (timestamped directories, request ids, session export bookkeeping); member measurement untouched"))
receipt$decision_verifier_failed <- receipt$decision
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), stage)
receipt$decision <- helpers$kitJourneyDecision(receipt, cases, bindings)
receipt$corrections <- c(receipt$corrections, list(receipt$cases[["control-2"]]$verifier_completion$correction))
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
cat("DECISION", receipt$decision$status, paste(receipt$decision$failing_rule, collapse = ","), "\n")
str(receipt$decision$noise)
