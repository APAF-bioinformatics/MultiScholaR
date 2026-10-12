# KIT-008 correction 1 (2026-10-12, verifier only): the gui member completed the screen journey in 303 s on the
# kit_008 build with both production dispatches on the managed worker, but the verifier never ran: the selector
# `kit:gui-008:candidate:1` was mapped to case id `candidate-1` while the stage's one case is `member`, so
# `verifyKitMember` aborted before writing a record (launch log kept: member/verification-launch.log). The member
# is not rerun: the verifier runs again from the corrected harness against the stage's frozen bindings and the
# decision is recomputed; the launch decision is kept in `decisions_before_correction`.
# Correction 2 (same day, user decision "1"): with every table, support entry and product equal to the reference, the
# recomputed decision was `rejected` on the sample-filter state's summary alone: the physical declaration records
# `core_utilisation = NA` where the native loop recorded its multidplyr worker spec and copied the support references
# into the args, and the two word the step description differently; every later state's summary carried that args
# section forward. The canonical form now drops those three bookkeeping fields and the decision recomputes both
# members' fingerprints from their immutable stores under it. Rerunning this script applies correction 2.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-008/member-verifier-correction.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, ns)
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage_root <- file.path(bindings$root, "kit-gui-008")
snapshot <- file.path(stage_root, "source")
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"), receipt$decision$status %in% c("invalid", "rejected"),
    unname(receipt$decision$failing_rule) %in% c("verifier_accepted", "reference_mismatch"))
owned <- file.path(stage_root, "member")
stopifnot(any(grepl("subscript out of bounds", readLines(file.path(owned, "verification-launch.log"), warn = FALSE), fixed = TRUE)))
run_dir <- file.path(owned, "verify-correction-1")
existing <- file.path(owned, "verification.rds")
# The verifier itself ran once in this correction; the decision recomputation then failed on a second harness
# defect (the runtime root field). That verification is reused, never a launch one (none was written).
verification <- if (file.exists(existing)) readRDS(existing) else NULL
if (is.null(verification) || !identical(verification$run_dir, run_dir) || !identical(verification$status, "complete")) {
    if (file.exists(existing)) unlink(existing)
    dir.create(run_dir, showWarnings = FALSE)
    verification <- helpers$verifyKitMember(run_dir, snapshot, "kit:gui-008:candidate:1", api)
}
stopifnot(identical(verification$status, "complete"))
correction_id <- length(receipt$corrections) + 1L
correction <- list(id = correction_id, harness_sha256 = digest::digest(file = file.path(repo, "tools/profiling/workflow_kit.R"), algo = "sha256"),
    reason = "kitSelectorStage mapped the gui selector to case id candidate-1 while the stage's case is member, so the verifier aborted before running; the decision read the reference member under a bindings field the real bindings do not carry; verifier rerun and decision recomputed only, member measurement untouched")
if (correction_id == 2L) correction$reason <- "canonical form drops core_utilisation, the args-embedded support copies and the step description (execution bookkeeping, not science); both members' fingerprints recomputed from their stores; decision recomputed only"
receipt$cases[["member"]]$verifier_completion <- list(accepted = TRUE, reproduced = verification$reproduced, reason = NULL, correction = correction)
receipt$cases[["member"]]$verifier_exit <- 0L
receipt$corrections <- c(receipt$corrections, list(correction))
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), "kit-gui-008")
receipt$decisions_before_correction <- c(receipt$decisions_before_correction, list(receipt$decision))
receipt$decision <- helpers$kitGuiMemberDecision(receipt, cases, bindings, member_dir = owned, api = api)
receipt$decision$protocol <- list(version = "1.3", sha256 = bindings$protocol_sha256,
    note = "recomputed after correction 1 (verifier only); the launch decision kept in decisions_before_correction")
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
d <- receipt$decision
cat("DECISION", d$status, paste(d$failing_rule, collapse = ","), "| elapsed", round(d$elapsed_seconds), "s | peak",
    round(d$peak_bytes / 1024^3, 3), "GiB | anon", round(d$peak_anon_bytes / 1024^3, 3), "GiB | resident end",
    round(d$resident_end_bytes / 1024^3, 3), "GiB | dispatched", paste(names(d$dispatched), vapply(d$dispatched, function(x) x$complete, 1L), collapse = ","),
    "| reference", d$reference, "| mismatches", d$reproduction$fingerprint_mismatches, "\n")
cat("verification: files", verification$files, "populations", length(verification$populations), "\n")
