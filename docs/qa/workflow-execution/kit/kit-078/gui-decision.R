# KIT-078: the `kit-gui-0NN` decision recomputed from a finished member, with the harness's current verifier, when the
# coordinator could not take it (the baseline's reference store had been renamed; the member itself is untouched).
# Usage (kit library on R_LIBS): Rscript docs/qa/workflow-execution/kit/kit-078/gui-decision.R <stage dir>
stage_root <- commandArgs(TRUE)[[1L]]
`%||%` <- function(x, y) if (is.null(x)) y else x
ns <- loadNamespace("MultiScholaR")
api <- function(name) get(name, envir = ns)
helpers <- new.env()
sys.source("tools/profiling/workflow_kit.R", envir = helpers)
bindings <- jsonlite::fromJSON("docs/qa/workflow-execution/kit/bindings.json", simplifyVector = FALSE)
receipt <- readRDS(file.path(stage_root, "stage.rds"))
stage <- basename(stage_root)
cases <- helpers$kitCases(bindings, stage)
decision <- helpers$kitGuiMemberDecision(receipt, cases, bindings, member_dir = file.path(stage_root, cases[[1L]]$id), api = api)
decision$note <- "decision recomputed from the finished member after the harness correction"
saveRDS(decision, file.path(stage_root, "decision.rds"))
receipt$decision <- decision
receipt$status <- "complete"
receipt$error <- NULL
saveRDS(receipt, file.path(stage_root, "stage.rds"))
cat("DECISION", decision$status, decision$failing_rule %||% "", "\n")
cat("elapsed", decision$elapsed_seconds, "s; anon", round(decision$peak_anon_bytes / 1024^3, 2), "GiB\n")
cat("dispatched:\n"); for (n in names(decision$dispatched)) cat(" ", n, decision$dispatched[[n]]$requests, "requests,", decision$dispatched[[n]]$complete, "complete\n")
cat("reproduction:", decision$reproduction$fingerprint_mismatches, "mismatches\n")
for (m in utils::head(decision$reproduction$mismatches, 12)) cat("  ", paste(format(unlist(m)), collapse = " | "), "\n")
