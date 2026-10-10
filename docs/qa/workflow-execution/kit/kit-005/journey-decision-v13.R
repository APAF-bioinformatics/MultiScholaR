# KIT-005: the kit-journey attempt 6 stage decision recomputed under protocol v1.3 (user decision 2026-10-11:
# a journey member tolerates isolated sampling gaps up to 500 ms; control-1 had one 280 ms gap in 66,823
# samples). No member is rerun; the receipts gain the count of gaps above 250 ms computed from their
# memory traces, the v1.2 decision is kept beside the new one.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-005/journey-decision-v13.R <library>
library_path <- commandArgs(TRUE)[[1L]]
.libPaths(c(library_path, .libPaths()))
loadNamespace("MultiScholaR")
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
stage_root <- file.path(bindings$root, "kit-journey")
snapshot <- file.path(stage_root, "source")
receipt_path <- file.path(stage_root, "stage.rds")
receipt <- readRDS(receipt_path)
stopifnot(identical(receipt$status, "complete"))
for (id in names(receipt$cases)) {
    trace <- utils::read.delim(file.path(stage_root, id, "memory-trace.tsv"), stringsAsFactors = FALSE)
    uptime <- suppressWarnings(as.numeric(trace$uptime))
    uptime <- uptime[is.finite(uptime)]
    gaps <- diff(uptime)
    receipt$cases[[id]]$operation$gaps_above_250ms <- sum(gaps > 0.25)
    receipt$cases[[id]]$operation$max_sample_interval <- max(gaps)
    operation_path <- file.path(stage_root, id, "operation-receipt.rds")
    operation <- readRDS(operation_path)
    operation$gaps_above_250ms <- receipt$cases[[id]]$operation$gaps_above_250ms
    saveRDS(operation, operation_path)
    cat(sprintf("%s: max gap %.3f s, gaps above 250 ms: %d\n", id, max(gaps), sum(gaps > 0.25)))
}
cases <- helpers$kitCases(helpers$kitBindings(file.path(snapshot, "docs/qa/workflow-execution/kit/bindings.json")), "kit-journey")
receipt$decision_protocol_v12 <- receipt$decision
receipt$decision <- helpers$kitJourneyDecision(receipt, cases, bindings)
receipt$decision$protocol <- list(version = "1.3", sha256 = bindings$protocol_sha256,
    note = "recomputed from the attempt 6 receipts under v1.3 (journey gap tolerance 500 ms); v1.2 decision kept as decision_protocol_v12")
saveRDS(receipt, receipt_path)
saveRDS(receipt$decision, file.path(stage_root, "decision.rds"))
cat("DECISION", receipt$decision$status, paste(receipt$decision$failing_rule, collapse = ","), "\n")
str(receipt$decision$noise)
str(receipt$decision$production_default)
