# KIT-070: compare a fresh route record with a stored oracle record step by step through the harness comparators
# (the stores behind the kit_065 oracles were overwritten by later recorder runs; the fingerprints in the RDS stand).
# usage: Rscript docs/qa/workflow-execution/kit/kit-070/compare-oracles.R <reference rds> <candidate rds>
args <- commandArgs(TRUE); reference <- readRDS(args[[1L]]); candidate <- readRDS(args[[2L]])
harness <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", harness)
ref <- reference$populations; cand <- candidate$populations
cat(sprintf("reference %s (%s, build %s): %d steps; candidate %s (build %s): %d steps\n", basename(args[[1L]]), reference$route %||% "?",
    substr(reference$build %||% "?", 1, 12), length(ref), basename(args[[2L]]), substr(candidate$build %||% "?", 1, 12), length(cand)))
mismatches <- harness$kitReproduceMismatches(ref, cand, ignore_products = harness$kitOperationIdentityProducts)
per_state <- unlist(lapply(names(ref), function(step) {
    if (is.null(cand[[step]])) return(paste0(step, "/missing"))
    m <- harness$kitOracleStateMismatches(cand[[step]]$fingerprint, ref[[step]]$fingerprint)
    if (length(m)) paste0(step, "/", m) else character()
}))
states_ok <- vapply(names(ref), function(step) identical(ref[[step]]$state, cand[[step]]$state), logical(1))
cat(sprintf("steps compared %d; state names identical %d; reproduce mismatches %d; per-state mismatches (class/rows/components/support/summary) %d\n",
    length(ref), sum(states_ok), length(mismatches), length(per_state)))
for (m in mismatches) cat("  MISMATCH", m, "\n"); for (m in per_state) cat("  STATE-MISMATCH", m, "\n")
additional <- attr(mismatches, "additional_products"); if (length(additional)) for (a in additional) cat("  additional product (reported, not judged):", a, "\n")
cat(if (!length(mismatches) && !length(per_state) && all(states_ok)) "RESULT identical\n" else "RESULT differs\n")
