# Tracked summary of the kit-journey stage (schema kit.journey_stage_summary.v1) regenerated from the stage
# receipt, so the tracked record carries the stage's current decision (KIT-005 audit PH-03).
# Usage: Rscript docs/qa/workflow-execution/kit/kit-005/journey-stage-summary.R [<output json>]
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
out <- commandArgs(TRUE)
out <- if (length(out)) out[[1L]] else file.path(repo, "docs/qa/workflow-execution/kit/kit-005/journey-attempt6.json")
receipt <- readRDS(file.path(bindings$root, "kit-journey", "stage.rds"))
gib <- function(x) if (is.numeric(x) && length(x) == 1L && is.finite(x)) round(x / 1024^3, 3) else NULL
member <- function(id) {
    op <- receipt$cases[[id]]$operation
    list(id = id, status = op$status, elapsed_seconds = round(op$elapsed_seconds), peak_gib = gib(op$peak_bytes),
        anon_peak_gib = gib(op$sampled_anon_bytes), resident_end_gib = gib(op$resident_end_bytes), swap_bytes = op$swap_bytes,
        max_sample_interval_s = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms, samples = op$samples,
        populations_complete = op$populations_complete, verifier_accepted = isTRUE(receipt$cases[[id]]$verifier_completion$accepted),
        reproduced = receipt$cases[[id]]$verifier_completion$reproduced, outer_guard = op$outer_guard,
        phase_seconds = lapply(op$phase_seconds, round), phase_peak_anon_gib = lapply(op$phase_peak_anon_bytes, gib),
        phase_status = op$phase_status)
}
decision <- receipt$decision
summary <- list(schema = "kit.journey_stage_summary.v1", stage = "kit-journey", attempt = 6L,
    control = sprintf("%s (%s, %s)", bindings$control, bindings$builds[[bindings$control]]$install_member,
        substr(bindings$builds[[bindings$control]]$source_commit, 1L, 8L)),
    decision = list(status = decision$status, protocol = decision$protocol, failing_rule = as.list(decision$failing_rule),
        noise = decision$noise, production_default = decision$production_default),
    decisions_superseded = list(protocol_v12 = receipt$decision_protocol_v12$status,
        v13_before_ph02 = receipt$decision_before_ph02$status),
    corrections = receipt$corrections, members = lapply(names(receipt$cases), member))
jsonlite::write_json(summary, out, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("summary written:", out, "decision", decision$status, "\n")
