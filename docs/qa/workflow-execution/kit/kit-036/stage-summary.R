# KIT-036: tracked summary of stage kit-journey-036 (the closure pair's decision, both members' measurements,
# verifications and corrections) from the stage receipt. Usage: Rscript docs/qa/workflow-execution/kit/kit-036/stage-summary.R <stage dir> <label>
args <- commandArgs(TRUE)
stage <- args[[1L]]
label <- if (length(args) > 1L) args[[2L]] else "kit-journey-036"
decision <- readRDS(file.path(stage, "decision.rds"))
receipt <- readRDS(file.path(stage, "stage.rds"))
member <- function(id) {
    item <- receipt$cases[[id]]
    op <- item$operation
    m <- decision$members[[id]]
    jpath <- file.path(stage, id, "journey-result.rds")
    journey <- if (file.exists(jpath)) readRDS(jpath) else NULL
    vpath <- file.path(stage, id, "verification.rds")
    verification <- if (file.exists(vpath)) readRDS(vpath) else NULL
    list(id = id, build = op$config$build, library = op$member$library, corpus = op$member$corpus, status = op$status,
        outcome = m$outcome, eligible = m$eligible, failing_rule = m$failing_rule, failed_step = journey$failed_step, error = journey$error,
        elapsed_seconds = m$elapsed_seconds, peak_bytes = m$peak_bytes, peak_anon_bytes = m$peak_anon_bytes,
        resident_end_bytes = m$resident_end_bytes, phase_seconds = m$phase_seconds, phase_peak_anon_bytes = m$phase_peak_anon_bytes,
        phase_status = m$phase_status, max_sample_interval = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms,
        swap_bytes = op$swap_bytes, slice_stopped = op$cleanup$slice_stopped, populations_complete = op$populations_complete,
        verification = list(accepted = item$verifier_completion$accepted, correction = item$verifier_completion$correction,
            status = verification$status, oracle_mismatches = verification$oracle_mismatches,
            additional_products = verification$additional_products, operation_identity_products = verification$operation_identity_products,
            canonical_form_differences = verification$canonical_form_differences))
}
summary <- list(schema = "kit.stage_summary.v1", ticket = "KIT-036", stage = label, date = format(Sys.Date()),
    decision = decision[setdiff(names(decision), "members")], corrections = receipt$corrections,
    decisions_before_correction = lapply(receipt$decisions_before_correction, function(d) d[c("status", "failing_rule")]),
    members = list(control = member("control-1"), candidate = member("candidate-1")))
jsonlite::write_json(summary, file.path("docs/qa/workflow-execution/kit/kit-036", paste0(label, "-stage.json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("written", label, decision$status, "\n")
