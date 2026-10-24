# KIT-057: tracked summary of a kit-journey-057 pair stage (both members' measurements and verification, the pair
# decision) from the stage receipt. Usage: Rscript docs/qa/workflow-execution/kit/kit-057/stage-summary.R <stage dir> <label>
args <- commandArgs(TRUE)
stage <- args[[1L]]
label <- if (length(args) > 1L) args[[2L]] else "kit-journey-057"
decision <- readRDS(file.path(stage, "decision.rds"))
receipt <- readRDS(file.path(stage, "stage.rds"))
member_summary <- function(id) {
    item <- receipt$cases[[id]]
    op <- item$operation
    verification <- tryCatch(readRDS(file.path(stage, id, "verification.rds")), error = function(e) NULL)
    m <- decision$members[[id]]
    list(id = id, role = m$role, corpus = op$member$corpus, build = op$config$build, status = op$status,
        eligible = m$eligible, failing_rule = m$failing_rule, outcome = m$outcome,
        elapsed_seconds = m$elapsed_seconds, peak_bytes = m$peak_bytes, peak_anon_bytes = m$peak_anon_bytes,
        resident_end_bytes = m$resident_end_bytes, phase_seconds = m$phase_seconds, phase_peak_anon_bytes = m$phase_peak_anon_bytes,
        phase_status = m$phase_status, max_sample_interval = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms,
        swap_bytes = op$swap_bytes, slice_stopped = op$cleanup$slice_stopped, populations_complete = op$populations_complete,
        verification = list(accepted = item$verifier_completion$accepted, status = verification$status, oracle = verification$oracle,
            populations = verification$populations, oracle_mismatches = verification$oracle_mismatches,
            canonical_form_differences = verification$canonical_form_differences))
}
summary <- list(schema = "kit.pair_stage_summary.v1", ticket = "KIT-057", stage = label, date = format(Sys.Date()),
    decision = list(status = decision$status, rule = decision$rule, failing_rule = decision$failing_rule,
        control_outcome = decision$control_outcome, steps_compared = decision$steps_compared, steps_guarded = decision$steps_guarded,
        step_ratios = decision$step_ratios, whole_slice_anon_ratio = decision$whole_slice_anon_ratio,
        resident_end_ratio = decision$resident_end_ratio, promotion_authorized = decision$promotion_authorized),
    members = lapply(stats::setNames(nm = intersect(c("control-1", "candidate-1"), names(receipt$cases))), member_summary))
jsonlite::write_json(summary, file.path("docs/qa/workflow-execution/kit/kit-057", paste0(label, "-stage.json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("written", label, decision$status, "\n")
