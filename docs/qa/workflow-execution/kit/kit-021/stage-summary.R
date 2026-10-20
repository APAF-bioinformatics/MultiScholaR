# KIT-021: tracked summary of a kit-journey-021 or kit-journey-021t stage (the single journey member's decision, measurements and
# verification) from the stage receipt. Usage: Rscript docs/qa/workflow-execution/kit/kit-021/stage-summary.R <stage dir> <label>
args <- commandArgs(TRUE)
stage <- args[[1L]]
label <- if (length(args) > 1L) args[[2L]] else "kit-journey-021"
decision <- readRDS(file.path(stage, "decision.rds"))
receipt <- readRDS(file.path(stage, "stage.rds"))
item <- receipt$cases[["candidate-1"]]
op <- item$operation
verification <- readRDS(file.path(stage, "candidate-1", "verification.rds"))
m <- decision$member
summary <- list(schema = "kit.stage_summary.v1", ticket = "KIT-021", stage = label, date = format(Sys.Date()),
    decision = list(status = decision$status, rule = decision$rule, failing_rule = decision$failing_rule,
        steps = decision$steps, stop_after = decision$stop_after, promotion_authorized = decision$promotion_authorized),
    member = list(id = "candidate-1", corpus = op$member$corpus, build = op$config$build, status = op$status, exact = m$exact,
        elapsed_seconds = m$elapsed_seconds, peak_bytes = m$peak_bytes, peak_anon_bytes = m$peak_anon_bytes,
        resident_end_bytes = m$resident_end_bytes, phase_seconds = m$phase_seconds, phase_peak_anon_bytes = m$phase_peak_anon_bytes,
        phase_status = m$phase_status, max_sample_interval = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms,
        swap_bytes = op$swap_bytes, slice_stopped = op$cleanup$slice_stopped, populations_complete = op$populations_complete),
    verification = list(accepted = item$verifier_completion$accepted, status = verification$status, oracle = verification$oracle,
        populations = verification$populations, oracle_mismatches = verification$oracle_mismatches,
        canonical_form_differences = verification$canonical_form_differences, files = verification$files))
jsonlite::write_json(summary, file.path("docs/qa/workflow-execution/kit/kit-021", paste0(label, "-stage.json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("written", label, decision$status, "\n")
