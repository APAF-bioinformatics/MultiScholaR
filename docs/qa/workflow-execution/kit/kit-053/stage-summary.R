# KIT-053: tracked summary of a kit-journey-053 stage (the single journey member's decision, measurements and
# verification) from the stage receipt. Usage: Rscript docs/qa/workflow-execution/kit/kit-053/stage-summary.R <stage dir> <label>
# KIT-068: a lifecycle stage (one case `member`, the lifecycle decision: status, failing rule, route, the in-memory
# and head-rule notes) is summarised by the same script; the ticket directory and id come from `KIT_LIMPA_TICKET_DIR`.
args <- commandArgs(TRUE)
stage <- args[[1L]]
label <- if (length(args) > 1L) args[[2L]] else "kit-journey-053"
# KIT-056: another ticket's single journey member writes under its own directory and ticket id (`KIT_LIMPA_TICKET_DIR`).
ticket_dir <- Sys.getenv("KIT_LIMPA_TICKET_DIR", "kit-053")
ticket_id <- paste0("KIT-", sub("^kit-", "", ticket_dir))
`%||%` <- function(x, y) if (is.null(x)) y else x
decision <- readRDS(file.path(stage, "decision.rds"))
receipt <- readRDS(file.path(stage, "stage.rds"))
case_id <- if (!is.null(receipt$cases[["candidate-1"]])) "candidate-1" else names(receipt$cases)[[1L]]
item <- receipt$cases[[case_id]]
op <- item$operation
verification_path <- file.path(stage, case_id, "verification.rds")
verification <- if (file.exists(verification_path)) readRDS(verification_path) else NULL
m <- decision$member %||% list()
steps <- names(op$phase_status)
summary <- list(schema = "kit.stage_summary.v1", ticket = ticket_id, stage = label, date = format(Sys.Date()),
    decision = list(status = decision$status, rule = decision$rule, failing_rule = decision$failing_rule,
        steps = decision$steps %||% steps, stop_after = decision$stop_after, promotion_authorized = decision$promotion_authorized,
        expect_failure = decision$expect_failure, route = decision$route %||% op$member$route, in_memory = decision$in_memory,
        head_rule = decision$head_rule),
    member = list(id = case_id, corpus = op$member$corpus, build = op$config$build, status = op$status, exact = m$exact,
        elapsed_seconds = m$elapsed_seconds %||% op$elapsed_seconds, peak_bytes = m$peak_bytes %||% op$peak_bytes,
        peak_anon_bytes = m$peak_anon_bytes %||% op$sampled_anon_bytes, resident_end_bytes = m$resident_end_bytes %||% op$resident_end_bytes,
        phase_seconds = m$phase_seconds %||% op$phase_seconds, phase_peak_anon_bytes = m$phase_peak_anon_bytes %||% op$phase_peak_anon_bytes,
        phase_status = m$phase_status %||% op$phase_status, max_sample_interval = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms,
        swap_bytes = op$swap_bytes, slice_stopped = op$cleanup$slice_stopped, populations_complete = op$populations_complete,
        guard = op$guard, failed_step = item$result_status, error = item$error),
    verification = list(accepted = item$verifier_completion$accepted, status = verification$status, oracle = verification$oracle,
        populations = verification$populations, oracle_mismatches = verification$oracle_mismatches,
        canonical_form_differences = verification$canonical_form_differences, files = verification$files, failure = verification$failure))
jsonlite::write_json(summary, file.path("docs/qa/workflow-execution/kit", ticket_dir, paste0(label, "-stage.json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("written", label, decision$status, "\n")
