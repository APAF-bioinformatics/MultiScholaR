# KIT-015: tracked summary of stage kit-015 (decision and per-member measurements) from the stage receipt.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-015/stage-summary.R <proof>/kit/runtime/kit-015
stage <- commandArgs(TRUE)[[1L]]
decision <- readRDS(file.path(stage, "decision.rds"))
receipt <- readRDS(file.path(stage, "stage.rds"))
members <- lapply(names(decision$members), function(id) {
    m <- decision$members[[id]]
    op <- receipt$cases[[id]]$operation
    list(id = id, operation = m$operation, variant = m$variant, population = m$population, oracle = m$oracle,
        status = op$status, exact = m$exact, failing_rule = m$failing_rule, elapsed_seconds = m$elapsed_seconds,
        peak_bytes = m$peak_bytes, peak_anon_bytes = m$peak_anon_bytes, batches = m$batches,
        max_sample_interval = op$max_sample_interval, swap_bytes = op$swap_bytes, slice_stopped = op$cleanup$slice_stopped)
})
summary <- list(schema = "kit.stage_summary.v1", ticket = "KIT-015", stage = "kit-015", date = format(Sys.Date()),
    decision = list(status = decision$status, rule = decision$rule, failing_rule = decision$failing_rule,
        promotion_authorized = decision$promotion_authorized), members = members)
jsonlite::write_json(summary, file.path("docs/qa/workflow-execution/kit/kit-015/kit-015-stage.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("written", decision$status, "\n")
