# KIT-009: tracked summary of stage kit-009 (decision, members, pair ratios, verification) from the stage
# receipt under the proof root. Usage: Rscript docs/qa/workflow-execution/kit/kit-009/stage-summary.R
repo <- normalizePath(".")
root <- "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/runtime/kit-009"
stage <- readRDS(file.path(root, "stage.rds"))
decision <- readRDS(file.path(root, "decision.rds"))
members <- lapply(stats::setNames(nm = c("control-1", "candidate-1")), function(id) {
    op <- stage$cases[[id]]$operation
    v <- readRDS(file.path(root, id, "verification.rds"))
    list(status = op$status, variant = op$science_execution$variant, elapsed_seconds = op$elapsed_seconds,
        peak_bytes = op$peak_bytes, anonymous_peak_bytes = op$sampled_anon_bytes, swap_bytes = op$swap_bytes,
        max_sample_interval = op$max_sample_interval, batches = op$science_execution$batches,
        groups = op$science_execution$groups, engine_limit_bytes = op$science_execution$engine_limit_bytes,
        admission = op$science_execution$admission, verification = list(status = v$status, oracle = v$oracle,
            rows = v$state$rows, components = v$state$components, support = v$state$support,
            products_reported = v$state$products_reported))
})
record <- list(schema = "kit.stage_summary.v1", ticket = "KIT-009", stage = "kit-009", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), decision = decision$status,
    failing_rule = decision$failing_rule, pair = decision$pairs[["1"]], members = members,
    attempts = list(list(name = "attempt1-harness", outcome = "both members refused by the adapters vehicle before any work (journey population bound by its state descriptor); harness corrected; user amendment: both members rerun")))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-009/kit-009-stage.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("decision", decision$status, "\n")
