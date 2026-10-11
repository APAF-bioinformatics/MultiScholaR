# Tracked summary of a ticket journey stage (schema kit.journey_stage_summary.v2) from the stage receipt.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-011/journey-stage-summary.R <stage> [<output json>]
args <- commandArgs(TRUE)
stage <- args[[1L]]
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
out <- if (length(args) > 1L) args[[2L]] else file.path(repo, "docs/qa/workflow-execution/kit", sub("^kit-journey-", "kit-", stage), paste0(stage, "-stage.json"))
receipt <- readRDS(file.path(bindings$root, stage, "stage.rds"))
gib <- function(x) if (is.numeric(x) && length(x) == 1L && is.finite(x)) round(x / 1024^3, 3) else NULL
member <- function(id) {
    item <- receipt$cases[[id]]
    op <- item$operation
    verification_path <- file.path(bindings$root, stage, id, "verification.rds")
    verification <- if (file.exists(verification_path)) readRDS(verification_path) else NULL
    list(id = id, status = op$status, build = op$config$build, record = isTRUE(op$member$record), elapsed_seconds = round(op$elapsed_seconds),
        peak_gib = gib(op$peak_bytes), anon_peak_gib = gib(op$sampled_anon_bytes), resident_end_gib = gib(op$resident_end_bytes),
        swap_bytes = op$swap_bytes, max_sample_interval_s = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms, samples = op$samples,
        outer_guard = op$outer_guard, populations_complete = op$populations_complete,
        verifier_accepted = isTRUE(item$verifier_completion$accepted), verification_status = verification$status,
        oracle_mismatches = verification$oracle_mismatches, operation_identity_products_reported = length(verification$operation_identity_products),
        canonical_form_differences = verification$canonical_form_differences,
        phase_seconds = lapply(op$phase_seconds, round), phase_peak_anon_gib = lapply(op$phase_peak_anon_bytes, gib), phase_status = op$phase_status,
        host_profile = op$host_at_start$profile, host_cpu_busy = round(op$host_at_start$cpu_busy, 3))
}
summary <- list(schema = "kit.journey_stage_summary.v2", stage = stage, status = receipt$status, control = bindings$control,
    input_bindings_sha256 = receipt$input_bindings_sha256, decision = receipt$decision,
    decisions_before_correction = lapply(receipt$decisions_before_correction, `[[`, "status"), corrections = receipt$corrections,
    members = lapply(names(receipt$cases), member))
jsonlite::write_json(summary, out, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("summary written:", out, "decision", receipt$decision$status, "\n")
