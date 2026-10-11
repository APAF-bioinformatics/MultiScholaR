# Tracked summary of a kit-0NN operation stage (schema kit.operation_stage_summary.v1) from the stage receipt.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-006/operation-stage-summary.R <stage> [<output json>]
args <- commandArgs(TRUE)
stage <- args[[1L]]
repo <- normalizePath(".")
helpers <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), helpers)
bindings <- helpers$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
out <- if (length(args) > 1L) args[[2L]] else file.path(repo, "docs/qa/workflow-execution/kit", stage, paste0(stage, "-stage.json"))
receipt <- readRDS(file.path(bindings$root, stage, "stage.rds"))
gib <- function(x) if (is.numeric(x) && length(x) == 1L && is.finite(x)) round(x / 1024^3, 3) else NULL
member <- function(id) {
    item <- receipt$cases[[id]]
    op <- item$operation
    verification_path <- file.path(bindings$root, stage, id, "verification.rds")
    verification <- if (file.exists(verification_path)) readRDS(verification_path) else NULL
    list(id = id, status = op$status, operation = op$member$operation, population = op$member$population,
        build = op$config$build, elapsed_seconds = round(op$elapsed_seconds), peak_gib = gib(op$peak_bytes),
        anon_peak_gib = gib(op$sampled_anon_bytes), swap_bytes = op$swap_bytes, max_sample_interval_s = op$max_sample_interval,
        samples = op$samples, outer_guard = op$outer_guard, verifier_accepted = isTRUE(item$verifier_completion$accepted),
        verification_status = verification$status,
        batches = op$science_execution$admission$batches %||% verification$execution$admission$batches %||%
            (if (length(op$science_execution$admission$count_batch_sizes)) length(op$science_execution$admission$count_batch_sizes) else NULL),
        batch_sizes = op$science_execution$batch_sizes %||% verification$execution$batch_sizes %||%
            op$science_execution$admission$count_batch_sizes %||% verification$execution$admission$count_batch_sizes,
        admission_estimated_gib = gib(op$science_execution$admission$estimated_bytes %||% verification$execution$admission$estimated_bytes),
        worker_cap_gib = gib(op$science_execution$admission$worker_cap %||% verification$execution$admission$worker_cap),
        phase_seconds = lapply(op$phase_seconds, round), phase_peak_anon_gib = lapply(op$phase_peak_anon_bytes, gib),
        host_profile = op$host_at_start$profile, host_cpu_busy = round(op$host_at_start$cpu_busy, 3))
}
`%||%` <- function(x, y) if (is.null(x)) y else x
summary <- list(schema = "kit.operation_stage_summary.v1", stage = stage, status = receipt$status,
    control = bindings$control, builds = receipt$build, input_bindings_sha256 = receipt$input_bindings_sha256,
    decision = receipt$decision, decisions_before_rerun = receipt$decisions_before_rerun, reruns = receipt$reruns,
    corrections = receipt$corrections, members = lapply(names(receipt$cases), member))
jsonlite::write_json(summary, out, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("summary written:", out, "decision", receipt$decision$status, "\n")
