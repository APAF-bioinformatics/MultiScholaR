# KIT-054: compare two browser members' publication files (the new lane on the owner path beside the legacy lane in
# memory) with the harness's file comparison (content hashes; PDFs and the session bookkeeping by path; paths with
# run ids and timestamps canonicalised). Reports counts and paths only (data handling rule: no file content is read
# here beyond the hashes the members recorded).
# Usage: Rscript docs/qa/workflow-execution/kit/kit-054/lane-compare.R <label> <new lane stage dir> <legacy lane stage dir>
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
label <- args[[1L]]; new_stage <- args[[2L]]; legacy_stage <- args[[3L]]
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
member <- function(stage) {
    jr <- readRDS(file.path(stage, "candidate-1", "journey-result.rds"))
    receipt <- readRDS(file.path(stage, "stage.rds"))
    list(stage = stage, files = jr$files, status = jr$status, completed = jr$completed,
        store = jr$store %||% NA_character_, decision = tryCatch(readRDS(file.path(stage, "decision.rds"))$status, error = function(e) NA_character_),
        elapsed = receipt$cases[["candidate-1"]]$operation$elapsed_seconds %||% NA_real_)
}
a <- member(new_stage); b <- member(legacy_stage)
# The project's store, scratch and logs are excluded by the harness; the new lane's store-backed exports differ by
# design (the legacy lane writes no store), so they are reported separately.
differing <- h$kitFileMismatches(b$files, a$files)
canonical <- function(path) { path <- gsub("[0-9]{8}_[0-9]{6}", "<time>", path); gsub("[0-9a-f]{13,}", "<id>", path) }
paths_a <- canonical(a$files$path); paths_b <- canonical(b$files$path)
common <- intersect(paths_a, paths_b)
only_new <- setdiff(paths_a, paths_b); only_legacy <- setdiff(paths_b, paths_a)
directory <- function(path) sub("^(([^/]+/){0,3}).*$", "\\1", path)
by_dir <- function(paths) if (!length(paths)) list() else as.list(table(vapply(paths, directory, "")))
diff_common <- intersect(canonical(differing), common)
summary <- list(schema = "kit.lane_compare.v1", ticket = "KIT-054", label = label, date = format(Sys.Date()),
    new_lane = list(stage = new_stage, status = a$status, decision = a$decision, files = nrow(a$files), elapsed_seconds = a$elapsed, store = !is.na(a$store)),
    legacy_lane = list(stage = legacy_stage, status = b$status, decision = b$decision, files = nrow(b$files), elapsed_seconds = b$elapsed, store = !is.na(b$store)),
    common_files = length(common), equal_common = length(common) - length(diff_common), differing_common = length(diff_common),
    only_new_lane = length(only_new), only_legacy_lane = length(only_legacy),
    differing_paths = diff_common, only_new_lane_by_directory = by_dir(only_new), only_legacy_lane_by_directory = by_dir(only_legacy),
    differing_by_directory = by_dir(diff_common))
jsonlite::write_json(summary, file.path("docs/qa/workflow-execution/kit/kit-054", paste0("lane-compare-", label, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("%s: %d common files, %d equal, %d differing; %d only on the new lane, %d only on the legacy lane\n", label,
    length(common), length(common) - length(diff_common), length(diff_common), length(only_new), length(only_legacy)))
if (length(diff_common)) cat("differing (common) files:\n", paste(" ", diff_common, collapse = "\n"), "\n")
