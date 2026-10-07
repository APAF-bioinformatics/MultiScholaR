# Reuse exact retained bytes and verified P14 controls. No rebuild or new timing.
path <- "docs/qa/workflow-execution/perf-108/closeout"
stopifnot(!file.exists(file.path(path, "bindings.rds")))
bindings <- readRDS("docs/qa/workflow-execution/perf-108/ordered-groups/bindings.rds")
rollup_path <- file.path(bindings$root, "ordered-rollup/stage.rds")
rollup <- readRDS(rollup_path)
source_bindings <- file.path(dirname(rollup_path), "source/docs/qa/workflow-execution/perf-108/ordered-groups/bindings.rds")
stopifnot(identical(digest::digest(file = source_bindings, algo = "sha256"),
    digest::digest(file = "docs/qa/workflow-execution/perf-108/ordered-groups/bindings.rds", algo = "sha256")))
for (role in c("accepted", "current")) {
    frozen <- bindings[[role]]
    stopifnot(all(vapply(names(frozen$files), function(file)
        identical(digest::digest(file = file, algo = "sha256"), frozen$files[[file]]), logical(1))))
    for (rep in 1:2) {
        item <- rollup$cases[[paste("p14", role, rep, sep = "-")]]
        stopifnot(isTRUE(item$passed), isTRUE(item$verifier_completion$accepted),
            identical(item$operation$config$build, frozen$build),
            isTRUE(item$operation$measurement_valid), isTRUE(item$operation$lifecycle_passed),
            identical(item$operation$swap_bytes, 0))
    }
}
for (rep in 1:2) {
    current <- rollup$cases[[paste0("p14-current-", rep)]]$operation
    accepted <- rollup$cases[[paste0("p14-accepted-", rep)]]$operation
    stopifnot(current$peak_bytes / accepted$peak_bytes <= 1.05,
        current$elapsed_seconds / accepted$elapsed_seconds <= 1.05)
}
bindings$reused_rollup <- list(path = rollup_path,
    sha256 = digest::digest(file = rollup_path, algo = "sha256"))
bindings$retained_rollup_eligible <- TRUE
bindings$root <- file.path(dirname(dirname(bindings$root)), "closeout/runtime")
bindings$library <- bindings$current$library
bindings$build <- bindings$current$build
bindings$candidate <- bindings$current
bindings$worker_options <- list()
bindings$finalist <- list(rollup = "P14", paired_state = FALSE, partition_validation = FALSE,
    persistent_reader = FALSE, native_output = FALSE, ordered_groups = FALSE)
bindings$protocol_sha256 <- digest::digest(file = file.path(path, "protocol.md"), algo = "sha256")
bindings$screen_members <- 0L
bindings$full_members <- 8L
saveRDS(bindings, file.path(path, "bindings.rds"))
jsonlite::write_json(list(library = bindings$library, root = bindings$root,
    roles = bindings[c("accepted", "current")], reused_rollup = bindings$reused_rollup,
    protocol_sha256 = bindings$protocol_sha256, experimental_options = bindings$finalist),
    file.path(path, "freeze.json"), auto_unbox = TRUE, pretty = TRUE)
