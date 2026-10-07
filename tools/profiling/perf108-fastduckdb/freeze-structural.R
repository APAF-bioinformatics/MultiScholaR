# Freeze existing accepted inputs and exact installed bytes; no package rebuild.
ordered <- identical(commandArgs(trailingOnly = TRUE), "ordered")
path <- file.path("docs/qa/workflow-execution/perf-108", if (ordered) "ordered-groups" else "structural")
stopifnot(!file.exists(file.path(path, "bindings.rds")))
old <- readRDS("docs/qa/workflow-execution/perf-108/cache-qualification/bindings.rds")
candidate_library <- if (ordered)
    "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/ordered-groups/browser-smoke.jEmE2F/library" else
    "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/structural/browser-smoke.IsZngp/library"
inspect <- function(library, native) callr::r(function(library, native) {
    .libPaths(c(library, .libPaths()))
    ns <- loadNamespace("MultiScholaR")
    artifact <- get(if (native) "diaRollupNativeInfo" else "diaRollupExtensionInfo", ns)()
    packages <- c("MultiScholaR", "fastduckdb", "duckdb", "arrow", "nanoarrow", "duckplyr")
    files <- unlist(lapply(packages, function(package) {
        root <- find.package(package)
        selected <- list.files(root, recursive = TRUE, full.names = TRUE)
        selected[grepl("[.](so|rdb|rdx)$|/DESCRIPTION$|/NAMESPACE$|manifest[.]rds$", selected)]
    }))
    list(library = normalizePath(library), build = get("workflowExecutionBuildIdentity", ns)(),
        exact_sum = artifact$manifest$binary_sha256,
        versions = stats::setNames(vapply(packages, function(package)
            as.character(utils::packageVersion(package)), character(1)), packages),
        files = stats::setNames(vapply(files, function(file) digest::digest(file = file,
            algo = "sha256"), character(1)), files))
}, args = list(library, native), libpath = c(library, .libPaths()),
    system_profile = FALSE, user_profile = FALSE)
accepted <- inspect(old$accepted$library, FALSE)
candidate <- inspect(candidate_library, TRUE)
current <- if (ordered) candidate else inspect(old$candidate$library, TRUE)
stopifnot(identical(accepted$build, old$accepted$build),
    identical(accepted$exact_sum, old$accepted$exact_sum))
if (!ordered) stopifnot(identical(current$build, old$candidate$build),
    identical(current$exact_sum, old$candidate$exact_sum))
bindings <- old
bindings$root <- file.path(dirname(dirname(candidate_library)), "runtime")
bindings$library <- candidate_library
bindings$build <- candidate$build
bindings$accepted <- accepted
bindings$current <- current
bindings$candidate <- candidate
bindings$protocol_sha256 <- digest::digest(file = file.path(path, "protocol.md"), algo = "sha256")
bindings$screen_members <- 4L
bindings$full_members <- 18L
bindings$thresholds <- list(memory_peak = 0.8, memory_elapsed = 1.15,
    time_elapsed = 0.9, time_peak = 1.05, non_regression = 1.05,
    reference_elapsed_variability = 0.15, reference_peak_variability = 0.10)
bindings$finalist <- list(rollup = "O", primary = "unchanged", wide = "unchanged",
    paired_state = FALSE, partition_validation = FALSE, persistent_reader = FALSE)
if (ordered) {
    bindings$finalist$rollup <- "S"
    bindings$finalist$native_output <- FALSE
    bindings$worker_options <- list(multischolar.perf108.ordered_groups = TRUE)
    bindings$historical_source_hashes <- bindings$source_hashes
    bindings$source_hashes <- NULL
}
saveRDS(bindings, file.path(path, "bindings.rds"))
stopifnot(file.copy(file.path(path, "protocol.md"), file.path(path, "confirmation-protocol-frozen.md")))
jsonlite::write_json(list(roles = bindings[c("accepted", "current", "candidate")],
    root = bindings$root, protocol_sha256 = bindings$protocol_sha256,
    finalist = bindings$finalist), file.path(path, "freeze.json"), auto_unbox = TRUE, pretty = TRUE)
