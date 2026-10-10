# KIT-006: physical-mode sample filter (adaptive complete-sample batches through the generic partitioner)
# checked against the unchanged native method on a screen population before the full members. The
# population is the KIT-001 screen journey's intensity_filtered state (8 samples; the state the sample
# filter consumes), materialised as a snapshot store by the harness so the evidence store is never written.
# Both modes run in-process from the working tree; the record is written beside this script.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-006/in-process-screen-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
screen <- readRDS(file.path(b$root, "kit-lifecycle-complete", "member", "journey-result.rds"))$populations$intensity
stopifnot(identical(screen$state, "intensity_filtered"))
input <- list(store = screen$store, head_sha256 = screen$head_sha256, state = screen$state, entity_id = screen$entity_id,
    head_entity = h$kitFindHeadEntity(screen$store, screen$head_sha256))
stopifnot(!is.null(input$head_entity))
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit006-screen-check")
unlink(work, recursive = TRUE)
snapshot <- h$kitPopulationStore(input, file.path(work, "population-store"), api)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
identity <- readSavedStoreIdentity(file.path(snapshot, "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(snapshot), identity$project_id, identity$workflow_id), reopen = TRUE)
on.exit(closeDIAStateAdapter(manager), add = TRUE, after = FALSE)
reference <- getDIAStateReference(manager)
head <- h$kitReadHead(snapshot, api)
stopifnot(identical(head$current_state, "intensity_filtered"), identical(head$states$intensity_filtered, screen$entity_id))
parameters <- list(peptides_per_sample_cutoff = b$journey$settings$sample[["sample_filter-min_peptides_per_sample"]], inclusion_list = NULL)
fingerprint <- function(value) {
    checkpoint <- value$checkpoint
    h$kitFingerprintFilterResult(manager, checkpoint, api)
}
timed <- function(mode, cap = 12 * 1024^3) {
    started <- proc.time()[["elapsed"]]
    value <- computeDIAReferenceSampleFilter(manager, reference, parameters, worker_cap = cap, output_mode = mode)
    list(mode = mode, seconds = round(proc.time()[["elapsed"]] - started, 1), fingerprint = fingerprint(value),
        summary = value$summary, admission = value$execution$admission, unit = value$execution$bulk_query$independence_unit,
        declared_unit = value$execution$bulk_query$declared_unit)
}
native <- timed("native")
physical <- timed("physical")
compare <- function(a, b) {
    fields <- c("page_hashes", "rows", "samples", "columns", "design_sha256", "summary_sha256", "support_sha256")
    vapply(fields, function(f) identical(a$fingerprint[[f]], b$fingerprint[[f]]), logical(1))
}
same <- compare(native, physical)
# A reduced cap forces more than the floor of three batches on the same input; the result must not move.
estimates <- physical$admission$estimates
cap <- (estimates[[4L]] + estimates[[3L]]) / 2
reduced <- timed("physical", cap)
same_reduced <- compare(native, reduced)
record <- list(schema = "kit.in_process_screen_check.v1", ticket = "KIT-006", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    population = list(source = "kit-lifecycle-complete member, intensity_filtered", rows = screen$rows, entity_id = screen$entity_id,
        head_entity = input$head_entity, samples = native$fingerprint$samples),
    parameters = parameters,
    native = native[c("seconds", "unit")], physical = physical[c("seconds", "unit", "declared_unit", "admission")],
    reduced_cap = list(cap = cap, seconds = reduced$seconds, batches = reduced$admission$batches,
        count_batch_sizes = reduced$admission$count_batch_sizes),
    identical_to_native = list(normal_cap = same, reduced_cap = same_reduced),
    retained = list(rows = native$fingerprint$rows, samples_retained = native$summary$retained_sample_count,
        samples_removed = native$summary$removed_sample_count))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-006/in-process-screen-check.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
cat("native", native$seconds, "s; physical", physical$seconds, "s, batches", physical$admission$batches,
    paste(physical$admission$count_batch_sizes, collapse = "/"), "; reduced cap batches", reduced$admission$batches, "\n")
cat("identical at the normal cap:", all(same), "; at the reduced cap:", all(same_reduced), "\n")
stopifnot(all(same), all(same_reduced), physical$admission$batches == 3L, reduced$admission$batches == 4L)
