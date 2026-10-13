# KIT-009 post-hoc audit: where the imputation step's memory goes. The partitioned variant runs in process on
# a snapshot of the KIT-005 `journey_replicate` population (300 samples) while a shell sampler records this
# process's anonymous resident memory every 0.5 s; the batch loop's exit and the compute's exit are stamped
# (base::trace on the loaded functions) so the kernel phase (batches) and the finishing phase (output
# ordering, writers, the shared finishing step) can be told apart. A second run lowers the admitted cap so
# the query engine's share is 256 MB instead of 1 GB, to show whether the engine budget is used.
# Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-009/memory-attribution-probe.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
input <- b$populations$journey_replicate
parameters <- list(proportion_missing_values = 0.5, exclusion_column = NULL)
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit009-memory-probe")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
uptime <- function() as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])
samples <- file.path(work, "samples.tsv")
sampler <- sprintf("while kill -0 %d 2>/dev/null; do u=$(cut -d' ' -f1 /proc/uptime); a=$(grep RssAnon /proc/%d/status | awk '{print $2}'); echo \"$u\t$a\" >> %s; sleep 0.5; done",
    Sys.getpid(), Sys.getpid(), shQuote(samples))
system2("bash", c("-c", shQuote(sampler)), wait = FALSE)
marks <- list()
mark <- function(label) { marks[[length(marks) + 1L]] <<- list(label = label, uptime = uptime()); cat(sprintf("%-28s %.1f\n", label, uptime())) }
trace(partitionByFeatureGroup, exit = quote(mark("batches_done")), print = FALSE, where = asNamespace("MultiScholaR"))
trace(computeDIAReferencePeptideGroups, exit = quote(mark("compute_done")), print = FALSE, where = asNamespace("MultiScholaR"))
run <- function(label, worker_cap) {
    snapshot <- h$kitPopulationStore(input, file.path(work, label, "population-store"), api)
    identity <- readSavedStoreIdentity(file.path(snapshot, "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(snapshot), identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    reference <- getDIAStateReference(manager)
    gc(full = TRUE)
    mark(paste0(label, ":start"))
    value <- computeDIAReferencePeptideGroups(manager, reference, "peptideMissingValueImputation", parameters,
        worker_cap = worker_cap, variant = "partitioned")
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
        value$section, value$counts)
    mark(paste0(label, ":finish_done"))
    closeDIAStateAdapter(manager)
    gc(full = TRUE)
    mark(paste0(label, ":closed"))
    value$execution$engine_limit_bytes
}
mark("loaded")
engines <- list(normal_cap = run("normal_cap", 12 * 1024^3), reduced_cap = run("reduced_cap", 2 * 1024^3))
Sys.sleep(1)
trace <- utils::read.delim(samples, header = FALSE, col.names = c("uptime", "anon_kb"))
trace$anon <- trace$anon_kb * 1024
phase_peak <- function(from, to) {
    rows <- trace$uptime >= from & trace$uptime <= to
    if (!any(rows)) return(NA_real_)
    max(trace$anon[rows])
}
at <- function(label) Filter(function(m) identical(m$label, label), marks)[[1L]]$uptime
batches_marks <- Filter(function(m) identical(m$label, "batches_done"), marks)
compute_marks <- Filter(function(m) identical(m$label, "compute_done"), marks)
phases <- list()
for (i in seq_along(engines)) {
    label <- names(engines)[[i]]
    start <- at(paste0(label, ":start"))
    # The compute exit mark may be absent (base::trace does not always fire on the traced exit); the output
    # phase then runs to the finishing step's own mark.
    compute_end <- if (length(compute_marks) >= i) compute_marks[[i]]$uptime else at(paste0(label, ":finish_done"))
    phases[[label]] <- list(engine_limit_bytes = engines[[label]],
        anon_at_start = trace$anon[which.min(abs(trace$uptime - start))],
        kernel_phase_peak = phase_peak(start, batches_marks[[i]]$uptime),
        kernel_phase_seconds = round(batches_marks[[i]]$uptime - start, 1),
        output_phase_peak = phase_peak(batches_marks[[i]]$uptime, compute_end),
        output_phase_seconds = round(compute_end - batches_marks[[i]]$uptime, 1),
        finish_phase_peak = phase_peak(compute_end, at(paste0(label, ":finish_done"))),
        finish_phase_seconds = round(at(paste0(label, ":finish_done")) - compute_end, 1),
        anon_after_close = trace$anon[which.min(abs(trace$uptime - at(paste0(label, ":closed"))))])
}
record <- list(schema = "kit.memory_attribution_probe.v1", ticket = "KIT-009", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    population = list(name = "journey_replicate", rows = input$rows), parameters = parameters,
    anon_after_load = trace$anon[which.min(abs(trace$uptime - at("loaded")))], phases = phases, marks = marks,
    samples = nrow(trace))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-009/memory-attribution-probe.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
for (label in names(phases)) with(phases[[label]], cat(sprintf(
    "%s: engine %.0f MB; start %.2f GiB; kernel phase peak %.2f GiB (%.0f s); output phase peak %.2f GiB (%.0f s); finish phase peak %.2f GiB (%.0f s); after close %.2f GiB\n",
    label, engine_limit_bytes / 1e6, anon_at_start / 2^30, kernel_phase_peak / 2^30, kernel_phase_seconds,
    output_phase_peak / 2^30, output_phase_seconds, finish_phase_peak / 2^30, finish_phase_seconds, anon_after_close / 2^30)))
cat("RESULT written\n")
