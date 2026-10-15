# KIT-036: the DA volcano page written by the two writers on the 300-sample project of the KIT-017 counted pair
# (candidate-1's store, head `da_published`), each in a fresh R process: the process's peak resident memory
# (VmHWM) and, for the pandoc writer, the peak resident memory of its child processes sampled at 50 ms. Not a
# guarded run; the record is written beside this script. Usage:
#   R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-036/page-writer-probe.R <member dir> <label>
args <- commandArgs(TRUE)
if (identical(args[[1L]], "--child")) {
    writer <- args[[2L]]; member <- args[[3L]]; out <- args[[4L]]
    # The working tree's package (the writer under test), with the proof library on R_LIBS for its dependencies.
    suppressMessages(devtools::load_all(quiet = TRUE))
    options(multischolar.widget_pandoc = identical(writer, "pandoc"))
    project <- file.path(member, "project", "dia_baseline_smoke")
    source_dir <- file.path(project, "scripts", "proteomics")
    roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
    identity <- MultiScholaR:::readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
    manager <- MultiScholaR:::newDIAStateAdapter(MultiScholaR:::newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
    manager$.context$session_format <- "native_v1"; manager$.context$reference_states <- TRUE
    hwm <- function() as.numeric(sub("^VmHWM:\\s*([0-9]+) kB", "\\1", grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE))) * 1024
    before <- hwm()
    directory <- tempfile("kit036-"); dir.create(directory)
    started <- proc.time()[["elapsed"]]
    files <- MultiScholaR:::withDIAStateScope(manager, function() {
        manifest <- MultiScholaR:::resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
        stopifnot(MultiScholaR:::isDIAResultManifest(manifest))
        MultiScholaR:::writeDIAInteractiveVolcanoOutputs(manifest, directory, data.frame(), 0.05, manager = NULL)
    })
    seconds <- round(proc.time()[["elapsed"]] - started, 1)
    html <- files[grepl("[.]html$", files)]
    record <- list(writer = writer, seconds = seconds, r_peak_before_bytes = before, r_peak_bytes = hwm(),
        page_bytes = file.size(html), page = basename(html), glimma_loaded = "Glimma" %in% loadedNamespaces())
    jsonlite::write_json(record, out, auto_unbox = TRUE)
    quit(status = 0)
}
member <- args[[1L]]; label <- args[[2L]]
poll_children <- function(pid) {
    max_child <- 0; names <- character()
    repeat {
        # A finished child stays a zombie until reaped, so completion is read from sys, not /proc.
        if (!is.na(sys::exec_status(pid, wait = FALSE))) break
        children <- system2("pgrep", c("-P", pid), stdout = TRUE, stderr = FALSE)
        for (child in children) {
            status <- tryCatch(readLines(file.path("/proc", child, "status"), warn = FALSE), error = function(e) character())
            rss <- as.numeric(sub("^VmRSS:\\s*([0-9]+) kB", "\\1", grep("^VmRSS:", status, value = TRUE))) * 1024
            name <- sub("^Name:\\s*", "", grep("^Name:", status, value = TRUE))
            if (length(rss) && rss > max_child) { max_child <- rss; names <- name }
        }
        Sys.sleep(0.05)
    }
    list(child_peak_bytes = max_child, child = names)
}
results <- list()
for (writer in c("inline", "pandoc")) {
    out <- tempfile(fileext = ".json")
    pid <- sys::exec_background(file.path(R.home("bin"), "Rscript"), c("--vanilla", "docs/qa/workflow-execution/kit/kit-036/page-writer-probe.R", "--child", writer, member, out),
        std_out = FALSE, std_err = FALSE)
    children <- poll_children(pid)  # reaps the child when it ends
    stopifnot(file.exists(out))
    results[[writer]] <- c(jsonlite::read_json(out), children)
}
record <- list(schema = "kit.page_writer_probe.v1", ticket = "KIT-036", label = label, date = format(Sys.Date()),
    library = dirname(system.file(package = "MultiScholaR")), member = member, writers = results)
jsonlite::write_json(record, paste0("docs/qa/workflow-execution/kit/kit-036/page-writer-probe-", label, ".json"), auto_unbox = TRUE, pretty = TRUE)
for (w in names(results)) { r <- results[[w]]; cat(sprintf("%-7s %5.1f s  page %.1f MB  R peak %.2f GiB (before %.2f)  child peak %.2f GiB (%s)  glimma %s\n", w, r$seconds, r$page_bytes / 1e6, r$r_peak_bytes / 2^30, r$r_peak_before_bytes / 2^30, r$child_peak_bytes / 2^30, paste(r$child, collapse = ","), r$glimma_loaded)) }
