# KIT-012: the memory the DA step's packages add to the application process when their namespaces load.
# The DA step is the first step of the DIA journey that uses limma, openxlsx, htmlwidgets and Glimma; Glimma
# imports edgeR, DESeq2 and SummarizedExperiment, whose namespaces (and their Bioconductor dependencies) load
# into the application process and stay for its life. This records that cost in isolation, in a fresh process
# with the package loaded as the probe loads it. Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-012/da-namespace-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
anon <- function() {
    line <- grep("^RssAnon:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
    as.numeric(strsplit(trimws(sub("^[^:]+:", "", line)), " +")[[1L]][[1L]]) * 1024
}
measure <- function(label) {
    collected <- gc(full = TRUE)
    list(label = label, anon_bytes = anon(), live_r_heap = sum(collected[c("Ncells", "Vcells"), "used"] * c(56, 8)),
        ncells = collected["Ncells", "used"], vcells = collected["Vcells", "used"], namespaces = length(loadedNamespaces()))
}
steps <- list(measure("after load_all"))
before <- loadedNamespaces()
for (package in c("limma", "openxlsx", "htmlwidgets", "Glimma")) {
    invisible(requireNamespace(package, quietly = TRUE))
    steps[[length(steps) + 1L]] <- measure(paste("after", package))
}
loaded <- sort(setdiff(loadedNamespaces(), before))
first <- steps[[1L]]
last <- steps[[length(steps)]]
record <- list(schema = "kit.da_namespace_check.v1", ticket = "KIT-012", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), steps = steps, namespaces_loaded = loaded,
    result = list(anon_delta_bytes = last$anon_bytes - first$anon_bytes, live_r_heap_delta = last$live_r_heap - first$live_r_heap,
        ncells_delta = last$ncells - first$ncells, vcells_delta = last$vcells - first$vcells, namespaces_added = length(loaded)))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-012/da-namespace-check.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
for (step in steps) cat(sprintf("%-18s anon=%.3f GiB live=%.3f GiB namespaces=%d\n", step$label, step$anon_bytes / 2^30,
    step$live_r_heap / 2^30, step$namespaces))
cat(sprintf("RESULT %d namespaces, anonymous +%.3f GiB, live R heap +%.3f GiB: %s\n", length(loaded),
    record$result$anon_delta_bytes / 2^30, record$result$live_r_heap_delta / 2^30, paste(loaded, collapse = " ")))
