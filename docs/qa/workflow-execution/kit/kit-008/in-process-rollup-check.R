# KIT-008: the GUI rollup's backend switch (legacy row-indexed path -> the declaration's P14 ordered-sum
# aggregate, the backend the managed worker runs) checked page-identical on the screen P14 population
# before the member. Both backends run in process from the working tree on a snapshot store of the bound
# `screen_p14` population; their output pages are compared with each other and with the bound P14
# screen oracle. Not a counted run; the record is written beside this script.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-008/in-process-rollup-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
input <- b$populations$screen_p14
stopifnot(is.null(h$kitPopulationInvalid(input)))
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit008-rollup-check")
unlink(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
run <- function(backend) {
    snapshot <- h$kitPopulationStore(input, file.path(work, backend, "population-store"), api)
    identity <- readSavedStoreIdentity(file.path(snapshot, "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(snapshot), identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    reference <- getDIAStateReference(manager)
    started <- proc.time()[["elapsed"]]
    value <- computeDIAReferencePrecursorRollup(manager, reference, backend)
    seconds <- round(proc.time()[["elapsed"]] - started, 1)
    output <- value$checkpoint$components$peptide_data
    pages <- h$kitFingerprintTablePages(manager$.context$store, output, api)
    scope <- newWorkflowTableReadScope(manager$.context$store, output)
    rows <- scope$metadata$rows
    columns <- scope$metadata$columns
    releaseWorkflowTableReadScope(scope)
    list(backend = backend, seconds = seconds, rows = rows, columns = columns, pages = length(pages),
        page_hashes = unname(pages), summary_sha256 = digest::digest(value$summary, algo = "sha256"),
        profile = value$execution$profile$name, query_owner = value$execution$query_owner,
        exact_sum_sha256 = value$execution$exact_sum_sha256)
}
legacy <- run("legacy")
p14 <- run("P14")
oracle <- unlist(b$oracles$p14$screen_page_hashes, use.names = FALSE)
identical_backends <- identical(legacy$page_hashes, p14$page_hashes) && identical(legacy$rows, p14$rows) &&
    identical(legacy$columns, p14$columns) && identical(legacy$summary_sha256, p14$summary_sha256)
record <- list(schema = "kit.in_process_rollup_check.v1", ticket = "KIT-008", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    population = list(name = "screen_p14", rows = input$rows, entity_id = input$entity_id),
    legacy = legacy[c("seconds", "rows", "pages", "summary_sha256", "profile", "query_owner")],
    p14 = p14[c("seconds", "rows", "pages", "summary_sha256", "profile", "query_owner", "exact_sum_sha256")],
    page_identical = identical_backends,
    p14_equals_oracle = identical(p14$page_hashes, oracle) && length(oracle) == b$oracles$p14$screen_pages,
    legacy_equals_oracle = identical(legacy$page_hashes, oracle))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-008/in-process-rollup-check.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat(sprintf("legacy %.1fs rows=%d pages=%d | P14 %.1fs rows=%d pages=%d (%s/%s) | identical=%s oracle=%s\n",
    legacy$seconds, legacy$rows, legacy$pages, p14$seconds, p14$rows, p14$pages, p14$profile, p14$query_owner,
    identical_backends, record$p14_equals_oracle))
