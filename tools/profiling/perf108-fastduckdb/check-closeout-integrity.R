# Read-only final integrity check after the fixed schedule has terminated.
directory <- "docs/qa/workflow-execution/perf-108/closeout"
bindings <- readRDS(file.path(directory, "bindings.rds"))
stage <- file.path(bindings$root, "closure-filters")
stopifnot(file.exists(file.path(stage, "decision.rds")))
receipt <- readRDS(file.path(stage, "stage.rds"))
stopifnot(length(receipt$cases) == 8L)
groups <- lapply(names(receipt$cases), function(id) {
    config <- readRDS(file.path(stage, id, "application.rds"))
    path <- file.path("/sys/fs/cgroup", sub("^/", "", config$control_group))
    stopifnot(!dir.exists(path))
    list(member = id, control_group = config$control_group, absent = TRUE)
})
inputs <- lapply(names(bindings$full), function(population) {
    input <- bindings$full[[population]]
    head <- file.path(input$store, "dia-state-head.rds")
    payload <- file.path(input$store, input$reference$payload$relative_path)
    head_hash <- digest::digest(file = head, algo = "sha256")
    payload_hash <- digest::digest(file = payload, algo = "sha256")
    stopifnot(identical(head_hash, input$head_sha256),
        identical(payload_hash, input$reference$payload$sha256))
    list(population = population, head_sha256 = head_hash, payload_sha256 = payload_hash,
        bytes = file.size(payload), unchanged = TRUE)
})
for (role in c("accepted", "current")) {
    frozen <- bindings[[role]]$files
    stopifnot(all(vapply(names(frozen), function(path)
        identical(digest::digest(file = path, algo = "sha256"), frozen[[path]]), logical(1))))
}
jsonlite::write_json(list(status = "complete", member_cgroups = groups,
    original_inputs = inputs, installed_artifacts_unchanged = TRUE),
    file.path(directory, "cleanup.json"), auto_unbox = TRUE, pretty = TRUE)
