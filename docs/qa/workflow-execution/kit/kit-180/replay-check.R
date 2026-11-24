# KIT-180 gate F2 (ticket section 8): does a label-free and a TMT record replay on this commit?
#
# Usage, from the package root on a run-tests runner (the commit under test checked out, full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-180/replay-check.R <out dir>
#   <out dir>  default docs/qa/workflow-execution/kit/kit-180/replay
#
# For prot_lfq and prot_tmt: the GUI-entry lane (tests/testthat/helper-gui-entry.R; helpers sourced into the global
# environment, no helper file edited) through its DA step, then one added step that writes the session's analysis record
# and copies the session's store and the uploaded inputs to scratch. The copy is cut to the record's import and design
# states (tools/profiling/workflow_kit.R `kitPrepareReplayProject`) and the record is replayed into it
# (`replayWorkflowAnalysisRecord`, inputs checked by sha256). Compared: the replayed state history with the record's,
# and each replayed state with the first run's by identical(). Written: <out dir>/replay.tsv (lane, commit, measure,
# value). A lane that stops is recorded with its condition class and message.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-180/replay"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-180 F2 at", commit, "\n")

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
suppressMessages(testthat::source_test_helpers("tests/testthat", env = globalenv()))
kit <- new.env(parent = globalenv())
sys.source(file.path("tools", "profiling", "workflow_kit.R"), kit)
api <- function(name) get(name, asNamespace("MultiScholaR"))

rows <- list()
put <- function(lane, measure, value) {
    value <- if (is.null(value) || !length(value)) NA_character_ else paste(as.character(value), collapse = ";")
    rows[[length(rows) + 1L]] <<- data.frame(lane = lane, commit = commit, measure = measure,
        value = gsub("[\t\r\n]+", " ", value), stringsAsFactors = FALSE)
    utils::write.table(do.call(rbind, rows), file.path(out_dir, "replay.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

scratch_root <- file.path(tempdir(), "kit180-replay")
dir.create(scratch_root, showWarnings = FALSE)
current <- new.env(parent = emptyenv())
original_steps <- guiEntrySteps
guiEntrySteps <- function(lane) {
    steps <- original_steps(lane)
    ids <- vapply(steps, `[[`, "", "id")
    c(steps[seq_len(which(ids == "da"))], list(guiEntryStep("kit180_keep", "da", function(api) {
        manager <- shiny::isolate(api$workflow_data$state_manager)
        keep <- file.path(scratch_root, lane$lane_id)
        dir.create(keep, recursive = TRUE, showWarnings = FALSE)
        json <- file.path(keep, "analysis_record.json")
        writeWorkflowAnalysisRecord(manager, json)
        # The session owner's own store (R0 correction: the store is not under the lane's root).
        store_root <- dirname(manager$.context$store$directory)
        current$store_count <- as.integer(dir.exists(file.path(store_root, ".workflow-entities")))
        current$origin <- file.path(keep, "origin")
        dir.create(current$origin, showWarnings = FALSE)
        stopifnot(current$store_count == 1L, system2("cp", c("-a", shQuote(file.path(store_root, ".")),
            shQuote(current$origin))) == 0L)
        current$inputs <- file.path(keep, "inputs")
        dir.create(current$inputs, showWarnings = FALSE)
        file.copy(list.files(api$inputs$scratch, full.names = TRUE), current$inputs, overwrite = TRUE)
        # The recorded inputs where the session read them (the import's uploads, a design file when one was imported).
        sources <- as.character(unlist(lapply(manager$.context$snapshot$file_dependencies, function(dependency) dependency$source)))
        sources <- unique(sources[!is.na(sources) & file.exists(sources)])
        if (length(sources)) file.copy(sources, current$inputs, overwrite = TRUE)
        current$json <- json
        current$history <- guiEntryHistory(manager)
        current$states <- lapply(stats::setNames(current$history, current$history), function(name)
            withWorkflowStateScope(manager, function() manager$getState(name)))
    })))
}

for (lane in c("prot_lfq", "prot_tmt")) {
    for (name in ls(current)) rm(list = name, envir = current)
    started <- Sys.time()
    where <- character()
    outcome <- tryCatch(withCallingHandlers({
        guiEntryRunLane(lane)
        record <- jsonlite::read_json(current$json, simplifyVector = FALSE)
        states <- vapply(record$steps, function(step) step$state, "")
        put(lane, "record_workflow", record$workflow)
        put(lane, "record_steps", paste(vapply(record$steps, function(step) step$operation, ""), collapse = ","))
        put(lane, "record_states", paste(states, collapse = ","))
        put(lane, "record_valid", validateWorkflowAnalysisRecord(record)$valid)
        put(lane, "store_count", current$store_count)
        preconditions <- states[seq_len(sum(cumprod(vapply(record$steps, function(step) step$operation, "") %in%
            c("import", "design_builder"))))]
        put(lane, "preconditions", paste(preconditions, collapse = ","))
        copy <- file.path(scratch_root, lane, "replay")
        dir.create(copy, showWarnings = FALSE)
        stopifnot(system2("cp", c("-a", shQuote(file.path(current$origin, ".")), shQuote(copy))) == 0L)
        prepared <- kit$kitPrepareReplayProject(copy, preconditions, api)
        inputs <- character()
        for (input in record$inputs) {
            if (is.null(input$role) || is.null(input$basename) || is.na(input$basename)) next
            path <- file.path(current$inputs, input$basename)
            if (file.exists(path)) inputs[[input$role]] <- path
        }
        put(lane, "inputs_found", paste(names(inputs), collapse = ","))
        paths <- list(source_dir = file.path(copy, "kit180-paths", "source"), results_dir = file.path(copy, "kit180-paths", "results"),
            da_output_dir = file.path(copy, "kit180-paths", "da"), publication_graphs_dir = file.path(copy, "kit180-paths", "pub"),
            protein_qc_dir = file.path(copy, "kit180-paths", "pqc"), time_dir = file.path(copy, "kit180-paths", "time"))
        for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
        replayed <- list()
        utils::capture.output(replay <- suppressMessages(suppressWarnings(replayWorkflowAnalysisRecord(record, prepared,
            inputs, paths, file.path(scratch_root, lane, "replay-scratch"), on_step = function(step, manager) {
                if (step$state %in% states) replayed[[step$state]] <<- withWorkflowStateScope(manager,
                    function() manager$getState(step$state))
            }))))
        put(lane, "replayed_states", paste(replay$states, collapse = ","))
        put(lane, "history_equal", identical(as.character(replay$states), states))
        later <- setdiff(states, preconditions)
        same <- vapply(later, function(state) identical(replayed[[state]], current$states[[state]]), logical(1))
        put(lane, "states_compared", length(later))
        put(lane, "states_identical", sum(same))
        put(lane, "states_not_identical", paste(later[!same], collapse = ","))
        "replayed"
    }, error = function(e) {
        # VAL-54 instrument run: record where a stop happened (the innermost calls), not only its message.
        calls <- vapply(utils::tail(sys.calls(), 8L), function(call) substr(paste(deparse(call, nlines = 1L), collapse = " "), 1L, 120L), "")
        where <<- calls
    }), error = function(e) {
        put(lane, "stopped", paste(class(e)[[1L]], conditionMessage(e), "| call:", paste(deparse(conditionCall(e), nlines = 1L), collapse = " ")))
        put(lane, "stopped_in", paste(where, collapse = " <- "))
        "stopped"
    })
    put(lane, "outcome", outcome)
    put(lane, "seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
cat("KIT-180 F2 done\n")
