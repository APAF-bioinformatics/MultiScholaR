# KIT-266 gate G0 (observation of the BASE; no product code, no stored file touched). Declaration: declaration.md.
# Lane: guiEntryRunLane("prot_dia_limpa") (the GUI-entry harness), observed through traces; then two fresh R processes
# measure the memory of building the combined figure (memory-build.R).
# Usage (run-tests scripts job, from the package root): Rscript docs/qa/workflow-execution/kit/kit-266/g0/g0.R

pkgroot <- normalizePath(getwd())
dir.create(file.path(pkgroot, "kit266-g0"), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(file.path(pkgroot, "kit266-g0"))          # absolute: a trace must never write to a relative path
memscript <- file.path(pkgroot, "docs", "qa", "workflow-execution", "kit", "kit-266", "g0", "memory-build.R")
imgdir <- file.path(tempdir(), "kit266-g0-images")
events_log <- file.path(out, "events.log")
say <- function(...) cat(..., "\n", sep = "")
writeUtf8 <- function(lines, path) writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
`%||%` <- function(a, b) if (!is.null(a)) a else b
one <- function(x) if (length(x)) paste(as.character(x), collapse = ",") else "<none>"

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
ns <- asNamespace("MultiScholaR")
writeUtf8(c(paste("R", R.version.string), paste("HEAD", paste(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), collapse = " ")),
    paste("BASE_BOUND 894b9cc82 (declared); R/ unchanged since it"), paste("out", out), paste("imgdir", imgdir)),
    file.path(out, "meta.txt"))
helperFn <- function(name) if (exists(name, envir = ns, inherits = FALSE)) get(name, envir = ns) else get(name, envir = globalenv())
runLane <- helperFn("guiEntryRunLane")

# --- Predictions (declared before the run; declaration.md) --------------------------------------------------------------
predictions <- data.frame(id = paste0("P", 1:7), predicted = c(
    "at step norm only post_norm_{correlation,density,pca,rle}.png exist among the pre/post/ruv/composite figures of protein_qc",
    "pre_norm_* files first exist at step ruv_covariates",
    "no composite_QC_figure file in the project at any step",
    "the final pre_normalization source entity equals exactly the correlation_filtered state's entity (matching states: correlation_filtered)",
    "copyToResultsSummary's returned list names Composite QC (PDF) and Composite QC (PNG) as missing",
    "the tab entry fires at step ruv_covariates with current_state == correlation_filtered",
    "Reset fallback (normalization_input_reference NULL) selects correlation_filtered"), stringsAsFactors = FALSE)

# --- The recorder and the traces ----------------------------------------------------------------------------------------
obs <- new.env(parent = emptyenv())
obs$events <- list()
obs$base_dir <- NULL
obs$end_done <- FALSE
obs$p4 <- NULL
obs$p7 <- NULL
fmt <- function(x) paste(vapply(names(x), function(n) paste0(n, "=", one(x[[n]])), ""), collapse = "; ")
rec <- function(type, ...) {
    ev <- c(list(type = type), list(...))
    obs$events[[length(obs$events) + 1L]] <- ev
    cat(sprintf("%s %s\n", format(Sys.time(), "%H:%M:%OS3"), fmt(ev)), file = events_log, append = TRUE)
    invisible(NULL)
}

endOfLane <- function(caller) {
    if (isTRUE(obs$end_done)) return(invisible(NULL))
    obs$end_done <- TRUE
    manager <- tryCatch(get("manager", envir = caller)(), error = function(e) NULL)
    if (is.null(manager)) { rec("end_of_lane", error = "no manager"); return(invisible(NULL)) }
    scope <- ns$withWorkflowStateScope
    history <- unique(helperFn("guiEntryHistory")(manager))
    # P4: the pre-normalisation source's entity against the entities of the active history's states.
    obs$p4 <- tryCatch(scope(manager, function() {
        src <- manager$.context$snapshot$qc_sources$pre_normalization
        src_id <- one(src$entity_id)
        ids <- vapply(history, function(s) one(ns$getWorkflowStateReference(manager, s)$entity_id), "")
        list(history = history, src_id = src_id, ids = ids, matching = history[ids == src_id])
    }), error = function(e) list(error = conditionMessage(e)))
    rec("p4_source", head = tryCatch(manager$current_state, error = function(e) "<error>"), history = history,
        src_id = obs$p4$src_id %||% NA, state_ids = if (!is.null(obs$p4$ids)) paste(history, obs$p4$ids, sep = ":", collapse = " | ") else NA,
        matching = if (!is.null(obs$p4$matching)) one(obs$p4$matching) else NA, error = obs$p4$error %||% "")
    # The stage images, kept for the fresh memory processes (the lane's project is deleted when the lane returns).
    qcdir <- file.path(obs$base_dir, "results", "proteomics", "protein_qc")
    dir.create(imgdir, recursive = TRUE, showWarnings = FALSE)
    pngs <- list.files(qcdir, pattern = "\\.png$", full.names = TRUE)
    file.copy(pngs, imgdir, overwrite = TRUE)
    rec("images", from = qcdir, copied = sort(list.files(imgdir)))
    # P7 (last: it moves the head): the Reset fallback with normalization_input_reference NULL.
    obs$p7 <- tryCatch({
        res <- scope(manager, function() ns$revertProtNormStateManagerToPreNormalization(
            list(state_manager = manager), list(normalization_input_reference = NULL), messageFn = function(...) NULL))
        list(previous = res$previousState %||% "<NULL>")
    }, error = function(e) list(previous = paste0("<error> ", class(e)[[1L]], ": ", conditionMessage(e))))
    rec("p7_reset", previous_state = obs$p7$previous)
    invisible(NULL)
}

onArtifacts <- function(base_dir, caller_unused) {
    # The capture frame is found by its variables (parent.frame() inside a trace's tracer is not the caller).
    caller <- NULL
    for (f in rev(sys.frames())) if (exists("step", envir = f, inherits = FALSE) && exists("refusal_before", envir = f, inherits = FALSE)) { caller <- f; break }
    if (is.null(caller)) { rec("capture", step = NA_character_, error = "no capture frame"); return(invisible(NULL)) }
    step <- tryCatch(as.character(get("step", envir = caller)$id)[[1L]], error = function(e) NA_character_)
    files <- list.files(base_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    qc <- sort(basename(files[grepl("^results/proteomics/protein_qc/[^/]+$", files)]))
    obs$base_dir <- base_dir
    rec("capture", step = step, n_qc = length(qc), qc = qc, composite_anywhere = grep("composite_QC_figure", files, value = TRUE))
    if (identical(step, "session_export")) endOfLane(caller)
    invisible(NULL)
}
onTab <- function(selectedTab, workflowData) {
    state <- tryCatch(shiny::isolate(workflowData$state_manager$current_state), error = function(e) paste("<error>", conditionMessage(e)))
    rec("tab_entry", selected = selectedTab %||% "<NULL>", current_state = state %||% "<NULL>")
}
onTabExit <- function(rv) { force(rv); rec("tab_entry_exit", returned = tryCatch(one(rv), error = function(e) "<unprintable>")) }
onCopy <- function(omic_type) rec("copy_entry", omic_type = omic_type)
onCopyExit <- function(rv) {
    force(rv)
    shown <- tryCatch(vapply(rv, function(x) as.character(x$display_name)[[1L]], ""), error = function(e) character())
    rec("copy_exit", returned_class = class(rv)[[1L]], n = length(rv), failed_copies = paste(shown, collapse = " | "))
}

traceStatus <- character()
safeTrace <- function(what, where, ...) {
    ok <- tryCatch({ suppressMessages(trace(what, where = where, print = FALSE, ...)); TRUE }, error = function(e) {
        say("trace failed for ", what, ": ", conditionMessage(e))
        FALSE })
    traceStatus[[what]] <<- if (ok) "installed" else "FAILED"
    invisible(ok)
}
safeTrace("guiEntryArtifacts", runLane, tracer = bquote(.(onArtifacts)(base_dir, NULL)))
safeTrace("runProtNormTabEntryWorkflow", ns, tracer = bquote(.(onTab)(selectedTab, workflowData)),
    exit = bquote(.(onTabExit)(returnValue("<no return>"))))
safeTrace("copyToResultsSummary", ns, tracer = bquote(.(onCopy)(omic_type)),
    exit = bquote(.(onCopyExit)(returnValue("<no return>"))))
safeTrace("generateProtNormPreNormalizationQcArtifacts", ns, tracer = bquote(.(rec)("pre_generator")))
say("traces: ", paste(names(traceStatus), traceStatus, sep = "=", collapse = "; "))

# --- Positive controls: every trace must fire and write before a reading is trusted -----------------------------------
control <- character()
ctrl <- function(name, ok, detail) {
    control[[name]] <<- sprintf("CONTROL %s: %s (%s)", name, if (isTRUE(ok)) "PASS" else "BLIND", detail)
    say(control[[name]])
}
eventsOf <- function(type) Filter(function(e) identical(e$type, type), obs$events)
logLines <- function() if (file.exists(events_log)) readLines(events_log) else character()

ctrlDir <- file.path(tempdir(), "kit266-control")
dir.create(file.path(ctrlDir, "results", "proteomics", "protein_qc"), recursive = TRUE)
for (f in c("post_norm_pca.png", "pre_norm_pca.png")) file.create(file.path(ctrlDir, "results", "proteomics", "protein_qc", f))
controlCapture <- function(step, refusal_before = NULL) helperFn("guiEntryArtifacts")(ctrlDir, character())
tryCatch(controlCapture(list(id = "control-step")), error = function(e) say("control capture error: ", conditionMessage(e)))
ev <- eventsOf("capture")
ctrl("artifacts_trace", length(ev) == 1L && identical(ev[[1L]]$step, "control-step") && identical(ev[[1L]]$qc, c("post_norm_pca.png", "pre_norm_pca.png")) &&
    any(grepl("type=capture; step=control-step", logLines(), fixed = TRUE)),
    paste("events:", length(ev), "; step:", if (length(ev)) ev[[1L]]$step else NA, "; log lines:", length(logLines())))

try(ns$runProtNormTabEntryWorkflow("other", list(), list(), function() NULL), silent = TRUE)
ev <- eventsOf("tab_entry"); ev2 <- eventsOf("tab_entry_exit")
ctrl("tab_entry_trace", length(ev) == 1L && identical(ev[[1L]]$selected, "other") && length(ev2) == 1L,
    paste("entry events:", length(ev), "; exit events:", length(ev2), "; exit returned:", if (length(ev2)) ev2[[1L]]$returned else NA))

try(ns$copyToResultsSummary(omic_type = "", experiment_label = "kit266-control"), silent = TRUE)
ev <- eventsOf("copy_entry"); ev2 <- eventsOf("copy_exit")
ctrl("copy_trace", length(ev) == 1L && length(ev2) == 1L,
    paste("entry events:", length(ev), "; exit events:", length(ev2), "; exit class:", if (length(ev2)) ev2[[1L]]$returned_class else NA))

ctrl("lane_harness", is.function(runLane) && is.function(helperFn("guiEntryComparable")) && file.exists(memscript) &&
    file.exists(file.path(pkgroot, "tests", "testdata", "gui-entry", "prot_dia_limpa.json")),
    paste("memscript:", file.exists(memscript)))

rscript <- file.path(R.home("bin"), "Rscript")
runMem <- function(mode) {
    lines <- suppressWarnings(system2(rscript, c(shQuote(memscript), mode, shQuote(imgdir)), stdout = TRUE, stderr = TRUE))
    writeUtf8(lines, file.path(out, paste0("memory-", mode, ".txt")))
    cat(lines, sep = "\n")
    lines
}
memCtl <- runMem("control")
ctrl("memory_measure", any(grepl("MEMORY_CONTROL .* verdict=PASS", memCtl)), paste(grep("MEMORY_CONTROL", memCtl, value = TRUE), collapse = " "))

writeUtf8(unname(control), file.path(out, "controls.txt"))
if (!all(grepl(": PASS ", control, fixed = TRUE)) || any(traceStatus[c("guiEntryArtifacts", "runProtNormTabEntryWorkflow", "copyToResultsSummary")] != "installed")) {
    say("INSTRUMENT BLIND: a positive control failed or a trace did not install; the lane is not run. Stop and report.")
    quit(save = "no", status = 1L)
}
obs$events <- list()                       # the controls' events are not readings
unlink(events_log)

# --- The lane ------------------------------------------------------------------------------------------------------------
started <- Sys.time()
record <- tryCatch(runLane("prot_dia_limpa"), error = function(e) e)
say(sprintf("lane finished in %.0f s: %s", as.numeric(difftime(Sys.time(), started, units = "secs")),
    if (inherits(record, "error")) paste("ERROR:", conditionMessage(record)) else "ok"))
laneOk <- !inherits(record, "error")
events <- obs$events

# --- Events by step (an event belongs to the capture that follows it) ---------------------------------------------------
pending <- list(); rows <- list(); stepQc <- list(); stepComposite <- list(); stepOrder <- character()
for (e in events) {
    if (identical(e$type, "capture")) {
        stepOrder <- c(stepOrder, e$step)
        stepQc[[e$step]] <- e$qc
        stepComposite[[e$step]] <- e$composite_anywhere
        for (p in pending) rows[[length(rows) + 1L]] <- data.frame(step = e$step, type = p$type, detail = fmt(p[names(p) != "type"]), stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(step = e$step, type = "capture", detail = fmt(e[c("n_qc", "composite_anywhere")]), stringsAsFactors = FALSE)
        pending <- list()
    } else pending[[length(pending) + 1L]] <- e
}
for (p in pending) rows[[length(rows) + 1L]] <- data.frame(step = "<after last capture>", type = p$type, detail = fmt(p[names(p) != "type"]), stringsAsFactors = FALSE)
eventTable <- if (length(rows)) do.call(rbind, rows) else data.frame(step = character(), type = character(), detail = character())
utils::write.table(eventTable, file.path(out, "events.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
steps <- data.frame(step = stepOrder, n_protein_qc_files = vapply(stepOrder, function(s) length(stepQc[[s]]), 0L),
    protein_qc_files = vapply(stepOrder, function(s) one(stepQc[[s]]), ""), composite_QC_figure_files = vapply(stepOrder, function(s) one(stepComposite[[s]]), ""),
    stringsAsFactors = FALSE)
utils::write.table(steps, file.path(out, "steps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("per-step protein_qc files:")
print(steps[, c("step", "n_protein_qc_files", "protein_qc_files")], row.names = FALSE, right = FALSE)

# --- P1 to P7 -------------------------------------------------------------------------------------------------------------
figureFiles <- function(s) grep("^(pre_norm|post_norm|ruv_corrected|ruv_covariates|.*composite_QC_figure)", stepQc[[s]], value = TRUE)
p1obs <- if ("norm" %in% stepOrder) figureFiles("norm") else "<no step norm>"
p1want <- paste0("post_norm_", c("correlation", "density", "pca", "rle"), ".png")
firstPre <- stepOrder[vapply(stepOrder, function(s) any(grepl("^pre_norm_", stepQc[[s]])), NA)]
p2obs <- if (length(firstPre)) firstPre[[1L]] else "<never>"
comp <- unlist(stepComposite)
p3obs <- if (length(comp)) paste(names(comp), comp, sep = ": ", collapse = " | ") else "none at any step"
p4 <- obs$p4
p4obs <- if (is.null(p4)) "<not observed>" else if (!is.null(p4$error)) paste("<error>", p4$error) else
    paste0("matching states: ", one(p4$matching), "; source entity ", p4$src_id, "; states ", paste(p4$history, p4$ids, sep = "=", collapse = ", "))
copyExit <- Filter(function(e) identical(e$type, "copy_exit"), events)
p5obs <- if (length(copyExit)) paste(vapply(copyExit, function(e) paste0("[", e$failed_copies, "] (", e$returned_class, ", n=", e$n, ")"), ""), collapse = " ; ") else "<copyToResultsSummary never returned>"
p5ok <- length(copyExit) > 0L && any(vapply(copyExit, function(e) grepl("Composite QC (PDF)", e$failed_copies, fixed = TRUE) && grepl("Composite QC (PNG)", e$failed_copies, fixed = TRUE), NA))
tabEv <- Filter(function(e) identical(e$type, "tab_entry"), events)
tabRows <- if (nrow(eventTable)) eventTable[eventTable$type == "tab_entry", , drop = FALSE] else eventTable
p6obs <- if (nrow(tabRows)) paste(tabRows$step, tabRows$detail, sep = " -> ", collapse = " | ") else "<the tab entry never fired>"
p6ok <- nrow(tabRows) > 0L && any(tabRows$step == "ruv_covariates" & grepl("selected=normalization; current_state=correlation_filtered", tabRows$detail, fixed = TRUE))
p7obs <- if (is.null(obs$p7)) "<not observed>" else obs$p7$previous
verdict <- function(ok) if (isTRUE(ok)) "as predicted" else "FALSE (stop rule 1: record and report, do not interpret)"
predictions$observed <- c(p1obs <- paste(one(p1obs), "(all steps' figure files in steps.tsv)"), p2obs, p3obs, p4obs, p5obs, p6obs, p7obs)
predictions$verdict <- c(
    verdict(identical(sort(figureFiles("norm")), p1want) && "norm" %in% stepOrder),
    verdict(identical(p2obs, "ruv_covariates")),
    verdict(!length(comp)),
    verdict(!is.null(p4) && is.null(p4$error) && identical(p4$matching, "correlation_filtered")),
    verdict(p5ok),
    verdict(p6ok),
    verdict(identical(p7obs, "correlation_filtered")))
utils::write.table(predictions, file.path(out, "predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("P1 to P7:")
for (i in seq_len(nrow(predictions))) say(predictions$id[[i]], " | predicted: ", predictions$predicted[[i]], " | observed: ", predictions$observed[[i]], " | ", predictions$verdict[[i]])

# --- Per-step record against the oracle: INFORMATION ONLY --------------------------------------------------------------
if (laneOk) {
    oracle <- jsonlite::read_json(file.path(pkgroot, "tests", "testdata", "gui-entry", "prot_dia_limpa.json"), simplifyVector = FALSE)
    got <- helperFn("guiEntryComparable")(record)
    cmp <- lapply(union(names(oracle$steps), names(got)), function(id) {
        a <- got[[id]]; b <- oracle$steps[[id]]
        fs <- union(names(a), names(b))
        diffFields <- fs[!vapply(fs, function(f) identical(a[[f]], b[[f]]), logical(1))]
        ks <- union(names(a$artifacts), names(b$artifacts))
        nKeys <- sum(!vapply(ks, function(k) identical(a$artifacts[[k]], b$artifacts[[k]]), logical(1)))
        data.frame(step = id, identical = identical(a, b), differing_fields = one(diffFields), differing_artifact_keys = nKeys, stringsAsFactors = FALSE)
    })
    cmpTable <- do.call(rbind, cmp)
    utils::write.table(cmpTable, file.path(out, "oracle-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    say("per-step record against tests/testdata/gui-entry/prot_dia_limpa.json (INFORMATION ONLY):")
    print(cmpTable, row.names = FALSE, right = FALSE)
} else say("per-step oracle comparison not possible: the lane errored")

# --- Memory: two fresh processes, one figure each -----------------------------------------------------------------------
memOut <- list()
if (dir.exists(imgdir) && length(list.files(imgdir))) {
    for (mode in c("three", "two")) memOut[[mode]] <- runMem(mode)
    readings <- unlist(lapply(memOut, function(l) grep("^MEMORY_(READING|RISE)", l, value = TRUE)))
    writeUtf8(readings, file.path(out, "memory-readings.txt"))
    say("MEMORY SUMMARY:")
    cat(readings, sep = "\n")
} else say("memory not measured: no stage images were kept at the end of the lane")

writeUtf8(c(paste("traces:", paste(names(traceStatus), traceStatus, sep = "=", collapse = "; ")), unname(control),
    paste("lane ok:", laneOk)), file.path(out, "lanes.txt"))
say("g0 finished")
