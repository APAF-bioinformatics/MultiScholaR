# KIT-264 (VAL-286): the transcriptomics analysis record lists a RUVIII integration export only while the export's source
# state is a state of the head's history (compared by entity, never by name); the report's integration-export figure
# follows the record; a normalisation run removes the earlier export's figure files from `publication_graphs_dir`.
# Written before the fix (a tests-only commit): on it every E row fails at the value stated below, every L and G row
# passes, and X1 only prints. Lines `KIT-264 <row> ...` are the evidence the ticket's records copy.
#
# PREDICTIONS (written before the first run; ticket r1.2, section 10.2). Base = this commit without the fix.
#   L1  true export (run k = 2, export M = group): the record lists `trx_ruviii_export`; `ruviii-text` is the "was applied
#       ... k = 2 ... group" sentence.                                                                  base: same; head: same
#   L2  never exported: no entry; text "No RUVIII integration export was recorded for this analysis."  base: same; head: same
#   E1  route A (run k = 2, export, run k = 1): no entry; text = the "No RUVIII ..." sentence.         base: listed, "k = 2"
#   E5  route A, a PNG at QC_figures/ruviii_before_and_after.png: the `ruviii-figure` chunk prints
#       "Figure not available: ruviii_before_and_after.png " and returns no image.                      base: returns include_graphics
#   E2  route B (run k = 2, export, run with the design none, panel option on): no entry; the "No RUVIII ..." sentence;
#       the last trx_ruvs step has ruv_applied FALSE.                                                   base: listed, "k = 2"
#   E3  route C (run, export, select trx_qc_complete), record built directly: no entry; the export's source entity is in
#       `snapshot$states` but in no state of the history.                                               base: listed
#   E6  route A with ruviii_before_and_after.png and .pdf in publication_graphs_dir: removed at the second run's start, the
#       log says so; `omicNormOwnFiles` returns five names.                                      base: both remain, three names
#   E7  run k = 2, export, `resetTrxNormalisation`, run with the design none: no entry (the entity was replaced under the
#       same name `trx_ruvs`; a name test would keep it).                                               base: listed, "k = 2"
#   E8  the owner's route B through the mounted tabs (run k = 2, DA, export; design none, run, DA; Save, Copy, Report): the
#       docx says "No RUVIII integration export was recorded for this analysis." and "No correction for unwanted variation
#       was applied", not "RUVIII (ruv::RUVIII) was applied"; QC_figures/ruviii_before_and_after.png absent. Control (the
#       same journey without the second run): the RUVIII sentence and the figure are present.  base: the opposite; control same
#   G1  true export, a figure file present: entry present; text equal to L1's; the chunk returns include_graphics.   same
#   G2  never exported, no figure file: no entry; text and "Figure not available: ..." as above.                      same
#   G5  true export: the record's entry equals `workflowRecordExportAnalysis` on the head's export, field by field, and
#       the literals (operation, source_state "trx_ruvs", k = 2, M group, entity ids).                                 same
#   G6  route E (run k = 2, export, run k = 1, export): the record lists the second export (k = 1).                    same
#   G7  after route A, `snapshot$exports$trx_ruviii_export` is still present with its first `source$entity_id`.        same
#   G8  route D alone (run, export, Reset): listed (the entity is still in the history).                               same
#   X1  characterisation, never a stop: the design none by KIT-207 G5's direct call, then the export, the record and the
#       `ruviii-text` chunk, inside tryCatch; one line `KIT-264 X1 ruv_applied=... export=... entry=... k=... text=...`;
#       the test ends with succeed(). It cannot pass or fail: it is not a result of the ticket.
# R1 (this commit): E1, E2, E3, E5, E6, E7, E8 fail by assertion; L1, L2, G1, G2, G5 to G8 and E8's control pass; X1 prints.

library(testthat)

k264None <- "No RUVIII integration export was recorded for this analysis."
k264Missing <- "Figure not available: ruviii_before_and_after.png "
k264Applied <- function(k, m = "group") paste0("RUVIII (ruv::RUVIII) was applied to the TMM logCPM, with the RUVs inferred ",
    "negative controls and k = ", k, ", for the integration export only, with the replicate structure M taken from the ", m,
    ". The differential expression results do not use it.")
k264Legacy2 <- list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10)
k264Legacy1 <- list(ruv_mode = "manual", ruv_k = 1L, ruv_percentage = 10)
k264Design <- list(ruv_mode = "automatic", ruv_design = "none")
k264FigureNames <- c("ruviii_before_and_after.png", "ruviii_before_and_after.pdf")

k264Row <- function(row, ...) cat(sprintf("KIT-264 %s %s\n", row, paste(..., collapse = " ")))

#' Copy of `k257TrxRun` (test-kit257-norm-own-files.R:235-252): trx's "Run Normalization Pipeline" in process, the
#' notifications recorded in `seen$notes`.
k264Run <- function(project, nd, inputs, seen) {
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) NULL, .package = "MultiScholaR")
    notify <- function(ui, ..., type = "default", duration = 5) {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    set_status <- function(tab, value) {
        status <- project$workflow$tab_status
        status[[tab]] <- value
        project$workflow$tab_status <- status
    }
    withr::local_seed(20261124L)
    kit172Quiet(utils::capture.output(shiny::isolate(runTrxNormalisation(list(managed = runTrxNormManaged,
        in_process = runTrxNormInProcess), transcriptomicsDeclarationBinding, project$workflow$state_manager,
        trxNormGridParameters(inputs), nd, project$workflow$config_list, project$paths, omicNormAddLog(nd), notify,
        set_status))))
    invisible(NULL)
}

k264Errors <- function(seen) grep("[[]error", seen$notes, value = TRUE)

k264Start <- function() {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen
}

k264Export <- function(owner, m = "group") kit172Quiet(runTrxReferenceRuviiiExport(owner, m))

k264Snapshot <- function(owner) owner$.context$snapshot

k264Entity <- function(owner, state) k264Snapshot(owner)$states[[state]]$entity_id

#' The record as Save writes it (`writeWorkflowAnalysisRecord`, no state scope) and as the template reads it.
k264Record <- function(owner) {
    path <- tempfile(fileext = ".json")
    on.exit(unlink(path))
    kit172Quiet(writeWorkflowAnalysisRecord(owner, path))
    jsonlite::read_json(path, simplifyVector = FALSE)
}

k264Listed <- function(record) "trx_ruviii_export" %in% names(record$analyses)

k264LastRuvs <- function(record) {
    steps <- Filter(function(step) identical(step$operation, "trx_ruvs"), record$steps)
    if (length(steps)) steps[[length(steps)]]$outcomes else list()
}

#' One chunk of the head template, extracted with knitr::purl (as test-kit227-report-ruv-sentences.R's `k227Exprs`).
k264Exprs <- function(label) {
    path <- system.file("reports", "transcriptomics", "transcriptomics_report.rmd", package = "MultiScholaR")
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

#' An environment holding the head template's own `%||%`, `summary_file`, `value` and `show_figure` (the function
#' assignments of its "Project Setup" chunk) and `export` as the template's line reads it.
k264Env <- function(record, results_dir) {
    env <- new.env(parent = globalenv())
    env$report_results_summary_dir <- results_dir
    env$record <- record
    for (e in k264Exprs("Project Setup")) {
        if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
            identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
    }
    env$export <- record$analyses$trx_ruviii_export
    env
}

#' A chunk's captured output and its last returned value (an `include_graphics` value is not printed by a loop).
k264Chunk <- function(env, label) {
    result <- NULL
    out <- utils::capture.output(for (e in k264Exprs(label)) result <- eval(e, env))
    list(out = paste(out, collapse = "\n"), value = result)
}

k264Plant <- function(results_dir) {
    file <- file.path(results_dir, "QC_figures", "ruviii_before_and_after.png")
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    writeBin(charToRaw("kit264"), file)
    file
}

k264Reports <- function(fixture) list.files(fixture$paths$results_summary_dir, pattern = "^transcriptomics_report.*\\.docx$",
    full.names = TRUE)

k264Squash <- function(x) gsub("\\s+", "", gsub("‘|’|&apos;|&#39;", "'", x))

# ------------------------------------------------------------------------------------------------
# L1, G1, G5: a true export.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 L1, G1, G5: a true export (run k = 2, export M = group) is listed, described and shown", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k264Run(project, omicNormReactiveState(), k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    expect_false(is.null(k264Snapshot(owner)$exports$trx_ruviii_export), info = "export precondition")
    results <- withr::local_tempdir()
    figure <- k264Plant(results)
    record <- k264Record(owner)
    env <- k264Env(record, results)
    text <- k264Chunk(env, "ruviii-text")$out
    shown <- k264Chunk(env, "ruviii-figure")
    k264Row("L1", "listed:", k264Listed(record), "text:", text)
    expect_true(k264Listed(record), info = "L1 record")
    expect_identical(text, k264Applied(2), info = "L1 text")
    k264Row("G1", "figure file:", file.exists(figure), "chunk output:", shown$out, "class:", paste(class(shown$value), collapse = ","))
    expect_identical(text, k264Applied(2), info = "G1 text")
    expect_true(inherits(shown$value, "knit_image_paths"), info = paste("G1 figure value:", paste(class(shown$value), collapse = ",")))
    expect_identical(shown$out, "", info = "G1 figure output")
    # G5: the entry as the (unchanged) function builds it from the head's export, and the literals.
    snapshot <- k264Snapshot(owner)
    entry <- buildWorkflowAnalysisRecord(owner)$analyses$trx_ruviii_export
    direct <- workflowRecordExportAnalysis(owner$.context, snapshot, snapshot$exports$trx_ruviii_export, readWorkflowShippedDefaults())
    k264Row("G5", "operation:", entry$operation, "source_state:", entry$source_state, "k:", entry$outcomes$k, "m:",
        entry$outcomes$ruviii_m_source, "controls:", entry$outcomes$negative_control_count)
    expect_identical(entry, direct, info = "G5 entry equals workflowRecordExportAnalysis")
    expect_identical(entry$operation, "trx_ruviii_export", info = "G5 operation")
    expect_identical(entry$source_state, "trx_ruvs", info = "G5 source_state")
    expect_identical(entry$source_entity_id, k264Entity(owner, "trx_ruvs"), info = "G5 source entity")
    expect_identical(entry$product_entity_id, snapshot$exports$trx_ruviii_export$reference$entity_id, info = "G5 product entity")
    expect_equal(entry$outcomes$k, 2, info = "G5 k")
    expect_identical(entry$outcomes$ruviii_m_source, "group", info = "G5 M source")
    expect_identical(vapply(entry$parameters, function(p) p$name, ""), "ruviii_m_source", info = "G5 parameter names")
    expect_identical(entry$parameters[[1L]]$provenance, direct$parameters[[1L]]$provenance, info = "G5 provenance")
})

# ------------------------------------------------------------------------------------------------
# L2, G2: never exported.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 L2, G2: never exported, the record lists nothing and the report says so", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k264Run(project, omicNormReactiveState(), k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    results <- withr::local_tempdir()
    record <- k264Record(owner)
    env <- k264Env(record, results)
    text <- k264Chunk(env, "ruviii-text")$out
    shown <- k264Chunk(env, "ruviii-figure")
    k264Row("L2", "listed:", k264Listed(record), "text:", text)
    k264Row("G2", "figure output:", shown$out, "value:", paste(class(shown$value), collapse = ","))
    expect_false(k264Listed(record), info = "L2 record")
    expect_identical(text, k264None, info = "L2 text")
    expect_identical(shown$out, k264Missing, info = "G2 figure output")
    expect_null(shown$value, info = "G2 figure value")
})

# ------------------------------------------------------------------------------------------------
# E1, E5, G7: route A.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 E1, E5, G7: route A (run k = 2, export, run k = 1): the export is not listed, its figure is not shown, the state keeps it", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k264Run(project, nd, k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 1 precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    first <- k264Snapshot(owner)$exports$trx_ruviii_export$source$entity_id
    k264Run(project, nd, k264Legacy1, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 2 precondition:", paste(seen$notes, collapse = " | ")))
    expect_identical(k264Errors(seen), character(), info = "run errors")
    expect_false(identical(first, k264Entity(owner, "trx_ruvs")), info = "precondition: the second run is a new trx_ruvs entity")
    results <- withr::local_tempdir()
    figure <- k264Plant(results)
    record <- k264Record(owner)
    env <- k264Env(record, results)
    text <- k264Chunk(env, "ruviii-text")$out
    shown <- k264Chunk(env, "ruviii-figure")
    k264Row("E1", "listed:", k264Listed(record), "text:", text)
    expect_false(k264Listed(record), info = "E1 record")
    expect_identical(text, k264None, info = "E1 text")
    k264Row("E5", "figure file:", file.exists(figure), "chunk output:", shown$out, "class:", paste(class(shown$value), collapse = ","))
    expect_identical(shown$out, k264Missing, info = "E5 figure output")
    expect_false(inherits(shown$value, "knit_image_paths"), info = "E5 figure value")
    k264Row("G7", "export kept:", !is.null(k264Snapshot(owner)$exports$trx_ruviii_export), "same source:",
        identical(k264Snapshot(owner)$exports$trx_ruviii_export$source$entity_id, first))
    expect_false(is.null(k264Snapshot(owner)$exports$trx_ruviii_export), info = "G7 export kept")
    expect_identical(k264Snapshot(owner)$exports$trx_ruviii_export$source$entity_id, first, info = "G7 first source entity")
})

# ------------------------------------------------------------------------------------------------
# E2: route B.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 E2: route B (run k = 2, export, run with the design none): no entry, the never-exported sentence", {
    skip_if_not_installed("edgeR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k264Run(project, nd, k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 1 precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    first <- k264Snapshot(owner)$exports$trx_ruviii_export$source$entity_id
    k264Run(project, nd, k264Design, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 2 precondition:", paste(seen$notes, collapse = " | ")))
    expect_identical(k264Errors(seen), character(), info = "run errors")
    expect_false(identical(first, k264Entity(owner, "trx_ruvs")), info = "precondition: the second run is a new trx_ruvs entity")
    record <- k264Record(owner)
    expect_identical(k264LastRuvs(record)$ruv_applied, FALSE, info = "precondition: the last trx_ruvs step applied no RUV")
    text <- k264Chunk(k264Env(record, withr::local_tempdir()), "ruviii-text")$out
    k264Row("E2", "listed:", k264Listed(record), "ruv_applied:", k264LastRuvs(record)$ruv_applied, "text:", text)
    expect_false(k264Listed(record), info = "E2 record")
    expect_identical(text, k264None, info = "E2 text")
})

# ------------------------------------------------------------------------------------------------
# E3: route C (a record built directly).
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 E3: route C (run, export, Revert to trx_qc_complete): the export's source is in the states but not in the history", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    k264Run(project, omicNormReactiveState(), k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    selectWorkflowStateReference(owner, "trx_qc_complete")
    expect_identical(owner$current_state, "trx_qc_complete", info = "head precondition")
    snapshot <- k264Snapshot(owner)
    history <- as.character(unlist(snapshot$state_history))
    expect_false("trx_ruvs" %in% history, info = "precondition: trx_ruvs left the history")
    expect_identical(snapshot$exports$trx_ruviii_export$source$entity_id, snapshot$states$trx_ruvs$entity_id,
        info = "precondition: the source is still a state of the snapshot")
    record <- buildWorkflowAnalysisRecord(owner)
    k264Row("E3", "history:", paste(history, collapse = ","), "listed:", "trx_ruviii_export" %in% names(record$analyses))
    expect_false("trx_ruviii_export" %in% names(record$analyses), info = "E3 record")
})

# ------------------------------------------------------------------------------------------------
# E6: the run's start removes the earlier export's figure files.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 E6: route A with the export's figure files in publication_graphs_dir: the second run removes them; omicNormOwnFiles names five", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k264Run(project, nd, k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 1 precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    files <- file.path(project$paths$publication_graphs_dir, k264FigureNames)
    for (file in files) writeBin(charToRaw("earlier"), file)
    expect_true(all(file.exists(files)), info = "precondition: the figure files are planted")
    own <- basename(omicNormOwnFiles(transcriptomicsDeclarationBinding, project$paths))
    k264Run(project, nd, k264Legacy1, seen)
    expect_identical(k264Errors(seen), character(), info = "run errors")
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 2 precondition:", paste(seen$notes, collapse = " | ")))
    log <- paste(shiny::isolate(nd$normalization_log), collapse = "\n")
    k264Row("E6", "own files:", paste(own, collapse = ","), "present after the second run:",
        paste(k264FigureNames[file.exists(files)], collapse = ","), "log removed line:",
        grepl("Removed the earlier run's files this run replaces:", log, fixed = TRUE))
    expect_identical(own, c("canonical_correlation_panels.png", "canonical_correlation_panels.pdf",
        "canonical_correlation_optimization_summary.xlsx", k264FigureNames), info = "E6 omicNormOwnFiles")
    expect_false(any(file.exists(files)), info = "E6 files removed")
    expect_true(grepl(paste(k264FigureNames, collapse = ", "), log, fixed = TRUE), info = paste("E6 log:", log))
})

# ------------------------------------------------------------------------------------------------
# E7, G8: Reset.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 E7: run k = 2, export, Reset, run with the design none: the entity under the name trx_ruvs is replaced, no entry", {
    skip_if_not_installed("edgeR")
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k264Run(project, nd, k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 1 precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    first <- k264Snapshot(owner)$exports$trx_ruviii_export$source$entity_id
    kit172Quiet(shiny::isolate(resetTrxNormalisation(owner, nd, project$workflow$config_list, function(...) NULL,
        function(...) NULL, function(...) NULL)))
    expect_identical(owner$current_state, "trx_reset", info = "reset precondition")
    k264Run(project, nd, k264Design, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 2 precondition:", paste(seen$notes, collapse = " | ")))
    expect_identical(k264Errors(seen), character(), info = "run errors")
    snapshot <- k264Snapshot(owner)
    expect_true("trx_ruvs" %in% as.character(unlist(snapshot$state_history)), info = "precondition: the name trx_ruvs is in the history")
    expect_false(identical(first, snapshot$states$trx_ruvs$entity_id), info = "precondition: the entity under that name is the second run's")
    record <- k264Record(owner)
    expect_identical(k264LastRuvs(record)$ruv_applied, FALSE, info = "precondition: the last trx_ruvs step applied no RUV")
    text <- k264Chunk(k264Env(record, withr::local_tempdir()), "ruviii-text")$out
    k264Row("E7", "listed:", k264Listed(record), "text:", text)
    expect_false(k264Listed(record), info = "E7 record")
    expect_identical(text, k264None, info = "E7 text")
})

test_that("KIT-264 G8: run, export, Reset alone: the export is still listed (its entity is in the history)", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k264Run(project, nd, k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    kit172Quiet(shiny::isolate(resetTrxNormalisation(owner, nd, project$workflow$config_list, function(...) NULL,
        function(...) NULL, function(...) NULL)))
    expect_identical(owner$current_state, "trx_reset", info = "reset precondition")
    snapshot <- k264Snapshot(owner)
    expect_identical(snapshot$exports$trx_ruviii_export$source$entity_id, snapshot$states$trx_ruvs$entity_id,
        info = "precondition: the source is the history's trx_ruvs")
    record <- k264Record(owner)
    k264Row("G8", "head:", owner$current_state, "listed:", k264Listed(record))
    expect_true(k264Listed(record), info = "G8 record")
    expect_equal(record$analyses$trx_ruviii_export$outcomes$k, 2, info = "G8 k")
})

# ------------------------------------------------------------------------------------------------
# G6: route E, a second export after the second run.
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 G6: run k = 2, export, run k = 1, export again: the record lists the second export", {
    skip_if_not_installed("edgeR")
    seen <- k264Start()
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    nd <- omicNormReactiveState()
    k264Run(project, nd, k264Legacy2, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 1 precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    k264Run(project, nd, k264Legacy1, seen)
    expect_identical(owner$current_state, "trx_ruvs", info = paste("run 2 precondition:", paste(seen$notes, collapse = " | ")))
    k264Export(owner)
    record <- k264Record(owner)
    text <- k264Chunk(k264Env(record, withr::local_tempdir()), "ruviii-text")$out
    k264Row("G6", "listed:", k264Listed(record), "k:", record$analyses$trx_ruviii_export$outcomes$k, "text:", text)
    expect_true(k264Listed(record), info = "G6 record")
    expect_equal(record$analyses$trx_ruviii_export$outcomes$k, 1, info = "G6 k")
    expect_identical(record$analyses$trx_ruviii_export$source_entity_id, k264Entity(owner, "trx_ruvs"), info = "G6 source entity")
    expect_identical(text, k264Applied(1), info = "G6 text")
})

# ------------------------------------------------------------------------------------------------
# E8: the owner's route B through the mounted tabs, and its control.
# ------------------------------------------------------------------------------------------------

#' KIT-186 T7's fixture (test-module-ci-trx-summary.R:149-156): the MCI-030 corpus and the `project_dirs` global the
#' template needs.
k264Fixture <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE, .local_envir = env)
    fixture <- module_ci_trx_project(module_ci_trx_corpus(), env = env)
    assign("project_dirs", fixture$project_dirs, envir = globalenv())
    withr::defer(if (exists("project_dirs", envir = globalenv())) rm("project_dirs", envir = globalenv()), envir = env)
    fixture
}

k264Docx <- function(fixture) {
    report <- k264Reports(fixture)
    expect_length(report, 1L)
    if (length(report)) module_ci_trx_docx_text(report[[1L]]) else ""
}

test_that("KIT-264 E8: the owner's route B through the mounted tabs: the report does not state the earlier export", {
    fixture <- k264Fixture()
    owner <- fixture$owner
    again <- function(session, output, out) {
        session$setInputs(`norm-ruv_mode` = "automatic", `norm-ruv_design` = "none")
        session$setInputs(`norm-run_normalization` = 2L)
        session$flushReact()
        session$setInputs(`de-run_da_analysis` = 2L)
        session$flushReact()
        module_ci_trx_summary_steps(session, output, out)
    }
    out <- module_ci_trx_tabs_journey(fixture, norm = list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10),
        ruviii = "group", summary = again)
    expect_identical(out$errors, character(), info = "precondition: no step failed")
    expect_false(is.null(k264Snapshot(owner)$exports$trx_ruviii_export), info = "precondition: the export is in the state")
    record <- buildWorkflowAnalysisRecord(owner)
    expect_identical(k264LastRuvs(record)$ruv_applied, FALSE, info = "precondition: the last trx_ruvs step applied no RUV")
    text <- k264Docx(fixture)
    figure <- file.path(fixture$paths$results_summary_dir, "QC_figures", "ruviii_before_and_after.png")
    has <- function(statement) grepl(k264Squash(statement), k264Squash(text), fixed = TRUE)
    k264Row("E8", "none-sentence:", has(k264None), "no-correction:", has("No correction for unwanted variation was applied"),
        "applied-sentence:", has("RUVIII (ruv::RUVIII) was applied"), "figure in the results summary:", file.exists(figure))
    expect_true(has(k264None), info = "E8 never-exported sentence")
    expect_true(has("No correction for unwanted variation was applied"), info = "E8 no-correction caption")
    expect_false(has("RUVIII (ruv::RUVIII) was applied"), info = "E8 the earlier export is not stated")
    expect_false(file.exists(figure), info = "E8 the earlier figure is not delivered")
})

test_that("KIT-264 E8 control: the same journey without the second run keeps the RUVIII sentence and the figure", {
    fixture <- k264Fixture()
    owner <- fixture$owner
    out <- module_ci_trx_tabs_journey(fixture, norm = list(ruv_mode = "manual", ruv_k = 2L, ruv_percentage = 10),
        ruviii = "group", summary = module_ci_trx_summary_steps)
    expect_identical(out$errors, character(), info = "precondition: no step failed")
    expect_false(is.null(k264Snapshot(owner)$exports$trx_ruviii_export), info = "precondition: the export is in the state")
    record <- buildWorkflowAnalysisRecord(owner)
    expect_identical(k264LastRuvs(record)$ruv_applied, TRUE, info = "precondition: the last trx_ruvs step applied RUV")
    text <- k264Docx(fixture)
    figure <- file.path(fixture$paths$results_summary_dir, "QC_figures", "ruviii_before_and_after.png")
    applied <- grepl(k264Squash("RUVIII (ruv::RUVIII) was applied"), k264Squash(text), fixed = TRUE)
    k264Row("E8 control", "applied-sentence:", applied, "figure in the results summary:", file.exists(figure))
    expect_true(applied, info = "E8 control: the RUVIII sentence")
    expect_true(file.exists(figure), info = "E8 control: the figure")
})

# ------------------------------------------------------------------------------------------------
# X1: characterisation (never a stop, never a pass or a failure of the ticket).
# ------------------------------------------------------------------------------------------------

test_that("KIT-264 X1 (characterisation only): an export made after a run without RUV is observed, not judged", {
    project <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- project$owner
    observed <- tryCatch({
        # KIT-207 G5's form (test-kit207-summary-head-guard.R:752): the state without RUV, no option needed.
        kit172Quiet(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none")))
        export <- tryCatch({
            k264Export(owner)
            "ok"
        }, error = function(e) paste0(paste(class(e), collapse = "/"), ": ", conditionMessage(e)))
        record <- k264Record(owner)
        env <- k264Env(record, withr::local_tempdir())
        text <- k264Chunk(env, "ruviii-text")$out
        sprintf("ruv_applied=%s export=%s entry=%s k=%s text=%s", format(k264LastRuvs(record)$ruv_applied), export,
            if (k264Listed(record)) "listed" else "absent", format(record$analyses$trx_ruviii_export$outcomes$k %||% NA),
            gsub("\n", " ", text, fixed = TRUE))
    }, error = function(e) paste0("observation failed: ", paste(class(e), collapse = "/"), ": ", conditionMessage(e)))
    cat(sprintf("KIT-264 X1 %s\n", observed))
    testthat::succeed()
})
