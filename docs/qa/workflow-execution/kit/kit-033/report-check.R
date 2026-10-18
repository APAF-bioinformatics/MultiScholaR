# KIT-033 report check on a real journey (the kit-journey-028 member's project, reflinked; kit_027 build):
# study_parameters.txt and analysis_record.json are written from the store (the record path), then the DIA-NN
# report is rendered three times through RenderReport: the pre-ticket template (baseline), this ticket's template
# with the record, and this ticket's template with the record removed. The documents are compared as text.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-033/report-check.R <baseline template> <label>
suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE); baseline_template <- normalizePath(args[[1L]]); label <- args[[2L]]
repo <- normalizePath(".")
h <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
journey <- readRDS(file.path(b$root, "kit-journey-028", "candidate-1", "journey-result.rds"))
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit033-report-", label))
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(journey$project, "."), work)) == 0L)
store_root <- file.path(work, substring(dirname(journey$store), nchar(journey$project) + 2L))
base <- file.path(work, "dia_baseline_smoke")
project_dirs <- setupDirectories(base_dir = base, omic_types = "proteomics", label = "dia_baseline_smoke", force = TRUE, reuse_existing = TRUE)
assign("project_dirs", project_dirs, envir = globalenv())
paths <- project_dirs[[1L]]
identity <- readSavedStoreIdentity(file.path(store_root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
manager$.context$reference_states <- TRUE
writeWorkflowAnalysisRecord(manager, file.path(paths$source_dir, "analysis_record.json"))
record <- jsonlite::read_json(file.path(paths$source_dir, "analysis_record.json"), simplifyVector = FALSE)
contrasts <- do.call(rbind, lapply(record$analyses$differential_abundance$contrasts, function(row) as.data.frame(lapply(row, unlist))))
invisible(capture.output(createWorkflowArgsFromConfig("dia_baseline_smoke", source_dir_path = paths$source_dir,
    final_s4_object = getDIAStateReference(manager), contrasts_tbl = contrasts,
    workflow_data = list(state_manager = manager, config_list = list()))))
closeDIAStateAdapter(manager)
study <- readLines(file.path(paths$source_dir, "study_parameters.txt"))
out <- list(schema = "kit.report_check.v1", ticket = "KIT-033", label = label, date = format(Sys.Date()),
    pandoc = as.character(rmarkdown::pandoc_version()),
    study_parameters = list(lines = length(study), truncation_marks = sum(grepl("omitted|\\[list object\\]|\\.\\.\\.", study)),
        analysis_record_section = "Analysis Record:" %in% study))
render <- function(template, name) {
    target <- file.path(base, "scripts", "proteomics", "DIANN_report.rmd")
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    file.copy(template, target, overwrite = TRUE)
    t0 <- Sys.time()
    path <- tryCatch(RenderReport(omic_type = "proteomics", experiment_label = "dia_baseline_smoke", rmd_filename = "DIANN_report.rmd"),
        error = function(e) e)
    if (inherits(path, "error")) return(list(error = conditionMessage(path)))
    keep <- file.path(work, paste0(name, ".docx")); file.copy(path, keep, overwrite = TRUE)
    xml <- utils::unzip(keep, files = "word/document.xml", exdir = file.path(work, name))
    text <- gsub("<[^>]+>", " ", paste(readLines(xml, warn = FALSE), collapse = ""))
    list(path = keep, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), text = gsub("[[:space:]]+", " ", text))
}
current <- file.path(repo, "inst/reports/proteomics/DIANN_report.rmd")
baseline <- render(baseline_template, "baseline")
with_record <- render(current, "with-record")
file.rename(file.path(paths$source_dir, "analysis_record.json"), file.path(work, "analysis_record.json.removed"))
without_record <- render(current, "without-record")
# A study_parameters.txt without the sample filter's section (as a session whose file lacks it): with the record
# the report states the recorded cutoff; without it the template default (the defect the ticket names).
study_path <- file.path(paths$source_dir, "study_parameters.txt")
start <- match("[filterMinNumPeptidesPerSample]", trimws(study))
end <- start + match(TRUE, !nzchar(trimws(study[-seq_len(start)])))
writeLines(study[-seq(start, end)], study_path)
missing_without <- render(current, "missing-without-record")
file.rename(file.path(work, "analysis_record.json.removed"), file.path(paths$source_dir, "analysis_record.json"))
missing_with <- render(current, "missing-with-record")
writeLines(study, study_path)
sentence <- function(text) regmatches(text, regexpr("Samples with fewer than [0-9]+ peptides", text))
out$missing_section <- list(removed = "[filterMinNumPeptidesPerSample]", without_record = sentence(missing_without$text %||% ""),
    with_record = sentence(missing_with$text %||% ""),
    recorded_cutoff = decodeWorkflowRecordValue(Filter(function(p) identical(p$name, "peptides_per_sample_cutoff"),
        Filter(function(step) identical(step$operation, "minimum_peptides_per_sample"), record$steps)[[1L]]$parameters)[[1L]]$value),
    errors = c(missing_without$error, missing_with$error))
mask <- function(text) gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}( [0-9]{2}:[0-9]{2}:[0-9]{2})?", "<time>", text)
for (name in c("baseline", "with_record", "without_record")) out[[name]] <- get(name)[c("error", "seconds")]
if (is.null(baseline$error) && is.null(without_record$error) && is.null(with_record$error)) {
    out$without_record_equals_baseline <- identical(mask(baseline$text), mask(without_record$text))
    out$with_record_differs <- !identical(mask(baseline$text), mask(with_record$text))
    out$with_record_note <- regmatches(with_record$text, regexpr("Parameter not recorded; default shown:[^.]*\\.", with_record$text))
    file.copy(with_record$path, file.path(repo, "docs/qa/workflow-execution/kit/kit-033", "DIANN_report-with-record.docx"), overwrite = TRUE)
}
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-033", paste0("report-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("STUDY lines", out$study_parameters$lines, "truncation marks", out$study_parameters$truncation_marks, "record section", out$study_parameters$analysis_record_section, "\n")
cat("RENDER baseline", baseline$error %||% "ok", "| with record", with_record$error %||% "ok", "| without", without_record$error %||% "ok", "\n")
cat("MISSING without record:", out$missing_section$without_record, "\n  with record:", out$missing_section$with_record, "\n")
cat("COMPARE without = baseline:", out$without_record_equals_baseline %||% NA, "| with differs:", out$with_record_differs %||% NA, "| note:", out$with_record_note %||% "", "\n")
