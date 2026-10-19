# KIT-039 render check on the kit-journey-028 member's project (reflinked). The record and study_parameters.txt
# are written by this build (the KIT-038 content), then:
#   DIA-NN report: the pre-ticket template (the KIT-038 commit's, the baseline) and this ticket's template with
#                  the record and with the record removed; compared as text (dates masked). With the record the
#                  Version Information section carries the record's version, source commit, build and R version
#                  and no "NA"; without it the section is the baseline's git fields.
#   RUV sentence:  the report's RUV text matches the record's outcomes$ruv_applied, k and percentage.
#   SI:            renderWorkflowAnalysisRecord to DOCX and HTML; every step's description sentence, every
#                  parameter's label and the legend present; the workbook's parameters sheet has a label column.
#   study:         study_parameters.txt identical to KIT-038's output modulo the Timestamp line.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-039/report-check.R <baseline template> <label>
suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE); baseline_template <- normalizePath(args[[1L]]); label <- args[[2L]]
repo <- normalizePath(".")
h <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
journey <- readRDS(file.path(b$root, "kit-journey-028", "candidate-1", "journey-result.rds"))
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit039-report-", label))
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
json <- file.path(paths$source_dir, "analysis_record.json")
writeWorkflowAnalysisRecord(manager, json)
record <- jsonlite::read_json(json, simplifyVector = FALSE)
contrasts <- do.call(rbind, lapply(record$analyses$differential_abundance$contrasts, function(row) as.data.frame(lapply(row, unlist))))
invisible(capture.output(createWorkflowArgsFromConfig("dia_baseline_smoke", source_dir_path = paths$source_dir,
    final_s4_object = getDIAStateReference(manager), contrasts_tbl = contrasts,
    workflow_data = list(state_manager = manager, config_list = list()))))
# The same writer on the same store, as KIT-038 ran it: the file must be identical modulo the Timestamp line.
study <- readLines(file.path(paths$source_dir, "study_parameters.txt"))
kit038_dir <- file.path(work, "kit038"); dir.create(kit038_dir)
invisible(capture.output(createWorkflowArgsFromConfig("dia_baseline_smoke", source_dir_path = kit038_dir,
    final_s4_object = getDIAStateReference(manager), contrasts_tbl = contrasts,
    workflow_data = list(state_manager = manager, config_list = list()))))
study_again <- readLines(file.path(kit038_dir, "study_parameters.txt"))
closeDIAStateAdapter(manager)
strip <- function(x) x[!grepl("^Timestamp:", x)]
out <- list(schema = "kit.report_check.v1", ticket = "KIT-039", label = label, date = format(Sys.Date()),
    pandoc = as.character(rmarkdown::pandoc_version()),
    study_parameters = list(lines = length(study), identical_modulo_timestamp = identical(strip(study), strip(study_again))))
docx_text <- function(path, name) {
    xml <- utils::unzip(path, files = "word/document.xml", exdir = file.path(work, name))
    text <- gsub("<[^>]+>", " ", paste(readLines(xml, warn = FALSE), collapse = ""))
    gsub("[[:space:]]+", " ", text)
}
render <- function(template, name) {
    target <- file.path(base, "scripts", "proteomics", "DIANN_report.rmd")
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    file.copy(template, target, overwrite = TRUE)
    t0 <- Sys.time()
    path <- tryCatch(RenderReport(omic_type = "proteomics", experiment_label = "dia_baseline_smoke", rmd_filename = "DIANN_report.rmd"),
        error = function(e) e)
    if (inherits(path, "error")) return(list(error = conditionMessage(path)))
    keep <- file.path(work, paste0(name, ".docx")); file.copy(path, keep, overwrite = TRUE)
    list(path = keep, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), text = docx_text(keep, name))
}
current <- file.path(repo, "inst/reports/proteomics/DIANN_report.rmd")
with_record <- render(current, "with-record")
# The baseline and the no-record render both run without the record, so only the template differs between them.
file.rename(json, file.path(work, "analysis_record.json.removed"))
baseline <- render(baseline_template, "baseline")
without_record <- render(current, "without-record")
file.rename(file.path(work, "analysis_record.json.removed"), json)
mask <- function(text) gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}( [0-9]{2}:[0-9]{2}:[0-9]{2})?", "<time>", text)
section <- function(text) regmatches(text, regexpr("Version Information.{0,400}", text, perl = TRUE))
for (name in c("baseline", "with_record", "without_record")) {
    out[[name]] <- get(name)[c("error", "seconds")]
    if (!is.null(get(name)$text)) writeLines(get(name)$text, file.path(repo, "docs/qa/workflow-execution/kit/kit-039", paste0("DIANN_report-", gsub("_", "-", name), "-text.txt")))
}
if (is.null(baseline$error) && is.null(without_record$error) && is.null(with_record$error)) {
    out$version_with_record <- section(with_record$text)
    out$version_without_record <- section(without_record$text)
    out$version_baseline <- section(baseline$text)
    out$with_record_names_software <- grepl(record$software$version, out$version_with_record, fixed = TRUE) &&
        grepl(substr(record$software$sha, 1, 12), out$version_with_record, fixed = TRUE) &&
        grepl(substr(record$software$build, 1, 12), out$version_with_record, fixed = TRUE) &&
        !grepl("\\bNA\\b", out$version_with_record)
    out$without_record_equals_baseline <- identical(mask(baseline$text), mask(without_record$text))
    if (!out$without_record_equals_baseline) {
        x <- strsplit(mask(baseline$text), "")[[1L]]; y <- strsplit(mask(without_record$text), "")[[1L]]
        n <- min(length(x), length(y)); first <- which(x[seq_len(n)] != y[seq_len(n)])
        at <- if (length(first)) first[[1L]] else n + 1L
        out$first_difference <- list(position = at, baseline = substr(mask(baseline$text), max(1L, at - 120L), at + 160L),
            without_record = substr(mask(without_record$text), max(1L, at - 120L), at + 160L))
    }
    ruv <- Filter(function(step) identical(step$operation, "normalization-ruv"), record$steps)[[1L]]$outcomes
    sentence <- regmatches(with_record$text, regexpr("(RUV-III batch correction was not applied[^.]*[.]|identified that using [^.]*[.] Based on this, [0-9]+ unwanted components)", with_record$text, perl = TRUE))
    out$ruv <- list(outcomes = ruv, sentence = sentence,
        matches = if (isTRUE(ruv$ruv_applied)) grepl(sprintf("using %s%% of features", ruv$negative_control_percentage), sentence, fixed = TRUE) &&
            grepl(sprintf("(%s proteins)", ruv$negative_control_count), sentence, fixed = TRUE) &&
            grepl(sprintf("%s unwanted components", ruv$k), sentence, fixed = TRUE) else grepl("not applied", sentence))
    file.copy(with_record$path, file.path(repo, "docs/qa/workflow-execution/kit/kit-039", "DIANN_report-with-record.docx"), overwrite = TRUE)
}
# The SI from the record: DOCX and HTML, the workbook.
si_dir <- file.path(work, "si"); dir.create(si_dir)
docx <- renderWorkflowAnalysisRecord(json, si_dir, "docx")
html <- renderWorkflowAnalysisRecord(json, si_dir, "html")
si_text <- docx_text(docx, "si-docx")
tables <- workflowAnalysisRecordTables(record)
descriptions <- Filter(function(d) !is.na(d) && nzchar(d), tables$steps$description)
labels <- unique(tables$parameters$label[tables$parameters$provenance != "unknown"])
out$si <- list(docx = basename(docx), html = basename(html),
    descriptions = length(descriptions), descriptions_present = all(vapply(descriptions, function(d) grepl(d, si_text, fixed = TRUE), logical(1))),
    labels = length(labels), labels_present = all(vapply(labels, function(l) grepl(l, si_text, fixed = TRUE), logical(1))),
    legend_present = grepl("Provenance legend", si_text, fixed = TRUE) && grepl("no_shipped_default", si_text, fixed = TRUE),
    outcome_lines = grepl("Outcome: rows", si_text, fixed = TRUE))
workbook <- file.path(si_dir, "analysis_record.xlsx")
writeWorkflowAnalysisRecordWorkbook(record, workbook)
sheet <- openxlsx::read.xlsx(workbook, sheet = "parameters")
out$workbook <- list(parameters_rows = nrow(sheet), label_column = "label" %in% names(sheet),
    empty_labels = sum(is.na(sheet$label) | !nzchar(sheet$label)))
file.copy(docx, file.path(repo, "docs/qa/workflow-execution/kit/kit-039", "analysis_record_SI-screen.docx"), overwrite = TRUE)
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-039", paste0("report-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("RENDER baseline", baseline$error %||% "ok", "| with record", with_record$error %||% "ok", "| without", without_record$error %||% "ok", "\n")
cat("VERSION with record:", out$version_with_record %||% "", "\n  names software:", out$with_record_names_software %||% NA, "| without = baseline:", out$without_record_equals_baseline %||% NA, "\n")
cat("RUV sentence:", out$ruv$sentence %||% "", "| matches:", out$ruv$matches %||% NA, "\n")
cat("SI descriptions", out$si$descriptions, "present", out$si$descriptions_present, "| labels", out$si$labels, "present", out$si$labels_present, "| legend", out$si$legend_present, "| outcomes", out$si$outcome_lines, "\n")
cat("WORKBOOK rows", out$workbook$parameters_rows, "label column", out$workbook$label_column, "empty labels", out$workbook$empty_labels, "\n")
cat("STUDY lines", out$study_parameters$lines, "identical modulo timestamp", out$study_parameters$identical_modulo_timestamp, "| pandoc", out$pandoc, "\n")
