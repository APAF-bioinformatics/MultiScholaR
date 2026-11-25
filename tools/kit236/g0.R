# KIT-236 gate G0 (ticket r1.2 section 4): observe each false sentence of section 1 on the real GUI lanes, with KIT-235's
# parser in, and each replacement's source. No package change. Nothing is asserted: every outcome is recorded beside the
# ticket's prediction; section 4 of the ticket decides from the record.
#
# Usage (run-tests scripts job, from the package root):  Rscript tools/kit236/g0.R
# Writes kit236-g0/:
#   meta.txt                      preconditions (pandoc, locale, versions, commits)
#   g0a/<lane>-report-text.txt    the lane's report (docx text: </w:p> to a line break, other tags removed)
#   g0a/<lane>-sentences.tsv      section 1's sentences: present or absent, the ticket's prediction, the paragraph
#   g0a/<lane>-record.txt         the record fields of section 1.4 (decoded), the DA model verbatim, n and k
#   g0a/<lane>-steps.tsv          each record step: dimensions, outcomes, parameters
#   g0a/<lane>-file.txt           study_parameters.txt facts (ITSD lines, DA UI block, Version Information)
#   g0a/<lane>-parsed.txt         the template's own parser, extractor and AnalyzeSamples (record overlay) on the files
#   g0a/<lane>-files.tsv          every non-store file of the lane's project (path, bytes)
#   g0a/<lane>-src-*              the lane's scripts/<omic>/ study_parameters.txt, analysis_record.json, design_matrix.tab
#   g0b/summary.tsv               per lane and variant V0 (unchanged copy) to V7: render outcome and each check
#   g0b/<lane>-<V>-text.txt, g0b/<lane>-<V>-parsed.txt, g0b/<lane>-<V>-edits.txt
#   g0c/summary.tsv, g0c/<case>.txt   chunk-level cases (lipid, metab, the four proteomics templates); plotting scale

out <- normalizePath("kit236-g0", mustWork = FALSE)
for (d in c("g0a", "g0b", "g0c")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
work <- file.path(tempdir(), "kit236-work")
dir.create(work, showWarnings = FALSE)

say <- function(...) cat(..., "\n", sep = "")
writeUtf8 <- function(lines, path) writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
toJson <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "string", force = TRUE,
    digits = NA))
strText <- function(x) paste(utils::capture.output(utils::str(x, give.attr = FALSE, vec.len = 6, nchar.max = 400)),
    collapse = "\n")
`%||%` <- function(a, b) if (!is.null(a)) a else b
tsv <- function(df, path) utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
oneLine <- function(x) gsub("[\t\r\n]+", " ", paste(as.character(x), collapse = " "))

# --- Preconditions -------------------------------------------------------------------------------------------------
gitOut <- function(args) tryCatch(system2("git", args, stdout = TRUE, stderr = TRUE), error = function(e) conditionMessage(e))
meta <- c(
    paste("R", R.version.string),
    paste("pandoc_available", rmarkdown::pandoc_available()),
    paste("pandoc_version", tryCatch(as.character(rmarkdown::pandoc_version()), error = function(e) conditionMessage(e))),
    paste("locale", Sys.getlocale()),
    paste("rmarkdown", as.character(utils::packageVersion("rmarkdown"))),
    paste("knitr", as.character(utils::packageVersion("knitr"))),
    paste("HEAD", paste(gitOut(c("rev-parse", "HEAD")), collapse = " ")),
    paste("log", paste(gitOut(c("log", "--oneline", "-6")), collapse = " | ")),
    paste("templates md5", paste(names(tools::md5sum(Sys.glob("inst/reports/*/*.rmd"))),
        unname(tools::md5sum(Sys.glob("inst/reports/*/*.rmd"))), collapse = "; ")))
writeUtf8(meta, file.path(out, "meta.txt"))
say(paste(meta, collapse = "\n"))

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
ns <- asNamespace("MultiScholaR")
say("loaded from the tree: ", as.character(utils::packageVersion("MultiScholaR")))
encodeValue <- get("encodeWorkflowRecordValue", envir = ns)
decodeValue <- get("decodeWorkflowRecordValue", envir = ns)
overlay <- get("workflowRecordReportParameters", envir = ns)

# --- Text helpers ----------------------------------------------------------------------------------------------------
normText <- function(x) {
    x <- paste(x, collapse = "\n")
    for (q in c("‘", "’")) x <- gsub(q, "'", x, fixed = TRUE)
    for (q in c("“", "”")) x <- gsub(q, "\"", x, fixed = TRUE)
    x <- gsub(" ", " ", x, fixed = TRUE)
    trimws(gsub("[[:space:]]+", " ", x))
}
docxXml <- function(path) {
    dir <- tempfile("docx-")
    utils::unzip(path, files = "word/document.xml", exdir = dir)
    paste(readLines(file.path(dir, "word", "document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
}
xmlText <- function(xml) {
    text <- gsub("</w:p>", "\n", xml, fixed = TRUE)
    text <- gsub("<[^>]+>", "", text)
    for (pair in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&apos;", "'"), c("&amp;", "&"))) {
        text <- gsub(pair[[1L]], pair[[2L]], text, fixed = TRUE)
    }
    text
}
countDrawings <- function(xml) sum(gregexpr("<w:drawing>", xml, fixed = TRUE)[[1L]] > 0L)
has <- function(norm, phrase) grepl(normText(phrase), norm, fixed = TRUE)
contextOf <- function(paragraphs, norm, phrase) {
    p <- normText(phrase)
    hit <- which(vapply(paragraphs, function(x) grepl(p, x, fixed = TRUE), logical(1)))
    x <- if (length(hit)) paragraphs[[hit[[1L]]]] else norm
    pos <- regexpr(p, x, fixed = TRUE)
    if (pos < 0) return("")
    if (nchar(x) > 700L) x <- substr(x, max(1L, pos - 250L), pos + 450L)
    x
}
readReport <- function(path) {
    xml <- docxXml(path)
    raw <- xmlText(xml)
    paragraphs <- vapply(strsplit(raw, "\n", fixed = TRUE)[[1L]], normText, "", USE.NAMES = FALSE)
    list(raw = raw, norm = normText(raw), paragraphs = paragraphs[nzchar(paragraphs)], drawings = countDrawings(xml))
}

# --- Section 1's sentences: fixed phrases (whitespace-normalised), and the ticket's G0a prediction per lane ----------
# Columns: id, lipid phrase, metab phrase (NA: not in that template), prediction lipid_canonical, prediction metab_lc,
# requires (id whose phrase must also be present: S8/S14 share their words and are told apart by their caption).
catalogue <- local({
    r <- function(id, lip, met = lip, pl, pm = pl, requires = "") data.frame(id = id, lip = lip, met = met,
        pred_lip = pl, pred_met = pm, requires = requires, stringsAsFactors = FALSE)
    rbind(
        r("S1", "The PCA plots (Figures 1a-c) demonstrate the effectiveness of batch effect removal", pl = "absent"),
        r("S2", "Initial log2 data showed notable batch effects", pl = "absent"),
        r("S3", "After cyclic loess normalization, batch effects remained in PC1 but were reduced in PC2", pl = "absent"),
        r("S4", "batches merged effectively, suggesting successful removal of unwanted variations", pl = "absent"),
        r("S5", "the distribution of the PC1 and PC2 values for each group are similar.", pl = "absent"),
        r("S6", "Successive rounds of normalization have reduced the variability in the PC1 and PC2 values", pl = "absent"),
        r("S7", "The RLE plots (Figures 1g-i) show progressive improvement in data quality", pl = "absent"),
        r("S8", "Smaller interquartile ranges (IQR) after normalization indicate reduced technical variation", pl = "absent",
            requires = "S7"),
        r("S9", "Final normalized data shows consistent median values near zero", pl = "absent"),
        r("S10", "biological samples was high and consistent across all normalization methods", pl = "absent"),
        r("S11", "After cyclic loess normalization, technical variation was reduced", pl = "present"),
        r("S12", "values for each group are similar after normalization.", pl = "present"),
        r("S13", "The RLE plots (Figures 1e-f) show improvement in data quality after normalization", pl = "present"),
        r("S14", "Smaller interquartile ranges (IQR) after normalization indicate reduced technical variation",
            pl = "present", requires = "S13"),
        r("S15", "Normalized data shows more consistent median values near zero", pl = "present"),
        r("S16", "biological samples was high and consistent across normalization methods", pl = "present"),
        r("S17", "All samples passed initial QC criteria", pl = "present"),
        r("S18", "Normalization successfully reduced technical variation", pl = "present"),
        r("S19", "Batch effects were substantially mitigated, though not completely eliminated", pl = "present"),
        r("S20", "High technical reproducibility achieved for majority of lipids",
            "High technical reproducibility achieved for majority of metabolites", pl = "present"),
        r("S21", "Correlation metrics indicate reliable quantification", pl = "present"),
        r("S22", "Sample processing met quality standards", pl = "present"),
        r("S23", "Residual batch effects present but within acceptable ranges", pl = "present"),
        r("S24", "Technical variation adequately controlled through normalization steps", pl = "present"),
        r("S25-lip", "The lipidomics data quality is sufficient for differential analysis.", NA, pl = "present", pm = "n/a"),
        r("S26-lip", "have effectively reduced technical noise while preserving biological signal", NA, pl = "present",
            pm = "n/a"),
        r("S27-lip", "Proceed with differential abundance analysis using the fully normalized dataset", NA,
            pl = "present", pm = "n/a"),
        r("S28-lip", "Consider RUV-adjusted data for downstream multivariate analyses", NA, pl = "present", pm = "n/a"),
        r("S29-lip", "Validate top significant lipids using targeted assays if possible", NA, pl = "present", pm = "n/a"),
        r("S25-met", NA, "The dataset demonstrates robust technical quality", pl = "n/a", pm = "present"),
        r("S26-met", NA, "Normalization strategy effectively reduced systematic biases", pl = "n/a", pm = "present"),
        r("S27-met", NA, "Proceed with downstream analysis using normalized dataset", pl = "n/a", pm = "present"),
        r("S28-met", NA, "Quality metrics support the reliability of the data", pl = "n/a", pm = "present"),
        r("S29-met", NA, "Technical reproducibility meets industry standards", pl = "n/a", pm = "present"),
        r("S30-met", NA, "Dataset is suitable for biological interpretation", pl = "n/a", pm = "present"),
        r("S30/S31", "biological samples with no technical replicates.", pl = "present"),
        r("M1", "normalization steps before integration for differential abundance analysis", pl = "present"),
        r("M2", "Missing values (zeros or non-detects) were coded as NA to distinguish true absence from missing measurements",
            pl = "present"),
        r("M3", "Between-sample normalization was performed using the 'cyclicloess' method from the 'limma' R package",
            pl = "present"),
        r("M4a", "Quality control metrics were evaluated separately for each assay to ensure data quality before integration",
            pl = "present"),
        r("M4b", "with cross-assay consistency verified in the combined analysis", pl = "present"),
        r("M5", "Figure 1 shows the QC metrics across three stages of data processing", pl = "absent"),
        r("M6", "Figure 1 shows the QC metrics across two stages of data processing", pl = "present"),
        r("M6-loess", "Cyclic loess normalization (rightmost column", pl = "present"),
        r("M7", "using the adjusted abundance matrix for comparing each pair of experimental groups", pl = "present"),
        r("M8", "A linear model for comparing each pair of conditions was fitted using the formula", pl = "present"),
        r("M9", "The false discovery rate correction was applied to the moderated p-values by calculating the q-values (Storey, 2002).",
            pl = "present"),
        r("M10", "Data analyses were performed based on SOP xxx.", pl = "present"),
        r("M11", "Comments on test specific conditions information (e.g. temperature), where relevant", pl = "present"),
        r("M12", "Sample identification (eg client #, sample name, Salesforce ID etc)", pl = "present"),
        r("M13", "Internal standard (ITSD) normalization was applied using the average_centered method with median aggregation",
            pl = "present"),
        r("M14", "The standard empirical Bayes method was used to test for any deviation from zero log-fold-change.",
            pl = "present"),
        r("M15", "Significant differentially abundant lipids were defined as those with q-values less than",
            "Significant differentially abundant metabolites were defined as those with q-values less than", pl = "present"),
        r("F1", "Volcano plots of differentially abundant lipids across all groups are shown in Figure 3",
            "Volcano plots of differentially abundant metabolites across all groups are shown in Figure 3", pl = "present"),
        r("F2-warning", "Volcano Plots directory not found", pl = "present"),
        r("F2-caption", "Figure 3: Volcano plots of differentially abundant lipids in alphabetical order.",
            "Figure 3: Volcano plots of differentially abundant metabolites in alphabetical order.", pl = "present"),
        r("F3", "DE_lipids_results.xlsx", "DE_metabolites_results.xlsx", pl = "present"),
        r("F4-warning", "NumSigDeMolecules directory not found", pl = "absent", pm = "present"),
        r("F4-caption", "Figure 2: Number of significant differentially abundant lipids for each comparison.",
            "Figure 2: Number of significant differentially abundant metabolites for each comparison.", pl = "present"),
        r("F5-warning", "Composite QC figure not found", pl = "absent"))
})

sentenceTable <- function(report, omic_short) {
    phrase_col <- if (identical(omic_short, "lip")) "lip" else "met"
    pred_col <- if (identical(omic_short, "lip")) "pred_lip" else "pred_met"
    phrase_of <- stats::setNames(catalogue[[phrase_col]], catalogue$id)
    rows <- lapply(seq_len(nrow(catalogue)), function(i) {
        phrase <- catalogue[[phrase_col]][[i]]
        if (is.na(phrase)) return(NULL)
        present <- has(report$norm, phrase)
        req <- catalogue$requires[[i]]
        if (nzchar(req)) present <- present && has(report$norm, phrase_of[[req]])
        observed <- if (present) "present" else "absent"
        predicted <- catalogue[[pred_col]][[i]]
        data.frame(id = catalogue$id[[i]], observed = observed, predicted = predicted,
            verdict = if (identical(observed, predicted)) "as predicted" else "DIFFERS", phrase = phrase,
            context = if (present) contextOf(report$paragraphs, report$norm, phrase) else "", stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# --- Template code: chunks read from the .rmd (no knitr), function definitions evaluated in a namespace child --------
templateChunks <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    starts <- grep("^```\\{r", lines)
    lapply(starts, function(s) {
        e <- s + which(grepl("^```[[:space:]]*$", lines[-seq_len(s)]))[[1L]]
        list(label = trimws(sub("^```\\{r[[:space:],]*([^,}]*).*$", "\\1", lines[[s]])), header = lines[[s]],
            code = if (e > s + 1L) lines[(s + 1L):(e - 1L)] else character())
    })
}
isFunctionDef <- function(e) is.call(e) && (identical(e[[1L]], as.name("<-")) || identical(e[[1L]], as.name("="))) &&
    is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))
templateEnv <- function(path, omic_type, source_dir, qc_dir = file.path(source_dir, "qc"),
    figures_dir = file.path(source_dir, "figures"), tables_dir = file.path(source_dir, "tables")) {
    env <- new.env(parent = ns)
    errors <- character()
    for (chunk in templateChunks(path)) {
        exprs <- tryCatch(parse(text = chunk$code, keep.source = FALSE), error = function(e) {
            errors <<- c(errors, paste0(chunk$label, ": parse: ", conditionMessage(e))); expression()
        })
        for (e in exprs) if (isFunctionDef(e)) tryCatch(eval(e, env), error = function(err) {
            errors <<- c(errors, paste0(chunk$label, ": ", conditionMessage(err)))
        })
    }
    env$`%||%` <- function(a, b) if (!is.null(a)) a else b
    env$d66_log <- function(...) invisible(NULL)
    env$params <- list(omic_type = omic_type, experiment_label = "kit236", workflow_name = "kit236",
        timestamp = "2026-11-25 00:00:00")
    env$report_source_dir <- source_dir
    env$report_results_summary_dir <- source_dir
    env$report_qc_figures_dir <- qc_dir
    env$report_publication_figures_dir <- figures_dir
    env$report_publication_tables_dir <- tables_dir
    env$.def_errors <- errors
    env
}
runChunk <- function(path, label, env) {
    chunk <- Filter(function(ch) identical(ch$label, label), templateChunks(path))
    if (!length(chunk)) return(paste("<no chunk", label, ">"))
    res <- tryCatch(utils::capture.output(for (e in parse(text = chunk[[1L]]$code, keep.source = FALSE)) {
        v <- withVisible(eval(e, env))
        if (v$visible) print(v$value)
    }), error = function(err) paste("<chunk error>", conditionMessage(err)))
    c(paste("header:", chunk[[1L]]$header), res)
}
fnName <- function(env, pattern) {
    hit <- sort(grep(pattern, ls(env), value = TRUE))
    if (length(hit)) hit[[1L]] else NA_character_
}
extractorOf <- function(env) fnName(env, "^Extract.*(StudyParameters|Parameters)")

fieldNames <- c("num_assays", "assay_names", "itsd_applied", "itsd_method_type", "itsd_aggregation", "itsd_removed",
    "itsd_features_per_assay", "log_transform_applied", "log_transform_offset", "normalisation_method", "ruv_applied",
    "s4_de_q_val_thresh", "s4_de_treat_lfc_cutoff", "ui_de_q_value_threshold", "ui_de_log_fold_change_cutoff",
    "ui_de_treat_enabled", "s4_formula_string", "da_block_reason", "da_model", "software_version", "build",
    "source_commit", "repository", "branch", "commit", "git_timestamp", "technical_replicates", "biological_samples",
    "total_samples", "contrasts", "record_present", "record_defaults_note")
pick <- function(x) stats::setNames(lapply(fieldNames, function(f) if (f %in% names(x)) x[[f]] else "<absent>"), fieldNames)

# The template's own parser and extractor on a source directory, then AnalyzeSamples (counts + extractor + overlay).
parseCase <- function(template, omic_type, source_dir) {
    env <- templateEnv(template, omic_type, source_dir)
    sp <- file.path(source_dir, "study_parameters.txt")
    study <- tryCatch(env$ReadStudyParamsNew(sp), error = function(e) structure(list(), parse_error = conditionMessage(e)))
    extractor <- extractorOf(env)
    extracted <- tryCatch(env[[extractor]](study), error = function(e) structure(list(), extract_error = conditionMessage(e)))
    design <- tryCatch(readr::read_tsv(file.path(source_dir, "design_matrix.tab"), show_col_types = FALSE),
        error = function(e) NULL)
    analyzed <- tryCatch(env$AnalyzeSamples(design, study), error = function(e) list(numbers = list(),
        error = conditionMessage(e)))
    lines <- c(
        paste("definition errors:", if (length(env$.def_errors)) paste(env$.def_errors, collapse = " | ") else "none"),
        paste("parse error:", attr(study, "parse_error") %||% "none"),
        paste("sections:", paste0("'", names(study), "'", collapse = ", ")),
        "== section 'ITSD Normalization' as parsed", strText(study[["ITSD Normalization"]]),
        "== section 'Log Transformation' as parsed", strText(study[["Log Transformation"]]),
        "== section 'Between-Sample Normalization' as parsed", strText(study[["Between-Sample Normalization"]]),
        "== section 'Version Information' as parsed", strText(study[["Version Information"]]),
        "== 'Workflow Parameters' names", strText(names(study[["Workflow Parameters"]])),
        "== 'Parameters from Final S4 Object' $deAnalysisParameters",
        strText(study[["Workflow Parameters"]][["Parameters from Final S4 Object"]][["deAnalysisParameters"]]),
        "== 'User Interface Parameters' as parsed",
        strText(study[["Workflow Parameters"]][["User Interface Parameters"]]),
        paste("extractor:", extractor, "; error:", attr(extracted, "extract_error") %||% "none", "; length:",
            length(extracted)),
        "== extractor fields (no record overlay)", toJson(pick(extracted)),
        paste("AnalyzeSamples numbers length:", length(analyzed$numbers)),
        "== AnalyzeSamples fields (record overlay when the record is present)", toJson(pick(analyzed$numbers)))
    list(env = env, study = study, extracted = extracted, numbers = analyzed$numbers, lines = lines)
}

# --- G0a: the real lanes ---------------------------------------------------------------------------------------------
helperFn <- function(name) {
    if (exists(name, envir = ns, inherits = FALSE)) get(name, envir = ns) else get(name, envir = globalenv())
}
runLane <- helperFn("guiEntryRunLane")
store_pattern <- "(^|/)\\.multischolar(/|$)|(^|/)\\.workflow|\\.duckdb|\\.qs2$|\\.wal$"
cap <- new.env(parent = emptyenv())
cap$lane <- NULL
cap$n <- 0L
cap$base <- list()
cap$dirs <- list()
kit236Capture <- function(base_dir) {
    lane <- cap$lane
    if (is.null(lane)) return(invisible(NULL))
    cap$n <- cap$n + 1L
    dest <- file.path(work, lane, "base")
    unlink(dest, recursive = TRUE)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    dirs <- list.dirs(base_dir, recursive = TRUE, full.names = FALSE)
    for (d in dirs[nzchar(dirs) & !grepl(store_pattern, dirs)]) dir.create(file.path(dest, d), recursive = TRUE,
        showWarnings = FALSE)
    files <- list.files(base_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    for (f in files[!grepl(store_pattern, files)]) file.copy(file.path(base_dir, f), file.path(dest, f), overwrite = TRUE)
    cap$base[[lane]] <- base_dir
    cap$dirs[[lane]] <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
    invisible(NULL)
}
trace_ok <- tryCatch({
    # As KIT-235's G0 instrument fix 1: the capture function embedded in the tracer call as an object.
    trace("guiEntryArtifacts", tracer = bquote(.(kit236Capture)(base_dir)), where = runLane, print = FALSE)
    TRUE
}, error = function(e) { say("trace failed: ", conditionMessage(e)); FALSE })
if (!trace_ok) {
    helper_env <- environment(runLane)
    original_artifacts <- get("guiEntryArtifacts", envir = helper_env)
    if (bindingIsLocked("guiEntryArtifacts", helper_env)) unlockBinding("guiEntryArtifacts", helper_env)
    assign("guiEntryArtifacts", function(base_dir, roots) {
        kit236Capture(base_dir)
        original_artifacts(base_dir, roots)
    }, envir = helper_env)
}
say("trace installed: ", trace_ok)

lanes <- list(
    lip = list(id = "lipid_canonical", omic = "lipidomics", prefix = "lipid", word = "lipids",
        template = "lipidomics_report.rmd"),
    met = list(id = "metab_lc", omic = "metabolomics", prefix = "metab", word = "metabolites",
        template = "metabolomics_report.rmd"))

relTo <- function(path, base) substring(path, nchar(base) + 2L)
mapPaths <- function(pd, old, new) rapply(pd, function(s) ifelse(startsWith(s, old), paste0(new, substring(s, nchar(old) + 1L)),
    s), classes = "character", how = "replace")
lanePaths <- function(short, new_base) {
    lane <- lanes[[short]]
    pd <- mapPaths(cap$dirs[[short]], cap$base[[short]], new_base)
    list(all = pd, paths = pd[[lane$omic]])
}
recordParam <- function(rec, operation, name) {
    steps <- Filter(function(s) identical(s$operation, operation), rec$steps)
    if (!length(steps)) return("<no step>")
    hit <- Filter(function(p) identical(p$name, name), steps[[length(steps)]]$parameters)
    if (!length(hit)) return("<no parameter>")
    decodeValue(hit[[1L]]$value)
}
sectionLines <- function(lines, heading) {
    i <- which(lines == heading)
    if (!length(i)) return("<section absent>")
    j <- i[[1L]] + 1L
    while (j <= length(lines) && nzchar(trimws(lines[[j]]))) j <- j + 1L
    lines[i[[1L]]:(j - 1L)]
}

lane_status <- list()
for (short in names(lanes)) {
    lane <- lanes[[short]]
    cap$lane <- short
    cap$n <- 0L
    started <- Sys.time()
    status <- tryCatch({ runLane(lane$id); "ran" }, error = function(e) paste("error:", conditionMessage(e)))
    lane_status[[short]] <- sprintf("%s %s: %s (%.0f s, %d captures)", short, lane$id, status,
        as.numeric(difftime(Sys.time(), started, units = "secs")), cap$n)
    say(lane_status[[short]])
}
cap$lane <- NULL
writeUtf8(c(paste("trace installed:", trace_ok), unlist(lane_status)), file.path(out, "g0a", "lanes.txt"))

laneInfo <- list()
for (short in names(lanes)) {
    lane <- lanes[[short]]
    tryCatch({
        base <- file.path(work, short, "base")
        if (is.null(cap$dirs[[short]])) stop("no capture for this lane")
        lp <- lanePaths(short, base)
        paths <- lp$paths
        src <- paths$source_dir
        files <- list.files(base, recursive = TRUE, all.files = TRUE, no.. = TRUE)
        tsv(data.frame(path = files, bytes = file.size(file.path(base, files))), file.path(out, "g0a", paste0(short, "-files.tsv")))
        all_paths <- unlist(rapply(cap$dirs[[short]], function(s) s, classes = "character", how = "unlist"))
        outside <- all_paths[!is.na(all_paths) & !startsWith(all_paths, cap$base[[short]])]
        for (f in c("study_parameters.txt", "analysis_record.json", "design_matrix.tab")) {
            if (file.exists(file.path(src, f))) file.copy(file.path(src, f), file.path(out, "g0a", paste0(short, "-src-", f)))
        }
        docx <- list.files(paths$results_summary_dir, pattern = "\\.docx$", full.names = TRUE)
        docx <- docx[!grepl("analysis_record", basename(docx), fixed = TRUE)]
        info <- list(base = base, paths = paths, docx = docx, outside = outside)
        if (length(docx)) {
            report <- readReport(docx[[1L]])
            info$report <- report
            writeUtf8(report$raw, file.path(out, "g0a", paste0(short, "-report-text.txt")))
            tab <- sentenceTable(report, short)
            tsv(tab, file.path(out, "g0a", paste0(short, "-sentences.tsv")))
            say(short, ": ", sum(tab$verdict == "DIFFERS"), " of ", nrow(tab), " sentences differ from the prediction")
        }
        # The record (section 1.4 and R1-N04, R1-N06).
        rec_path <- file.path(src, "analysis_record.json")
        rl <- c(paste("report docx:", paste(basename(docx), collapse = ", ")),
            paste("project paths outside base_dir:", paste(outside, collapse = ", ")))
        if (file.exists(rec_path)) {
            rec <- jsonlite::read_json(rec_path, simplifyVector = FALSE)
            da <- rec$analyses$differential_abundance
            da_params <- stats::setNames(lapply(da$parameters, function(p) decodeValue(p$value)),
                vapply(da$parameters, function(p) p$name, ""))
            steps <- rec$steps
            tsv(do.call(rbind, lapply(steps, function(s) data.frame(order = s$order %||% NA, state = s$state %||% NA,
                operation = s$operation %||% NA, dimensions = toJson(s$dimensions), outcomes = toJson(s$outcomes),
                parameters = toJson(stats::setNames(lapply(s$parameters, function(p) decodeValue(p$value)),
                    vapply(s$parameters, function(p) p$name, ""))), stringsAsFactors = FALSE))),
                file.path(out, "g0a", paste0(short, "-steps.tsv")))
            design_rows <- rec$design$rows
            trid <- rec$design$technical_replicate_id %||% "<absent>"
            trvals <- vapply(design_rows, function(row) as.character(row[[trid]] %||% NA), "")
            design_file <- tryCatch(utils::read.delim(file.path(src, "design_matrix.tab"), stringsAsFactors = FALSE,
                check.names = FALSE), error = function(e) NULL)
            pearson_op <- paste0(lane$prefix, "_pearson_selection")
            pearson <- Filter(function(s) identical(s$operation, pearson_op), steps)
            rl <- c(rl,
                paste("schema:", rec$schema %||% "<absent>"),
                paste("normalisation method (record):", toJson(recordParam(rec, paste0(lane$prefix,
                    "_between_sample_normalisation"), "normalisation_method"))),
                paste("log2 offset (record):", toJson(recordParam(rec, paste0(lane$prefix, "_log2"), "offset"))),
                paste("ruv step present:", any(vapply(steps, function(s) identical(s$operation, paste0(lane$prefix, "_ruv")),
                    logical(1)))),
                paste("itsd step outcomes:", toJson(Filter(Negate(is.null), lapply(steps, function(s)
                    if (identical(s$operation, paste0(lane$prefix, "_itsd_normalisation"))) s$outcomes else NULL)))),
                paste("DA operation:", da$operation %||% "<absent>"),
                paste("DA parameters (decoded):", toJson(da_params)),
                paste("DA contrasts:", toJson(da$contrasts)),
                paste("DA model (verbatim):", toJson(da$model)),
                paste("DA results:", toJson(da$results)),
                paste("DA failures:", toJson(da$failures)),
                paste("software version:", toJson(rec$software$version), "; build:", toJson(rec$software$build),
                    "; sha:", toJson(rec$software$sha)),
                paste("design: sample_id =", toJson(rec$design$sample_id), "; group_id =", toJson(rec$design$group_id),
                    "; technical_replicate_id =", toJson(trid)),
                paste("record design rows (k):", length(design_rows), "; distinct technical_replicate_id values:",
                    length(unique(trvals)), "; values:", paste(trvals, collapse = ",")),
                paste("design_matrix.tab rows (n):", if (is.null(design_file)) "<unreadable>" else nrow(design_file),
                    "; columns:", if (is.null(design_file)) "" else paste(names(design_file), collapse = ",")),
                paste("pearson step:", pearson_op, "present =", length(pearson) > 0L, "; threshold =",
                    toJson(recordParam(rec, pearson_op, "min_pearson_correlation_threshold"))),
                paste("n of k: n =", if (is.null(design_file)) NA else nrow(design_file), "; k =", length(design_rows)),
                "== step dimensions",
                vapply(steps, function(s) sprintf("%s %s %s: dimensions %s", s$order %||% "", s$state %||% "",
                    s$operation %||% "", toJson(s$dimensions)), ""))
        } else {
            rl <- c(rl, "analysis_record.json absent in scripts/<omic>/")
        }
        writeUtf8(rl, file.path(out, "g0a", paste0(short, "-record.txt")))
        # The file's facts.
        sp_path <- file.path(src, "study_parameters.txt")
        if (file.exists(sp_path)) {
            sp <- readLines(sp_path, warn = FALSE, encoding = "UTF-8")
            fl <- c("== ITSD Normalization section", sectionLines(sp, "ITSD Normalization:"),
                paste("'ITSD Features Per Assay' block present:", any(grepl("ITSD Features Per Assay", sp, fixed = TRUE))),
                paste("'[Differential Abundance UI Parameters]' present:", any(grepl("Differential Abundance UI Parameters", sp,
                    fixed = TRUE))),
                paste("'[Differential Expression UI Parameters]' present:", any(grepl("Differential Expression UI Parameters", sp,
                    fixed = TRUE))),
                paste("'User Interface Parameters' present:", any(grepl("User Interface Parameters", sp, fixed = TRUE))),
                "== Version Information section", sectionLines(sp, "Version Information:"),
                "== Log Transformation section", sectionLines(sp, "Log Transformation:"),
                "== Between-Sample Normalization section", sectionLines(sp, "Between-Sample Normalization:"),
                "== deAnalysisParameters lines", sp[grepl("deAnalysisParameters|treat_lfc_cutoff|formula_string|lfc \\[", sp)])
            writeUtf8(fl, file.path(out, "g0a", paste0(short, "-file.txt")))
            parsed <- parseCase(file.path(src, lane$template), lane$omic, src)
            writeUtf8(parsed$lines, file.path(out, "g0a", paste0(short, "-parsed.txt")))
        }
        laneInfo[[short]] <- info
    }, error = function(e) say("G0a ", short, ": harness error: ", conditionMessage(e)))
}

# --- G0b: planted variants through the real render (RenderReport as the Summary calls it) ----------------------------
replaceGroup <- function(text, pattern, value) {
    m <- gregexpr(pattern, text, perl = TRUE)[[1L]]
    if (m[[1L]] == -1L) return(list(text = text, n = 0L))
    st <- attr(m, "capture.start")[, "v"]
    ln <- attr(m, "capture.length")[, "v"]
    for (i in rev(seq_along(st))) text <- paste0(substr(text, 1L, st[[i]] - 1L), value, substr(text, st[[i]] + ln[[i]],
        nchar(text)))
    list(text = text, n = length(st))
}
editRecord <- function(path, pattern, value, label) {
    text <- readChar(path, file.size(path), useBytes = TRUE)
    res <- replaceGroup(text, pattern, value)
    writeLines(res$text, path, sep = "", useBytes = TRUE)
    sprintf("record %s: %d replacement(s) with %s", label, res$n, value)
}
editFile <- function(path, fn, label) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    new <- fn(lines)
    writeUtf8(new, path)
    sprintf("file %s: %d line(s) changed or added", label, sum(!new %in% lines))
}
pngAt <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(path, width = 200, height = 200)
    graphics::plot.new()
    graphics::text(0.5, 0.5, "KIT-236 G0 placeholder")
    grDevices::dev.off()
    sprintf("image written: %s (%s)", basename(path), file.exists(path))
}

variants <- list(
    V0 = list(what = "unchanged copy (control: the re-render reproduces the lane's report)", edit = function(p, short) "none",
        checks = list()),
    V1 = list(what = "normalisation 'none' in the record and the file", edit = function(p, short) c(
        editRecord(file.path(p$source_dir, "analysis_record.json"),
            '"name": "normalisation_method",\\s*"source": "[^"]*",\\s*"value": \\{\\s*"type": "character",\\s*"value": (?<v>"[^"]*")',
            '"none"', "normalisation_method"),
        editFile(file.path(p$source_dir, "study_parameters.txt"), function(l) sub("^\\* Method: .*$", "* Method: none", l),
            "* Method: none"),
        paste("encodeWorkflowRecordValue('none') =", toJson(encodeValue("none")))),
        checks = list(c("'none' method from the 'limma' R package", "present"), c("Cyclic loess normalization (rightmost column",
            "present"))),
    V2 = list(what = "log offset 1 in the record and the file", edit = function(p, short) c(
        editRecord(file.path(p$source_dir, "analysis_record.json"),
            '"name": "offset",\\s*"source": "[^"]*",\\s*"value": \\{\\s*"type": "double",\\s*"value": (?<v>[-0-9.eE]+)', "1",
            "offset"),
        editFile(file.path(p$source_dir, "study_parameters.txt"), function(l) sub("^\\* Transform: log2$", "* Offset: 1", l),
            "* Offset: 1")),
        checks = list(c("transformed as log2(x + 1)", "present"), c("were coded as NA", "present"))),
    V3 = list(what = "record design rows: two runs of one group share a technical_replicate_id value", edit = function(p, short) {
        rec_path <- file.path(p$source_dir, "analysis_record.json")
        rec <- jsonlite::read_json(rec_path, simplifyVector = FALSE)
        trid <- rec$design$technical_replicate_id
        rows <- rec$design$rows
        groups <- vapply(rows, function(r) as.character(r[[rec$design$group_id]]), "")
        g <- groups[[1L]]
        idx <- which(groups == g)
        first <- as.character(rows[[idx[[1L]]]][[trid]])
        second <- as.character(rows[[idx[[2L]]]][[trid]])
        c(editRecord(rec_path, sprintf('"%s": (?<v>"%s")', trid, second), sprintf('"%s"', first), paste(trid, second)),
            sprintf("%s of run 2 of group %s set from %s to %s", trid, g, second, first))
    }, checks = list(c("with no technical replicates", "present"))),
    V4 = list(what = "an image in Publication_figures/Volcano_Plots and NumSigDeMolecules", edit = function(p, short) c(
        pngAt(file.path(p$publication_figures_dir, "Volcano_Plots", "kit236_g0_placeholder.png")),
        pngAt(file.path(p$publication_figures_dir, "NumSigDeMolecules", "kit236_g0_placeholder.png"))),
        checks = list(c("across all groups are shown in Figure 3", "present"), c("Figure 3: Volcano plots", "present"),
            c("Figure 2: Number of significant", "present"), c("Volcano Plots directory not found", "absent"),
            c("NumSigDeMolecules directory not found", "absent"), c("<more drawings than V0>", "present"))),
    V5 = list(what = "DA_results_<omic>.xlsx removed", edit = function(p, short) {
        f <- file.path(p$publication_tables_dir, paste0("DA_results_", lanes[[short]]$omic, ".xlsx"))
        existed <- file.exists(f)
        unlink(f)
        sprintf("%s existed = %s; exists now = %s", basename(f), existed, file.exists(f))
    }, checks = list(c("_results.xlsx within the Publication_tables folder", "present"))),
    V6 = list(what = "the file's ITSD block lists one assay", edit = function(p, short) {
        sp <- file.path(p$source_dir, "study_parameters.txt")
        lines <- readLines(sp, warn = FALSE, encoding = "UTF-8")
        assays <- trimws(strsplit(sub("^\\* Assay Names: ", "", grep("^\\* Assay Names: ", lines, value = TRUE)[[1L]]),
            ",")[[1L]])
        at <- grep("^\\* ITSD Removed After Normalization:", lines)[[1L]]
        if (at < length(lines) && grepl("^\\* Timestamp:", lines[[at + 1L]])) at <- at + 1L
        block <- c("* ITSD Features Per Assay:", sprintf("    - %s (2 features): IS_kit236_a, IS_kit236_b", assays[[1L]]))
        writeUtf8(append(lines, block, after = at), sp)
        c(sprintf("block inserted after line %d: %s", at, paste(block, collapse = " / ")),
            paste("assays:", paste(assays, collapse = ", ")))
    }, checks = list(c("was applied using the average_centered method", "present"),
        c("The ITSD features used per assay were:", "present"), c("<ITSD paragraph names the second assay>", "absent"))),
    V7 = list(what = "record DA lfc 0.5 and the file's treat_lfc_cutoff 0.5", edit = function(p, short) c(
        editRecord(file.path(p$source_dir, "analysis_record.json"),
            '"name": "lfc",\\s*"source": "[^"]*",\\s*"value": \\{\\s*"type": "double",\\s*"value": (?<v>[-0-9.eE]+)', "0.5", "lfc"),
        editFile(file.path(p$source_dir, "study_parameters.txt"), function(l) {
            at <- grep("^\\[deAnalysisParameters\\]$", l)
            l <- sub("^(  lfc \\[[a-z_]+\\] = )0$", "\\10.5", l)
            if (length(at)) append(l, "  treat_lfc_cutoff = 0.5", after = at[[1L]]) else l
        }, "treat_lfc_cutoff = 0.5")),
        checks = list(c("The standard empirical Bayes method was used to test for any deviation from zero log-fold-change",
            "present"), c("absolute log2 fold change", "absent"), c("TREAT", "absent"))))

g0b_rows <- list()
v0_drawings <- list()
for (short in names(laneInfo)) {
    lane <- lanes[[short]]
    for (vid in names(variants)) {
        v <- variants[[vid]]
        row <- tryCatch({
            root <- file.path(tempdir(), "kit236-g0b", short, vid)
            unlink(root, recursive = TRUE)
            dir.create(root, recursive = TRUE)
            file.copy(laneInfo[[short]]$base, root, recursive = TRUE)
            new_base <- file.path(root, "base")
            lp <- lanePaths(short, new_base)
            p <- lp$paths
            for (old in list.files(p$results_summary_dir, pattern = "\\.docx$", full.names = TRUE)) {
                if (!grepl("analysis_record", basename(old), fixed = TRUE)) unlink(old)
            }
            edits <- v$edit(p, short)
            writeUtf8(c(v$what, edits), file.path(out, "g0b", paste0(short, "-", vid, "-edits.txt")))
            assign("project_dirs", lp$all, envir = globalenv())
            started <- Sys.time()
            rendered <- tryCatch(RenderReport(omic_type = lane$omic, experiment_label = lane$id, rmd_filename = lane$template),
                error = function(e) structure(NA_character_, render_error = conditionMessage(e)))
            secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
            parsed <- parseCase(file.path(p$source_dir, lane$template), lane$omic, p$source_dir)
            writeUtf8(parsed$lines, file.path(out, "g0b", paste0(short, "-", vid, "-parsed.txt")))
            ok <- !is.null(rendered) && !is.na(rendered) && file.exists(rendered)
            check_text <- "render failed"
            if (ok) {
                report <- readReport(rendered)
                writeUtf8(report$raw, file.path(out, "g0b", paste0(short, "-", vid, "-text.txt")))
                if (identical(vid, "V0")) {
                    v0_drawings[[short]] <- report$drawings
                    tab <- sentenceTable(report, short)
                    tsv(tab, file.path(out, "g0b", paste0(short, "-V0-sentences.tsv")))
                    check_text <- sprintf("V0 sentences: %d of %d as the G0a prediction", sum(tab$verdict == "as predicted"),
                        nrow(tab))
                }
                results <- vapply(v$checks, function(ch) {
                    observed <- if (identical(ch[[1L]], "<more drawings than V0>")) {
                        report$drawings > (v0_drawings[[short]] %||% Inf)
                    } else if (identical(ch[[1L]], "<ITSD paragraph names the second assay>")) {
                        assays <- trimws(strsplit(as.character(parsed$numbers$assay_names %||% ""), ",")[[1L]])
                        ctx <- contextOf(report$paragraphs, report$norm, "Internal standard (ITSD) normalization was applied")
                        length(assays) >= 2L && grepl(assays[[2L]], ctx, fixed = TRUE)
                    } else has(report$norm, ch[[1L]])
                    obs <- if (isTRUE(observed)) "present" else "absent"
                    sprintf("[%s] '%s': %s (predicted %s)", if (identical(obs, ch[[2L]])) "as predicted" else "DIFFERS",
                        ch[[1L]], obs, ch[[2L]])
                }, "")
                if (length(results)) check_text <- paste(results, collapse = " | ")
                contexts <- vapply(c("normalization was performed using", "Missing values", "technical replicates",
                    "Figure 3", "results table", "Internal standard (ITSD)", "empirical Bayes method"), function(k) {
                    paste0(k, " => ", contextOf(report$paragraphs, report$norm, k))
                }, "")
                writeUtf8(c(sprintf("drawings: %d", report$drawings), contexts),
                    file.path(out, "g0b", paste0(short, "-", vid, "-contexts.txt")))
            }
            data.frame(lane = short, variant = vid, rendered = ok, seconds = round(secs),
                render_error = attr(rendered, "render_error") %||% "", drawings = if (ok) report$drawings else NA,
                checks = check_text, edits = oneLine(edits), stringsAsFactors = FALSE)
        }, error = function(e) data.frame(lane = short, variant = vid, rendered = FALSE, seconds = NA,
            render_error = paste("harness:", conditionMessage(e)), drawings = NA, checks = "", edits = "",
            stringsAsFactors = FALSE))
        say("G0b ", short, " ", vid, ": rendered = ", row$rendered, "; ", row$checks)
        g0b_rows[[length(g0b_rows) + 1L]] <- row
    }
}
if (length(g0b_rows)) tsv(do.call(rbind, g0b_rows), file.path(out, "g0b", "summary.tsv"))

# --- G0c: chunk-level cases ------------------------------------------------------------------------------------------
g0c_rows <- list()
g0cRecord <- function(case, template, check, observed, predicted, text) {
    obs <- if (isTRUE(observed)) "present" else "absent"
    g0c_rows[[length(g0c_rows) + 1L]] <<- data.frame(case = case, template = basename(template), check = check,
        observed = obs, predicted = predicted, verdict = if (identical(obs, predicted)) "as predicted" else
            if (identical(predicted, "record")) "recorded" else "DIFFERS", stringsAsFactors = FALSE)
    writeUtf8(c(paste("template:", template), paste("check:", check), paste("observed:", obs), "== text", text),
        file.path(out, "g0c", paste0(case, ".txt")))
}
plantRecord <- function(dir, lfc = NULL, block_reason = NULL, treat_applied = NULL, contrasts = NULL, fdr = NULL,
    steps = list()) {
    da <- list(operation = "kit236_planted", parameters = list(), contrasts = contrasts %||% list(),
        model = list(model = "limma"), results = list())
    if (!is.null(lfc)) da$parameters <- list(list(name = "lfc", source = "manifest", value = encodeValue(lfc),
        provenance = "declared"))
    if (!is.null(block_reason)) da$model$block_reason <- block_reason
    if (!is.null(treat_applied)) da$model$treat_applied <- treat_applied
    if (!is.null(fdr)) da$results <- list(list(contrast = "groupA-groupB", assay = "assay_1", entity_id = "e1", rows = 10L,
        fdr_method = fdr, fdr_fallback_reason = "planted"))
    rec <- list(schema = "multischolar.analysis_record.v1", software = list(package = "MultiScholaR", version = "0.9.0",
        build = "b1"), design = list(rows = list()), steps = steps, analyses = list(differential_abundance = da))
    jsonlite::write_json(rec, file.path(dir, "analysis_record.json"), auto_unbox = TRUE, null = "null", digits = NA,
        pretty = TRUE)
    invisible(file.path(dir, "analysis_record.json"))
}
caseDir <- function(name) {
    d <- file.path(tempdir(), "kit236-g0c", name)
    unlink(d, recursive = TRUE)
    dir.create(d, recursive = TRUE)
    d
}
numbersFor <- function(env, design, study = list()) {
    res <- tryCatch(env$AnalyzeSamples(design, study), error = function(e) list(numbers = list()))
    res$numbers
}
design2 <- data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
design3 <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B", "C"), each = 2), stringsAsFactors = FALSE)
designTr <- data.frame(Run = c("S1", "S2", "S3", "S4", "S5", "S6"), group = c("A", "A", "A", "B", "B", "B"),
    replicates = c(1L, 1L, 2L, 1L, 1L, 2L), tech_rep_group = c("A_1", "A_1", "A_2", "B_1", "B_1", "B_2"),
    stringsAsFactors = FALSE)

omic_templates <- list(lipid = list(path = "inst/reports/lipidomics/lipidomics_report.rmd", omic = "lipidomics"),
    metab = list(path = "inst/reports/metabolomics/metabolomics_report.rmd", omic = "metabolomics"))
for (o in names(omic_templates)) {
    tp <- omic_templates[[o]]
    tryCatch({
        # RUV applied: S1 to S10 and M5.
        env <- templateEnv(tp$path, tp$omic, caseDir(paste0(o, "-ruv")))
        txt <- tryCatch(env$FormatResults(NULL, ruv_applied = TRUE), error = function(e) paste("<error>", conditionMessage(e)))
        norm <- normText(txt)
        for (id in c(paste0("S", 1:10), "M5")) {
            g0cRecord(paste0(o, "-ruv-", id), tp$path, catalogue$lip[catalogue$id == id], has(norm,
                catalogue$lip[catalogue$id == id]), "present", txt)
        }
        # Three groups, one comparison in the record's contrasts: M7, M8.
        d <- caseDir(paste0(o, "-3groups"))
        plantRecord(d, contrasts = list(list(contrasts = "groupA-groupB", friendly_names = "A_vs_B",
            full_format = "A_vs_B=groupA-groupB")))
        env <- templateEnv(tp$path, tp$omic, d)
        numbers <- numbersFor(env, design3)
        txt <- tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
        g0cRecord(paste0(o, "-3groups-M7"), tp$path, "comparing each pair of experimental groups",
            has(normText(txt), "comparing each pair of experimental groups"), "present",
            c(txt, paste("numbers length:", length(numbers)), paste("contrasts:", toJson(numbers$contrasts))))
        # block_reason through the overlay: the single and the joined string.
        for (br in c("technical_replicates", "no_duplicates; technical_replicates")) {
            d <- caseDir(paste0(o, "-block-", gsub("[^a-z]", "", br)))
            plantRecord(d, block_reason = br)
            env <- templateEnv(tp$path, tp$omic, d)
            numbers <- numbersFor(env, design2)
            txt <- tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
            g0cRecord(paste0(o, "-block-", if (grepl(";", br)) "joined" else "single"), tp$path, "duplicateCorrelation",
                has(normText(txt), "duplicateCorrelation"), "absent",
                c(txt, paste("da_block_reason through the overlay:", toJson(numbers$da_block_reason))))
        }
        # The composite absent: the QC_metrics chunk.
        d <- caseDir(paste0(o, "-nocomposite"))
        env <- templateEnv(tp$path, tp$omic, d, qc_dir = file.path(d, "empty_qc"))
        dir.create(file.path(d, "empty_qc"))
        chunk_out <- runChunk(tp$path, "QC_metrics", env)
        g0cRecord(paste0(o, "-nocomposite-F5"), tp$path, "Composite QC figure not found",
            any(grepl("Composite QC figure not found", chunk_out, fixed = TRUE)), "present", chunk_out)
    }, error = function(e) say("G0c ", o, ": harness error: ", conditionMessage(e)))
}

# The plotting scale of the composite's left column (post_filter data): does plotPca/plotRle/plotDensity/plotPearson
# transform the data? (1) the methods' and their helpers' source lines naming log; (2) plotRle and plotPca on a fixture
# whose intensities double across samples, compared with the raw-scale and log2-scale values computed here.
scale_lines <- character()
fx <- data.frame(database_identifier = c("F1", "F2", "F3", "F4"), annotation = c("A1", "A2", "A3", "A4"),
    S1 = c(100, 300, 1000, 50), S2 = c(200, 500, 1500, 80), S3 = c(400, 700, 4000, 160), S4 = c(800, 1100, 9000, 400),
    check.names = FALSE, stringsAsFactors = FALSE)
fx_design <- data.frame(Sample_ID = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
fx_mat <- as.matrix(fx[, c("S1", "S2", "S3", "S4")])
rle_raw <- sort(as.vector(fx_mat - apply(fx_mat, 1L, stats::median)))
rle_log <- sort(as.vector(log2(fx_mat) - apply(log2(fx_mat), 1L, stats::median)))
scale_lines <- c(scale_lines, paste("expected RLE values, raw scale:", paste(signif(rle_raw, 5), collapse = ",")),
    paste("expected RLE values, log2 scale:", paste(signif(rle_log, 5), collapse = ",")))
makeFixture <- function(o) {
    args <- list(ITSDNormalization = list(applied = FALSE), log_transformed = FALSE, normalisation_method = "none")
    if (identical(o, "lipid")) {
        d <- fx
        names(d)[[2L]] <- "lipid_identification"
        createLipidomicsAssayData(lipid_data = list(assay_1 = d), design_matrix = fx_design,
            lipid_id_column = "database_identifier", annotation_id_column = "lipid_identification", sample_id = "Sample_ID",
            group_id = "group", database_identifier_type = "LipidMaps", internal_standard_regex = "^LIPID_ISTD", args = args)
    } else {
        d <- fx
        names(d)[[2L]] <- "metabolite_identification"
        createMetaboliteAssayData(metabolite_data = list(assay_1 = d), design_matrix = fx_design,
            metabolite_id_column = "database_identifier", annotation_id_column = "metabolite_identification",
            sample_id = "Sample_ID", group_id = "group", database_identifier_type = "HMDB", internal_standard_regex = "^ISTD",
            args = args)
    }
}
matchScale <- function(values, raw, logv) {
    v <- sort(as.numeric(values[is.finite(values)]))
    if (length(v) == length(raw) && isTRUE(all.equal(v, raw, tolerance = 1e-6))) return("raw scale")
    if (length(v) == length(logv) && isTRUE(all.equal(v, logv, tolerance = 1e-6))) return("log2 scale")
    "neither"
}
for (o in c("lipid", "metab")) {
    cls <- if (identical(o, "lipid")) "LipidomicsAssayData" else "MetaboliteAssayData"
    for (f in c("plotPca", "plotRle", "plotDensity", "plotPearson")) {
        src <- tryCatch(deparse(methods::getMethod(f, cls)), error = function(e) paste("<no method>", conditionMessage(e)))
        helpers <- unique(unlist(regmatches(src, gregexpr("[A-Za-z0-9_.]+Helper[A-Za-z0-9_]*(?=\\()", src, perl = TRUE))))
        log_lines <- grep("log2|log\\(|log1p|log10", src, value = TRUE)
        scale_lines <- c(scale_lines, sprintf("%s %s: %d source lines; helpers called: %s; lines naming log: %s", o, f,
            length(src), paste(helpers, collapse = ", "), if (length(log_lines)) paste(trimws(log_lines), collapse = " || ") else
            "none"))
        for (h in helpers) {
            hsrc <- tryCatch(deparse(get(h, envir = ns)), error = function(e) character())
            hlog <- grep("log2|log\\(|log1p|log10", hsrc, value = TRUE)
            scale_lines <- c(scale_lines, sprintf("  helper %s: %d lines; lines naming log: %s", h, length(hsrc),
                if (length(hlog)) paste(trimws(hlog), collapse = " || ") else "none"))
        }
    }
    fixture <- tryCatch(makeFixture(o), error = function(e) structure(NA, err = conditionMessage(e)))
    if (!isS4(fixture)) {
        scale_lines <- c(scale_lines, paste(o, "fixture failed:", attr(fixture, "err")))
        next
    }
    rle <- tryCatch(plotRle(fixture, grouping_variable = "group", yaxis_limit = c(-6, 6)), error = function(e) conditionMessage(e))
    p <- if (is.list(rle) && !inherits(rle, "ggplot")) rle[[1L]] else rle
    if (inherits(p, "ggplot")) {
        nums <- Filter(is.numeric, as.list(p$data))
        for (nm in names(nums)) scale_lines <- c(scale_lines, sprintf("%s plotRle data column %s: range %s; matches %s", o, nm,
            paste(signif(range(nums[[nm]], na.rm = TRUE), 5), collapse = " to "), matchScale(nums[[nm]], rle_raw, rle_log)))
        built <- tryCatch(ggplot2::ggplot_build(p)$data[[1L]], error = function(e) NULL)
        if (!is.null(built)) scale_lines <- c(scale_lines, sprintf("%s plotRle built layer 1 middle: %s", o,
            paste(signif(built$middle %||% NA, 5), collapse = ",")))
    } else scale_lines <- c(scale_lines, paste(o, "plotRle returned:", oneLine(utils::capture.output(print(class(rle))))))
    pca <- tryCatch(plotPca(fixture, grouping_variable = "group", label_column = "", shape_variable = NULL, title = "",
        font_size = 8), error = function(e) conditionMessage(e))
    q <- if (is.list(pca) && !inherits(pca, "ggplot")) pca[[1L]] else pca
    if (inherits(q, "ggplot")) {
        pc_raw <- sort(abs(stats::prcomp(t(fx_mat))$x[, 1L]))
        pc_log <- sort(abs(stats::prcomp(t(log2(fx_mat)))$x[, 1L]))
        nums <- Filter(is.numeric, as.list(q$data))
        for (nm in names(nums)) scale_lines <- c(scale_lines, sprintf("%s plotPca data column %s: %s; |PC1| matches %s", o, nm,
            paste(signif(nums[[nm]], 5), collapse = ","), matchScale(abs(nums[[nm]]), pc_raw, pc_log)))
    } else scale_lines <- c(scale_lines, paste(o, "plotPca returned:", oneLine(if (is.character(pca)) pca else class(pca))))
}
writeUtf8(scale_lines, file.path(out, "g0c", "plotting-scale.txt"))
say(paste(scale_lines, collapse = "\n"))

# Part B (E-3): the four proteomics templates.
prot <- c(LFQ = "inst/reports/proteomics/LFQ_report.rmd", TMT = "inst/reports/proteomics/TMT_report.rmd",
    DIANN = "inst/reports/proteomics/DIANN_report.rmd", DIANN_limpa = "inst/reports/proteomics/DIANN_limpa_report.rmd")
for (t in names(prot)) {
    tp <- prot[[t]]
    tryCatch({
        methodFn <- function(env) fnName(env, "^FormatMethodDetails")
        # Default case (no record, two groups): the text part B must leave unchanged.
        d <- caseDir(paste0(t, "-default"))
        env <- templateEnv(tp, "proteomics", d)
        numbers <- numbersFor(env, design2)
        texts <- c(sample = tryCatch(env$FormatSampleDetails(numbers), error = function(e) paste("<error>", conditionMessage(e))),
            method = tryCatch(env[[methodFn(env)]](numbers), error = function(e) paste("<error>", conditionMessage(e))),
            de = tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e))))
        g0cRecord(paste0(t, "-default"), tp, "default text (recorded)", TRUE, "record",
            c(paste("definition errors:", paste(env$.def_errors, collapse = " | ")), paste("numbers length:", length(numbers)),
                texts))
        # B1: technical replicates declared.
        env <- templateEnv(tp, "proteomics", caseDir(paste0(t, "-B1")))
        numbers <- numbersFor(env, designTr)
        txt <- tryCatch(env$FormatSampleDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
        g0cRecord(paste0(t, "-B1"), tp, "with no technical replicates", has(normText(txt), "with no technical replicates"),
            "present", txt)
        # B2: normalisation 'none' (record step parameter).
        d <- caseDir(paste0(t, "-B2"))
        plantRecord(d, steps = list(list(order = 1L, operation = "normalization-ruv", parameters = list(list(name = "norm_method",
            source = "manifest", value = encodeValue("none"), provenance = "declared")))))
        env <- templateEnv(tp, "proteomics", d)
        numbers <- numbersFor(env, design2)
        txt <- tryCatch(env[[methodFn(env)]](numbers), error = function(e) paste("<error>", conditionMessage(e)))
        g0cRecord(paste0(t, "-B2"), tp, "'none' method from the 'limma' R package", has(normText(txt),
            "'none' method from the 'limma' R package"), "present",
            c(paste("normalisation_method:", toJson(numbers$normalisation_method)), txt))
        # B3: a record whose DA result carries fdr_method BH.
        d <- caseDir(paste0(t, "-B3"))
        plantRecord(d, fdr = "BH")
        env <- templateEnv(tp, "proteomics", d)
        numbers <- numbersFor(env, design2)
        txt <- tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
        g0cRecord(paste0(t, "-B3"), tp, "Benjamini", has(normText(txt), "Benjamini"), "absent", txt)
        # B4: record lfc 0.5 (no treat_applied).
        d <- caseDir(paste0(t, "-B4"))
        plantRecord(d, lfc = 0.5)
        env <- templateEnv(tp, "proteomics", d)
        numbers <- numbersFor(env, design2)
        txt <- tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
        g0cRecord(paste0(t, "-B4"), tp, "test for any deviation from zero log-fold-change",
            has(normText(txt), "test for any deviation from zero log-fold-change"), "present",
            c(paste("ui_de_log_fold_change_cutoff:", toJson(numbers$ui_de_log_fold_change_cutoff), "; s4_de_treat_lfc_cutoff:",
                toJson(numbers$s4_de_treat_lfc_cutoff), "; ui_de_treat_enabled:", toJson(numbers$ui_de_treat_enabled)), txt))
        # B5: record block_reason technical_replicates.
        d <- caseDir(paste0(t, "-B5"))
        plantRecord(d, block_reason = "technical_replicates")
        env <- templateEnv(tp, "proteomics", d)
        numbers <- numbersFor(env, design2)
        txt <- tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
        g0cRecord(paste0(t, "-B5"), tp, "duplicateCorrelation", has(normText(txt), "duplicateCorrelation"), "absent",
            c(paste("da_block_reason:", toJson(numbers$da_block_reason)), txt))
        if (identical(t, "DIANN_limpa")) {
            # B6: the non-limpa branch (no record), and with block_reason.
            env <- templateEnv(tp, "proteomics", caseDir("DIANN_limpa-B6"))
            numbers <- numbersFor(env, designTr)
            txt <- tryCatch(env[[methodFn(env)]](numbers), error = function(e) paste("<error>", conditionMessage(e)))
            g0cRecord("DIANN_limpa-B6", tp, "Technical replicate averaging", has(normText(txt), "Technical replicate averaging"),
                "present", txt)
            d <- caseDir("DIANN_limpa-B6-block")
            plantRecord(d, block_reason = "technical_replicates")
            env <- templateEnv(tp, "proteomics", d)
            numbers <- numbersFor(env, designTr)
            txt <- tryCatch(env[[methodFn(env)]](numbers), error = function(e) paste("<error>", conditionMessage(e)))
            g0cRecord("DIANN_limpa-B6-block", tp, "Technical replicate averaging", has(normText(txt),
                "Technical replicate averaging"), "present", c(paste("da_block_reason:", toJson(numbers$da_block_reason)), txt))
            # R1-N03: the non-limpa branch with a file cutoff of 0.5 and no record; (a) deAnalysisParameters only,
            # (b) with the writer's DIA UI block too.
            for (form in c("a", "b")) {
                d <- caseDir(paste0("DIANN_limpa-N03", form))
                sp <- c("Workflow Parameters:", "-------------------", "", "  Parameters from Final S4 Object:",
                    "  --------------------------------", "[deAnalysisParameters]", "  formula_string = ~ 0 + group",
                    "  treat_lfc_cutoff = 0.5")
                if (identical(form, "b")) sp <- c(sp, "  User Interface Parameters:", "  -------------------------",
                    "[Differential Expression UI Parameters]", "  q_value_threshold = 0.05", "  log_fold_change_cutoff = 0.5",
                    "  treat_enabled = TRUE", "")
                writeUtf8(sp, file.path(d, "study_parameters.txt"))
                env <- templateEnv(tp, "proteomics", d)
                study <- tryCatch(env$ReadStudyParamsNew(file.path(d, "study_parameters.txt")), error = function(e) list())
                numbers <- numbersFor(env, design2, study)
                txt <- tryCatch(env$FormatDEDetails(numbers), error = function(e) paste("<error>", conditionMessage(e)))
                treat_line <- if (has(normText(txt), "The TREAT method was used")) "TREAT sentence" else
                    if (has(normText(txt), "standard empirical Bayes method")) "standard empirical Bayes sentence" else "neither"
                g0cRecord(paste0("DIANN_limpa-N03", form), tp, paste("TREAT sentence printed:", treat_line),
                    has(normText(txt), "The TREAT method was used"), "record",
                    c(paste("parsed deAnalysisParameters:", toJson(study[["Workflow Parameters"]][["Parameters from Final S4 Object"]][["deAnalysisParameters"]])),
                        paste("ui_de_log_fold_change_cutoff:", toJson(numbers$ui_de_log_fold_change_cutoff),
                            "; s4_de_treat_lfc_cutoff:", toJson(numbers$s4_de_treat_lfc_cutoff), "; ui_de_treat_enabled:",
                            toJson(numbers$ui_de_treat_enabled)), txt))
            }
        }
    }, error = function(e) say("G0c ", t, ": harness error: ", conditionMessage(e)))
}
if (length(g0c_rows)) tsv(do.call(rbind, g0c_rows), file.path(out, "g0c", "summary.tsv"))
say("G0 script finished")
