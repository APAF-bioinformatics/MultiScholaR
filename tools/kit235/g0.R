# KIT-235 gate G0 (ticket section 3): observe the lipid/metab report defect before any code. No package change.
#
# Usage (run-tests scripts job, from the package root after package_commit is installed):
#   Rscript tools/kit235/g0.R <staging commit> <main commit> <kit-227 commit>
# Writes kit235-g0/:
#   meta.txt                 preconditions (pandoc, locale, versions, commits)
#   g0a/summary.tsv          per case x template version x bullet shape: the Main outcome, extractor and AnalyzeSamples
#   g0a/<case>-<ver>-<shape>.txt   sample_details, method_details, parsed sections and keys, the fields of section 1.3
#   g0b/                     the GUI lanes' study_parameters.txt, analysis_record.json, design_matrix.tab, report copies
#                            (one copy per changed content), the report's docx text, the lane records, files.tsv
#   g0c/summary.tsv, g0c/<lane>-<ver>-<shape>.txt   the captured GUI files through the staging and kit-227 chunks
# Nothing is asserted: every outcome is recorded; section 3 of the ticket decides from the record.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
refs <- list(staging = arguments[[1L]], main = arguments[[2L]], kit227 = arguments[[3L]])
out <- normalizePath("kit235-g0", mustWork = FALSE)
dir.create(file.path(out, "g0a"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "g0b"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "g0c"), recursive = TRUE, showWarnings = FALSE)
bullet <- "\u2022"

say <- function(...) cat(..., "\n", sep = "")
writeUtf8 <- function(lines, path) writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
asText <- function(x) if (is.null(x)) "<NULL>" else paste(as.character(x), collapse = "")
toJson <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "string", force = TRUE,
    digits = NA))

# --- Preconditions (P10) ------------------------------------------------------------------------------------------
gitOut <- function(args) tryCatch(system2("git", args, stdout = TRUE, stderr = TRUE), error = function(e) conditionMessage(e))
meta <- c(
    paste("R", R.version.string),
    paste("pandoc_available", rmarkdown::pandoc_available()),
    paste("pandoc_version", tryCatch(as.character(rmarkdown::pandoc_version()), error = function(e) conditionMessage(e))),
    paste("locale", Sys.getlocale()),
    paste("rmarkdown", as.character(utils::packageVersion("rmarkdown"))),
    paste("knitr", as.character(utils::packageVersion("knitr"))),
    paste("xml2", tryCatch(as.character(utils::packageVersion("xml2")), error = function(e) "absent")),
    paste("MultiScholaR installed", as.character(utils::packageVersion("MultiScholaR"))),
    paste("HEAD", paste(gitOut(c("rev-parse", "HEAD")), collapse = " ")),
    paste("HEAD^", paste(gitOut(c("rev-parse", "HEAD^")), collapse = " ")),
    paste("origin/v0.9-staging", paste(gitOut(c("rev-parse", "origin/v0.9-staging")), collapse = " ")),
    sprintf("ref %s = %s", names(refs), vapply(refs, function(r) paste(gitOut(c("rev-parse", r)), collapse = " "), "")))
writeUtf8(meta, file.path(out, "meta.txt"))
say(paste(meta, collapse = "\n"))

# --- The templates of three versions, read by commit -----------------------------------------------------------------
omics <- list(
    lipid = list(file = "inst/reports/lipidomics/lipidomics_report.rmd", omic = "lipidomics",
        extract = "ExtractLipidStudyParameters"),
    metab = list(file = "inst/reports/metabolomics/metabolomics_report.rmd", omic = "metabolomics",
        extract = "ExtractStudyParametersNew"))
tdir <- file.path(tempdir(), "kit235-templates")
dir.create(tdir, showWarnings = FALSE)
templatePath <- function(ver, omic) {
    path <- file.path(tdir, paste0(ver, "-", omic, ".rmd"))
    if (!file.exists(path)) {
        status <- system2("git", c("show", paste0(refs[[ver]], ":", omics[[omic]]$file)), stdout = path, stderr = FALSE)
        if (!identical(as.integer(status), 0L)) stop("git show failed for ", ver, " ", omic)
    }
    path
}
for (ver in names(refs)) for (omic in names(omics)) say("template ", ver, " ", omic, ": ",
    tools::md5sum(templatePath(ver, omic)))

chunkLabels <- c("Hierarchical Parser", "Format Sample Details", "Method Section", "Result Section", "Comment Section",
    "Opinions Section", "Differential Abundance Section", "Pathway Enrichment Section", "Sample Analysis Functions",
    "Main Execution Function")

chunkExprs <- function(path, label) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    if (length(start) != 1L) stop("chunk '", label, "' found ", length(start), " times in ", path)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

#' The chunks evaluated as RenderReport evaluates the template: new.env(parent = asNamespace("MultiScholaR"))
#' (R/func_general_filemgmt_report_render_helpers.R:190), with the Project Setup chunk's paths, %||% and a d66_log stub.
#' The Main chunk's own execution (Main() under its handler) runs last and leaves results_data.
runTemplate <- function(ver, omic, source_dir) {
    env <- new.env(parent = asNamespace("MultiScholaR"))
    env$report_source_dir <- source_dir
    env$report_publication_tables_dir <- file.path(source_dir, "publication_tables")
    env$report_results_summary_dir <- source_dir
    env$params <- list(omic_type = omics[[omic]]$omic, experiment_label = "kit235", workflow_name = "kit235")
    env$`%||%` <- function(a, b) if (!is.null(a)) a else b
    env$d66_log <- function(...) invisible(NULL)
    chunk_errors <- character()
    for (label in chunkLabels) {
        for (e in chunkExprs(templatePath(ver, omic), label)) {
            tryCatch(suppressWarnings(eval(e, env)), error = function(err) {
                chunk_errors <<- c(chunk_errors, paste0(label, ": ", conditionMessage(err)))
            })
        }
    }
    env$.chunk_errors <- chunk_errors
    env
}

fieldNames <- c("num_assays", "assay_names", "itsd_applied", "itsd_method_type", "itsd_aggregation", "itsd_removed",
    "itsd_features_per_assay", "log_transform_applied", "log_transform_offset", "normalisation_method", "ruv_applied",
    "ruv_grouping_variable", "ruv_assay_list", "best_k", "num_neg_ctrl", "best_percentage", "ruv_number_k",
    "ruv_per_assay_details", "total_metabolites", "ruv_design_methods")

#' One file through one template version: the template's own Main under its handler, then (outside it) the parser, the
#' extractor alone and AnalyzeSamples' steps one by one, so that an error AnalyzeSamples would swallow is recorded.
runCase <- function(ver, omic, source_dir) {
    env <- runTemplate(ver, omic, source_dir)
    rd <- env$results_data
    study <- tryCatch(env$ReadStudyParamsNew(file.path(source_dir, "study_parameters.txt")),
        error = function(e) structure(list(), parse_error = conditionMessage(e)))
    extracted <- tryCatch(env[[omics[[omic]]$extract]](study), error = function(e) structure(list(),
        extract_error = conditionMessage(e)))
    design <- tryCatch(readr::read_tsv(file.path(source_dir, "design_matrix.tab"), show_col_types = FALSE),
        error = function(e) NULL)
    analyzed <- tryCatch(env$AnalyzeSamples(design, study), error = function(e) list(numbers = list()))
    swallowed <- tryCatch({
        group_summary <- env$AnalyzeGroups(design)
        numbers <- env$CalculateSampleCounts(group_summary)
        numbers <- env$CalculatePatientMetrics(numbers, group_summary)
        params_extracted <- env[[omics[[omic]]$extract]](study)
        if (any(grepl("workflowRecordReportParameters", deparse(env$AnalyzeSamples), fixed = TRUE))) {
            MultiScholaR:::workflowRecordReportParameters(file.path(source_dir, "analysis_record.json"), params_extracted,
                study)
        }
        "none"
    }, error = function(e) conditionMessage(e))
    plain <- function(v) if (is.list(v)) lapply(v, plain) else if (is.atomic(v)) structure(as.vector(v), names = names(v)) else
        paste(class(v), collapse = "/")
    pick <- function(x) stats::setNames(lapply(fieldNames, function(f) if (f %in% names(x)) plain(x[[f]]) else "<absent>"),
        fieldNames)
    list(env = env, results = rd, study = study, extracted = extracted, analyzed = analyzed, swallowed = swallowed,
        extractor_fields = pick(extracted), analyzed_fields = pick(analyzed$numbers))
}

caseText <- function(res) {
    sections <- vapply(names(res$study), function(s) {
        keys <- if (is.list(res$study[[s]])) names(res$study[[s]]) else character()
        paste0("  [", s, "] keys: ", paste0("'", keys, "'", collapse = ", "))
    }, "")
    c("== chunk evaluation errors", if (length(res$env$.chunk_errors)) res$env$.chunk_errors else "none",
        "== sample_details", asText(res$results$content$sample_details),
        "== method_details", asText(res$results$content$method_details),
        "== results (first 400 characters)", substr(asText(res$results$content$results), 1L, 400L),
        "== length(results_data$numbers)", length(res$results$numbers),
        "== parsed sections and keys", sections,
        "== extractor called alone: length", length(res$extracted),
        "== extractor error", asText(attr(res$extracted, "extract_error")),
        "== AnalyzeSamples: length(numbers)", length(res$analyzed$numbers),
        "== error AnalyzeSamples would swallow", res$swallowed,
        "== fields from the extractor (no record overlay)", toJson(res$extractor_fields),
        "== fields after AnalyzeSamples (record overlay when a record is present)", toJson(res$analyzed_fields))
}

summaryRow <- function(label, ver, shape, res) {
    sd <- asText(res$results$content$sample_details)
    md <- asText(res$results$content$method_details)
    data.frame(case = label, version = ver, shape = shape,
        main_error = grepl("Error running Main function", sd, fixed = TRUE),
        sample_details_head = substr(gsub("[\t\r\n]+", " ", sd), 1L, 120L),
        method_details_chars = nchar(md),
        method_has_METHOD_DETAILS = grepl("METHOD DETAILS", md, fixed = TRUE),
        extractor_length = length(res$extracted), analyze_numbers_length = length(res$analyzed$numbers),
        swallowed_error = gsub("[\t\r\n]+", " ", res$swallowed), chunk_errors = length(res$env$.chunk_errors),
        num_assays = toJson(res$extractor_fields$num_assays), ruv_applied = toJson(res$extractor_fields$ruv_applied),
        ruv_assay_list_names = paste(names(res$extractor_fields$ruv_assay_list), collapse = "|"),
        stringsAsFactors = FALSE)
}

toBulletShape <- function(lines) sub("^\\* ", paste0(bullet, " "), lines)

# --- G0a: the writer's own files through three template versions -------------------------------------------------------
ns <- asNamespace("MultiScholaR")
nsGet <- function(name) get(name, envir = ns)

# Fixtures of tests/testthat/test-general-filemgmt-s4-contracts.R:12-107, copied (ruv = FALSE drops the RUV args: S).
makeMetab <- function(ruv = TRUE) {
    design_matrix <- data.frame(Sample_ID = c("S1", "S2", "S3"), group = c("A", "A", "B"), stringsAsFactors = FALSE)
    assay_pos <- data.frame(database_identifier = c("M1", "M2"), metabolite_identification = c("Met1", "Met2"),
        S1 = c(10, 20), S2 = c(11, 21), S3 = c(12, 22), check.names = FALSE)
    assay_neg <- data.frame(database_identifier = c("M3", "M4", "M5"), metabolite_identification = c("Met3", "Met4", "Met5"),
        S1 = c(30, 40, 50), S2 = c(31, 41, 51), S3 = c(32, 42, 52), check.names = FALSE)
    args <- list(
        ITSDNormalization = list(applied = TRUE, method_type = "median", itsd_aggregation = "sum",
            itsd_pattern_columns = c("ISTD_1", "ISTD_2"), removed_itsd = TRUE,
            timestamp = as.POSIXct("2026-04-18 09:00:00", tz = "UTC"),
            itsd_features_per_assay = list(assay_pos = c("ISTD_A1", "ISTD_A2"), assay_neg = c("ISTD_B1")),
            itsd_counts_per_assay = list(assay_pos = 2L, assay_neg = 1L)),
        log_transformed = TRUE, log_transform_offset = 1, normalisation_method = "PQN")
    if (ruv) args <- c(args, list(ruv_grouping_variable = "group", ruv_number_k = list(assay_pos = 2L, assay_neg = 3L),
        ctrl = list(assay_pos = c(TRUE, FALSE), assay_neg = c(TRUE, TRUE, FALSE))))
    nsGet("createMetaboliteAssayData")(metabolite_data = list(assay_pos = assay_pos, assay_neg = assay_neg),
        design_matrix = design_matrix, metabolite_id_column = "database_identifier",
        annotation_id_column = "metabolite_identification", sample_id = "Sample_ID", group_id = "group",
        database_identifier_type = "HMDB", internal_standard_regex = "^ISTD", args = args)
}
makeLipid <- function(ruv = c("single", "per_assay", "none")) {
    ruv <- match.arg(ruv)
    design_matrix <- data.frame(Sample_ID = c("S1", "S2", "S3"), group = c("A", "A", "B"), stringsAsFactors = FALSE)
    assay_data <- data.frame(database_identifier = c("L1", "L2", "L3"), lipid_identification = c("Lipid1", "Lipid2", "Lipid3"),
        S1 = c(100, 200, 300), S2 = c(101, 201, 301), S3 = c(102, 202, 302), check.names = FALSE)
    args <- list(ITSDNormalization = list(applied = FALSE), log_transformed = FALSE, normalisation_method = "median")
    if (identical(ruv, "single")) args <- c(args, list(ruv_grouping_variable = "group", ruv_number_k = 4L,
        ctrl = c(TRUE, FALSE, TRUE)))
    if (identical(ruv, "per_assay")) args <- c(args, list(ruv_grouping_variable = "group",
        ruv_number_k = list(assay_main = 4L), ctrl = list(assay_main = c(TRUE, FALSE, TRUE))))
    nsGet("createLipidomicsAssayData")(lipid_data = list(assay_main = assay_data), design_matrix = design_matrix,
        lipid_id_column = "database_identifier", annotation_id_column = "lipid_identification", sample_id = "Sample_ID",
        group_id = "group", database_identifier_type = "LipidMaps", internal_standard_regex = "^LIPID_ISTD", args = args)
}

assign("config_list", list(), envir = .GlobalEnv)
writeCaseFile <- function(label, object, workflow_name) {
    dir <- file.path(tempdir(), "kit235-g0a", label)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    invisible(utils::capture.output(nsGet("createWorkflowArgsFromConfig")(workflow_name = workflow_name,
        source_dir_path = dir, final_s4_object = object,
        contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)), type = "output"))
    file.path(dir, "study_parameters.txt")
}
designLines <- c("Run\tgroup", "S1\tA", "S2\tA", "S3\tB")

cases <- list(
    M = list(omic = "metab", make = function() makeMetab(TRUE), wf = "metabolomics"),
    L1 = list(omic = "lipid", make = function() makeLipid("single"), wf = "lipidomics"),
    L2 = list(omic = "lipid", make = function() makeLipid("per_assay"), wf = "lipidomics"),
    `S-lip` = list(omic = "lipid", make = function() makeLipid("none"), wf = "lipidomics"),
    `S-met` = list(omic = "metab", make = function() makeMetab(FALSE), wf = "metabolomics"))

g0a_rows <- list()
for (label in names(cases)) {
    case <- cases[[label]]
    written <- tryCatch(writeCaseFile(label, case$make(), case$wf), error = function(e) {
        say("G0a ", label, ": writer failed: ", conditionMessage(e)); NULL
    })
    if (is.null(written)) {
        g0a_rows[[length(g0a_rows) + 1L]] <- data.frame(case = label, version = "-", shape = "-", main_error = NA,
            sample_details_head = "writer failed", method_details_chars = NA, method_has_METHOD_DETAILS = NA,
            extractor_length = NA, analyze_numbers_length = NA, swallowed_error = NA, chunk_errors = NA, num_assays = NA,
            ruv_applied = NA, ruv_assay_list_names = NA, stringsAsFactors = FALSE)
        next
    }
    star <- readLines(written, warn = FALSE, encoding = "UTF-8")
    writeUtf8(star, file.path(out, "g0a", paste0(label, "-study_parameters.txt")))
    for (shape in c("star", "bullet")) {
        for (ver in c("main", "staging", "kit227")) {
            dir <- file.path(tempdir(), "kit235-g0a-run", label, shape, ver)
            dir.create(dir, recursive = TRUE, showWarnings = FALSE)
            writeUtf8(if (identical(shape, "star")) star else toBulletShape(star), file.path(dir, "study_parameters.txt"))
            writeUtf8(designLines, file.path(dir, "design_matrix.tab"))
            res <- tryCatch(runCase(ver, case$omic, dir), error = function(e) {
                say("G0a ", label, " ", ver, " ", shape, ": harness error: ", conditionMessage(e)); NULL
            })
            if (is.null(res)) next
            writeUtf8(caseText(res), file.path(out, "g0a", paste0(label, "-", ver, "-", shape, ".txt")))
            g0a_rows[[length(g0a_rows) + 1L]] <- summaryRow(label, ver, shape, res)
        }
    }
}
g0a <- do.call(rbind, g0a_rows)
utils::write.table(g0a, file.path(out, "g0a", "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(g0a[, c("case", "version", "shape", "main_error", "method_has_METHOD_DETAILS", "extractor_length",
    "analyze_numbers_length", "swallowed_error")])

# --- G0b: the application's real report (the GUI-entry lanes) ----------------------------------------------------------
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
say("loaded from the tree: ", as.character(utils::packageVersion("MultiScholaR")))
# The test helpers are sourced into the namespace (pkgload) or onto the search path: found either way.
helperFn <- function(name) {
    dev_ns <- asNamespace("MultiScholaR")
    if (exists(name, envir = dev_ns, inherits = FALSE)) get(name, envir = dev_ns) else get(name, envir = globalenv())
}
runLane <- helperFn("guiEntryRunLane")
say("guiEntryRunLane found in ", environmentName(environment(runLane)))
docxText <- function(path) {
    dir <- tempfile("docx-")
    utils::unzip(path, files = "word/document.xml", exdir = dir)
    xml <- paste(readLines(file.path(dir, "word", "document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
    text <- gsub("</w:p>", "\n", xml, fixed = TRUE)
    text <- gsub("<[^>]+>", "", text)
    for (pair in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&apos;", "'"), c("&amp;", "&"))) {
        text <- gsub(pair[[1L]], pair[[2L]], text, fixed = TRUE)
    }
    text
}

capture_state <- new.env(parent = emptyenv())
capture_state$lane <- NULL
capture_state$n <- 0L
capture_state$seen <- character()
capture_state$files <- data.frame(lane = character(), capture = integer(), source = character(), saved = character(),
    md5 = character(), stringsAsFactors = FALSE)
capture_state$last <- list()
kit235Capture <- function(base_dir) {
    lane <- capture_state$lane
    if (is.null(lane)) return(invisible(NULL))
    capture_state$n <- capture_state$n + 1L
    files <- list.files(base_dir, recursive = TRUE, all.files = TRUE)
    files <- files[!grepl("(^|/)\\.multischolar(/|$)|(^|/)\\.workflow|\\.duckdb|\\.qs2$|\\.wal$", files)]
    keep <- files[grepl("(^|/)(study_parameters\\.txt|analysis_record\\.json|design_matrix\\.tab)$|\\.docx$|_report\\.rmd$",
        files)]
    for (f in keep) {
        src <- file.path(base_dir, f)
        md5 <- unname(tools::md5sum(src))
        key <- paste(lane, f, md5)
        if (key %in% capture_state$seen) next
        capture_state$seen <- c(capture_state$seen, key)
        # The report reads scripts/<omic>/ (report_source_dir); results_summary/<omic>/Study_report/ holds copies.
        short <- if (grepl("\\.docx$", f)) {
            if (grepl("analysis_record", basename(f), fixed = TRUE)) "record.docx" else "report.docx"
        } else if (grepl("_report\\.rmd$", f)) "template.rmd" else {
            paste0(if (grepl("^scripts/", f)) "src-" else "sr-", basename(f))
        }
        saved <- sprintf("%s-c%02d-%s", lane, capture_state$n, short)
        file.copy(src, file.path(out, "g0b", saved), overwrite = TRUE)
        capture_state$files[nrow(capture_state$files) + 1L, ] <- list(lane, capture_state$n, f, saved, md5)
        capture_state$last[[paste(lane, short)]] <- file.path(out, "g0b", saved)
    }
    invisible(NULL)
}
assign("kit235Capture", kit235Capture, envir = .GlobalEnv)
trace_ok <- tryCatch({
    trace("guiEntryArtifacts", tracer = quote(kit235Capture(base_dir)), where = runLane, print = FALSE)
    TRUE
}, error = function(e) { say("trace failed: ", conditionMessage(e)); FALSE })
if (!trace_ok) {
    # Fallback with the same effect: the helper's binding wrapped in its own environment.
    helper_env <- environment(runLane)
    original_artifacts <- get("guiEntryArtifacts", envir = helper_env)
    if (bindingIsLocked("guiEntryArtifacts", helper_env)) unlockBinding("guiEntryArtifacts", helper_env)
    assign("guiEntryArtifacts", function(base_dir, roots) {
        kit235Capture(base_dir)
        original_artifacts(base_dir, roots)
    }, envir = helper_env)
    say("guiEntryArtifacts wrapped in ", environmentName(helper_env))
}
say("trace installed: ", trace_ok)

record_dir <- file.path(out, "g0b", "records")
lanes <- c(lip = "lipid_canonical", met = "metab_lc")
lane_status <- character()
for (short in names(lanes)) {
    capture_state$lane <- short
    capture_state$n <- 0L
    started <- Sys.time()
    status <- tryCatch({
        withr::with_envvar(c(MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = record_dir), runLane(lanes[[short]]))
        "ran"
    }, error = function(e) paste("error:", conditionMessage(e)))
    lane_status[[short]] <- sprintf("%s %s: %s (%.0f s, %d captures)", short, lanes[[short]], status,
        as.numeric(difftime(Sys.time(), started, units = "secs")), capture_state$n)
    say(lane_status[[short]])
}
capture_state$lane <- NULL
utils::write.table(capture_state$files, file.path(out, "g0b", "files.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

g0b_lines <- c(paste("trace installed:", trace_ok), unname(lane_status))
for (short in names(lanes)) {
    docx <- capture_state$last[[paste(short, "report.docx")]]
    if (is.null(docx)) {
        g0b_lines <- c(g0b_lines, paste(short, "no report docx captured"))
        next
    }
    text <- tryCatch(docxText(docx), error = function(e) paste("<docx read error>", conditionMessage(e)))
    writeUtf8(text, file.path(out, "g0b", paste0(short, "-report-text.txt")))
    g0b_lines <- c(g0b_lines,
        sprintf("%s report %s: contains 'Error running Main function' = %s; contains 'METHOD DETAILS' = %s", short,
            basename(docx), grepl("Error running Main function", text, fixed = TRUE),
            grepl("METHOD DETAILS", text, fixed = TRUE)))
}
writeUtf8(g0b_lines, file.path(out, "g0b", "summary.txt"))
say(paste(g0b_lines, collapse = "\n"))

# --- G0c: the GUI files through the fixed reading (the bullet shape) ---------------------------------------------------
g0c_rows <- list()
g0c_checks <- character()
for (short in names(lanes)) {
    sp <- capture_state$last[[paste(short, "src-study_parameters.txt")]]
    if (is.null(sp)) {
        g0c_checks <- c(g0c_checks, paste(short, "no study_parameters.txt captured"))
        next
    }
    raw <- readLines(sp, warn = FALSE, encoding = "UTF-8")
    omic <- if (identical(short, "lip")) "lipid" else "metab"
    shape_facts <- c(
        sprintf("%s file %s: %d lines; lines starting '* ' = %d; starting '%s' = %d; starting '- ' = %d", short, basename(sp),
            length(raw), sum(grepl("^\\* ", raw)), bullet, sum(startsWith(raw, bullet)), sum(grepl("^- ", raw))),
        sprintf("%s ITSD Features Per Assay block present = %s", short, any(grepl("ITSD Features Per Assay", raw, fixed = TRUE))),
        sprintf("%s RUV form: single ('k value:' line) = %s; per-assay ('Per-Assay Parameters') = %s; 'Status: Not Applied' = %s; 'Status: Applied' = %s",
            short, any(grepl("^[*-] k value:|^\u2022 k value:", raw)), any(grepl("Per-Assay Parameters", raw, fixed = TRUE)),
            any(grepl("Status: Not Applied", raw, fixed = TRUE)), any(grepl("Status: Applied$", raw))),
        sprintf("%s Automatic RUV Optimization Results section = %s", short,
            any(grepl("Automatic RUV Optimization Results", raw, fixed = TRUE))))
    g0c_checks <- c(g0c_checks, shape_facts)
    for (shape in c("star", "bullet")) {
        for (ver in c("staging", "kit227")) {
            dir <- file.path(tempdir(), "kit235-g0c", short, shape, ver)
            dir.create(dir, recursive = TRUE, showWarnings = FALSE)
            writeUtf8(if (identical(shape, "star")) raw else toBulletShape(raw), file.path(dir, "study_parameters.txt"))
            for (other in c("analysis_record.json", "design_matrix.tab")) {
                src <- capture_state$last[[paste(short, paste0("src-", other))]]
                if (!is.null(src)) file.copy(src, file.path(dir, other), overwrite = TRUE)
            }
            res <- tryCatch(runCase(ver, omic, dir), error = function(e) {
                say("G0c ", short, " ", ver, " ", shape, ": harness error: ", conditionMessage(e)); NULL
            })
            if (is.null(res)) next
            writeUtf8(caseText(res), file.path(out, "g0c", paste0(short, "-", ver, "-", shape, ".txt")))
            g0c_rows[[length(g0c_rows) + 1L]] <- summaryRow(short, ver, shape, res)
            if (identical(shape, "bullet")) {
                md <- asText(res$results$content$method_details)
                f <- res$analyzed_fields
                g0c_checks <- c(g0c_checks, sprintf(paste(
                    "%s %s bullet: num_assays=%s; assay_names=%s; log_transform_applied=%s; normalisation_method=%s;",
                    "itsd_applied=%s; itsd_method_type=%s; itsd_aggregation=%s; ruv_applied=%s;",
                    "Methods say ITSD 'was applied' = %s; 'Internal standard normalization was not applied' = %s;",
                    "'using the NA method' = %s; RUV 'was not applied' = %s; 'without log transformation' = %s"),
                    short, ver, toJson(f$num_assays), toJson(f$assay_names), toJson(f$log_transform_applied),
                    toJson(f$normalisation_method), toJson(f$itsd_applied), toJson(f$itsd_method_type),
                    toJson(f$itsd_aggregation), toJson(f$ruv_applied),
                    grepl("Internal standard (ITSD) normalization was applied", md, fixed = TRUE),
                    grepl("Internal standard normalization was not applied", md, fixed = TRUE),
                    grepl("using the NA method", md, fixed = TRUE),
                    grepl("RUV-III batch correction was not applied", md, fixed = TRUE),
                    grepl("without log transformation", md, fixed = TRUE)))
            }
        }
    }
}
if (length(g0c_rows)) utils::write.table(do.call(rbind, g0c_rows), file.path(out, "g0c", "summary.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
writeUtf8(g0c_checks, file.path(out, "g0c", "checks.txt"))
say(paste(g0c_checks, collapse = "\n"))
say("G0 script finished")
