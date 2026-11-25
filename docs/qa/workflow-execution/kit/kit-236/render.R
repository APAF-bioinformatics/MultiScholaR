# KIT-236 correction pass, I04: render the six reports from the head on the GUI lanes' real outputs, and read what
# they print. Lanes: lipid_canonical, metab_lc, prot_lfq, prot_tmt, prot_dia, prot_dia_limpa (guiEntryRunLane, as
# the GUI-entry tests run them). For lipidomics and metabolomics, the lane's project is copied and rendered again with
# one condition changed (V1 to V7, as G0b). Every report's text is written out (docx: </w:p> to a line break, other
# tags removed); each phrase of the prediction table is looked up (whitespace collapsed, curly quotes made straight).
# Nothing is asserted: the table records predicted against observed.
#
# Usage (run-tests scripts job, from the package root): Rscript docs/qa/workflow-execution/kit/kit-236/render.R

out <- normalizePath("kit236-render", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
work <- file.path(tempdir(), "kit236-render-work")
dir.create(work, showWarnings = FALSE)
say <- function(...) cat(..., "\n", sep = "")
writeUtf8 <- function(lines, path) writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
ns <- asNamespace("MultiScholaR")
writeUtf8(c(paste("R", R.version.string), paste("pandoc", tryCatch(as.character(rmarkdown::pandoc_version()),
    error = function(e) "none")), paste("HEAD", paste(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), collapse = " "))),
    file.path(out, "meta.txt"))

normText <- function(x) {
    x <- paste(x, collapse = "\n")
    for (q in c("‘", "’")) x <- gsub(q, "'", x, fixed = TRUE)
    for (q in c("“", "”")) x <- gsub(q, "\"", x, fixed = TRUE)
    x <- gsub(" ", " ", x, fixed = TRUE)
    trimws(gsub("[[:space:]]+", " ", x))
}
readReport <- function(path) {
    dir <- tempfile("docx-")
    utils::unzip(path, files = "word/document.xml", exdir = dir)
    xml <- paste(readLines(file.path(dir, "word", "document.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
    text <- gsub("</w:p>", "\n", xml, fixed = TRUE)
    text <- gsub("<[^>]+>", "", text)
    for (pair in list(c("&lt;", "<"), c("&gt;", ">"), c("&quot;", "\""), c("&apos;", "'"), c("&amp;", "&"))) {
        text <- gsub(pair[[1L]], pair[[2L]], text, fixed = TRUE)
    }
    list(raw = text, norm = normText(text), drawings = sum(gregexpr("<w:drawing>", xml, fixed = TRUE)[[1L]] > 0L))
}

# --- The predictions (written before the run; KIT-236 outcome.md, "Render predictions") -----------------------------
P <- function(...) list(phrase = paste0(...), predicted = "present")
A <- function(...) list(phrase = paste0(...), predicted = "absent")
omicPredictions <- function(word, omic, features, count) c(list(
    P("The analysis included 6 samples, each a separate biological sample (no technical replicates were declared in the design)."),
    A("with no technical replicates."),
    P("Data analyses were performed with MultiScholaR 0.9.0 (build "),
    P("SOP reference: [to be completed by the analyst]."),
    A("SOP xxx"), A("Salesforce"),
    P("Each assay was processed and tested for differential abundance separately; the results of the 2 assays are reported together."),
    A("before integration"), A("cross-assay consistency"),
    P("Intensities were log2-transformed; a value of zero was read as not detected"),
    A("were coded as NA"),
    P("Internal standard (ITSD) normalization was applied using the average_centered method with median aggregation"),
    P("Between-sample normalization was performed using the 'cyclicloess' method from the 'limma' R package"),
    P("Per-assay quality control: the quality-control figures were drawn separately for each assay. The assays were required to contain the same samples."),
    P("Test conditions, deviations from the method and critical reagent details: [to be completed by the analyst]."),
    P("one column per stage: after filtering, on the intensities before log2 transformation (left); after cyclic loess normalisation (normalizeCyclicLoess) (right). RUV-III correction was not applied to this dataset."),
    P("this report draws no conclusion from them."),
    A("Cyclic loess normalization (rightmost column"), A("high and consistent"),
    P("Features remaining after each step (total across the 2 assays): import and design ", features[[1L]],
        "; intensity filter ", features[[2L]], "; duplicate resolution ", features[[2L]], "; quality-control finalisation ",
        features[[2L]], "; internal-standard normalisation ", features[[2L]], "; log2 transformation ", features[[2L]],
        "; between-sample normalisation ", features[[2L]], "; sample-correlation filter ", features[[2L]], "."),
    P("6 of 6 samples were kept by the sample-correlation filter (threshold 0.5)."),
    A("norm complete"),
    P("[To be completed by the analyst: quality-control assessment from Figure 1"),
    A("All samples passed initial QC criteria"),
    P("[To be completed by the analyst: interpretation of the results"),
    P("Differential abundance of ", word, " was tested separately within each assay, using the adjusted abundance matrix, for the comparisons defined in the contrast table (1 comparison: KO_vs_WT)."),
    P("p-values were calculated using the 'eBayes' function."),
    P("False discovery rates were estimated as q-values (Storey, 2002); where they cannot be computed, Benjamini-Hochberg adjusted p-values are used."),
    P("Each comparison was tested for any difference from zero log2 fold change."),
    P("Differentially abundant ", word, " were those with a q-value below 0.05."),
    A("Volcano plots of differentially abundant ", word, " are shown in Figure 3"),
    P("Full details of all ", word, " are provided in the supplementary results table DA_results_", omic, ".xlsx"),
    A("DE_", if (identical(word, "lipids")) "lipids" else "metabolites", "_results.xlsx"),
    A("Figure 3: Volcano plots"),
    A("Volcano Plots directory not found"), A("NumSigDeMolecules directory not found"), A("Composite QC figure not found"),
    A("test could not be computed")),
    if (count) list(P("The number of significant differentially abundant ", word, " is summarised in Figure 2 (q-value < 0.05). For the comparison KO_vs_WT in LCMS_Pos, there were"),
        P("Figure 2: Number of significant differentially abundant ", word), A("groupKO-groupWT")) else
        list(A("is summarised in Figure 2"), A("Figure 2: Number of significant")))
protPredictions <- function(limpa) c(list(
    A("biological samples with no technical replicates."),
    A("The design is not available to this report."),
    P("False discovery rates were estimated as q-values (Storey, 2002)"),
    A("The standard empirical Bayes method was used to test for any deviation from zero log-fold-change."),
    A("Technical replicate averaging"),
    A("'none' method from the 'limma'")),
    if (!limpa) list(P("Each comparison was tested for any difference from zero log2 fold change.")))
lanes <- list(
    lip = list(id = "lipid_canonical", omic = "lipidomics", template = "lipidomics_report.rmd",
        predictions = omicPredictions("lipids", "lipidomics", c(20, 17), TRUE)),
    met = list(id = "metab_lc", omic = "metabolomics", template = "metabolomics_report.rmd",
        predictions = omicPredictions("metabolites", "metabolomics", c(6, 6), FALSE)),
    lfq = list(id = "prot_lfq", omic = "proteomics", predictions = protPredictions(FALSE)),
    tmt = list(id = "prot_tmt", omic = "proteomics", predictions = protPredictions(FALSE)),
    dia = list(id = "prot_dia", omic = "proteomics", predictions = protPredictions(FALSE)),
    dial = list(id = "prot_dia_limpa", omic = "proteomics", predictions = protPredictions(TRUE)))

checkRows <- list()
check <- function(lane, variant, report, predictions) {
    for (p in predictions) {
        observed <- if (grepl(normText(p$phrase), report$norm, fixed = TRUE)) "present" else "absent"
        checkRows[[length(checkRows) + 1L]] <<- data.frame(lane = lane, variant = variant, predicted = p$predicted,
            observed = observed, verdict = if (identical(observed, p$predicted)) "as predicted" else "DIFFERS",
            phrase = p$phrase, stringsAsFactors = FALSE)
    }
}

# --- The lanes, captured at each step (the project copied when a report exists) -------------------------------------
helperFn <- function(name) if (exists(name, envir = ns, inherits = FALSE)) get(name, envir = ns) else get(name, envir = globalenv())
runLane <- helperFn("guiEntryRunLane")
store_pattern <- "(^|/)\\.multischolar(/|$)|(^|/)\\.workflow|\\.duckdb|\\.qs2$|\\.wal$"
cap <- new.env(parent = emptyenv())
cap$lane <- NULL
cap$base <- list()
cap$dirs <- list()
kit236RenderCapture <- function(base_dir) {
    lane <- cap$lane
    if (is.null(lane)) return(invisible(NULL))
    files <- list.files(base_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    files <- files[!grepl(store_pattern, files)]
    if (!any(grepl("^results_summary/.*\\.docx$", files) & !grepl("analysis_record", files))) return(invisible(NULL))
    dest <- file.path(work, lane, "base")
    unlink(dest, recursive = TRUE)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    dirs <- list.dirs(base_dir, recursive = TRUE, full.names = FALSE)
    for (d in dirs[nzchar(dirs) & !grepl(store_pattern, dirs)]) dir.create(file.path(dest, d), recursive = TRUE, showWarnings = FALSE)
    for (f in files) file.copy(file.path(base_dir, f), file.path(dest, f), overwrite = TRUE)
    cap$base[[lane]] <- base_dir
    cap$dirs[[lane]] <- get0("project_dirs", envir = globalenv(), inherits = FALSE)
    invisible(NULL)
}
trace_ok <- tryCatch({
    trace("guiEntryArtifacts", tracer = bquote(.(kit236RenderCapture)(base_dir)), where = runLane, print = FALSE)
    TRUE
}, error = function(e) { say("trace failed: ", conditionMessage(e)); FALSE })
say("trace installed: ", trace_ok)

status <- character()
for (short in names(lanes)) {
    cap$lane <- short
    started <- Sys.time()
    result <- tryCatch({ runLane(lanes[[short]]$id); "ran" }, error = function(e) paste("error:", conditionMessage(e)))
    status[[short]] <- sprintf("%s %s: %s (%.0f s)", short, lanes[[short]]$id, result,
        as.numeric(difftime(Sys.time(), started, units = "secs")))
    say(status[[short]])
}
cap$lane <- NULL

mapPaths <- function(pd, old, new) rapply(pd, function(s) ifelse(startsWith(s, old), paste0(new, substring(s, nchar(old) + 1L)),
    s), classes = "character", how = "replace")
reportOf <- function(dir) {
    docx <- list.files(dir, pattern = "\\.docx$", full.names = TRUE)
    docx[!grepl("analysis_record", basename(docx), fixed = TRUE)]
}
for (short in names(lanes)) {
    lane <- lanes[[short]]
    tryCatch({
        if (is.null(cap$dirs[[short]])) stop("no report captured")
        paths <- mapPaths(cap$dirs[[short]], cap$base[[short]], file.path(work, short, "base"))[[lane$omic]]
        docx <- reportOf(paths$results_summary_dir)
        report <- readReport(docx[[1L]])
        writeUtf8(c(paste("report:", basename(docx[[1L]])), paste("drawings:", report$drawings), report$raw),
            file.path(out, paste0(short, "-text.txt")))
        check(short, "lane", report, lane$predictions)
    }, error = function(e) {
        status[[short]] <<- paste(status[[short]], "| read:", conditionMessage(e))
        say(short, ": ", conditionMessage(e))
    })
}

# --- Lipidomics and metabolomics variants through RenderReport ------------------------------------------------------
replaceGroup <- function(text, pattern, value) {
    m <- gregexpr(pattern, text, perl = TRUE)[[1L]]
    if (m[[1L]] == -1L) return(list(text = text, n = 0L))
    st <- attr(m, "capture.start")[, "v"]
    ln <- attr(m, "capture.length")[, "v"]
    for (i in rev(seq_along(st))) text <- paste0(substr(text, 1L, st[[i]] - 1L), value, substr(text, st[[i]] + ln[[i]], nchar(text)))
    list(text = text, n = length(st))
}
editRecord <- function(path, pattern, value) {
    text <- readChar(path, file.size(path), useBytes = TRUE)
    res <- replaceGroup(text, pattern, value)
    writeLines(res$text, path, sep = "", useBytes = TRUE)
    res$n
}
editFile <- function(path, fn) writeUtf8(fn(readLines(path, warn = FALSE, encoding = "UTF-8")), path)
touchPng <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(path, width = 200, height = 200)
    graphics::plot.new()
    grDevices::dev.off()
}
words <- list(lip = "lipids", met = "metabolites")
variants <- list(
    V1 = list(edit = function(p, short) {
        editRecord(file.path(p$source_dir, "analysis_record.json"),
            '"name": "normalisation_method",\\s*"source": "[^"]*",\\s*"value": \\{\\s*"type": "character",\\s*"value": (?<v>"[^"]*")', '"none"')
        editFile(file.path(p$source_dir, "study_parameters.txt"), function(l) sub("^\\* Method: .*$", "* Method: none", l))
    }, predictions = function(word) list(P("No between-sample normalisation was applied."),
        P("after the normalisation step, at which no normalisation was applied (right)"), A("'none' method from"))),
    V2 = list(edit = function(p, short) {
        editRecord(file.path(p$source_dir, "analysis_record.json"),
            '"name": "offset",\\s*"source": "[^"]*",\\s*"value": \\{\\s*"type": "double",\\s*"value": (?<v>[-0-9.eE]+)', "1")
        editFile(file.path(p$source_dir, "study_parameters.txt"), function(l) sub("^\\* Transform: log2$", "* Offset: 1", l))
    }, predictions = function(word) list(P("intensities were transformed as log2(x + 1). Negative values were set to zero before the transformation, zeros were transformed with the offset, and cells with no value stayed missing."),
        A("were coded as NA"))),
    V3 = list(edit = function(p, short) editRecord(file.path(p$source_dir, "analysis_record.json"),
        '"tech_rep_group": (?<v>"KO_2")', '"KO_1"'),
        predictions = function(word) list(P("The analysis included 6 measured runs from 5 biological samples; 1 biological sample was measured more than once (technical replicates), giving 1 additional runs."))),
    V4 = list(edit = function(p, short) {
        touchPng(file.path(p$publication_figures_dir, "Volcano_Plots", "kit236_placeholder.png"))
        touchPng(file.path(p$publication_figures_dir, "NumSigDeMolecules", "kit236_placeholder.png"))
    }, predictions = function(word) list(P("Volcano plots of differentially abundant ", word, " are shown in Figure 3, with a threshold of q-value < 0.05."),
        P("Figure 3: Volcano plots of differentially abundant ", word), P("Figure 2: Number of significant differentially abundant ", word))),
    V5 = list(edit = function(p, short) unlink(list.files(p$publication_tables_dir, pattern = "^DA_results_.*\\.xlsx$", full.names = TRUE)),
        predictions = function(word) list(A("supplementary results table"))),
    V6 = list(edit = function(p, short) editFile(file.path(p$source_dir, "study_parameters.txt"), function(l) {
        at <- grep("^\\* ITSD Removed After Normalization:", l)[[1L]]
        if (at < length(l) && grepl("^\\* Timestamp:", l[[at + 1L]])) at <- at + 1L
        append(l, c("* ITSD Features Per Assay:", "    - LCMS_Pos (2 features): IS_kit236_a, IS_kit236_b"), after = at)
    }), predictions = function(word) list(P("on the intensities before internal-standard normalisation and log2 transformation (left)"),
        P("Internal standard (ITSD) normalization was applied using the average_centered method"))),
    V7 = list(edit = function(p, short) {
        rec <- file.path(p$source_dir, "analysis_record.json")
        editRecord(rec, '"name": "lfc",\\s*"source": "[^"]*",\\s*"value": \\{\\s*"type": "double",\\s*"value": (?<v>[-0-9.eE]+)', "0.5")
        editRecord(rec, '"treat_applied": (?<v>false)', "true")
    }, predictions = function(word) list(P("Each comparison was tested against a log2 fold-change threshold of 0.5 with limma's TREAT method (McCarthy and Smyth, 2009)."),
        P("p-values were calculated with limma's TREAT test against the fold-change threshold"),
        P("Differentially abundant ", word, " were those with a q-value below 0.05 and an absolute log2 fold change of at least 0.5."),
        A("p-values were calculated using the 'eBayes' function"))))

for (short in c("lip", "met")) {
    lane <- lanes[[short]]
    if (is.null(cap$dirs[[short]])) next
    for (vid in names(variants)) {
        tryCatch({
            root <- file.path(tempdir(), "kit236-variants", short, vid)
            unlink(root, recursive = TRUE)
            dir.create(root, recursive = TRUE)
            file.copy(file.path(work, short, "base"), root, recursive = TRUE)
            pd <- mapPaths(cap$dirs[[short]], cap$base[[short]], file.path(root, "base"))
            p <- pd[[lane$omic]]
            for (old in reportOf(p$results_summary_dir)) unlink(old)
            variants[[vid]]$edit(p, short)
            assign("project_dirs", pd, envir = globalenv())
            rendered <- RenderReport(omic_type = lane$omic, experiment_label = lane$id, rmd_filename = lane$template)
            if (is.null(rendered) || !file.exists(rendered)) stop("render failed")
            report <- readReport(rendered)
            writeUtf8(c(paste("drawings:", report$drawings), report$raw), file.path(out, paste0(short, "-", vid, "-text.txt")))
            check(short, vid, report, variants[[vid]]$predictions(words[[short]]))
        }, error = function(e) {
            checkRows[[length(checkRows) + 1L]] <<- data.frame(lane = short, variant = vid, predicted = "render",
                observed = paste("error:", conditionMessage(e)), verdict = "DIFFERS", phrase = "", stringsAsFactors = FALSE)
        })
    }
}

writeUtf8(c(paste("trace installed:", trace_ok), unname(status)), file.path(out, "lanes.txt"))
if (length(checkRows)) {
    table <- do.call(rbind, checkRows)
    utils::write.table(table, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    say(sum(table$verdict == "DIFFERS"), " of ", nrow(table), " checks differ from the prediction")
}
say("render script finished")
