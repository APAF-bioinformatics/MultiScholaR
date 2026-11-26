# KIT-259 third pass (VAL-307; AUDIT-KIT-259-impl-20261126 closure check C.4): render the four proteomics reports from
# the GUI lanes' real projects (prot_lfq, prot_tmt, prot_dia, prot_dia_limpa), with KIT-236's harness
# (docs/qa/workflow-execution/kit/kit-236/render.R: guiEntryRunLane, the project copied when its report exists,
# RenderReport, the docx text). Variants on a copy of each lane's project:
#   V0: the lane's project as run, rendered with the head's template and with staging 82a5df02's template: the texts
#       (normText) are equal; drawings 3, 3, 3, 1.
#   V1: the combined figure's list of replaced panels written by the head's writer (writeProtNormCompositePlaceholders)
#       from markers on pre_norm_pca and post_norm_correlation (and ruv_corrected_cancor where the lane ran RUV), and
#       carried into QC_figures as Copy does: exactly one added sentence (LFQ, TMT, DIA: letters a and h, the lanes skip
#       RUV); DIA with limpa: no combined figure, so no list and no sentence.
#   V2: a list in QC_figures with no combined figure: no sentence.
# Nothing is asserted: checks.tsv records predicted against observed.
# Usage (run-tests scripts job, package root): Rscript docs/qa/workflow-execution/kit/kit-259/render3.R

out <- normalizePath("kit259-render", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
work <- file.path(tempdir(), "kit259-render-work")
dir.create(work, showWarnings = FALSE)
say <- function(...) cat(..., "\n", sep = "")
writeUtf8 <- function(lines, path) writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
`%||%` <- function(a, b) if (!is.null(a)) a else b
BASE <- "82a5df02d"

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
ns <- asNamespace("MultiScholaR")
writeUtf8(c(paste("R", R.version.string), paste("HEAD", paste(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    collapse = " ")), paste("BASE", BASE)), file.path(out, "meta.txt"))

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

# --- Predictions (written before the run) ------------------------------------------------------------------------
SENTENCE <- paste("In Figure 1, panels a (principal component analysis after log2 transformation) and h (Pearson",
    "correlation after between-sample normalisation) could not be drawn by the application; each shows a note giving",
    "the reason instead of the plot.")
lanes <- list(
    lfq = list(id = "prot_lfq", template = "LFQ_report.rmd", ruv = "skip", drawings = 3L, sentence = TRUE),
    tmt = list(id = "prot_tmt", template = "TMT_report.rmd", ruv = "skip", drawings = 3L, sentence = TRUE),
    dia = list(id = "prot_dia", template = "DIANN_report.rmd", ruv = "skip", drawings = 3L, sentence = TRUE),
    dial = list(id = "prot_dia_limpa", template = "DIANN_limpa_report.rmd", ruv = "covariates", drawings = 1L,
        sentence = FALSE))

rows <- list()
record <- function(lane, variant, check, predicted, observed) {
    rows[[length(rows) + 1L]] <<- data.frame(lane = lane, variant = variant, check = check, predicted = predicted,
        observed = observed, verdict = if (identical(as.character(predicted), as.character(observed))) "as predicted" else
            "DIFFERS", stringsAsFactors = FALSE)
}

# --- The lanes, captured when the report exists (KIT-236's mechanism) ---------------------------------------------
helperFn <- function(name) if (exists(name, envir = ns, inherits = FALSE)) get(name, envir = ns) else get(name, envir = globalenv())
runLane <- helperFn("guiEntryRunLane")
store_pattern <- "(^|/)\\.multischolar(/|$)|(^|/)\\.workflow|\\.duckdb|\\.qs2$|\\.wal$"
cap <- new.env(parent = emptyenv())
cap$lane <- NULL
cap$base <- list()
cap$dirs <- list()
kit259RenderCapture <- function(base_dir) {
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
    trace("guiEntryArtifacts", tracer = bquote(.(kit259RenderCapture)(base_dir)), where = runLane, print = FALSE)
    TRUE
}, error = function(e) { say("trace failed: ", conditionMessage(e)); FALSE })
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
templateText <- function(file, side) {
    if (identical(side, "head")) return(readLines(file.path("inst", "reports", "proteomics", file), warn = FALSE, encoding = "UTF-8"))
    system2("git", c("show", paste0(BASE, ":inst/reports/proteomics/", file)), stdout = TRUE)
}

#' Copy the lane's captured project, let `edit(p)` change it, render with the side's template; the report's text.
renderVariant <- function(short, vid, side, edit = function(p) NULL) {
    lane <- lanes[[short]]
    root <- file.path(tempdir(), "kit259-variants", short, vid, side)
    unlink(root, recursive = TRUE)
    dir.create(root, recursive = TRUE)
    file.copy(file.path(work, short, "base"), root, recursive = TRUE)
    pd <- mapPaths(cap$dirs[[short]], cap$base[[short]], file.path(root, "base"))
    p <- pd[["proteomics"]]
    for (old in list.files(p$results_summary_dir, pattern = "\\.docx$", full.names = TRUE)) {
        if (!grepl("analysis_record", basename(old), fixed = TRUE)) unlink(old)
    }
    writeUtf8(templateText(lane$template, side), file.path(p$base_dir, "scripts", "proteomics", lane$template))
    edit(p)
    assign("project_dirs", pd, envir = globalenv())
    rendered <- RenderReport(omic_type = "proteomics", experiment_label = lane$id, rmd_filename = lane$template)
    if (is.null(rendered) || !file.exists(rendered)) stop("render failed")
    report <- readReport(rendered)
    writeUtf8(c(paste("drawings:", report$drawings), report$raw), file.path(out, paste0(short, "-", vid, "-", side, ".txt")))
    report
}

markAndList <- function(short) function(p) {
    lane <- lanes[[short]]
    figures <- c("pre_norm_pca", "post_norm_correlation", if (!identical(lane$ruv, "skip")) "ruv_corrected_cancor")
    for (f in figures) writeLines("injected", file.path(p$feature_qc_dir, paste0(f, ".png.placeholder")))
    written <- writeProtNormCompositePlaceholders(lane$ruv, p$feature_qc_dir, "proteomics")
    if (!is.null(written)) file.copy(written, file.path(p$qc_figures_dir, basename(written)), overwrite = TRUE)
}

for (short in names(lanes)) {
    lane <- lanes[[short]]
    if (is.null(cap$dirs[[short]])) {
        record(short, "lane", "captured", TRUE, FALSE)
        next
    }
    tryCatch({
        head0 <- renderVariant(short, "V0", "head")
        base0 <- renderVariant(short, "V0", "base")
        record(short, "V0", "text equal to staging's template", TRUE, identical(head0$norm, base0$norm))
        record(short, "V0", "drawings", lane$drawings, head0$drawings)
        record(short, "V0", "no placeholder sentence", FALSE, grepl("could not be drawn by the application", head0$norm, fixed = TRUE))
        v1 <- renderVariant(short, "V1", "head", markAndList(short))
        count <- lengths(regmatches(v1$norm, gregexpr("could not be drawn by the application", v1$norm, fixed = TRUE)))
        record(short, "V1", "sentences", if (lane$sentence) 1L else 0L, count)
        if (lane$sentence) {
            record(short, "V1", "the predicted sentence (letters a and h)", TRUE, grepl(normText(SENTENCE), v1$norm, fixed = TRUE))
            record(short, "V1", "text otherwise equal to V0", TRUE,
                identical(normText(sub(normText(SENTENCE), "", v1$norm, fixed = TRUE)), head0$norm))
        } else {
            record(short, "V1", "text equal to V0", TRUE, identical(v1$norm, head0$norm))
        }
        v2 <- renderVariant(short, "V2", "head", function(p) {
            writeLines("h\tpost_norm_correlation", file.path(p$qc_figures_dir, "proteomics_composite_QC_figure.placeholders"))
            unlink(file.path(p$qc_figures_dir, paste0("proteomics_composite_QC_figure.", c("png", "pdf"))))
        })
        record(short, "V2", "no placeholder sentence", FALSE, grepl("could not be drawn by the application", v2$norm, fixed = TRUE))
    }, error = function(e) record(short, "render", "error", "none", paste("error:", conditionMessage(e))))
}

writeUtf8(c(paste("trace installed:", trace_ok), unname(status)), file.path(out, "lanes.txt"))
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say(sum(table$verdict == "DIFFERS"), " of ", nrow(table), " checks differ from the prediction")
