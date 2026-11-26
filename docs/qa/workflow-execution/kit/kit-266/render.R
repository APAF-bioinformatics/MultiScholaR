out <- normalizePath("kit266-render", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
work <- file.path(tempdir(), "kit266-render-work")
dir.create(work, showWarnings = FALSE)
say <- function(...) cat(..., "\n", sep = "")
writeUtf8 <- function(lines, path) writeLines(enc2utf8(as.character(lines)), path, useBytes = TRUE)
`%||%` <- function(a, b) if (!is.null(a)) a else b
BASE <- "36e3180b0"  # staging with KIT-272: the V0 reference for every template (restated; the KIT-236 render holds the old bullets)

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


# --- Predictions (written before the run; ticket section 11 R3, restated on the merged tree) -----------------------
lanes <- list(
    lfq = list(id = "prot_lfq", template = "LFQ_report.rmd", drawings = 3L),
    tmt = list(id = "prot_tmt", template = "TMT_report.rmd", drawings = 3L),
    dia = list(id = "prot_dia", template = "DIANN_report.rmd", drawings = 3L),
    dial = list(id = "prot_dia_limpa", template = "DIANN_limpa_report.rmd", drawings = 2L))
CAPTION <- paste("Figure 1: QC metrics across normalisation stages. On DIA with limpa the first two columns are before the",
    "sample-correlation filter, and a third column, when the figure has one, is after it.")
LEFT <- "DPC-Quant's protein estimates (a complete log2 matrix; model-based where peptides were not detected)"
THREE <- paste0("Figure 1 shows four quality-control views at three stages of processing: ", LEFT,
    " (left column: panels a, d, g, j); after the peptide count filter, with no between-sample normalisation applied, before the ",
    "sample-correlation filter (centre column: b, e, h, k); and for the RUV-III-adjusted matrix after the sample-correlation filter ",
    "(drawn for quality control only; the differential abundance model carried the RUV-III factors as covariates) (right column: c, f, i, l).")
TWO <- paste0("Figure 1 shows four quality-control views at two stages of processing: ", LEFT,
    " (left column: panels a, c, e, g) and after the peptide count filter, with no between-sample normalisation applied, before the ",
    "sample-correlation filter (right column: b, d, f, h). RUV-III correction was not applied to this dataset.")
LIST_SENTENCE <- paste("In Figure 1, panels a (principal component analysis of DPC-Quant's protein estimates), f (principal component",
    "scores by group of the RUV-III-adjusted matrix after the sample-correlation filter) and m (RUV-III canonical correlation) could",
    "not be drawn by the application; each shows a note giving the reason instead of the plot.")
# KIT-272's bullets as the staging template prints them (the limpa branches copy them verbatim)
stagingBullets <- function() {
    lines <- system2("git", c("show", paste0(BASE, ":inst/reports/proteomics/DIANN_limpa_report.rmd")), stdout = TRUE)
    c(grep("^• Principal component analysis \\(panels a-c\\)", lines, value = TRUE),
      grep("^• Principal component scores by group \\(panels d-f\\)", lines, value = TRUE))
}

rows <- list()
record <- function(lane, variant, check, predicted, observed) {
    rows[[length(rows) + 1L]] <<- data.frame(lane = lane, variant = variant, check = check, predicted = as.character(predicted),
        observed = as.character(observed), verdict = if (identical(as.character(predicted), as.character(observed))) "as predicted" else
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


pngSize <- function(path) if (file.exists(path)) readBin(readBin(path, "raw", 24L)[17:24], "integer", n = 2L, size = 4L, endian = "big") else c(NA, NA)
limpaResolver <- function(mode, qc) ns$resolveDiaLimpaCompositeFigureInputs(mode, qc, "none")

for (short in names(lanes)) {
    lane <- lanes[[short]]
    if (is.null(cap$dirs[[short]])) {
        record(short, "lane", "captured", TRUE, FALSE)
        next
    }
    tryCatch({
        head0 <- renderVariant(short, "V0", "head")
        base0 <- renderVariant(short, "V0", "base")
        record(short, "V0", "drawings", lane$drawings, head0$drawings)
        if (!identical(short, "dial")) {
            record(short, "V0", "text equal to the staging template's (normText)", TRUE, identical(head0$norm, base0$norm))
        } else {
            record(short, "V0", "caption line present", TRUE, grepl(normText(CAPTION), head0$norm, fixed = TRUE))
            record(short, "V0", "Results paragraph is the three-stage covariates sentence (method none)", TRUE,
                grepl(normText(THREE), head0$norm, fixed = TRUE))
            record(short, "V0", "KIT-272's two bullets present", TRUE, all(vapply(stagingBullets(), function(b)
                grepl(normText(b), head0$norm, fixed = TRUE), NA)))
            record(short, "V0", "no 'protein intensities per sample' and no 'Density plots'", TRUE,
                !grepl("protein intensities per sample|Density plots", head0$norm))
            record(short, "V0", "staging's text had no Figure 1 caption line", FALSE, grepl(normText(CAPTION), base0$norm, fixed = TRUE))
            p <- mapPaths(cap$dirs[[short]], cap$base[[short]], file.path(work, short, "base"))[["proteomics"]]
            record(short, "V0", "QC_figures holds the composite PNG and PDF", TRUE,
                all(file.exists(file.path(p$qc_figures_dir, paste0("proteomics_composite_QC_figure.", c("png", "pdf"))))))
            file.copy(file.path(p$qc_figures_dir, "proteomics_composite_QC_figure.png"), file.path(out, "dial-V0-figure.png"), overwrite = TRUE)
            size <- pngSize(file.path(p$qc_figures_dir, "proteomics_composite_QC_figure.png"))
            record(short, "V0", "PNG width x height (three columns with m)", "1950x3900", paste0(size[[1L]], "x", size[[2L]]))
            # V1: markers on a, f and m, the list written through the limpa layout and carried into QC_figures as Copy does.
            v1 <- renderVariant(short, "V1", "head", function(p) {
                for (f in c("pre_norm_pca", "ruv_corrected_density", "ruv_covariates_cancor"))
                    writeLines("injected", file.path(p$feature_qc_dir, paste0(f, ".png.placeholder")))
                written <- ns$writeProtNormCompositePlaceholders("covariates", p$feature_qc_dir, "proteomics", resolveInputsFn = limpaResolver)
                if (!is.null(written)) file.copy(written, file.path(p$qc_figures_dir, basename(written)), overwrite = TRUE)
            })
            count <- lengths(regmatches(v1$norm, gregexpr("could not be drawn by the application", v1$norm, fixed = TRUE)))
            record(short, "V1", "one added sentence", 1L, count)
            record(short, "V1", "the predicted sentence (panels a, f, m; limpa words)", TRUE, grepl(normText(LIST_SENTENCE), v1$norm, fixed = TRUE))
            record(short, "V1", "text otherwise equal to V0", TRUE, identical(normText(sub(normText(LIST_SENTENCE), "", v1$norm, fixed = TRUE)), head0$norm))
            # V2: the project at the pass-through: no covariates step in the record, no RUV figures, the two-column figure.
            v2 <- renderVariant(short, "V2", "head", function(p) {
                unlink(file.path(p$feature_qc_dir, c(paste0("ruv_corrected_", c("pca", "rle", "density", "correlation"), ".png"), "ruv_covariates_cancor.png")))
                ns$generateProtNormCompositeQcFigure("skip", p$feature_qc_dir, "proteomics", resolveInputsFn = limpaResolver)
                for (ext in c("png", "pdf")) file.copy(file.path(p$feature_qc_dir, paste0("proteomics_composite_QC_figure.", ext)),
                    file.path(p$qc_figures_dir, paste0("proteomics_composite_QC_figure.", ext)), overwrite = TRUE)
                record_path <- file.path(p$base_dir, "scripts", "proteomics", "analysis_record.json")
                json <- jsonlite::read_json(record_path)
                json$steps <- Filter(function(s) !identical(s$operation, "protein_ruv_covariates"), json$steps)
                jsonlite::write_json(json, record_path, auto_unbox = TRUE, null = "null", digits = NA)
                file.copy(file.path(p$qc_figures_dir, "proteomics_composite_QC_figure.png"), file.path(out, "dial-V2-figure.png"), overwrite = TRUE)
                assign("v2size", pngSize(file.path(p$qc_figures_dir, "proteomics_composite_QC_figure.png")), envir = globalenv())
            })
            record(short, "V2", "PNG 1500 x 3300", "1500x3300", paste0(v2size[[1L]], "x", v2size[[2L]]))
            record(short, "V2", "the two-stage sentence", TRUE, grepl(normText(TWO), v2$norm, fixed = TRUE))
        }
    }, error = function(e) record(short, "render", "error", "none", paste("error:", conditionMessage(e))))
}

writeUtf8(c(paste("trace installed:", trace_ok), unname(status)), file.path(out, "lanes.txt"))
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say(sum(table$verdict == "DIFFERS"), " of ", nrow(table), " checks differ from the prediction")
