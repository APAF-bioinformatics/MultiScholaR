# KIT-131 diagnostic (run-win-k131; build-time verification, not package code): run as a plain worker would be, by the
# package's own workflowLaunchProcess (same Rscript, arguments, environment, TMPDIR and working directory as the real
# worker of the same request), before that request's real worker starts. Prints PROBE lines:
#  - the process's environment, temporary folder, working directory, capabilities, bitmap type, devices available;
#  - ggsave of a lipid QC figure name at the path withWorkflowShortFigureFolder gives a worker (tempfile("wfx-")), and
#    the same in the controller's short temporary folder, with lengths, existence and conditions;
#  - for a lipid qc_figures request: the request's handler replayed (computeWorkflowOmicsAction) with the figure store
#    switched off (nothing saved, so the real worker is not served from the store), every figure traced
#    (figure-trace.R), then the replay's request publication folder is removed, so the real worker finds it as before.
# Usage (by lane-worker.R): Rscript --vanilla -e 'source(Sys.getenv("KIT131_PROBE"))' <request file> <hash>
`%||%` <- function(a, b) if (is.null(a)) b else a
arguments <- commandArgs(trailingOnly = TRUE)
request_path <- arguments[[1L]]
out <- function(...) cat("PROBE", sprintf(...), "\n")
out("pid %d; R %s; R.home %s", Sys.getpid(), as.character(getRversion()), R.home())
for (name in c("TMPDIR", "TMP", "TEMP", "R_LIBS", "MULTISCHOLAR_WORKFLOW_WORKER", "LANG", "R_USER", "HOME"))
    out("env %s = %s (%d)", name, Sys.getenv(name), nchar(Sys.getenv(name)))
out("tempdir %s (%d characters)", tempdir(), nchar(tempdir()))
out("working directory %s (%d characters)", getwd(), nchar(getwd()))
out("request %s (%d characters)", request_path, nchar(request_path))
out(".libPaths %s", paste(.libPaths(), collapse = " ; "))
capability <- capabilities()
out("capabilities %s", paste(names(capability), capability, sep = "=", collapse = " "))
out("bitmapType %s", getOption("bitmapType") %||% "NULL")
for (package in c("ragg", "Cairo", "png", "ggplot2", "svglite", "systemfonts", "textshaping"))
    out("package %s: %s", package, if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else "not available")
if (.Platform$OS.type == "windows") {
    value <- tryCatch(utils::readRegistry("SYSTEM\\CurrentControlSet\\Control\\FileSystem", "HLM")$LongPathsEnabled,
        error = function(error) conditionMessage(error))
    out("LongPathsEnabled as this process reads it: %s", paste(value, collapse = " "))
}
figure <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
device <- tryCatch(paste(deparse(ggplot2:::plot_dev("png", "x.png", dpi = 150))[1:3], collapse = " "), error = function(error) conditionMessage(error))
out("ggsave's png device: %s", substr(gsub("\\s+", " ", device), 1L, 300L))
attempt <- function(label, path, writer) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    message <- ""
    withCallingHandlers(tryCatch(writer(path), error = function(error) message <<- paste("error:", conditionMessage(error))),
        warning = function(warning) {
            message <<- paste(message, "warning:", conditionMessage(warning))
            invokeRestart("muffleWarning")
        })
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
    listed <- basename(path) %in% tryCatch(list.files(dirname(path)), error = function(error) character())
    out("write %s: %d characters, listed %s, bytes %s, folder exists %s; condition: %s; path %s", label, nchar(path), listed,
        suppressWarnings(file.size(path)), dir.exists(dirname(path)), substr(gsub("[\r\n\t]+", " ", trimws(message)), 1L, 300L), path)
}
writers <- list(
    ggsave = function(path) ggplot2::ggsave(path, figure, width = 8, height = 6, dpi = 150),
    agg_png = function(path) { ragg::agg_png(path); plot(1:3); grDevices::dev.off() },
    png_default = function(path) { grDevices::png(path); plot(1:3); grDevices::dev.off() },
    writeLines = function(path) writeLines("x", path))
scratch <- tempfile("wfx-")
for (writer in names(writers)) {
    extension <- if (writer == "writeLines") ".txt" else ".png"
    attempt(paste(writer, "in the worker's short folder (tempfile wfx-)"), file.path(scratch, "lipid_qc_dir",
        paste0("lcms_pos_pre_norm_correlation", extension)), writers[[writer]])
    short <- Sys.getenv("KIT131_SHORT")
    if (nzchar(short)) attempt(paste(writer, "in the controller's temporary folder"), file.path(short, "kit131-probe",
        "lipid_qc_dir", paste0("lcms_pos_pre_norm_correlation", extension)), writers[[writer]])
}

request <- tryCatch(readRDS(request_path), error = function(error) NULL)
operation <- request$operation %||% ""
out("operation %s; request folder %s (%d characters)", operation, request$directory %||% "", nchar(request$directory %||% ""))
if (grepl("^lipid_qc_figures$", operation)) {
    suppressPackageStartupMessages(library(MultiScholaR))
    ns <- asNamespace("MultiScholaR")
    source(file.path(Sys.getenv("KIT131_CHECKOUT"), "packaging", "verify", "figure-trace.R"))
    kit134FigureTraceInstall()
    .kit134FigureTrace$item <- "replay"
    utils::assignInNamespace("workflowOmicsFigureOwner", function(omic_type) NULL, "MultiScholaR")
    logged <- character()
    trace("log_warn", tracer = quote(cat("PROBE log_warn:", substr(paste(..., collapse = " "), 1L, 400L), "\n")),
        where = asNamespace("logger"), print = FALSE)
    replay <- tryCatch(withCallingHandlers(ns$computeWorkflowOmicsAction(request), warning = function(warning) {
        out("replay warning: %s", substr(gsub("[\r\n\t]+", " ", conditionMessage(warning)), 1L, 300L))
        invokeRestart("muffleWarning")
    }), error = function(error) error)
    out("replay: %s", if (inherits(replay, "error")) paste("error:", class(replay)[[1L]], substr(gsub("[\r\n\t]+", " ",
        conditionMessage(replay)), 1L, 700L)) else "returned")
    published <- file.path(request$directory, "p")
    files <- list.files(published, recursive = TRUE, full.names = TRUE)
    out("replay's request publication folder %s: %d files", published, length(files))
    for (file in files) out("  published %d characters %s bytes: %s", nchar(file), file.size(file), file)
    kit134FigureTraceReport(file.path(Sys.getenv("KIT131_OUT"), "probe-figures.tsv"), "replay")
    unlink(published, recursive = TRUE)
    out("replay's publication folder removed: %s", !dir.exists(published))
}
out("done")
