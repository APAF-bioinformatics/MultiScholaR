# KIT-266 G0, memory measurement: ONE fresh R process builds ONE combined QC figure with the existing builder
# (generateProtNormCompositeQcFigure) from the lane's stage images, and prints /proc/self/status readings verbatim.
# The rise is defined as KIT-122 measured it: VmHWM after the build minus VmRSS just before it.
# Usage: Rscript memory-build.R <three|two|control> <image dir>
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1L]]
imgdir <- args[[2L]]
statusLines <- function() {
    lines <- readLines("/proc/self/status")
    lines[grepl("^(VmPeak|VmSize|VmHWM|VmRSS):", lines)]
}
fieldKb <- function(lines, field) as.numeric(gsub("[^0-9]", "", lines[startsWith(lines, paste0(field, ":"))]))
pngDim <- function(path) {
    b <- as.integer(readBin(path, "raw", 24L))
    w <- c(16777216, 65536, 256, 1)
    c(sum(b[17:20] * w), sum(b[21:24] * w))
}
cat("mode: ", mode, "\n", sep = "")

if (identical(mode, "control")) {
    # The measure itself: a 240 MB vector, touched, must show as a rise of about 0.22 GiB.
    before <- statusLines()
    x <- rep(1, 3e7)
    keep <- sum(x)
    after <- statusLines()
    cat("BEFORE\n", paste(before, collapse = "\n"), "\nAFTER\n", paste(after, collapse = "\n"), "\n", sep = "")
    rise <- (fieldKb(after, "VmHWM") - fieldKb(before, "VmRSS")) / 1024^2
    cat(sprintf("MEMORY_CONTROL rise_GiB=%.3f (allocated 3e7 doubles = 0.224 GiB) verdict=%s\n", rise,
        if (rise >= 0.18 && rise <= 0.60) "PASS" else "BLIND"))
    quit(save = "no", status = 0L)
}

suppressMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")
work <- file.path(tempdir(), "qc")
dir.create(work, recursive = TRUE)
file.copy(list.files(imgdir, full.names = TRUE), work)
expected <- c(paste0(rep(c("pre_norm", "post_norm", "ruv_corrected"), each = 4), "_", c("pca", "density", "rle", "correlation"), ".png"),
    "ruv_covariates_cancor.png")
cat("images present: ", paste(sort(list.files(work)), collapse = " "), "\n", sep = "")
cat("images missing: ", paste(setdiff(expected, list.files(work)), collapse = " "), "\n", sep = "")
msgs <- character()
msgFn <- function(...) {
    msgs <<- c(msgs, paste(as.character(c(...)), collapse = " "))
    invisible(NULL)
}
# The two-line local resolver of the three-column build: the default resolver, with the limpa cancor file name.
limpaResolver <- function(ruvMode, qcDir) {
    r <- ns$resolveProtNormCompositeFigureInputs(ruvMode, qcDir)
    r$plotFiles <- sub("ruv_corrected_cancor.png", "ruv_covariates_cancor.png", r$plotFiles, fixed = TRUE)
    r
}
ruvMode <- if (identical(mode, "three")) "covariates" else "skip"
resolver <- if (identical(mode, "three")) limpaResolver else ns$resolveProtNormCompositeFigureInputs
before <- statusLines()
started <- Sys.time()
ns$generateProtNormCompositeQcFigure(ruvMode = ruvMode, qcDir = work, omicType = "proteomics",
    resolveInputsFn = resolver, messageFn = msgFn)
secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
after <- statusLines()
cat("BEFORE\n", paste(before, collapse = "\n"), "\nAFTER\n", paste(after, collapse = "\n"), "\n", sep = "")
cat("builder messages:\n", paste(msgs, collapse = "\n"), "\n", sep = "")
png <- file.path(work, "proteomics_composite_QC_figure.png")
cat(sprintf("output files: %s\n", paste(sort(grep("composite", list.files(work), value = TRUE)), collapse = " ")))
if (file.exists(png)) cat(sprintf("PNG %d x %d px\n", pngDim(png)[[1L]], pngDim(png)[[2L]])) else cat("PNG missing\n")
failed <- any(grepl("Could not generate|Error", msgs)) || !file.exists(png)
rise <- (fieldKb(after, "VmHWM") - fieldKb(before, "VmRSS")) / 1024^2
cat(sprintf("MEMORY_READING mode=%s VmRSS_before_kB=%.0f VmHWM_before_kB=%.0f VmRSS_after_kB=%.0f VmHWM_after_kB=%.0f build_seconds=%.1f\n",
    mode, fieldKb(before, "VmRSS"), fieldKb(before, "VmHWM"), fieldKb(after, "VmRSS"), fieldKb(after, "VmHWM"), secs))
cat(sprintf("MEMORY_RISE mode=%s rise_GiB=%.3f (VmHWM after minus VmRSS before) build_ok=%s\n", mode, rise, !failed))
