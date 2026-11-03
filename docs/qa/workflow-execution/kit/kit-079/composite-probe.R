# One process: the composite rendered from a member's published per-assay files with the shipped generator, reading the
# PNGs as doubles or as native rasters. Usage: composite-probe.R <lipid|metab> <qc dir> <out dir> <double|native>
suppressMessages(library(MultiScholaR))
a <- commandArgs(TRUE); omic <- a[[1]]; qc <- a[[2]]; out <- a[[3]]; mode <- a[[4]]
ns <- asNamespace("MultiScholaR"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
hwm <- function() as.numeric(gsub("[^0-9]", "", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) / 1024^2
reader <- if (identical(mode, "native")) function(path) png::readPNG(path, native = TRUE) else png::readPNG
assays <- unique(sub("_(pre_norm|post_norm|ruv_corrected)_(pca|density|rle|correlation)[.]png$", "", grep("_pca[.]png$", list.files(qc), value = TRUE)))
files <- character(); labels <- list(); i <- 1L
for (assay in assays) for (type in c("pca", "density", "rle", "correlation")) {
  files <- c(files, file.path(qc, sprintf("%s_%s_%s.png", assay, c("pre_norm", "post_norm", "ruv_corrected"), type)))
  labels[[paste(assay, type, sep = "_")]] <- sprintf("%s)", letters[i + 0:2]); i <- i + 3L
}
base <- hwm()
if (omic == "lipid") {
  g <- ns$buildLipidNormCompositeFromFilesGenerator(readPngFn = reader)
} else {
  g <- function(...) ns$generateMetabNormCompositeFromFiles(..., loadImageAsPlotFn = function(path, ...) ns$loadMetabNormImageAsPlot(path, readPngFn = reader))
}
res <- g(plot_files = files, ncol = 3, row_labels = labels, column_labels = c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected"))
built <- hwm()
ns$savePlot(res$plot, out, "composite", width = res$width, height = res$height, dpi = 150, limitsize = FALSE)
done <- hwm()
cat(sprintf("%s %s: %d files; peak GiB after load %.2f, after build %.2f, after render %.2f; png %s %.1f MB; pdf %.1f MB; rds %.1f MB\n", omic, mode, length(files), base, built, done,
  substr(digest::digest(file = file.path(out, "composite.png"), algo = "sha256"), 1, 12), file.size(file.path(out, "composite.png"))/1e6, file.size(file.path(out, "composite.pdf"))/1e6, file.size(file.path(out, "composite.rds"))/1e6))
