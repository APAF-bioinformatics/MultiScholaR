# KIT-077 Simulation Setup: the same plot rendered twice on the proof host, bytes compared. A ggplot with no randomness,
# each device the omics figures use. Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-077/figure-determinism.R
suppressMessages(library(ggplot2))
plot <- ggplot(data.frame(x = 1:20, y = (1:20)^2), aes(x, y)) + geom_point() + theme_bw()
same <- function(ext, render) {
    a <- tempfile(fileext = ext); b <- tempfile(fileext = ext)
    render(a); Sys.sleep(2); render(b)
    identical(readBin(a, "raw", file.size(a)), readBin(b, "raw", file.size(b)))
}
cat("ggsave png, 150 dpi, identical:", same(".png", function(p) ggsave(p, plot, width = 8, height = 6, dpi = 150)), "\n")
cat("ggsave png, 300 dpi, identical:", same(".png", function(p) ggsave(p, plot, width = 8, height = 6, dpi = 300)), "\n")
cat("ggsave pdf, identical:", same(".pdf", function(p) ggsave(p, plot, width = 8, height = 6)), "\n")
cat("png() device, identical:", same(".png", function(p) { grDevices::png(p, width = 10, height = 8, units = "in", res = 300); print(plot); grDevices::dev.off() }), "\n")
