# Shared-floor probe (2026-10-09, read-only; not a counted run).
# Measures the anonymous memory one fresh R process carries after loading a set of namespaces and one gc()
# Each namespace is loaded by its own top-level call, as the managed processes do (`Rscript -e`);
# the same loads inside a `for` loop measured about 130 MiB higher per process.
# Every process in a managed slice (controller, monitor, worker) loads the full package namespace
# before any workflow work starts, so three times the full figure is the slice floor.
# Usage: CORE_LIBRARY=<installed library> Rscript floor-probe.R   (from the repository root)
lib <- Sys.getenv("CORE_LIBRARY")
measure <- function(pkgs, reps = 3L) {
    script <- tempfile(fileext = ".R")
    writeLines(c(sprintf(".libPaths(c(%s, .libPaths()))", deparse(lib)),
        sprintf("invisible(suppressPackageStartupMessages(loadNamespace(%s)))", vapply(pkgs, deparse, "")),
        "invisible(gc()); cat(as.numeric(gsub('[^0-9]', '', grep('^RssAnon', readLines('/proc/self/status'), value = TRUE))) / 1024)"), script)
    on.exit(unlink(script))
    stats::median(vapply(seq_len(reps), function(i) as.numeric(system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", script), stdout = TRUE, stderr = FALSE)), numeric(1)))
}
d <- read.dcf(file.path(lib, "MultiScholaR", "DESCRIPTION"))
imports <- trimws(gsub("\\(.*\\)", "", strsplit(paste(d[, "Imports"], d[, "Depends"], sep = ","), ",")[[1]]))
imports <- setdiff(imports[nzchar(imports)], c("R", "methods", "stats", "utils", "graphics", "grDevices"))
app_only <- c("Glimma", "cowplot", "ComplexHeatmap", "ggraph", "Biobase", "seqinr", "git2r", "devtools", "circlize", "GGally", "ggpubr")
worker_set <- c("arrow", "duckdb", "DBI", "digest", "cli", "jsonlite", "dplyr", "fastduckdb", "rlang", "tibble", "data.table", "later")
m <- c(bare_r = measure(character(0)), full_namespace = measure("MultiScholaR"),
    imports_without_app_only = measure(setdiff(imports, app_only)), worker_style_set = measure(worker_set))
cat(sprintf("imports=%d app_only=%d (median of 3 fresh processes each)\n", length(imports), length(app_only)))
cat("RESULT anonymous MiB per process:\n"); print(round(m))
cat(sprintf("RESULT2 three-process floor GiB: full=%.2f without_app_only=%.2f controller_full_plus_two_worker_style=%.2f\n",
    3 * m[["full_namespace"]] / 1024, 3 * m[["imports_without_app_only"]] / 1024, (m[["full_namespace"]] + 2 * m[["worker_style_set"]]) / 1024))
