# KIT-174 A3: the proteomics normalisation tab's options column (which holds the RUV section) rendered
# with the base commit's proteomics configuration and with the head commit's, and their line diff.
# The renderer (R/mod_omic_norm_panels.R) must be identical on base and head; the base configuration
# is R/utils_omic_gui_config.R at the base commit, sourced over the installed head namespace.
# Usage: Rscript norm-panel.R <base sha> <out dir>
args <- commandArgs(trailingOnly = TRUE)
base <- args[[1L]]
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
renderer_same <- system2("git", c("diff", "--quiet", base, "HEAD", "--", "R/mod_omic_norm_panels.R"))
if (!identical(renderer_same, 0L)) stop("R/mod_omic_norm_panels.R differs between base and head")
render <- function(config) {
    column <- MultiScholaR:::omicNormOptionsColumn(shiny::NS("normalization"), config$norm)
    strsplit(as.character(htmltools::renderTags(column)$html), "\n", fixed = TRUE)[[1L]]
}
head_lines <- render(MultiScholaR:::omicGuiConfig("proteomics"))
source_file <- tempfile(fileext = ".R")
writeLines(system2("git", c("show", paste0(base, ":R/utils_omic_gui_config.R")), stdout = TRUE),
    source_file)
base_env <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(source_file, envir = base_env)
base_lines <- render(base_env$omicGuiConfig("proteomics"))
base_file <- file.path(out, "norm-panel-base.html")
head_file <- file.path(out, "norm-panel-head.html")
writeLines(base_lines, base_file)
writeLines(head_lines, head_file)
difference <- suppressWarnings(system2("diff", c(base_file, head_file), stdout = TRUE))
writeLines(difference, file.path(out, "norm-panel.diff"))
changed <- grep("^[<>]", difference, value = TRUE)
cat(sprintf("base lines %d, head lines %d, changed lines %d\n", length(base_lines), length(head_lines),
    length(changed)))
cat(changed, sep = "\n")
