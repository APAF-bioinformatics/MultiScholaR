# char-diff.R: where the current characterisation differs from the KIT-094 snapshot, after the KIT-095 exception (diagnostic).
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
source("tests/testthat/test-workflow-omics-binding-steps.R", local = TRUE, echo = FALSE, max.deparse.length = 0) |> tryCatch(error = function(e) NULL)
expected <- readRDS("tests/testthat/fixtures/kit094-characterisation.rds")
current <- suppressWarnings(suppressMessages(kit094Characterisation(tempfile())))
stripped <- kit095WithoutTranscriptomics(current, expected)
walk <- function(a, b, path) {
  if (identical(a, b)) return(invisible())
  if (is.list(a) && is.list(b) && identical(names(a), names(b)) && length(a)) { for (i in seq_along(a)) walk(a[[i]], b[[i]], paste0(path, "/", names(a)[i] %||% i)); return(invisible()) }
  cat("DIFF", path, ":", substr(paste(deparse(a), collapse = ""), 1, 150), "vs", substr(paste(deparse(b), collapse = ""), 1, 150), "\n")
}
for (part in names(expected)) walk(stripped[[part]], expected[[part]], part)
