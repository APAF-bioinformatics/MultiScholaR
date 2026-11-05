# KIT-094: the characterisation snapshot of the lists the bindings take over, for what the KIT-073 frozen dump does not
# cover (the report sources, the type allowlists and routes, the replay's input roles, the fixture records' steps,
# parameters and outcomes, and the registry fields this ticket touches). Run on the base commit before any edit; the
# functions are the ones the test calls (tests/testthat/test-workflow-omics-binding-steps.R), read from that file so
# there is one definition.
#
# Usage: Rscript characterise.R <package root> <output .rds>
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
suppressMessages(pkgload::load_all(args[[1L]], quiet = TRUE, export_all = TRUE, helpers = TRUE))
test_file <- file.path(args[[1L]], "tests", "testthat", "test-workflow-omics-binding-steps.R")
definitions <- new.env(parent = asNamespace("MultiScholaR"))
# Only the file's top-level assignments: its `test_that` blocks are not run here.
for (expression in parse(test_file, keep.source = FALSE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-"))) eval(expression, definitions)
}
root <- tempfile("kit094-characterise-")
dir.create(root)
snapshot <- suppressWarnings(suppressMessages(definitions$kit094Characterisation(root)))
saveRDS(snapshot, args[[2L]])
cat(sprintf("KIT094_CHARACTERISATION %s parts=%s sha256=%s\n", args[[2L]], paste(names(snapshot), collapse = ","),
    digest::digest(file = args[[2L]], algo = "sha256")))
