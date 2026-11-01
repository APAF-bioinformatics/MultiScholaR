# KIT-073 Simulation Setup: freeze the observable surface of the lipidomics and metabolomics declarations on a tree
# (the base commit before the shared core): the registry entries and declaration fields, and both omics' wrapper
# journeys by content (tests/testthat/helper-workflow-omics-declarations.R's dump).
# Usage: R_LIBS=<kit library> Rscript docs/qa/workflow-execution/kit/kit-073/freeze-dump.R <tree> <dump helper> <out rds>
args <- commandArgs(TRUE); tree <- args[[1L]]; helper <- args[[2L]]; out <- args[[3L]]
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
sys.source(helper, envir = environment())
root <- tempfile("kit073-freeze-"); dir.create(root)
dump <- suppressWarnings(suppressMessages(omicsDeclarationsDump(root)))
saveRDS(dump, out)
cat(sprintf("frozen: %d lipid and %d metab registry entries; lipid states %d, metab states %d -> %s\n",
  length(dump$registry$lipidomics), length(dump$registry$metabolomics), length(dump$lipidomics$history),
  length(dump$metabolomics$history), out))
