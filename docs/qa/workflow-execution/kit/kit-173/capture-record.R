# KIT-173 T7(a): a record made on the base code, for the replay test T7(b).
#
# Usage, from the package root on a run-tests runner, after the commit under test is installed (in R1, the tests-only
# commit, whose R/ equals the base a60c4e9a):
#   Rscript docs/qa/workflow-execution/kit/kit-173/capture-record.R <out dir>
#
# One proteomics project on the clusterProfiler route (taxon 999999) with the real kernel (`processEnrichments`), on the
# committed E. coli fixture (tests/testdata/go-propagation/t7-annotation-ecoli.tsv, 70 proteins), the annotation table
# made by `getUniprotAnnotations` with the download stubbed (no service call). Writes to <out dir>:
#   project.tar     the project directory (store, source/uniprot_dat_cln.RDS, results), uncompressed tar (one file:
#                   the store's own paths would exceed the 100-byte path limit under tests/testdata)
#   record.json     its analysis record
#   enricher.tsv    per enricher call: SHA-256 of TERM2GENE and TERM2NAME as passed (helper goPropTableHash)
#   evidence.rds    the stored result tables, the plots' data digests, the pathway files' SHA-256 (TSV) or existence
#   product.rds     the enrichment product's parameters, service and contrasts
#   capture.tsv     commit, versions, declared names, and the files above with their SHA-256
# The outputs are committed under tests/testdata/go-propagation/replay-base/ before the fix commit.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[[1L]] else "k173-capture-out"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(out)
suppressPackageStartupMessages(library(MultiScholaR))
Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)

work <- file.path(tempdir(), "k173-capture")
dir.create(work, recursive = TRUE)
annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
root <- file.path(work, "project")
acc <- new.env()
restore <- h$goPropRecordEnricher(acc)
utils::capture.output(result <- suppressMessages(suppressWarnings(h$goPropRunProject(root, annotations))))
restore()

evidence <- h$goPropRunEvidence(result$state, result$paths$pathway_dir)
saveRDS(evidence, file.path(out, "evidence.rds"))
saveRDS(list(parameters = result$product$parameters, service = result$product$service,
    contrasts = result$product$contrasts, provenance_names = names(result$provenance$values)),
    file.path(out, "product.rds"))
calls <- do.call(rbind, lapply(seq_along(acc$calls), function(i) data.frame(call = i,
    term2gene_rows = nrow(acc$calls[[i]]$term2gene), term2gene_sha256 = acc$calls[[i]]$term2gene_sha256,
    term2name_sha256 = acc$calls[[i]]$term2name_sha256)))
utils::write.table(calls, file.path(out, "enricher.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
file.copy(result$record, file.path(out, "record.json"))
owd <- setwd(root)
utils::tar(file.path(out, "project.tar"), files = list.files(".", recursive = TRUE, all.files = TRUE),
    compression = "none", tar = "internal")
setwd(owd)

record <- jsonlite::read_json(result$record, simplifyVector = FALSE)
step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
sha <- function(f) digest::digest(file = file.path(out, f), algo = "sha256")
values <- c(
    commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]], error = function(e) NA),
    multischolar = as.character(utils::packageVersion("MultiScholaR")),
    clusterProfiler = as.character(utils::packageVersion("clusterProfiler")),
    go_db = as.character(utils::packageVersion("GO.db")),
    declared_names = paste(vapply(step$parameters, `[[`, "", "name"), collapse = ","),
    outcome_names = paste(names(step$outcomes), collapse = ","),
    enricher_calls = length(acc$calls),
    pathway_files = paste(evidence$files$file, collapse = ","),
    annotation_attribute = paste(names(attributes(annotations)), collapse = ","),
    sha256_project_tar = sha("project.tar"), sha256_record_json = sha("record.json"),
    sha256_evidence_rds = sha("evidence.rds"), sha256_product_rds = sha("product.rds"),
    sha256_enricher_tsv = sha("enricher.tsv"))
utils::write.table(data.frame(measure = names(values), value = unname(values)), file.path(out, "capture.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
print(values)
print(calls)
