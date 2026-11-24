# KIT-203 G0 fixtures: the base commit's outputs that S3, S4 and S6 compare head against. Evidence code under docs/qa:
# no package code. Runs on a GitHub-hosted runner (run-tests `scripts`) on the base code (installed), from the package
# root, after gate-k203.R g0 has written the E0 accessions. Each capture runs inside testthat::test_that() (edition 3,
# the test context of the files that read them; KIT-199's form). No service is called (the UniProt download stubbed).
#
# Usage:
#   Rscript docs/qa/workflow-execution/kit/kit-203/capture-k203.R <e0 accessions tsv> <out dir>
#
# Writes to <out dir> (committed under tests/testdata/enrichment-significance/ before the tests):
#   e0-accessions.tsv            the 126 accessions of KIT-173's null draw 5 (copied)
#   s3-none-up.tsv, .rds         that list through processEnrichments() on the full E. coli fixture, go_propagation
#                                "none": the written up TSV, and list(result = @result, plot_inputs)
#   t7-<value>-<up|down>.tsv, .rds   KIT-173's T7 annotation fixture with T12's table (p_term up, d_term down), for
#                                <value> in none, is_a_part_of_regulates (value "iapor" in the file name)
#   s6-identity.rds              the T7 project (helper goPropRunProject, real kernel, the tab's default propagation):
#                                the product identity, declared names and values, the outcome names
#   capture.tsv                  commit, versions, checks, and every file's SHA-256
# The S3 and T7 captures run with the kernel's plot files stubbed (no assertion reads them); the project capture without.

args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
e0_path <- args[[1L]]
out <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(MultiScholaR))
suppressPackageStartupMessages(library(testthat))
testthat::local_edition(3)
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
# The fixture directory from the package root, whatever testthat's test_path() would resolve in this context.
h$goPropFixtureDir <- function(...) file.path("tests", "testdata", "go-propagation", ...)
restore_widget <- h$goPropReplaceBinding("MultiScholaR", "saveSelfContainedWidget", function(...) invisible(NULL))
restore_ggsave <- h$goPropReplaceBinding("ggplot2", "ggsave", function(...) invisible(NULL))
checks <- list()

#' The S3 table: every entry of the E. coli fixture, the list at q = 0.01 and log2FC 1, the rest 0.5 and 0
#' (KIT-189's runRoute); the same function is in test-enrichment-significance.R.
s3Run <- function(list_ids, pathway_dir) {
    raw <- h$goPropReadTsv("uniprot-ecoli-83333.tsv.gz")
    annotations <- h$goPropConvert(raw)
    entries <- sort(unique(raw$Entry), method = "radix")
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = ifelse(entries %in% list_ids, 0.01, 0.5),
        log2FC = ifelse(entries %in% list_ids, 1, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "effect"),
        da_data = list(effect = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    processEnrichments(da, "999999", q_cutoff = 0.05, pathway_dir = pathway_dir, go_annotations = annotations,
        exclude_iea = FALSE, protein_id_column = "Protein.IDs", contrast_names = "effect", go_propagation = "none")
}

#' The T7 table (KIT-173 T12): the T7 annotation fixture, p_term proteins 3 up, d_term 3 down, the rest 0
t7Run <- function(value, pathway_dir) {
    roles <- h$goPropReadTsv("t7-fixture-entries.tsv")
    direct <- h$goPropConvert(h$goPropReadTsv("t7-annotation-ecoli.tsv"))
    entries <- sort(direct$Entry)
    table <- data.frame(Protein.IDs = entries, fdr_qvalue = 0.01,
        log2FC = ifelse(entries %in% roles$entry[roles$role == "p_term"], 3,
            ifelse(entries %in% roles$entry[roles$role == "d_term"], -3, 0)), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "t7"),
        da_data = list(t7 = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    processEnrichments(da, "999999", q_cutoff = 0.05, pathway_dir = pathway_dir, go_annotations = direct,
        exclude_iea = FALSE, protein_id_column = "Protein.IDs", contrast_names = "t7", go_propagation = value)
}

keep <- function(result, tsv_from, stem) {
    file.copy(tsv_from, file.path(out, paste0(stem, ".tsv")), overwrite = TRUE)
    saveRDS(list(result = result@result, plot_inputs = attr(result, "enrichment_plot_inputs", exact = TRUE)),
        file.path(out, paste0(stem, ".rds")), version = 3L)
}

testthat::test_that("KIT-203 G0 capture: S3 list, go_propagation none", {
    testthat::local_test_context()
    list_ids <- utils::read.delim(e0_path, colClasses = "character")$entry
    file.copy(e0_path, file.path(out, "e0-accessions.tsv"), overwrite = TRUE)
    pathway_dir <- withr::local_tempdir()
    results <- NULL
    utils::capture.output(results <- suppressMessages(suppressWarnings(s3Run(list_ids, pathway_dir))))
    up <- results@enrichment_data$effect$up
    keep(up, file.path(pathway_dir, "effect_up_enrichment_results.tsv"), "s3-none-up")
    tsv <- utils::read.delim(file.path(out, "s3-none-up.tsv"), colClasses = "character")
    checks$s3_list_size <<- length(list_ids)
    checks$s3_tested <<- nrow(tsv)
    checks$s3_significant <<- sum(tsv$significant %in% "TRUE")
    checks$s3_m <<- length(up@geneSets)
    testthat::expect_identical(checks$s3_significant, 120L)
})

testthat::test_that("KIT-203 G0 capture: T7 fixture, both propagation values", {
    testthat::local_test_context()
    for (value in c("none", "is_a_part_of_regulates")) {
        pathway_dir <- withr::local_tempdir()
        results <- NULL
        utils::capture.output(results <- suppressMessages(suppressWarnings(t7Run(value, pathway_dir))))
        short <- if (identical(value, "none")) "none" else "iapor"
        for (direction in c("up", "down")) {
            result <- results@enrichment_data$t7[[direction]]
            testthat::expect_false(is.null(result), info = paste(value, direction))
            keep(result, file.path(pathway_dir, paste0("t7_", direction, "_enrichment_results.tsv")),
                paste("t7", short, direction, sep = "-"))
            checks[[paste("t7", short, direction, "rows", sep = "_")]] <<- nrow(result@result)
        }
    }
})

testthat::test_that("KIT-203 G0 capture: the T7 project's identity (S6)", {
    testthat::local_test_context()
    # The project publishes the DA's interactive volcano before the enrichment: the plot-file stubs of the two captures
    # above are lifted first (VAL-54 instrument repair of run 45326277295: they made the publication refuse).
    restore_widget()
    restore_ggsave()
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    work <- withr::local_tempdir()
    annotations <- h$goPropAnnotate(h$goPropT7Raw(), file.path(work, "cache"))
    result <- NULL
    utils::capture.output(result <- suppressMessages(suppressWarnings(h$goPropRunProject(file.path(work, "project"),
        annotations))))
    record <- jsonlite::read_json(result$record, simplifyVector = FALSE)
    step <- Filter(function(s) identical(s$operation, "enrichment"), record$steps)[[1L]]
    saveRDS(list(identity = result$run$reference$identity, declared_names = vapply(step$parameters, `[[`, "", "name"),
        parameters = result$run$parameters, outcome_names = names(step$outcomes)),
        file.path(out, "s6-identity.rds"), version = 3L)
    checks$s6_declared_names <<- paste(vapply(step$parameters, `[[`, "", "name"), collapse = ",")
    testthat::expect_false(is.null(result$run$reference$identity))
})

restore_widget()
restore_ggsave()
files <- sort(setdiff(list.files(out), "capture.tsv"))
values <- c(commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)[[1L]], error = function(e) NA),
    r = R.version.string,
    vapply(c("MultiScholaR", "clusterProfiler", "enrichit", "qvalue", "GO.db", "testthat"),
        function(p) as.character(utils::packageVersion(p)), ""),
    vapply(checks, function(x) paste(as.character(x), collapse = ","), ""),
    stats::setNames(vapply(files, function(f) digest::digest(file = file.path(out, f), algo = "sha256"), ""),
        paste0("sha256_", files)))
utils::write.table(data.frame(measure = names(values), value = unname(values)), file.path(out, "capture.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
print(values)
