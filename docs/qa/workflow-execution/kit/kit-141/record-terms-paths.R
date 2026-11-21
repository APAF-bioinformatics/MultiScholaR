# KIT-141, granted investigation after run 44861344601: the new workbook test found no "Term or pathway" column in the
# enrichment terms sheet written from the built (in-memory) record. This script builds the same enriched project as
# that test (the test file's own fixtures: the tier-1 corpus, the stand-in trainer, the STRING double, no network) and
# prints, for the built record and for the JSON round-tripped record (as integrationWriteRecordBundle reads it):
# str() of the first two rows of block$enrichment[[1]]$terms, the class of the hydrated value$terms that
# integrationRecordEnrichmentBlock reads, and the terms table integrationRecordEnrichmentTables returns; then the
# terms sheet of a workbook written from the built record, read back with openxlsx's default column names and with
# sep.names = " ".
#
# Usage (from the package root, head installed): Rscript record-terms-paths.R <out dir>
# Writes <out dir>/record-terms-paths.txt; exits 1 when the two paths give different terms tables.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("MultiScholaR", quietly = TRUE), requireNamespace("testthat", quietly = TRUE))
out_file <- file.path(out_dir, "record-terms-paths.txt")
sink(out_file, split = TRUE)

# The test file's fixtures, evaluated as testthat would (package namespace as parent); its tests are not run.
ns <- asNamespace("MultiScholaR")
env <- new.env(parent = ns)
# Each helper on its own: one that defines S4 classes cannot be sourced into a child of the locked namespace (run
# 44863659582); it is skipped and named, and the fixtures used here do not need it.
for (helper in list.files("tests/testthat", pattern = "^helper.*\\.[rR]$", full.names = TRUE)) {
    tryCatch(sys.source(helper, envir = env), error = function(error)
        cat(sprintf("helper skipped: %s (%s)\n", basename(helper), conditionMessage(error))))
}
exprs <- parse("tests/testthat/test-integration-enrichment.R", keep.source = FALSE)
for (e in exprs) {
    if (is.call(e) && identical(e[[1L]], as.name("test_that"))) next
    eval(e, env)
}

same <- local({
    paths <- env$intEnrProject()
    fixture <- env$intEnrModel(paths, close = FALSE)
    manager <- fixture$manager
    local({
        env$intEnrLocalString()
        ns$runIntegrationEnrichment(manager, list(factor = "Factor1", view = "proteomics", database = "string"),
            api_key = env$SENTINEL_KEY, query = env$intEnrQuery)
    })
    export <- manager$.context$snapshot$exports$integration_mofa_model
    entry <- export$enrichment[[1L]]
    value <- ns$withWorkflowStateScope(manager, function() ns$hydrateSavedWorkflowEntity(manager$.context$store,
        entry$reference)$value)
    cat("== class of value$terms (hydrated in the owner's scope)\n")
    print(class(value$terms))
    unscoped <- tryCatch(class(ns$hydrateSavedWorkflowEntity(manager$.context$store, entry$reference)$value$terms),
        error = function(error) paste("error:", conditionMessage(error)))
    cat("== class of value$terms (hydrated as integrationRecordEnrichmentBlock does, no scope)\n")
    print(unscoped)
    built <- ns$buildWorkflowAnalysisRecord(manager)
    json <- file.path(paths$base_dir, "record.json")
    jsonlite::write_json(built, json, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE)
    read <- jsonlite::read_json(json, simplifyVector = FALSE)
    tables <- list()
    for (side in c("built", "json")) {
        record <- if (side == "built") built else read
        block <- record$analyses$integration_mofa_model$integration
        cat(sprintf("\n==== %s record\n== enrichment entries: %d; complete: %s\n", side, length(block$enrichment),
            paste(vapply(block$enrichment, \(x) paste(unlist(x$complete), collapse = ""), ""), collapse = ",")))
        cat("== str(block$enrichment[[1]]$terms[1:2])\n")
        str(utils::head(block$enrichment[[1L]]$terms, 2L))
        terms <- ns$integrationRecordEnrichmentTables(block)$terms
        cat("== integrationRecordEnrichmentTables(block)$terms: names, then rows\n")
        print(names(terms))
        print(utils::head(terms, 4L))
        tables[[side]] <- terms
    }
    cat("\n== identical(built terms table, json terms table):", identical(tables$built, tables$json), "\n")
    xlsx <- file.path(paths$base_dir, "built.xlsx")
    ns$writeWorkflowAnalysisRecordWorkbook(built, xlsx)
    cat("== terms sheet of the workbook from the built record, read.xlsx(check.names = FALSE): names\n")
    print(names(openxlsx::read.xlsx(xlsx, sheet = "integration_enrichment_terms", check.names = FALSE)))
    cat("== the same with sep.names = \" \": names, then column 'Term or pathway'\n")
    sheet <- openxlsx::read.xlsx(xlsx, sheet = "integration_enrichment_terms", check.names = FALSE, sep.names = " ")
    print(names(sheet))
    print(sheet[["Term or pathway"]])
    ns$closeWorkflowStateOwner(manager)
    identical(tables$built, tables$json)
})
sink()
if (!isTRUE(same)) quit(status = 1L)
