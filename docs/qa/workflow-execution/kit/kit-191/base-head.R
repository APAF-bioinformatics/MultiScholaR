# KIT-191 U7 (b) and U8's branch-point store (ticket 10; VAL-64 (2): built and used in one runner job, never kept).
# Base (the branch point) and head are each loaded from their own tree in their own fresh process. Each process:
#   - T1's fixture (KIT-183's `intEtPlatformTerms()`): camera's and gsea's terms;
#   - the corpus store (metabolomics_GCMS, the KEGG fixture table of `intEtCorpusSets()`, `intEtTrainer()`): camera's
#     and gsea's terms through `runIntegrationEnrichment()`.
# The base process keeps its store (camera and gsea entries) for the head process, which then (U8) replays it with
# `integrationReplayEnrichment()` (each entry with its own test, terms identical) and (U7 (b), AUDIT-KIT-183-impl I03's
# form) computes the camera product's identity at head on the base store's source with the recorded parameters: equal
# to the recorded identity in every component but `software`. The fixtures are KIT-183's, read from the head tree's
# tests/testthat/test-integration-enrichment-tests.R (functions named intEt*; nothing else of the file is run).
#
# Usage (from the head package root, full git history): Rscript base-head.R <base commit> <out dir>
# Writes <out dir>/base-head.tsv, replay.tsv, identity.tsv, verdict.txt and each process's log.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
work <- file.path(tempdir(), "k191-base-head")
dir.create(work)
rscript <- file.path(R.home("bin"), "Rscript")
base_tree <- file.path(work, "base")
stopifnot(system2("git", c("worktree", "add", "--detach", base_tree, base)) == 0L)
fixtures <- file.path(head, "tests/testthat/test-integration-enrichment-tests.R")

child <- file.path(work, "child.R")
writeLines(c(
    "a <- commandArgs(trailingOnly = TRUE)",
    "tree <- a[[1L]]; side <- a[[2L]]; fixtures <- a[[3L]]; out <- a[[4L]]; base_store <- a[[5L]]",
    "suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = FALSE))",
    "setwd(tree)",
    "for (e in parse(fixtures, keep.source = FALSE)) {",
    "    if (is.call(e) && identical(e[[1L]], as.name('<-')) && is.name(e[[2L]]) && startsWith(as.character(e[[2L]]), 'intEt')) eval(e, globalenv())",
    "}",
    "t1 <- intEtPlatformTerms()",
    "paths <- intEtProject(.local_envir = globalenv())",
    "manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0('store-k191-', side)))",
    "scan <- integrationScanProject(paths$integration_dir)",
    "integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))",
    "model <- runIntegrationModel(manager, list(), trainer = intEtTrainer())",
    "mapping <- intEtMapping(intEtCorpusSets(manager, model$value))",
    "request <- list(factor = 'Factor1', view = 'metabolomics_GCMS', database = 'kegg', id_column = 'kegg_id')",
    "camera <- runIntegrationEnrichment(manager, request, mapping = mapping)",
    "gsea <- runIntegrationEnrichment(manager, c(request, list(test = 'gsea')), mapping = mapping)",
    "saveRDS(list(t1 = t1, camera = as.data.frame(camera$run$value$terms), gsea = as.data.frame(gsea$run$value$terms),",
    "    camera_test = camera$entry$settings$test), file.path(out, paste0(side, '-terms.rds')))",
    "if (side == 'base') {",
    "    store <- dirname(manager$.context$store$directory)",
    "    closeWorkflowStateOwner(manager)",
    "    stopifnot(file.copy(store, base_store, recursive = TRUE))",
    "} else {",
    "    copied <- file.path(base_store, basename(list.files(base_store)[[1L]]))",
    "    identity <- readSavedStoreIdentity(file.path(copied, '.workflow-entities', 'store.json'))",
    "    old <- newWorkflowStateOwner(newWorkflowEntityStore(copied, identity$project_id, identity$workflow_id), reopen = TRUE)",
    "    recorded <- integrationRecordedEnrichment(old)",
    "    tests_recorded <- vapply(recorded, function(item) as.character(item$entry$settings$test), character(1))",
    "    replay <- integrationReplayEnrichment(old, recorded)",
    "    replay$test_recorded <- tests_recorded[replay$enrichment]",
    "    entries <- old$.context$snapshot$exports$integration_mofa_model$enrichment",
    "    replay$test_after <- vapply(replay$enrichment, function(key) as.character(entries[[key]]$settings$test), character(1))",
    "    utils::write.table(replay, file.path(out, 'replay.tsv'), sep = '\\t', quote = FALSE, row.names = FALSE)",
    "    item <- recorded[['Factor1|metabolomics_GCMS|kegg|camera']]",
    "    entry <- item$entry; s <- entry$settings",
    "    export <- old$.context$snapshot$exports$integration_mofa_model",
    "    declaration <- integrationEnrichmentCorrelatedOperation()",
    "    parameters <- list(ranking = entry$ranking$entity_id, database = s$database, id_column = s$id_column,",
    "        mapping = entry$stored$sha256, organism = s$kegg_organism, min_set_size = s$min_set_size,",
    "        max_set_size = s$max_set_size, test = 'camera', test_settings = INTEGRATION_ENRICHMENT_TEST_SETTINGS[['camera']],",
    "        values_slot = 'views', design_covariates = 'none')",
    "    at_head <- withWorkflowStateScope(old, function() workflowComponentScientificIdentity(old, export$source,",
    "        declaration$name, parameters[declaration$parameters]))",
    "    recorded_identity <- withWorkflowStateScope(old, function() describeSavedWorkflowEntity(old$.context$store,",
    "        entry$reference$entity_id)$identity)",
    "    keys <- union(names(at_head), names(recorded_identity))",
    "    identity <- data.frame(component = keys, equal = vapply(keys, function(k) identical(at_head[[k]], recorded_identity[[k]]),",
    "        logical(1)))",
    "    utils::write.table(identity, file.path(out, 'identity.tsv'), sep = '\\t', quote = FALSE, row.names = FALSE)",
    "    closeWorkflowStateOwner(old)",
    "    closeWorkflowStateOwner(manager)",
    "}"), child)

base_store <- file.path(work, "base-store")
dir.create(base_store)
run <- function(tree, side) {
    log <- file.path(out_dir, paste0(side, ".log"))
    code <- system2(rscript, c("--vanilla", child, tree, side, fixtures, work, base_store), stdout = log, stderr = log)
    cat(side, "exit", code, "\n")
    code
}
stopifnot(run(base_tree, "base") == 0L)
stopifnot(length(list.files(base_store)) == 1L)
stopifnot(run(head, "head") == 0L)

b <- readRDS(file.path(work, "base-terms.rds"))
h <- readRDS(file.path(work, "head-terms.rds"))
sha <- function(x) digest::digest(x, algo = "sha256")
compare <- data.frame(check = c("t1_terms_camera_and_gsea", "corpus_camera_terms", "corpus_gsea_terms", "corpus_default_test"),
    base_sha256 = c(sha(b$t1), sha(b$camera), sha(b$gsea), b$camera_test),
    head_sha256 = c(sha(h$t1), sha(h$camera), sha(h$gsea), h$camera_test),
    identical = c(identical(b$t1, h$t1), identical(b$camera, h$camera), identical(b$gsea, h$gsea),
        identical(b$camera_test, h$camera_test)))
utils::write.table(compare, file.path(out_dir, "base-head.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
replay <- utils::read.delim(file.path(work, "replay.tsv"), stringsAsFactors = FALSE)
identity <- utils::read.delim(file.path(work, "identity.tsv"), stringsAsFactors = FALSE)
file.copy(file.path(work, c("replay.tsv", "identity.tsv")), out_dir)
checks <- c(
    terms_identical = all(compare$identical),
    replay_two_entries = nrow(replay) == 2L,
    replay_terms_identical = all(replay$terms_identical),
    replay_tests_kept = all(replay$test_recorded == replay$test_after) && setequal(replay$test_recorded, c("camera", "gsea")),
    identity_equal_but_software = all(identity$equal[identity$component != "software"]) &&
        !isTRUE(identity$equal[identity$component == "software"]))
writeLines(sprintf("%s\t%s", names(checks), checks), file.path(out_dir, "verdict.txt"))
print(compare)
print(replay)
print(identity)
print(checks)
stopifnot(all(checks))
