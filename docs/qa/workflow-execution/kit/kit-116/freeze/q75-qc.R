# KIT-116 W2, quarantine row Q-75: the frozen values of three blocks of test-omics-qc-core.R, computed from the metab
# and trx QC functions KIT-104 W6 deleted (analyzeMetabQcItsdData, readMetabQcItsdState, validateMetabQcS4FinalizeState,
# trxQcFinalizeCollaborators), at the commit 827f2328, where the comparison passed. Run from a checkout of 827f2328:
#
#   Rscript q75-qc.R <checkout of 827f2328> <output dir>
#
# The script loads the package and test helpers of that checkout and the test file's own helper definitions, runs each
# block's old side on the block's inputs, and checks each old value against the core of the same commit as the block
# did. Values are literals, or sha256 of the value the block compared with expect_identical. Output:
# <output dir>/q75-frozen.R and q75-old-values.rds.

args <- commandArgs(trailingOnly = TRUE)
green <- normalizePath(args[[1]])
out_dir <- args[[2]]
stopifnot(identical(system2("git", c("-C", green, "rev-parse", "--short=8", "HEAD"), stdout = TRUE), "827f2328"))
pkgload::load_all(green, quiet = TRUE)
withr::local_dir(green)
# The helpers defer their temporary owners to testthat's teardown environment, which test_file initialises.
testthat:::local_teardown_env()
testthat::source_test_helpers(file.path(green, "tests", "testthat"), env = globalenv())
library(testthat)
stopifnot(workflowOwnerCapable())

exprs <- parse(file.path(green, "tests", "testthat", "test-omics-qc-core.R"), keep.source = FALSE)
for (e in exprs) if (!(is.call(e) && identical(e[[1L]], as.name("test_that")))) eval(e, globalenv())

digestOf <- function(x) digest::digest(x, algo = "sha256")
frozen <- list()
old_values <- list()
unfaithful <- character()
record <- function(group, key, value, core_value) {
    if (is.null(frozen[[group]])) frozen[[group]] <<- list()
    frozen[[group]][key] <<- list(value)
    if (!identical(value, core_value)) unfaithful <<- c(unfaithful, paste(group, key))
}
frozenCode <- function(frozen, prefix) {
    vapply(names(frozen), function(group) {
        values <- frozen[[group]]
        keys <- names(values)
        keys <- ifelse(make.names(keys) == keys, keys, paste0("`", keys, "`"))
        lines <- vapply(seq_along(values), function(i) {
            paste0("    ", keys[[i]], " = ", paste(deparse(values[[i]], width.cutoff = 500L), collapse = " "))
        }, character(1))
        code <- paste0(prefix, group, " <- list(\n", paste(lines, collapse = ",\n"), ")")
        stopifnot(identical(eval(parse(text = code)), values))
        code
    }, character(1))
}

# Metab's ITSD analysis (analyzeMetabQcItsdData on readMetabQcItsdState), per pattern ("NULL" for a NULL pattern).
(function() {
    itsd <- omicGuiConfig("metabolomics")$qc$itsd
    quiet <- function(...) invisible(NULL)
    for (pattern in list("", "^IS_", "^NO_MATCH$", NULL)) {
        key <- if (is.null(pattern)) "NULL" else pattern
        workflow <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
        manager <- workflow$state_manager
        today <- analyzeMetabQcItsdData(readMetabQcItsdState(manager, pattern), inputPattern = pattern, logInfoFn = quiet)
        core <- analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, pattern, itsd)
        old_values$itsd_metab[[paste0("pattern", key)]] <<- today
        record("itsd_metab", paste0(key, " metrics"), digestOf(today$metrics), digestOf(core$metrics))
        record("itsd_metab", paste0(key, " long"), digestOf(today$longData), digestOf(core$long_data))
        record("itsd_metab", paste0(key, " text"), today$resultText, core$result_text)
        record("itsd_metab", paste0(key, " total"), today$nIsTotal, core$total)
        record("itsd_metab", paste0(key, " pattern"), today$pattern, core$pattern)
        record("itsd_metab", paste0(key, " by_assay"), digestOf(today$metricsByAssay), digestOf(core$metrics_by_assay))
    }
    cases <- list(no_standard = list(object = module_ci_metab_qc_object(layout = "gc", include_itsd = FALSE),
        pattern = "^IS_"), no_pattern = list(object = module_ci_metab_qc_object(layout = "gc",
        internal_standard_regex = NA_character_), pattern = ""))
    for (name in names(cases)) {
        case <- cases[[name]]
        manager <- module_ci_metab_qc_workflow(case$object)$state_manager
        today <- qcCoreError(analyzeMetabQcItsdData(readMetabQcItsdState(manager, case$pattern),
            inputPattern = case$pattern, logInfoFn = quiet))
        stopifnot(!is.null(today))
        record("itsd_metab", paste0("refusal_", name), today,
            qcCoreError(analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, case$pattern, itsd)))
    }
    manager <- module_ci_lipid_qc_workflow()$state_manager
    today <- qcCoreError(analyzeMetabQcItsdData(readMetabQcItsdState(manager, ""), inputPattern = ""))
    stopifnot(!is.null(today))
    record("itsd_metab", "refusal_other_omic", today,
        qcCoreError(analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, "", itsd)))
})()

# Metab's finalise validation (validateMetabQcS4FinalizeState).
(function() {
    qc <- qcCoreConfig("metabolomics")
    aligned <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager, measures = "features")
    # A pass returned its input unchanged (the shell carries its owner's per-run entity reference, so the frozen value
    # is that identity, not a digest of the shell).
    passed <- validateMetabQcS4FinalizeState(aligned)
    record("finalise_metab", "aligned_returned_unchanged", identical(passed, aligned),
        identical(validateWorkflowOmicsQcFinalize(aligned, qc, metabolomicsDeclarationBinding), passed))
    zero_rows <- aligned
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    passed <- validateMetabQcS4FinalizeState(zero_rows)
    record("finalise_metab", "zero_rows_returned_unchanged", identical(passed, zero_rows),
        identical(validateWorkflowOmicsQcFinalize(zero_rows, qc, metabolomicsDeclarationBinding), passed))
    no_design <- aligned
    no_design@design_matrix <- no_design@design_matrix[0L, , drop = FALSE]
    no_assays <- aligned
    no_assays@metabolite_data <- list()
    attr(no_assays, "assay_rows") <- integer()
    other_class <- readLipidWorkflowShell(module_ci_lipid_qc_workflow()$state_manager)
    cases <- list(no_design = no_design, no_assays = no_assays, other_class = other_class)
    for (name in names(cases)) {
        today <- qcCoreError(validateMetabQcS4FinalizeState(cases[[name]]))
        stopifnot(!is.null(today))
        record("finalise_metab", name, today,
            qcCoreError(validateWorkflowOmicsQcFinalize(cases[[name]], qc, metabolomicsDeclarationBinding)))
    }
    # Parity (D11): a blank group passed metab's finalise (no core counterpart; the core refuses it).
    misaligned <- readMetabWorkflowShell(module_ci_metab_qc_workflow(qcCoreBlankGroup(module_ci_metab_qc_object()))$state_manager)
    record("finalise_metab", "blank_group_passed_class", class(validateMetabQcS4FinalizeState(misaligned))[[1L]],
        "MetaboliteAssayData")
})()

# Trx's finalise checks (trxQcFinalizeCollaborators()$validateFinalizeStateFn on the metadata).
(function() {
    fixture <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    shell <- fixture$object
    shell@count_data <- lapply(shell@count_data, function(table) table[0L, , drop = FALSE])
    attr(shell, "assay_rows") <- vapply(fixture$object@count_data, nrow, integer(1))
    qc <- qcCoreConfig("transcriptomics")
    today <- trxQcFinalizeCollaborators()$validateFinalizeStateFn
    metadata <- readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, fixture$owner)
    record("finalise_trx", "metadata_returned_unchanged", identical(today(metadata), metadata), TRUE)
    zero_rows <- shell
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    metadata_zero <- metadata
    metadata_zero$assay_rows[[1L]] <- 0L
    record("finalise_trx", "zero_rows", qcCoreError(today(metadata_zero)),
        qcCoreError(validateWorkflowOmicsQcFinalize(zero_rows, qc, transcriptomicsDeclarationBinding)))
})()

if (length(unfaithful)) stop("the old value is not the green core's for: ", paste(unfaithful, collapse = ", "))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(old_values, file.path(out_dir, "q75-old-values.rds"))
writeLines(frozenCode(frozen, "qcCoreFrozen_"), file.path(out_dir, "q75-frozen.R"))
cat("frozen groups:", paste(names(frozen), lengths(frozen), sep = "=", collapse = " "), "\n")
