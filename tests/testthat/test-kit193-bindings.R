# KIT-193 section 8.2: the binding tripwire. Each filled row of tests/testdata/kit193/bindings.tsv names an object the
# Lean models were checked against, with its digest (a package function: the SHA-256 of deparse(body(fn)); a constant:
# of deparse(x)) or its version (a library) as recorded when its paired test first ran on merged code. This test
# recomputes each one and fails when the object has changed: the model must then be re-reviewed (who may change a row:
# ticket 8.2). Source references are removed before deparsing, so a package loaded with them (load_all) digests as the
# installed package does. Rows still `pending` are those whose pair waits for a product ticket (listed below).

kit193BoundNames <- c(
    "exceedsRoundedPercent", "roundedPercentNumerator", "allowedByRoundedPercent", "updateMissingValueParameters",
    "removeRowsWithMissingValuesPercentHelper", "workflowReplayRuleMarkers", "workflowReplayRule",
    "workflowFixedInForce", "workflowReplayRuleOptions", "buildWorkflowAnalysisRecord", "goAncestorPairs",
    "propagateGoAnnotations", "buildUniprotGoGeneSets", "goPropagationRelations", "goAspectRoots", "processEnrichments",
    "qvalue", "stats", "integrationEnrichmentPermutation", "INTEGRATION_ENRICHMENT_PERMUTATION_TIE", "ruv", "RUVIIIC",
    "restoreRuvFeatureLevels", "computeRuvCovariates", "computeOmicsRuvAssayFactors", "computeProtRuvModelFactors",
    "ruvExactContrastsFit", "limma", "edgeR", "integrationEnrichmentSetT", "integrationEnrichmentFloorMinimum",
    "ruvIII_C_Varying,ProteinQuantitativeData", "protRuvFirstPassControls", "omicsRuvFirstPassControls",
    "trxRuvDesignControls", "OMICS_RUV_FIRST_PASS_DROPPED_SHARE", "TRX_FIRST_PASS_DROPPED_SHARE")

#' Rows not yet recorded: target 2's (KIT-190)
kit193PendingNames <- c("workflowReplayRuleMarkers", "workflowReplayRule", "workflowFixedInForce",
    "workflowReplayRuleOptions", "buildWorkflowAnalysisRecord")

kit193Bindings <- function() {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", "bindings.tsv"), sep = "\t", quote = "",
        colClasses = "character", na.strings = character())
}

#' The value of one row on the loaded package
kit193BindingValue <- function(kind, name) {
    ns <- asNamespace("MultiScholaR")
    if (identical(kind, "library")) return(format(utils::packageVersion(name)))
    if (identical(kind, "package constant")) return(digest::digest(deparse(get(name, envir = ns)), algo = "sha256"))
    fn <- if (grepl(",", name, fixed = TRUE)) {
        parts <- strsplit(name, ",", fixed = TRUE)[[1L]]
        methods::getMethod(parts[[1L]], parts[[2L]], where = ns)@.Data
    } else {
        get(name, envir = ns)
    }
    digest::digest(deparse(body(utils::removeSource(fn))), algo = "sha256")
}

test_that("KIT-193 bindings: the table names exactly the bound objects; only the listed rows are pending", {
    table <- kit193Bindings()
    expect_setequal(unique(table$name), kit193BoundNames)
    expect_setequal(table$name[table$digest_or_version == "pending"], kit193PendingNames)
})

test_that("KIT-193 bindings: every filled row's object is unchanged since its Lean model was checked against it", {
    table <- kit193Bindings()
    table <- table[table$digest_or_version != "pending", , drop = FALSE]
    for (i in seq_len(nrow(table))) {
        row <- table[i, ]
        now <- kit193BindingValue(row$kind, row$name)
        testthat::expect(identical(now, row$digest_or_version), sprintf(paste("KIT-193 binding %s: %s has changed since",
            "the Lean model %s was checked against it; re-review the model (KIT-193 section 8.2)"), row$item, row$name,
            row$spec_name))
    }
})

test_that("KIT-193 bindings, planted defect: a recorded digest that differs from the object fails", {
    row <- kit193Bindings()
    row <- row[row$name == "computeRuvCovariates", , drop = FALSE]
    now <- kit193BindingValue(row$kind, row$name)
    testthat::expect_failure(testthat::expect(identical(now, strrep("0", 64L)), "changed"))
})
