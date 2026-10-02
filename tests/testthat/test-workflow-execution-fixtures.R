# WFX-001 audit repairs: current-package tests, with no external evidence library.
# Direct S4 slots follow the existing fixture convention for these legacy classes.

invokeWorkflowFilter <- function(x = makeWorkflowExecutionFilterFixture(), threshold = x$threshold, inclusion = x$inclusion) {
    MultiScholaR::filterMinNumPeptidesPerSampleHelper(x$input,
        peptides_per_sample_cutoff = threshold, sample_id_column = "Run",
        protein_id_column = "Protein.Group", peptide_sequence_column = "Stripped.Sequence",
        core_utilisation = NA, inclusion_list = inclusion, design_matrix = x$design,
        design_sample_id_column = "Run", return_filter_result = TRUE)
}

testthat::test_that("mixed retention, native types and independent expected ledgers", {
    fixture <- makeWorkflowExecutionFilterFixture()
    for (n in c(1L, 1000L)) {
        x <- makeWorkflowExecutionFilterFixture(n)
        a <- invokeWorkflowFilter(x)
        testthat::expect_identical(a$data$source_ordinal, x$expected$retained_source_ordinal)
        testthat::expect_identical(unique(a$data$Run), x$expected$retained)
        testthat::expect_type(a$data$Run, "character")
        support <- a$support_table
        testthat::expect_identical(support$distinct_peptide_identity_count,
            rep(c(2L, 1L, 1L, 4L), each = n))
        testthat::expect_identical(support$decision_reason,
            rep(c("meets_distinct_identity_threshold", "explicit_inclusion_override",
                "below_distinct_identity_threshold", "meets_distinct_identity_threshold"), each = n))
        testthat::expect_identical(a$removal_ledger$Run, x$expected$removed)
        testthat::expect_identical(a$removal_ledger$distinct_peptide_identity_count, rep(1L, n))
        testthat::expect_identical(a$removal_ledger$required_distinct_peptide_identities, rep(2L, n))
        testthat::expect_identical(a$summary$inclusion_entries_absent_from_data, "absent")
        testthat::expect_equal(a$summary$retained_sample_count, 3L * n)
        testthat::expect_equal(a$summary$removed_sample_count, n)
    }
    a <- invokeWorkflowFilter()
    for (chunk in c(1L, 2L, 3L, 7L)) {
        x <- fixture
        x$input <- do.call(rbind, split(x$input, ceiling(seq_len(nrow(x$input)) / chunk)))
        rownames(x$input) <- NULL
        testthat::expect_identical(invokeWorkflowFilter(x), a)
    }
    testthat::expect_false(identical(invokeWorkflowFilter(inclusion = character())$data$source_ordinal,
        fixture$expected$retained_source_ordinal))
    mutant <- fixture
    mutant$input <- subset(mutant$input, !is.na(Quantity) & Quantity > 0)
    testthat::expect_false(identical(invokeWorkflowFilter(mutant)$data$source_ordinal,
        fixture$expected$retained_source_ordinal))
})
testthat::test_that("empty/all/none and wide annotations", {
    fixture <- makeWorkflowExecutionFilterFixture()
    testthat::expect_equal(nrow(invokeWorkflowFilter(threshold = 1L)$data), 10L)
    testthat::expect_equal(nrow(invokeWorkflowFilter(threshold = 5L, inclusion = character())$data), 0L)
    x <- fixture
    x$input <- x$input[0, , drop = FALSE]
    x$design <- x$design[0, , drop = FALSE]
    testthat::expect_equal(nrow(invokeWorkflowFilter(x)$data), 0L)
    wide <- fixture
    wide$input$Protein.Ids <- strrep("A", 2048L)
    testthat::expect_identical(invokeWorkflowFilter(wide)$support_table, invokeWorkflowFilter()$support_table)
})
testthat::test_that("invalid inputs preserve native rejection behavior", {
    fixture <- makeWorkflowExecutionFilterFixture()
    for (k in list(0, -1, 1.5, NA_real_, Inf))
        testthat::expect_error(invokeWorkflowFilter(threshold = k), "positive integer")
    for (key in c("Run", "Protein.Group", "Stripped.Sequence")) {
        for (value in c(NA_character_, " ")) {
            x <- fixture
            x$input[[key]] <- as.character(x$input[[key]])
            x$input[[key]][1] <- value
            testthat::expect_error(invokeWorkflowFilter(x), "identity is missing")
        }
    }
    for (kind in c("duplicate", "missing", "blank", "unmapped")) {
        x <- fixture
        if (kind == "duplicate") x$design <- rbind(x$design, x$design[1, ])
        if (kind == "missing") x$design$Run[1] <- NA_character_
        if (kind == "blank") x$design$Run[1] <- " "
        if (kind == "unmapped") x$design <- x$design[-1, ]
        testthat::expect_error(invokeWorkflowFilter(x), "duplicate run|identity is missing|Unmapped input")
    }
})
testthat::test_that("native S4 design cleanup and argument ledgers", {
    fixture <- makeWorkflowExecutionFilterFixture()
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@peptide_data <- fixture$input
    object@sample_id <- "Run"
    object@protein_id_column <- "Protein.Group"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@design_matrix <- rbind(fixture$design[4:1, ], data.frame(Run = "extra", group = "case"))
    result <- MultiScholaR::filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = 2L, core_utilisation = NA, inclusion_list = fixture$inclusion)
    testthat::expect_identical(result@design_matrix$Run, fixture$expected$retained)
    testthat::expect_identical(result@args$filterMinNumPeptidesPerSample$support_table,
        invokeWorkflowFilter()$support_table)
    testthat::expect_identical(result@peptide_data, invokeWorkflowFilter()$data)
})

testthat::test_that("both identity components and abundance-independent counts are necessary", {
    fixture <- makeWorkflowExecutionFilterFixture()
    correct <- invokeWorkflowFilter(fixture)
    # Mutate counting only. Retain the original table for the final sample join.
    mutantRetained <- function(keys, observed_only = FALSE, inclusion = fixture$inclusion) {
        rows <- fixture$input
        if (observed_only) rows <- rows[!is.na(rows$Quantity), , drop = FALSE]
        counts <- table(unique(rows[c("Run", keys)])$Run)
        included <- unique(trimws(inclusion[!is.na(inclusion)]))
        samples <- names(counts)[counts >= fixture$threshold | names(counts) %in% included]
        fixture$input$source_ordinal[as.character(fixture$input$Run) %in% samples]
    }
    testthat::expect_identical(correct$data$source_ordinal, fixture$expected$retained_source_ordinal)
    for (keys in list("Protein.Group", "Stripped.Sequence")) {
        testthat::expect_false(identical(mutantRetained(keys), correct$data$source_ordinal))
    }
    testthat::expect_false(identical(mutantRetained(c("Protein.Group", "Stripped.Sequence"),
        observed_only = TRUE), correct$data$source_ordinal))
    testthat::expect_false(identical(mutantRetained(c("Protein.Group", "Stripped.Sequence"),
        inclusion = character()), correct$data$source_ordinal))
})
