# KIT-005 (user decision 2026-10-10): the automatic RUV percentage sweep prepares the canonical-correlation
# input (NIPALS imputation) once; every percentage's plot equals the single-call ruvCancor() result.

test_that("the percentage sweep prepares the canonical-correlation input once and reproduces ruvCancor per percentage", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    normalized <- suppressMessages(normaliseBetweenSamples(object, normalisation_method = "quantile"))
    # A few missing values so the NIPALS imputation branch is the one under test.
    samples <- normalized@design_matrix[[normalized@sample_id]]
    normalized@protein_quant_table[c(2L, 7L, 11L), samples[c(1L, 3L)]] <- NA_real_
    expect_true(anyNA(normalized@protein_quant_table))
    built <- function(plot) ggplot2::ggplot_build(plot)$data
    prepared <- suppressMessages(prepareRuvCancorInput(normalized, 5))
    expect_identical(dim(prepared), c(nrow(normalized@design_matrix), nrow(normalized@protein_quant_table)))
    expect_false(anyNA(prepared))
    for (percentage in c(20, 50)) {
        controls <- suppressMessages(getNegCtrlProtAnova(normalized, ruv_grouping_variable = "group",
            percentage_as_neg_ctrl = percentage, ruv_qval_cutoff = 0.05, ruv_fdr_method = "BH"))
        single <- suppressMessages(ruvCancor(normalized, ctrl = controls, num_components_to_impute = 5, ruv_grouping_variable = "group"))
        shared <- suppressMessages(ruvCancorFromInput(prepared, normalized, controls, "group"))
        expect_equal(built(shared), built(single), tolerance = 1e-12, info = percentage)
    }
    prepared_calls <- 0L
    original <- prepareRuvCancorInput
    sweep <- testthat::with_mocked_bindings(suppressMessages(findBestNegCtrlPercentage(normalized,
        percentage_range = c(20, 35, 50), num_components_to_impute = 5, ruv_grouping_variable = "group",
        ruv_fdr_method = "BH", separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, verbose = FALSE)),
        prepareRuvCancorInput = function(...) { prepared_calls <<- prepared_calls + 1L; original(...) })
    expect_identical(prepared_calls, 1L)
    expect_identical(nrow(sweep$optimization_results), 3L)
    expect_true(all(sweep$optimization_results$status == "ok"))
    direct <- suppressMessages(ruvCancor(normalized, ctrl = sweep$best_control_genes_index,
        num_components_to_impute = 5, ruv_grouping_variable = "group"))
    expect_equal(built(sweep$best_cancor_plot), built(direct), tolerance = 1e-12)
    expect_error(suppressMessages(findBestNegCtrlPercentage(normalized, percentage_range = 20,
        ruv_grouping_variable = "no_such_column", verbose = FALSE)), "not a column")
})
