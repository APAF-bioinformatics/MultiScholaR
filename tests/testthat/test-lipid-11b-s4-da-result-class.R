library(testthat)

test_that("LipidomicsDifferentialAbundanceResults keeps its slot contract and defaults", {
    lipid_object <- createLipidomicsAssayData(
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L2"),
                S1 = c(10, 20),
                S2 = c(11, 21),
                check.names = FALSE
            )
        ),
        design_matrix = data.frame(
            Run = c("S1", "S2"),
            group = c("A", "B"),
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "lipid_id",
        sample_id = "Run",
        group_id = "group"
    )

    result_object <- methods::new(
        "LipidomicsDifferentialAbundanceResults",
        theObject = lipid_object
    )

    expect_s4_class(result_object, "LipidomicsDifferentialAbundanceResults")
    expect_identical(methods::slotNames(result_object), c(
        "theObject",
        "fit.eb",
        "contrasts_results_table",
        "num_sig_diff_exp_bar_plot",
        "num_sig_diff_table",
        "volcano_plot",
        "interactive_volcano_plot",
        "p_value_dist_plot",
        "results_table_long",
        "results_table_wide"
    ))
    expect_s4_class(result_object@theObject, "LipidomicsAssayData")
    expect_null(result_object@fit.eb)
    expect_identical(result_object@contrasts_results_table, list())
    expect_identical(result_object@num_sig_diff_exp_bar_plot, list())
    expect_equal(nrow(result_object@num_sig_diff_table), 0)
    expect_identical(result_object@volcano_plot, list())
    expect_identical(result_object@interactive_volcano_plot, list())
    expect_identical(result_object@p_value_dist_plot, list())
    expect_equal(nrow(result_object@results_table_long), 0)
    expect_equal(nrow(result_object@results_table_wide), 0)
    expect_true(isTRUE(methods::validObject(result_object, test = TRUE)))
})

test_that("LipidomicsDifferentialAbundanceResults still requires a LipidomicsAssayData object", {
    expect_error(
        methods::new(
            "LipidomicsDifferentialAbundanceResults",
            theObject = list()
        ),
        regexp = "invalid object for slot \"theObject\""
    )
})
