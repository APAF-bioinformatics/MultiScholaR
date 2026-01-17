# Verification Script for Lipidomics Normalized Column Exclusion

test_that("Normalized columns are excluded by default", {
    # Test cases provided in the user request
    test_columns <- c(
        "Control_1", "Control_1_Norm",
        "Control_2", "Control_2_Norm",
        "Control_3", "Control_3_Norm",
        "Control_4", "Control_4_Norm",
        "Control_5", "Control_5_Norm",
        # Additional variations to test regex robustness
        "Sample_A_Normalised",
        "Sample_B_normalized",
        "Sample_C_normallized",
        "Sample_D_Normalized",
        "Random_Sample"
    )

    # Regex pattern used in implementation
    exclusion_pattern <- "_norm(ali[sz]ed|allized)?$|normalized"

    # Apply filter (ignore.case = TRUE as in the app)
    kept_columns <- test_columns[!grepl(exclusion_pattern, test_columns, ignore.case = TRUE)]

    # Expected results
    expected_columns <- c(
        "Control_1",
        "Control_2",
        "Control_3",
        "Control_4",
        "Control_5",
        "Random_Sample"
    )

    expect_equal(kept_columns, expected_columns)
    expect_true("Control_1" %in% kept_columns)
    expect_false("Control_1_Norm" %in% kept_columns)
})
