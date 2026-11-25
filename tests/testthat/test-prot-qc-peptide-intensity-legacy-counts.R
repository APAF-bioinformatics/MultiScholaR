# KIT-187 (VAL-71, D6): the peptide intensity filter's deprecated percentage-only call reads each percentage as the
# three-decimal rounding of an exact share, as the protein filter does; the count call (the GUI's) is unchanged.

k187PeptideCall <- function(...) {
    design <- data.frame(Run = paste0("R", 1:6), group = rep(c("G1", "G2"), each = 3L), stringsAsFactors = FALSE)
    table <- data.frame(Protein.Ids = c("P1", "P2"), Stripped.Sequence = c("FULL", "TWOOFTHREE"),
        R1 = c(20, 10), R2 = c(20, 20), R3 = c(20, 20), R4 = c(20, 10), R5 = c(20, 20), R6 = c(20, 20),
        stringsAsFactors = FALSE)
    peptideIntensityFilteringHelper(table, design, min_peptide_intensity_threshold = 15, sample_id_column = "Run",
        grouping_variable = "group", protein_id_column = "Protein.Ids", peptide_sequence_column = "Stripped.Sequence",
        peptide_quantity_column = "Peptide.Normalised", core_utilisation = NA, ...)
}

test_that("T10: a peptide observed in two of three runs of each group passes the percentage call at 33.333 and 0", {
    expect_warning(result <- k187PeptideCall(groupwise_percentage_cutoff = 33.333, max_groups_percentage_cutoff = 0),
        "deprecated")
    expect_setequal(result$Stripped.Sequence, c("FULL", "TWOOFTHREE"))
})

test_that("T10: the count call (the GUI's path) keeps the same peptide", {
    result <- k187PeptideCall(min_reps_per_group = 2, min_groups = 2)
    expect_setequal(result$Stripped.Sequence, c("FULL", "TWOOFTHREE"))
})
