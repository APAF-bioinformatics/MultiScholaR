# The per-group reference loop hands each protein group the complete design; a group lacking runs
# where it was not observed must still filter exactly as the whole table does (KIT-001, 2026-10-09).

test_that("a protein-group batch with absent design runs filters exactly as the whole table", {
    design <- data.frame(Run = paste0("R", 1:4), group = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
    table <- data.frame(
        Protein.Group = c(rep("P1", 4), rep("P2", 3)),
        Stripped.Sequence = c(rep("AAA", 4), rep("BBB", 3)),
        Run = c("R1", "R2", "R3", "R4", "R1", "R2", "R3"),
        Peptide.Normalised = c(10, 12, 11, 13, 20, 21, 22), stringsAsFactors = FALSE)
    call <- function(part, complete) peptideIntensityFilteringHelper(part, design,
        min_peptide_intensity_threshold = 5, sample_id_column = "Run", grouping_variable = "group",
        protein_id_column = "Protein.Group", peptide_sequence_column = "Stripped.Sequence",
        peptide_quantity_column = "Peptide.Normalised", core_utilisation = NA,
        min_reps_per_group = 2, min_groups = 2, strict_mode = FALSE, return_filter_result = TRUE,
        require_complete_design_runs = complete)
    whole <- call(table, TRUE)
    expect_error(call(table[table$Protein.Group == "P2", ], TRUE), "Design run\\(s\\) absent from input: R4")
    parts <- lapply(split(table, table$Protein.Group), call, complete = FALSE)
    per_group <- do.call(rbind, lapply(parts, `[[`, "data"))
    rownames(per_group) <- NULL
    expect_identical(per_group, whole$data)
    support <- do.call(rbind, lapply(parts, `[[`, "support_table"))
    rownames(support) <- NULL
    expect_identical(support, whole$support_table)
    expect_identical(unique(support$design_run_count[support$group == "B"]), 2L)
    expect_error(call(transform(table, Run = sub("R4", "R9", Run)), FALSE), "Unmapped input run")
})
