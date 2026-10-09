# RUV-III replicate structure (user decision 2026-10-09): technical replicates when the design has
# them, otherwise the RUV grouping variable.

test_that("technical replicates become the RUV-III replicate sets when the design has them", {
    design <- data.frame(Run = paste0("S", 1:6), group = rep(c("A", "B"), each = 3),
        replicates = c(1L, 1L, 2L, 1L, 1L, 2L), stringsAsFactors = FALSE)
    blocks <- resolveProtRuvReplicateBlocks(design, "Run", "group", "replicates", "group")
    expect_identical(blocks$structure, "technical_replicates")
    expect_identical(blocks$blocks, 4L)
    expect_identical(blocks$design$.ruv_replicate_block, c("A_1", "A_1", "A_2", "B_1", "B_1", "B_2"))
    M <- getRuvIIIReplicateMatrixHelper(blocks$design, Run, .ruv_replicate_block)
    expect_equal(dim(M), c(6L, 4L))
    expect_equal(unname(M["S1", "A_1"]), 1)
    expect_equal(unname(M["S3", "A_1"]), 0)
    expect_equal(unname(rowSums(M)), rep(1, 6))
})

test_that("without technical replicates every sample of the grouping variable is a replicate", {
    design <- data.frame(Run = paste0("S", 1:4), group = c("A", "A", "B", "B"),
        replicates = 1:4, batch = c("x", "y", "x", "y"), stringsAsFactors = FALSE)
    blocks <- resolveProtRuvReplicateBlocks(design, "Run", "group", "replicates", "group")
    expect_identical(blocks$structure, "group")
    expect_identical(blocks$blocks, 2L)
    expect_identical(blocks$design$.ruv_replicate_block, c("A", "A", "B", "B"))
    by_batch <- resolveProtRuvReplicateBlocks(design, "Run", "group", "replicates", "batch")
    expect_identical(by_batch$design$.ruv_replicate_block, c("x", "y", "x", "y"))
    no_column <- resolveProtRuvReplicateBlocks(design[c("Run", "group")], "Run", "group", "replicates", "group")
    expect_identical(no_column$structure, "group")
    expect_error(resolveProtRuvReplicateBlocks(design[c("Run", "group")], "Run", "group", "replicates", "batch"),
        "absent from the design matrix")
})

test_that("ruvIII_C_Varying records the replicate structure it used", {
    object <- module_ci_prot_norm_object()
    design <- object@design_matrix
    expected <- resolveProtRuvReplicateBlocks(design, object@sample_id, object@group_id,
        object@technical_replicate_id, "group")$structure
    ctrl <- stats::setNames(rep(TRUE, nrow(object@protein_quant_table)), object@protein_quant_table[[object@protein_id_column]])
    corrected <- ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = 1L, ctrl = ctrl)
    expect_identical(corrected@args$ruvIII_C_Varying$ruv_replicate_structure, expected)
    expect_identical(corrected@args$ruvIII_C_Varying$ruv_grouping_variable, "group")
})
