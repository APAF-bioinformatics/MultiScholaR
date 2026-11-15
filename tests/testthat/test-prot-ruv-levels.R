# KIT-123 R2: an RUV-III correction keeps each protein's level. RUV-III estimates W from the controls' uncentred
# values, so their mean level entered W and every protein was shifted by a constant (parked row 14).

test_that("restoreRuvFeatureLevels removes a constant shift per feature and nothing else", {
    original <- matrix(c(20, 21, 22, 23, 30, NA, 32, 33), nrow = 2L, byrow = TRUE, dimnames = list(c("a", "b"), paste0("s", 1:4)))
    corrected <- original + c(5, -7) + matrix(c(0.1, -0.1, 0.2, -0.2, 0.3, NA, -0.1, -0.2), nrow = 2L, byrow = TRUE)
    corrected["a", "s4"] <- NA # not corrected in this sample: left out of both means
    restored <- restoreRuvFeatureLevels(corrected, original)
    expect_equal(mean(restored["a", 1:3]), mean(original["a", 1:3]))
    expect_equal(mean(restored["b", c(1, 3, 4)]), mean(original["b", c(1, 3, 4)]))
    expect_identical(is.na(restored), is.na(corrected))
    # differences between samples within a feature are the correction's own
    expect_equal(restored["a", 1:3] - restored["a", 1L], corrected["a", 1:3] - corrected["a", 1L])
    expect_equal(restored["b", c(3, 4)] - restored["b", 1L], corrected["b", c(3, 4)] - corrected["b", 1L])
    # a feature with no value in common is left as it is
    none <- restoreRuvFeatureLevels(matrix(c(NA, 1), 1L, dimnames = list("a", c("s1", "s2"))),
        matrix(c(1, NA), 1L, dimnames = list("a", c("s1", "s2"))))
    expect_identical(unname(none[1L, ]), c(NA, 1))
    expect_error(restoreRuvFeatureLevels(corrected, original[, 4:1]))
})

test_that("ruvIII_C_Varying keeps every protein's mean and still corrects within proteins", {
    object <- module_ci_prot_norm_object()
    samples <- as.character(object@design_matrix[[object@sample_id]])
    ids <- object@protein_quant_table[[object@protein_id_column]]
    before <- as.matrix(object@protein_quant_table[, samples]); rownames(before) <- ids
    ctrl <- stats::setNames(rep(TRUE, length(ids)), ids)
    corrected <- ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = 1L, ctrl = ctrl)
    # proteins and samples the correction gave no value to are dropped by the method
    samples <- intersect(samples, names(corrected@protein_quant_table))
    after <- as.matrix(corrected@protein_quant_table[, samples]); rownames(after) <- corrected@protein_quant_table[[corrected@protein_id_column]]
    kept <- intersect(rownames(before), rownames(after))
    expect_gt(length(kept), 0L)
    before_all <- before
    before <- before[, samples, drop = FALSE]
    both <- is.finite(before[kept, , drop = FALSE]) & is.finite(after[kept, , drop = FALSE])
    mean_in <- function(m) rowSums(replace(m[kept, , drop = FALSE], !both, 0)) / rowSums(both)
    expect_equal(mean_in(after), mean_in(before), tolerance = 1e-9)
    # the correction itself is still applied: values move within proteins
    centred <- function(m) (m[kept, , drop = FALSE] - mean_in(m))[both]
    expect_gt(max(abs(centred(after) - centred(before))), 1e-6)
    # and it is RUV-III-C's own output up to a constant per protein
    Y <- t(before_all)
    blocks <- resolveProtRuvReplicateBlocks(object@design_matrix, object@sample_id, object@group_id, object@technical_replicate_id, "group")
    M <- getRuvIIIReplicateMatrixHelper(blocks$design, !!rlang::sym(object@sample_id), .ruv_replicate_block)
    raw <- t(RUVIIIC::RUVIII_C_Varying(k = 1L, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = names(ctrl), progress = FALSE))
    difference <- after[kept, samples, drop = FALSE] - raw[kept, samples, drop = FALSE]
    spread <- apply(difference, 1L, function(row) if (sum(is.finite(row)) > 1L) diff(range(row, na.rm = TRUE)) else 0)
    expect_lt(max(spread), 1e-8)
    expect_gt(max(abs(rowMeans(raw[kept, samples, drop = FALSE], na.rm = TRUE) - rowMeans(before[kept, samples, drop = FALSE], na.rm = TRUE))), 0.01)
})

# KIT-124 W1: the reference match on complete data and on data with missing values, k of 1 and 2 (the test above is one
# four-sample fixture at k 1). The reference's replicate matrix is built here (`ruv::replicate.matrix`), not by the
# package's helper.
protRuvReferenceInput <- function(missing) {
    set.seed(20261116)
    samples <- sprintf("S%02d", 1:8); ids <- sprintf("P%02d", 1:30)
    group <- rep(c("A", "B"), each = 4L)
    values <- matrix(stats::rnorm(30L, 22, 2), 30L, 8L, dimnames = list(ids, samples)) + outer(stats::rnorm(30L), stats::rnorm(8L)) +
        matrix(stats::rnorm(240L, sd = 0.2), 30L, 8L)
    values[1:6, group == "B"] <- values[1:6, group == "B"] + 1.5
    # seven missing values, three of them in control proteins (P14, P20, P27)
    if (missing) values[cbind(c(2L, 5L, 9L, 9L, 14L, 20L, 27L), c(1L, 6L, 2L, 7L, 4L, 8L, 3L))] <- NA
    table <- data.frame(Protein.Ids = ids, values, check.names = FALSE, row.names = NULL)
    object <- ProteinQuantitativeData(protein_quant_table = table, protein_id_column = "Protein.Ids", protein_id_table = table["Protein.Ids"],
        design_matrix = data.frame(Run = samples, group = group, replicates = rep(1:4, 2L)), sample_id = "Run", group_id = "group",
        technical_replicate_id = "replicates", args = list())
    list(object = object, values = values, group = group, ctrl = stats::setNames(seq_along(ids) > 12L, ids))
}

test_that("(KIT-124 W1) ruvIII_C_Varying is RUVIIIC::RUVIII_C_Varying up to a constant per protein, with each protein's mean kept", {
    for (missing in c(FALSE, TRUE)) for (k in 1:2) {
        case <- sprintf("missing values %s, k %d", missing, k)
        input <- protRuvReferenceInput(missing)
        corrected <- suppressMessages(ruvIII_C_Varying(input$object, ruv_grouping_variable = "group", ruv_number_k = k, ctrl = input$ctrl))
        after <- as.matrix(corrected@protein_quant_table[, colnames(input$values)])
        rownames(after) <- corrected@protein_quant_table$Protein.Ids
        expect_identical(dimnames(after), dimnames(input$values), info = case)
        Y <- t(input$values)
        M <- ruv::replicate.matrix(input$group)
        rownames(M) <- rownames(Y)
        raw <- t(RUVIIIC::RUVIII_C_Varying(k = k, Y = Y, M = M, toCorrect = colnames(Y), potentialControls = names(input$ctrl)[input$ctrl],
            progress = FALSE))
        expect_identical(is.na(after), is.na(raw), info = case)
        expect_identical(sum(is.na(after)), if (missing) 7L else 0L, info = case)
        # equal up to a constant per protein
        spread <- apply(after - raw, 1L, function(row) diff(range(row, na.rm = TRUE)))
        expect_lt(max(spread), 1e-8, label = paste("the largest spread of package minus reference within a protein,", case))
        # the constant is the level: the package's mean is the mean before, the reference's is not
        both <- is.finite(after) & is.finite(input$values)
        mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
        expect_equal(mean_in(after), mean_in(input$values), tolerance = 1e-9, info = case)
        expect_gt(max(abs(mean_in(raw) - mean_in(input$values))), 1, label = paste("the reference's level shift,", case))
    }
})
