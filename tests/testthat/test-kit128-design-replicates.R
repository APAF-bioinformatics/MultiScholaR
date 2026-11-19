# KIT-128 audit K128-2: the shared design builder (lipid, metab, trx) gives each run `tech_rep_group = <group>_<replicates>`,
# and declared technical replicates are now blocked (lipid, metab) or summed (trx). Runs of a group assigned in two
# batches, each starting at replicate 1, would share numbers and so become technical replicates by accident. The save
# refuses a replicate number shared inside a group unless every run sharing it was marked by the technical-replicate
# action (its `tech_reps`).
library(testthat)

kit128BuilderDesign <- function() {
  tibble::tibble(Run = paste0("S", 1:6), group = NA_character_, factor1 = NA_character_, factor2 = NA_character_,
    factor3 = NA_character_, batch = NA_character_, replicates = NA_integer_, tech_reps = NA_integer_)
}

kit128BuilderAssays <- function() {
  values <- as.data.frame(matrix(1:12, nrow = 2L, dimnames = list(NULL, paste0("S", 1:6))))
  list(LCMS_Pos = cbind(data.frame(LipidName = c("L1", "L2"), LipidClass = "PC"), values))
}

kit128BuilderSave <- function(design) {
  buildOmicsDesignSavePayload(design, character(), kit128BuilderAssays(), list(lipid_id_col = "LipidName",
    annotation_col = "LipidClass"), NULL, list(), "~ 0 + group", "lipid_id_col")
}

test_that("KIT-128 K128-2: runs of one group assigned in two batches from replicate 1 are refused with the runs named", {
  design <- applyOmicsDesignMetadataAssignment(kit128BuilderDesign(), c("S1", "S2"), "A", "", "", 1L)
  design <- applyOmicsDesignMetadataAssignment(design, c("S3", "S4"), "A", "", "", 1L)
  design <- applyOmicsDesignMetadataAssignment(design, c("S5", "S6"), "B", "", "", 1L)
  refused <- kit128BuilderSave(design)
  expect_type(refused, "character")
  expect_match(refused, "S1, S2, S3, S4", fixed = TRUE)
  expect_match(refused, "Assign Technical Replicates", fixed = TRUE)
  # Numbered on from the first batch: accepted, every run its own biological replicate.
  design <- applyOmicsDesignMetadataAssignment(design, c("S3", "S4"), "A", "", "", 3L)
  payload <- kit128BuilderSave(design)
  expect_type(payload, "list")
  expect_false(anyDuplicated(payload$design_matrix$tech_rep_group) > 0L)
})

test_that("KIT-128 K128-2: runs marked by the technical-replicate action share their number and are accepted", {
  design <- applyOmicsDesignMetadataAssignment(kit128BuilderDesign(), c("S1", "S2", "S3", "S4"), "A", "", "", 1L)
  design <- applyOmicsDesignMetadataAssignment(design, c("S5", "S6"), "B", "", "", 1L)
  design <- applyOmicsDesignTechRepAssignment(design, c("S1", "S2"), "lowest")$designMatrix
  payload <- kit128BuilderSave(design)
  expect_type(payload, "list")
  expect_identical(payload$design_matrix$tech_rep_group[1:2], c("A_1", "A_1"))
})
