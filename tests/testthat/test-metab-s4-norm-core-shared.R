# fidelity-coverage-compare: shared
library(methods)
library(testthat)

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      annotation_id_column = "character",
      internal_standard_regex = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      args = "list",
      technical_replicate_id = "character",
      group_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      annotation_id_column = "annotation_id",
      internal_standard_regex = "ISTD",
      design_matrix = data.frame(),
      sample_id = "Run",
      args = list(),
      technical_replicate_id = "TechRep",
      group_id = "group"
    )
  )
}

if (!methods::isGeneric("logTransformAssays")) {
  methods::setGeneric(
    "logTransformAssays",
    function(theObject, offset = 1, ...) {
      standardGeneric("logTransformAssays")
    }
  )
}

if (!methods::isGeneric("normaliseUntransformedData")) {
  methods::setGeneric(
    "normaliseUntransformedData",
    function(theObject,
             method = "ITSD",
             itsd_feature_ids = NULL,
             itsd_pattern_columns = NULL,
             itsd_aggregation = "sum",
             remove_itsd_after_norm = TRUE,
             ...) {
      standardGeneric("normaliseUntransformedData")
    }
  )
}

if (!methods::isGeneric("normaliseBetweenSamples")) {
  methods::setGeneric(
    "normaliseBetweenSamples",
    function(theObject, normalisation_method = NULL) {
      standardGeneric("normaliseBetweenSamples")
    }
  )
}

if (!methods::isGeneric("cleanDesignMatrix")) {
  methods::setGeneric(
    "cleanDesignMatrix",
    function(theObject) {
      standardGeneric("cleanDesignMatrix")
    }
  )
}

newMetabNormCoreObject <- function(assay_list, design_matrix, args = list()) {
  sample_ids <- unique(unlist(lapply(assay_list, function(tbl) grep("^Sample_", colnames(tbl), value = TRUE))))
  valid_design_matrix <- design_matrix[design_matrix$Run %in% sample_ids, , drop = FALSE]

  if (nrow(valid_design_matrix) == 0) {
    valid_design_matrix <- data.frame(
      Run = sample_ids,
      group = rep("group", length(sample_ids)),
      TechRep = rep("rep", length(sample_ids)),
      stringsAsFactors = FALSE
    )
  }

  object <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = assay_list,
    metabolite_id_column = "Name",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "ISTD",
    design_matrix = valid_design_matrix,
    sample_id = "Run",
    technical_replicate_id = "TechRep",
    group_id = "group",
    args = args
  )

  object@design_matrix <- design_matrix
  object
}

newMetabNormCoreEmptyObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "Name",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "ISTD",
    design_matrix = data.frame(
      Run = "Sample_1",
      group = "A",
      TechRep = "rep1",
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    technical_replicate_id = "TechRep",
    group_id = "group",
    args = list()
  )
}

test_that("logTransformAssays preserves negative-value cleanup, unnamed assays, and validation branches", {
  the_object <- newMetabNormCoreObject(
    assay_list = list(
      tibble::tibble(
        Name = c("M1", "M2", "ISTD1"),
        annotation_id = c("alpha", "beta", "ISTD spike"),
        Sample_1 = c(-1, 10, 5),
        Sample_2 = c(4, 20, 10),
        Sample_3 = c(8, 40, 20),
        Sample_4 = c(16, 80, 40)
      )
    ),
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      group = c("A", "A", "B", "B"),
      TechRep = c("TR1", "TR1", "TR2", "TR2"),
      stringsAsFactors = FALSE
    )
  )
  names(the_object@metabolite_data) <- ""

  transformed <- suppressWarnings(logTransformAssays(the_object, offset = 1))
  transformed_assay <- transformed@metabolite_data$Assay_1

  expect_named(transformed@metabolite_data, "Assay_1")
  expect_equal(transformed_assay$Sample_1[1], 0)
  expect_equal(transformed_assay$Sample_4[2], log2(81))
  expect_true(isTRUE(transformed@args$log_transformed))
  expect_identical(transformed@args$log_transform_offset, 1)

  expect_warning(
    logTransformAssays(newMetabNormCoreEmptyObject(), offset = 1),
    "no assays",
    ignore.case = TRUE
  )

  expect_error(
    logTransformAssays(the_object, offset = 0),
    "single positive numeric"
  )
})

test_that("normaliseUntransformedData preserves regex and manual ITSD normalization flows", {
  assay_tbl <- tibble::tibble(
    Name = c("M1", "M2", "ISTD1"),
    annotation_id = c("alpha", "beta", "ISTD spike"),
    Sample_1 = c(-1, 10, 5),
    Sample_2 = c(4, 20, 10),
    Sample_3 = c(8, 40, 20),
    Sample_4 = c(16, 80, 40)
  )

  the_object <- newMetabNormCoreObject(
    assay_list = list(LCMS_Pos = assay_tbl),
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      group = c("A", "A", "B", "B"),
      TechRep = c("TR1", "TR1", "TR2", "TR2"),
      stringsAsFactors = FALSE
    )
  )

  regex_norm <- suppressWarnings(normaliseUntransformedData(the_object, method = "ITSD"))
  regex_assay <- regex_norm@metabolite_data$LCMS_Pos

  expect_equal(regex_assay$Name, c("M1", "M2"))
  expect_equal(regex_assay$Sample_2[1], 7.5)
  expect_equal(regex_assay$Sample_4[2], 37.5)
  expect_true(isTRUE(regex_norm@args$ITSDNormalization$applied))
  expect_identical(regex_norm@args$ITSDNormalization$itsd_counts_per_assay$LCMS_Pos, 1L)

  manual_norm <- suppressWarnings(
    normaliseUntransformedData(
      the_object,
      method = "ITSD",
      itsd_feature_ids = list(LCMS_Pos = "ISTD1"),
      itsd_aggregation = "mean",
      remove_itsd_after_norm = FALSE
    )
  )
  manual_assay <- manual_norm@metabolite_data$LCMS_Pos

  expect_equal(sort(manual_assay$Name), c("ISTD1", "M1", "M2"))
  expect_identical(manual_norm@args$ITSDNormalization$itsd_aggregation, "mean")
  expect_false(isTRUE(manual_norm@args$ITSDNormalization$removed_itsd))

  expect_error(
    normaliseUntransformedData(the_object, method = "PQN"),
    "only supports method = 'ITSD'"
  )

  expect_error(
    normaliseUntransformedData(the_object, method = "ITSD", itsd_aggregation = "mode"),
    "must be one of"
  )
})

test_that("normaliseBetweenSamples and cleanDesignMatrix preserve normalization dispatch and design cleanup", {
  reordered_object <- newMetabNormCoreObject(
    assay_list = list(
      LCMS_Pos = tibble::tibble(
        Name = c("M1", "M2"),
        annotation_id = c("alpha", "beta"),
        Sample_2 = c(20, 40),
        Sample_1 = c(10, 30)
      )
    ),
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_3", "Sample_2"),
      Batch = c("A", "B", "C"),
      stringsAsFactors = FALSE
    )
  )

  cleaned <- cleanDesignMatrix(reordered_object)
  expect_identical(cleaned@design_matrix$Run, c("Sample_2", "Sample_1"))
  expect_identical(cleaned@design_matrix$Batch, c("C", "A"))

  unmatched_object <- newMetabNormCoreObject(
    assay_list = list(
      LCMS_Pos = tibble::tibble(
        Name = c("M1", "M2"),
        annotation_id = c("alpha", "beta"),
        Sample_1 = c(10, 20),
        Sample_2 = c(30, 40)
      )
    ),
    design_matrix = data.frame(
      Run = c("Missing_1", "Missing_2"),
      Batch = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )

  expect_warning(
    cleanDesignMatrix(unmatched_object),
    "No sample columns identified",
    ignore.case = TRUE
  )

  quantile_object <- newMetabNormCoreObject(
    assay_list = list(
      LCMS_Pos = tibble::tibble(
        Name = c("M1", "M2", "M3"),
        annotation_id = c("a", "b", "c"),
        Sample_1 = c(10, 20, 30),
        Sample_2 = c(11, 21, 31),
        Sample_3 = c(12, 22, 32),
        Sample_4 = c(13, 23, 33)
      )
    ),
    design_matrix = data.frame(
      Run = c("Sample_4", "Sample_2", "Sample_1", "Sample_3"),
      stringsAsFactors = FALSE
    ),
    args = list(
      normaliseBetweenSamples = list(
        normalisation_method = "quantile"
      )
    )
  )

  quantile_norm <- suppressWarnings(normaliseBetweenSamples(quantile_object))
  quantile_assay <- quantile_norm@metabolite_data$LCMS_Pos

  expect_equal(quantile_assay$Sample_1, c(11.5, 21.5, 31.5))
  expect_equal(quantile_assay$Sample_4, c(11.5, 21.5, 31.5))
  expect_identical(quantile_norm@design_matrix$Run, c("Sample_1", "Sample_2", "Sample_3", "Sample_4"))

  none_object <- newMetabNormCoreObject(
    assay_list = list(
      LCMS_Pos = tibble::tibble(
        Name = c("M1", "M2", "M3"),
        annotation_id = c("a", "b", "c"),
        Sample_1 = c(10, 20, 30),
        Sample_2 = c(11, 21, 31),
        Sample_3 = c(12, 22, 32),
        Sample_4 = c(13, 23, 33)
      )
    ),
    design_matrix = data.frame(
      Run = c("Sample_4", "Sample_2", "Sample_1", "Sample_3"),
      stringsAsFactors = FALSE
    ),
    args = list(
      normaliseBetweenSamples = list(
        normalisation_method = "none"
      )
    )
  )

  none_norm <- suppressWarnings(normaliseBetweenSamples(none_object))
  expect_equal(none_norm@metabolite_data$LCMS_Pos$Sample_1, c(10, 20, 30))
  expect_equal(none_norm@metabolite_data$LCMS_Pos$Sample_4, c(13, 23, 33))
})
