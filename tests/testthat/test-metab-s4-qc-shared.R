# fidelity-coverage-compare: shared
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
      technical_replicate_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "feature",
      annotation_id_column = "annotation_id",
      internal_standard_regex = "^ISTD",
      design_matrix = data.frame(),
      sample_id = "Run",
      args = list(),
      technical_replicate_id = "TechRep"
    )
  )
}

if (!methods::isGeneric("metaboliteIntensityFiltering")) {
  methods::setGeneric(
    "metaboliteIntensityFiltering",
    function(
      theObject,
      metabolites_intensity_cutoff_percentile = NULL,
      metabolites_proportion_of_samples_below_cutoff = NULL
    ) {
      standardGeneric("metaboliteIntensityFiltering")
    }
  )
}

if (!methods::isGeneric("resolveDuplicateFeatures")) {
  methods::setGeneric(
    "resolveDuplicateFeatures",
    function(theObject, itsd_pattern_columns = NULL) {
      standardGeneric("resolveDuplicateFeatures")
    }
  )
}

if (!methods::isGeneric("filterSamplesByMetaboliteCorrelationThreshold")) {
  methods::setGeneric(
    "filterSamplesByMetaboliteCorrelationThreshold",
    function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = 0.5) {
      standardGeneric("filterSamplesByMetaboliteCorrelationThreshold")
    }
  )
}

if (!methods::isGeneric("pearsonCorForSamplePairs")) {
  methods::setGeneric(
    "pearsonCorForSamplePairs",
    function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
      standardGeneric("pearsonCorForSamplePairs")
    }
  )
}

if (!methods::isGeneric("plotPearson")) {
  methods::setGeneric(
    "plotPearson",
    function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {
      standardGeneric("plotPearson")
    }
  )
}

newMetabS4QcObject <- function(args = list()) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        feature = c("M1", "M2", "M3"),
        annotation_id = c("alpha", "beta", "gamma"),
        Sample_1 = c(1, 2, 3),
        Sample_2 = c(1, 2, 3),
        Sample_3 = c(3, 2, 1),
        Sample_4 = c(1, 2, 3)
      ),
      LCMS_Neg = tibble::tibble(
        feature = c("N1", "N2", "N3"),
        annotation_id = c("delta", "epsilon", "zeta"),
        Sample_1 = c(5, 6, 7),
        Sample_2 = c(5, 6, 7),
        Sample_3 = c(10, 9, 8),
        Sample_4 = c(10, 9, 8)
      )
    ),
    metabolite_id_column = "feature",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "^ISTD",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      TechRep = c("pair_a", "pair_a", "pair_b", "pair_b"),
      AltGroup = c("grp1", "grp1", "grp2", "grp2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    technical_replicate_id = "TechRep",
    args = args
  )
}

newMetabIntensityObject <- function(args = list()) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        feature = c("M1", "M2", "M3"),
        annotation_id = c("alpha", "beta", "gamma"),
        Sample_1 = c(10, 100, 60),
        Sample_2 = c(20, 200, NA)
      )
    ),
    metabolite_id_column = "feature",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "^ISTD",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      TechRep = c("pair_a", "pair_a"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    technical_replicate_id = "TechRep",
    args = args
  )
}

newMetabResolveObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        database_identifier = c("ISTD_LOW", "ISTD_HIGH", "DB_DUP", "DB_DUP", "DB_KEEP"),
        metabolite = c("ISTD Mix", "ISTD Mix", "Non-ITSD low", "Non-ITSD high", "Unique"),
        annotation_id = c("ISTD_alpha", "ISTD_alpha", "feature_low", "feature_high", "feature_unique"),
        Sample_1 = c(10, 30, 5, 50, 1),
        Sample_2 = c(20, 40, 15, 60, 2)
      )
    ),
    metabolite_id_column = "database_identifier",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "^ISTD",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    technical_replicate_id = "TechRep",
    args = list()
  )
}

newMetabEmptyObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "^ISTD",
    design_matrix = data.frame(
      Run = "Sample_1",
      TechRep = "pair_a",
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    technical_replicate_id = "TechRep",
    args = list()
  )
}

test_that("metaboliteIntensityFiltering preserves config cleanup, helper filtering, and validation branches", {
  filtered <- suppressWarnings(
    metaboliteIntensityFiltering(
      newMetabIntensityObject(
        args = list(
          metaboliteIntensityFiltering = list(
            metabolites_intensity_cutoff_percentile = "50 # config percentile",
            metabolites_proportion_of_samples_below_cutoff = "0.6 # config cutoff"
          )
        )
      )
    )
  )

  filtered_assay <- filtered@metabolite_data$LCMS_Pos

  expect_s3_class(filtered_assay, "tbl_df")
  expect_identical(filtered_assay$feature, c("M2", "M3"))
  expect_equal(filtered_assay$Sample_1, c(100, 60))
  expect_equal(filtered_assay$Sample_2, c(200, NA))
  expect_null(filtered@args$metaboliteIntensityFiltering$metabolites_intensity_cutoff_percentile)
  expect_null(filtered@args$metaboliteIntensityFiltering$metabolites_proportion_of_samples_below_cutoff)

  expect_warning(
    metaboliteIntensityFiltering(
      newMetabEmptyObject(),
      metabolites_intensity_cutoff_percentile = 50,
      metabolites_proportion_of_samples_below_cutoff = 0.6
    ),
    "no assays",
    ignore.case = TRUE
  )

  expect_error(
    suppressWarnings(
      metaboliteIntensityFiltering(
        newMetabIntensityObject(
          args = list(
            metaboliteIntensityFiltering = list(
              metabolites_intensity_cutoff_percentile = "bad_value # config percentile",
              metabolites_proportion_of_samples_below_cutoff = "0.6"
            )
          )
        )
      )
    ),
    "Failed to convert cleaned metabolites_intensity_cutoff_percentile"
  )

  helper_filtered <- metaboliteIntensityFilteringHelper(
    tibble::tibble(
      feature = c("M1", "M2"),
      Sample_1 = c(1, 10),
      Sample_2 = c(2, 20)
    ),
    min_metabolite_intensity_threshold = 5,
    metabolites_proportion_of_samples_below_cutoff = 0.6,
    metabolite_id_column = "feature"
  )

  expect_identical(helper_filtered$feature, "M2")
})

test_that("resolveDuplicateFeatures preserves ITSD renaming, non-ITSD deduplication, and duplicate discovery", {
  input_object <- newMetabResolveObject()

  duplicates <- suppressMessages(findMetabDuplicateFeatureIDs(input_object))
  expect_named(duplicates, "LCMS_Pos")
  expect_identical(duplicates$LCMS_Pos$database_identifier, "DB_DUP")
  expect_identical(duplicates$LCMS_Pos$count, 2L)

  resolved <- suppressWarnings(suppressMessages(resolveDuplicateFeatures(input_object)))
  resolved_assay <- resolved@metabolite_data$LCMS_Pos

  expect_equal(nrow(resolved_assay), 4L)
  expect_false(".original_row_id" %in% names(resolved_assay))
  expect_equal(
    sort(resolved_assay$database_identifier),
    c("DB_DUP", "DB_KEEP", "ISTD Mix_1", "ISTD Mix_2")
  )

  resolved_by_intensity <- resolveDuplicateFeaturesByIntensity(
    tibble::tibble(
      id = c("A", "A", "B"),
      Sample_1 = c(1, 5, 2),
      Sample_2 = c(1, 6, 2)
    ),
    id_col = "id",
    sample_cols = c("Sample_1", "Sample_2")
  )

  expect_identical(resolved_by_intensity$id, c("A", "B"))
  expect_equal(resolved_by_intensity$Sample_1, c(5, 2))
})

test_that("pearsonCorForSamplePairs preserves grouping overrides, regex filtering, and helper fallbacks", {
  the_object <- newMetabS4QcObject(args = list(pearsonCorForSamplePairs = list()))

  results_default <- suppressWarnings(
    suppressMessages(pearsonCorForSamplePairs(the_object, tech_rep_remove_regex = "pool"))
  )

  expect_named(results_default, c("LCMS_Pos", "LCMS_Neg"))
  expect_equal(results_default$LCMS_Pos$pearson_correlation, c(1, -1))
  expect_identical(results_default$LCMS_Pos$TechRep, c("pair_a", "pair_b"))

  results_override <- suppressWarnings(
    suppressMessages(pearsonCorForSamplePairs(the_object, correlation_group = "AltGroup"))
  )
  expect_identical(results_override$LCMS_Pos$AltGroup, c("grp1", "grp2"))
  expect_equal(results_override$LCMS_Pos$pearson_correlation, c(1, -1))

  expect_error(
    suppressWarnings(suppressMessages(pearsonCorForSamplePairs(the_object, correlation_group = "MissingGroup"))),
    "MissingGroup"
  )

  pair_table <- tibble::tibble(
    feature = c("M1", "M1", "M2", "M2"),
    Run = c("S1", "S2", "S1", "S2"),
    abundance = c(1, 1, 2, 2)
  )
  expect_equal(calculateMetabolitePairCorrelation(pair_table, "feature", "Run", "abundance"), 1)

  invalid_pair_table <- tibble::tibble(
    feature = c("M1", "M1", "M1"),
    Run = c("S1", "S2", "S3"),
    abundance = c(1, 1, 1)
  )
  expect_warning(
    invalid_correlation <- calculateMetabolitePairCorrelation(
      invalid_pair_table,
      "feature",
      "Run",
      "abundance"
    ),
    "does not contain exactly two samples"
  )
  expect_true(is.na(invalid_correlation))
})

test_that("filterSamplesByMetaboliteCorrelationThreshold and plotPearson preserve filtering and histogram output", {
  the_object <- newMetabS4QcObject(args = list(pearsonCorForSamplePairs = list()))

  pearson_results <- suppressWarnings(
    suppressMessages(pearsonCorForSamplePairs(the_object))
  )

  filtered <- suppressWarnings(
    suppressMessages(
      filterSamplesByMetaboliteCorrelationThreshold(
        the_object,
        pearson_correlation_per_pair = pearson_results,
        min_pearson_correlation_threshold = 0.5
      )
    )
  )

  expect_identical(filtered@design_matrix$Run, c("Sample_1", "Sample_2"))
  expect_identical(
    colnames(filtered@metabolite_data$LCMS_Pos),
    c("feature", "annotation_id", "Sample_1", "Sample_2")
  )

  expect_error(
    filterSamplesByMetaboliteCorrelationThreshold(
      the_object,
      pearson_correlation_per_pair = NULL
    ),
    "must be a list"
  )

  pearson_plots <- suppressWarnings(
    suppressMessages(plotPearson(the_object))
  )

  expect_named(pearson_plots, c("LCMS_Pos", "LCMS_Neg"))
  expect_s3_class(pearson_plots$LCMS_Pos, "ggplot")
  expect_s3_class(pearson_plots$LCMS_Neg, "ggplot")
})
