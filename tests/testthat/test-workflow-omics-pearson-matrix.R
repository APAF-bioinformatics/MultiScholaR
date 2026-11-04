# The lipid and metabolite Pearson methods on one matrix call: 1:1 with the per-pair kernel they replaced.
library(testthat)

perPairReference <- function(assay, id_column, pairs, pair_kernel) {
  sample_cols <- setdiff(colnames(assay), id_column)
  long <- tidyr::pivot_longer(assay, tidyr::all_of(sample_cols), names_to = "Run", values_to = "abundance")
  suppressWarnings(purrr::map2_dbl(pairs$Run.x, pairs$Run.y, function(x, y) {
    pair_kernel(dplyr::filter(long, .data$Run %in% c(x, y)), id_column, "Run", "abundance")
  }))
}

smallAssay <- function(id_column, n_features = 500L, n_samples = 24L, seed = 84L) {
  set.seed(seed)
  base <- rnorm(n_features, 10, 2)
  values <- matrix(rep(base, n_samples) + rnorm(n_features * n_samples, 0, 0.5), nrow = n_features)
  values[sample(length(values), length(values) %/% 10L)] <- NA_real_
  wide <- seq_len(min(20L, n_features))
  values[wide, 3L] <- values[wide, 3L] * 1e6 # a wide scale
  values[, 7L] <- 5 # no variance
  values[, 9L] <- NA_real_ # one observation
  values[1L, 9L] <- 1
  values[, 11L] <- NA_real_ # none
  colnames(values) <- sprintf("S%02d", seq_len(n_samples))
  assay <- tibble::as_tibble(values)
  assay[[id_column]] <- sprintf("F%04d", seq_len(n_features))
  assay[c(id_column, colnames(values))]
}

smallDesign <- function(n_samples = 24L) {
  data.frame(Run = sprintf("S%02d", seq_len(n_samples)), group = rep(c("A", "B"), length.out = n_samples),
    stringsAsFactors = FALSE)
}

metabObject <- function(assays) {
  methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = smallDesign(), sample_id = "Run",
    metabolite_id_column = "feature", technical_replicate_id = "group", args = list())
}

lipidObject <- function(assays) {
  methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = smallDesign(), sample_id = "Run",
    lipid_id_column = "feature", technical_replicate_id = "group", args = list())
}

quietly <- function(expr) suppressWarnings(suppressMessages(expr))

cases <- list(
  metabolomics = list(object = metabObject, kernel = calculateMetabolitePairCorrelation),
  lipidomics = list(object = lipidObject, kernel = calculateLipidPairCorrelation)
)

for (omic in names(cases)) {
  test_that(paste(omic, "Pearson equals the per-pair kernel, bit for bit, on every pair"), {
    assays <- list(
      first = smallAssay("feature"),
      second = smallAssay("feature", n_features = 300L, seed = 85L),
      duplicated = local({
        assay <- smallAssay("feature", n_features = 50L, seed = 86L)
        assay$feature[[2L]] <- assay$feature[[1L]]
        assay
      }),
      one_feature = smallAssay("feature", n_features = 1L, seed = 87L)
    )
    results <- quietly(pearsonCorForSamplePairs(cases[[omic]]$object(assays)))
    expect_named(results, names(assays))
    for (assay_name in names(assays)) {
      observed <- results[[assay_name]]
      reference <- perPairReference(assays[[assay_name]], "feature", observed, cases[[omic]]$kernel)
      expect_identical(observed$pearson_correlation, reference, info = assay_name)
    }
    expect_gt(sum(!is.na(results$first$pearson_correlation)), 100L)
    expect_gt(sum(is.na(results$first$pearson_correlation)), 0L)
    expect_true(all(is.na(results$duplicated$pearson_correlation)))
    expect_true(all(is.na(results$one_feature$pearson_correlation)))
  })
}

test_that("neither method calls the per-pair kernel", {
  fail <- function(...) stop("per-pair kernel called")
  testthat::local_mocked_bindings(calculateMetabolitePairCorrelation = fail, calculateLipidPairCorrelation = fail,
    .package = "MultiScholaR")
  assays <- list(first = smallAssay("feature", n_features = 50L))
  for (omic in names(cases)) {
    results <- quietly(pearsonCorForSamplePairs(cases[[omic]]$object(assays)))
    expect_gt(sum(!is.na(results$first$pearson_correlation)), 50L)
  }
})
