# The lipid and metabolite RUV sweeps and corrections carry proteomics' rule (KIT-005): RUV-III cannot estimate more factors
# than it has controls observed in every sample, so a sweep skips such a percentage and the correction refuses such a k.
library(testthat)

controlDesign <- function(n_samples = 24L) {
  biological <- rep(sprintf("B%02d", seq_len(n_samples %/% 2L)), each = 2L)
  data.frame(Run = sprintf("S%02d", seq_len(n_samples)), group = rep(c("A", "B"), each = n_samples %/% 2L), tech_rep_group = biological,
    stringsAsFactors = FALSE)
}

controlAssay <- function(id_column, seed, complete_features, n_features = 300L, n_samples = 24L) {
  set.seed(seed)
  values <- matrix(rnorm(n_features * n_samples, 15, 2), nrow = n_features) + rep(rnorm(n_samples, 0, 0.3), each = n_features)
  colnames(values) <- sprintf("S%02d", seq_len(n_samples))
  gaps <- seq_len(n_features)
  gaps <- gaps[gaps > complete_features]
  values[cbind(gaps, (gaps %% n_samples) + 1L)] <- NA_real_
  assay <- tibble::as_tibble(values)
  assay[[id_column]] <- sprintf("F%04d", seq_len(n_features))
  assay[c(id_column, colnames(values))]
}

controlCases <- function() {
  list(
    metab = list(sweep = MultiScholaR:::runPerAssayRuvOptimization, make = function(assays) methods::new("MetaboliteAssayData",
      metabolite_data = assays, design_matrix = controlDesign(), sample_id = "Run", group_id = "group",
      metabolite_id_column = "feature", technical_replicate_id = "tech_rep_group", args = list())),
    lipid = list(sweep = MultiScholaR:::runLipidPerAssayRuvOptimization, make = function(assays) methods::new("LipidomicsAssayData",
      lipid_data = assays, design_matrix = controlDesign(), sample_id = "Run", group_id = "group",
      lipid_id_column = "feature", technical_replicate_id = "tech_rep_group", args = list())))
}

controlSweep <- function(case, object) {
  suppressWarnings(suppressMessages(case$sweep(object, ruv_mode = "automatic",
    params = list(percentage_min = 40, percentage_max = 42, ruv_grouping_variable = "group"))))
}

for (omic in names(controlCases())) {
  test_that(paste(omic, "sweep skips a percentage whose elbow k exceeds its complete controls"), {
    case <- controlCases()[[omic]]
    object <- case$make(list(one = controlAssay("feature", 21L, complete_features = 4L), two = controlAssay("feature", 22L, complete_features = 250L)))
    local_mocked_bindings(findBestKElbow = function(...) 5L, .package = "MultiScholaR")
    result <- controlSweep(case, object)
    one <- result$one$optimization_results
    expect_true(all(one$status == "skipped_k_exceeds_complete_controls"))
    expect_false(isTRUE(result$one$success))
    # The second assay has many complete controls: the same k is viable there, assay by assay.
    expect_true(all(result$two$optimization_results$status == "ok"))
    expect_true(isTRUE(result$two$success))
  })

  test_that(paste(omic, "correction refuses a k above its complete controls and accepts one within them"), {
    case <- controlCases()[[omic]]
    object <- case$make(list(one = controlAssay("feature", 21L, complete_features = 4L)))
    ctrl <- list(one = stats::setNames(rep(TRUE, 300L), sprintf("F%04d", seq_len(300L))))
    expect_error(suppressWarnings(suppressMessages(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = 5L, ctrl = ctrl))),
      class = "omics_ruv_insufficient_controls")
    corrected <- suppressWarnings(suppressMessages(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = 1L, ctrl = ctrl)))
    expect_s4_class(corrected, class(object))
  })

  test_that(paste(omic, "rule changes nothing on complete data"), {
    case <- controlCases()[[omic]]
    object <- case$make(list(one = controlAssay("feature", 21L, complete_features = 300L)))
    result <- controlSweep(case, object)
    expect_true(all(result$one$optimization_results$status == "ok"))
    expect_true(isTRUE(result$one$success))
  })
}
