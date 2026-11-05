# The lipid and metabolite RUV sweeps impute the canonical-correlation input once per assay, as proteomics does (KIT-005):
# same results as imputing at every percentage, including on data with missing values.
library(testthat)

ruvAssay <- function(id_column, seed, n_features = 300L, n_samples = 24L, missing = 0.05) {
  set.seed(seed)
  values <- matrix(rnorm(n_features * n_samples, 15, 2), nrow = n_features) + rep(rnorm(n_samples, 0, 0.3), each = n_features)
  values[seq_len(30L), seq(1L, n_samples, by = 2L)] <- values[seq_len(30L), seq(1L, n_samples, by = 2L)] + 3
  if (missing > 0) values[sample(length(values), floor(missing * length(values)))] <- NA_real_
  colnames(values) <- sprintf("S%02d", seq_len(n_samples))
  assay <- tibble::as_tibble(values)
  assay[[id_column]] <- sprintf("F%04d", seq_len(n_features))
  assay[c(id_column, colnames(values))]
}

ruvDesign <- function() {
  data.frame(Run = sprintf("S%02d", seq_len(24L)), group = rep(c("A", "B"), length.out = 24L), stringsAsFactors = FALSE)
}

ruvMetabObject <- function(assays, args = list()) {
  methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = ruvDesign(), sample_id = "Run", group_id = "group",
    metabolite_id_column = "feature", technical_replicate_id = "group", args = args)
}

ruvLipidObject <- function(assays, args = list()) {
  methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = ruvDesign(), sample_id = "Run", group_id = "group",
    lipid_id_column = "feature", technical_replicate_id = "group", args = args)
}

ruvCases <- function() {
  assays <- function(missing) list(one = ruvAssay("feature", 11L, missing = missing), two = ruvAssay("feature", 12L, missing = missing))
  list(
    metab = list(make = ruvMetabObject, sweep = MultiScholaR:::runPerAssayRuvOptimization, assays = assays),
    lipid = list(make = ruvLipidObject, sweep = MultiScholaR:::runLipidPerAssayRuvOptimization, assays = assays))
}

ruvSweep <- function(case, object) {
  suppressWarnings(suppressMessages(case$sweep(object, ruv_mode = "automatic",
    params = list(percentage_min = 2, percentage_max = 6, ruv_grouping_variable = "group"))))
}

# A sweep result without its plot objects, plus the data each plot was drawn from.
ruvComparable <- function(result) {
  lapply(result, function(assay) {
    lapply(assay, function(value) {
      if (inherits(value, "ggplot")) return(value$data)
      if (is.list(value) && !is.data.frame(value)) {
        return(lapply(value, function(inner) if (inherits(inner, "ggplot")) inner$data else inner))
      }
      value
    })
  })
}

for (omic in names(ruvCases())) {
  test_that(paste(omic, "sweep imputes each assay once"), {
    case <- ruvCases()[[omic]]
    calls <- 0L
    original <- MultiScholaR:::imputeNipalsMatrix
    local_mocked_bindings(imputeNipalsMatrix = function(...) { calls <<- calls + 1L; original(...) }, .package = "MultiScholaR")
    ruvSweep(case, case$make(case$assays(0.05)))
    expect_identical(calls, 2L)
  })

  test_that(paste(omic, "sweep is 1:1 with imputing at every percentage, with missing values"), {
    case <- ruvCases()[[omic]]
    object <- case$make(case$assays(0.05))
    new <- ruvSweep(case, object)
    local_mocked_bindings(prepareOmicsRuvCancorObject = function(theObject, ...) theObject, .package = "MultiScholaR")
    old <- ruvSweep(case, object)
    expect_identical(names(new), names(old))
    expect_identical(ruvComparable(new), ruvComparable(old))
    expect_true(all(vapply(new, function(assay) isTRUE(assay$success), logical(1))))
  })

  test_that(paste(omic, "preparation leaves complete data untouched and honours the object's components"), {
    case <- ruvCases()[[omic]]
    complete <- case$make(case$assays(0))
    slot_name <- if (omic == "metab") "metabolite_data" else "lipid_data"
    id_slot <- if (omic == "metab") "metabolite_id_column" else "lipid_id_column"
    expect_identical(MultiScholaR:::prepareOmicsRuvCancorObject(complete, slot_name, id_slot), complete)
    seen_ncomp <- NULL
    original <- MultiScholaR:::imputeNipalsMatrix
    local_mocked_bindings(imputeNipalsMatrix = function(X, ncomp, ...) { seen_ncomp <<- ncomp; original(X, ncomp, ...) }, .package = "MultiScholaR")
    with_args <- case$make(case$assays(0.05), args = list(ruvCancor = list(num_components_to_impute = 3)))
    prepared <- MultiScholaR:::prepareOmicsRuvCancorObject(with_args, slot_name, id_slot)
    expect_identical(seen_ncomp, 3)
    expect_false(anyNA(methods::slot(prepared, slot_name)[["one"]]))
    observed <- !is.na(as.matrix(methods::slot(with_args, slot_name)[["one"]][-1L]))
    expect_equal(as.matrix(methods::slot(prepared, slot_name)[["one"]][-1L])[observed],
      as.matrix(methods::slot(with_args, slot_name)[["one"]][-1L])[observed])
  })
}
