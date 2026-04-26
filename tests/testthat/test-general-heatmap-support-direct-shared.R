# fidelity-coverage-compare: shared
library(testthat)

localBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

test_that("heatmap support helpers preserve annotation filtering and file outputs", {
  localBinding(
    parent.env(environment(getSamplesCorrelationHeatMap)),
    "rowAnnotation",
    ComplexHeatmap::rowAnnotation
  )
  localBinding(
    parent.env(environment(getSamplesCorrelationHeatMap)),
    "Legend",
    ComplexHeatmap::Legend
  )
  localBinding(
    parent.env(environment(getProteinsHeatMap)),
    "Heatmap",
    function(...) {
      args <- list(...)
      args$core_utilisation_columns <- NULL
      do.call(ComplexHeatmap::Heatmap, args)
    }
  )

  correlation_matrix <- matrix(
    c(1, 0.8, 0.8, 1),
    nrow = 2,
    dimnames = list(c("S1", "S2"), c("S1", "S2"))
  )
  protein_matrix <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S2", "S1"))
  )
  metadata_tbl <- data.frame(
    Run = c("S1", "S2", "S3"),
    is_HEK = c(FALSE, FALSE, TRUE),
    condition = c("A", "B", "A"),
    batch = c("X", "Y", "X"),
    stringsAsFactors = FALSE
  )
  colour_rules <- list(
    Condition = c(A = "#1b9e77", B = "#d95f02"),
    Batch = c(X = "#7570b3", Y = "#e7298a")
  )

  sample_heatmap <- suppressMessages(suppressWarnings(getSamplesCorrelationHeatMap(
    correlation_matrix = correlation_matrix,
    metadata_tbl = metadata_tbl,
    is_HEK_column = is_HEK,
    metadata_column_labels = c("Condition", "Batch"),
    metadata_column_selected = c("condition", "batch"),
    colour_rules = colour_rules,
    columns_to_exclude = character(),
    sample_id_column = Run,
    raster_device = "CairoPNG"
  )))

  protein_heatmap <- suppressMessages(suppressWarnings(getProteinsHeatMap(
    protein_matrix = protein_matrix,
    metadata_tbl = metadata_tbl,
    is_HEK_column = is_HEK,
    metadata_column_selected = c("condition", "batch"),
    metadata_column_labels = c("Condition", "Batch"),
    colour_rules = colour_rules,
    columns_to_exclude = character(),
    sort_by_sample_id = TRUE,
    sample_id_column = Run
  )))

  output_dir <- tempfile("heatmap-support-")
  dir.create(output_dir)

  saved_with_clusters <- suppressMessages(suppressWarnings(save_heatmap_products(
    heatmap_obj = protein_heatmap$heatmap,
    row_clusters = c(P1 = 1, P2 = 2),
    params_list = list(cluster_rows = TRUE, selected_assay = "Combined"),
    output_dir = output_dir,
    file_prefix = "demo"
  )))

  saved_without_clusters <- suppressMessages(suppressWarnings(save_heatmap_products(
    heatmap_obj = sample_heatmap$heatmap,
    row_clusters = NULL,
    params_list = list(cluster_rows = FALSE, selected_assay = "Correlation"),
    output_dir = output_dir,
    file_prefix = "demo_no_clusters"
  )))

  expect_true(is.list(sample_heatmap))
  expect_s4_class(sample_heatmap$heatmap, "Heatmap")
  expect_length(sample_heatmap$legend, 2L)

  expect_true(is.list(protein_heatmap))
  expect_s4_class(protein_heatmap$heatmap, "Heatmap")
  expect_length(protein_heatmap$legend, 2L)

  expect_true(is.list(saved_with_clusters))
  expect_true(is.list(saved_without_clusters))
  expect_true(file.exists(file.path(output_dir, "Heatmap", "demo_heatmap_methods.md")))
  expect_true(file.exists(file.path(output_dir, "Heatmap", "demo_clusters.csv")))
  expect_true(file.exists(file.path(output_dir, "Heatmap", "demo_no_clusters_heatmap_methods.md")))
  expect_false(file.exists(file.path(output_dir, "Heatmap", "demo_no_clusters_clusters.csv")))
})
