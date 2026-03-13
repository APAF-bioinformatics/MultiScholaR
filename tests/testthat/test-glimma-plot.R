library(testthat)
library(MultiScholaR)

test_that("generateProtDAVolcanoPlotGlimma works with captured session data", {
  fixture_dir <- testthat::test_path("..", "testdata", "glimma_fixtures")
  
  if (!dir.exists(fixture_dir)) {
    skip(paste0("Glimma test fixtures directory not found at: ", fixture_dir))
  }
  
  snapshots <- list.files(fixture_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(snapshots) == 0) {
    skip("No Glimma snapshots found. Run the app with start_glimma_capture() to create one.")
  }
  
  # Use the most recent snapshot
  latest_snapshot <- rev(sort(snapshots))[1]
  message("Testing with snapshot: ", latest_snapshot)
  
  data <- readRDS(latest_snapshot)
  
  # Run the function
  p <- do.call(generateProtDAVolcanoPlotGlimma, data)
  
  # Assertions
  expect_true(!is.null(p), info = "Function returned NULL")
  expect_s3_class(p, "htmlwidget")
  expect_s3_class(p, "glimmaXY")
})
