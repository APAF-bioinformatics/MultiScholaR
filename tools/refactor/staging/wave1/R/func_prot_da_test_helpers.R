#' Start Glimma Data Capture
#' @export
start_glimma_capture <- function() {
  options(multischolar.capture_glimma_data = TRUE)
  message("MultiScholaR: Glimma data capture ENABLED.")
  message("Snapshots will be saved to: tests/testdata/glimma_fixtures/")
}

#' Stop Glimma Data Capture
#' @export
stop_glimma_capture <- function() {
  options(multischolar.capture_glimma_data = FALSE)
  message("MultiScholaR: Glimma data capture DISABLED.")
}

