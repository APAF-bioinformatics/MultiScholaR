# Load the test-owned browser contract for standalone CI runners.

source(
  file.path("tests", "testthat", "helper-00-e2e-browser-preflight.R"),
  local = TRUE
)
