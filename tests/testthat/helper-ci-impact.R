`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

ci_impact_script <- function() {
  testthat::test_path("..", "..", "tools", "ci", "detect-impact.R")
}

ci_impact_run <- function(changed_files, output = tempfile(fileext = ".json"), summary = NULL, extra_args = character()) {
  args <- c(
    ci_impact_script(),
    "--output",
    output
  )
  if (length(changed_files) > 0L) {
    args <- c(args, "--changed-files", paste(changed_files, collapse = ","))
  }
  if (!is.null(summary)) {
    args <- c(args, "--summary", summary)
  }
  args <- c(args, extra_args)
  result <- system2(file.path(R.home("bin"), "Rscript"), args, stdout = TRUE, stderr = TRUE)
  testthat::expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
  jsonlite::read_json(output, simplifyVector = FALSE)
}

ci_impact_module_rows <- function(payload) {
  payload$module_matrix$include %||% list()
}

ci_impact_e2e_rows <- function(payload) {
  payload$e2e_matrix$include %||% list()
}

ci_impact_browser_rows <- function(payload) {
  payload$browser_matrix$include %||% list()
}

ci_impact_modules <- function(payload) {
  vapply(ci_impact_module_rows(payload), function(row) paste(row$omic, row$module, sep = "/"), character(1))
}

ci_impact_e2e_filters <- function(payload) {
  vapply(ci_impact_e2e_rows(payload), `[[`, character(1), "filter")
}

ci_impact_rule_ids <- function(payload) {
  unique(vapply(payload$matched_rules, `[[`, character(1), "id"))
}

ci_impact_workflow_text <- function() {
  paste(readLines(testthat::test_path("..", "..", ".github", "workflows", "module-ci.yml"), warn = FALSE), collapse = "\n")
}

ci_impact_map <- function() {
  jsonlite::read_json(testthat::test_path("..", "..", "tools", "ci", "impact-map.json"), simplifyVector = FALSE)
}
