library(testthat)

test_that("ci impact source map owns every tracked R source path", {
  impact_map <- ci_impact_map()
  repo_root <- testthat::test_path("..", "..")
  r_files <- system2("git", c("-C", repo_root, "ls-files", "R/*.R"), stdout = TRUE)
  expect_gt(length(r_files), 0L)

  unmatched <- vapply(r_files, function(path) {
    !any(vapply(impact_map$rules, function(rule) {
      any(vapply(unlist(rule$paths, use.names = FALSE), function(pattern) {
        grepl(pattern, path, perl = TRUE)
      }, logical(1)))
    }, logical(1)))
  }, logical(1))

  expect_false(any(unmatched), info = paste(r_files[unmatched], collapse = "\n"))
})

test_that("ci impact source map routes canonical source families", {
  cases <- list(
    lipid_design = list(path = "R/mod_lipid_design_builder_helpers.R", module = "lipidomics/design", rule = "lipidomics-design"),
    prot_limpa = list(path = "R/func_prot_limpa.R", module = "proteomics/normalization", rule = "proteomics-normalization-limpa"),
    metab_da = list(path = "R/mod_metab_da_server.R", module = "metabolomics/differential_abundance", rule = "metabolomics-da"),
    shared_state = list(path = "R/utils_workflow_state.R", module = "proteomics/all", rule = "shared-state-and-workflow")
  )

  for (case in cases) {
    payload <- ci_impact_run(case$path)
    expect_true(case$module %in% ci_impact_modules(payload), info = case$path)
    expect_true(case$rule %in% ci_impact_rule_ids(payload), info = case$path)
  }
})

test_that("ci impact source map merges narrow consequences and broad escalation safely", {
  payload <- ci_impact_run(c("R/mod_lipid_da.R", "R/utils_workflow_state.R"))

  expect_identical(payload$impact_level, "broad")
  expect_true("lipidomics/differential_abundance" %in% ci_impact_modules(payload))
  expect_true("proteomics/all" %in% ci_impact_modules(payload))
  expect_true(isTRUE(payload$run_cross_omic))
})
