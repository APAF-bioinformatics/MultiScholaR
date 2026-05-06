library(testthat)

test_that("MCI-024.1 sample identity sentinel follows declared drops and renames", {
  catalog <- read_module_ci_fixture_packs()
  pack <- module_ci_fixture_packs(catalog, omic = "proteomics", import_format = "dia", fixture_class = "happy_path")[[1]]
  oracle <- module_ci_read_pack_oracle(pack)
  samples <- module_ci_as_character(oracle$sample_order)

  import_digest <- module_ci_sentinel_digest(samples, kind = "sample", metadata = list(stage = "import"))
  design_digest <- module_ci_sentinel_digest(samples, kind = "sample", metadata = list(stage = "design"))
  expect_no_error(module_ci_sentinel_assert_identity(import_digest, design_digest, label = "sample"))

  filtered_samples <- samples[-1]
  expect_no_error(module_ci_sentinel_assert_identity(
    import_digest,
    filtered_samples,
    expected_dropped = samples[[1]],
    label = "sample"
  ))
  expect_failure(module_ci_sentinel_assert_identity(import_digest, filtered_samples, label = "sample"))

  renamed_samples <- samples
  renamed_samples[[2]] <- paste0(samples[[2]], "_RENAMED")
  expect_no_error(module_ci_sentinel_assert_identity(
    import_digest,
    renamed_samples,
    expected_renamed = stats::setNames(renamed_samples[[2]], samples[[2]]),
    label = "sample"
  ))
  expect_failure(module_ci_sentinel_assert_identity(import_digest, rev(samples), label = "sample"))
})

test_that("MCI-024.2 feature identity sentinel covers proteins, metabolites, and lipids", {
  catalog <- read_module_ci_fixture_packs()
  cases <- list(
    proteomics = module_ci_fixture_packs(catalog, omic = "proteomics", import_format = "dia", fixture_class = "happy_path")[[1]],
    metabolomics = module_ci_fixture_packs(catalog, omic = "metabolomics", import_format = "lc", fixture_class = "happy_path")[[1]],
    lipidomics = module_ci_fixture_packs(catalog, omic = "lipidomics", import_format = "lc", fixture_class = "happy_path")[[1]]
  )

  for (omic in names(cases)) {
    pack <- cases[[omic]]
    oracle <- module_ci_read_pack_oracle(pack)
    data <- module_ci_read_pack_table(pack)
    feature_col <- oracle$schema$feature_id_col
    extracted <- module_ci_sentinel_feature_ids(data, feature_col = feature_col)
    expected <- module_ci_as_character(oracle$feature_ids)

    expect_identical(extracted, as.character(data[[feature_col]]), info = omic)
    expect_no_error(module_ci_sentinel_assert_identity(
      module_ci_sentinel_digest(expected, kind = "feature"),
      extracted,
      label = sprintf("%s feature", omic)
    ))

    dropped <- extracted[-length(extracted)]
    expect_no_error(module_ci_sentinel_assert_identity(
      extracted,
      dropped,
      expected_dropped = extracted[[length(extracted)]],
      label = sprintf("%s feature", omic)
    ))
    expect_failure(module_ci_sentinel_assert_identity(extracted, c(dropped, "UNDECLARED_FEATURE"), label = omic))
  }
})

test_that("MCI-024.3 assay provenance sentinel spans S4 objects, session exports, DA tables, and report params", {
  lipid_object <- module_ci_lipid_da_object(layout = "pos_gcms")
  lipid_session <- list(assay_names = c("LCMS_Pos", "GCMS"))
  lipid_da_table <- data.frame(
    lipid_id = c("L1", "L2"),
    assay = c("LCMS_Pos", "GCMS"),
    comparison = "KO_vs_WT",
    stringsAsFactors = FALSE
  )
  lipid_report_params <- list(assay_names = c("LCMS_Pos", "GCMS"))

  expect_no_error(module_ci_sentinel_assert_assays(lipid_object, lipid_session))
  expect_no_error(module_ci_sentinel_assert_assays(lipid_session, lipid_da_table))
  expect_no_error(module_ci_sentinel_assert_assays(lipid_da_table, lipid_report_params))
  expect_failure(module_ci_sentinel_assert_assays(lipid_object, list(assay_names = "LCMS_Pos")))

  metab_object <- module_ci_metab_da_object(layout = "combined")
  metab_da_table <- data.frame(assay = c("LCMS_Pos", "GCMS"), metabolite_id = c("M1", "M2"))
  expect_no_error(module_ci_sentinel_assert_assays(metab_object, metab_da_table))
})

test_that("MCI-024.4 parameter fidelity sentinel compares live state and serialized payloads", {
  live_state <- list(
    config_list = list(globalParameters = list(workflow_type = "lipidomics", report_template = "lipidomics_report.rmd")),
    design = list(formula = "~ 0 + group + batch"),
    normalization = list(method = "quantile", log_offset = 2),
    ruv = list(mode = "automatic", best_k = 2L),
    itsd = list(method = "regex_global", aggregation = "sum"),
    da = list(q_cutoff = 0.05, logfc_cutoff = 1),
    enrichment = list(backend = "clusterprofiler")
  )
  serialized_state <- live_state

  paths <- c(
    "config_list.globalParameters.workflow_type",
    "config_list.globalParameters.report_template",
    "design.formula",
    "normalization.method",
    "normalization.log_offset",
    "ruv.mode",
    "ruv.best_k",
    "itsd.method",
    "da.q_cutoff",
    "enrichment.backend"
  )

  expect_no_error(module_ci_sentinel_assert_parameters(live_state, serialized_state, paths))

  drifted <- serialized_state
  drifted$da$q_cutoff <- 0.1
  expect_failure(module_ci_sentinel_assert_parameters(live_state, drifted, paths))

  missing <- serialized_state
  missing$ruv$best_k <- NULL
  expect_failure(module_ci_sentinel_assert_parameters(live_state, missing, paths))
})

test_that("MCI-024.5 shared-state sentinels detect pairwise and triple-omic collisions", {
  env <- new.env(parent = emptyenv())
  env$project_dirs <- list(
    proteomics = list(base_dir = tempfile("prot-"), source_dir = tempfile("prot-source-")),
    metabolomics = list(base_dir = tempfile("metab-"), source_dir = tempfile("metab-source-")),
    lipidomics = list(base_dir = tempfile("lipid-"), source_dir = tempfile("lipid-source-"))
  )
  env$config_list <- list(globalParameters = list(workflow_type = "proteomics"))
  env$selected_omic <- "proteomics"

  before <- module_ci_sentinel_env_snapshot(env, c("project_dirs", "config_list", "selected_omic"))
  env$selected_omic <- "metabolomics"
  after_allowed <- module_ci_sentinel_env_snapshot(env, c("project_dirs", "config_list", "selected_omic"))
  expect_no_error(module_ci_sentinel_assert_env_unchanged(before, after_allowed, allowed_changed = "selected_omic"))

  env$project_dirs$lipidomics$source_dir <- env$project_dirs$metabolomics$source_dir
  collision_snapshot <- module_ci_sentinel_env_snapshot(env, c("project_dirs", "config_list", "selected_omic"))
  expect_failure(module_ci_sentinel_assert_env_unchanged(before, collision_snapshot, allowed_changed = "selected_omic"))
  expect_failure(module_ci_sentinel_assert_project_dirs_isolated(env$project_dirs))

  env$project_dirs$lipidomics$source_dir <- tempfile("lipid-source-fixed-")
  expect_no_error(module_ci_sentinel_assert_project_dirs_isolated(env$project_dirs))
})

test_that("MCI-024.6 artifact schema sentinel validates TSV, XLSX, RDS, and report files", {
  artifact_dir <- tempfile("module-ci-sentinel-artifacts-")
  dir.create(artifact_dir, recursive = TRUE)

  tsv_path <- file.path(artifact_dir, "da_results.tsv")
  utils::write.table(
    data.frame(feature_id = c("F1", "F2"), sample = c("WT_1", "KO_1"), q_value = c(0.01, 0.2)),
    tsv_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  expect_no_error(module_ci_sentinel_assert_artifact_schema(list(
    path = tsv_path,
    type = "tsv",
    required_columns = c("feature_id", "sample", "q_value"),
    min_rows = 2L
  )))
  expect_failure(module_ci_sentinel_assert_artifact_schema(list(
    path = tsv_path,
    type = "tsv",
    required_columns = c("feature_id", "missing_column"),
    min_rows = 1L
  )))

  skip_if_not_installed("openxlsx")
  xlsx_path <- file.path(artifact_dir, "summary.xlsx")
  openxlsx::write.xlsx(data.frame(assay = "LCMS_Pos", n = 2L), xlsx_path)
  expect_no_error(module_ci_sentinel_assert_artifact_schema(list(
    path = xlsx_path,
    type = "xlsx",
    required_columns = c("assay", "n"),
    min_rows = 1L
  )))

  rds_path <- file.path(artifact_dir, "session_state.RDS")
  saveRDS(structure(list(omic_type = "lipidomics"), class = "module_ci_session_state"), rds_path)
  expect_no_error(module_ci_sentinel_assert_artifact_schema(list(
    path = rds_path,
    type = "rds",
    class = "module_ci_session_state"
  )))
  expect_failure(module_ci_sentinel_assert_artifact_schema(list(
    path = rds_path,
    type = "rds",
    class = "wrong_state_class"
  )))

  report_path <- file.path(artifact_dir, "report.html")
  writeLines("<html><body>lipidomics_report.rmd</body></html>", report_path)
  expect_no_error(module_ci_sentinel_assert_artifact_schema(list(
    path = report_path,
    type = "html",
    contains = "lipidomics_report.rmd"
  )))
  expect_failure(module_ci_sentinel_assert_artifact_schema(list(
    path = report_path,
    type = "html",
    contains = "proteomics_report.rmd"
  )))
})
