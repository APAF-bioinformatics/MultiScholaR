# fidelity-coverage-compare: shared
library(testthat)

downloadReportTemplate <- get(
  "downloadReportTemplate",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

write_results <- get(
  "write_results",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

writeMockWorkbook <- function(path) {
  workbook <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(workbook, "Sheet1")
  openxlsx::writeData(workbook, "Sheet1", data.frame(x = 1))
  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
}

writeMockBinary <- function(path, text = basename(path)) {
  writeBin(charToRaw(text), path)
}

localSharedGraphicsReset <- function(.local_envir = parent.frame()) {
  old_device <- getOption("device")

  while (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }

  options(device = grDevices::pdf)

  withr::defer({
    while (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
    options(device = old_device)
  }, envir = .local_envir)
}

makeSharedReportPaths <- function(root) {
  list(
    results_dir = file.path(root, "results"),
    results_summary_dir = file.path(root, "results_summary", "proteomics_Project42"),
    publication_graphs_dir = file.path(root, "publication_graphs"),
    time_dir = file.path(root, "time"),
    qc_dir = file.path(root, "qc"),
    da_output_dir = file.path(root, "da"),
    pathway_dir = file.path(root, "pathway"),
    source_dir = file.path(root, "source"),
    feature_qc_dir = file.path(root, "feature_qc"),
    subfeature_qc_dir = file.path(root, "peptide_qc")
  )
}

test_that("copyToResultsSummary validates path resolution failures", {
  copy_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        stop("mock missing project paths")
      }
    )
  )

  expect_error(
    copy_helper("proteomics", "Project42"),
    "Failed to get project paths: mock missing project paths"
  )

  expect_error(
    copy_helper("", "Project42"),
    "`omic_type` must be a single non-empty character string."
  )
})

test_that("copyToResultsSummary validates experiment labels and essential path sets", {
  partial_paths_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        list(results_dir = tempfile("results-only-"))
      }
    )
  )

  expect_error(
    partial_paths_helper("proteomics", ""),
    "`experiment_label` must be a single non-empty character string."
  )

  expect_error(
    partial_paths_helper("proteomics", "Project42"),
    "Essential paths missing from project_dirs"
  )
})

test_that("copyToResultsSummary backs up and copies proteomics report artifacts", {
  root <- tempfile("general-filemgmt-report-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  paths <- makeSharedReportPaths(root)
  dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$publication_graphs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$da_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$pathway_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$feature_qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$subfeature_qc_dir, recursive = TRUE, showWarnings = FALSE)

  dir.create(file.path(paths$publication_graphs_dir, "Interactive_Volcano_Plots"), recursive = TRUE)
  dir.create(file.path(paths$publication_graphs_dir, "NumSigDaMolecules"), recursive = TRUE)
  dir.create(file.path(paths$publication_graphs_dir, "Volcano_Plots"), recursive = TRUE)
  dir.create(file.path(paths$publication_graphs_dir, "Heatmap"), recursive = TRUE)

  writeMockBinary(file.path(paths$results_summary_dir, "existing.txt"), "already-here")
  writeMockBinary(file.path(paths$time_dir, "12_correlation_filtered_combined_plots.png"))
  writeMockBinary(file.path(paths$feature_qc_dir, "composite_QC_figure.pdf"))
  writeMockBinary(file.path(paths$feature_qc_dir, "composite_QC_figure.png"))
  writeMockBinary(file.path(paths$publication_graphs_dir, "Interactive_Volcano_Plots", "widget.html"))
  writeMockBinary(file.path(paths$publication_graphs_dir, "NumSigDaMolecules", "bars.png"))
  writeLines(
    "sig\tcount\nA\t1",
    file.path(paths$publication_graphs_dir, "NumSigDaMolecules", "study_num_sig_da_molecules.tab")
  )
  writeMockBinary(file.path(paths$publication_graphs_dir, "Volcano_Plots", "volcano.png"))
  writeMockBinary(file.path(paths$publication_graphs_dir, "Heatmap", "heatmap.html"))
  writeLines("pathway\tvalue\nset1\t0.1", file.path(paths$pathway_dir, "contrast_up_enrichment_results.tsv"))
  writeMockBinary(file.path(paths$pathway_dir, "pathway-plot.png"))
  writeLines(
    c("Workflow Name: Shared Workflow", "Timestamp: 2026-04-23 18:00:00"),
    file.path(paths$source_dir, "study_parameters.txt")
  )
  writeLines("protein\tvalue\nA\t10", file.path(paths$feature_qc_dir, "ruv_normalised_results_cln_with_replicates.tsv"))
  saveRDS(data.frame(protein = "A", value = 10), file.path(paths$feature_qc_dir, "ruv_normalised_results_cln_with_replicates.RDS"))

  writeMockWorkbook(file.path(paths$da_output_dir, "da_example_long_annot.xlsx"))

  current_rmd <- file.path(root, "current-report.Rmd")
  writeLines("---\ntitle: report\n---", current_rmd)

  contrasts_tbl <- data.frame(lhs = "A", rhs = "B", stringsAsFactors = FALSE)
  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  assign("contrasts_tbl", contrasts_tbl, envir = .GlobalEnv)
  assign("design_matrix", design_matrix, envir = .GlobalEnv)
  withr::defer(rm(list = c("contrasts_tbl", "design_matrix"), envir = .GlobalEnv))

  copy_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        paths
      }
    )
  )

  failures <- NULL
  suppressWarnings(suppressMessages(capture.output(
    failures <- copy_helper(
      omic_type = "proteomics",
      experiment_label = "Project42",
      contrasts_tbl = contrasts_tbl,
      design_matrix = design_matrix,
      force = TRUE,
      current_rmd = current_rmd
    )
  )))

  expect_length(failures, 0)
  expect_true(dir.exists(file.path(paths$results_summary_dir, "QC_figures")))
  expect_true(dir.exists(file.path(paths$results_summary_dir, "Publication_figures", "Interactive_Volcano_Plots")))
  expect_true(dir.exists(file.path(paths$results_summary_dir, "Publication_tables")))
  expect_true(dir.exists(file.path(paths$results_summary_dir, "Study_report")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "12_correlation_filtered_combined_plots.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "proteomics_composite_QC_figure.pdf")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "da_proteomics_num_sig_da_molecules.tab")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "RUV_normalised_results.tsv")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "ruv_normalised_results.RDS")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "DA_results_proteomics.xlsx")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "Pathway_enrichment_results_proteomics.xlsx")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Study_report", "contrasts_tbl.tab")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Study_report", "design_matrix.tab")))
  expect_true(file.exists(file.path(paths$source_dir, basename(current_rmd))))

  backup_dirs <- list.dirs(file.path(root, "results_summary"), full.names = TRUE, recursive = FALSE)
  backup_dirs <- backup_dirs[grepl("proteomics_Project42_backup_", basename(backup_dirs), fixed = TRUE)]
  expect_length(backup_dirs, 1)
  expect_true(file.exists(file.path(backup_dirs[[1]], "existing.txt")))
  expect_true(file.exists(file.path(backup_dirs[[1]], "backup_info.txt")))
})

test_that("copyToResultsSummary loads fallback source tables and reports copy failures", {
  root <- tempfile("general-filemgmt-fallback-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  paths <- makeSharedReportPaths(root)
  dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$publication_graphs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$da_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$pathway_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$feature_qc_dir, recursive = TRUE, showWarnings = FALSE)

  writeLines("lhs\trhs\nA\tB", file.path(paths$source_dir, "contrast_strings.tab"))
  writeLines("Sample_ID\tgroup\nS1\tA\nS2\tB", file.path(paths$source_dir, "design_matrix.tab"))
  writeLines(
    "protein\tvalue\nA\t10",
    file.path(paths$feature_qc_dir, "normalised_results_cln_with_replicates.tsv")
  )
  saveRDS(
    data.frame(protein = "A", value = 10),
    file.path(paths$feature_qc_dir, "normalised_results_cln_with_replicates.RDS")
  )

  current_rmd <- file.path(root, "fallback-report.Rmd")
  writeLines("---\ntitle: fallback\n---", current_rmd)

  old_contrasts <- if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
    get("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  had_contrasts <- exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
  old_design <- if (exists("design_matrix", envir = .GlobalEnv, inherits = FALSE)) {
    get("design_matrix", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  had_design <- exists("design_matrix", envir = .GlobalEnv, inherits = FALSE)

  if (had_contrasts) {
    rm(list = "contrasts_tbl", envir = .GlobalEnv)
  }
  if (had_design) {
    rm(list = "design_matrix", envir = .GlobalEnv)
  }

  withr::defer({
    if (had_contrasts) {
      assign("contrasts_tbl", old_contrasts, envir = .GlobalEnv)
    } else if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "contrasts_tbl", envir = .GlobalEnv)
    }

    if (had_design) {
      assign("design_matrix", old_design, envir = .GlobalEnv)
    } else if (exists("design_matrix", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "design_matrix", envir = .GlobalEnv)
    }
  })

  copy_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        paths
      }
    )
  )

  failures <- NULL
  suppressMessages(capture.output(
    failures <- copy_helper(
      omic_type = "proteomics",
      experiment_label = "Project42",
      force = TRUE,
      current_rmd = current_rmd
    )
  ))

  expect_true(file.exists(file.path(paths$results_summary_dir, "Study_report", "contrasts_tbl.tab")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Study_report", "design_matrix.tab")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "normalised_results.tsv")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "normalised_results.RDS")))
  expect_length(failures, 0)
})

test_that("copyToResultsSummary preserves cancellation and backup failure branches", {
  root <- tempfile("general-filemgmt-backup-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  base_paths <- makeSharedReportPaths(root)
  dir.create(base_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$publication_graphs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$time_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$da_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$pathway_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$feature_qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(base_paths$subfeature_qc_dir, recursive = TRUE, showWarnings = FALSE)
  writeMockBinary(file.path(base_paths$results_summary_dir, "existing.txt"), "already-here")

  cancel_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        base_paths
      },
      readline = function(prompt) "n"
    )
  )

  cancelled <- suppressMessages(suppressWarnings(capture.output(
    cancel_helper("proteomics", "Project42", force = FALSE)
  )))
  cancelled_result <- cancel_helper("proteomics", "Project42", force = FALSE)
  expect_identical(cancelled_result$status, "cancelled")
  expect_true(file.exists(file.path(base_paths$results_summary_dir, "existing.txt")))

  backup_dir_fail_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        base_paths
      },
      dir.create = function(path, recursive = FALSE, showWarnings = TRUE, mode = "0777") {
        if (grepl("_backup_", path, fixed = TRUE)) {
          FALSE
        } else {
          base::dir.create(path, recursive = recursive, showWarnings = showWarnings, mode = mode)
        }
      }
    )
  )

  backup_failures <- NULL
  expect_warning(
    suppressMessages(capture.output(
      backup_failures <- backup_dir_fail_helper("proteomics", "Project42", force = TRUE)
    )),
    "failed to copy correctly"
  )
  expect_true(any(vapply(backup_failures, `[[`, character(1), "type") == "backup_dir_creation"))

  recreate_attempts <- 0L
  recreate_fail_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        base_paths
      },
      dir.create = function(path, recursive = FALSE, showWarnings = TRUE, mode = "0777") {
        if (identical(path, base_paths$results_summary_dir)) {
          recreate_attempts <<- recreate_attempts + 1L
          return(FALSE)
        }
        base::dir.create(path, recursive = recursive, showWarnings = showWarnings, mode = mode)
      }
    )
  )

  recreate_failures <- NULL
  expect_warning(
    suppressMessages(capture.output(
      recreate_failures <- recreate_fail_helper("proteomics", "Project42", force = TRUE)
    )),
    "failed to copy correctly"
  )
  expect_true(any(vapply(recreate_failures, `[[`, character(1), "type") == "critical_dir_recreation"))

  unlink_fail_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        base_paths
      },
      unlink = function(...) {
        stop("mock unlink failure")
      }
    )
  )

  unlink_failures <- NULL
  expect_warning(
    suppressMessages(capture.output(
      unlink_failures <- unlink_fail_helper("proteomics", "Project42", force = TRUE)
    )),
    "failed to copy correctly"
  )
  expect_true(any(vapply(unlink_failures, `[[`, character(1), "type") == "dir_clear_failure"))
})

test_that("copyToResultsSummary preserves metabolomics-specific result exports", {
  root <- tempfile("general-filemgmt-metab-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  paths <- makeSharedReportPaths(root)
  dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$publication_graphs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$da_output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$pathway_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$feature_qc_dir, recursive = TRUE, showWarnings = FALSE)

  dir.create(file.path(paths$publication_graphs_dir, "Heatmaps"), recursive = TRUE, showWarnings = FALSE)
  writeMockBinary(file.path(paths$feature_qc_dir, "composite_QC_figure.png"))
  writeMockBinary(file.path(paths$feature_qc_dir, "lcms_neg_pre_norm_pca.png"))
  writeMockBinary(file.path(paths$feature_qc_dir, "lcms_pos_post_norm_density.png"))
  writeMockBinary(file.path(paths$feature_qc_dir, "batch_assay_metrics.png"))
  writeMockBinary(file.path(paths$feature_qc_dir, "itsd_normalization_plot.png"))
  writeMockBinary(file.path(paths$feature_qc_dir, "plate_correlation_plot.png"))
  saveRDS(data.frame(metab = "M1", value = 1), file.path(paths$feature_qc_dir, "normalised_results.RDS"))
  writeMockBinary(file.path(paths$publication_graphs_dir, "Heatmaps", "heatmap.png"))

  contrasts_tbl <- data.frame(lhs = "A", rhs = "B", stringsAsFactors = FALSE)
  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  copy_helper <- makeFunctionWithOverrides(
    copyToResultsSummary,
    list(
      getProjectPaths = function(...) {
        paths
      }
    )
  )

  failures <- NULL
  suppressMessages(capture.output(
    failures <- copy_helper(
      omic_type = "metabolomics",
      experiment_label = "Project42",
      contrasts_tbl = contrasts_tbl,
      design_matrix = design_matrix,
      force = TRUE
    )
  ))

  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "metabolomics_composite_QC_figure.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "lcms_neg_pre_norm_pca.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "lcms_pos_post_norm_density.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "batch_assay_metrics.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "itsd_normalization_plot.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "plate_correlation_plot.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_figures", "Heatmaps", "heatmap.png")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "normalised_results.RDS")))
  expect_length(failures, 0)
})

test_that("downloadReportTemplate reuses fresh cache and refreshes stale cache", {
  cache_root <- tempfile("report-template-cache-")
  dir.create(cache_root, recursive = TRUE)
  withr::defer(unlink(cache_root, recursive = TRUE, force = TRUE))

  cached_path <- file.path(
    cache_root, "MultiScholaR_cache", "report_templates",
    "proteomics", "report", "DIANN_report.rmd"
  )
  dir.create(dirname(cached_path), recursive = TRUE, showWarnings = FALSE)
  writeLines("fresh-cache", cached_path)

  template_helper <- makeFunctionWithOverrides(
    downloadReportTemplate,
    list(
      requireNamespace = function(...) FALSE,
      tempdir = function(...) cache_root
    )
  )

  expect_identical(
    template_helper("proteomics", "DIANN_report.rmd"),
    cached_path
  )
  expect_identical(readLines(cached_path), "fresh-cache")

  stale_path <- file.path(
    cache_root, "MultiScholaR_cache", "report_templates",
    "lipidomics", "report", "lipidomics_report.rmd"
  )
  dir.create(dirname(stale_path), recursive = TRUE, showWarnings = FALSE)
  writeLines("stale-cache", stale_path)
  Sys.setFileTime(stale_path, Sys.time() - (10 * 24 * 60 * 60))

  download_call <- NULL
  stale_helper <- makeFunctionWithOverrides(
    downloadReportTemplate,
    list(
      requireNamespace = function(...) FALSE,
      tempdir = function(...) cache_root,
      download.file = function(url, destfile, mode, quiet) {
        download_call <<- list(url = url, destfile = destfile, mode = mode, quiet = quiet)
        writeLines("downloaded-template", destfile)
        invisible(NULL)
      }
    )
  )

  downloaded_path <- stale_helper("lipidomics", "lipidomics_report.rmd")

  expect_identical(downloaded_path, stale_path)
  expect_match(download_call$url, "Workbooks/lipidomics/report/lipidomics_report.rmd", fixed = TRUE)
  expect_identical(download_call$destfile, stale_path)
  expect_identical(readLines(downloaded_path), "downloaded-template")
})

test_that("RenderReport resolves template paths and forwards render parameters", {
  root <- tempfile("render-report-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  current_paths <- list(
    base_dir = root,
    results_summary_dir = file.path(root, "results_summary", "proteomics_Project42"),
    source_dir = file.path(root, "source")
  )
  dir.create(file.path(root, "scripts", "proteomics"), recursive = TRUE, showWarnings = FALSE)
  dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(current_paths$source_dir, recursive = TRUE, showWarnings = FALSE)

  render_helper <- makeFunctionWithOverrides(
    RenderReport,
    list(
      getProjectPaths = function(...) {
        current_paths
      }
    )
  )

  expect_error(
    render_helper(
      omic_type = "proteomics",
      experiment_label = "Project42",
      rmd_filename = "DIANN_report.rmd"
    ),
    "R Markdown template file not found"
  )

  template_path <- file.path(root, "scripts", "proteomics", "DIANN_report.rmd")
  writeLines("---\ntitle: demo\n---", template_path)
  writeLines(
    c("Workflow Name: Demo Workflow", "Timestamp: 2026-04-23 21:00:00"),
    file.path(current_paths$source_dir, "study_parameters.txt")
  )

  render_ns <- asNamespace("rmarkdown")
  render_call <- NULL
  localNamespaceBinding(
    render_ns,
    "render",
    function(input, params, output_file, output_format, envir) {
      render_call <<- list(
        input = input,
        params = params,
        output_file = output_file,
        output_format = output_format,
        envir_parent = parent.env(envir)
      )
      writeLines("rendered", output_file)
      output_file
    }
  )

  output_path <- suppressMessages(
    render_helper(
      omic_type = "proteomics",
      experiment_label = "Project42",
      rmd_filename = "DIANN_report.rmd",
      output_format = "html_document"
    )
  )

  expect_identical(
    output_path,
    file.path(current_paths$results_summary_dir, "DIANN_report_proteomics_Project42.html")
  )
  expect_identical(render_call$input, template_path)
  expect_identical(render_call$output_format, "html_document")
  expect_identical(render_call$params$omic_type, "proteomics")
  expect_identical(render_call$params$experiment_label, "Project42")
  expect_identical(render_call$params$workflow_name, "Demo Workflow")
  expect_identical(render_call$params$timestamp, "2026-04-23 21:00:00")
  expect_identical(render_call$envir_parent, globalenv())
  expect_true(file.exists(output_path))
})

test_that("RenderReport validates inputs and preserves alternate output branches", {
  expect_error(
    RenderReport("", "Project42"),
    "`omic_type` must be a single non-empty character string."
  )
  expect_error(
    RenderReport("proteomics", ""),
    "`experiment_label` must be a single non-empty character string."
  )
  expect_error(
    RenderReport("proteomics", "Project42", rmd_filename = ""),
    "`rmd_filename` must be a single non-empty character string."
  )

  null_paths_helper <- makeFunctionWithOverrides(
    RenderReport,
    list(
      getProjectPaths = function(...) {
        NULL
      }
    )
  )
  expect_error(
    null_paths_helper("proteomics", "Project42"),
    "Essential paths \\(base_dir, results_summary_dir\\) missing from project_dirs"
  )

  root <- tempfile("render-report-alt-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  current_paths <- list(
    base_dir = root,
    results_summary_dir = file.path(root, "results_summary", "proteomics_Project42"),
    source_dir = file.path(root, "source")
  )
  dir.create(file.path(root, "scripts", "proteomics"), recursive = TRUE, showWarnings = FALSE)
  dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(current_paths$source_dir, recursive = TRUE, showWarnings = FALSE)
  template_path <- file.path(root, "scripts", "proteomics", "DIANN_report.rmd")
  writeLines("---\ntitle: demo\n---", template_path)

  alt_helper <- makeFunctionWithOverrides(
    RenderReport,
    list(
      getProjectPaths = function(...) {
        current_paths
      }
    )
  )

  render_calls <- list()
  localNamespaceBinding(
    asNamespace("rmarkdown"),
    "render",
    function(input, params, output_file, output_format, envir) {
      render_calls[[length(render_calls) + 1L]] <<- list(
        input = input,
        params = params,
        output_file = output_file,
        output_format = output_format,
        envir_parent = parent.env(envir)
      )
      if (grepl("\\.pdf$", output_file)) {
        writeLines("rendered", output_file)
      }
      output_file
    }
  )

  word_output <- suppressMessages(
    alt_helper(
      omic_type = "proteomics",
      experiment_label = "Project42",
      rmd_filename = "DIANN_report.rmd",
      output_format = "word_document"
    )
  )

  pdf_output <- suppressMessages(
    alt_helper(
      omic_type = "proteomics",
      experiment_label = "Project42",
      rmd_filename = "DIANN_report.rmd",
      output_format = "pdf_document"
    )
  )

  expect_identical(word_output, file.path(current_paths$results_summary_dir, "DIANN_report_proteomics_Project42.docx"))
  expect_identical(pdf_output, file.path(current_paths$results_summary_dir, "DIANN_report_proteomics_Project42.pdf"))
  expect_identical(render_calls[[1]]$params$workflow_name, "Unknown Workflow")
  expect_true(grepl("^\\d{4}-\\d{2}-\\d{2} ", render_calls[[1]]$params$timestamp))
  expect_identical(render_calls[[1]]$output_format, "word_document")
  expect_identical(render_calls[[2]]$output_format, "pdf_document")
  expect_identical(render_calls[[2]]$envir_parent, globalenv())
  expect_true(file.exists(pdf_output))
  expect_false(file.exists(word_output))
})

test_that("report helper utility functions write expected artifacts", {
  output_dir <- tempfile("report-utils-")
  dir.create(output_dir, recursive = TRUE)
  withr::defer(unlink(output_dir, recursive = TRUE, force = TRUE))
  localSharedGraphicsReset()

  plot_one <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  plot_two <- ggplot2::ggplot(mtcars, ggplot2::aes(hp, mpg)) + ggplot2::geom_point()

  pdf_calls <- new.env(parent = emptyenv())
  pdf_calls$printed <- list()
  save_list_wrapper <- makeFunctionWithOverrides(
    saveListOfPdfs,
    list(
      cairo_pdf = function(filename) {
        pdf_calls$opened <- filename
        invisible(NULL)
      },
      print = function(plot) {
        pdf_calls$printed[[length(pdf_calls$printed) + 1L]] <<- plot
        invisible(NULL)
      },
      dev.off = function() {
        pdf_calls$closed <- TRUE
        1L
      }
    )
  )

  pdf_path <- file.path(output_dir, "combined.pdf")
  expect_invisible(save_list_wrapper(list(plot_one, plot_two), pdf_path))
  expect_identical(pdf_calls$opened, pdf_path)
  expect_length(pdf_calls$printed, 2L)
  expect_true(isTRUE(pdf_calls$closed))

  simple_out <- file.path(output_dir, "simple.txt")
  source_out <- file.path(output_dir, "source.txt")
  withr::local_envvar(c(
    SHARED_RMD_SIMPLE_OUT = simple_out,
    SHARED_RMD_SOURCE_OUT = source_out
  ))

  source_simple_wrapper <- makeFunctionWithOverrides(
    sourceRmdFileSimple,
    list(
      purl = function(x, output) {
        writeLines("writeLines('simple', Sys.getenv('SHARED_RMD_SIMPLE_OUT'))", output)
        output
      }
    )
  )

  expect_invisible(source_simple_wrapper("ignored.Rmd"))
  expect_identical(readLines(simple_out), "simple")

  localNamespaceBinding(
    asNamespace("knitr"),
    "purl",
    function(file, output) {
      writeLines("writeLines('sourced', Sys.getenv('SHARED_RMD_SOURCE_OUT'))", output)
      invisible(output)
    }
  )

  old_device <- getOption("device")
  expect_invisible(sourceRmdFile("ignored.Rmd", skip_plots = TRUE))
  expect_identical(readLines(source_out), "sourced")
  expect_identical(getOption("device"), old_device)

  save_plot_calls <- new.env(parent = emptyenv())
  save_plot_calls$ggsave <- list()
  save_plot_wrapper <- makeFunctionWithOverrides(
    savePlot,
    list(
      saveRDS = function(object, file) {
        save_plot_calls$rds_paths <- c(save_plot_calls$rds_paths, file)
        base::saveRDS(list(kind = "mock"), file)
      },
      ggsave = function(filename, plot, device, width, height, ...) {
        save_plot_calls$ggsave[[length(save_plot_calls$ggsave) + 1L]] <<- list(
          filename = filename,
          device = device,
          width = width,
          height = height
        )
        writeLines("saved-plot", filename)
        invisible(NULL)
      }
    )
  )

  save_plot_wrapper(plot_one, output_dir, "single_plot", formats = c("pdf", "png"), width = 4, height = 4)
  expect_true(file.exists(file.path(output_dir, "single_plot.rds")))
  expect_true(file.exists(file.path(output_dir, "single_plot.pdf")))
  expect_true(file.exists(file.path(output_dir, "single_plot.png")))
  expect_identical(save_plot_calls$ggsave[[1]]$filename, file.path(output_dir, "single_plot.pdf"))
  expect_true(is.function(save_plot_calls$ggsave[[1]]$device))
  expect_identical(save_plot_calls$ggsave[[2]]$device, "png")

  named_plots <- list(first = plot_one, second = plot_two)
  save_plot_wrapper(named_plots, output_dir, "named_plot_list", formats = "png", width = 4, height = 4)
  expect_true(file.exists(file.path(output_dir, "named_plot_list.rds")))
  expect_true(file.exists(file.path(output_dir, "named_plot_list_first.png")))
  expect_true(file.exists(file.path(output_dir, "named_plot_list_second.png")))

  unnamed_plots <- unname(list(plot_one, plot_two))
  save_plot_wrapper(unnamed_plots, output_dir, "unnamed_plot_list", formats = "png", width = 4, height = 4)
  expect_true(file.exists(file.path(output_dir, "unnamed_plot_list_plot_1.png")))
  expect_true(file.exists(file.path(output_dir, "unnamed_plot_list_plot_2.png")))

  alias_call <- NULL
  alias_wrapper <- makeFunctionWithOverrides(
    save_plot,
    list(
      savePlot = function(plot, base_path, plot_name, formats, width, height, ...) {
        alias_call <<- list(
          plot = plot,
          base_path = base_path,
          plot_name = plot_name,
          formats = formats,
          width = width,
          height = height
        )
        invisible(NULL)
      }
    )
  )

  expect_invisible(alias_wrapper(plot_two, output_dir, "alias_plot", formats = "png", width = 4, height = 4))
  expect_identical(alias_call$base_path, output_dir)
  expect_identical(alias_call$plot_name, "alias_plot")
  expect_identical(alias_call$formats, "png")
  expect_identical(alias_call$width, 4)
  expect_identical(alias_call$height, 4)
})

test_that("write_results writes into the protein_qc directory", {
  output_dir <- tempfile("write-results-")
  dir.create(file.path(output_dir, "protein_qc"), recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(output_dir, recursive = TRUE, force = TRUE))

  writer <- makeFunctionWithOverrides(
    write_results,
    list(results_dir = output_dir)
  )

  data <- data.frame(sample = c("S1", "S2"), value = c(1, 2), stringsAsFactors = FALSE)
  writer(data, "written.tsv")

  written_path <- file.path(output_dir, "protein_qc", "written.tsv")
  expect_true(file.exists(written_path))
  written_lines <- readLines(written_path)
  expect_true(any(grepl("sample", written_lines, fixed = TRUE)))
  expect_true(any(grepl("S1", written_lines, fixed = TRUE)))
})
