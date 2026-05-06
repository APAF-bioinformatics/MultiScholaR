module_ci_lipid_summary_config <- function(layout = c("pos_neg", "pos_gcms", "single_pos"),
                                           normalization_method = "quantile",
                                           ruv_mode = c("automatic", "manual", "skip"),
                                           itsd_method = "regex_global",
                                           da_q_cutoff = 0.05) {
  layout <- match.arg(layout)
  ruv_mode <- match.arg(ruv_mode)
  list(
    layout = layout,
    workflow_type = "lipidomics",
    report_template = "lipidomics_report.rmd",
    normalization_method = normalization_method,
    ruv_mode = ruv_mode,
    itsd_method = itsd_method,
    da_q_cutoff = da_q_cutoff,
    da_logfc_cutoff = 1,
    formula_string = if (identical(layout, "single_pos")) "~ 0 + group" else "~ 0 + group + batch"
  )
}

module_ci_lipid_summary_assay_names <- function(layout) {
  names(module_ci_lipid_da_assays(layout = layout))
}

module_ci_lipid_summary_args <- function(config = module_ci_lipid_summary_config()) {
  assay_names <- module_ci_lipid_summary_assay_names(config$layout)
  ruv_args <- if (identical(config$ruv_mode, "skip")) {
    list(
      ruv_skipped = TRUE,
      ruv_mode = "skip"
    )
  } else {
    list(
      ruv_grouping_variable = "group",
      ruv_mode = config$ruv_mode,
      ruv_number_k = stats::setNames(as.list(rep(2L, length(assay_names))), assay_names),
      ctrl = stats::setNames(lapply(assay_names, function(...) c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE)), assay_names)
    )
  }

  c(
    list(
      globalParameters = list(
        workflow_type = config$workflow_type,
        report_template = config$report_template,
        assay_layout = config$layout
      ),
      ITSDNormalization = list(
        applied = !identical(config$itsd_method, "skip"),
        method_type = config$itsd_method,
        itsd_aggregation = "sum",
        itsd_pattern_columns = c("lipid", "LipidClass"),
        removed_itsd = TRUE,
        itsd_features_per_assay = stats::setNames(lapply(assay_names, function(assay) paste0("IS_", assay)), assay_names),
        itsd_counts_per_assay = stats::setNames(as.list(rep(1L, length(assay_names))), assay_names),
        timestamp = as.POSIXct("2026-05-06 10:00:00", tz = "UTC")
      ),
      log_transformed = TRUE,
      log_transform_offset = 2,
      normalisation_method = config$normalization_method,
      daAnalysisParameters = list(
        formula_string = config$formula_string,
        q_cutoff = config$da_q_cutoff,
        logfc_cutoff = config$da_logfc_cutoff,
        contrasts = module_ci_lipid_da_contrasts("multi_group")$contrasts,
        friendly_names = module_ci_lipid_da_contrasts("multi_group")$friendly_names
      )
    ),
    ruv_args
  )
}

module_ci_lipid_summary_object <- function(config = module_ci_lipid_summary_config()) {
  object <- module_ci_lipid_da_object(
    layout = config$layout,
    formula_string = config$formula_string
  )
  object@args <- module_ci_lipid_summary_args(config)
  object
}

module_ci_lipid_summary_state_manager <- function(object,
                                                  state = "lipid_correlation_filtered",
                                                  state_config = list(ruv_optimization_result = list(best_k = 2L, best_percentage = 40))) {
  manager <- new.env(parent = emptyenv())
  manager$current_state <- state
  manager$state_history <- unique(c("lipid_qc_complete", "lipid_normalized", state))
  manager$getHistory <- function() manager$state_history
  manager$getState <- function(state_name = manager$current_state) {
    if (state_name %in% manager$state_history) {
      return(object)
    }
    NULL
  }
  manager$getStateConfig <- function(state_name) state_config
  manager
}

module_ci_lipid_summary_workflow <- function(config = module_ci_lipid_summary_config()) {
  object <- module_ci_lipid_summary_object(config)
  ruv_result <- if (identical(config$ruv_mode, "skip")) {
    list(ruv_skipped = TRUE, reason = "module-ci")
  } else {
    list(mode = config$ruv_mode, best_k = 2L, best_percentage = 40)
  }
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- module_ci_lipid_summary_state_manager(
    object,
    state = if (identical(config$ruv_mode, "skip")) "lipid_normalized" else "lipid_correlation_filtered",
    state_config = list(ruv_optimization_result = ruv_result)
  )
  workflow$config_list <- list(globalParameters = module_ci_lipid_summary_args(config)$globalParameters)
  workflow$contrasts_tbl <- module_ci_lipid_da_contrasts("multi_group")
  workflow$design_matrix <- object@design_matrix
  workflow$da_ui_params <- list(
    q_cutoff = config$da_q_cutoff,
    logfc_cutoff = config$da_logfc_cutoff,
    formula_string = config$formula_string
  )
  workflow$normalization_ui_params <- list(
    normalisation_method = config$normalization_method,
    log_offset = 2
  )
  workflow$itsd_ui_params <- module_ci_lipid_summary_args(config)$ITSDNormalization
  workflow$ruv_optimization_result <- ruv_result
  workflow
}

module_ci_lipid_summary_paths <- function() {
  root <- tempfile("module-ci-lipid-summary-")
  paths <- list(
    base_dir = root,
    results_dir = file.path(root, "results"),
    results_summary_dir = file.path(root, "results_summary"),
    publication_graphs_dir = file.path(root, "publication_graphs"),
    time_dir = file.path(root, "time"),
    qc_dir = file.path(root, "qc"),
    da_output_dir = file.path(root, "da_output"),
    pathway_dir = file.path(root, "pathway_enrichment"),
    source_dir = file.path(root, "source"),
    feature_qc_dir = file.path(root, "feature_qc"),
    integration_dir = file.path(root, "integration")
  )
  lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)
  list(lipidomics = paths)
}

module_ci_lipid_summary_write_template <- function(project_dirs) {
  template_dir <- file.path(project_dirs$lipidomics$base_dir, "scripts", "lipidomics")
  dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("---", "title: module ci lipidomics report", "---"), file.path(template_dir, "lipidomics_report.rmd"))
  invisible(template_dir)
}

module_ci_lipid_summary_write_publication_inputs <- function(project_dirs,
                                                             config = module_ci_lipid_summary_config(),
                                                             include_optional = TRUE,
                                                             ruv_applied = !identical(config$ruv_mode, "skip")) {
  paths <- project_dirs$lipidomics
  object <- module_ci_lipid_summary_object(config)
  da_long <- module_ci_lipid_da_results_long("mixed")
  da_long <- da_long[da_long$assay %in% names(object@lipid_data), , drop = FALSE]
  normalised_results <- data.frame(
    lipid_id = module_ci_lipid_da_feature_ids()[seq_len(3)],
    assay = rep(names(object@lipid_data), length.out = 3),
    Run = module_ci_lipid_da_samples()[seq_len(3)],
    normalised_intensity = c(10.5, 20.5, 30.5),
    stringsAsFactors = FALSE
  )

  for (assay_name in unique(da_long$assay)) {
    assay_mode <- if (grepl("pos", assay_name, ignore.case = TRUE)) {
      "posmode"
    } else if (grepl("neg", assay_name, ignore.case = TRUE)) {
      "negmode"
    } else {
      gsub("[^A-Za-z0-9]", "", tolower(assay_name))
    }
    assay_rows <- da_long[da_long$assay == assay_name & da_long$comparison == "groupKO-groupWT", , drop = FALSE]
    file_base <- sprintf("de_%s_lipids_groupKO-groupWT_long_annot", assay_mode)
    readr::write_tsv(assay_rows, file.path(paths$da_output_dir, paste0(file_base, ".tsv")))
    openxlsx::write.xlsx(assay_rows, file.path(paths$da_output_dir, paste0(file_base, ".xlsx")))
  }

  normalized_stem <- if (isTRUE(ruv_applied)) "ruv_normalised_results" else "normalised_results"
  readr::write_tsv(normalised_results, file.path(paths$feature_qc_dir, paste0(normalized_stem, ".tsv")))
  saveRDS(normalised_results, file.path(paths$feature_qc_dir, paste0(normalized_stem, ".RDS")))
  readr::write_tsv(object@design_matrix, file.path(paths$source_dir, "design_matrix.tab"))
  readr::write_tsv(module_ci_lipid_da_contrasts("multi_group"), file.path(paths$source_dir, "contrasts_tbl.tab"))
  writeLines("study parameters", file.path(paths$source_dir, "study_parameters.txt"))

  if (isTRUE(include_optional)) {
    writeBin(charToRaw("%PDF-1.4 module-ci"), file.path(paths$feature_qc_dir, "composite_QC_figure.pdf"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$feature_qc_dir, "composite_QC_figure.png"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$feature_qc_dir, "lcms_pos_pre_norm_pca.png"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$feature_qc_dir, "lcms_neg_ruv_corrected_rle.png"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$feature_qc_dir, "itsd_lcms_pos.png"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$time_dir, "12_correlation_filtered_combined_plots.png"))

    for (subdir in c(
      "Interactive_Volcano_Plots",
      "NumSigDaMolecules",
      "Volcano_Plots",
      "Heatmap",
      "Heatmaps",
      "NumSigDeMolecules"
    )) {
      dir.create(file.path(paths$publication_graphs_dir, subdir), recursive = TRUE, showWarnings = FALSE)
      writeLines("artifact", file.path(paths$publication_graphs_dir, subdir, "artifact.txt"))
    }
    writeLines(
      c("mode\tcomparison\tup\tdown", "posmode\tgroupKO-groupWT\t2\t1"),
      file.path(paths$publication_graphs_dir, "NumSigDaMolecules", "lipids_num_sig_da_molecules.tab")
    )
    writeLines(
      c("mode\tcomparison\tup\tdown", "posmode\tgroupKO-groupWT\t2\t1"),
      file.path(paths$publication_graphs_dir, "NumSigDeMolecules", "lipids_num_sig_de_molecules.tab")
    )
  }

  invisible(project_dirs)
}

module_ci_lipid_summary_render_stub <- function(output_dir, template = "lipidomics_report.rmd") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rendered <- file.path(output_dir, sub("\\.rmd$", ".html", template, ignore.case = TRUE))
  writeLines(sprintf("<html><body>%s</body></html>", template), rendered)
  rendered
}

module_ci_lipid_summary_scorecard <- function(project_dirs, session_export_path = NULL, report_path = NULL) {
  paths <- project_dirs$lipidomics
  artifacts <- c(
    study_parameters = file.path(paths$source_dir, "study_parameters.txt"),
    da_workbook = file.path(paths$results_summary_dir, "Publication_tables", "DA_results_lipidomics.xlsx"),
    normalized_tsv = file.path(paths$results_summary_dir, "Publication_tables", "RUV_normalised_results.tsv"),
    normalized_rds = file.path(paths$results_summary_dir, "Publication_tables", "ruv_normalised_results.RDS"),
    qc_figure = file.path(paths$results_summary_dir, "QC_figures", "lcms_pos_pre_norm_pca.png"),
    num_sig = file.path(paths$results_summary_dir, "Publication_tables", "da_lipidomics_num_sig_de_molecules.tab"),
    session_state = session_export_path,
    report_file = report_path
  )
  data.frame(
    artifact = names(artifacts),
    path = unname(artifacts),
    exists = vapply(artifacts, function(path) !is.null(path) && !is.na(path) && file.exists(path), logical(1)),
    stringsAsFactors = FALSE
  )
}

module_ci_lipid_summary_with_project_dirs <- function(project_dirs, code) {
  had_project_dirs <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
  old_project_dirs <- if (had_project_dirs) get("project_dirs", envir = .GlobalEnv, inherits = FALSE) else NULL
  assign("project_dirs", project_dirs, envir = .GlobalEnv)
  on.exit({
    if (had_project_dirs) {
      assign("project_dirs", old_project_dirs, envir = .GlobalEnv)
    } else if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
      rm("project_dirs", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(code)
}
