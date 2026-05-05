module_ci_prot_summary_branch_config <- function(branch = "DIA",
                                                 enrichment = TRUE,
                                                 ruv = "run",
                                                 da_q_cutoff = 0.05) {
  branch <- match.arg(branch, c("DIA", "DIA_limpa", "TMT", "MaxQuant_LFQ", "FragPipe_LFQ"))
  workflow_type <- switch(branch,
    DIA = "DIA",
    DIA_limpa = "DIA",
    TMT = "TMT_PD",
    MaxQuant_LFQ = "LFQ",
    FragPipe_LFQ = "LFQ"
  )
  report_template <- switch(branch,
    DIA = "DIANN_report.rmd",
    DIA_limpa = "DIANN_limpa_report.rmd",
    TMT = "TMT_report.rmd",
    MaxQuant_LFQ = "LFQ_report.rmd",
    FragPipe_LFQ = "LFQ_report.rmd"
  )

  list(
    branch = branch,
    workflow_type = workflow_type,
    report_template = report_template,
    use_limpa = identical(branch, "DIA_limpa"),
    lfq_engine = if (branch %in% c("MaxQuant_LFQ", "FragPipe_LFQ")) branch else NULL,
    enrichment = isTRUE(enrichment),
    ruv = ruv,
    da_q_cutoff = da_q_cutoff,
    da_logfc_cutoff = 1
  )
}

module_ci_prot_summary_args <- function(config = module_ci_prot_summary_branch_config()) {
  args <- list(
    globalParameters = list(
      workflow_type = config$workflow_type,
      use_limpa = config$use_limpa,
      report_template = config$report_template,
      lfq_engine = config$lfq_engine,
      branch = config$branch
    ),
    deAnalysisParameters = list(
      q_cutoff = config$da_q_cutoff,
      logfc_cutoff = config$da_logfc_cutoff,
      contrasts = "groupB-groupA"
    ),
    normalization = list(
      ruv_mode = config$ruv,
      normalisation_method = "cyclicloess"
    )
  )

  if (isTRUE(config$use_limpa)) {
    args$proteinMissingValueImputationLimpa <- list(method = "limpa_dpc_quant")
    args$limpa_dpc_quant_results <- list(dpc_method = "limpa_dpc_quant")
  }

  if (identical(config$ruv, "run")) {
    args$ruv <- list(best_k = 2L, best_percentage = 10)
  } else {
    args$ruv <- list(ruv_skipped = TRUE, reason = "module-ci")
  }

  if (isTRUE(config$enrichment)) {
    args$enrichmentAnalysisUI <- list(
      database_source = if (identical(config$report_template, "LFQ_report.rmd")) "clusterprofiler" else "gprofiler2",
      organism_selected = if (identical(config$report_template, "LFQ_report.rmd")) "999999" else "9606",
      selected_contrast = "B_vs_A",
      q_value_cutoff = 0.05
    )
  }

  args
}

module_ci_prot_summary_object <- function(config = module_ci_prot_summary_branch_config()) {
  module_ci_prot_da_object(
    workflow_type = config$workflow_type,
    use_limpa = config$use_limpa,
    report_template = config$report_template,
    args = module_ci_prot_summary_args(config)
  )
}

module_ci_prot_summary_state_manager <- function(object,
                                                 state = "correlation_filtered",
                                                 ruv_result = list(best_k = 2L, best_percentage = 10)) {
  manager <- new.env(parent = emptyenv())
  manager$current_state <- state
  manager$state_history <- state
  manager$getHistory <- function() manager$state_history
  manager$getState <- function(state_name) {
    if (identical(state_name, state)) {
      return(object)
    }
    NULL
  }
  manager$getStateConfig <- function(state_name) {
    list(ruv_optimization_result = ruv_result)
  }
  manager
}

module_ci_prot_summary_workflow <- function(config = module_ci_prot_summary_branch_config()) {
  object <- module_ci_prot_summary_object(config)
  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- module_ci_prot_summary_state_manager(
    object,
    state = if (identical(config$ruv, "run")) "ruv_corrected" else "correlation_filtered",
    ruv_result = if (identical(config$ruv, "run")) list(best_k = 2L, best_percentage = 10) else list(ruv_skipped = TRUE)
  )
  workflow$config_list <- list(globalParameters = module_ci_prot_summary_args(config)$globalParameters)
  workflow$contrasts_tbl <- module_ci_prot_enrich_contrasts(raw = "groupB-groupA", friendly = "B_vs_A")
  workflow$design_matrix <- module_ci_prot_da_design()
  workflow$da_ui_params <- list(q_cutoff = config$da_q_cutoff, logfc_cutoff = config$da_logfc_cutoff)
  workflow$enrichment_ui_params <- if (isTRUE(config$enrichment)) {
    module_ci_prot_summary_args(config)$enrichmentAnalysisUI
  } else {
    NULL
  }
  workflow$ruv_optimization_result <- if (identical(config$ruv, "run")) {
    list(best_k = 2L, best_percentage = 10)
  } else {
    list(ruv_skipped = TRUE, reason = "module-ci")
  }
  workflow
}

module_ci_prot_summary_paths <- function() {
  root <- tempfile("module-ci-prot-summary-")
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
  list(proteomics = paths)
}

module_ci_prot_summary_write_template_set <- function(project_dirs) {
  template_dir <- file.path(project_dirs$proteomics$base_dir, "scripts", "proteomics")
  dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)
  for (template in c("DIANN_report.rmd", "DIANN_limpa_report.rmd", "TMT_report.rmd", "LFQ_report.rmd")) {
    writeLines(c("---", "title: module ci report", "---"), file.path(template_dir, template))
  }
  invisible(template_dir)
}

module_ci_prot_summary_write_publication_inputs <- function(project_dirs,
                                                            backend = "gprofiler2",
                                                            include_optional = TRUE,
                                                            ruv_applied = TRUE) {
  paths <- project_dirs$proteomics
  normalised_results <- data.frame(
    Protein.Group = c("P001", "P002"),
    Run = c("S1", "S2"),
    normalised_intensity = c(1.25, 2.5)
  )
  readr::write_tsv(module_ci_prot_da_result_table(), file.path(paths$da_output_dir, "da_groupB-groupA_long_annot.tsv"))
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(module_ci_prot_da_result_table(), file.path(paths$da_output_dir, "da_groupB-groupA_long_annot.xlsx"))
  }
  readr::write_tsv(module_ci_prot_enrich_gprofiler_table("positive"), file.path(paths$pathway_dir, "B_vs_A_up_enrichment_results.tsv"))
  readr::write_tsv(module_ci_prot_enrich_gprofiler_table("negative"), file.path(paths$pathway_dir, "B_vs_A_down_enrichment_results.tsv"))
  normalised_stem <- if (isTRUE(ruv_applied)) {
    "ruv_normalised_results_cln_with_replicates"
  } else {
    "normalised_results_cln_with_replicates"
  }
  readr::write_tsv(normalised_results, file.path(paths$feature_qc_dir, paste0(normalised_stem, ".tsv")))
  saveRDS(normalised_results, file.path(paths$feature_qc_dir, paste0(normalised_stem, ".RDS")))
  readr::write_tsv(module_ci_prot_da_design(), file.path(paths$source_dir, "design_matrix.tab"))
  readr::write_tsv(module_ci_prot_enrich_contrasts(), file.path(paths$source_dir, "contrasts_tbl.tab"))
  writeLines("study parameters", file.path(paths$source_dir, "study_parameters.txt"))

  if (isTRUE(include_optional)) {
    writeBin(charToRaw("%PDF-1.4 module-ci"), file.path(paths$feature_qc_dir, "composite_QC_figure.pdf"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$feature_qc_dir, "composite_QC_figure.png"))
    writeBin(as.raw(c(0x89, 0x50, 0x4e, 0x47)), file.path(paths$time_dir, "12_correlation_filtered_combined_plots.png"))
    for (subdir in c("Interactive_Volcano_Plots", "NumSigDaMolecules", "Volcano_Plots", "Heatmap")) {
      dir.create(file.path(paths$publication_graphs_dir, subdir), recursive = TRUE, showWarnings = FALSE)
      writeLines("artifact", file.path(paths$publication_graphs_dir, subdir, "artifact.txt"))
    }
  }

  invisible(project_dirs)
}

module_ci_prot_summary_render_stub <- function(output_dir, template) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rendered <- file.path(output_dir, sub("\\.rmd$", ".html", template, ignore.case = TRUE))
  writeLines(sprintf("<html><body>%s</body></html>", template), rendered)
  rendered
}

module_ci_prot_summary_scorecard <- function(project_dirs, session_export_path = NULL, report_path = NULL) {
  paths <- project_dirs$proteomics
  artifacts <- c(
    study_parameters = file.path(paths$source_dir, "study_parameters.txt"),
    da_workbook = file.path(paths$results_summary_dir, "Publication_tables", "DA_results_proteomics.xlsx"),
    enrichment_workbook = file.path(paths$results_summary_dir, "Publication_tables", "Pathway_enrichment_results_proteomics.xlsx"),
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
