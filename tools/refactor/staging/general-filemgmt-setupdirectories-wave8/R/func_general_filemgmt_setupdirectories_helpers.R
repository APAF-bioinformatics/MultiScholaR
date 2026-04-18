# ----------------------------------------------------------------------------
# initializeSetupDirectoriesOmicPaths
# ----------------------------------------------------------------------------
initializeSetupDirectoriesOmicPaths <- function(
    base_dir,
    current_omic_type,
    results_path,
    results_summary_path,
    scripts_path,
    omic_config,
    timestamp
) {
    publication_graphs_dir_base <- file.path(results_path, "publication_graphs")
    # For omics types like 'integration', 'publication_graphs' might not be a
    # primary subdir, but the timestamped QC path still nests underneath it.
    qc_dir_base <- file.path(results_path, "publication_graphs", "filtering_qc")

    # Create the omic-specific data directory (for example data/proteomics).
    omic_specific_data_dir <- file.path(base_dir, "data", current_omic_type)
    dir.create(omic_specific_data_dir, recursive = TRUE, showWarnings = FALSE)
    logger::log_info("Ensured omic-specific data directory exists: {omic_specific_data_dir}")

    current_omic_paths_def <- list(
        results_base = results_path,
        results_summary_base = results_summary_path,
        data_dir = omic_specific_data_dir,
        scripts_dest_dir = scripts_path,
        scripts_source_dir = file.path(base_dir, "scripts", omic_config$scripts_source_leaf),
        publication_graphs_dir = publication_graphs_dir_base,
        qc_dir = qc_dir_base,
        time_dir = file.path(qc_dir_base, timestamp)
    )

    # Ensure the session timestamp directory exists before downstream
    # materialization references it.
    dir.create(current_omic_paths_def$time_dir, recursive = TRUE, showWarnings = FALSE)
    logger::log_info("Ensured timestamped directory exists: {current_omic_paths_def$time_dir}")

    current_omic_paths_def
}

