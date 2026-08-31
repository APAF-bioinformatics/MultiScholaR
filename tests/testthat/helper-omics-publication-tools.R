publicationToolPath <- function(name) {
    normalizePath(
        testthat::test_path("..", "..", "tools", "profiling", name),
        mustWork = TRUE
    )
}

publication_tool_files <- c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_lock.R",
    "omics_publication_builds.R",
    "omics_publication_repository_inputs.R",
    "omics_publication_remote_installs.R",
    "omics_publication_restore_reproducibility.R",
    "omics_publication_comparator_builds.R",
    "omics_publication_comparator_build_reproducibility.R",
    "omics_publication_comparator_evidence.R",
    "omics_publication_comparator_cleanup.R",
    "omics_publication_comparator_envelopes.R",
    "omics_publication_measure_spec.R",
    "omics_publication_linux_resources.R",
    "omics_publication_retained_resources.R",
    "omics_publication_host_safety.R",
    "omics_publication_schedule.R",
    "omics_publication_campaign_state.R",
    "omics_publication_statistics.R",
    "omics_dia_commit_repair.R"
)
publication_tool_paths <- testthat::test_path(
    "..",
    "..",
    "tools",
    "profiling",
    publication_tool_files
)
publication_tool_file <- NULL
if (all(file.exists(publication_tool_paths))) {
    for (publication_tool_file in publication_tool_files) {
        sys.source(
            publicationToolPath(publication_tool_file),
            envir = globalenv()
        )
    }
}

rm(publication_tool_file, publication_tool_files, publication_tool_paths)

publicationScheduleFixtureProjects <- function() {
    list(
        "proteomics.diann.peptide.dia.v1" = list(
            list(project_id = "project-1", source_id = "source-1"),
            list(project_id = "project-2", source_id = "source-2"),
            list(project_id = "project-3", source_id = "source-3")
        )
    )
}

publicationScheduleFixtureWorkBindings <- function() {
    subject_id <- "proteomics.diann.peptide.dia.v1"
    digest <- paste(rep("a", 64L), collapse = "")
    records <- lapply(seq_len(3L), \(index) {
        list(
            work_unit_id = "validated_input_byte",
            work_count = 1000 + index,
            exact_input_bytes = 1000 + index,
            input_sha256 = digest
        )
    })
    names(records) <- vapply(seq_len(3L), \(index) {
        paste(
            subject_id,
            paste0("project-", index),
            "estimand-1",
            sep = "::"
        )
    }, character(1))
    records
}
publicationTestSafetyEvidence <- function() {
    builder <- get(
        "publicationSafetyTraceEvidence",
        envir = globalenv(),
        inherits = FALSE
    )
    builder(
        list(list(
            safe = TRUE,
            safety_mode = "dynamic_observed_v1",
            checks = list(
                frequency_telemetry = TRUE,
                governor_stability = TRUE
            ),
            frequency_drift_fraction = NULL,
            observation = list(
                governors = list("powersave"),
                frequencies = list(
                    available = TRUE,
                    minimum_hz = 600000000,
                    median_hz = 2400000000,
                    maximum_hz = 4000000000
                ),
                load_average = list(0, 0, 0),
                maximum_thermal_celsius = 60,
                available_memory_bytes = 32 * 1024^3,
                available_disk_bytes = 512 * 1024^3
            ),
            reason = NULL,
            elapsed_seconds = 0
        )),
        timeout_seconds = 10,
        interval_seconds = 1
    )
}
