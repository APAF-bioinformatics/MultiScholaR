# E2E fixture manifest reader and validator.
# Loads tests/testdata/e2e/manifest.json and verifies all paths exist.

.E2E_REQUIRED_FIELDS <- c(
    "lane_id"
    , "omic_type"
    , "workflow_type"
    , "import_tool"
    , "use_limpa"
    , "fixture_dir"
    , "seed_file"
    , "assays"
    , "expected_contrasts"
    , "report_template"
    , "enrichment_backend"
    , "expected_exports"
    , "sample_count"
    , "group_count"
    , "groups"
    , "test_filter"
    , "module_families"
    , "touchpoints"
    , "critical_shared"
    , "report_export"
)

.E2E_VALID_MODULE_FAMILIES <- c(
    "import"
    , "design"
    , "qc"
    , "qc_peptide"
    , "qc_protein"
    , "normalization"
    , "differential_abundance"
    , "enrichment"
    , "summary_report"
)

.E2E_VALID_TOUCHPOINTS <- c(
    "browser"
    , "import"
    , "design"
    , "qc"
    , "qc_bypass"
    , "limpa_imputation"
    , "normalization"
    , "session_export"
    , "da"
    , "enrichment"
    , "summary_report"
    , "report_export"
)

.e2e_manifest_path <- function() {
    testthat::test_path("..", "testdata", "e2e", "manifest.json")
}

.e2e_fixture_root <- function() {
    testthat::test_path("..", "testdata", "e2e")
}

.validate_lane_fields <- function(lane) {
    missing_fields <- setdiff(.E2E_REQUIRED_FIELDS, names(lane))
    if (length(missing_fields) > 0) {
        rlang::abort(paste0(
            "Lane '"
            , lane[["lane_id"]] %||% "<unknown>"
            , "' is missing required fields: "
            , paste(missing_fields, collapse = ", ")
        ))
    }

    if (!is.character(lane$test_filter) || !nzchar(lane$test_filter)) {
        rlang::abort(paste0("Lane '", lane$lane_id, "' has invalid test_filter"))
    }
    tryCatch(grepl(lane$test_filter, "e2e-probe"), error = function(err) {
        rlang::abort(paste0(
            "Lane '", lane$lane_id, "' has invalid test_filter regex: ",
            conditionMessage(err)
        ))
    })

    module_families <- unlist(lane$module_families, use.names = FALSE)
    unknown_modules <- setdiff(module_families, .E2E_VALID_MODULE_FAMILIES)
    if (length(module_families) == 0L || length(unknown_modules) > 0L) {
        rlang::abort(paste0(
            "Lane '", lane$lane_id, "' has unsupported module_families: ",
            paste(unknown_modules, collapse = ", ")
        ))
    }

    touchpoints <- unlist(lane$touchpoints, use.names = FALSE)
    unknown_touchpoints <- setdiff(touchpoints, .E2E_VALID_TOUCHPOINTS)
    if (length(touchpoints) == 0L || length(unknown_touchpoints) > 0L) {
        rlang::abort(paste0(
            "Lane '", lane$lane_id, "' has unsupported touchpoints: ",
            paste(unknown_touchpoints, collapse = ", ")
        ))
    }

    if (!is.logical(lane$critical_shared) || length(lane$critical_shared) != 1L) {
        rlang::abort(paste0("Lane '", lane$lane_id, "' has invalid critical_shared"))
    }
    if (!is.logical(lane$report_export) || length(lane$report_export) != 1L) {
        rlang::abort(paste0("Lane '", lane$lane_id, "' has invalid report_export"))
    }
}

.validate_lane_paths <- function(lane, fixture_root) {
    lane_dir <- file.path(fixture_root, lane$fixture_dir)
    if (!dir.exists(lane_dir)) {
        rlang::abort(paste0(
            "Lane '"
            , lane$lane_id
            , "': fixture_dir not found: "
            , lane_dir
        ))
    }

    seed_path <- file.path(lane_dir, lane$seed_file)
    if (!file.exists(seed_path)) {
        rlang::abort(paste0(
            "Lane '"
            , lane$lane_id
            , "': seed_file not found: "
            , seed_path
        ))
    }

    if (!is.null(lane$fasta_file)) {
        fasta_path <- file.path(lane_dir, lane$fasta_file)
        if (!file.exists(fasta_path)) {
            rlang::abort(paste0(
                "Lane '"
                , lane$lane_id
                , "': fasta_file not found: "
                , fasta_path
            ))
        }
    }

    if (!is.null(lane$assay2_file)) {
        assay2_path <- file.path(lane_dir, lane$assay2_file)
        if (!file.exists(assay2_path)) {
            rlang::abort(paste0(
                "Lane '"
                , lane$lane_id
                , "': assay2_file not found: "
                , assay2_path
            ))
        }
    }
}

#' Read and validate the E2E fixture manifest
#'
#' Parses `tests/testdata/e2e/manifest.json` and validates that every lane
#' has the required fields, and that all `fixture_dir` and `seed_file` paths
#' resolve to real filesystem locations.
#'
#' @param manifest_path Path to the manifest JSON. Defaults to the standard
#'   E2E manifest location resolved via `testthat::test_path()`.
#' @param fixture_root Base directory for resolving `fixture_dir` values.
#'   Defaults to `tests/testdata/e2e/`.
#' @return Named list of lane objects, keyed by `lane_id`.
#' @export
read_e2e_manifest <- function(
    manifest_path = .e2e_manifest_path()
    , fixture_root = .e2e_fixture_root()
) {
    if (!file.exists(manifest_path)) {
        rlang::abort(paste0("E2E manifest not found: ", manifest_path))
    }

    raw <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)

    if (!is.list(raw) || is.null(raw[["lanes"]])) {
        rlang::abort(paste0(
            "Invalid manifest structure — expected top-level 'lanes' array in: "
            , manifest_path
        ))
    }

    lanes <- raw[["lanes"]]

    for (lane in lanes) {
        .validate_lane_fields(lane)
        .validate_lane_paths(lane, fixture_root)
    }

    stats::setNames(lanes, vapply(lanes, `[[`, character(1), "lane_id"))
}
