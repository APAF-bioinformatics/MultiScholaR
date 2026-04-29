# E2E helpers for cross-omic fixture packs used in multi-omic integration tests.
# Packs bundle multiple single-lane fixtures by reference — no data duplication.

.e2e_cross_omic_packs_path <- function() {
    testthat::test_path("..", "testdata", "e2e", "cross_omic_packs.json")
}

.e2e_fixture_base_dir <- function() {
    testthat::test_path("..", "testdata", "e2e")
}

#' Read and validate cross-omic pack definitions
#'
#' Reads `tests/testdata/e2e/cross_omic_packs.json` and validates that all
#' required fields are present for each pack.
#'
#' @return A list with `schema_version` (character) and `packs` (list of pack
#'   definitions, each with pack_id, lanes, integration_type, expected_factors,
#'   expected_views).
read_cross_omic_packs <- function() {
    packs_path <- .e2e_cross_omic_packs_path()

    if (!file.exists(packs_path)) {
        rlang::abort(paste0("cross_omic_packs.json not found: ", packs_path))
    }

    packs <- jsonlite::fromJSON(packs_path, simplifyVector = FALSE)

    required_top <- c("schema_version", "packs")
    missing_top <- setdiff(required_top, names(packs))
    if (length(missing_top) > 0) {
        rlang::abort(paste0(
            "cross_omic_packs.json missing required fields: "
            , paste(missing_top, collapse = ", ")
        ))
    }

    pack_required <- c("pack_id", "lanes", "integration_type", "expected_factors", "expected_views")
    purrr::walk(packs$packs, function(pack) {
        missing_fields <- setdiff(pack_required, names(pack))
        if (length(missing_fields) > 0) {
            rlang::abort(paste0(
                "Pack '", pack$pack_id, "' missing required fields: "
                , paste(missing_fields, collapse = ", ")
            ))
        }
        if (length(pack$lanes) < 2) {
            rlang::abort(paste0(
                "Pack '", pack$pack_id, "' must reference at least 2 lanes, got "
                , length(pack$lanes)
            ))
        }
    })

    packs
}

#' Load a combined cross-omic fixture from individual lane fixtures
#'
#' Given a pack ID, looks up all referenced lanes in the manifest and assembles
#' a combined fixture structure suitable for multi-omic module initialization.
#' All referenced lanes must exist in the manifest; fixture file paths are
#' resolved relative to `tests/testdata/e2e/`.
#'
#' No data duplication: each lane's fixture is referenced by its existing seed
#' file path rather than copied.
#'
#' @param pack_id Character. Pack ID from `cross_omic_packs.json`.
#' @param manifest Named list of lane definitions as returned by
#'   `read_e2e_manifest()` (i.e. a list keyed by lane_id).
#'
#' @return A named list with:
#'   - `pack_id` (character)
#'   - `integration_type` (character, e.g. "mofa")
#'   - `expected_factors` (integer)
#'   - `expected_views` (character vector)
#'   - `lanes` (named list keyed by lane_id; each entry is the lane's manifest
#'     fields plus `fixture_path` pointing to the seed file)
load_cross_omic_fixture <- function(pack_id, manifest) {
    packs <- read_cross_omic_packs()

    pack <- purrr::detect(packs$packs, ~ .x$pack_id == pack_id)
    if (is.null(pack)) {
        rlang::abort(paste0(
            "Pack '", pack_id, "' not found in cross_omic_packs.json. "
            , "Available: "
            , paste(purrr::map_chr(packs$packs, "pack_id"), collapse = ", ")
        ))
    }

    # Accept named list from read_e2e_manifest() — already keyed by lane_id
    lane_index <- manifest

    missing_lanes <- setdiff(unlist(pack$lanes), names(lane_index))
    if (length(missing_lanes) > 0) {
        rlang::abort(paste0(
            "Pack '", pack_id, "' references lane(s) not in manifest: "
            , paste(missing_lanes, collapse = ", ")
        ))
    }

    base_dir <- .e2e_fixture_base_dir()

    lane_fixtures <- purrr::set_names(
        purrr::map(unlist(pack$lanes), function(lid) {
            lane <- lane_index[[lid]]
            fixture_path <- file.path(base_dir, lane$fixture_dir, lane$seed_file)
            c(lane, list(fixture_path = fixture_path))
        })
        , unlist(pack$lanes)
    )

    list(
        pack_id = pack$pack_id
        , integration_type = pack$integration_type
        , expected_factors = pack$expected_factors
        , expected_views = unlist(pack$expected_views)
        , lanes = lane_fixtures
    )
}
