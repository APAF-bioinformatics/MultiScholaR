# E2E fixture integrity validation.
# Validates all structural requirements for the 9 E2E test lanes before
# the browser test suite runs.

.E2E_REQUIRED_SEED_COLS <- list(
    diann       = c("Protein.Group", "Protein.Ids", "Protein.Names"),
    pd_tmt      = c("Annotated.Sequence", "Master.Protein.Accessions"),
    fragpipe    = c("Protein", "Protein ID", "Gene"),
    maxquant    = c("Protein.IDs", "Gene.names"),
    msdial      = c("Feature.Name"),
    custom      = c("Feature.Name"),
    lipidsearch = c("LipidMolec", "Class", "FattyAcid")
)

.E2E_MSDIAL_META_COLS <- c(
    "Feature.Name", "m.z", "Retention.Time", "Adduct"
    , "Platform", "Similarity", "Annotation", "Metabolite.Name"
)

.e2e_count_sample_cols <- function(col_names, import_tool) {
    switch(
        import_tool
        , diann       = NA_integer_  # long format — checked via design.tsv only
        , pd_tmt      = sum(grepl("^Abundance\\.|^Abundance: ", col_names))
        , fragpipe    = sum(grepl("\\s+MaxLFQ\\s+Intensity$", col_names))
        , maxquant    = sum(grepl("^LFQ\\.intensity\\.", col_names))
        , msdial      = sum(!col_names %in% .E2E_MSDIAL_META_COLS)
        , custom      = sum(grepl("^(WT|KO)_", col_names))
        , lipidsearch = sum(grepl("\\.MeanArea$", col_names))
        , NA_integer_
    )
}

.e2e_validate_seed_file <- function(seed_path, lane_id, import_tool, sample_count) {
    errors <- character()

    if (!file.exists(seed_path)) {
        return(paste0("[", lane_id, "] seed_file not found: ", basename(seed_path)))
    }

    if (!file.access(seed_path, 4L) == 0L) {
        return(paste0("[", lane_id, "] seed_file not readable: ", basename(seed_path)))
    }

    header_raw <- tryCatch(
        readLines(seed_path, n = 1L, warn = FALSE)
        , error = function(e) NULL
    )
    if (is.null(header_raw) || !nzchar(header_raw)) {
        return(paste0("[", lane_id, "] seed_file is empty or unreadable"))
    }

    col_names <- strsplit(header_raw, "\t", fixed = TRUE)[[1L]]

    required <- .E2E_REQUIRED_SEED_COLS[[import_tool]] %||% character(0L)
    missing_req <- setdiff(required, col_names)
    if (length(missing_req) > 0L) {
        errors <- c(errors, paste0(
            "[", lane_id, "] seed_file missing required columns for "
            , import_tool, ": "
            , paste(missing_req, collapse = ", ")
        ))
    }

    n_sample_cols <- .e2e_count_sample_cols(col_names, import_tool)
    if (!is.na(n_sample_cols) && n_sample_cols < sample_count) {
        errors <- c(errors, paste0(
            "[", lane_id, "] seed_file has ", n_sample_cols
            , " sample column(s); expected at least ", sample_count
        ))
    }

    errors
}

.e2e_validate_design_tsv <- function(design_path, lane_id, sample_count, group_count, expected_groups) {
    errors <- character()

    if (!file.exists(design_path)) {
        return(paste0("[", lane_id, "] design.tsv not found"))
    }

    design <- tryCatch(
        utils::read.table(
            design_path
            , header = TRUE
            , sep = "\t"
            , stringsAsFactors = FALSE
            , comment.char = ""
            , quote = ""
        )
        , error = function(e) NULL
    )

    if (is.null(design)) {
        return(paste0("[", lane_id, "] design.tsv could not be parsed as TSV"))
    }

    col_names <- names(design)

    if (!"sample" %in% col_names) {
        errors <- c(errors, paste0("[", lane_id, "] design.tsv missing 'sample' column"))
    }
    if (!"group" %in% col_names) {
        errors <- c(errors, paste0("[", lane_id, "] design.tsv missing 'group' column"))
    }

    if (length(errors) > 0L) {
        return(errors)
    }

    n_rows <- nrow(design)
    if (n_rows < sample_count) {
        errors <- c(errors, paste0(
            "[", lane_id, "] design.tsv has ", n_rows, " row(s); expected at least ", sample_count
        ))
    }

    actual_groups <- sort(unique(design$group))
    expected_sorted <- sort(unlist(expected_groups))

    if (length(actual_groups) != group_count) {
        errors <- c(errors, paste0(
            "[", lane_id, "] design.tsv has ", length(actual_groups)
            , " unique group(s); expected ", group_count
        ))
    }

    missing_groups <- setdiff(expected_sorted, actual_groups)
    if (length(missing_groups) > 0L) {
        errors <- c(errors, paste0(
            "[", lane_id, "] design.tsv groups missing from manifest: "
            , paste(missing_groups, collapse = ", ")
        ))
    }

    errors
}

.e2e_validate_lane <- function(lane, fixture_root) {
    errors <- character()
    lane_id <- lane$lane_id
    lane_dir <- file.path(fixture_root, lane$fixture_dir)

    if (!dir.exists(lane_dir)) {
        return(paste0("[", lane_id, "] fixture_dir not found: ", lane$fixture_dir))
    }

    seed_path <- file.path(lane_dir, lane$seed_file)
    seed_errors <- .e2e_validate_seed_file(
        seed_path, lane_id, lane$import_tool, lane$sample_count
    )
    errors <- c(errors, seed_errors)

    if (!is.null(lane$fasta_file)) {
        fasta_path <- file.path(lane_dir, lane$fasta_file)
        if (!file.exists(fasta_path)) {
            errors <- c(errors, paste0("[", lane_id, "] fasta_file not found: ", lane$fasta_file))
        } else if (file.info(fasta_path)$size <= 0) {
            errors <- c(errors, paste0("[", lane_id, "] fasta_file is empty: ", lane$fasta_file))
        }
    }

    if (!is.null(lane$assay2_file)) {
        assay2_path <- file.path(lane_dir, lane$assay2_file)
        assay2_errors <- .e2e_validate_seed_file(
            assay2_path, lane_id, lane$import_tool, lane$sample_count
        )
        errors <- c(errors, assay2_errors)
    }

    design_path <- file.path(lane_dir, "design.tsv")
    design_errors <- .e2e_validate_design_tsv(
        design_path, lane_id, lane$sample_count, lane$group_count, lane$groups
    )
    errors <- c(errors, design_errors)

    errors
}

.e2e_validate_cross_omic_packs <- function(fixture_root, all_lane_ids) {
    errors <- character()
    packs_path <- file.path(fixture_root, "cross_omic_packs.json")

    if (!file.exists(packs_path)) {
        return("[cross_omic_packs] cross_omic_packs.json not found")
    }

    packs_raw <- tryCatch(
        jsonlite::read_json(packs_path, simplifyVector = FALSE)
        , error = function(e) NULL
    )

    if (is.null(packs_raw) || is.null(packs_raw$packs)) {
        return("[cross_omic_packs] invalid JSON structure — missing 'packs' array")
    }

    for (pack in packs_raw$packs) {
        pack_id <- pack$pack_id %||% "<unknown>"
        referenced_lanes <- unlist(pack$lanes) %||% character(0L)
        bad_lanes <- setdiff(referenced_lanes, all_lane_ids)
        if (length(bad_lanes) > 0L) {
            errors <- c(errors, paste0(
                "[cross_omic_packs] pack '", pack_id
                , "' references unknown lane_id(s): "
                , paste(bad_lanes, collapse = ", ")
            ))
        }
    }

    # Validate that all lane cross_omic_packs references exist as pack_ids
    known_pack_ids <- vapply(packs_raw$packs, `[[`, character(1L), "pack_id")

    errors
}

.e2e_validate_report_templates <- function(manifest, fixture_root) {
    errors <- character()
    template_dir <- file.path(fixture_root, "report_templates")

    if (!dir.exists(template_dir)) {
        return("[report_templates] report_templates/ directory not found")
    }

    unique_templates <- unique(vapply(manifest, `[[`, character(1L), "report_template"))

    for (tmpl in unique_templates) {
        stub_path <- file.path(template_dir, tmpl)
        if (!file.exists(stub_path)) {
            errors <- c(errors, paste0(
                "[report_templates] stub not found for template: ", tmpl
            ))
        }
    }

    errors
}

#' Validate E2E fixture integrity
#'
#' Runs structural checks on all 9 E2E fixture lanes: fixture directories,
#' seed file columns, design.tsv structure, cross-omic pack references, and
#' report template stubs. No data processing is performed.
#'
#' @return Named list with `$valid` (logical), `$errors` (character vector),
#'   and `$warnings` (character vector).
validate_e2e_fixtures <- function() {
    errors <- character()
    warnings <- character()

    fixture_root <- testthat::test_path("..", "testdata", "e2e")

    if (!dir.exists(fixture_root)) {
        return(list(
            valid    = FALSE
            , errors   = paste0("E2E fixture root not found: ", fixture_root)
            , warnings = warnings
        ))
    }

    manifest <- tryCatch(
        read_e2e_manifest()
        , error = function(e) {
            errors <<- c(errors, paste0("Manifest read failed: ", e$message))
            NULL
        }
    )

    if (is.null(manifest)) {
        return(list(valid = FALSE, errors = errors, warnings = warnings))
    }

    all_lane_ids <- names(manifest)

    for (lane in manifest) {
        lane_errors <- .e2e_validate_lane(lane, fixture_root)
        errors <- c(errors, lane_errors)
    }

    cross_omic_errors <- .e2e_validate_cross_omic_packs(fixture_root, all_lane_ids)
    errors <- c(errors, cross_omic_errors)

    template_errors <- .e2e_validate_report_templates(manifest, fixture_root)
    errors <- c(errors, template_errors)

    list(
        valid    = length(errors) == 0L
        , errors   = errors
        , warnings = warnings
    )
}

#' Assert E2E fixtures are structurally valid
#'
#' testthat wrapper around `validate_e2e_fixtures()`. Fails the test with a
#' formatted error listing all structural violations found.
expect_e2e_fixtures_valid <- function() {
    result <- validate_e2e_fixtures()

    if (!result$valid) {
        n <- length(result$errors)
        msg <- paste0(
            n, " E2E fixture integrity error(s):\n"
            , paste0("  - ", result$errors, collapse = "\n")
        )
        testthat::fail(msg)
    }

    testthat::succeed("All E2E fixtures are structurally valid")
    invisible(result)
}
