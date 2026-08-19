.make_minimal_lane <- function(
    lane_id = "prot_dia"
    , fixture_dir = "prot_dia"
    , seed_file = "report.tsv"
) {
    list(
        lane_id = lane_id
        , omic_type = "proteomics"
        , workflow_type = "DIA"
        , import_tool = "diann"
        , use_limpa = FALSE
        , fixture_dir = fixture_dir
        , seed_file = seed_file
        , assays = list(NULL)
        , expected_contrasts = list("KO_vs_WT")
        , report_template = "DIANN_report.rmd"
        , enrichment_backend = "gprofiler"
        , expected_exports = list("differential_abundance.xlsx")
        , sample_count = 6L
        , group_count = 2L
        , groups = list("WT", "KO")
        , cross_omic_packs = list()
        , test_filter = "^e2e-proteomics-dia$"
        , module_families = list("import", "design", "qc", "normalization", "differential_abundance", "enrichment", "summary_report")
        , touchpoints = list("browser", "import", "design", "qc", "normalization", "da", "enrichment", "summary_report", "report_export")
        , critical_shared = TRUE
        , report_export = TRUE
    )
}

.write_tmp_manifest <- function(lanes) {
    tmp <- tempfile(fileext = ".json")
    jsonlite::write_json(list(lanes = lanes), tmp, auto_unbox = TRUE)
    tmp
}

test_that("read_e2e_manifest returns named list of 9 lanes", {
    manifest <- read_e2e_manifest()

    expect_type(manifest, "list")
    expect_equal(length(manifest), 9L)
    expect_setequal(names(manifest), c(
        "prot_dia"
        , "prot_dia_limpa"
        , "prot_tmt"
        , "prot_lfq"
        , "prot_lfq_fragpipe"
        , "metab_lc"
        , "metab_gc"
        , "metab_combined"
        , "lipid_canonical"
    ))
    for (lane_id in names(manifest)) {
        expect_equal(manifest[[lane_id]][["lane_id"]], lane_id)
    }
})

test_that("read_e2e_manifest lanes have all required fields", {
    required <- c(
        "lane_id", "omic_type", "workflow_type", "import_tool", "use_limpa"
        , "fixture_dir", "seed_file", "assays", "expected_contrasts"
        , "report_template", "enrichment_backend", "expected_exports"
        , "sample_count", "group_count", "groups"
        , "test_filter", "module_families", "touchpoints"
        , "critical_shared", "report_export"
    )
    manifest <- read_e2e_manifest()

    for (lane in manifest) {
        missing_fields <- setdiff(required, names(lane))
        expect_equal(
            length(missing_fields)
            , 0L
            , info = paste("Lane", lane[["lane_id"]], "missing:", paste(missing_fields, collapse = ", "))
        )
    }
})

test_that("read_e2e_manifest lanes expose impact routing metadata", {
    manifest <- read_e2e_manifest()

    expected_filters <- c(
        prot_dia = "^e2e-proteomics-dia$",
        prot_dia_limpa = "^e2e-proteomics-dia-limpa$",
        prot_tmt = "^e2e-proteomics-tmt-lfq$",
        prot_lfq = "^e2e-proteomics-tmt-lfq$",
        prot_lfq_fragpipe = "^e2e-proteomics-tmt-lfq$",
        metab_lc = "^e2e-metabolomics-lc-gc$",
        metab_gc = "^e2e-metabolomics-lc-gc$",
        metab_combined = "^e2e-metabolomics-lc-gc$",
        lipid_canonical = "^e2e-lipidomics-canonical$"
    )

    expect_setequal(names(manifest), names(expected_filters))
    for (lane_id in names(expected_filters)) {
        lane <- manifest[[lane_id]]
        expect_identical(lane$test_filter, unname(expected_filters[[lane_id]]))
        expect_true(length(lane$module_families) > 0L, info = lane_id)
        expect_true("import" %in% unlist(lane$module_families), info = lane_id)
        expect_true("report_export" %in% unlist(lane$touchpoints), info = lane_id)
        expect_type(lane$critical_shared, "logical")
        expect_type(lane$report_export, "logical")
    }
})

test_that("read_e2e_manifest validates optional proteomics FASTA fixtures", {
    manifest <- read_e2e_manifest()
    proteomics_lanes <- manifest[vapply(
        manifest,
        \(lane) identical(lane$omic_type, "proteomics"),
        logical(1L)
    )]

    expect_gt(length(proteomics_lanes), 0L)
    for (lane in proteomics_lanes) {
        expect_false(is.null(lane$fasta_file), info = lane$lane_id)
        expect_true(
            file.exists(file.path(.e2e_fixture_root(), lane$fixture_dir, lane$fasta_file)),
            info = lane$lane_id
        )
    }
})

test_that("read_e2e_manifest errors on missing fixture_dir", {
    fixture_root <- testthat::test_path("..", "testdata", "e2e")
    lane <- .make_minimal_lane(fixture_dir = "nonexistent_dir_xyz")
    tmp <- .write_tmp_manifest(list(lane))
    on.exit(unlink(tmp))

    expect_error(
        read_e2e_manifest(manifest_path = tmp, fixture_root = fixture_root)
        , regexp = "fixture_dir not found"
    )
})

test_that("read_e2e_manifest errors on missing optional fasta_file", {
    fixture_root <- testthat::test_path("..", "testdata", "e2e")
    lane <- .make_minimal_lane()
    lane$fasta_file <- "nonexistent_fasta_xyz.fasta"
    tmp <- .write_tmp_manifest(list(lane))
    on.exit(unlink(tmp))

    expect_error(
        read_e2e_manifest(manifest_path = tmp, fixture_root = fixture_root)
        , regexp = "fasta_file not found"
    )
})

test_that("read_e2e_manifest errors on missing seed_file", {
    fixture_root <- testthat::test_path("..", "testdata", "e2e")
    lane <- .make_minimal_lane(seed_file = "nonexistent_seed_xyz.tsv")
    tmp <- .write_tmp_manifest(list(lane))
    on.exit(unlink(tmp))

    expect_error(
        read_e2e_manifest(manifest_path = tmp, fixture_root = fixture_root)
        , regexp = "seed_file not found"
    )
})

test_that("read_e2e_manifest errors on missing required field", {
    fixture_root <- testthat::test_path("..", "testdata", "e2e")
    lane <- .make_minimal_lane()
    lane[["workflow_type"]] <- NULL  # drop a required field
    tmp <- .write_tmp_manifest(list(lane))
    on.exit(unlink(tmp))

    expect_error(
        read_e2e_manifest(manifest_path = tmp, fixture_root = fixture_root)
        , regexp = "missing required fields"
    )
})

test_that("all e2e fixture directories exist", {
    fixture_root <- testthat::test_path("..", "testdata", "e2e")
    expected_dirs <- c(
        "prot_dia", "prot_dia_limpa", "prot_tmt", "prot_lfq"
        , "metab_lc", "metab_gc", "metab_combined", "lipid_canonical"
    )
    for (d in expected_dirs) {
        expect_true(dir.exists(file.path(fixture_root, d)), info = paste("Missing dir:", d))
    }
})

test_that("all 9 e2e seed files exist", {
    manifest <- read_e2e_manifest()
    fixture_root <- testthat::test_path("..", "testdata", "e2e")

    for (lane in manifest) {
        seed_path <- file.path(fixture_root, lane$fixture_dir, lane$seed_file)
        expect_true(file.exists(seed_path), info = paste("Missing seed:", seed_path))
    }
})

test_that("DIA limpa lane uses the intended MNAR fixture", {
    lane <- read_e2e_manifest()[["prot_dia_limpa"]]
    seed_path <- file.path(.e2e_fixture_root(), lane$fixture_dir, lane$seed_file)
    fixture <- utils::read.delim(seed_path, check.names = FALSE)

    expect_identical(lane$seed_file, "seed_report.tsv")
    expect_equal(nrow(fixture), 72L)
    expect_setequal(unique(fixture$Protein.Group), sprintf("P%05d", 1:6))
    expect_equal(sum(is.na(fixture$Precursor.Normalised)), 20L)
})
