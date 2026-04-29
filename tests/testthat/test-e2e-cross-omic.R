library(testthat)

# ---- read_cross_omic_packs() ----

test_that("read_cross_omic_packs returns 3 pack definitions", {
    packs <- read_cross_omic_packs()

    expect_type(packs, "list")
    expect_true(!is.null(packs$schema_version))
    expect_equal(length(packs$packs), 3L)
})

test_that("read_cross_omic_packs pack_ids are correct", {
    packs <- read_cross_omic_packs()
    pack_ids <- vapply(packs$packs, `[[`, character(1L), "pack_id")

    expect_setequal(pack_ids, c("prot_metab", "prot_lipid", "triple_omic"))
})

test_that("read_cross_omic_packs all packs have required fields", {
    packs <- read_cross_omic_packs()
    required <- c("pack_id", "lanes", "integration_type", "expected_factors", "expected_views")

    for (pack in packs$packs) {
        missing_fields <- setdiff(required, names(pack))
        expect_equal(
            length(missing_fields)
            , 0L
            , info = paste("Pack", pack$pack_id, "missing:", paste(missing_fields, collapse = ", "))
        )
    }
})

test_that("read_cross_omic_packs all packs reference at least 2 lanes", {
    packs <- read_cross_omic_packs()

    for (pack in packs$packs) {
        expect_gte(
            length(pack$lanes)
            , 2L
            , label = paste("Pack", pack$pack_id, "lane count")
        )
    }
})

test_that("read_cross_omic_packs integration_type is mofa for all packs", {
    packs <- read_cross_omic_packs()

    for (pack in packs$packs) {
        expect_equal(
            pack$integration_type
            , "mofa"
            , info = paste("Pack", pack$pack_id)
        )
    }
})

test_that("read_cross_omic_packs prot_metab pack is correct", {
    packs <- read_cross_omic_packs()
    pack <- purrr::detect(packs$packs, ~ .x$pack_id == "prot_metab")

    expect_false(is.null(pack))
    expect_setequal(unlist(pack$lanes), c("prot_dia", "metab_lc"))
    expect_equal(pack$expected_factors, 2L)
    expect_setequal(unlist(pack$expected_views), c("proteomics", "metabolomics"))
})

test_that("read_cross_omic_packs prot_lipid pack is correct", {
    packs <- read_cross_omic_packs()
    pack <- purrr::detect(packs$packs, ~ .x$pack_id == "prot_lipid")

    expect_false(is.null(pack))
    expect_setequal(unlist(pack$lanes), c("prot_dia", "lipid_canonical"))
    expect_equal(pack$expected_factors, 2L)
    expect_setequal(unlist(pack$expected_views), c("proteomics", "lipidomics"))
})

test_that("read_cross_omic_packs triple_omic pack is correct", {
    packs <- read_cross_omic_packs()
    pack <- purrr::detect(packs$packs, ~ .x$pack_id == "triple_omic")

    expect_false(is.null(pack))
    expect_setequal(unlist(pack$lanes), c("prot_dia", "metab_lc", "lipid_canonical"))
    expect_equal(pack$expected_factors, 3L)
    expect_setequal(unlist(pack$expected_views), c("proteomics", "metabolomics", "lipidomics"))
})

test_that("read_cross_omic_packs errors on missing file", {
    local_mocked_bindings(
        file.exists = function(path) FALSE
        , .package = "base"
    )

    expect_error(
        read_cross_omic_packs()
        , regexp = "cross_omic_packs.json not found"
    )
})

# ---- load_cross_omic_fixture() ----

test_that("load_cross_omic_fixture returns correct structure for prot_metab", {
    manifest <- read_e2e_manifest()
    fixture <- load_cross_omic_fixture("prot_metab", manifest)

    expect_type(fixture, "list")
    expect_equal(fixture$pack_id, "prot_metab")
    expect_equal(fixture$integration_type, "mofa")
    expect_equal(fixture$expected_factors, 2L)
    expect_setequal(fixture$expected_views, c("proteomics", "metabolomics"))
    expect_true(!is.null(fixture$lanes))
})

test_that("load_cross_omic_fixture lanes keyed by lane_id for prot_metab", {
    manifest <- read_e2e_manifest()
    fixture <- load_cross_omic_fixture("prot_metab", manifest)

    expect_setequal(names(fixture$lanes), c("prot_dia", "metab_lc"))
})

test_that("load_cross_omic_fixture each lane entry has fixture_path", {
    manifest <- read_e2e_manifest()
    fixture <- load_cross_omic_fixture("prot_metab", manifest)

    for (lid in names(fixture$lanes)) {
        expect_true(
            !is.null(fixture$lanes[[lid]]$fixture_path)
            , info = paste("Lane", lid, "missing fixture_path")
        )
    }
})

test_that("load_cross_omic_fixture fixture_paths point to real files", {
    manifest <- read_e2e_manifest()
    fixture <- load_cross_omic_fixture("prot_metab", manifest)

    for (lid in names(fixture$lanes)) {
        expect_true(
            file.exists(fixture$lanes[[lid]]$fixture_path)
            , info = paste("fixture_path does not exist for lane", lid)
        )
    }
})

test_that("load_cross_omic_fixture triple_omic has 3 lanes", {
    manifest <- read_e2e_manifest()
    fixture <- load_cross_omic_fixture("triple_omic", manifest)

    expect_equal(length(fixture$lanes), 3L)
    expect_setequal(names(fixture$lanes), c("prot_dia", "metab_lc", "lipid_canonical"))
    expect_setequal(fixture$expected_views, c("proteomics", "metabolomics", "lipidomics"))
})

test_that("load_cross_omic_fixture errors on unknown pack_id", {
    manifest <- read_e2e_manifest()

    expect_error(
        load_cross_omic_fixture("nonexistent_pack", manifest)
        , regexp = "not found in cross_omic_packs.json"
    )
})

test_that("load_cross_omic_fixture errors when manifest missing a referenced lane", {
    manifest <- read_e2e_manifest()
    manifest[["prot_dia"]] <- NULL  # drop required lane

    expect_error(
        load_cross_omic_fixture("prot_metab", manifest)
        , regexp = "references lane\\(s\\) not in manifest"
    )
})

test_that("load_cross_omic_fixture preserves lane metadata fields", {
    manifest <- read_e2e_manifest()
    fixture <- load_cross_omic_fixture("prot_metab", manifest)

    prot_lane <- fixture$lanes[["prot_dia"]]
    expect_equal(prot_lane$omic_type, "proteomics")
    expect_equal(prot_lane$import_tool, "diann")
    expect_equal(prot_lane$sample_count, 6L)
})

# ---- Cross-reference: manifest cross_omic_packs fields ----

test_that("manifest lanes reference only existing pack_ids", {
    manifest <- read_e2e_manifest()
    packs <- read_cross_omic_packs()
    known_pack_ids <- vapply(packs$packs, `[[`, character(1L), "pack_id")

    for (lane in manifest) {
        bad_refs <- setdiff(unlist(lane$cross_omic_packs), known_pack_ids)
        expect_equal(
            length(bad_refs)
            , 0L
            , info = paste(
                "Lane", lane$lane_id
                , "references unknown pack(s):", paste(bad_refs, collapse = ", ")
            )
        )
    }
})

test_that("prot_dia participates in all 3 packs", {
    manifest <- read_e2e_manifest()
    prot_dia_packs <- unlist(manifest[["prot_dia"]]$cross_omic_packs)

    expect_setequal(prot_dia_packs, c("prot_metab", "prot_lipid", "triple_omic"))
})

test_that("metab_lc participates in prot_metab and triple_omic", {
    manifest <- read_e2e_manifest()
    metab_lc_packs <- unlist(manifest[["metab_lc"]]$cross_omic_packs)

    expect_setequal(metab_lc_packs, c("prot_metab", "triple_omic"))
})

test_that("lipid_canonical participates in prot_lipid and triple_omic", {
    manifest <- read_e2e_manifest()
    lipid_packs <- unlist(manifest[["lipid_canonical"]]$cross_omic_packs)

    expect_setequal(lipid_packs, c("prot_lipid", "triple_omic"))
})
