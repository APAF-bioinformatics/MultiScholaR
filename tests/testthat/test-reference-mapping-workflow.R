# KIT-167 section 10 and 7: identity, store round trip, replay offline, compute-again, record sheet, Methods text,
# a curated mapping without network. The engine is the stub of helper-reference-mapping.R.

k167Map <- function(strain, ...) {
    stub <- k167StubDiamond()
    k167UseStubEngine(stub)
    mapStrainToReference(strain$path, "UP000000750", ...)
}

test_that("a curated mapping runs end to end with no network call and pairs each protein with its source", {
    k167NoNetwork()
    strain <- k167StrainFromReference()
    result <- k167Map(strain)
    expect_identical(nrow(result$mapping), 60L)
    expect_identical(result$counts$by_status$one_to_one, 60L)
    expect_identical(result$mapping$reference_accession_used, strain$reference_accessions)
    expect_true(all(grepl("^SPy_[0-9]{4}", result$mapping$reference_locus[!is.na(result$mapping$reference_locus)])))
    expect_true(all(is.na(result$mapping$reference_string_id)))
    expect_identical(result$provenance$reference_upid, "UP000000750")
    expect_identical(result$provenance$reference_release, "2026_04")
    expect_equal(result$identity_summary$median, 99.5)
    expect_false(result$identity_summary$caution)
})

test_that("the identity changes with each setting and not with the thread count", {
    strain <- k167StrainFromReference()
    base <- k167Map(strain)$provenance
    sha <- function(p) referenceMappingIdentity(p)$sha256
    expect_identical(sha(base), sha(k167Map(strain, threads = 2L)$provenance))
    for (changed in list(list(coverage = 70), list(evalue = 1e-6), list(min_identity = 40),
        list(header_classes = c("refseq", "plain")))) {
        expect_false(identical(sha(base), sha(do.call(k167Map, c(list(strain), changed))$provenance)),
            info = names(changed))
    }
    other_engine <- base
    other_engine$engine_sha256 <- strrep("b", 64)
    expect_false(identical(sha(base), sha(other_engine)))
})

test_that("min_identity NULL gives the table without a floor; a floor changes the table and the identity", {
    strain <- k167StrainFromReference()
    none <- k167Map(strain)
    expect_identical(none$mapping, k167Map(strain, min_identity = NULL)$mapping)
    floored <- k167Map(strain, min_identity = 99.9)
    expect_identical(floored$counts$by_status$one_to_one, 0L)
    expect_false(identical(none$provenance$mapping_sha256, floored$provenance$mapping_sha256))
})

test_that("the store round trip, replay offline and the reuse rule", {
    manager <- k167Owner()
    strain <- k167StrainFromReference()
    result <- k167Map(strain)
    storeReferenceMapping(manager, result, inputs = c(reference_mapping_strain_fasta = strain$path))
    stored <- readReferenceMapping(manager)
    expect_identical(canonicalReferenceMappingText(stored$mapping), canonicalReferenceMappingText(result$mapping))
    expect_identical(manager$.context$snapshot$file_dependencies$reference_mapping_strain_fasta$sha256,
        digest::digest(file = strain$path, algo = "sha256"))
    withr::local_envvar(MSR_DIAMOND = file.path(tempdir(), "no-diamond-here"))
    replayed <- replayReferenceMapping(manager)
    expect_identical(canonicalReferenceMappingSha256(replayed), result$provenance$mapping_sha256)
    expect_true(referenceMappingReusable(stored, referenceMappingIdentity(result$provenance)$sha256))
    expect_false(referenceMappingReusable(stored, strrep("c", 64)))
})

test_that("compute-again is refused in words without the engine and equals the stored table with it", {
    manager <- k167Owner()
    strain <- k167StrainFromReference()
    storeReferenceMapping(manager, k167Map(strain))
    withr::local_envvar(MSR_DIAMOND = file.path(tempdir(), "no-diamond-here"))
    withr::local_options(MultiScholaR.diamond_path = NULL)
    expect_error(computeReferenceMappingAgain(manager, strain$path, "UP000000750"),
        class = "reference_mapping_engine_missing")
    k167UseStubEngine(k167StubDiamond())
    again <- computeReferenceMappingAgain(manager, strain$path, "UP000000750")
    expect_true(again$identical)
})

test_that("the record sheet carries what A14 lists and the Methods paragraph states the applicability", {
    sheet <- referenceMappingRecordSheet(k167Map(k167StrainFromReference()))
    for (field in c("reference", "uniprot_upid", "uniprot_release", "strain_fasta_sha256", "reference_fasta_sha256",
        "reference_xref_sha256", "string_species", "string_status", "kegg_code", "engine", "engine_version",
        "engine_sha256", "coverage", "evalue", "min_identity", "count_one_to_one", "share_one_way", "identity_median",
        "identity_p05", "identity_share_below_40", "applicability", "mapping_sha256", "identity_sha256")) {
        expect_true(field %in% sheet$field, info = field)
    }
    expect_identical(sheet$value[sheet$field == "string_status"], "unavailable")
    methods <- referenceMappingMethodsText(sheet)
    expect_match(methods, referenceMappingTexts$applicability, fixed = TRUE)
    expect_match(methods, "UP000000750", fixed = TRUE)
    expect_match(methods, "2026_04", fixed = TRUE)
    expect_match(methods, "60 of 60 strain proteins", fixed = TRUE)
})

test_that("the gate is off by default", {
    withr::local_options(MultiScholaR.reference_mapping = NULL)
    expect_false(referenceMappingEnabled())
})

test_that("the transcriptomics features are unchanged without a gene table and gain the gene columns with one", {
    features <- data.frame(feature_id = c("T_0001", "T_0002", "X_9"), display_name = c("thrL", "T_0002", "X_9"),
        stringsAsFactors = FALSE)
    expect_identical(trxReferenceMappingFeatures(features, NULL), features)
    genes <- data.frame(gene_id = c("T_0001", "T_0002"), locus_tag = c("T_0001", "T_0002"),
        reference_accession_used = c("P1", NA), map_status = c("one_to_one", "no_protein"), gene_length = c(300, 401),
        stringsAsFactors = FALSE)
    with <- trxReferenceMappingFeatures(features, genes)
    expect_identical(with[names(features)], features)
    expect_identical(with$reference_accession_used, c("P1", NA, NA))
    expect_identical(with$map_status, c("one_to_one", "no_protein", NA))
    expect_identical(with$gene_length, c(300, 401, NA))
})
