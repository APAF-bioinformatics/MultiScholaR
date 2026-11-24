# KIT-167 section 10 and 5.7: the shipped reference list and the multi-valued cross-reference parsing, on rows copied
# from the captured cross-reference tables (UniProt release 2026_04, CC BY 4.0).

k167Xref <- function(...) {
    rows <- list(...)
    data.frame(Entry = vapply(rows, `[[`, "", 1), `Gene Names (primary)` = vapply(rows, `[[`, "", 2),
        `Gene Names (ordered locus)` = vapply(rows, `[[`, "", 3), STRING = vapply(rows, `[[`, "", 4),
        KEGG = vapply(rows, `[[`, "", 5), check.names = FALSE, stringsAsFactors = FALSE)
}

test_that("the shipped list names the four references with their service identifiers", {
    organisms <- readReferenceOrganisms()
    expect_identical(organisms$uniprot_upid, c("UP000000625", "UP000008816", "UP000000265", "UP000000750"))
    expect_identical(organisms$uniprot_taxon, c(83333L, 93061L, 272620L, 301447L))
    expect_identical(organisms$kegg_code, c("eco", "sao", "kpn", "spy"))
    expect_identical(organisms$string_status[4], "unavailable")
    expect_true(is.na(organisms$string_species[4]))
    expect_identical(unique(organisms$uniprot_release), "2026_04")
})

test_that("the shipped files match the SHA-256 recorded in the reference list", {
    organisms <- readReferenceOrganisms()
    for (i in seq_len(nrow(organisms))) {
        files <- referenceOrganismFiles(organisms$uniprot_upid[i])
        expect_identical(digest::digest(file = files$fasta, algo = "sha256"), organisms$fasta_sha256[i])
        expect_identical(digest::digest(file = files$xref, algo = "sha256"), organisms$xref_sha256[i])
        expect_identical(digest::digest(file = files$go, algo = "sha256"), organisms$go_sha256[i])
    }
})

test_that("K-12: the eco KEGG entry, its locus, and the STRING entry are taken from multi-valued fields", {
    xref <- k167Xref(list("P0AE06", "acrA", "b0463 JW0452", "511145.b0463;", "ecj:JW0452;eco:b0463;ecoc:C3026_02270;"),
        list("P75837", "ycaO", "b1109 JW1095", "", "ecj:JW1095;eco:b1109;ecoc:C3026_06690;"))
    parsed <- parseReferenceCrossReferences(xref, kegg_code = "eco", string_species = "511145",
        locus_pattern = "^b[0-9]{4}$")
    expect_identical(parsed$reference_kegg_id, c("eco:b0463", "eco:b1109"))
    expect_identical(parsed$reference_locus, c("b0463", "b1109"))
    expect_identical(parsed$reference_string_id, c("511145.b0463", NA))
    expect_identical(parsed$reference_gene, c("acrA", "ycaO"))
})

test_that("K-12: an IS element with several eco entries keeps them all, joined in C-locale order, counted", {
    xref <- k167Xref(list("A0A385XJ53", "insA9", "b4709", "", "eco:b3444;eco:b0022;eco:b1894;ecoc:C3026_00105;"))
    parsed <- parseReferenceCrossReferences(xref, "eco", "511145", "^b[0-9]{4}$")
    expect_identical(parsed$reference_kegg_id, "eco:b0022;eco:b1894;eco:b3444")
    expect_identical(parsed$reference_locus, "b0022;b1894;b3444")
    expect_identical(parsed$reference_locus_values, 3L)
})

test_that("SF370: entries shared with another strain keep only the reference's own locus", {
    xref <- k167Xref(list("P0C0J1", "speB", "SPy_2039 M5005_Spy1735", "", "spy:SPy_2039;spz:M5005_Spy1735;"))
    parsed <- parseReferenceCrossReferences(xref, "spy", NA_character_, "^SPy_[0-9]{4}$")
    expect_identical(parsed$reference_kegg_id, "spy:SPy_2039")
    expect_identical(parsed$reference_locus, "SPy_2039")
    expect_true(is.na(parsed$reference_string_id))
})

test_that("MGH 78578: the KEGG locus wins over a gene_oln of another form; no gene_oln still gives a locus", {
    xref <- k167Xref(list("A6T923", "hpxO", "KPN78578_16330", "272620.KPN_01663;", "kpn:KPN_01663;"),
        list("A6TGC5", "priA", "", "272620.KPN_04230;", "kpn:KPN_04230;"))
    parsed <- parseReferenceCrossReferences(xref, "kpn", "272620", "^KPN_[0-9]{5}$")
    expect_identical(parsed$reference_locus, c("KPN_01663", "KPN_04230"))
    expect_identical(parsed$reference_string_id, c("272620.KPN_01663", "272620.KPN_04230"))
})

test_that("without a KEGG entry the gene_oln token matching the ordered-locus pattern is the locus", {
    xref <- k167Xref(list("Q2G0X1", "", "SAOUHSC_00001 other", "", ""))
    parsed <- parseReferenceCrossReferences(xref, "sao", "93061", "^SAOUHSC_[0-9]{5}$")
    expect_identical(parsed$reference_locus, "SAOUHSC_00001")
    expect_true(is.na(parsed$reference_kegg_id))
    expect_true(is.na(parsed$reference_gene))
})

test_that("the rows of the shipped capture parse as the hand-written rows do", {
    organisms <- readReferenceOrganisms()
    xref <- readReferenceCrossReferenceTable(referenceOrganismFiles("UP000000265")$xref)
    parsed <- parseReferenceCrossReferences(xref[xref$Entry %in% c("A6T923", "A6TGC5"), ], "kpn", "272620",
        organisms$locus_pattern[organisms$uniprot_upid == "UP000000265"])
    expect_identical(parsed$reference_locus, c("KPN_01663", "KPN_04230"))
})
