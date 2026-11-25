# R1's record-phase UniProt wiring (KIT-159 10.4 R1, r1.12; 5.1): valRecordUniprot on a committed
# accession list writes one allow-listed `uniprot.tsv` (taxon, accessions, rows, release,
# release_source, sha256); the cache itself stays private; refusals pass through.

accFile <- function(dir) {
    path <- file.path(dir, "uniprot-accessions.tsv")
    writeLines(c("accession\ttaxon\tgene", "P69905\t9606\tHBA1", "P68871\t9606\tHBB"), path)
    path
}

test_that("T-R1-1: the record writes uniprot.tsv with the release and hash; the cache stays private", {
    dir <- withr::local_tempdir()
    out <- file.path(dir, "out")
    cache <- file.path(dir, "cache")
    got <- new.env()
    annotate <- function(data_tbl, protein_id_column, taxon_id, cache_dir) {
        got$ids <- data_tbl[[protein_id_column]]
        got$taxon <- taxon_id
        data.frame(Entry = c("P69905", "P68871"), Gene = c("HBA1", "HBB"))
    }
    res <- valRecordUniprot(accFile(dir), out, cache, annotate = annotate,
                            release = function() "2026_04")
    expect_setequal(got$ids, c("P69905", "P68871"))
    expect_identical(as.character(got$taxon), "9606")
    tab <- utils::read.delim(file.path(out, "uniprot.tsv"), colClasses = "character")
    expect_identical(names(tab), msrval:::EVIDENCE_COLUMNS[["uniprot.tsv"]])
    expect_identical(tab$release, "2026_04")
    expect_identical(tab$accessions, "2")
    expect_true(grepl("^[0-9a-f]{64}$", tab$sha256))
    expect_identical(list.files(out), "uniprot.tsv")
})

test_that("T-R1-2: uniprot.tsv is on the allow-list with its columns; an extra column is refused", {
    m <- licenceManifest(inputRepo())
    ev <- evidenceDir(m)
    cols <- msrval:::EVIDENCE_COLUMNS[["uniprot.tsv"]]
    utils::write.table(as.data.frame(stats::setNames(as.list(rep("x", length(cols))), cols)),
                       file.path(ev, "uniprot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    expect_true(valCheckEvidence(ev, m))
    utils::write.table(data.frame(taxon = "9606", extra = "1"), file.path(ev, "uniprot.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    expect_error(valCheckEvidence(ev, m), class = "val_publish_forbidden")
})

test_that("T-R1-3: an incomplete UniProt answer refuses the record and writes nothing", {
    dir <- withr::local_tempdir()
    out <- file.path(dir, "out")
    annotate <- function(...) {
        warning(structure(class = c("prot_annotation_uniprot_incomplete", "warning", "condition"),
                          list(message = "mocked", call = NULL)))
        data.frame(Entry = "P69905")
    }
    expect_error(suppressWarnings(valRecordUniprot(accFile(dir), out, file.path(dir, "c"),
                                                   annotate = annotate,
                                                   release = function() "2026_04")),
                 class = "val_annotation_incomplete")
    expect_false(file.exists(file.path(out, "uniprot.tsv")))
})
