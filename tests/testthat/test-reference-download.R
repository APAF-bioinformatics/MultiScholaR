# KIT-167 sections 4 and 8.4 (A10): the explicit "current UniProt release" choice. The service is replaced by a stub
# request function; nothing here calls UniProt.

k167Fasta <- function(n, ox = 301447L) {
    paste0(sprintf(">sp|P%05d|PRT%d_STRPY Protein %d OS=Streptococcus pyogenes OX=%d PE=1 SV=1\nMKTAYIAKQRQISFVKSHFSRQ%03d",
        seq_len(n), seq_len(n), seq_len(n), ox, seq_len(n)), collapse = "\n")
}
k167Tsv <- function(header, n, row) paste(c(header, vapply(rev(seq_len(n)), row, "")), collapse = "\n")
k167Replies <- function(n = 60L, release = c("2026_05", "2026_05", "2026_05"), status = c(200L, 200L, 200L),
                        xref_rows = n) {
    replies <- list(
        fasta = list(status = status[1], release = release[1], body = k167Fasta(n)),
        xref = list(status = status[2], release = release[2], body = k167Tsv(
            "Entry\tGene Names (primary)\tGene Names (ordered locus)\tSTRING\tKEGG", xref_rows,
            function(i) sprintf("P%05d\tg%d\tSPy_%04d\t\tspy:SPy_%04d;", i, i, i, i))),
        go = list(status = status[3], release = release[3], body = k167Tsv(paste("Entry\tEntry Name\tProtein names",
            "Gene Names\tOrganism\tLength\tGene Ontology IDs\tReviewed\tProtein existence\tAnnotation", sep = "\t"), n,
            function(i) sprintf("P%05d\tPRT%d_STRPY\tProtein %d\tg%d\tStreptococcus pyogenes\t25\tGO:0005737\treviewed\tEvidence at protein level\t5.0",
                i, i, i, i))))
    function(url) {
        part <- if (grepl("format=fasta", url, fixed = TRUE)) "fasta" else if (grepl("go_id", url, fixed = TRUE)) "go" else "xref"
        reply <- replies[[part]]
        list(status = reply$status, headers = if (is.na(reply$release)) list() else
            list(`x-uniprot-release` = reply$release), body = reply$body)
    }
}

test_that("the current release of a shipped reference is downloaded, sorted, hashed and labelled", {
    dest <- tempfile("k167-download-")
    urls <- character()
    stub <- k167Replies()
    reference <- downloadReferenceProteome("UP000000750", dest, request_fn = function(url) {
        urls <<- c(urls, url)
        stub(url)
    })
    expect_identical(reference$release, "2026_05")
    expect_identical(reference$upid, "UP000000750")
    expect_true(all(file.exists(c(reference$fasta, reference$xref, reference$go))))
    expect_identical(dirname(reference$fasta), file.path(normalizePath(dest), "UP000000750_2026_05"))
    xref <- readReferenceCrossReferenceTable(reference$xref)
    expect_identical(xref$Entry, sort(xref$Entry, method = "radix"))
    expect_match(reference$name, "not the validated reference", fixed = TRUE)
    expect_identical(reference$kegg_code, "spy")
    expect_true(is.na(reference$string_species))
    expect_identical(reference$uniprot_taxon, 301447L)
    expect_identical(reference$fasta_sha256, digest::digest(file = reference$fasta, algo = "sha256"))
    expect_true(all(grepl("^https://rest.uniprot.org/uniprotkb/stream\\?query=proteome:UP000000750", urls)))
    resolved <- resolveMappingReference(reference)
    expect_true(resolved$supplied)
    expect_identical(resolved$release, "2026_05")
})

test_that("an absent httr2 is refused in words", {
    local_mocked_bindings(optionalPackageAvailable = function(package) FALSE)
    expect_error(downloadReferenceProteome("UP000000750", tempfile(), request_fn = k167Replies()),
        class = "reference_mapping_download_unavailable")
})

test_that("a failed reply, a missing release header, differing releases or an implausible count are refused", {
    expect_error(downloadReferenceProteome("UP000000750", tempfile(),
        request_fn = k167Replies(status = c(200L, 503L, 200L))), class = "reference_mapping_download_failed")
    expect_error(downloadReferenceProteome("UP000000750", tempfile(),
        request_fn = k167Replies(release = c("2026_05", NA, "2026_05"))), class = "reference_mapping_download_failed")
    expect_error(downloadReferenceProteome("UP000000750", tempfile(),
        request_fn = k167Replies(release = c("2026_05", "2026_05", "2026_06"))), class = "reference_mapping_download_failed")
    expect_error(downloadReferenceProteome("UP000000750", tempfile(), request_fn = k167Replies(xref_rows = 59L)),
        class = "reference_mapping_download_failed")
    expect_error(downloadReferenceProteome("UP000000750", tempfile(), request_fn = k167Replies(n = 10L)),
        class = "reference_mapping_download_failed")
    expect_error(downloadReferenceProteome("UP000000750", tempfile(), request_fn = function(url) stop("timeout")),
        class = "reference_mapping_download_failed")
})

test_that("a proteome outside the shipped list needs its service identifiers", {
    expect_error(downloadReferenceProteome("UP000999999", tempfile(), request_fn = k167Replies()),
        class = "reference_mapping_reference_invalid")
    reference <- downloadReferenceProteome("UP000999999", tempfile(), request_fn = k167Replies(),
        name = "Other strain", string_species = NA_character_, kegg_code = "xyz")
    expect_identical(reference$kegg_code, "xyz")
})
