# Input manifest, retrieval, the record phase, the UniProt cache refusal and phase pre
# (KIT-159 5.1; section 10.3: T-IN-1 to T-IN-9, T-ANN-1, T-PRE-1).

test_that("T-IN-1 to T-IN-3: a missing column, a bad URL and a bad hash format are refused", {
    dir <- inputRepo()
    row <- manifestRow(dir, "open", "open.tsv")
    m1 <- writeManifest(row[, setdiff(names(row), "sha256")], file.path(dir, "m1.tsv"))
    expect_error(valReadManifest(m1), class = "val_manifest_invalid")
    bad_url <- row
    bad_url$url <- "ftp://example.org/open.tsv"
    expect_error(valReadManifest(writeManifest(bad_url, file.path(dir, "m2.tsv"))),
                 class = "val_manifest_invalid")
    bad_hash <- row
    bad_hash$sha256 <- "abc"
    expect_error(valReadManifest(writeManifest(bad_hash, file.path(dir, "m3.tsv"))),
                 class = "val_manifest_invalid")
    dup <- rbind(row, row)
    expect_error(valReadManifest(writeManifest(dup, file.path(dir, "m4.tsv"))),
                 class = "val_manifest_invalid")
    no_cite <- row
    no_cite$citation <- ""
    expect_error(valReadManifest(writeManifest(no_cite, file.path(dir, "m5.tsv"))),
                 class = "val_manifest_invalid")
})

test_that("T-IN-4 and T-IN-5: size and hash mismatches are refused and the file deleted", {
    dir <- inputRepo()
    dest <- withr::local_tempdir()
    row <- manifestRow(dir, "open", "open.tsv")
    size <- row
    size$bytes <- row$bytes + 1
    m <- writeManifest(size, file.path(dir, "m.tsv"))
    expect_error(valFetchInputs(m, dest, repo_root = dir), class = "val_input_size_mismatch")
    expect_length(list.files(dest, recursive = TRUE), 0L)
    hash <- row
    hash$sha256 <- strrep("a", 64)
    m <- writeManifest(hash, file.path(dir, "m.tsv"))
    expect_error(valFetchInputs(m, dest, repo_root = dir), class = "val_input_hash_mismatch")
    expect_length(list.files(dest, recursive = TRUE), 0L)
})

test_that("T-IN-6: a verified input is reached by role; an undeclared role is refused", {
    dir <- inputRepo()
    dest <- withr::local_tempdir()
    m <- writeManifest(manifestRow(dir, "open", "open.tsv"), file.path(dir, "m.tsv"))
    got <- valFetchInputs(m, dest, repo_root = dir)
    expect_identical(got$status, "verified")
    expect_identical(sha256File(valInputPath("open")), sha256File(file.path(dir, "data/open.tsv")))
    expect_error(valInputPath("nope"), class = "val_input_undeclared")
})

test_that("T-IN-7: 2.5 times the bytes to fetch above the free space is refused", {
    dir <- inputRepo()
    m <- writeManifest(manifestRow(dir, "open", "open.tsv"), file.path(dir, "m.tsv"))
    expect_error(valFetchInputs(m, withr::local_tempdir(), repo_root = dir, free_bytes = 10),
                 class = "val_disk_insufficient")
})

test_that("T-IN-8: three members of one archive from one download, each verified", {
    dir <- inputRepo()
    src <- file.path(withr::local_tempdir(), "pack")
    dir.create(src)
    for (f in c("a.tsv", "b.tsv", "c.tsv")) writeLines(c("x", f), file.path(src, f))
    archive <- file.path(dir, "data", "pack.tar.gz")
    withr::with_dir(src, utils::tar(archive, files = c("a.tsv", "b.tsv", "c.tsv"),
                                    compression = "gzip", tar = "internal"))
    members <- data.frame(member = c("a.tsv", "b.tsv", "c.tsv"),
                          sha256 = vapply(file.path(src, c("a.tsv", "b.tsv", "c.tsv")),
                                          sha256File, character(1)),
                          bytes = file.size(file.path(src, c("a.tsv", "b.tsv", "c.tsv"))))
    utils::write.table(members, file.path(dir, "members-pack.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    m <- writeManifest(manifestRow(dir, "pack", "pack.tar.gz", members = "members-pack.tsv"),
                       file.path(dir, "m.tsv"))
    got <- valFetchInputs(m, withr::local_tempdir(), repo_root = dir)
    expect_identical(sum(got$member %in% c("a.tsv", "b.tsv", "c.tsv")), 3L)
    expect_true(all(got$status == "verified"))
    expect_identical(readLines(valInputPath("pack", member = "b.tsv")), c("x", "b.tsv"))
    bad <- members
    bad$sha256[2] <- strrep("b", 64)
    utils::write.table(bad, file.path(dir, "members-pack.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    expect_error(valFetchInputs(m, withr::local_tempdir(), repo_root = dir),
                 class = "val_input_hash_mismatch")
})

test_that("T-IN-9: an optional role that cannot be fetched is recorded, a required one refused", {
    dir <- inputRepo()
    fail <- function(url, destfile) stop("HTTP 503")
    opt <- manifestRow(dir, "remote", "remote.tsv", url = "https://example.org/remote.tsv",
                       optional = TRUE)
    m <- writeManifest(rbind(manifestRow(dir, "open", "open.tsv"), opt), file.path(dir, "m.tsv"))
    got <- valFetchInputs(m, withr::local_tempdir(), repo_root = dir, download = fail,
                          backoff = 0)
    expect_identical(got$status[got$role == "remote"], "input_unavailable")
    expect_identical(got$status[got$role == "open"], "verified")
    req <- opt
    req$optional <- FALSE
    m2 <- writeManifest(req, file.path(dir, "m2.tsv"))
    expect_error(valFetchInputs(m2, withr::local_tempdir(), repo_root = dir, download = fail,
                                backoff = 0),
                 class = "val_input_unavailable")
})

test_that("the record phase hashes inputs and members; two equal records compare equal", {
    dir <- inputRepo()
    draft <- data.frame(dataset = "toy", role = "open", url = "repo://data/open.tsv", members = "",
                        licence = "CC BY 4.0", kind = "data", redistribution = "permitted",
                        citation = "Example, X1", source_note = "fixture", truth = FALSE)
    a <- withr::local_tempdir()
    b <- withr::local_tempdir()
    valRecordInputs(draft, withr::local_tempdir(), a, repo_root = dir)
    valRecordInputs(draft, withr::local_tempdir(), b, repo_root = dir)
    rec <- utils::read.delim(file.path(a, "inputs.tsv"), stringsAsFactors = FALSE)
    expect_identical(rec$sha256, sha256File(file.path(dir, "data/open.tsv")))
    expect_true(valCompareRecords(a, b)$identical)
    rec$sha256 <- strrep("c", 64)
    utils::write.table(rec, file.path(b, "inputs.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    expect_false(valCompareRecords(a, b)$identical)
})

test_that("T-ANN-1: an incomplete or empty UniProt answer refuses the record", {
    tbl <- data.frame(Protein.Ids = c("P69905", "P68871"))
    signal <- function(cls) {
        function(...) {
            warning(structure(class = c(cls, "warning", "condition"),
                              list(message = "mocked", call = NULL)))
            data.frame(Entry = "P69905")
        }
    }
    release <- function() "2026_03"
    expect_error(valRecordUniprotCache(tbl, "Protein.Ids", 9606, withr::local_tempdir(),
                                       annotate = signal("prot_annotation_uniprot_incomplete"),
                                       release = release),
                 class = "val_annotation_incomplete")
    expect_error(valRecordUniprotCache(tbl, "Protein.Ids", 9606, withr::local_tempdir(),
                                       annotate = signal("prot_annotation_uniprot_empty"),
                                       release = release),
                 class = "val_annotation_empty")
    expect_error(valRecordUniprotCache(tbl, "Protein.Ids", 9606, withr::local_tempdir(),
                                       annotate = function(...) data.frame(Entry = character()),
                                       release = release),
                 class = "val_annotation_empty")
    ok <- valRecordUniprotCache(tbl, "Protein.Ids", 9606, withr::local_tempdir(),
                                annotate = function(...) data.frame(Entry = c("P68871", "P69905")),
                                release = release)
    expect_identical(ok$release, "2026_03")
    expect_match(ok$sha256, "^[0-9a-f]{64}$")
})

test_that("T-PRE-1: phase pre admits listed aggregates only and refuses truth reads", {
    dir <- inputRepo()
    truth_row <- manifestRow(dir, "closed", "closed.tsv")
    m <- writeManifest(rbind(manifestRow(dir, "open", "open.tsv"), truth_row),
                       file.path(dir, "m.tsv"))
    valFetchInputs(m, withr::local_tempdir(), repo_root = dir, truth_roles = "closed")
    spec <- data.frame(file = c("counts.tsv", "counts.tsv"), column = c("group", "n"),
                       type = c("label", "count"), labels = c("all", ""))
    spec_path <- file.path(dir, "pre-outputs.tsv")
    utils::write.table(spec, spec_path, sep = "\t", quote = FALSE, row.names = FALSE)
    script <- function(body) {
        p <- tempfile(fileext = ".R")
        writeLines(body, p)
        p
    }
    good <- script(c("x <- utils::read.delim(valInputPath('open'))",
                     "utils::write.table(data.frame(group = 'all', n = nrow(x)),",
                     "  file.path(out_dir, 'counts.tsv'), sep = '\\t', quote = FALSE,",
                     "  row.names = FALSE)"))
    out <- withr::local_tempdir()
    expect_true(valRunPre(good, spec_path, out))
    extra <- script(c("utils::write.table(data.frame(group = 'all', n = 1, extra = 2),",
                      "  file.path(out_dir, 'counts.tsv'), sep = '\\t', quote = FALSE,",
                      "  row.names = FALSE)"))
    expect_error(valRunPre(extra, spec_path, withr::local_tempdir()),
                 "counts.tsv: columns differ from the declaration", class = "val_publish_forbidden",
                 fixed = TRUE)
    text <- script("writeLines('report', file.path(out_dir, 'feasibility.md'))")
    expect_error(valRunPre(text, spec_path, withr::local_tempdir()),
                 "wrote an unlisted file feasibility.md", class = "val_publish_forbidden",
                 fixed = TRUE)
    truth <- script("x <- valInputPath('closed')")
    expect_error(valRunPre(truth, spec_path, withr::local_tempdir()), class = "val_spec_invalid")
})
