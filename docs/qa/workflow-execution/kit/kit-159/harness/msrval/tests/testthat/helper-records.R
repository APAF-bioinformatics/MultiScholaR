# Fixtures for the records, licence and declaration tests (KIT-159 sections 5.1 to 5.3, 10.3).

sha256File <- function(path) digest::digest(file = path, algo = "sha256")

# A repository folder with input files and a manifest whose rows point at them by repo://.
inputRepo <- function(dir = withr::local_tempdir(.local_envir = parent.frame())) {
    dir.create(file.path(dir, "data"), showWarnings = FALSE)
    writeLines(c("id\tvalue", "a\t1", "b\t2"), file.path(dir, "data", "open.tsv"))
    writeLines(c("id\tvalue", "s1\t9"), file.path(dir, "data", "closed.tsv"))
    dir
}

manifestRow <- function(dir, role, file, kind = "data", redistribution = "permitted",
                        optional = FALSE, commit = FALSE, members = "", url = NULL,
                        citation = "Example data set, accession X1, doi:10.0/x") {
    path <- file.path(dir, "data", file)
    data.frame(
        dataset = "toy", role = role, file = file,
        url = url %||% paste0("repo://data/", file), accession = "X1",
        licence = if (redistribution == "permitted") "CC BY 4.0" else "none stated",
        kind = kind, redistribution = redistribution, citation = citation,
        sha256 = if (file.exists(path)) sha256File(path) else strrep("0", 64),
        bytes = if (file.exists(path)) file.size(path) else 1,
        members = members, keep_archive = FALSE, optional = optional, commit = commit,
        retrieved = "2026-11-23", source_note = "test fixture", stringsAsFactors = FALSE
    )
}

writeManifest <- function(rows, path) {
    utils::write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    path
}

# A small git repository for the declaration tests: a product file, a harness folder and an
# evidence folder; returns the path. Commits are made with a fixed identity.
gitRepo <- function(dir = withr::local_tempdir(.local_envir = parent.frame())) {
    git <- function(...) system2("git", shQuote(c("-C", dir, ...)), stdout = TRUE, stderr = TRUE)
    git("init", "-q")
    git("config", "user.email", "test@example.invalid")
    git("config", "user.name", "msrval test")
    dir.create(file.path(dir, "R"))
    writeLines("f <- function() 1", file.path(dir, "R", "f.R"))
    h <- file.path(dir, "docs/qa/workflow-execution/kit/kit-159/harness")
    dir.create(h, recursive = TRUE)
    writeLines("harness v1", file.path(h, "VERSION"))
    git("add", "-A")
    git("commit", "-q", "-m", "base")
    dir
}

gitDo <- function(dir, ...) {
    system2("git", shQuote(c("-C", dir, ...)), stdout = TRUE, stderr = TRUE)
}

gitCommitFiles <- function(dir, files, message = "change") {
    for (f in names(files)) {
        path <- file.path(dir, f)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        writeLines(files[[f]], path)
    }
    gitDo(dir, "add", "-A")
    gitDo(dir, "commit", "-q", "-m", message)
    gitDo(dir, "rev-parse", "HEAD")
}

# Licence-guard fixtures (shared by test-licence.R and test-evidence.R).
writeTsv <- function(df, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    path
}

emptyTable <- function(name, n = 1L) {
    cols <- EVIDENCE_COLUMNS[[name]]
    as.data.frame(stats::setNames(rep(list(rep("x", n)), length(cols)), cols))
}

# A licence fixture: one permitted and one not-permitted data input, and a GPL program.
licenceManifest <- function(dir) {
    rbind(manifestRow(dir, "open", "open.tsv", commit = TRUE),
          manifestRow(dir, "closed", "closed.tsv", redistribution = "not_permitted"),
          manifestRow(dir, "tool", "open.tsv", kind = "program", citation = "tool 1.0, GPL-3"))
}

evidenceDir <- function(manifest) {
    ev <- withr::local_tempdir(.local_envir = parent.frame())
    writeTsv(manifest[, c("role", "kind", "licence", "redistribution", "citation")],
             file.path(ev, "attribution.tsv"))
    writeTsv(emptyTable("cells.tsv"), file.path(ev, "cells.tsv"))
    ev
}

