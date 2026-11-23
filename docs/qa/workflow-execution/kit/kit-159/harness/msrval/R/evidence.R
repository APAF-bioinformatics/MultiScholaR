# The evidence writer's hashes and phase `score` (KIT-159 sections 5.7, 7.8).

#' Write `outputs.sha256` over every file of a folder (written last)
#'
#' @param dir The folder.
#' @return The path of `outputs.sha256`, invisibly.
#' @export
valWriteOutputsSha <- function(dir) {
    files <- sort(setdiff(list.files(dir, recursive = TRUE), "outputs.sha256"))
    lines <- sprintf("%s  %s", vapply(file.path(dir, files), valSha256, character(1)), files)
    path <- file.path(dir, "outputs.sha256")
    writeLines(lines, path)
    invisible(path)
}

#' Verify a folder against its `outputs.sha256`
#' @noRd
valVerifyOutputsSha <- function(dir) {
    path <- file.path(dir, "outputs.sha256")
    if (!file.exists(path)) return(FALSE)
    lines <- readLines(path, warn = FALSE)
    sha <- sub("  .*$", "", lines)
    files <- sub("^[0-9a-f]{64}  ", "", lines)
    present <- setdiff(list.files(dir, recursive = TRUE), "outputs.sha256")
    setequal(files, present) &&
        all(vapply(seq_along(files), function(i) {
            identical(valSha256(file.path(dir, files[[i]])), sha[[i]])
        }, logical(1)))
}

#' Phase score: combine the per-unit rows of several requests of one counted run
#'
#' Each `runs/<label>/` folder is verified against its `outputs.sha256` and its provenance (the same
#' `package_commit`, harness tree and declaration commit); every expected unit x arm cell must be
#' present exactly once across the folders' `cells.tsv`.
#'
#' @param dir The ticket's evidence folder (holding `runs/`).
#' @param labels The requests' labels.
#' @param expected_cells Data frame of the declared cells (`unit`, `arm`).
#' @return A list: `metrics` and `cells` combined, `package_commit`.
#' @export
valScoreRuns <- function(dir, labels, expected_cells) {
    changed <- function(what) valAbort("phase score: {what}", "val_declaration_changed")
    provs <- lapply(labels, function(l) {
        d <- file.path(dir, "runs", l)
        if (!valVerifyOutputsSha(d)) changed(paste(l, "does not match its outputs.sha256"))
        jsonlite::fromJSON(file.path(d, "provenance.json"), simplifyVector = FALSE)
    })
    ids <- vapply(provs, function(p) {
        paste(p$code$package_commit %||% "", p$code$harness$tree %||% "",
              p$code$declaration$commit %||% "")
    }, character(1))
    if (length(unique(ids)) != 1L) changed("the requests differ in commit, harness or declaration")
    read <- function(l, f) {
        utils::read.delim(file.path(dir, "runs", l, f), colClasses = "character",
                          na.strings = character(), stringsAsFactors = FALSE)
    }
    cells <- do.call(rbind, lapply(labels, function(l) cbind(read(l, "cells.tsv"), request = l)))
    metrics <- do.call(rbind, lapply(labels, function(l) {
        cbind(read(l, "metrics.tsv"), request = l)
    }))
    key_c <- paste(cells$unit, cells$arm)
    key_e <- paste(expected_cells$unit, expected_cells$arm)
    if (!setequal(key_c, key_e) || anyDuplicated(key_c)) {
        valAbort("phase score: the cells are not each declared cell exactly once",
                 "val_spec_invalid")
    }
    list(metrics = metrics, cells = cells, package_commit = provs[[1L]]$code$package_commit)
}
