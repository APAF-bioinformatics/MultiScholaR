# S0, the evidence run on generated inputs (KIT-159 10.4, r1.14; labelled "evidence run; not a
# lifecycle"; gates nothing, VAL-41): the accepted scripts of the ten routes, G1(b)'s alternatives
# that still exist after KIT-180 (VAL-70) and the two-species route, each DA table compared with
# its accepted twin. Its job's first script runs the new pieces' tests; valRunS0 refuses before any
# route unless they all ran and passed (VAL-54: an instrument fault ends the job).

#' Columns of one comparison of an alternative's DA table with the accepted script's
#' @noRd
ALT_COMPARE_COLUMNS <- c("digest_a", "digest_b", "identical", "rows_a", "rows_b",
                         "keys_one_side", "max_abs_diff", "na_mismatch", "called_a", "called_b",
                         "called_changed")

#' Compare an alternative's DA table with the accepted script's (G1(b): masked digest, then the
#' largest value difference and the change in called keys)
#'
#' Keys are assay x contrast x feature. Over the keys both tables hold: the largest absolute
#' difference of `log2fc`, `p` and `s_own` where both are present, the count of values present on
#' one side only, and the keys whose `called_own` differs.
#'
#' @param a The accepted script's `da.tsv` rows.
#' @param b The alternative's `da.tsv` rows.
#' @return A one-row data frame with the columns of `ALT_COMPARE_COLUMNS`.
#' @export
valCompareAlternative <- function(a, b) {
    key <- function(x) paste(x$assay, x$contrast, x$feature_id, sep = "\r")
    ka <- key(a)
    kb <- key(b)
    shared <- intersect(ka, kb)
    ia <- match(shared, ka)
    ib <- match(shared, kb)
    diffs <- 0
    na_mismatch <- 0L
    for (col in c("log2fc", "p", "s_own")) {
        x <- as.numeric(a[[col]][ia])
        y <- as.numeric(b[[col]][ib])
        na_mismatch <- na_mismatch + sum(is.na(x) != is.na(y))
        both <- !is.na(x) & !is.na(y)
        diffs <- c(diffs, abs(x[both] - y[both]))
    }
    called_a <- as.logical(a$called_own)
    called_b <- as.logical(b$called_own)
    digest_a <- valMaskedDigest(a)
    digest_b <- valMaskedDigest(b)
    data.frame(digest_a = digest_a, digest_b = digest_b, identical = identical(digest_a, digest_b),
               rows_a = nrow(a), rows_b = nrow(b),
               keys_one_side = length(setdiff(ka, kb)) + length(setdiff(kb, ka)),
               max_abs_diff = max(diffs), na_mismatch = as.integer(na_mismatch),
               called_a = sum(called_a, na.rm = TRUE), called_b = sum(called_b, na.rm = TRUE),
               called_changed = sum(called_a[ia] != called_b[ib], na.rm = TRUE),
               stringsAsFactors = FALSE)
}

#' The S0 plan (10.4, r1.14): fifteen route runs
#'
#' The ten routes' accepted scripts; export-then-load on `prot_dia`, `lipid_canonical` and
#' `metab_lc`; DIA_LIMPA without the covariates button; `prot_dia` on the two-species input
#' (compared with `prot_dia`: the same report, a FASTA of two organisms and the box ticked).
#'
#' @return A data frame: `cell`, `route`, `corpus`, `alternative`, `compare_with`.
#' @export
valS0Plan <- function() {
    routes <- valRouteTable()$route
    alt_routes <- c("prot_dia", "lipid_canonical", "metab_lc", "prot_dia_limpa")
    alts <- c("export_load", "export_load", "export_load", "no_covariates")
    data.frame(
        cell = c(routes, paste0(alt_routes, "-", alts), PARITY_TWO_SPECIES),
        route = c(routes, alt_routes, "prot_dia"),
        corpus = c(vapply(c(routes, alt_routes), valCorpusFor, character(1), USE.NAMES = FALSE),
                   PARITY_TWO_SPECIES),
        alternative = c(rep(NA_character_, length(routes)), alts, NA_character_),
        compare_with = c(rep(NA_character_, length(routes)), alt_routes, "prot_dia"),
        stringsAsFactors = FALSE
    )
}

#' Did every named test file run, with no test failed, errored or skipped?
#'
#' @param path The `harness-tests.tsv` of the job's first script.
#' @param files The test files required (names after `test-`, without `.R`).
#' @return `TRUE` or `FALSE`, with attribute `reason`.
#' @export
valS0TestsPassed <- function(path, files) {
    no <- function(reason) structure(FALSE, reason = reason)
    if (is.null(path) || !file.exists(path)) return(no("no test table"))
    t <- utils::read.delim(path, stringsAsFactors = FALSE)
    if (!nrow(t)) return(no("no test ran"))
    missing <- setdiff(paste0("test-", files, ".R"), t$file)
    if (length(missing)) return(no(paste("not run:", paste(missing, collapse = ", "))))
    bad <- t$failed > 0 | t$error > 0 | t$skipped > 0
    if (any(bad)) return(no(sprintf("%d of %d tests did not pass", sum(bad), nrow(t))))
    structure(TRUE, reason = "passed")
}

#' Run S0 (10.4, r1.14) into `out_dir`: corpus.tsv, routes.tsv, steps.tsv, modal.tsv,
#' alternatives.tsv, the ledgers and outputs.sha256
#'
#' Refuses (`val_spec_invalid`) before writing a corpus or running a route unless the new pieces'
#' tests all ran and passed. Each route run is the accepted (or alternative) script under
#' [valRunWorkflow()] with the shipped arm; an error in one run is recorded and the next runs.
#'
#' @param out_dir The evidence folder.
#' @param pkg_root The package checkout.
#' @param tests `list(path, files)` for [valS0TestsPassed()].
#' @param work_dir A private work folder (corpora, project folders, tables).
#' @param cells Plan cells to run (default all).
#' @return The routes table, invisibly.
#' @export
valRunS0 <- function(out_dir, pkg_root, tests, work_dir, cells = NULL) {
    gate <- valS0TestsPassed(tests$path, tests$files)
    if (!isTRUE(gate)) {
        reason <- attr(gate, "reason")  # nolint: object_usage_linter (read by cli's glue)
        valAbort(paste("S0 stops before any route: the new pieces' tests did not pass ({reason});",
                       "an instrument fault (VAL-54)"), "val_spec_invalid")
    }
    plan <- valS0Plan()
    if (!is.null(cells)) plan <- plan[plan$cell %in% cells, , drop = FALSE]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    corpora <- list()
    for (id in unique(plan$corpus)) {
        corpora[[id]] <- valWriteCorpus(id, file.path(work_dir, "corpus", id), pkg_root)
    }
    write <- function(x, name) {
        if (is.null(x)) x <- data.frame(none = character())
        utils::write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
                           row.names = FALSE, na = "NA")
    }
    write(do.call(rbind, lapply(corpora, function(x) {
        cbind(x$manifest, generator = x$generator, arguments = x$arguments)
    })), "corpus.tsv")
    helpers <- valGuiHelpers(pkg_root)
    results <- lapply(seq_len(nrow(plan)), function(i) {
        valRunS0Cell(plan[i, , drop = FALSE], corpora[[plan$corpus[[i]]]], helpers, pkg_root,
                     work_dir, out_dir)
    })
    names(results) <- plan$cell
    routes <- do.call(rbind, lapply(results, `[[`, "route"))
    write(routes, "routes.tsv")
    write(do.call(rbind, lapply(results, `[[`, "steps")), "steps.tsv")
    write(do.call(rbind, lapply(results, `[[`, "modal")), "modal.tsv")
    write(valS0Alternatives(plan, results), "alternatives.tsv")
    valWriteOutputsSha(out_dir)
    invisible(routes)
}

#' One S0 route run and its rows (route, steps, modal inputs) and DA table
#' @noRd
valRunS0Cell <- function(row, corpus, helpers, pkg_root, work_dir, out_dir) {
    alternative <- if (is.na(row$alternative)) NULL else row$alternative
    cell <- valCorpusCell(row$route, corpus, alternative)
    cell_dir <- file.path(work_dir, "cells", row$cell)
    tables <- file.path(cell_dir, "tables")
    valStatus(phase = "s0", arm = row$cell, status = "started")
    started <- Sys.time()
    res <- tryCatch(
        valRunWorkflow(cell$script, files = cell$files, design = cell$design, helpers = helpers,
                       forced_inputs = cell$forced_inputs, uniprot_cache = cell$uniprot_cache,
                       annotation = cell$annotation, tables_dir = tables,
                       ledger_path = file.path(out_dir, "ledgers", paste0(row$cell, ".tsv")),
                       pkg_root = pkg_root, work_dir = file.path(cell_dir, "work")),
        error = function(e) e
    )
    seconds <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
    out <- data.frame(cell = row$cell, route = row$route, corpus = row$corpus,
                      alternative = row$alternative, exit = "error", reason = NA_character_,
                      steps_ok = NA_integer_, steps = NA_integer_, fallback = NA_character_,
                      da_rows = NA_integer_, da_called = NA_integer_, da_digest = NA_character_,
                      ledger_sha256 = NA_character_, seconds = seconds, message = NA_character_,
                      stringsAsFactors = FALSE)
    if (inherits(res, "error")) {
        out$reason <- class(res)[[1L]]
        out$message <- substr(gsub("[\t\r\n]+", " ", valMaskMessage(conditionMessage(res))), 1L,
                              300L)
        valStatus(phase = "s0", arm = row$cell, status = "error", class = out$reason)
        return(list(route = out, steps = NULL, modal = NULL, da = NULL))
    }
    da_path <- file.path(tables, "da.tsv")
    da <- if (identical(res$exit, "completed") && file.exists(da_path)) {
        utils::read.delim(da_path, stringsAsFactors = FALSE)
    }
    out$exit <- res$exit
    out$reason <- res$reason
    out$steps <- NROW(res$steps)
    out$steps_ok <- sum(res$steps$outcome == "ok")
    out$fallback <- res$fallback
    out$ledger_sha256 <- res$ledger_sha256
    if (!is.null(da)) {
        out$da_rows <- nrow(da)
        out$da_called <- sum(as.logical(da$called_own), na.rm = TRUE)
        out$da_digest <- valMaskedDigest(da)
    }
    steps <- if (!is.null(res$steps)) {
        cbind(cell = row$cell, res$steps[, c("step", "outcome", "class", "reason", "head",
                                             "message")], stringsAsFactors = FALSE)
    }
    led <- res$ledger
    pick <- led$location %in% "modal" | endsWith(led$input_id, "-mixed_species_fasta")
    modal <- led[pick, c("step", "action", "input_id", "value", "source", "location"), drop = FALSE]
    modal <- modal[!duplicated(modal$input_id, fromLast = TRUE), , drop = FALSE]
    modal <- if (nrow(modal)) cbind(cell = row$cell, modal, stringsAsFactors = FALSE)
    valStatus(phase = "s0", arm = row$cell, status = out$exit)
    list(route = out, steps = steps, modal = modal, da = da)
}

#' The alternatives table: each planned comparison, or `not_compared` when a side has no DA table
#' @noRd
valS0Alternatives <- function(plan, results) {
    rows <- lapply(which(!is.na(plan$compare_with)), function(i) {
        base <- data.frame(cell = plan$cell[[i]], compare_with = plan$compare_with[[i]],
                           stringsAsFactors = FALSE)
        a <- results[[plan$compare_with[[i]]]]$da
        b <- results[[plan$cell[[i]]]]$da
        if (is.null(a) || is.null(b)) {
            empty <- as.data.frame(as.list(stats::setNames(rep(NA, length(ALT_COMPARE_COLUMNS)),
                                                           ALT_COMPARE_COLUMNS)))
            return(cbind(base, status = "not_compared", empty, stringsAsFactors = FALSE))
        }
        cbind(base, status = "compared", valCompareAlternative(a, b), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
