# The run folder of a cell phase (KIT-159 5.2, 5.6, 5.7; audit I07, I08): the cells' results and
# values, calls.tsv, coverage.tsv, ledgers, the table hashes, scored.tsv.gz and provenance blocks.

#' The cells' results: every declared cell, filled from the chunks' rows
#' @noRd
valCellResults <- function(cells, ran) {
    res <- cells[, c("unit", "run", "arm")]
    res$exit <- "refused"
    res$stop_reason <- NA_character_
    res$step <- "identity_gate"
    res$class <- "val_kernel_identity_failed"
    res$message_sha256 <- NA_character_
    res$seconds <- NA_real_
    res$chunk <- NA_character_
    res$random_state_changed <- NA
    res$vmhwm_mib <- NA_real_
    if (!is.null(ran)) {
        j <- match(res$run, ran$run)
        hit <- !is.na(j)
        for (col in c("exit", "stop_reason", "step", "class", "message_sha256", "seconds",
                      "chunk", "random_state_changed", "vmhwm_mib")) {
            if (!is.null(ran[[col]])) res[[col]][hit] <- ran[[col]][j[hit]]
        }
    }
    rownames(res) <- NULL
    res
}

#' Each run's returned value, read back from its chunk's `<run>.rds`
#' @noRd
valCellValues <- function(res, chunk_root) {
    runs <- unique(res$run)
    lapply(stats::setNames(runs, runs), function(r) {
        ch <- res$chunk[match(r, res$run)]
        if (is.na(ch)) return(NULL)
        p <- file.path(chunk_root, paste0("chunk-", ch), paste0(r, ".rds"))
        if (!file.exists(p)) return(NULL)
        readRDS(p)$value
    })
}

#' A value as one line of compact JSON, or NA
#' @noRd
valOneLineJson <- function(x) {
    if (is.null(x)) return(NA_character_)
    txt <- tryCatch(as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")),
                    error = function(e) NA_character_)
    gsub("[\t\r\n]+", " ", txt)
}

#' cells.tsv's fallback, options at start, RUV k and controls from each cell's value and ruv.tsv
#' @noRd
valCellSettings <- function(res, values) {
    one <- function(r, f) {
        v <- values[[r]]
        f(v)
    }
    res$fallback <- vapply(res$run, one, character(1), f = function(v) {
        x <- v$fallback
        if (is.null(x) || length(x) != 1L) NA_character_ else as.character(x)
    })
    res$options_at_start <- vapply(res$run, one, character(1), f = function(v) {
        valOneLineJson(v$options_at_start)
    })
    ruv <- lapply(res$run, function(r) {
        p <- values[[r]]$tables[["ruv"]]
        if (is.null(p) || !file.exists(p)) return(c(NA_integer_, NA_integer_))
        x <- utils::read.delim(p, stringsAsFactors = FALSE)
        c(suppressWarnings(as.integer(x$k[[1L]])), suppressWarnings(as.integer(x$controls_n[[1L]])))
    })
    res$ruv_k <- vapply(ruv, `[[`, integer(1), 1L)
    res$ruv_controls <- vapply(ruv, `[[`, integer(1), 2L)
    res
}

#' A unit's truth tables: the data set's, plus its `per_unit` truths built by the ticket's
#' function on the unit's inputs (4.1, KIT-163 N5). They are built in the scoring process, not as
#' the unit's first cell; outcome.md records the divergence.
#' @noRd
valUnitTruth <- function(decl, dataset, env, prep, unit) {
    ds <- valDeclDataset(decl, dataset)
    tr <- prep$truth
    for (t1 in ds[["truth"]] %||% list()) {
        if (!is.list(t1) || is.null(t1[["per_unit"]])) next
        fn <- env[[t1[["per_unit"]]]]
        if (!is.function(fn)) {
            valAbort("per-unit truth {t1$id}: the ticket's code defines no {t1$per_unit}()",
                     "val_spec_invalid")
        }
        got <- fn(unit, prep$units[[unit]])
        valCheckTruth(got, t1[["id"]])
        tr[[t1[["id"]]]] <- got
    }
    tr
}

#' Units whose call counts stay private in the lifecycle: every scored unit of a data set with a
#' primary metric, except `not_under_test` units (4.8, A28)
#' @noRd
valPrimaryUnits <- function(decl, dataset, unit_ids) {
    has_primary <- any(vapply(decl[["metrics"]] %||% list(), function(m) {
        !is.null(m[["primary_family"]]) || isTRUE(m[["primary"]])
    }, logical(1)))
    if (!has_primary) return(character())
    nut <- unlist(valDeclDataset(decl, dataset)[["not_under_test"]][["units"]])
    setdiff(unit_ids, nut)
}

#' calls.tsv: per completed cell and contrast the call statistic and the q-value fallbacks
#' (`qvalue_bh_fallback` where the product's DA result flags the contrast; the count is NA where the
#' product's DA object was not read: off proteomics, custom and kernel arms)
#' @noRd
valWriteCalls <- function(res, values, cell_root, run_dir) {
    rows <- list()
    for (i in which(res$exit == "completed")) {
        p <- file.path(valCellDir(cell_root, res$run[[i]]), "da.tsv")
        if (!file.exists(p)) next
        da <- utils::read.delim(p, stringsAsFactors = FALSE)
        fb <- values[[res$run[[i]]]]$qvalue_fallback
        for (ct in unique(da$contrast)) {
            cs <- unique(stats::na.omit(da$call_statistic[da$contrast == ct]))
            cs <- if (length(cs)) as.character(cs[[1L]]) else NA_character_
            if (!is.null(fb) && ct %in% fb) cs <- "qvalue_bh_fallback"
            rows[[length(rows) + 1L]] <- data.frame(
                unit = res$unit[[i]], arm = res$arm[[i]], contrast = ct, call_statistic = cs,
                fallback_count = if (is.null(fb)) NA_integer_ else length(fb),
                stringsAsFactors = FALSE
            )
        }
    }
    out <- as.data.frame(stats::setNames(replicate(5L, character(), simplify = FALSE),
                                         EVIDENCE_COLUMNS[["calls.tsv"]]))
    if (length(rows)) out <- do.call(rbind, rows)
    utils::write.table(out, file.path(run_dir, "calls.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
    invisible(out)
}

#' coverage.tsv: M11 per truth, unit, arm and contrast in the own view (4.1, 4.3); in the
#' lifecycle the call counts of primary units are not written (A28)
#' @noRd
valWriteCoverageTable <- function(decl, done, cell_root, truths, universe, env, run_dir,
                                  blind_seed, primary_units) {
    rows <- list()
    for (u in unique(done$unit)) {
        r <- done[done$unit == u, , drop = FALSE]
        paths <- file.path(valCellDir(cell_root, r$run), "da.tsv")
        keep <- file.exists(paths)
        if (!any(keep)) next
        tabs <- valArmTables(stats::setNames(as.list(paths[keep]), r$arm[keep]), env)
        for (tid in names(truths[[u]])) {
            tr <- truths[[u]][[tid]]
            if (!is.null(blind_seed)) tr <- valBlindTruth(tr, as.integer(blind_seed))
            s <- valScoreTable(tabs, tr[, c("key", "class", "expected", "truth_status")],
                               universe, unit = u, truth_id = tid)
            for (a in names(tabs)) {
                mat <- valCellMatrix(file.path(dirname(paths[keep][r$arm[keep] == a][[1L]]),
                                               "matrix.tsv"))
                for (ct in unique(s$contrast[s$arm == a])) {
                    cv <- valCoverage(s[s$arm == a & s$contrast == ct, , drop = FALSE], mat)
                    if (is.null(cv) || nrow(cv) == 0L) next
                    if (!is.null(blind_seed) && u %in% primary_units) {
                        cv$called_own <- NA_integer_
                        cv$called_bh <- NA_integer_
                    }
                    id <- data.frame(truth_id = tid, unit = u, arm = a, contrast = ct,
                                     stringsAsFactors = FALSE)
                    rows[[length(rows) + 1L]] <- cbind(id, cv)
                }
            }
        }
    }
    out <- as.data.frame(stats::setNames(replicate(11L, character(), simplify = FALSE),
                                         EVIDENCE_COLUMNS[["coverage.tsv"]]))
    if (length(rows)) out <- do.call(rbind, rows)
    utils::write.table(out[, EVIDENCE_COLUMNS[["coverage.tsv"]]],
                       file.path(run_dir, "coverage.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
    invisible(out)
}

#' Copy each workflow cell's ledger as `ledger-<route>-<arm>-<unit>.tsv` (5.2 "always committed";
#' the arm and unit in the name keep two arms on one route apart)
#' @noRd
valCopyLedgers <- function(res, cells, cell_root, run_dir) {
    safe <- function(x) gsub("[^A-Za-z0-9._]", "_", x)
    for (i in which(!duplicated(res$run) & cells$kind == "workflow")) {
        route <- cells$route[[i]]
        src <- file.path(valCellDir(cell_root, res$run[[i]]), "work",
                         paste0("ledger-", route, ".tsv"))
        if (!file.exists(src)) next
        file.copy(src, file.path(run_dir, sprintf("ledger-%s-%s-%s.tsv", route, safe(res$arm[[i]]),
                                                  safe(res$unit[[i]]))), overwrite = TRUE)
    }
    invisible(TRUE)
}

#' scored.tsv.gz under 5.2's rule; returns the `not-committed.sha256` line when withheld
#' (lifecycle, a not-permitted source, or over 5 MB)
#' @noRd
valWriteScored <- function(done, cell_root, truths, universe, env, run_dir, lifecycle, private) {
    parts <- list()
    for (u in unique(done$unit)) {
        r <- done[done$unit == u, , drop = FALSE]
        paths <- file.path(valCellDir(cell_root, r$run), "da.tsv")
        keep <- file.exists(paths)
        if (!any(keep)) next
        tabs <- valArmTables(stats::setNames(as.list(paths[keep]), r$arm[keep]), env)
        for (tid in names(truths[[u]])) {
            tr <- truths[[u]][[tid]][, c("key", "class", "expected", "truth_status")]
            parts[[length(parts) + 1L]] <- valScoreTable(tabs, tr, universe, unit = u,
                                                         truth_id = tid)
        }
    }
    if (!length(parts)) return(character())
    s <- do.call(rbind, parts)[, SCORED_COLUMNS]
    held <- lifecycle || private
    target <- if (held) file.path(cell_root, "scored.tsv.gz") else
        file.path(run_dir, "scored.tsv.gz")
    con <- gzfile(target, "w")
    utils::write.table(s, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    close(con)
    if (!held && file.size(target) <= MAX_FILE_BYTES) return(character())
    line <- sprintf("%s  scored.tsv.gz  rows=%d", valSha256(target), nrow(s))
    if (!held) unlink(target)
    line
}

#' tables.sha256 (every cell's private standard tables) and not-committed.sha256 (feature-level
#' files withheld by the licence rule); returns the table rows for provenance
#' @noRd
valWriteTableHashes <- function(res, values, run_dir, withheld, prepared) {
    rows <- list()
    for (r in unique(res$run[res$exit == "completed"])) {
        v <- values[[r]]
        tabs <- v$tables %||% (if (!is.null(v$da)) c(da = v$da))
        for (nm in names(tabs)) {
            p <- tabs[[nm]]
            if (!file.exists(p)) next
            rows[[length(rows) + 1L]] <- data.frame(
                run = r, name = nm, sha256 = valSha256(p),
                rows = max(0L, length(readLines(p, warn = FALSE)) - 1L), stringsAsFactors = FALSE
            )
        }
    }
    tab <- if (length(rows)) do.call(rbind, rows) else NULL
    lines <- character()
    if (!is.null(tab)) {
        lines <- sprintf("%s  %s/%s.tsv  rows=%d", tab$sha256,
                         gsub("[^A-Za-z0-9._@-]", "_", tab$run), tab$name, tab$rows)
    }
    writeLines(lines, file.path(run_dir, "tables.sha256"))
    nc <- c(withheld, if (!is.null(prepared) && nrow(prepared)) {
        sprintf("%s  prepared/%s", prepared$sha256, basename(prepared$path))
    })
    writeLines(as.character(nc), file.path(run_dir, "not-committed.sha256"))
    tab %||% list()
}

#' The provenance `truth` block: per truth id its SHA-256 (of the key-sorted table), rows and
#' counts by truth status; and the blinding
#' @noRd
valTruthProvenance <- function(truth, work_dir, blind_seed) {
    out <- lapply(truth, function(t1) {
        p <- tempfile(tmpdir = work_dir, fileext = ".tsv")
        utils::write.table(t1[order(as.character(t1$key), method = "radix"), , drop = FALSE], p,
                           sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        sh <- valSha256(p)
        unlink(p)
        st <- table(as.character(t1$truth_status))
        list(sha256 = sh, rows = nrow(t1),
             by_truth_status = as.list(stats::setNames(as.integer(st), names(st))))
    })
    c(out, list(blind = list(seed = blind_seed %||% NA_integer_, applied = !is.null(blind_seed))))
}

#' The provenance `cells` block, one entry per unit x arm (5.6)
#' @noRd
valCellsProvenance <- function(res, cells, values, run_dir, cell_root, hook_declared) {
    lapply(seq_len(nrow(res)), function(i) {
        v <- values[[res$run[[i]]]]
        led <- file.path(valCellDir(cell_root, res$run[[i]]), "work",
                         paste0("ledger-", cells$route[[i]], ".tsv"))
        tabs <- v$tables %||% (if (!is.null(v$da)) c(da = v$da))
        hook_none <- if (hook_declared && res$exit[[i]] == "completed") "not_run" else "none"
        hook <- v$hook %||% list(status = hook_none)
        list(unit = res$unit[[i]], arm = res$arm[[i]], kind = cells$kind[[i]],
             chunk = res$chunk[[i]], seed = as.integer(cells$seed[[i]]),
             rng_kind = "Mersenne-Twister/Inversion/Rejection",
             random_state_changed = res$random_state_changed[[i]], seconds = res$seconds[[i]],
             exit = res$exit[[i]],
             stop = list(reason = res$stop_reason[[i]], class = res$class[[i]],
                         step = res$step[[i]], message_sha256 = res$message_sha256[[i]]),
             options_at_start = valOneLineJson(v$options_at_start), fallback = v$fallback,
             qvalue_fallback = v$qvalue_fallback,
             ruv = list(k = res$ruv_k[[i]], controls = res$ruv_controls[[i]]),
             ledger_sha256 = if (file.exists(led)) valSha256(led) else NA_character_,
             tables = lapply(names(tabs), function(nm) {
                 sha <- if (file.exists(tabs[[nm]])) valSha256(tabs[[nm]]) else NA_character_
                 list(name = nm, sha256 = sha)
             }),
             # TODO(KIT-159 audit I07, G5): the declared `replay` cells are not replayed yet;
             # replayWorkflowAnalysisRecord needs the project inputs G5 records.
             replay = list(status = "not_built"), hook = hook, vmhwm_mib = res$vmhwm_mib[[i]])
    })
}

#' Refuse any file a job wrote under `val-out/` outside `runs/` (5.2 "on val-out/", I08)
#' @noRd
valCheckOutsideRuns <- function(out_dir) {
    files <- list.files(out_dir, recursive = TRUE, all.files = TRUE)
    stray <- files[!startsWith(files, "runs/")]
    if (length(stray)) {
        valAbort("evidence: {stray[[1]]} was written outside runs/", "val_publish_forbidden")
    }
    invisible(TRUE)
}

#' The row bound of metrics.tsv: per metric, comparisons' arms x views x units
#' @noRd
valMetricsRowBound <- function(decl, cells) {
    n_units <- length(unique(cells$unit))
    sum(vapply(decl[["metrics"]] %||% list(), function(m) {
        length(unlist(m[["views"]])) * sum(lengths(lapply(m[["comparisons"]], unlist))) * n_units
    }, numeric(1)))
}

#' The chunk backstop (I17): the chunk's cells x `cell_minutes` + 10
#' @noRd
valChunkMinutes <- function(n_runs, chunks, cell_minutes) {
    ceiling(n_runs / min(chunks, n_runs)) * cell_minutes + 10
}

#' The provenance `code` block of a run made outside `val-run.R` (tests): every field present
#' @noRd
valNoCode <- function(stage) {
    list(package_commit = "none", release_commit = "none", product_tree_identical = NA,
         harness = list(commit = "none", tree = "none"),
         declaration = list(stage = stage, commit = "none", sha256_md = "none",
                            sha256_json = "none"),
         helpers_extracted = list())
}
