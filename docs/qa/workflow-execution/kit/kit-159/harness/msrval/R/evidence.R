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
#' @return A list: `metrics`, `cells` and the within-run comparison rows (`within`) combined,
#'   `package_commit`.
#' @export
valScoreRuns <- function(dir, labels, expected_cells) {
    changed <- function(what) valAbort("phase score: {what}", "val_declaration_changed")
    provs <- lapply(labels, function(l) {
        d <- file.path(dir, "runs", l)
        if (!valVerifyOutputsSha(d)) changed(paste(l, "does not match its outputs.sha256"))
        jsonlite::fromJSON(file.path(d, "provenance.json"), simplifyVector = FALSE)
    })
    ids <- vapply(provs, function(p) {
        paste(p[["code"]][["package_commit"]] %||% "", p[["code"]][["harness"]][["tree"]] %||% "",
              p[["code"]][["declaration"]][["commit"]] %||% "")
    }, character(1))
    if (length(unique(ids)) != 1L) changed("the requests differ in commit, harness or declaration")
    read <- function(l, f) {
        utils::read.delim(file.path(dir, "runs", l, f), colClasses = "character",
                          na.strings = character(), stringsAsFactors = FALSE)
    }
    cells <- do.call(rbind, lapply(labels, function(l) cbind(read(l, "cells.tsv"), request = l)))
    cmp <- do.call(rbind, lapply(labels, function(l) {
        cbind(read(l, "comparisons.tsv"), request = l)
    }))
    within <- cmp[cmp$scope == "within_run", , drop = FALSE]
    for (col in c("diff", "ci_low", "ci_high", "conf", "margin")) {
        within[[col]] <- suppressWarnings(as.numeric(within[[col]]))
    }
    metrics <- do.call(rbind, lapply(labels, function(l) {
        cbind(read(l, "metrics.tsv"), request = l)
    }))
    key_c <- paste(cells$unit, cells$arm)
    key_e <- paste(expected_cells$unit, expected_cells$arm)
    if (!setequal(key_c, key_e) || anyDuplicated(key_c)) {
        valAbort("phase score: the cells are not each declared cell exactly once",
                 "val_spec_invalid")
    }
    list(metrics = metrics, cells = cells, within = if (nrow(within)) within,
         package_commit = provs[[1L]]$code$package_commit)
}

#' Widen a run-level interval in phase score to the requests' within-run intervals (4.4, A38;
#' revision 1.8 item 9, T-SCORE-2)
#'
#' Below 6 units the widening is mandatory; at 6 or more it is the default and may be declared off
#' (`widen = FALSE`), the within-run interval then reported beside the verdict.
#'
#' @param interval The run-level interval.
#' @param comparisons The requests' `comparisons.tsv` rows combined.
#' @param arm,comparator,metric The comparison.
#' @param widen The declaration's `widen_to_within_run` (default `TRUE`).
#' @param view Only rows of this view (`NULL`: any).
#' @param unit Only the rows of this unit (`unit=<id>;` in `flags`): the jackknife's `full` unit.
#' @return The interval, widened bound by bound or with `within_run_reported`.
#' @export
valWidenFromRequests <- function(interval, comparisons, arm, comparator, metric, widen = NULL,
                                 view = NULL, unit = NULL) {
    sel <- comparisons[["scope"]] == "within_run" & comparisons[["arm"]] == arm &
        comparisons[["comparator"]] == comparator & comparisons[["metric"]] == metric
    if (!is.null(view) && !is.null(comparisons[["view"]])) {
        sel <- sel & comparisons[["view"]] == view
    }
    if (!is.null(unit)) {
        sel <- sel & grepl(paste0("unit=", unit, ";"), comparisons[["flags"]], fixed = TRUE)
    }
    rows <- comparisons[sel %in% TRUE, , drop = FALSE]
    if (nrow(rows) == 0L) return(interval)
    lo <- min(as.numeric(rows[["ci_low"]]))
    hi <- max(as.numeric(rows[["ci_high"]]))
    n <- interval[["n_units"]]
    if (n < 6L && isFALSE(widen)) {
        valAbort("widening to the within-run interval is mandatory below 6 units",
                 "val_spec_invalid")
    }
    if (n < 6L || !isFALSE(widen)) {
        interval[["ci_low"]] <- min(interval[["ci_low"]], lo)
        interval[["ci_high"]] <- max(interval[["ci_high"]], hi)
        interval[["flags"]] <- c(interval[["flags"]], "widened_to_within_run")
    } else {
        interval[["within_run_reported"]] <- c(lo, hi)
        interval[["flags"]] <- c(interval[["flags"]],
                                 sprintf("within_run_reported=[%s,%s]", signif(lo, 6),
                                         signif(hi, 6)))
    }
    interval
}

HARNESS_REGIME_NSIM <- 2000L

#' Write `results.md`: every verdict with its regime's measured coverage and Monte Carlo SE, each
#' primary family's nominal level beside the bound from measured misses, and the attribution
#' (sections 4.4, 5.2; I05; O1 and O5 of VAL-109)
#'
#' The family bound sums the measured miss rates (1 - coverage) of the family's classified
#' components; a regime with no observed miss enters at 3 / simulations (the rule of three), and
#' the note says so. A verdict whose regime's worst-condition miss rate is at least twice nominal is
#' flagged.
#'
#' @param path Output path.
#' @param comparisons Verdict rows with `member`, `metric`, `arm`, `comparator`, `outcome`, `conf`,
#'   `regime`, `coverage`, `mc_se`, `statistic`, `family`, `m` (and optionally `nsim`,
#'   `regime_status`, `miss_twice_nominal`).
#' @param attribution `attribution.tsv` rows.
#' @return `path`, invisibly.
#' @export
valWriteResults <- function(path, comparisons, attribution) {
    cmp <- comparisons
    lines <- c("# Results", "")
    fams <- unique(cmp[["family"]][!is.na(cmp[["family"]])])
    if (length(fams) > 0L) lines <- c(lines, "## Primary families", "")
    for (f in fams) {
        r <- cmp[cmp[["family"]] %in% f, , drop = FALSE]
        m <- r[["m"]][[1L]]
        conf <- r[["conf"]][is.finite(r[["conf"]])][1L]
        nominal <- 1 - m * (1 - conf)
        k <- r[is.finite(r[["coverage"]]), , drop = FALSE]
        nsim <- if (is.null(k[["nsim"]])) rep(NA_real_, nrow(k)) else as.numeric(k[["nsim"]])
        nsim[!is.finite(nsim)] <- HARNESS_REGIME_NSIM
        miss <- 1 - k[["coverage"]]
        zero <- miss <= 0
        miss[zero] <- 3 / nsim[zero]
        zero_note <- ""
        if (any(zero)) {
            zero_note <- sprintf(
                " (%d component(s) with no observed miss entered at 3/simulations)", sum(zero)
            )
        }
        lines <- c(lines, sprintf(
            "- %s: m = %d, nominal family level %s, family error bound from measured misses %s%s",
            f, as.integer(m), format(round(nominal, 6)), format(round(sum(miss), 4)), zero_note
        ))
    }
    lines <- c(lines, "", "## Verdicts", "")
    for (i in seq_len(nrow(cmp))) {
        r <- cmp[i, ]
        twice <- !is.null(r[["miss_twice_nominal"]]) && isTRUE(r[["miss_twice_nominal"]])
        st <- if (is.null(r[["regime_status"]])) NA_character_ else r[["regime_status"]]
        st2 <- if (is.null(r[["regime_status_r113"]])) NA_character_ else r[["regime_status_r113"]]
        # O5 (VAL-109): both labels where T2's band and r1.13's miss-rate rule differ
        label <- ""
        if (!is.na(st) && !is.na(st2) && st != st2) {
            label <- sprintf(" [%s under T2's band; %s under r1.13's miss-rate rule]", st, st2)
        } else if (!is.na(st)) {
            label <- sprintf(" [%s]", st)
        }
        note <- if (is.null(r[["within_note"]])) NA_character_ else r[["within_note"]]
        note_txt <- ""
        if (!is.na(note) && grepl("6_or_more_open", note, fixed = TRUE)) {
            note_txt <- paste("; within-run intervals reported beside the verdict, not widened:",
                              "at 6 or more independent units the widening rule is open (O8)")
        } else if (!is.na(note) && grepl("below_6", note, fixed = TRUE)) {
            note_txt <- paste("; within-run intervals reported beside the verdict, not widened:",
                              "independent units below 6 (O8, VAL-109)")
        }
        lines <- c(lines, sprintf(
            paste("- %s (%s: %s - %s): %s; regime %s%s, measured coverage %s (Monte Carlo SE %s)",
                  "at level %s, simulated on %s%s%s"),
            r[["member"]], r[["metric"]], r[["arm"]], r[["comparator"]], r[["outcome"]],
            r[["regime"]], label, format(r[["coverage"]]), format(r[["mc_se"]]),
            format(r[["conf"]]), r[["statistic"]],
            if (twice) "; FLAG: worst-condition miss rate at least twice nominal" else "",
            note_txt
        ))
    }
    lines <- c(lines, "", "## Sources and attribution", "")
    for (i in seq_len(nrow(attribution))) {
        a <- attribution[i, ]
        lines <- c(lines, sprintf("- %s (%s): %s; %s; %s", a[["role"]], a[["kind"]],
                                  a[["licence"]], a[["redistribution"]], a[["citation"]]))
    }
    writeLines(lines, path)
    invisible(path)
}

#' Write `handoff.sha256` over every file of a hand-off folder
#'
#' @param root The hand-off folder.
#' @return Its path, invisibly.
#' @export
valWriteHandoffSha <- function(root) {
    files <- sort(setdiff(list.files(root, recursive = TRUE), "handoff.sha256"))
    path <- file.path(root, "handoff.sha256")
    writeLines(sprintf("%s  %s", vapply(file.path(root, files), valSha256, character(1)), files),
               path)
    invisible(path)
}
