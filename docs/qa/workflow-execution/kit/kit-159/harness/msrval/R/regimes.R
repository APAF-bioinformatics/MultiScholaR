# Ticket coverage regimes into the classifier, and phase score (KIT-159 4.4, 5.7, 7.8; audit I03,
# I05, I06).

#' Run a declaration's `coverage_regimes` (4.4; `val-tests.R --ticket --declaration`)
#'
#' Each entry's `generator` (a file of the stage folder) defines `simulate(nsim, seed, levels)`
#' returning one row per simulation and level with `level`, `covered`, `miss_low` and `miss_high`
#' (the truth below the lower or above the upper bound; C66), and optionally `cell`, `achieved` and
#' `critical`; an attribute `cells` (`cell`, `parameter`, `value`) gives the generator's values.
#' The pass rule is the harness band of `kind`, and a declared stricter rule
#' (`pass = list(bound, max_miss)`) is applied with it; never a looser one.
#'
#' @param decl The declaration (a list).
#' @param stage_dir The stage folder the generator paths are relative to.
#' @return `NULL` without regimes, else a list: `coverage` (the columns of `ticket-coverage.tsv`)
#'   and `cells` (`ticket-coverage-cells.tsv`).
#' @export
valTicketCoverage <- function(decl, stage_dir) {
    regs <- decl[["coverage_regimes"]] %||% list()
    if (!length(regs)) return(NULL)
    cov <- list()
    cells <- list()
    for (r in regs) {
        env <- new.env(parent = globalenv())
        sys.source(file.path(stage_dir, r[["generator"]]), envir = env)
        if (!is.function(env$simulate) || is.null(r[["seed"]])) {
            valAbort("coverage regime {r$id}: no simulate() in its generator, or no seed",
                     "val_spec_invalid")
        }
        spec <- list(id = r[["id"]], method = r[["method"]], class = r[["class"]] %||% NA,
                     form = r[["form"]] %||% NA, g = r[["g"]],
                     levels = as.numeric(unlist(r[["levels"]])),
                     kind = r[["kind"]] %||% "approx_run", pass = r[["pass"]] %||% "harness",
                     simulate = env$simulate)
        got <- valTicketRegimeTable(spec, as.integer(r[["simulations"]]), as.integer(r[["seed"]]))
        cov[[length(cov) + 1L]] <- got$coverage
        cells[[length(cells) + 1L]] <- got$cells
    }
    list(coverage = do.call(rbind, cov), cells = do.call(rbind, cells))
}

#' One ticket regime's rows of `ticket-coverage.tsv`: per cell and level the two-sided coverage
#' and each bound's one-sided miss, the declared stricter rule applied with the harness band
#' @noRd
valTicketRegimeTable <- function(regime, nsim, seed) {
    checkmate::assert_function(regime$simulate)
    if (nsim < 2000L) {
        valAbort("a ticket regime needs at least 2,000 simulations", "val_spec_invalid")
    }
    strict <- regime[["pass"]] %||% "harness"
    if (!identical(strict, "harness") &&
            !(is.list(strict) && isTRUE(strict[["bound"]] %in% c("two_sided", "lower", "upper")) &&
                  is.numeric(strict[["max_miss"]]))) {
        valAbort("regime {regime$id}: pass is harness or list(bound, max_miss)", "val_spec_invalid")
    }
    res <- regime$simulate(nsim, seed, regime$levels)
    cells_attr <- attr(res, "cells")
    checkmate::assert_names(names(res), must.include = c("level", "covered", "miss_low",
                                                         "miss_high"))
    if (is.null(res[["cell"]])) res$cell <- "all"
    kind <- regime[["kind"]] %||% "approx_run"
    out <- list()
    for (cell in unique(res$cell)) {
        for (lv in regime$levels) {
            r <- res[res$cell == cell & abs(res$level - lv) < 5e-5, , drop = FALSE]
            n <- nrow(r)
            if (n < 2000L) {
                valAbort("regime {regime$id}: fewer than 2,000 simulations at a level",
                         "val_spec_invalid")
            }
            cov <- mean(as.logical(r$covered))
            lo <- mean(as.logical(r$miss_low))
            hi <- mean(as.logical(r$miss_high))
            ach <- if (is.null(r[["achieved"]])) NA_real_ else mean(r[["achieved"]])
            crit <- if (is.null(r[["critical"]])) NA_real_ else mean(r[["critical"]])
            ok <- valCoverageStatus(kind, cov, lv, ach, n) %in% c("pass", "conservative")
            if (is.list(strict)) {
                miss <- switch(strict[["bound"]], two_sided = 1 - cov, lower = lo, upper = hi)
                ok <- ok && miss <= strict[["max_miss"]]
            }
            se <- function(p) sqrt(p * (1 - p) / n)
            out[[length(out) + 1L]] <- data.frame(
                regime = regime$id, cell = cell, method = regime$method,
                class = as.character(regime$class), form = as.character(regime$form),
                g = as.integer(regime$g), level = lv, simulations = n,
                bound = c("two_sided", "lower", "upper"), coverage = c(cov, 1 - lo, 1 - hi),
                mc_se = c(se(cov), se(lo), se(hi)), miss = c(1 - cov, lo, hi), critical = crit,
                pass = ok, stringsAsFactors = FALSE)
        }
    }
    cells_tab <- cells_attr %||% data.frame(cell = "all", parameter = "simulations", value = nsim)
    list(coverage = do.call(rbind, out)[, EVIDENCE_COLUMNS[["ticket-coverage.tsv"]]],
         cells = data.frame(regime = regime$id, cell = cells_tab$cell,
                            parameter = cells_tab$parameter, value = cells_tab$value,
                            stringsAsFactors = FALSE))
}

#' Ticket regimes as classifier rows: one per regime and level, `pass` only when every cell
#' passes, coverage the worst cell's
#' @noRd
valTicketRegimeRows <- function(tc) {
    two <- tc[tc$bound == "two_sided", , drop = FALSE]
    keys <- unique(two[, c("regime", "method", "class", "form", "g", "level")])
    do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
        k <- keys[i, ]
        r <- two[two$regime == k$regime & abs(two$level - k$level) < 5e-5, , drop = FALSE]
        data.frame(regime = k$regime, method = k$method, class = as.character(k$class),
                   form = as.character(k$form), n = as.integer(k$g), level = k$level,
                   status = if (all(as.logical(r$pass) %in% TRUE)) "pass" else "fail",
                   coverage = min(r$coverage), mc_se = max(r$mc_se), nsim = min(r$simulations),
                   source = "ticket", statistic = paste("ticket regime", k$regime),
                   coverage_min = min(r$coverage), coverage_max = max(r$coverage),
                   stringsAsFactors = FALSE)
    }))
}

#' Bind two regime tables over the union of their columns
#' @noRd
valBindRegimes <- function(a, b) {
    if (is.null(b) || !nrow(b)) return(a)
    cols <- union(names(a), names(b))
    for (cl in setdiff(cols, names(a))) a[[cl]] <- rep(NA, nrow(a))
    for (cl in setdiff(cols, names(b))) b[[cl]] <- rep(NA, nrow(b))
    rbind(a[, cols], b[, cols])
}

#' The classifier's regimes: the harness's table and the committed ticket regimes (4.4, I03)
#'
#' The declaration pins the committed `ticket-coverage.tsv` by hash (`ticket_coverage`:
#' `path`, `sha256`); a declaration with `coverage_regimes` and no pin, a missing file or another
#' hash is refused.
#'
#' @param decl The declaration (a list).
#' @param evidence_dir The ticket's evidence folder.
#' @return The regime table, with attribute `ticket_coverage` (the committed rows) when present.
#' @export
valClassifierRegimes <- function(decl, evidence_dir) {
    harness <- valHarnessRegimes()
    if (!length(decl[["coverage_regimes"]])) return(harness)
    pin <- decl[["ticket_coverage"]]
    if (is.null(pin[["sha256"]])) {
        valAbort("coverage_regimes are declared but ticket-coverage.tsv is not pinned by hash",
                 "val_spec_invalid")
    }
    path <- file.path(evidence_dir, pin[["path"]] %||% "ticket-coverage.tsv")
    if (!file.exists(path)) valAbort("the pinned ticket-coverage.tsv is absent", "val_spec_invalid")
    if (!identical(valSha256(path), pin[["sha256"]])) {
        valAbort("ticket-coverage.tsv differs from its declared hash", "val_declaration_changed")
    }
    tc <- utils::read.delim(path, stringsAsFactors = FALSE)
    out <- valBindRegimes(harness, valTicketRegimeRows(tc))
    attr(out, "ticket_coverage") <- tc
    out
}

#' Phase score: the committed runs of several requests scored at run level (7.8; I06)
#'
#' Reads each request's verified folder ([valScoreRuns()]): its metric rows, its cells' exits (for
#' the refusal rule) and its within-run rows (for the widening); classifies with the given regimes
#' (harness and ticket); writes `comparisons.tsv`, `attribution.tsv`, `results.md` (outside the
#' lifecycle), provenance (`combined_runs`) and `outputs.sha256` into `runs/<label>/`.
#'
#' @param decl The declaration.
#' @param dataset Data set id.
#' @param units The unit table.
#' @param evidence_dir The evidence folder holding `runs/<request>/`.
#' @param labels The requests' labels.
#' @param label This run's label.
#' @param out_dir `val-out/`.
#' @param manifest The input manifest.
#' @param stage The declaration's stage.
#' @param code The provenance `code` block.
#' @param regimes The classifier's regimes.
#' @param started Start time.
#' @return The run folder, invisibly.
#' @export
valScorePhase <- function(decl, dataset, units, evidence_dir, labels, label, out_dir, manifest,
                          stage = "counted", code = NULL, regimes = NULL,
                          started = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) {
    units <- valUnits(units)
    valCheckDeclarationSchema(decl, units)
    code <- code %||% valNoCode(stage)
    code$code_files <- code$code_files %||% list()
    cells <- valPhaseCells(decl, dataset, units, "counted")
    got <- valScoreRuns(evidence_dir, labels, cells[, c("unit", "arm")])
    metrics <- got$metrics
    metrics$value <- as.numeric(metrics$value)
    out <- file.path(out_dir, "runs", label)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    cmp <- valRunLevel(decl, metrics, units, regimes = regimes, phase = stage,
                       within = got$within, cells = got$cells[, c("unit", "arm", "exit")])
    verdicts <- attr(cmp, "verdicts")
    utils::write.table(cmp, file.path(out, "comparisons.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
    att <- manifest[, c("role", "kind", "licence", "redistribution", "citation")]
    utils::write.table(att, file.path(out, "attribution.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    if (!identical(stage, "lifecycle") && !is.null(verdicts)) {
        valWriteResults(file.path(out, "results.md"), verdicts, att)
    }
    prov <- valProvenance(
        ticket = decl[["ticket"]] %||% "unknown", stage = stage, label = label, phase = "score",
        status = "completed", blinded = identical(stage, "lifecycle"), started = started,
        code = code,
        scoring = list(
            comparisons_sha256 = valSha256(file.path(out, "comparisons.tsv")),
            combined_runs = lapply(labels, function(l) {
                list(label = l, outputs_sha256 = valSha256(file.path(evidence_dir, "runs", l,
                                                                     "outputs.sha256")))
            }),
            intervals = verdicts %||% list(),
            ticket_coverage = attr(regimes, "ticket_coverage") %||% list()))
    valWriteProvenance(prov, file.path(out, "provenance.json"))
    valCheckProvenance(file.path(out, "provenance.json"))
    valCheckEvidence(out, manifest, phase = stage)
    valWriteOutputsSha(out)
    invisible(out)
}
