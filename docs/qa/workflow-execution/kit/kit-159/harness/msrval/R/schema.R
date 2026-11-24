# The declaration's statistical schema (KIT-159 4.3 to 4.7, 5.3; audit I11, I16, I20).

#' Check a declaration's statistical schema before any cell runs
#'
#' Refuses (`val_spec_invalid`): more primaries than data sets, a primary metric with more than
#' one contrast or comparison, a `primary_family` with m above 4 or members off the level
#' 1 - 0.05/m; a compared metric without `margin` or `contrast`; `feature_interval_decides`
#' without `feature_interval_reason`; `accept_flagged` without `label` and `ledger_row`;
#' `widen_to_within_run = false` below 6 units; a `worse_when` group naming a guard or an
#' undeclared component; an interval method that does not fit the unit roles (`wilson` and
#' `newcombe_paired` are never run-level methods; `clopper_pearson` is one-arm; `jackknife` needs
#' `full` and `deletion` units); a within-run B below 50 / tail or without a seed; a data set
#' without a seed, or no `blind.seed`; `valPartialAuc` without `fpr_max`; class `accuracy` in the
#' union view; an unknown `metric_input`.
#'
#' @param decl The declaration (a list).
#' @param units The unit table (data frame), the unit ids, or `NULL`.
#' @return `TRUE`, invisibly.
#' @export
valCheckDeclarationSchema <- function(decl, units = NULL) {
    bad <- function(what) valAbort("declaration: {what}", "val_spec_invalid")
    n_units <- if (is.data.frame(units)) nrow(units) else if (is.character(units)) {
        length(units)
    } else {
        NA_integer_
    }
    roles <- if (is.data.frame(units)) units$unit_role else NULL
    for (ds in decl[["datasets"]] %||% list()) {
        if (is.null(ds[["seed"]])) bad(paste("data set", ds[["id"]], "declares no seed"))
    }
    if (length(decl[["datasets"]]) && is.null(decl[["blind"]][["seed"]])) bad("no blind.seed")
    metrics <- decl[["metrics"]] %||% list()
    components <- character()
    primaries <- 0L
    for (m in metrics) {
        id <- m[["id"]] %||% "?"
        cmps <- m[["comparisons"]] %||% list()
        method <- (m[["interval"]] %||% list())[["method"]] %||% "paired_t"
        for (cmp in cmps) {
            components <- c(components, paste0(id, ":", paste(unlist(cmp), collapse = "-")))
        }
        if (isTRUE(m[["primary"]])) {
            primaries <- primaries + 1L
            if (length(cmps) != 1L || length(unlist(m[["contrast"]])) != 1L) {
                bad(paste("primary metric", id, "needs one contrast and one comparison"))
            }
        }
        if (!length(cmps)) next
        if (is.null(m[["margin"]]) || is.null(m[["contrast"]])) {
            bad(paste("compared metric", id, "lacks a margin or a contrast"))
        }
        if (isTRUE(m[["feature_interval_decides"]]) &&
                !nzchar(m[["feature_interval_reason"]] %||% "")) {
            bad(paste("metric", id, "sets feature_interval_decides without a reason"))
        }
        af <- (m[["interval"]] %||% list())[["accept_flagged"]]
        if (!is.null(af) &&
                (is.null(af[["label"]]) || !grepl("^VAL-", af[["ledger_row"]] %||% ""))) {
            bad(paste("metric", id, "accept_flagged must name a label and a ledger row"))
        }
        if (isFALSE(m[["widen_to_within_run"]]) && is.finite(n_units) && n_units < 6L) {
            bad(paste("metric", id, "cannot decline the widening below 6 units"))
        }
        if (method %in% c("wilson", "newcombe_paired")) {
            bad(paste("metric", id, "uses", method, "as a run-level method"))
        }
        if (identical(method, "clopper_pearson") && any(lengths(lapply(cmps, unlist)) != 1L)) {
            bad(paste("metric", id, "clopper_pearson is one-arm"))
        }
        if (!is.null(roles)) {
            ok <- switch(method,
                         jackknife = all(c("full", "deletion") %in% roles),
                         bootstrap_within_run = TRUE,
                         any(roles %in% c("independent", "draw")))
            if (!ok) bad(paste("metric", id, "method", method, "does not fit the unit roles"))
        }
        wr <- m[["within_run"]]
        if (!is.null(wr) || identical(method, "bootstrap_within_run")) {
            conf <- m[["conf"]] %||% 0.95
            if (is.null(wr[["seed"]])) bad(paste("metric", id, "within_run declares no seed"))
            if (!is.null(wr[["B"]]) && wr[["B"]] < valBootstrapMinB(conf)) {
                bad(paste("metric", id, "within_run B below 50 / tail"))
            }
        }
        if (identical(m[["function"]], "valPartialAuc") && is.null(m[["args"]][["fpr_max"]])) {
            bad(paste("metric", id, "valPartialAuc needs fpr_max"))
        }
        if (identical(m[["class"]], "accuracy") && "union" %in% unlist(m[["views"]])) {
            bad(paste("metric", id, "accuracy is never scored in the union view"))
        }
        if (!(m[["metric_input"]] %||% "scored") %in% c("scored", "confusion", "matrix")) {
            bad(paste("metric", id, "metric_input must be scored, confusion or matrix"))
        }
    }
    n_ds <- max(1L, length(decl[["datasets"]]))
    if (primaries + length(decl[["primary_family"]]) > n_ds) {
        bad("more than one primary per data set")
    }
    for (f in names(decl[["primary_family"]])) {
        pf <- decl[["primary_family"]][[f]]
        mem <- pf[["members"]] %||% list()
        m_decl <- as.integer(pf[["m"]] %||% length(mem))
        if (m_decl > 4L || m_decl != length(mem)) {
            bad(paste("family", f, "m must be <= 4 and count its members"))
        }
        lvl <- 1 - 0.05 / m_decl
        # Re-audit 2 R2-02: every metric a member names (all_of, worse_when) is held to the level,
        # whether or not it declares primary_family.
        named <- unique(sub(":.*$", "", unlist(lapply(mem, function(x) {
            c(unlist(x[["all_of"]]), unlist(x[["worse_when"]]))
        }))))
        for (m in metrics) {
            in_family <- identical(m[["primary_family"]], f) || isTRUE(m[["id"]] %in% named)
            if (in_family && abs((m[["conf"]] %||% 0.95) - lvl) > 1e-6) {
                bad(paste("family", f, "metric", m[["id"]], "is not at 1 - 0.05/m"))
            }
        }
        for (nm in names(mem)) {
            guards <- vapply(mem[[nm]][["guards"]] %||% list(), function(g) g[["component"]],
                             character(1))
            ww <- unlist(mem[[nm]][["worse_when"]])
            if (any(ww %in% guards)) bad(paste("member", nm, "worse_when names a guard"))
            if (length(components) && !all(c(unlist(mem[[nm]][["all_of"]]), ww) %in% components)) {
                bad(paste("member", nm, "names an undeclared component"))
            }
        }
    }
    invisible(TRUE)
}
