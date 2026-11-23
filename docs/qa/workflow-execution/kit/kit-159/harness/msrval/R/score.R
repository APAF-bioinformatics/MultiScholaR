# The scored table, the matched BH rule, views and blinding (KIT-159 sections 4.1, 4.3, 7.8).

SCORED_COLUMNS <- c(
    "truth_id", "unit", "arm", "assay", "contrast", "key", "key_status", "class", "truth_status",
    "expected", "estimate", "p", "s_own", "s_bh", "called_own", "called_bh", "quantified",
    "testable"
)

#' Build the scored table of one truth and one unit
#'
#' One row per arm x contrast x (assay, key) of the key universe. A universe key absent from an
#' arm's table is a row with no estimate and no p (a miss in the union view). Two or more features
#' of one arm and one assay mapping to one key make that key `duplicate` for the arm. The matched
#' rule `s_bh` is Benjamini-Hochberg of the raw p within arm x unit x contrast x assay over the
#' testable keys (A34).
#'
#' @param arms Named list of arm tables with columns `feature_id`, `key`, `key_status`, `assay`,
#'   `contrast`, `estimate`, `p`, `s_own` (`key_status` and `assay` optional).
#' @param truth Truth table: `key`, `class`, `expected`, `truth_status`.
#' @param universe Key universe: `key` and optionally `assay`.
#' @param unit Unit id.
#' @param truth_id Truth table id.
#' @param alpha Call threshold for `s_own` and `s_bh`.
#' @param lambda Fold-change threshold of the arm's call rule (0: none).
#' @return The scored table (columns of section 4.1); attribute `not_in_universe` counts arm rows
#'   whose key is outside the universe.
#' @export
valScoreTable <- function(arms, truth, universe, unit = "all", truth_id = "primary", alpha = 0.05,
                          lambda = 0) {
    checkmate::assert_list(arms, names = "unique", min.len = 1L)
    checkmate::assert_data_frame(truth)
    checkmate::assert_names(names(truth), must.include = c("key", "class", "expected",
                                                           "truth_status"))
    checkmate::assert_data_frame(universe)
    checkmate::assert_names(names(universe), must.include = "key")
    checkmate::assert_number(alpha, lower = 0, upper = 1)
    checkmate::assert_number(lambda, lower = 0)
    if (anyDuplicated(truth$key)) {
        valAbort("the truth table holds a key more than once", "val_truth_invalid")
    }
    if (is.null(universe$assay)) universe$assay <- "all"
    universe <- unique(data.frame(assay = as.character(universe$assay),
                                  key = as.character(universe$key), stringsAsFactors = FALSE))
    outside <- 0L
    parts <- lapply(names(arms), function(arm) {
        res <- valScoreOneArm(arms[[arm]], arm, universe)
        outside <<- outside + res$outside
        res$rows
    })
    s <- do.call(rbind, parts)
    s$unit <- unit
    s$truth_id <- truth_id
    m <- match(s$key, as.character(truth$key))
    s$class <- as.character(truth$class)[m]
    s$expected <- as.numeric(truth$expected)[m]
    s$truth_status <- ifelse(is.na(m), "excluded", as.character(truth$truth_status)[m])
    ok <- s$key_status == "ok"
    s$quantified <- ok & is.finite(s$estimate)
    s$testable <- ok & is.finite(s$p)
    family <- paste(s$arm, s$unit, s$contrast, s$assay, sep = "\r")
    s$s_bh <- NA_real_
    for (f in unique(family[s$testable])) {
        i <- which(family == f & s$testable)
        s$s_bh[i] <- stats::p.adjust(s$p[i], method = "BH")
    }
    fc_ok <- if (lambda > 0) is.finite(s$estimate) & abs(s$estimate) >= lambda else TRUE
    s$called_own <- s$testable & is.finite(s$s_own) & s$s_own < alpha & fc_ok
    s$called_bh <- s$testable & is.finite(s$s_bh) & s$s_bh < alpha & fc_ok
    s <- s[, SCORED_COLUMNS]
    rownames(s) <- NULL
    attr(s, "not_in_universe") <- outside
    s
}

#' Expand one arm's table over the universe for each of its contrasts
#' @noRd
valScoreOneArm <- function(tab, arm, universe) {
    checkmate::assert_data_frame(tab)
    checkmate::assert_names(names(tab), must.include = c("key", "contrast", "estimate", "p"))
    tab$key <- as.character(tab$key)
    if (is.null(tab$assay)) tab$assay <- "all"
    if (is.null(tab$key_status)) tab$key_status <- "ok"
    if (is.null(tab$s_own)) tab$s_own <- NA_real_
    if (is.null(tab$feature_id)) tab$feature_id <- tab$key
    id <- paste(tab$assay, tab$key, sep = "\r")
    u_id <- paste(universe$assay, universe$key, sep = "\r")
    outside <- sum(!is.na(tab$key) & !(id %in% u_id))
    rows <- lapply(unique(tab$contrast), function(ct) {
        t1 <- tab[tab$contrast == ct & id %in% u_id, , drop = FALSE]
        t1 <- valCollapseDuplicates(t1)
        t1_id <- paste(t1$assay, t1$key, sep = "\r")
        m <- match(u_id, t1_id)
        data.frame(
            arm = arm, assay = universe$assay, contrast = ct, key = universe$key,
            key_status = ifelse(is.na(m), "ok", t1$key_status[m]),
            estimate = as.numeric(t1$estimate)[m], p = as.numeric(t1$p)[m],
            s_own = as.numeric(t1$s_own)[m], stringsAsFactors = FALSE
        )
    })
    list(rows = do.call(rbind, rows), outside = outside)
}

#' Mark keys reached by two or more features of one arm and assay as `duplicate` (one row each)
#' @noRd
valCollapseDuplicates <- function(t1) {
    if (nrow(t1) == 0L) return(t1)
    id <- paste(t1$assay, t1$key, sep = "\r")
    ok_id <- id[t1$key_status == "ok"]
    dup <- unique(ok_id[duplicated(ok_id)])
    if (length(dup) == 0L) return(t1[!duplicated(id), , drop = FALSE])
    t1$key_status[id %in% dup] <- "duplicate"
    t1$estimate[id %in% dup] <- NA_real_
    t1$p[id %in% dup] <- NA_real_
    t1$s_own[id %in% dup] <- NA_real_
    t1[!duplicated(id), , drop = FALSE]
}

#' Restrict a scored table to a view (section 4.3)
#'
#' `own`: the keys each arm quantified (accuracy) or tested (detection), arm by arm.
#' `intersection`: keys every arm in `arms` quantified or tested. `union`: the universe restricted
#' to scored truth keys. A key whose status is not `ok` in any arm of the comparison is excluded
#' for all its arms; attribute `excluded` counts the excluded keys by reason.
#'
#' @param scored A scored table.
#' @param arms Arms of the comparison.
#' @param view One of `own`, `intersection`, `union`.
#' @param metric_type `accuracy` (quantified keys) or `detection` (tested keys).
#' @return The rows of the view.
#' @export
valViews <- function(scored, arms, view = c("own", "intersection", "union"),
                     metric_type = c("accuracy", "detection")) {
    view <- match.arg(view)
    metric_type <- match.arg(metric_type)
    checkmate::assert_character(arms, min.len = 1L)
    if (!all(arms %in% scored$arm)) {
        valAbort("an arm of the view is not in the scored table", "val_spec_invalid")
    }
    s <- scored[scored$arm %in% arms, , drop = FALSE]
    id <- paste(s$contrast, s$assay, s$key, sep = "\r")
    flag <- if (metric_type == "accuracy") s$quantified else s$testable
    bad <- s$key_status != "ok"
    first_bad <- s[bad, c("key_status"), drop = FALSE]
    first_bad$id <- id[bad]
    first_bad <- first_bad[!duplicated(first_bad$id), , drop = FALSE]
    excluded <- if (nrow(first_bad) == 0L) {
        data.frame(reason = character(), n = integer())
    } else {
        tab <- table(first_bad$key_status)
        data.frame(reason = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
    }
    keep <- !(id %in% first_bad$id)
    if (view == "own") {
        keep <- flag & s$key_status == "ok"
    } else if (view == "intersection") {
        all_flag <- tapply(flag, id, function(f) length(f) == length(arms) && all(f))
        keep <- keep & all_flag[id]
    } else {
        keep <- keep & s$truth_status == "scored"
    }
    out <- s[keep, , drop = FALSE]
    attr(out, "excluded") <- excluded
    out
}

#' Blind a truth table for the lifecycle (section 7.8, step 7)
#'
#' The (class, expected) labels are sorted and dealt to the sorted keys in a seeded order, so the
#' result depends only on the key set, the multiset of labels and the seed.
#'
#' @param truth Truth table: `key`, `class`, `expected`, `truth_status`.
#' @param seed The declaration's blind seed.
#' @return The blinded truth table, ordered by key, with `blinded = TRUE`.
#' @export
valBlindTruth <- function(truth, seed) {
    checkmate::assert_data_frame(truth)
    checkmate::assert_names(names(truth), must.include = c("key", "class", "expected",
                                                           "truth_status"))
    checkmate::assert_int(seed)
    if (anyDuplicated(truth$key)) {
        valAbort("the truth table holds a key more than once", "val_truth_invalid")
    }
    keys <- sort(as.character(truth$key), method = "radix")
    lab <- truth[order(as.character(truth$class), as.numeric(truth$expected), method = "radix"),
                 c("class", "expected")]
    ord <- valWithSeed(seed, sample.int(length(keys)))
    class_b <- character(length(keys))
    expected_b <- numeric(length(keys))
    class_b[ord] <- as.character(lab$class)
    expected_b[ord] <- as.numeric(lab$expected)
    m <- match(keys, as.character(truth$key))
    data.frame(key = keys, class = class_b, expected = expected_b,
               truth_status = as.character(truth$truth_status)[m], blinded = TRUE,
               stringsAsFactors = FALSE)
}
