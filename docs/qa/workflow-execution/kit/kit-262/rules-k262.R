# KIT-262 (VAL-280, VAL-302): the script-local rules. No package code. Sourced by m1-k262.R and m2-k262.R after the
# package is loaded (pkgload::load_all, export_all = TRUE). Evidence code only; decides nothing.
#
# Rules (ticket section 3), per assay and percentage, k the elbow of the whole curve, c the controls complete in every
# sample, cap = calculateAdaptiveMaxK(samples):
#   P   (shipped)  usable when status is "ok" (k <= c), ranked by the composite on k, the cap only a penalty.
#   C1             usable when min(k, cap) <= c (R5 judged at the capped k); ranking unchanged (composite on the full k).
#   C2             (proteomics' rule, KIT-123 R1) the curve is cut at K <= cap before the elbow and the separation; R5
#                  and the ranking on that k_b. Reported for coverage and k only (VAL-302, O-1 (ii)).

K262_R5 <- "skipped_k_exceeds_complete_controls"
K262_NOT_REACHED <- c("error_neg_ctrl_selection", "error_cancor", "skipped_insufficient_controls")
K262_REASON <- "^k = ([0-9]+) exceeds (NA|[0-9]+) controls complete in every sample$"

# R5's reason -> list(ok, k, c). c is NA_integer_ when the reason says NA.
k262ParseReason <- function(reason) {
  m <- regmatches(reason, regexec(K262_REASON, reason))[[1L]]
  if (length(m) != 3L) return(list(ok = FALSE, k = NA_integer_, c = NA_integer_))
  list(ok = TRUE, k = as.integer(m[[2L]]), c = if (identical(m[[3L]], "NA")) NA_integer_ else as.integer(m[[3L]]))
}

# The long curve data (featureset, K, cc) of a canonical-correlation plot, or NULL.
k262CurveData <- function(cancorplot) {
  if (inherits(cancorplot, "ggplot")) cancorplot$data else if (is.list(cancorplot)) cancorplot$data else NULL
}

# Runs code() with every call of findBestKElbow in the package's namespace recorded: the curve data of its argument,
# in call order. The recorder never signals (a signal inside the elbow's tryCatch would silently become R4) and
# returns the package function's own value. `trace` is not used; testthat's mocking is the test file's own method.
k262WithCurveCapture <- function(code) {
  store <- new.env(parent = emptyenv())
  store$curves <- list()
  store$record_failures <- 0L
  original <- get("findBestKElbow", envir = asNamespace("MultiScholaR"))
  recorder <- function(cancorplot, ...) {
    tryCatch({
      store$curves[length(store$curves) + 1L] <- list(k262CurveData(cancorplot))
    }, error = function(e) {
      store$record_failures <- store$record_failures + 1L
      store$curves[length(store$curves) + 1L] <- list(NULL)
    })
    original(cancorplot, ...)
  }
  value <- testthat::with_mocked_bindings(code(), findBestKElbow = recorder, .package = "MultiScholaR")
  list(value = value, curves = store$curves, record_failures = store$record_failures)
}

# Assigns the recorded curves to the trace rows of each assay: one capture per row whose status is not R1, R2 or R3
# (R4 rows included), assays in the result's order, rows in percentage order (the sweep's imap over assays, map over
# percentages). Returns a named list (per assay) of lists aligned with the trace rows (NULL where none), and whether the
# number of captures matched.
k262AlignCurves <- function(result, curves) {
  position <- 0L
  aligned <- list()
  for (assay in names(result)) {
    trace <- result[[assay]]$optimization_results
    rows <- vector("list", if (is.data.frame(trace)) nrow(trace) else 0L)
    if (length(rows)) for (i in seq_along(rows)) {
      if (!trace$status[[i]] %in% K262_NOT_REACHED) {
        position <- position + 1L
        if (position <= length(curves)) rows[i] <- list(curves[[position]])
      }
    }
    aligned[[assay]] <- rows
  }
  list(curves = aligned, matched = identical(position, length(curves)))
}

# The curve cut at K <= cap, as proteomics does (R/func_prot_norm_optimization_helpers.R:905-908).
k262Bounded <- function(data, cap) {
  if (is.data.frame(data) && "K" %in% names(data)) data[data$K <= cap, , drop = FALSE] else data
}

k262Winner <- function(percentage, composite, k, separation, usable) {
  idx <- which(usable & is.finite(composite))
  if (!length(idx)) return(NA_integer_)
  idx[order(-composite[idx], k[idx], -separation[idx], percentage[idx])[[1L]]]
}

# The three rules on one assay's trace.
#   trace     the assay's optimization_results (data.frame)
#   curves    list aligned with the trace rows (curve data or NULL); NULL when no curve is available (then C2 is NA)
#   complete  integer vector, controls complete in every sample per row; NA where not known (rows with status "ok"
#             need none: k <= c holds there)
#   cap       calculateAdaptiveMaxK(samples)
# Returns list(rows = data.frame, winner = list(P, C1, C2) each list(row, percentage, k_used)).
k262Rules <- function(trace, curves, complete, cap, k_penalty_weight = 0.5, metric = "max_difference") {
  n <- nrow(trace)
  reached <- !trace$status %in% K262_NOT_REACHED
  k <- trace$best_k
  sep <- trace$separation_score
  is_ok <- trace$status == "ok"
  is_r5 <- trace$status == K262_R5
  composite <- vapply(seq_len(n), function(i) {
    if (!reached[[i]] || is.na(k[[i]])) NA_real_ else
      tryCatch(calculateCompositeScore(sep[[i]], k[[i]], k_penalty_weight, cap), error = function(e) NA_real_)
  }, numeric(1L))
  # P: the shipped status, and the composite recomputed from the trace's own separation and k
  p_usable <- is_ok & is.finite(trace$composite_score)
  # C1: R5 judged at min(k, cap); a row with k and a finite composite that passes it
  k_cap <- ifelse(is.na(k), NA_integer_, pmin(k, cap))
  c1_pass <- reached & !is.na(k) & (is_ok | !vapply(seq_len(n), function(i) isTRUE(k_cap[[i]] > complete[[i]]), logical(1L)))
  c1_usable <- c1_pass & is.finite(composite)
  # C2: the truncated curve
  k_b <- rep(NA_integer_, n)
  sep_b <- rep(NA_real_, n)
  composite_b <- rep(NA_real_, n)
  finite_at_cap <- rep(NA, n)
  if (!is.null(curves)) for (i in seq_len(n)) {
    if (!reached[[i]]) next
    bounded <- list(data = k262Bounded(curves[[i]], cap))
    sep_b[[i]] <- tryCatch(calculateSeparationScore(bounded, metric), error = function(e) NA_real_)
    k_b[[i]] <- tryCatch(as.integer(findBestKElbow(bounded)), error = function(e) NA_integer_)
    finite_at_cap[[i]] <- is.finite(sep_b[[i]])
    composite_b[[i]] <- if (is.na(k_b[[i]])) NA_real_ else
      tryCatch(calculateCompositeScore(sep_b[[i]], k_b[[i]], k_penalty_weight, cap), error = function(e) NA_real_)
  }
  c2_pass <- reached & !is.na(k_b) & !vapply(seq_len(n), function(i) isTRUE(k_b[[i]] > complete[[i]]), logical(1L))
  c2_usable <- c2_pass & is.finite(composite_b)
  rows <- data.frame(percentage = trace$percentage_requested, status = trace$status, k = k, c = complete,
    separation = sep, composite_shipped = trace$composite_score, composite_recomputed = composite,
    p_usable = p_usable, c1_usable = c1_usable, k_b = k_b, separation_b = sep_b, composite_b = composite_b,
    c2_usable = c2_usable, finite_at_cap = finite_at_cap, stringsAsFactors = FALSE)
  pick <- function(usable, comp, kk, ss) {
    idx <- k262Winner(trace$percentage_requested, comp, kk, ss, usable)
    list(row = idx, percentage = if (is.na(idx)) NA_real_ else trace$percentage_requested[[idx]],
      k_used = if (is.na(idx)) NA_integer_ else as.integer(kk[[idx]]))
  }
  # P and C1 are ranked on the shipped composite (the same number); the k they use is min(k, cap) after the search
  winner_p <- pick(p_usable, trace$composite_score, k, sep)
  winner_c1 <- pick(c1_usable, composite, k, sep)
  winner_c2 <- pick(c2_usable, composite_b, k_b, sep_b)
  if (!is.na(winner_p$row)) winner_p$k_used <- as.integer(min(k[[winner_p$row]], cap))
  if (!is.na(winner_c1$row)) winner_c1$k_used <- as.integer(min(k[[winner_c1$row]], cap))
  list(rows = rows, winner = list(P = winner_p, C1 = winner_c1, C2 = winner_c2), cap = cap)
}

# The cause of each R5 rejection: "cap_side" when min(k, cap) <= c (the cap would bring k within the complete controls,
# C1 admits it) and "completeness" when c < min(k, cap) (too few complete controls even at the capped k). NA for rows
# that are not R5.
k262Cause <- function(rows, cap) {
  r5 <- rows$status == K262_R5
  cause <- rep(NA_character_, nrow(rows))
  cause[r5] <- ifelse(is.na(rows$c[r5]) | pmin(rows$k[r5], cap) <= rows$c[r5], "cap_side", "completeness")
  cause
}
