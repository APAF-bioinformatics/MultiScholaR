# KIT-262 E5b (ticket 9.3): rules-k262.R on constructed traces and curves, the winners computed by hand
# (e5b-hand-k262.py works the same numbers independently of the package). Evidence code: no package change.
#   Rscript docs/qa/workflow-execution/kit/kit-262/e5b-k262.R <out dir relative to the package root>
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(file.path(getwd(), args[[1L]]), mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("docs/qa/workflow-execution/kit/kit-262/rules-k262.R")

longCurve <- function(delta) {
  data.frame(featureset = rep(c("All", "Control"), each = length(delta)), K = rep(seq_along(delta), 2L),
    cc = c(0.5 + delta, rep(0.5, length(delta))), stringsAsFactors = FALSE)
}
rowOf <- function(pct, delta, status, k, composite) data.frame(percentage_requested = pct, status = status, best_k = as.integer(k),
  separation_score = max(delta), composite_score = composite, error_reason = "", stringsAsFactors = FALSE)
checks <- data.frame(check = character(), pass = logical(), detail = character(), stringsAsFactors = FALSE)
add <- function(check, pass, detail = "") checks[nrow(checks) + 1L, ] <<- list(check, isTRUE(pass), as.character(detail))
same <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-3))
winnerIs <- function(rules, rule, pct, k) {
  w <- rules$winner[[rule]]
  identical(w$percentage, as.numeric(pct)) && identical(w$k_used, as.integer(k))
}

# A: a re-admitted percentage outranks P's winner under C1 (P02's example: k = 6 scores 0.0055, k = 3 with c = 2 scores 0.0552)
dA5 <- c(.1, .2, .3, .4, .5, .6, .6, .6); dA10 <- c(.1, .2, .3, .3, .3, .3, .3, .3)
tA <- rbind(rowOf(5, dA5, "ok", 6, 0.005494692), rowOf(10, dA10, K262_R5, 3, NA_real_))
rA <- k262Rules(tA, list(longCurve(dA5), longCurve(dA10)), c(8L, 2L), 2L)
add("A_composites", same(rA$rows$composite_recomputed, c(0.005495, 0.055182)), paste(rA$rows$composite_recomputed, collapse = " "))
add("A_P_winner_5_k2", winnerIs(rA, "P", 5, 2))
add("A_C1_winner_10_k2", winnerIs(rA, "C1", 10, 2))
add("A_C2_winner_5_k2_by_tie_rule", winnerIs(rA, "C2", 5, 2) && same(rA$rows$composite_b, c(0.1, 0.1)))
add("A_usable_sets", identical(rA$rows$p_usable, c(TRUE, FALSE)) && identical(rA$rows$c1_usable, c(TRUE, TRUE)) && identical(rA$rows$c2_usable, c(TRUE, TRUE)))
# B: a re-admitted percentage that does not outrank
dB5 <- c(.3, .6, .9, .9, .9, .9); dB10 <- c(0, .05, .1, .15, .2, .2)
tB <- rbind(rowOf(5, dB5, "ok", 3, 0.165546), rowOf(10, dB10, K262_R5, 5, NA_real_))
rB <- k262Rules(tB, list(longCurve(dB5), longCurve(dB10)), c(8L, 2L), 2L)
add("B_all_winners_5", winnerIs(rB, "P", 5, 2) && winnerIs(rB, "C1", 5, 2) && winnerIs(rB, "C2", 5, 2))
add("B_C1_readmits_10", identical(rB$rows$c1_usable, c(TRUE, TRUE)) && same(rB$rows$composite_recomputed[[2L]], 0.004979))
# C: a curve whose truncated elbow is below min(k, cap) (cap 3: k = 4, k_b = 1)
dC <- c(.5, .5, .5, 1, 1, 1)
tC <- rowOf(5, dC, "ok", 4, 0.18394)
rC <- k262Rules(tC, list(longCurve(dC)), 6L, 3L)
add("C_kb_below_min_k_cap", identical(rC$rows$k_b, 1L) && winnerIs(rC, "P", 5, 3) && winnerIs(rC, "C1", 5, 3) && winnerIs(rC, "C2", 5, 1))
# D: a truncated separation <= 0 (composite 0, usable): k_b = 1 (a weak-effect curve), separation_b = -0.05
dD <- c(-.1, -.05, .4, .4, .4, .4)
tD <- rowOf(5, dD, "ok", 3, 0.073576)
rD <- k262Rules(tD, list(longCurve(dD)), 6L, 2L)
add("D_nonpositive_truncated_separation_usable", same(rD$rows$separation_b, -0.05) && same(rD$rows$composite_b, 0) && rD$rows$c2_usable &&
  winnerIs(rD, "C2", 5, 1))
# E: an unalignable control set (c = NA) passes R5 in every rule; an R5-status row with NA c is admitted by C1 and C2
tE <- rbind(rowOf(5, dA10, "ok", 3, 0.055182), rowOf(10, dA10, K262_R5, 3, NA_real_))
rE <- k262Rules(tE, list(longCurve(dA10), longCurve(dA10)), c(NA_integer_, NA_integer_), 2L)
add("E_NA_complete_passes", identical(rE$rows$p_usable, c(TRUE, FALSE)) && identical(rE$rows$c1_usable, c(TRUE, TRUE)) && identical(rE$rows$c2_usable, c(TRUE, TRUE)))
# the cause split: an R5 row with c below min(k, cap) is completeness, with c at or above it cap_side
rows_cause <- data.frame(status = c(K262_R5, K262_R5, "ok"), k = c(5L, 5L, 2L), c = c(1L, 3L, NA_integer_), stringsAsFactors = FALSE)
cs <- k262Cause(rows_cause, 3L)
add("cause_split", identical(cs, c("completeness", "cap_side", NA_character_)), paste(cs, collapse = " "))
# the reason parser
p <- k262ParseReason("k = 5 exceeds 2 controls complete in every sample")
add("reason_parser", p$ok && p$k == 5L && p$c == 2L && !k262ParseReason("something else")$ok && is.na(k262ParseReason("k = 3 exceeds NA controls complete in every sample")$c))
utils::write.table(checks, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
quit(status = if (all(checks$pass)) 0L else 1L)
