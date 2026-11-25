# KIT-262 M1 (ticket section 4): the rejection reason of each percentage on the cancor-once fixture, lipidomics and
# metabolomics, read from an UNHOOKED sweep; then the same sweep with the curve recorder (E1'), the rules P, C1 and C2
# (rules-k262.R), and the checks E1, E1', E2, E3, E7. Evidence code: no package change.
#
# Usage, from the package root on a run-tests runner:
#   Rscript docs/qa/workflow-execution/kit/kit-262/m1-k262.R <out dir relative to the package root>
suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = FALSE, quiet = TRUE))
logger::log_threshold(logger::FATAL)
args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(file.path(getwd(), args[[1L]]), mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("docs/qa/workflow-execution/kit/kit-262/rules-k262.R")

# ---- copy: tests/testthat/test-workflow-omics-ruv-cancor-once.R:5-28, :31-40 and :43-54, verbatim ----
ruvAssay <- function(id_column, seed, n_features = 300L, n_samples = 24L, missing = 0.05) {
  set.seed(seed)
  values <- matrix(rnorm(n_features * n_samples, 15, 2), nrow = n_features) + rep(rnorm(n_samples, 0, 0.3), each = n_features)
  values[seq_len(30L), seq(1L, n_samples, by = 2L)] <- values[seq_len(30L), seq(1L, n_samples, by = 2L)] + 3
  if (missing > 0) values[sample(length(values), floor(missing * length(values)))] <- NA_real_
  colnames(values) <- sprintf("S%02d", seq_len(n_samples))
  assay <- tibble::as_tibble(values)
  assay[[id_column]] <- sprintf("F%04d", seq_len(n_features))
  assay[c(id_column, colnames(values))]
}

ruvDesign <- function() {
  data.frame(Run = sprintf("S%02d", seq_len(24L)), group = rep(c("A", "B"), length.out = 24L), stringsAsFactors = FALSE)
}

ruvMetabObject <- function(assays, args = list()) {
  methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = ruvDesign(), sample_id = "Run", group_id = "group",
    metabolite_id_column = "feature", technical_replicate_id = "group", args = args)
}

ruvLipidObject <- function(assays, args = list()) {
  methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = ruvDesign(), sample_id = "Run", group_id = "group",
    lipid_id_column = "feature", technical_replicate_id = "group", args = args)
}

ruvCases <- function() {
  assays <- function(missing) list(one = ruvAssay("feature", 11L, missing = missing), two = ruvAssay("feature", 12L, missing = missing))
  list(
    metab = list(make = ruvMetabObject, sweep = MultiScholaR:::runPerAssayRuvOptimization, assays = assays),
    lipid = list(make = ruvLipidObject, sweep = MultiScholaR:::runLipidPerAssayRuvOptimization, assays = assays))
}

ruvSweep <- function(case, object) {
  suppressWarnings(suppressMessages(case$sweep(object, ruv_mode = "automatic",
    params = list(percentage_min = 2, percentage_max = 6, ruv_grouping_variable = "group"))))
}

# A sweep result without its plot objects, plus the data each plot was drawn from.
ruvComparable <- function(result) {
  lapply(result, function(assay) {
    lapply(assay, function(value) {
      if (inherits(value, "ggplot")) return(value$data)
      if (is.list(value) && !is.data.frame(value)) {
        return(lapply(value, function(inner) if (inherits(inner, "ggplot")) inner$data else inner))
      }
      value
    })
  })
}
# ---- end of the copy ----

# Complete controls recomputed for every percentage with the package's own helpers (the sweep's calls).
m1CompleteControls <- function(object, omic, percentages) {
  slots <- if (omic == "lipid") c("lipid_data", "lipid_id_column") else c("metabolite_data", "metabolite_id_column")
  assays <- methods::slot(object, slots[[1L]])
  id_column <- methods::slot(object, slots[[2L]])
  design_samples <- as.character(object@design_matrix[[object@sample_id]])
  result <- list()
  for (pct in percentages) {
    ctrl <- suppressWarnings(suppressMessages(MultiScholaR:::getNegCtrlMetabAnova(object, ruv_grouping_variable = "group",
      percentage_as_neg_ctrl = pct)))
    for (assay in names(assays)) {
      complete_rows <- MultiScholaR:::omicsRuvCompleteRows(assays[[assay]], id_column, design_samples)
      ctrl_assay <- if (is.list(ctrl) && assay %in% names(ctrl)) ctrl[[assay]] else ctrl
      result[[paste(assay, pct)]] <- list(c = MultiScholaR:::omicsRuvCompleteControlCount(complete_rows, ctrl_assay),
        controls = sum(ctrl_assay, na.rm = TRUE), complete_features = sum(complete_rows))
    }
  }
  result
}

checks <- data.frame(check = character(), omic = character(), assay = character(), pass = logical(), detail = character(),
  stringsAsFactors = FALSE)
addCheck <- function(check, omic, assay, pass, detail = "") {
  checks[nrow(checks) + 1L, ] <<- list(check, omic, assay, isTRUE(pass), as.character(detail))
}
traces <- list()
rules <- list()
winners <- list()
observation <- character()
verdicts <- list()

for (omic in c("lipid", "metab")) {
  case <- ruvCases()[[omic]]
  object <- case$make(case$assays(0.05))
  # 1. the unhooked sweep: the verdict is read from it
  unhooked <- ruvSweep(case, object)
  # 2. the hooked sweep (E1')
  captured <- k262WithCurveCapture(function() ruvSweep(case, object))
  hooked <- captured$value
  addCheck("E1prime", omic, "both", identical(ruvComparable(unhooked), ruvComparable(hooked)),
    "optimization_results and per-assay results, all fields, unhooked against hooked")
  addCheck("E3_recorder_clean", omic, "both", identical(captured$record_failures, 0L), paste("failures", captured$record_failures))
  aligned <- k262AlignCurves(hooked, captured$curves)
  addCheck("E3_count", omic, "both", aligned$matched, paste("captures", length(captured$curves)))
  for (assay in names(unhooked)) {
    trace <- unhooked[[assay]]$optimization_results
    cc <- m1CompleteControls(object, omic, trace$percentage_requested)
    complete <- vapply(trace$percentage_requested, function(p) as.integer(cc[[paste(assay, p)]]$c), integer(1L))
    cap <- MultiScholaR:::calculateAdaptiveMaxK(trace$sample_size[[1L]])
    ru <- k262Rules(trace, aligned$curves[[assay]], complete, cap)
    cause <- k262Cause(ru$rows, cap)
    parsed <- lapply(trace$error_reason, k262ParseReason)
    is_r5 <- trace$status == K262_R5
    # E2: every R5 reason parses, and its c equals the recomputed one (and its k the trace's)
    addCheck("E2", omic, assay, all(vapply(parsed[is_r5], function(p) p$ok, logical(1L))) &&
      all(vapply(which(is_r5), function(i) identical(parsed[[i]]$c, complete[[i]]) && identical(parsed[[i]]$k, trace$best_k[[i]]), logical(1L))),
      paste(sum(is_r5), "R5 rows"))
    # E1: P's usable set, composites and winner against the shipped trace and result
    shipped_usable <- trace$status == "ok"
    same_usable <- identical(ru$rows$p_usable, shipped_usable)
    same_comp <- isTRUE(all.equal(ru$rows$composite_recomputed[shipped_usable], trace$composite_score[shipped_usable], tolerance = 0))
    r <- unhooked[[assay]]
    same_winner <- if (isTRUE(r$success)) identical(ru$winner$P$percentage, r$best_percentage) && identical(ru$winner$P$k_used, min(as.integer(r$best_k), cap)) &&
      identical(trace$best_k[[ru$winner$P$row]], as.integer(r$best_k)) else is.na(ru$winner$P$row)
    # P is unusable on an R5 row exactly when k exceeds the recomputed c, and a row with status ok has k <= c
    same_r5 <- all(!is_r5 | (!is.na(complete) & trace$best_k > complete)) &&
      all(trace$status != "ok" | (!is.na(complete) & trace$best_k <= complete))
    addCheck("E1", omic, assay, same_usable && same_comp && same_winner && same_r5,
      sprintf("usable %s composite %s winner %s r5 %s", same_usable, same_comp, same_winner, same_r5))
    # E3: the recorded curves reproduce the trace's k
    k_rec <- vapply(seq_len(nrow(trace)), function(i) {
      cu <- aligned$curves[[assay]][[i]]
      if (trace$status[[i]] %in% K262_NOT_REACHED) return(NA_integer_)
      as.integer(findBestKElbow(list(data = cu)))
    }, integer(1L))
    reached <- !trace$status %in% K262_NOT_REACHED
    addCheck("E3", omic, assay, identical(k_rec[reached], trace$best_k[reached]), paste(sum(reached), "rows reached the elbow"))
    # E7: P within C1 within C2, k_b <= min(k, cap), where the truncated curve has finite rows
    usable_nested <- all(!ru$rows$p_usable | ru$rows$c1_usable) && all(!ru$rows$c1_usable | ru$rows$c2_usable)
    kb_ok <- all(is.na(ru$rows$k_b) | is.na(ru$rows$k) | ru$rows$k_b <= pmin(ru$rows$k, cap))
    addCheck("E7", omic, assay, usable_nested && kb_ok, sprintf("nested %s k_b %s", usable_nested, kb_ok))
    traces[[paste(omic, assay)]] <- data.frame(omic = omic, assay = assay, percentage = trace$percentage_requested,
      status = trace$status, error_reason = trace$error_reason, realized_controls = trace$realized_num_controls,
      complete_controls = complete, k = trace$best_k, separation = trace$separation_score,
      composite_shipped = trace$composite_score, sample_size = trace$sample_size, cap = cap,
      complete_features = vapply(trace$percentage_requested, function(p) cc[[paste(assay, p)]]$complete_features, integer(1L)),
      cause = cause, stringsAsFactors = FALSE)
    rules[[paste(omic, assay)]] <- cbind(omic = omic, assay = assay, ru$rows[c("percentage", "p_usable", "c1_usable", "k_b",
      "separation_b", "composite_b", "c2_usable", "finite_at_cap")], stringsAsFactors = FALSE)
    winners[[paste(omic, assay)]] <- data.frame(omic = omic, assay = assay, rule = c("P", "C1", "C2"),
      percentage = unname(vapply(ru$winner, function(w) w$percentage, numeric(1L))),
      k_used = unname(vapply(ru$winner, function(w) w$k_used, integer(1L))), cap = cap, stringsAsFactors = FALSE)
    observation <- c(observation, sprintf("KIT-190 observation %s %s: %s%% %s k=%s reason=%s", omic, assay,
      trace$percentage_requested, trace$status, trace$best_k, trace$error_reason),
      if (!isTRUE(unhooked[[assay]]$success)) sprintf("KIT-190 observation %s %s: %s", omic, assay, unhooked[[assay]]$error))
  }
  # the verdict of section 4, from the unhooked trace of assay two
  two <- unhooked$two$optimization_results
  one <- unhooked$one$optimization_results
  no_usable <- !isTRUE(unhooked$two$success)
  has_r5 <- any(two$status == K262_R5)
  all_r5 <- all(two$status == K262_R5)
  verdicts[[omic]] <- if (!has_r5 || !no_usable) "refuted" else if (all_r5) "confirmed" else "confirmed_in_part"
  addCheck("assay_one_has_ok_row", omic, "one", any(one$status == "ok"), "prediction, not a gate")
}

overall <- if (any(unlist(verdicts) == "refuted")) "refuted" else if (all(unlist(verdicts) == "confirmed")) "confirmed" else "confirmed_in_part"
write.table(do.call(rbind, traces), file.path(out, "trace.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, rules), file.path(out, "rules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, winners), file.path(out, "winners.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(checks, file.path(out, "checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(observation, file.path(out, "observation.txt"))
writeLines(c(sprintf("%s\t%s", names(verdicts), unlist(verdicts)), sprintf("overall\t%s", overall)), file.path(out, "verdict.tsv"))
cat(observation, sep = "\n")
cat(sprintf("\nverdict: %s\n", overall))
print(checks)
cat("\nwinners\n")
print(do.call(rbind, winners))
