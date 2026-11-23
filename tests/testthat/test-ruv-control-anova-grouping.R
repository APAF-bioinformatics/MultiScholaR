# KIT-190: lipidomics and metabolomics legacy RUV choose their control features by the one-way ANOVA on the RUV
# grouping column (ticket section 2.6), not by the row order of the assay table. Tests T1 to T6 and T7 (ii) of the
# ticket's section 9; T7 (i) and (iii) to (vi) need KIT-187's replay-rule table and are added after it merges.
#
# Objects use the package's own classes, built as the characterisation files build them:
# `newLipidRuvObject` (test-lipid-12c-ruv-s4-characterization.R:247-293) and `newNegCtrlMetabObject`
# (test-metab-02y-get-neg-ctrl-metab-anova-characterization.R:190-215); T6's protein object as
# test-kit128-no-biological-replication-prot.R:39-44 builds one, its peptide object as
# test-dia-limpa-ruv-rules.R:53-63 does. `k190LocalBinding` and `k190MethodEnv` are copies of `localBinding` and
# `getTargetMethodEnv` of test-lipid-12c-ruv-s4-characterization.R:202-245.

k190Samples <- paste0("S", 1:8)
k190Groups <- rep(c("A", "B"), each = 4L)
k190Within <- c(-0.15, -0.05, 0.05, 0.15)
k190GroupingClass <- "multischolar_ruv_grouping_unusable"

k190MethodEnv <- function(generic_name, signature_name, fallback_env = environment()) {
  if (methods::hasMethod(generic_name, signature_name)) {
    return(environment(methods::getMethod(generic_name, signature_name)))
  }
  fallback_env
}

k190LocalBinding <- function(env, name, value, .local_envir = parent.frame()) {
  target_env <- env
  while (!exists(name, envir = target_env, inherits = FALSE) && environmentIsLocked(target_env)) {
    parent_env <- parent.env(target_env)
    if (identical(parent_env, emptyenv())) {
      break
    }
    target_env <- parent_env
  }
  had_binding <- exists(name, envir = target_env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = target_env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, target_env)
  if (was_locked) {
    unlockBinding(name, target_env)
  }
  assign(name, value, envir = target_env)
  if (was_locked) {
    lockBinding(name, target_env)
  }
  withr::defer({
    if (exists(name, envir = target_env, inherits = FALSE) && bindingIsLocked(name, target_env)) {
      unlockBinding(name, target_env)
    }
    if (had_binding) {
      assign(name, old_value, envir = target_env)
    } else if (exists(name, envir = target_env, inherits = FALSE)) {
      rm(list = name, envir = target_env)
    }
    if (was_locked && exists(name, envir = target_env, inherits = FALSE)) {
      lockBinding(name, target_env)
    }
  }, envir = .local_envir)
}

# M40 (section 9): F01-F10 group B = A + 4 log2; F11-F40 group B = A + 0.01 * (j - 10), so the oracle's p falls with j.
k190M40 <- function() {
  shift <- c(rep(4, 10L), 0.01 * (11:40 - 10))
  m <- t(vapply(shift, function(d) c(20 + k190Within, 20 + d + k190Within), numeric(8L)))
  dimnames(m) <- list(sprintf("F%02d", 1:40), k190Samples)
  m
}

# M45 (T5): M40 followed by F41-F45, group A as M40's F11, group B missing; `reversed` puts them as F45 ... F41.
k190M45 <- function(reversed = FALSE) {
  one_group <- matrix(rep(c(20 + k190Within, rep(NA_real_, 4L)), 5L), nrow = 5L, byrow = TRUE,
    dimnames = list(sprintf("F%02d", 41:45), k190Samples))
  if (reversed) {
    one_group <- one_group[rev(rownames(one_group)), , drop = FALSE]
  }
  rbind(k190M40(), one_group)
}

k190Design <- function(groups = k190Groups, samples = k190Samples) {
  data.frame(Run = samples, group = groups, factor1 = groups, stringsAsFactors = FALSE)
}

k190Object <- function(omic, m, design = k190Design(samples = colnames(m))) {
  assay <- tibble::as_tibble(data.frame(database_identifier = rownames(m), m, check.names = FALSE,
    stringsAsFactors = FALSE))
  if (identical(omic, "lipid")) {
    return(methods::new("LipidomicsAssayData", lipid_data = list(LCMS_Pos = assay),
      lipid_id_column = "database_identifier", design_matrix = design, sample_id = "Run", group_id = "group",
      args = list()))
  }
  methods::new("MetaboliteAssayData", metabolite_data = list(LCMS_Pos = assay),
    metabolite_id_column = "database_identifier", design_matrix = design, sample_id = "Run", group_id = "group",
    args = list())
}

k190Quiet <- function(expr) {
  value <- NULL
  utils::capture.output(value <- suppressWarnings(suppressMessages(expr)), type = "output")
  value
}

# The rule of section 2.6 with base R only: valid rows have at least two values (N_v, the count's base); one-way
# ANOVA per feature on the columns `anova_cols`; p = 1 when fewer than two groups are observed or F is NaN; a feature
# with zero residual degrees of freedom is excluded before BH; BH as qvalue(pi0 = 1) computes it; q > cutoff; ranked
# by q descending, ties in table order; round(N_v * pct / 100) features.
k190Oracle <- function(m, groups, pct, anova_cols = colnames(m), cutoff = 0.05) {
  m <- m[rowSums(!is.na(m)) > 1, , drop = FALSE]
  count <- round(nrow(m) * pct / 100)
  keep <- colnames(m) %in% anova_cols
  g_all <- groups[keep]
  p <- vapply(rownames(m), function(id) {
    x <- m[id, keep]
    observed <- !is.na(x)
    g <- factor(g_all[observed])
    if (nlevels(g) < 2L) {
      return(1)
    }
    fit <- stats::aov(x[observed] ~ g)
    if (fit$df.residual == 0L) {
      return(NA_real_)
    }
    value <- summary(fit)[[1]][["Pr(>F)"]][1]
    if (is.na(value)) 1 else value
  }, numeric(1L))
  p <- p[!is.na(p)]
  n_tested <- length(p)
  o <- order(p)
  q <- numeric(n_tested)
  q[o] <- pmin(1, rev(cummin(rev(p[o] * n_tested / seq_len(n_tested)))))
  names(q) <- names(p)
  candidates <- names(q)[q > cutoff]
  ranked <- candidates[order(-q[candidates])]
  chosen <- utils::head(ranked, max(min(count, length(ranked)), 1L))
  list(selected = rownames(m)[rownames(m) %in% chosen], p = p, q = q, count = count)
}

# Trace of the helper's p-values (section 9, "p trace"): applied first; `untrace` deferred before tracing.
k190Trace <- function(frame = parent.frame()) {
  k190 <- new.env(parent = emptyenv())
  ns <- asNamespace("MultiScholaR")
  withr::defer(suppressMessages(untrace("getNegCtrlProtAnovaHelper", where = ns)), envir = frame)
  suppressMessages(trace("getNegCtrlProtAnovaHelper", exit = bquote(assign("ps", ps, envir = .(k190))),
    where = ns, print = FALSE))
  k190
}

# Pass-through double (T1): records the design and the grouping the helper will read, asserts nothing, returns the
# traced helper's value with every argument passed on.
k190PassThrough <- function(signature_name, frame = parent.frame()) {
  rec <- new.env(parent = emptyenv())
  rec$calls <- 0L
  traced <- get("getNegCtrlProtAnovaHelper", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  double <- function(data_matrix, design_matrix, grouping_variable = "group", ...) {
    rec$calls <- rec$calls + 1L
    rec$design <- design_matrix
    row <- design_matrix[colnames(data_matrix), grouping_variable]
    rec$row_is_null <- is.null(row)
    traced(data_matrix, design_matrix = design_matrix, grouping_variable = grouping_variable, ...)
  }
  k190LocalBinding(k190MethodEnv("getNegCtrlMetabAnova", signature_name), "getNegCtrlProtAnovaHelper", double,
    .local_envir = frame)
  rec
}

k190Signature <- c(lipid = "LipidomicsAssayData", metab = "MetaboliteAssayData")

k190Selected <- function(result) names(result)[result]

test_that("KIT-190 self-check: the trace records the helper's p-values", {
  k190 <- k190Trace()
  m <- k190M40()
  design <- data.frame(group = k190Groups, row.names = k190Samples)
  # Called through the namespace binding the trace modified (the S4 methods' enclosure), not the attached copy that a
  # name looked up from the test environment finds.
  traced <- get("getNegCtrlProtAnovaHelper", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  k190Quiet(traced(m, design, grouping_variable = "group", percentage_as_neg_ctrl = 10, ruv_fdr_method = "BH"))
  expect_length(unlist(k190$ps), 40L)
})

for (omic in names(k190Signature)) {
  test_that(sprintf("KIT-190 T1 (%s): controls are chosen by the ANOVA on the grouping, not by row order", omic), {
    k190 <- k190Trace()
    rec <- k190PassThrough(k190Signature[[omic]])
    m <- k190M40()
    oracle <- k190Oracle(m, k190Groups, 10)
    result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "group",
      percentage_as_neg_ctrl = 10))

    # Instrument checks: one helper call, the assay returned with 40 named values, 40 traced p-values.
    expect_identical(rec$calls, 1L)
    expect_length(result, 1L)
    expect_true(is.logical(result[[1]]) && length(result[[1]]) == 40L && !is.null(names(result[[1]])))
    expect_length(unlist(k190$ps), 40L)

    # Values (base: F01-F04, every p = 1, the design row NULL).
    expect_identical(k190Selected(result[[1]]), oracle$selected)
    expect_false(any(k190Selected(result[[1]]) %in% sprintf("F%02d", 1:10)))
    expect_equal(unlist(k190$ps)[names(oracle$p)], oracle$p, tolerance = 1e-10)
    expect_false(isTRUE(rec$row_is_null))
  })

  test_that(sprintf("KIT-190 T2 (%s): another grouping column is unchanged", omic), {
    k190 <- k190Trace()
    m <- k190M40()
    oracle <- k190Oracle(m, k190Groups, 10)
    result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "factor1",
      percentage_as_neg_ctrl = 10))
    expect_identical(k190Selected(result[[1]]), oracle$selected)
    expect_equal(unlist(k190$ps)[names(oracle$p)], oracle$p, tolerance = 1e-10)
  })

  test_that(sprintf("KIT-190 T3 (%s): an unusable grouping is refused with a classed error", omic), {
    m <- k190M40()
    # (ii) the method with a grouping column absent from the design.
    expect_error(k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "absent",
      percentage_as_neg_ctrl = 10)), class = k190GroupingClass)
    # (iii) a design whose only other group is a pool.
    pool_design <- k190Design(groups = rep(c("A", "Pool"), each = 4L))
    err <- tryCatch(k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m, pool_design), ruv_grouping_variable = "group",
      percentage_as_neg_ctrl = 10)), error = function(e) e)
    expect_s3_class(err, k190GroupingClass)
    if (inherits(err, "condition")) {
      # cli wraps the message to the console width: words are matched across any white space.
      expect_match(conditionMessage(err), "after\\s+leaving\\s+out\\s+the\\s+pooled\\s+samples,\\s+fewer\\s+than\\s+two\\s+groups")
      expect_match(conditionMessage(err$parent), "After excluding Pool samples, only 1 group(s) remain", fixed = TRUE)
    }
  })

  test_that(sprintf("KIT-190 T3 (iv) (%s): the sweeps pass the refusal to the caller", omic), {
    local_mocked_bindings(getNegCtrlMetabAnova = function(...) {
      stop(errorCondition("grouping unusable (test double)", class = k190GroupingClass))
    }, .package = "MultiScholaR")
    object <- k190Object(omic, k190M40())
    sweep <- if (identical(omic, "lipid")) runLipidPerAssayRuvOptimization else runPerAssayRuvOptimization
    expect_error(k190Quiet(sweep(object, "automatic", params = list(ruv_grouping_variable = "group",
      percentage_min = 1, percentage_max = 2))), class = k190GroupingClass)
    expect_error(k190Quiet(sweep(object, "manual", params = list(ruv_grouping_variable = "group", manual_k = 1L,
      manual_percentage = 10))), class = k190GroupingClass)
  })

  test_that(sprintf("KIT-190 T4 (%s): pool samples leave the ANOVA", omic), {
    samples <- paste0("S", 1:12)
    groups <- rep(c("A", "B", "QC_pool"), each = 4L)
    pool_shift <- c(rep(0, 10L), rep(3, 10L), rep(0, 20L))
    m <- cbind(k190M40(), t(vapply(pool_shift, function(d) 20 + d + k190Within, numeric(4L))))
    colnames(m) <- samples
    without_pool <- k190Oracle(m, groups, 10, anova_cols = samples[1:8])
    with_pool <- k190Oracle(m, groups, 10)
    if (identical(without_pool$selected, with_pool$selected)) {
      stop("KIT-190 T4 is vacuous: the pool samples do not move the oracle's selection")
    }
    result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m, k190Design(groups, samples)),
      ruv_grouping_variable = "group", percentage_as_neg_ctrl = 10))
    expect_identical(k190Selected(result[[1]]), without_pool$selected)
  })

  test_that(sprintf("KIT-190 T5 (%s): one-group features score p = 1 and ties keep the table order", omic), {
    k190 <- k190Trace()
    selections <- lapply(c(forward = FALSE, reversed = TRUE), function(reversed) {
      m <- k190M45(reversed)
      oracle <- k190Oracle(m, k190Groups, 8)
      result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "group",
        percentage_as_neg_ctrl = 8))
      ps <- unlist(k190$ps)
      expect_identical(oracle$count, 4)
      expect_equal(unname(ps[sprintf("F%02d", 41:45)]), rep(1, 5L))
      expect_equal(unname(oracle$q[c("F11", sprintf("F%02d", 41:45))]), rep(1, 6L))
      expect_identical(k190Selected(result[[1]]), oracle$selected)
      oracle$selected
    })
    if (identical(sort(selections$forward), sort(selections$reversed))) {
      stop("KIT-190 T5 is vacuous: the table-order tie-break is not exercised")
    }
    expect_setequal(selections$forward, c("F11", "F41", "F42", "F43"))
    expect_setequal(selections$reversed, c("F11", "F45", "F44", "F43"))
  })
}

test_that("KIT-190 T3 (i): the helper refuses a design without its grouping column", {
  design <- data.frame(factor2 = k190Groups, row.names = k190Samples)
  expect_error(k190Quiet(getNegCtrlProtAnovaHelper(k190M40(), design, grouping_variable = "group",
    percentage_as_neg_ctrl = 10, ruv_fdr_method = "BH")), class = k190GroupingClass)
})

test_that("KIT-190 T6: proteomics and peptide selections follow the ANOVA (unaffected)", {
  k190 <- k190Trace()
  m <- k190M40()
  oracle <- k190Oracle(m, k190Groups, 10)
  design <- data.frame(Run = k190Samples, group = k190Groups, stringsAsFactors = FALSE)

  protein <- ProteinQuantitativeData(
    protein_quant_table = data.frame(Protein.Ids = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE,
      row.names = NULL),
    protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(Protein.Ids = rownames(m), stringsAsFactors = FALSE),
    design_matrix = design, sample_id = "Run", group_id = "group", technical_replicate_id = "replicates",
    args = list())
  protein_result <- k190Quiet(getNegCtrlProtAnova(protein, ruv_grouping_variable = "group",
    percentage_as_neg_ctrl = 10, ruv_fdr_method = "BH"))
  expect_identical(k190Selected(protein_result), oracle$selected)
  expect_false(any(k190Selected(protein_result) %in% sprintf("F%02d", 1:10)))
  expect_equal(unlist(k190$ps)[names(oracle$p)], oracle$p, tolerance = 1e-10)

  rm(list = ls(k190), envir = k190)
  peptide <- methods::new("PeptideQuantitativeData")
  peptide@peptide_data <- data.frame(Protein.Ids = rownames(m), Stripped.Sequence = rownames(m),
    stringsAsFactors = FALSE)
  peptide@peptide_matrix <- m
  peptide@protein_id_column <- "Protein.Ids"
  peptide@peptide_sequence_column <- "Stripped.Sequence"
  peptide@is_logged_data <- TRUE
  peptide@design_matrix <- design
  peptide@sample_id <- "Run"
  peptide@group_id <- "group"
  peptide@args <- list()
  peptide_result <- k190Quiet(getNegCtrlProtAnovaPeptides(peptide, ruv_grouping_variable = "group",
    percentage_as_neg_ctrl = 10, ruv_fdr_method = "BH"))
  expect_identical(k190Selected(peptide_result), oracle$selected)
  expect_false(any(k190Selected(peptide_result) %in% sprintf("F%02d", 1:10)))
  expect_equal(unlist(k190$ps)[names(oracle$p)], oracle$p, tolerance = 1e-10)
})

test_that("KIT-190 T7 (ii): under the earlier rule the selection is the row-order one", {
  withr::local_options(MultiScholaR.ruv_control_grouping = "dropped")
  k190 <- k190Trace()
  for (omic in names(k190Signature)) {
    rm(list = ls(k190), envir = k190)
    m <- k190M40()
    result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "group",
      percentage_as_neg_ctrl = 10))
    expect_identical(k190Selected(result[[1]]), sprintf("F%02d", 1:4))
    expect_equal(unname(unlist(k190$ps)), rep(1, 40L))
  }
})
