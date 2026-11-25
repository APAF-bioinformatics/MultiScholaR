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
  one_group <- rownames(m)[apply(!is.na(m[, keep, drop = FALSE]), 1L, function(x) length(unique(g_all[x])) < 2L)]
  list(selected = rownames(m)[rownames(m) %in% chosen], p = p, q = q, count = count, one_group = one_group,
    rows = rownames(m))
}

# Rule (d) of section 2.7a (D2, VAL-110; lipid and metab, VAL-126) on an oracle's p and q: candidates q > 0.05 less the
# features observed in one group (which stay in BH's family), ranked by p descending, ties by identifier in C-locale
# byte order, round(N_v * pct / 100) of them (all if fewer).
k190RuleD <- function(oracle, cutoff = 0.05) {
  candidates <- setdiff(names(oracle$q)[oracle$q > cutoff], oracle$one_group)
  ranked <- candidates[order(-oracle$p[candidates], candidates, method = "radix")]
  chosen <- utils::head(ranked, max(min(oracle$count, length(ranked)), 1L))
  oracle$rows[oracle$rows %in% chosen]
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

  # Under D2's rule (d) (VAL-110, VAL-126) the one-group features F41 to F45 still score p = 1 but are not candidates,
  # so both orders select the four highest p among the two-group features (section 9: "under another D2 option the
  # expected values follow that option").
  test_that(sprintf("KIT-190 T5 (%s): one-group features score p = 1 and are never chosen, in either order", omic), {
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
      chosen <- k190RuleD(oracle)
      expect_identical(k190Selected(result[[1]]), chosen)
      chosen
    })
    expect_setequal(selections$forward, c("F11", "F12", "F13", "F14"))
    expect_setequal(selections$reversed, c("F11", "F12", "F13", "F14"))
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

# ---- T7: records and replay (D4 (a)); KIT-187's shared table `workflowReplayRuleMarkers` -----------------------------

k190Marker <- function(record, state) {
  step <- Filter(function(s) identical(s$state, state), record$steps)
  if (!length(step)) return(NULL)
  rows <- Filter(function(p) identical(p$name, "ruv_control_grouping"), step[[1L]]$parameters)
  if (!length(rows)) return(NULL)
  list(source = rows[[1L]]$source, value = decodeWorkflowRecordValue(rows[[1L]]$value))
}

k190FixedNames <- function(record, state) {
  step <- Filter(function(s) identical(s$state, state), record$steps)
  vapply(Filter(function(p) identical(p$source, "fixed"), step[[1L]]$parameters), function(p) p$name, "")
}

# As test-kit172-reopen.R:19-23.
k190Paths <- function(root, class) {
  qc <- if (identical(class, "lipid")) "lipid_qc_dir" else "metabolite_qc_dir"
  stats::setNames(list(file.path(root, "source"), file.path(root, "qc"), file.path(root, "da"), file.path(root, "pub")),
    c("source_dir", qc, "da_output_dir", "publication_graphs_dir"))
}

# A KIT-172 journey (helper-kit172.R:110-130) to the normalised state, then its RUV step: `legacy = TRUE` computes it
# under the earlier rule, so its provenance has no marker, as a build before this change wrote it. Rows: F172 seed
# 172101 with planted effects, the planted features first, so the row-order and the ANOVA selections differ.
k190Owner <- function(root, class, ruv, legacy = FALSE) {
  sim <- kit172Values(172101L, planted = TRUE)
  values <- sim$values[order(sim$effects == 0), , drop = FALSE]
  journey <- kit172Journey(root, kit172Object(values, class), class, ruv = NULL, da = FALSE)
  run <- function() kit172Quiet(runWorkflowOmicsReferenceRuv(journey$binding, journey$owner, ruv,
    journey$fixture$workflow$config_list))
  if (legacy) withr::with_options(list(MultiScholaR.ruv_control_grouping = "dropped"), run()) else run()
  journey
}

# Write the record, close the owner, prepare a copy cut to the design state and replay it (as test-kit172-reopen.R:
# 44-56 does); `during` wraps the replay (spies). Returns the record, the replayed project's record, the mismatches and
# the corrected matrix of the origin.
k190Replay <- function(journey, root, class, during = function(expr) expr) {
  prefix <- journey$prefix
  corrected <- kit172StateMatrix(journey, paste0(prefix, "_ruv_corrected"))
  json <- file.path(root, "analysis_record.json")
  writeWorkflowAnalysisRecord(journey$owner, json)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  closeDIAStateAdapter(journey$owner)
  source_dir <- file.path(root, "source")
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE,
    full.names = TRUE)[[1L]]))
  copy <- tempfile("k190-copy-")
  dir.create(copy)
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  harness <- new.env(parent = globalenv())
  sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
  api <- function(name) get(name, asNamespace("MultiScholaR"))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", paste0(prefix, "_raw_data_s4")), api)
  inputs <- c(stats::setNames(file.path(source_dir, "import_sources", "LCMS_Pos_export.txt"), paste0(prefix, "_assay_1")),
    design_file = file.path(source_dir, "design_matrix.tab"))
  replay <- during(kit172Quiet(harness$kitReplayJourney(record, prepared, store_root, inputs, k190Paths(root, class),
    tempfile("k190-scratch-"), api)))
  identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
  replayed <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
  replayed_record <- buildWorkflowAnalysisRecord(replayed)
  closeDIAStateAdapter(replayed)
  list(record = record, replayed_record = replayed_record, mismatches = as.character(replay$mismatches),
    corrected = corrected, state = paste0(prefix, "_ruv_corrected"))
}

test_that("KIT-190 T7 (i): the replay rule of the control grouping follows the step's marker", {
  markers <- get0("workflowReplayRuleMarkers", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(!is.null(markers[["ruv_control_grouping"]]))
  if (!is.null(markers[["ruv_control_grouping"]])) {
    entry <- function(source) list(name = "ruv_control_grouping", source = source, value = encodeWorkflowRecordValue("kept"),
      provenance = "fixed")
    for (operation in c("lipid_ruv", "metab_ruv")) {
      step <- function(parameters) list(state = "s", operation = operation, parameters = parameters)
      expect_identical(workflowReplayRule(step(list(entry("fixed"))), "ruv_control_grouping"), "kept", info = operation)
      expect_identical(workflowReplayRule(step(list()), "ruv_control_grouping"), "dropped", info = operation)
      expect_identical(workflowReplayRule(step(list(entry("registry"))), "ruv_control_grouping"), "dropped", info = operation)
    }
  }
})

test_that("KIT-190 T7 (iii): a legacy RUV step made now lists the marker; one made under the earlier rule lists none", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  journey <- k190Owner(root, "lipid", kit172LegacyRuv("automatic"))
  withr::defer(closeDIAStateAdapter(journey$owner))
  expect_identical(k190Marker(buildWorkflowAnalysisRecord(journey$owner), "lipid_ruv_corrected"),
    list(source = "fixed", value = "kept"))
  # The builder does not back-fill the marker into a step whose provenance has no fixed constants (KIT-187).
  other <- k190Owner(withr::local_tempdir(), "lipid", kit172LegacyRuv("manual"), legacy = TRUE)
  withr::defer(closeDIAStateAdapter(other$owner))
  expect_null(k190Marker(buildWorkflowAnalysisRecord(other$owner), "lipid_ruv_corrected"))
})

for (class in c("lipid", "metab")) {
  test_that(sprintf("KIT-190 T7 (iv) (%s): records made before and after the change replay with their own selection", class), {
    skip_if_not(workflowOwnerCapable())
    root_a <- withr::local_tempdir()
    root_b <- withr::local_tempdir()
    a_owner <- k190Owner(root_a, class, kit172LegacyRuv("manual"), legacy = TRUE)
    b_owner <- k190Owner(root_b, class, kit172LegacyRuv("manual"))
    a <- k190Replay(a_owner, root_a, class)
    b <- k190Replay(b_owner, root_b, class)
    if (identical(a$corrected, b$corrected)) {
      stop("KIT-190 T7 (iv) is vacuous: the two selections give the same corrected state")
    }
    expect_identical(a$mismatches, character())
    expect_identical(b$mismatches, character())
    expect_null(k190Marker(a$record, a$state))
    expect_null(k190Marker(a$replayed_record, a$state))
    expect_identical(k190Marker(b$record, b$state), list(source = "fixed", value = "kept"))
    expect_identical(k190Marker(b$replayed_record, b$state), list(source = "fixed", value = "kept"))
  })
}

test_that("KIT-190 T7 (v): a new-design RUV step replays unchanged and never reaches the control selection", {
  skip_if_not(workflowOwnerCapable())
  skip_if_not_installed("sva")
  withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
  root <- withr::local_tempdir()
  ruv <- list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = 1L, ruv_k_rule = "user",
    ruv_k_seed = NULL, ruv_control_ids = NULL)
  journey <- k190Owner(root, "lipid", ruv)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  ns <- asNamespace("MultiScholaR")
  original_method <- get("getNegCtrlMetabAnova", envir = ns)
  original_helper <- get("getNegCtrlProtAnovaHelper", envir = ns)
  spied <- function(expr) {
    testthat::local_mocked_bindings(
      getNegCtrlMetabAnova = function(...) { calls$n <- calls$n + 1L; original_method(...) },
      getNegCtrlProtAnovaHelper = function(...) { calls$n <- calls$n + 1L; original_helper(...) },
      .package = "MultiScholaR")
    expr
  }
  result <- k190Replay(journey, root, "lipid", during = spied)
  expect_identical(result$mismatches, character())
  expect_identical(calls$n, 0L)
  for (record in list(result$record, result$replayed_record)) {
    fixed <- k190FixedNames(record, result$state)
    expect_true(all(c("replicate_sets", "test") %in% fixed))
    expect_false("ruv_control_grouping" %in% fixed)
  }
})

# Copy of kit129OmicsObject (test-kit129-run-path.R:14-45), lipidomics, no technical replicates.
k190Kit129Object <- function(seed = 129L, p = 50L) {
  set.seed(seed)
  design <- data.frame(Sample_ID = sprintf("S%02d", 1:12), group = rep(c("A", "B"), each = 6L), stringsAsFactors = FALSE)
  f <- stats::residuals(stats::lm(stats::rnorm(12L) ~ factor(design$group)))
  f <- f / sqrt(mean(f^2))
  tables <- lapply(c(pos = "pos", neg = "neg"), function(assay) {
    effect <- c(rep(1, 5L), rep(-1, 5L), rep(0, p - 10L))
    y <- stats::rnorm(p, 20, 2) + outer(effect, as.numeric(design$group == "B")) + outer(stats::rnorm(p, sd = 0.6), f) +
      matrix(stats::rnorm(p * 12L, sd = 0.5), p)
    seen <- y >= stats::quantile(y, 0.08)
    kept <- rowSums(seen[, design$group == "A"]) >= 2L & rowSums(seen[, design$group == "B"]) >= 2L
    seen[!kept, ] <- TRUE
    y[!seen] <- NA_real_
    dimnames(y) <- list(sprintf("%s_f%03d", assay, seq_len(p)), design$Sample_ID)
    table <- tibble::as_tibble(as.data.frame(y, check.names = FALSE))
    table <- tibble::add_column(table, database_identifier = rownames(y), .before = 1L)
    tibble::add_column(table, name = rownames(y), .after = 1L)
  })
  methods::new("LipidomicsAssayData", lipid_data = tables, lipid_id_column = "database_identifier",
    annotation_id_column = "name", design_matrix = design, sample_id = "Sample_ID", group_id = "group", args = list())
}

test_that("KIT-190 T7 (vi): a managed legacy RUV step records the marker; a managed new-design step does not", {
  skip_if_not(workflowOwnerCapable())
  skip_if_not_installed("sva")
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
  managed <- function(parameters, env = parent.frame()) {
    # The object test-kit129-run-path.R runs through the worker (its kit129OmicsObject, lipid, copied below): the
    # declarations fixture's default object holds duplicate ids, which the declared designs refuse.
    fixture <- lipidDeclarationsOwner(withr::local_tempdir(.local_envir = env), object = k190Kit129Object())
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    worker <- omicsWorkerRun(fixture, "lipid_ruv", parameters)
    expect_identical(worker$status, "complete", info = worker$error)
    k190Quiet(completeWorkflowAction(list(manager = fixture$owner, accepting = TRUE), "lipid_ruv", worker$result,
      context = list(state_name = "lipid_ruv_corrected", config = list(), description = "RUV (managed)"),
      values = parameters))
    buildWorkflowAnalysisRecord(fixture$owner)
  }
  legacy <- managed(lipidDeclarationsRuvParameters())
  expect_identical(k190Marker(legacy, "lipid_ruv_corrected"), list(source = "fixed", value = "kept"))
  withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
  declared <- managed(list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = 1L,
    ruv_k_rule = "user", ruv_k_seed = NULL, ruv_control_ids = NULL))
  expect_null(k190Marker(declared, "lipid_ruv_corrected"))
})

# ---- T8: rule (d) for lipid and metab (VAL-110, VAL-126); proteomics and peptide unchanged ---------------------------

# M40 followed by three exact copies of F11 named so that table order (T_c, T_a, T_b) and identifier order (T_a, T_b,
# T_c) differ: a tie block of four equal p-values at the top (F11 and its copies); at 5 % of 43 rows, 2 controls.
k190TieMatrix <- function() {
  copies <- k190M40()[rep("F11", 3L), , drop = FALSE]
  rownames(copies) <- c("T_c", "T_a", "T_b")
  rbind(k190M40(), copies)
}

for (omic in names(k190Signature)) {
  test_that(sprintf("KIT-190 T8 (i) (%s): a tie block is decided by the feature identifier, not the row order", omic), {
    m <- k190TieMatrix()
    oracle <- k190Oracle(m, k190Groups, 5)
    expect_identical(k190RuleD(oracle), c("F11", "T_a"))
    # Non-vacuity: the shipped (table-order) rule would take F11 and T_c.
    if (identical(oracle$selected, k190RuleD(oracle))) stop("KIT-190 T8 (i) is vacuous: both rules select the same")
    result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "group",
      percentage_as_neg_ctrl = 5))
    expect_setequal(k190Selected(result[[1]]), c("F11", "T_a"))
  })

  test_that(sprintf("KIT-190 T8 (ii) (%s): a feature observed in one group is never chosen", omic), {
    m <- k190M45()
    result <- k190Quiet(getNegCtrlMetabAnova(k190Object(omic, m), ruv_grouping_variable = "group",
      percentage_as_neg_ctrl = 8))
    expect_false(any(k190Selected(result[[1]]) %in% sprintf("F%02d", 41:45)))
    expect_setequal(k190Selected(result[[1]]), k190RuleD(k190Oracle(m, k190Groups, 8)))
  })

  test_that(sprintf("KIT-190 T8 (iii) (%s): a permuted table gives the same controls", omic), {
    m <- rbind(k190TieMatrix(), k190M45()[sprintf("F%02d", 41:45), , drop = FALSE])
    permuted <- m[c(48L, 3L, 44L, 41L, 27L, 12L, setdiff(seq_len(nrow(m)), c(48L, 3L, 44L, 41L, 27L, 12L))), , drop = FALSE]
    select <- function(x) sort(k190Selected(k190Quiet(getNegCtrlMetabAnova(k190Object(omic, x),
      ruv_grouping_variable = "group", percentage_as_neg_ctrl = 10))[[1]]))
    expect_identical(select(permuted), select(m))
  })
}

test_that("KIT-190 T8 (iv): proteomics and peptide keep the shipped rule (q descending, table order)", {
  m <- k190TieMatrix()
  oracle <- k190Oracle(m, k190Groups, 5)
  if (identical(oracle$selected, k190RuleD(oracle))) stop("KIT-190 T8 (iv) is vacuous: both rules select the same")
  design <- data.frame(Run = k190Samples, group = k190Groups, stringsAsFactors = FALSE)
  protein <- ProteinQuantitativeData(
    protein_quant_table = data.frame(Protein.Ids = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE,
      row.names = NULL),
    protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(Protein.Ids = rownames(m), stringsAsFactors = FALSE),
    design_matrix = design, sample_id = "Run", group_id = "group", technical_replicate_id = "replicates", args = list())
  expect_identical(k190Selected(k190Quiet(getNegCtrlProtAnova(protein, ruv_grouping_variable = "group",
    percentage_as_neg_ctrl = 5, ruv_fdr_method = "BH"))), oracle$selected)
  peptide <- methods::new("PeptideQuantitativeData")
  peptide@peptide_data <- data.frame(Protein.Ids = rownames(m), Stripped.Sequence = rownames(m), stringsAsFactors = FALSE)
  peptide@peptide_matrix <- m
  peptide@protein_id_column <- "Protein.Ids"
  peptide@peptide_sequence_column <- "Stripped.Sequence"
  peptide@is_logged_data <- TRUE
  peptide@design_matrix <- design
  peptide@sample_id <- "Run"
  peptide@group_id <- "group"
  peptide@args <- list()
  expect_identical(k190Selected(k190Quiet(getNegCtrlProtAnovaPeptides(peptide, ruv_grouping_variable = "group",
    percentage_as_neg_ctrl = 5, ruv_fdr_method = "BH"))), oracle$selected)
  # The helper called as every other caller calls it (no rule argument) selects the same.
  design_helper <- data.frame(group = k190Groups, row.names = k190Samples)
  expect_identical(k190Selected(k190Quiet(getNegCtrlProtAnovaHelper(m, design_helper, grouping_variable = "group",
    percentage_as_neg_ctrl = 5, ruv_fdr_method = "BH"))), oracle$selected)
})

# ---- T9: the bound on the automatic k for lipidomics and metabolomics (VAL-182, VAL-235) ------------------------------
# The untouched sweep, then k = min(k_sweep, calculateAdaptiveMaxK(samples of the assay)) in legacy Automatic only;
# k_sweep, k_used and the cap recorded; a plain-language notice when k is reduced; a warning (k unchanged) when a Manual k
# is above the cap; per-assay RUV status recorded (KIT-242 addendum); a replay marker so earlier records replay unbounded.

k190RuvStep <- function(record) Filter(function(s) identical(s$state, "lipid_ruv_corrected"), record$steps)[[1L]]

k190BoundJourney <- function(ruv, options = list()) {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  # F172 (seed 172001), whose automatic search chose a k above 2 under rule (d) in KIT-190's measurement (`ed/`).
  journey <- withr::with_options(options, kit172Quiet(kit172Journey(root, kit172Object(kit172Values(172001L)$values, "lipid"),
    "lipid", ruv = ruv, da = FALSE)))
  withr::defer(closeDIAStateAdapter(journey$owner), envir = parent.frame())
  journey
}

test_that("KIT-190 T9 (i): Automatic k is the sweep's k bounded by the sample-size cap, recorded with a notice", {
  skip_if_not(workflowOwnerCapable())
  journey <- k190BoundJourney(kit172LegacyRuv("automatic"))
  value <- journey$ruv$value
  expect_false(is.null(value$k_sweep))
  k_sweep <- unlist(value$k_sweep)
  expect_true(all(k_sweep > 2L))
  expect_identical(unlist(value$k_cap), c(LCMS_Pos = 2L))
  expect_identical(unlist(value$best_k), c(LCMS_Pos = 2L))
  notice <- value$k_notice[["LCMS_Pos"]]
  expect_true(is.character(notice) && grepl("LCMS_Pos", notice, fixed = TRUE) &&
    grepl(as.character(k_sweep[[1L]]), notice, fixed = TRUE) && grepl("12 samples", notice, fixed = TRUE))
  outcomes <- k190RuvStep(buildWorkflowAnalysisRecord(journey$owner))$outcomes
  expect_identical(unlist(outcomes$k_per_assay), c(LCMS_Pos = 2L))
  expect_identical(unlist(outcomes$k_sweep_per_assay), k_sweep)
  expect_identical(unlist(outcomes$k_cap_per_assay), c(LCMS_Pos = 2L))
  expect_identical(unlist(outcomes$ruv_status_per_assay), c(LCMS_Pos = "corrected"))
})

test_that("KIT-190 T9 (ii): a Manual k above the cap is kept, with a classed warning naming the cap", {
  skip_if_not(workflowOwnerCapable())
  ruv <- kit172LegacyRuv("manual")
  ruv$manual_k <- 3L
  root <- withr::local_tempdir()
  object <- kit172Object(kit172Values(172001L)$values, "lipid")
  journey <- kit172Quiet(kit172Journey(root, object, "lipid", ruv = NULL, da = FALSE))
  withr::defer(closeDIAStateAdapter(journey$owner))
  warned <- NULL
  run <- withCallingHandlers(suppressMessages(runWorkflowOmicsReferenceRuv(journey$binding, journey$owner, ruv,
    journey$fixture$workflow$config_list)), warning = function(w) {
      if (inherits(w, "multischolar_ruv_k_notice")) warned <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    })
  expect_identical(unlist(run$value$best_k), c(LCMS_Pos = 3L))
  expect_true(is.character(warned) && grepl("LCMS_Pos", warned, fixed = TRUE) && grepl("3", warned, fixed = TRUE) &&
    grepl("2", warned, fixed = TRUE))
})

test_that("KIT-190 T9 (iii): the bound's replay marker; under the earlier rule the Automatic k is the sweep's", {
  skip_if_not(workflowOwnerCapable())
  markers <- get0("workflowReplayRuleMarkers", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_identical(markers$ruv_k_bound$current, "sample_size")
  expect_identical(markers$ruv_k_bound$earlier, "none")
  expect_identical(workflowOperationEntry("lipid_ruv")$declaration$fixed$ruv_k_bound, "sample_size")
  journey <- k190BoundJourney(kit172LegacyRuv("automatic"), options = list(MultiScholaR.ruv_k_bound = "none"))
  value <- journey$ruv$value
  expect_identical(unlist(value$best_k), unlist(value$k_sweep))
  step <- k190RuvStep(buildWorkflowAnalysisRecord(journey$owner))
  names_fixed <- vapply(Filter(function(p) identical(p$source, "fixed"), step$parameters), function(p) p$name, "")
  expect_false("ruv_k_bound" %in% names_fixed)
})

test_that("KIT-190 T9 (iv): the tab and the export show the k removed, and the notices reach the tab on both routes", {
  # A capped Automatic step's value (as `omicsRuvOperation` returns it) and its sweep result for one assay.
  value <- list(best_k = list(LCMS_Pos = 2L), k_sweep = list(LCMS_Pos = 8L), k_cap = list(LCMS_Pos = 2L),
    k_samples = list(LCMS_Pos = 12L), k_notice = list(LCMS_Pos = "LCMS_Pos: the automatic search found 8 unwanted factors."))
  results <- list(LCMS_Pos = list(success = TRUE, best_k = 8L, best_percentage = 11, separation_score = 0.5,
    control_genes_index = c(a = TRUE, b = TRUE, c = FALSE)))
  annotate <- get0("omicsRuvAnnotateResults", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(is.function(annotate))
  if (is.function(annotate)) {
    annotated <- annotate(results, value)
    summary <- omicNormRuvSummary(annotated$LCMS_Pos)
    expect_match(summary, "Factors removed: 2 (the search found 8; at most 2 with 12 samples)", fixed = TRUE)
    expect_false(grepl("Best k: 8", summary, fixed = TRUE))
    # Not capped: the summary is as before.
    plain <- annotate(list(LCMS_Pos = utils::modifyList(results$LCMS_Pos, list(best_k = 2L))),
      utils::modifyList(value, list(k_sweep = list(LCMS_Pos = 2L), k_notice = list(LCMS_Pos = NULL))))
    expect_match(omicNormRuvSummary(plain$LCMS_Pos), "^Best k: 2\n")
  }
  report <- get0("omicsRuvReportNotices", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(is.function(report))
  if (is.function(report)) {
    logged <- character()
    warned <- character()
    report(value, function(text) logged <<- c(logged, text), function(text) warned <<- c(warned, text))
    expect_identical(warned, value$k_notice$LCMS_Pos)
    expect_true(any(grepl(value$k_notice$LCMS_Pos, logged, fixed = TRUE)))
  }
  # Both routes deliver the annotation and the notices (in process: the tab's notify; managed: showNotification).
  for (route in c("runWorkflowOmicsInProcessNormalisation", "runWorkflowOmicsManagedNormalisation")) {
    body_text <- paste(deparse(body(get(route, envir = asNamespace("MultiScholaR")))), collapse = "\n")
    expect_true(grepl("omicsRuvAnnotateResults", body_text, fixed = TRUE), info = route)
    expect_true(grepl("omicsRuvReportNotices", body_text, fixed = TRUE), info = route)
  }
  # The session-export summary prints the k removed with the k found.
  line_fn <- get0("omicsRuvExportSummaryLine", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  expect_true(is.function(line_fn))
  if (is.function(line_fn) && is.function(annotate)) {
    expect_match(line_fn("LCMS_Pos", annotate(results, value)$LCMS_Pos), "k=2 (found 8, at most 2 with 12 samples)",
      fixed = TRUE)
  }
})
