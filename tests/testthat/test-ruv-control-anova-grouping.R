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

test_that("KIT-190 T7 (vi): a managed legacy RUV step records the marker; a managed new-design step does not", {
  skip_if_not(workflowOwnerCapable())
  skip_if_not_installed("sva")
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
    assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
  managed <- function(parameters, env = parent.frame()) {
    # F172 on the log2 scale without missing cells: unique feature ids (the declarations fixture's default object holds
    # duplicates, which the declared designs refuse before the duplicate step) and values the declared correction needs.
    fixture <- lipidDeclarationsOwner(withr::local_tempdir(.local_envir = env),
      object = kit172Object(log2(kit172Values(172001L, missing = FALSE)$values), "lipid"))
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
