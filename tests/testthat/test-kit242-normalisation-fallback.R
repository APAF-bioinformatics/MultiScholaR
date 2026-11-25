# KIT-242 (ticket r1.4, section 7): lipidomics and metabolomics corrections that do not happen are not recorded or
# reported as done. This file holds the option-free rows, written and run failing-first before the fix: E2 (a method
# failure stops), E3 (a malformed assay stops), E10 (the record names the FDR method per comparison and assay), E16 (an
# invalid k stops), and the negative record assertions NA1 to NA3 (the record never says a step was applied to an assay
# the step left as it was, and never names an assay the RUV state lacks). The rows that follow the owner's decided
# options are added with the fix. None of these rows runs the correlation filter (KIT-255's kernels): they assume
# nothing about which sample columns exist after it.
#
# Synthetic objects (ticket section 4): two assays A1, A2 of 40 features, 2 groups x 6 samples (KIT-172's design),
# levels about 20 log2, one loading per sample N(0, 0.5^2), noise N(0, 0.25^2), 5 standards `IS_<assay>_<i>` whose
# annotation cell holds their own name. Row k draws with set.seed(242000 + k) before any package function is called.

kit242Runs <- function() kit172Design()$Run

kit242Draw <- function(k, n = c(A1 = 40L, A2 = 40L), n_is = c(A1 = 5L, A2 = 5L)) {
  set.seed(242000L + k, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  runs <- kit242Runs()
  loading <- stats::rnorm(length(runs), 0, 0.5)
  out <- list()
  for (a in names(n)) {
    level <- stats::rnorm(n[[a]], 20, 1.5)
    noise <- matrix(stats::rnorm(n[[a]] * length(runs), 0, 0.25), n[[a]])
    y <- outer(level, rep(1, length(runs))) + outer(rep(1, n[[a]]), loading) + noise
    ids <- sprintf("%s_F%03d", a, seq_len(n[[a]]))
    if (n_is[[a]] > 0L) {
      s <- matrix(16 + stats::rnorm(n_is[[a]] * length(runs), 0, 0.1), n_is[[a]]) + outer(rep(1, n_is[[a]]), loading)
      y <- rbind(y, s)
      ids <- c(ids, sprintf("IS_%s_%d", a, seq_len(n_is[[a]])))
    }
    m <- 2^y
    dimnames(m) <- list(ids, runs)
    out[[a]] <- m
  }
  out
}

kit242Object <- function(values, class, annotations = list()) {
  tables <- lapply(stats::setNames(names(values), names(values)), function(a) {
    m <- values[[a]]
    ids <- rownames(m)
    ann <- annotations[[a]] %||% ifelse(startsWith(ids, "IS_"), ids, "PC")
    table <- if (identical(class, "lipid")) data.frame(LipidName = ids, LipidClass = ann, stringsAsFactors = FALSE) else
      data.frame(Metabolite = ids, Ontology = ann, stringsAsFactors = FALSE)
    cbind(table, stats::setNames(as.data.frame(unname(m)), colnames(m)))
  })
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(class, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = tables, design_matrix = kit172Design(), lipid_id_column = "LipidName",
      annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = tables, design_matrix = kit172Design(), metabolite_id_column = "Metabolite",
      annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  }
}

#' An assay's sample-value matrix, rows by feature id (sorted), columns the design's runs the table holds.
kit242Matrix <- function(object, binding, assay) {
  table <- methods::slot(object, binding$data_slot)[[assay]]
  if (is.null(table)) return(NULL)
  table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(table[, intersect(kit242Runs(), names(table)), drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(table[[methods::slot(object, binding$id_slot)]])
  m[order(rownames(m)), , drop = FALSE]
}

#' Equal (ticket section 7: all.equal, tolerance 0) on the non-standard rows and the columns both hold.
kit242Equal <- function(before, after) {
  rows <- setdiff(intersect(rownames(before), rownames(after)), grep("^IS_", rownames(before), value = TRUE))
  cols <- intersect(colnames(before), colnames(after))
  isTRUE(all.equal(unname(before[rows, cols, drop = FALSE]), unname(after[rows, cols, drop = FALSE]), tolerance = 0))
}

kit242Owner <- function(object, class, env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  fixture <- kit172Quiet(kit172Owner(root, object, class))
  withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
  fixture
}

kit242Outcomes <- function(owner, state) {
  record <- kit172Quiet(buildWorkflowAnalysisRecord(owner))
  for (step in record$steps) if (identical(step$state, state)) return(step$outcomes)
  stop("no record step for state ", state)
}

kit242Prepare <- function(binding, owner, config) {
  kit172Quiet({
    commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), config, "QC complete")
    runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
  })
}

for (class in c("lipid", "metab")) {
  binding <- kit172Binding(class)
  prefix <- binding$prefix

  # ---- E2: a method failure stops, under every option (G0 C2: today a warning and the assay unnormalised) -----------
  test_that(paste("KIT-242 E2:", class, "quantile on a sample holding one value stops with omics_normalisation_failed"), {
    v <- kit242Draw(2L)
    keep <- v$A2[1L, "B_6"]
    v$A2[, "B_6"] <- NA_real_
    v$A2[1L, "B_6"] <- keep
    object <- kit242Object(v, class)
    err <- expect_error(suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "quantile"))),
      class = "omics_normalisation_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$method, "quantile")
    expect_identical(err$reason, "normaliser_error")
    expect_s3_class(err$parent, "error")
    expect_match(conditionMessage(err$parent), "two non-NA values")
    # Through the declaration, in process: the same refusal, and no normalised state committed.
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(object, class)
    kit242Prepare(binding, fixture$owner, fixture$workflow$config_list)
    expect_error(suppressWarnings(suppressMessages(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, fixture$owner,
      "quantile", fixture$workflow$config_list))), class = "omics_normalisation_failed")
    expect_false(paste0(prefix, "_normalized") %in% unlist(fixture$owner$getHistory()))
  })

  # ---- E3: a malformed assay stops (UVR-082: a wrong argument stops); today a warning and the assay unchanged ---------
  test_that(paste("KIT-242 E3:", class, "an assay with no sample columns stops with omics_normalisation_failed (malformed_assay)"), {
    object <- kit242Object(kit242Draw(3L), class)
    tables <- methods::slot(object, binding$data_slot)
    tables$A2 <- tables$A2[, c(methods::slot(object, binding$id_slot), object@annotation_id_column)]
    methods::slot(object, binding$data_slot, check = FALSE) <- tables
    err <- expect_error(suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"))),
      class = "omics_normalisation_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$reason, "malformed_assay")
  })

  # ---- E16: an invalid k stops (both Q6.2 options); today the assay is dropped with a logger line ---------------------
  test_that(paste("KIT-242 E16:", class, "legacy RUV-III-C with k = -1 for A2 stops with omics_ruv_failed (invalid_k)"), {
    values <- lapply(kit242Draw(16L, n_is = c(A1 = 0L, A2 = 0L)), log2)
    object <- kit242Object(values, class)
    controls <- list(A1 = sprintf("A1_F%03d", 1:8), A2 = sprintf("A2_F%03d", 1:8))
    err <- expect_error(suppressWarnings(suppressMessages(ruvIII_C_Varying(object, ruv_grouping_variable = "group",
      ruv_number_k = list(A1 = 1L, A2 = -1L), ctrl = controls))), class = "omics_ruv_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$reason, "invalid_k")
  })

  # ---- E10: the record names the FDR method per comparison and assay (outcome 7; proteomics' rule) -------------------
  test_that(paste("KIT-242 E10:", class, "the record's DA entries carry fdr_method, BH with its reason after the q-value fallback"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(7L, n = c(A1 = 40L, A2 = 10L), n_is = c(A1 = 0L, A2 = 0L))
    v$A2[1:8, startsWith(kit242Runs(), "B_")] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    kit242Prepare(binding, owner, fixture$workflow$config_list)
    kit172Quiet(runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner),
      list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts()))
    manifest <- resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results)
    entries <- kit172Quiet(buildWorkflowAnalysisRecord(owner))$analyses$differential_abundance$results
    expect_length(entries, 2L)
    for (entry in entries) {
      fell_back <- length(manifest$results[[entry$contrast]]$qvalue_warnings[[entry$assay]]) > 0L
      expect_identical(entry$fdr_method, if (fell_back) "BH" else "qvalue", info = paste(entry$contrast, entry$assay))
      if (fell_back) expect_true(is.character(entry$fdr_fallback_reason) && nzchar(entry$fdr_fallback_reason), info = entry$assay)
    }
    a2 <- Filter(function(entry) identical(entry$assay, "A2"), entries)
    expect_identical(a2[[1L]]$fdr_method, "BH") # A2 has 2 valid p-values: the fallback is certain
  })

  # ---- NA1: never "normalised" for an assay the between-sample step left as it was ---------------------------------
  test_that(paste("KIT-242 NA1:", class, "the record states between-sample normalisation per assay, FALSE where the assay is unchanged"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(1L)
    v$A2[, "B_6"] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    kit242Prepare(binding, owner, fixture$workflow$config_list)
    try(suppressWarnings(suppressMessages(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess",
      fixture$workflow$config_list))), silent = TRUE)
    before <- owner$getState(paste0(prefix, "_log2"))
    after <- owner$getState(owner$current_state)
    skip_if(identical(owner$current_state, paste0(prefix, "_log2")), "the step refused: nothing was recorded as applied")
    outcomes <- kit242Outcomes(owner, owner$current_state)
    for (assay in c("A1", "A2")) {
      applied <- outcomes$normalisation_applied_per_assay[[assay]]
      expect_true(is.logical(applied) && length(applied) == 1L, info = paste(assay, "per-assay outcome present"))
      if (kit242Equal(kit242Matrix(before, binding, assay), kit242Matrix(after, binding, assay))) {
        expect_false(isTRUE(applied), info = paste(assay, "unchanged, so not recorded as normalised"))
      }
    }
  })

  # ---- NA2: never "ITSD applied" for an assay the ITSD step left as it was ------------------------------------------
  test_that(paste("KIT-242 NA2:", class, "the record states ITSD per assay, FALSE where no standard was matched"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(4L)
    unmatched <- list(A2 = ifelse(startsWith(rownames(v$A2), "IS_"), "ISTD", "PC"))
    fixture <- kit242Owner(kit242Object(v, class, unmatched), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet(commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"))
    try(suppressWarnings(suppressMessages(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config))),
      silent = TRUE)
    skip_if(identical(owner$current_state, paste0(prefix, "_qc_complete")), "the step refused: nothing was recorded as applied")
    before <- owner$getState(paste0(prefix, "_qc_complete"))
    after <- owner$getState(paste0(prefix, "_itsd_norm"))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_itsd_norm"))
    for (assay in c("A1", "A2")) {
      applied <- outcomes$itsd_applied_per_assay[[assay]]
      expect_true(is.logical(applied) && length(applied) == 1L, info = paste(assay, "per-assay outcome present"))
      if (kit242Equal(kit242Matrix(before, binding, assay), kit242Matrix(after, binding, assay))) {
        expect_false(isTRUE(applied), info = paste(assay, "unchanged, so not recorded as ITSD-normalised"))
      }
    }
  })

  # ---- NA3: the RUV record never names an assay the RUV state lacks, nor "applied" for one left as it was ----------
  test_that(paste("KIT-242 NA3:", class, "every assay the legacy RUV record names is in the RUV state, FALSE where uncorrected"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(12L)
    runs <- kit242Runs()
    for (j in seq_len(nrow(v$A2))) v$A2[j, runs[((j - 1L) %% length(runs)) + 1L]] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242Prepare(binding, owner, config)
    kit172Quiet(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config))
    try(suppressWarnings(suppressMessages(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config))),
      silent = TRUE)
    skip_if(identical(owner$current_state, paste0(prefix, "_normalized")), "the step refused: nothing was recorded as applied")
    before <- owner$getState(paste0(prefix, "_normalized"))
    after <- owner$getState(paste0(prefix, "_ruv_corrected"))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_ruv_corrected"))
    named <- names(outcomes$k_per_assay)
    expect_setequal(named, c("A1", "A2"))
    for (assay in named) {
      expect_false(is.null(kit242Matrix(after, binding, assay)), info = paste(assay, "named by the record, present in the RUV state"))
      applied <- outcomes$ruv_applied_per_assay[[assay]]
      expect_true(is.logical(applied) && length(applied) == 1L, info = paste(assay, "per-assay outcome present"))
      corrected <- kit242Matrix(after, binding, assay)
      if (!is.null(corrected) && kit242Equal(kit242Matrix(before, binding, assay), corrected)) {
        expect_false(isTRUE(applied), info = paste(assay, "unchanged, so not recorded as RUV-corrected"))
      }
    }
  })
}
