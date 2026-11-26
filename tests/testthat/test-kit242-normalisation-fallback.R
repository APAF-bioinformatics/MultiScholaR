# KIT-242 (ticket r1.5, section 7), Part A: lipidomics and metabolomics normalisation that does not happen is not
# recorded or reported as done (VAL-242). Rows of Part A, written and run failing-first before its fix: E1 (a sample
# with no value in one assay: the assay normalised without it, the sample empty there only; OD-1 (c)), E2 (a method
# failure stops), E3 (a malformed assay stops), E4 and E5 (ITSD with no standard matched: not applied to that assay,
# recorded with the count of standards used; OD-2 (b)), E6 (a sample without standard signal: missing in that assay
# only, the assay normalised without it; OD-3 (d)), E7 (the record per assay), E8 (the tab's notice, managed and in
# process), E9 (the session summary's ITSD line), E10 (the record names the FDR method), E11 (older records: the
# normalisation markers), E13 (KIT-145's guard counts the other samples), E17 (Part C: a proteomics sample with no
# value leaves the design), E26 (a sample with no value in any assay leaves the design), NA1 and NA2 (the record never
# says applied for an assay left unchanged), and the guards G1, G2, G3 and G7 (runs without a fallback are unchanged:
# the earlier rule and the current one give identical results). Part B's rows (E14 to E16, E18 to E20, E22 to E25, E25e,
# NA3, E11c, the row of ruling (7), G4, G5, G8, G9; E16 and NA3 were run failing-first in R1a) follow at the end of the
# file, written on a base holding KIT-224 and KIT-190; item vii (E21, G6) moved to KIT-255 (VAL-283). None of these rows
# runs the correlation filter (KIT-255's kernels): they assume nothing about which sample columns exist after it.
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


kit242Rule <- function(rule) list(MultiScholaR.normalisation_fallback = rule)

#' The assay's matrix in the table's own row order (the order the method receives).
kit242TableMatrix <- function(object, binding, assay, runs = kit242Runs()) {
  table <- as.data.frame(methods::slot(object, binding$data_slot)[[assay]], check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(table[, intersect(runs, names(table)), drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(table[[methods::slot(object, binding$id_slot)]])
  m
}

#' Evaluate `expr`, keeping the conditions of the given classes (muffled) and the plain warnings' texts.
kit242Conditions <- function(expr, classes) {
  kept <- list()
  plain <- character()
  value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    hit <- intersect(class(w), classes)
    if (length(hit)) kept[[length(kept) + 1L]] <<- w else plain <<- c(plain, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, kept = kept, plain = plain)
}

kit242Normaliser <- function(method) switch(method, cyclicloess = limma::normalizeCyclicLoess, quantile = limma::normalizeQuantiles,
  scale = limma::normalizeMedianAbsValues)

kit242ProteinObject <- function(m) {
  table <- data.frame(Protein.Ids = rownames(m), as.data.frame(m, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(table) <- NULL
  methods::new("ProteinQuantitativeData", protein_quant_table = table, protein_id_column = "Protein.Ids",
    design_matrix = kit172Design(), sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group", args = list())
}

kit242ProteinValues <- function(k) {
  set.seed(242000L + k, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  runs <- kit242Runs()
  y <- outer(stats::rnorm(60L, 20, 1.5), rep(1, length(runs))) + outer(rep(1, 60L), stats::rnorm(length(runs), 0, 0.5)) +
    matrix(stats::rnorm(60L * length(runs), 0, 0.25), 60L)
  dimnames(y) <- list(sprintf("P%03d", 1:60), runs)
  y
}

#' The empty-sample object of E1, E7, E8, E11, E13 and NA1: A2's sample B_6 has no value (log2 scale).
kit242EmptySampleObject <- function(class) {
  v <- kit242Draw(1L)
  v$A2[, "B_6"] <- NA_real_
  kit242Object(lapply(v, log2), class)
}

#' The normalisation tab (helper-workflow-omics-worker.R's omicsTabServer), managed or in process, on `object`.
kit242Tab <- function(class, object, managed, inputs) {
  omic <- if (identical(class, "lipid")) "lipidomics" else "metabolomics"
  binding <- kit172Binding(class)
  fixture <- kit242Owner(object, class)
  owner <- fixture$owner
  probe <- new.env()
  notes <- new.env()
  notes$updates <- list(); notes$notes <- list(); notes$modals <- integer(); notes$choices <- list()
  res <- new.env()
  local({
    guiEntryInstallBrowserDouble(notes)
    if (managed) omicsManagedChainMocks(fixture, probe)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(omicsTabServer(fixture, omic, probe, managed), {
      session$userData$probe_owner <- owner
      base <- omicsTabInputs()
      base[names(inputs)] <- inputs
      do.call(session$setInputs, base)
      commitWorkflowOmicsStoreStep(owner, paste0(binding$prefix, "_qc_finalize"), paste0(binding$prefix, "_qc_complete"), list(),
        "QC processing complete - ready for normalization")
      session$setInputs(run_normalization = 1)
      for (i in seq_len(10L)) {
        session$elapse(200)
        session$flushReact()
      }
      res$normalization_complete <- shiny::isolate(session$userData$probe_norm$normalization_complete)
      res$ruv_complete <- shiny::isolate(session$userData$probe_norm$ruv_complete)
      res$ruv_results <- shiny::isolate(session$userData$probe_norm$ruv_optimization_results)
    }))))
  })
  list(log = as.character(unlist(probe$log)), notes = notes$notes, normalization_complete = res$normalization_complete,
    ruv_complete = res$ruv_complete, ruv_results = res$ruv_results)
}

for (class in c("lipid", "metab")) {
  binding <- kit172Binding(class)
  prefix <- binding$prefix

  # ---- E1: a sample with no value in one assay (OD-1 (c)); today the whole assay is left unnormalised ----------------
  for (method in c("cyclicloess", "quantile", "scale")) {
    test_that(paste("KIT-242 E1:", class, method, "normalises A2 without the empty sample and leaves it empty in A2 only"), {
      object <- kit242EmptySampleObject(class)
      run <- kit242Conditions(normaliseBetweenSamples(object, normalisation_method = method), "omics_normalisation_sample_excluded")
      expect_length(run$kept, 1L)
      if (length(run$kept)) {
        expect_identical(run$kept[[1L]]$assay, "A2")
        expect_identical(run$kept[[1L]]$sample, "B_6")
        expect_identical(run$kept[[1L]]$method, method)
      }
      out <- kit242TableMatrix(run$value, binding, "A2")
      input <- kit242TableMatrix(object, binding, "A2")
      expect_true(all(is.na(out[, "B_6"])))
      others <- setdiff(colnames(input), "B_6")
      expect_equal(unname(out[, others]), unname(kit242Normaliser(method)(input[, others, drop = FALSE])), tolerance = 0)
      expect_equal(unname(kit242TableMatrix(run$value, binding, "A1")), unname(kit242Normaliser(method)(kit242TableMatrix(object, binding, "A1"))),
        tolerance = 0)
      expect_true("B_6" %in% as.character(run$value@design_matrix$Run))
    })
  }

  # ---- E2: a method failure stops, under every option (G0 C2: today a warning and the assay unnormalised) -----------
  test_that(paste("KIT-242 E2:", class, "quantile on a sample holding one value stops with omics_normalisation_failed"), {
    v <- kit242Draw(2L)
    keep <- v$A2[1L, "B_6"]
    v$A2[, "B_6"] <- NA_real_
    v$A2[1L, "B_6"] <- keep
    object <- kit242Object(v, class)
    err <- tryCatch(suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "quantile"))),
      error = function(e) e)
    expect_s3_class(err, "omics_normalisation_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$method, "quantile")
    expect_identical(err$reason, "normaliser_error")
    expect_s3_class(err$parent, "error")
    expect_match(if (inherits(err$parent, "condition")) conditionMessage(err$parent) else "", "two non-NA values")
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
    err <- tryCatch(suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"))),
      error = function(e) e)
    expect_s3_class(err, "omics_normalisation_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$reason, "malformed_assay")
  })

  # ---- E4, E5: ITSD with no standard matched (OD-2 (b)); today applied TRUE and nothing per assay ------------------
  test_that(paste("KIT-242 E4:", class, "ITSD with no standard matched in A2: not applied there, recorded with 0 standards"), {
    v <- kit242Draw(4L)
    object <- kit242Object(v, class, list(A2 = ifelse(startsWith(rownames(v$A2), "IS_"), "ISTD", "PC")))
    run <- kit242Conditions(normaliseUntransformedData(object, method = "ITSD", itsd_aggregation = "median"), "omics_itsd_not_applied")
    expect_length(run$kept, 1L)
    if (length(run$kept)) {
      expect_identical(run$kept[[1L]]$assay, "A2")
      expect_identical(run$kept[[1L]]$reason, "no_standard_matched")
    }
    itsd <- run$value@args$ITSDNormalization
    expect_true(isTRUE(itsd$applied))
    expect_identical(itsd$applied_per_assay, list(A1 = TRUE, A2 = FALSE))
    expect_identical(itsd$not_applied_reason_per_assay, list(A2 = "no_standard_matched"))
    expect_identical(lapply(itsd$itsd_counts_per_assay, as.integer), list(A1 = 5L, A2 = 0L))
    expect_identical(names(itsd$itsd_features_per_assay), "A1")
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(object, class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet(commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"))
    step <- kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config))
    expect_identical(step$value$itsd$not_applied_reason_per_assay, list(A2 = "no_standard_matched"))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_itsd_norm"))
    expect_true(isTRUE(outcomes$itsd_applied))
    expect_identical(outcomes$itsd_applied_per_assay, list(A1 = TRUE, A2 = FALSE))
    expect_identical(lapply(outcomes$itsd_counts_per_assay, as.integer), list(A1 = 5L, A2 = 0L))
  })

  test_that(paste("KIT-242 E5:", class, "ITSD with no standard in any assay: not applied, recorded FALSE"), {
    v <- kit242Draw(5L)
    unmatched <- lapply(v, function(m) ifelse(startsWith(rownames(m), "IS_"), "ISTD", "PC"))
    object <- kit242Object(v, class, unmatched)
    run <- kit242Conditions(normaliseUntransformedData(object, method = "ITSD", itsd_aggregation = "median"), "omics_itsd_not_applied")
    expect_length(run$kept, 2L)
    itsd <- run$value@args$ITSDNormalization
    expect_false(isTRUE(itsd$applied))
    expect_identical(itsd$applied_per_assay, list(A1 = FALSE, A2 = FALSE))
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(object, class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet(commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"))
    kit172Quiet(runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_itsd_norm"))
    expect_false(isTRUE(outcomes$itsd_applied))
    expect_identical(outcomes$itsd_applied_per_assay, list(A1 = FALSE, A2 = FALSE))
  })

  # ---- E6: a sample without standard signal (OD-3 (d)); today the whole assay then loses its normalisation ----------
  test_that(paste("KIT-242 E6:", class, "a sample without standard signal is missing in A2 only and A2 is normalised without it"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(6L)
    v$A2[startsWith(rownames(v$A2), "IS_"), "B_6"] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet({
      commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
      itsd <- runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)
      runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
      runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
    })
    expect_identical(itsd$value$itsd$samples_without_standard_per_assay, list(A2 = "B_6"))
    itsd_outcomes <- kit242Outcomes(owner, paste0(prefix, "_itsd_norm"))
    expect_identical(itsd_outcomes$itsd_samples_without_standard_per_assay, list(A2 = "B_6"))
    before <- kit242TableMatrix(owner$getState(paste0(prefix, "_log2")), binding, "A2")
    after <- kit242TableMatrix(owner$getState(paste0(prefix, "_normalized")), binding, "A2")
    others <- setdiff(colnames(before), "B_6")
    expect_true(all(is.na(after[, "B_6"])))
    expect_equal(unname(after[, others]), unname(limma::normalizeCyclicLoess(before[, others, drop = FALSE])), tolerance = 0)
    a1 <- kit242TableMatrix(owner$getState(paste0(prefix, "_normalized")), binding, "A1")
    expect_identical(sum(is.finite(a1[, "B_6"])), nrow(a1))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_normalized"))
    expect_identical(outcomes$normalisation_samples_excluded_per_assay, list(A2 = "B_6"))
  })

  # ---- E7: the declaration's product and the record per assay (OD-1 (c)) ------------------------------------------
  test_that(paste("KIT-242 E7:", class, "the between-sample step records per assay what it normalised and left out"), {
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(kit242EmptySampleObject(class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet(commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete"))
    run <- kit172Quiet(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config))
    expect_identical(run$value$per_assay$samples_excluded, list(A2 = "B_6"))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_normalized"))
    expect_identical(outcomes$normalisation_applied_per_assay, list(A1 = TRUE, A2 = TRUE))
    expect_identical(outcomes$normalisation_samples_excluded_per_assay, list(A2 = "B_6"))
  })

  # ---- E8: the tab names the assay and the sample (managed and in process) -----------------------------------------
  for (managed in c(TRUE, FALSE)) {
    test_that(paste("KIT-242 E8:", class, if (managed) "managed" else "in process", "tab names A2 and B_6 and warns"), {
      skip_if_not(workflowOwnerCapable())
      tab <- kit242Tab(class, kit242EmptySampleObject(class), managed,
        list(apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "skip", log_offset = 0))
      expect_true(isTRUE(tab$normalization_complete))
      expect_true(any(grepl("Assay A2 was normalised without samples B_6", tab$log, fixed = TRUE)), info = paste(tab$log, collapse = " | "))
      types <- vapply(tab$notes, function(note) as.character(note$type), "")
      messages <- vapply(tab$notes, function(note) as.character(note$message), "")
      expect_true(any(types == "warning" & grepl("B_6", messages, fixed = TRUE)))
    })
  }

  # ---- E9: the session summary's ITSD line per assay (outcome 6) ----------------------------------------------------
  test_that(paste("KIT-242 E9:", class, "the session summary names ITSD per assay; the session field keeps the checkbox"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(9L)
    fixture <- kit242Owner(kit242Object(v, class, list(A2 = ifelse(startsWith(rownames(v$A2), "IS_"), "ISTD", "PC"))), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit172Quiet({
      commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
      runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)
    })
    omic_type <- if (identical(class, "lipid")) "lipidomics" else "metabolomics"
    export <- omicGuiConfig(omic_type)$norm$export
    input <- list(apply_itsd = TRUE, itsd_aggregation = "median", norm_method = "cyclicloess", ruv_mode = "skip",
      min_pearson_correlation_threshold = 0.5, ruv_grouping_variable = "group", log_offset = 0)
    session <- kit172Quiet(collectWorkflowOmicsNormSession(binding, fixture$workflow, list(), input, "kit242", export))
    expect_true(isTRUE(session$itsd_applied))
    dir <- withr::local_tempdir()
    written <- new.env()
    kit172Quiet(saveWorkflowOmicsNormSessionArtifacts(binding, session, dir, list(request_id = "dia-export-kit242", timestamp = "20261125_120000"),
      export, writeLinesFn = function(text, path) {
        written$text <- text
        writeLines(text, path)
      }))
    expect_match(paste(written$text, collapse = "\n"), "- ITSD applied: A1 yes (5 standards); A2 no (no internal standard matched)", fixed = TRUE)
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

  # ---- E11: older records (OD-4 (c)): the normalisation marker and its earlier rule -----------------------------------
  test_that(paste("KIT-242 E11:", class, "the normalisation marker: current 'exclude', earlier 'warn' reproduces the earlier result, and replays"), {
    marker <- workflowReplayRuleMarkers$normalisation_fallback
    expect_identical(marker$current, "exclude")
    expect_identical(marker$earlier, "warn")
    expect_identical(workflowOperationEntry(paste0(prefix, "_between_sample_normalisation"))$declaration$fixed$normalisation_fallback, "exclude")
    expect_identical(workflowReplayRuleMarkers$protein_normalisation_fallback$current, "remove")
    # (ii) under the earlier token the earlier result: A2 unchanged, with the earlier warning.
    object <- kit242EmptySampleObject(class)
    earlier <- withr::with_options(kit242Rule("warn"), kit242Conditions(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"),
      "omics_normalisation_sample_excluded"))
    expect_length(earlier$kept, 0L)
    expect_true(any(grepl("At least one sample column contains only NA values. Skipping normalization.", earlier$plain, fixed = TRUE)))
    expect_equal(unname(kit242TableMatrix(earlier$value, binding, "A2")), unname(kit242TableMatrix(object, binding, "A2")), tolerance = 0)
    # (iii) a record made under each rule replays with no mismatch (the earlier one without the marker).
    skip_if_not(workflowOwnerCapable())
    harness <- new.env(parent = globalenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    for (rule in c("warn", "exclude")) {
      root <- withr::local_tempdir()
      fixture <- kit172Quiet(kit172Owner(root, object, class))
      owner <- fixture$owner
      config <- fixture$workflow$config_list
      withr::with_options(kit242Rule(if (identical(rule, "warn")) "warn" else NULL), kit172Quiet({
        kit242Prepare(binding, owner, config)
        runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
      }))
      json <- file.path(root, "record.json")
      kit172Quiet(writeWorkflowAnalysisRecord(owner, json))
      history <- unlist(owner$getHistory())
      closeDIAStateAdapter(owner)
      record <- jsonlite::read_json(json, simplifyVector = FALSE)
      step <- Filter(function(s) identical(s$state, paste0(prefix, "_normalized")), record$steps)[[1L]]
      has_marker <- any(vapply(step$parameters, function(p) identical(p$name, "normalisation_fallback"), logical(1)))
      expect_identical(has_marker, identical(rule, "exclude"), info = rule)
      source_dir <- fixture$source_dir
      copy <- withr::local_tempdir()
      stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
      store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
      copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
      operations <- vapply(record$steps, function(s) s$operation, "")
      states <- vapply(record$steps, function(s) s$state, "")
      prepared <- harness$kitPrepareReplayProject(copy_root, states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))], api)
      files <- c(list.files(file.path(source_dir, "import_sources"), full.names = TRUE), file.path(source_dir, "design_matrix.tab"))
      sums <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
      inputs <- character()
      for (input in record$inputs) {
        hit <- files[sums == (input$sha256 %||% "")]
        if (length(hit)) inputs[[input$role]] <- hit[[1L]]
      }
      replay <- suppressWarnings(suppressMessages(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths,
        withr::local_tempdir(), api)))
      expect_identical(replay$states, history, info = rule)
      expect_identical(as.character(replay$mismatches), character(), info = rule)
    }
  })

  # ---- E3b (r1.6, the coordinator's ruling on item 2): duplicate feature ids are a malformed assay and stop ----------
  test_that(paste("KIT-242 E3b:", class, "an assay with duplicate feature ids stops, naming the cause and 'Resolve Duplicates'"), {
    v <- kit242Draw(32L)
    object <- kit242Object(lapply(v, log2), class)
    tables <- methods::slot(object, binding$data_slot)
    tables$A2 <- rbind(tables$A2, tables$A2[1:2, , drop = FALSE])
    methods::slot(object, binding$data_slot, check = FALSE) <- tables
    err <- tryCatch(suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"))),
      error = function(e) e)
    expect_s3_class(err, "omics_normalisation_failed")
    if (inherits(err, "omics_normalisation_failed")) {
      expect_identical(err$assay, "A2")
      expect_identical(err$reason, "malformed_assay")
      expect_match(conditionMessage(err), "2 duplicated feature identifiers \\(A2_F001, A2_F002\\)")
      expect_match(conditionMessage(err), "Run 'Resolve Duplicates' in the Quality Control tab", fixed = TRUE)
    }
  })

  # ---- E11b (r1.6, the coordinator's ruling on item 1): the ITSD record's marker ------------------------------------
  test_that(paste("KIT-242 E11b:", class, "an older record replays the ITSD step exactly; a new one records the per-assay fields"), {
    marker <- workflowReplayRuleMarkers$itsd_outcome
    expect_identical(marker$current, "per_assay")
    expect_identical(marker$earlier, "applied_flag")
    expect_identical(workflowOperationEntry(paste0(prefix, "_itsd_normalisation"))$declaration$fixed$itsd_outcome, "per_assay")
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(33L)
    object <- kit242Object(v, class, list(A2 = ifelse(startsWith(rownames(v$A2), "IS_"), "ISTD", "PC")))
    harness <- new.env(parent = globalenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    earlier <- list(MultiScholaR.itsd_outcome = "applied_flag", MultiScholaR.normalisation_fallback = "warn")
    for (rule in c("earlier", "current")) {
      root <- withr::local_tempdir()
      fixture <- kit172Quiet(kit172Owner(root, object, class))
      owner <- fixture$owner
      config <- fixture$workflow$config_list
      withr::with_options(if (identical(rule, "earlier")) earlier else list(MultiScholaR.itsd_outcome = NULL), kit172Quiet({
        commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
        runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)
        runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
        runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
      }))
      itsd <- owner$getState(paste0(prefix, "_itsd_norm"))@args$ITSDNormalization
      outcomes <- kit242Outcomes(owner, paste0(prefix, "_itsd_norm"))
      if (identical(rule, "earlier")) {
        expect_null(itsd$applied_per_assay)
        expect_true(isTRUE(itsd$applied))
        expect_true(isTRUE(outcomes$itsd_applied))
        expect_null(outcomes$itsd_applied_per_assay)
      } else {
        expect_identical(itsd$applied_per_assay, list(A1 = TRUE, A2 = FALSE))
        expect_identical(outcomes$itsd_applied_per_assay, list(A1 = TRUE, A2 = FALSE))
      }
      json <- file.path(root, "record.json")
      kit172Quiet(writeWorkflowAnalysisRecord(owner, json))
      history <- unlist(owner$getHistory())
      closeDIAStateAdapter(owner)
      record <- jsonlite::read_json(json, simplifyVector = FALSE)
      step <- Filter(function(s) identical(s$state, paste0(prefix, "_itsd_norm")), record$steps)[[1L]]
      has_marker <- any(vapply(step$parameters, function(p) identical(p$name, "itsd_outcome"), logical(1)))
      expect_identical(has_marker, identical(rule, "current"), info = rule)
      source_dir <- fixture$source_dir
      copy <- withr::local_tempdir()
      stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
      store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
      copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
      operations <- vapply(record$steps, function(s) s$operation, "")
      states <- vapply(record$steps, function(s) s$state, "")
      prepared <- harness$kitPrepareReplayProject(copy_root, states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))], api)
      files <- c(list.files(file.path(source_dir, "import_sources"), full.names = TRUE), file.path(source_dir, "design_matrix.tab"))
      sums <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
      inputs <- character()
      for (input in record$inputs) {
        hit <- files[sums == (input$sha256 %||% "")]
        if (length(hit)) inputs[[input$role]] <- hit[[1L]]
      }
      replay <- suppressWarnings(suppressMessages(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths,
        withr::local_tempdir(), api)))
      expect_identical(replay$states, history, info = rule)
      expect_identical(as.character(replay$mismatches), character(), info = rule)
    }
  })

  # ---- E13: KIT-145's guard counts the other samples (item viii); today it refuses with fewest 0 --------------------
  test_that(paste("KIT-242 E13:", class, "with the cyclic-loess guard on, the empty sample is left out before the guard counts"), {
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    object <- kit242EmptySampleObject(class)
    run <- tryCatch(kit242Conditions(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"), "omics_normalisation_sample_excluded"),
      error = function(e) e)
    expect_false(inherits(run, "omics_cyclic_loess_too_few_features"))
    expect_false(inherits(run, "error"))
    if (!inherits(run, "error")) expect_length(run$kept, 1L)
  })

  # ---- E26: a sample with no value in any assay leaves the design (item 1a), with a notice ------------------------------
  test_that(paste("KIT-242 E26:", class, "a sample with no value in any assay is removed from the design and every assay"), {
    v <- kit242Draw(26L)
    v$A1[, "B_6"] <- NA_real_
    v$A2[, "B_6"] <- NA_real_
    object <- kit242Object(lapply(v, log2), class)
    run <- kit242Conditions(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"),
      c("omics_normalisation_sample_excluded", "omics_normalisation_sample_removed"))
    removed <- Filter(function(w) inherits(w, "omics_normalisation_sample_removed"), run$kept)
    expect_length(removed, 1L)
    if (length(removed)) expect_identical(removed[[1L]]$sample, "B_6")
    expect_false("B_6" %in% as.character(run$value@design_matrix$Run))
    for (assay in c("A1", "A2")) expect_false("B_6" %in% names(methods::slot(run$value, binding$data_slot)[[assay]]), info = assay)
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
    # r1.5 (K242-RC4-N06): under OD-1 (c) the step completes.
    expect_identical(owner$current_state, paste0(prefix, "_normalized"))
    before <- owner$getState(paste0(prefix, "_log2"))
    after <- owner$getState(owner$current_state)
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
    # r1.5 (K242-RC4-N06): under OD-2 (b) the step completes.
    expect_identical(owner$current_state, paste0(prefix, "_itsd_norm"))
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

  # ---- G1: complete objects where ITSD acts: the earlier and the current rule give identical results ---------------
  test_that(paste("KIT-242 G1:", class, "with no fallback, ITSD, normalisation and the DA are identical under the earlier and the current rule"), {
    skip_if_not(workflowOwnerCapable())
    object <- kit242Object(kit242Draw(41L), class)
    journey <- function(rule) {
      fixture <- kit242Owner(object, class, env = parent.frame())
      owner <- fixture$owner
      config <- fixture$workflow$config_list
      da <- withr::with_options(kit242Rule(rule), kit172Quiet({
        commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
        runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)
        runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
        runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
        runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner),
          list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts())
      }))
      itsd <- owner$getState(paste0(prefix, "_itsd_norm"))
      list(itsd_rows = vapply(c("A1", "A2"), function(a) sum(startsWith(rownames(kit242Matrix(itsd, binding, a)), "IS_")), integer(1)),
        counts = lapply(itsd@args$ITSDNormalization$itsd_counts_per_assay, as.integer),
        normalised = lapply(c(A1 = "A1", A2 = "A2"), function(a) kit242Matrix(owner$getState(paste0(prefix, "_normalized")), binding, a)),
        da = as.data.frame(da$results[[binding$long_table]]))
    }
    current <- journey(NULL)
    earlier <- journey("warn")
    expect_identical(unname(current$itsd_rows), c(0L, 0L)) # ITSD acted: no standard left
    expect_identical(current$counts, list(A1 = 5L, A2 = 5L))
    expect_identical(current$normalised, earlier$normalised)
    expect_identical(current$da, earlier$da)
  })

  # ---- G2: method none ----------------------------------------------------------------------------------------------
  test_that(paste("KIT-242 G2:", class, "method none leaves every assay as it is, with no KIT-242 condition"), {
    object <- kit242EmptySampleObject(class)
    run <- kit242Conditions(normaliseBetweenSamples(object, normalisation_method = "none"),
      c("omics_normalisation_sample_excluded", "omics_normalisation_sample_removed"))
    expect_length(run$kept, 0L)
    for (assay in c("A1", "A2")) {
      expect_equal(unname(kit242TableMatrix(run$value, binding, assay)), unname(kit242TableMatrix(object, binding, assay)), tolerance = 0)
    }
  })

  # ---- G7: KIT-145's guard on, no empty sample: unchanged ---------------------------------------------------------
  test_that(paste("KIT-242 G7:", class, "with the guard on and no empty sample, the result equals the guard off"), {
    object <- kit242Object(lapply(kit242Draw(47L), log2), class)
    off <- suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "cyclicloess")))
    on <- withr::with_options(list(MultiScholaR.cyclic_loess_minimum = TRUE),
      suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = "cyclicloess"))))
    for (assay in c("A1", "A2")) {
      expect_identical(kit242TableMatrix(on, binding, assay), kit242TableMatrix(off, binding, assay), info = assay)
    }
  })
}

# ---- E17 and G3: Part C, proteomics normalisation with a sample with no value -------------------------------------------
for (method in c("cyclicloess", "quantile", "scale")) {
  test_that(paste("KIT-242 E17: proteomics", method, "removes a sample with no value from the table and the design, the others normalised without it"), {
    y <- kit242ProteinValues(17L)
    y[, "B_6"] <- NA_real_
    run <- kit242Conditions(normaliseBetweenSamples(kit242ProteinObject(y), normalisation_method = method),
      "omics_normalisation_sample_removed")
    expect_length(run$kept, 1L)
    if (length(run$kept)) expect_identical(run$kept[[1L]]$sample, "B_6")
    table <- run$value@protein_quant_table
    expect_false("B_6" %in% names(table))
    expect_false("B_6" %in% as.character(run$value@design_matrix$Run))
    others <- setdiff(colnames(y), "B_6")
    out <- as.matrix(table[match(rownames(y), table$Protein.Ids), others])
    expect_equal(unname(out), unname(kit242Normaliser(method)(y[, others])), tolerance = 0)
  })

  test_that(paste("KIT-242 G3: proteomics", method, "on a table with no empty sample is identical under the earlier and the current rule"), {
    object <- kit242ProteinObject(kit242ProteinValues(43L))
    current <- suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = method)))
    earlier <- withr::with_options(list(MultiScholaR.protein_normalisation_fallback = "keep"),
      suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = method))))
    expect_identical(current@protein_quant_table, earlier@protein_quant_table)
    expect_identical(current@design_matrix, earlier@design_matrix)
  })
}

test_that("KIT-242 E17: proteomics quantile on a sample holding one value stops with omics_normalisation_failed", {
  y <- kit242ProteinValues(18L)
  keep <- y[1L, "B_6"]
  y[, "B_6"] <- NA_real_
  y[1L, "B_6"] <- keep
  err <- tryCatch(suppressWarnings(suppressMessages(normaliseBetweenSamples(kit242ProteinObject(y), normalisation_method = "quantile"))),
    error = function(e) e)
  expect_s3_class(err, "omics_normalisation_failed")
  expect_identical(err$reason, "normaliser_error")
})

# ---- (r1.6, the coordinator's ruling on item 4): the proteomics normalisation text is byte-identical when nothing is removed
test_that("KIT-242 Part C: the protein normalisation text equals the earlier text when no sample is removed", {
  earlier <- function(skipped, protein_count, samples, method) paste(
    if (skipped) "Between-Sample Normalisation Skipped (none)\n" else "Between-Sample Normalisation Applied Successfully\n",
    "================================================\n",
    sprintf("Proteins remaining: %d\n", protein_count), sprintf("Samples: %d\n", samples),
    sprintf("Method: %s\n", method),
    if (skipped) "The quantified log2 matrix is carried unchanged (the limpa script's choice).\n" else "",
    "State saved as: 'protein_normalised'\n")
  expect_identical(diaProteinNormalisationText(FALSE, 60L, 12L, "cyclicloess"), earlier(FALSE, 60L, 12L, "cyclicloess"))
  expect_identical(diaProteinNormalisationText(TRUE, 60L, 12L, "none"), earlier(TRUE, 60L, 12L, "none"))
  expect_match(diaProteinNormalisationText(FALSE, 60L, 11L, "scale", "B_6"), "Warning: samples B_6 have no value and were removed from the design.",
    fixed = TRUE)
})

# ==== Part B (ticket r1.6; the coordinator's rulings of 2026-11-26): legacy RUV-III-C and the declared design ===========
# Rows E14 to E16, E18 to E20, E22 to E25, E25e (VAL-290), NA3, the marker row E11c, the row of ruling (7), and the guards
# G4, G5, G8 and G9. Rules (VAL-242 OD-1 (c), OD-3 (d), Q6.1, Q6.2 (c); VAL-290): in each assay a sample or a feature with
# no value is left out of the search and the correction and put back as missing; a feature the control vector does not
# name is not a control; an assay with no usable controls, or for which no number of factors is supported by the data,
# keeps its values without RUV, with a notice and a record; a failed computation and an invalid argument stop
# (`omics_ruv_failed`); no assay is dropped. Seeds 242000 + k as above.

kit242RuvRule <- function(rule) list(MultiScholaR.ruv_fallback = rule)

#' The object's tables as tibbles (the legacy RUV method's input type; the owner's states are tibbles).
kit242Tibbles <- function(object, binding) {
  methods::slot(object, binding$data_slot) <- lapply(methods::slot(object, binding$data_slot), tibble::as_tibble)
  object
}

#' Intensities without standards, A2's sample B_6 with no value (the object of E18 to E20; G0 C4 after ITSD).
kit242RuvEmptySampleValues <- function(k) {
  v <- kit242Draw(k, n_is = c(A1 = 0L, A2 = 0L))
  v$A2[, "B_6"] <- NA_real_
  v
}

#' C12 a-decl's object: every A2 feature misses one value, so no A2 feature is observed in every sample.
kit242RuvIncompleteValues <- function(k) {
  v <- kit242Draw(k, n_is = c(A1 = 0L, A2 = 0L))
  runs <- kit242Runs()
  for (j in seq_len(nrow(v$A2))) v$A2[j, runs[((j - 1L) %% length(runs)) + 1L]] <- NA_real_
  v
}

#' QC, log2 and between-sample normalisation through the declarations.
kit242RuvChain <- function(binding, owner, config, method = "cyclicloess") {
  kit242Prepare(binding, owner, config)
  kit172Quiet(runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, method, config))
}

#' Legacy RUV on an object as the declaration's kernel computes it with nothing to leave out: the search, the bound on k,
#' then `ruvIII_C_Varying` (the reference "legacy RUV on a one-assay object of A2 without B_6" of E18 and E19).
kit242LegacyRuvReference <- function(binding, object, parameters) {
  results <- binding$optimise_ruv(object, parameters$ruv_mode, parameters[setdiff(names(parameters), "ruv_mode")])
  k <- omicsRuvKBound(object, results, parameters$ruv_mode, binding$extract_best_k(results))$best_k
  list(k = k, object = ruvIII_C_Varying(object, ruv_grouping_variable = parameters$ruv_grouping_variable, ruv_number_k = k,
    ctrl = binding$extract_controls(results)))
}

#' One assay of an object without some sample columns (the object's other assays removed).
kit242OneAssayWithout <- function(object, binding, assay, samples) {
  tables <- methods::slot(object, binding$data_slot)
  methods::slot(object, binding$data_slot) <- stats::setNames(list(tables[[assay]][, setdiff(names(tables[[assay]]), samples)]), assay)
  object
}

kit242ManualRuv <- function() utils::modifyList(kit172LegacyRuv("manual"), list(manual_k = 1L, manual_percentage = 10))

#' Write the owner's record, close it, and replay the record on a copy of its project (E11's procedure).
kit242ReplayRecord <- function(fixture, root, harness, api) {
  owner <- fixture$owner
  json <- file.path(root, "record.json")
  kit172Quiet(writeWorkflowAnalysisRecord(owner, json))
  history <- unlist(owner$getHistory())
  closeDIAStateAdapter(owner)
  record <- jsonlite::read_json(json, simplifyVector = FALSE)
  source_dir <- fixture$source_dir
  copy <- withr::local_tempdir(.local_envir = parent.frame())
  stopifnot(system2("cp", c("-a", file.path(source_dir, "."), copy)) == 0L)
  store_root <- dirname(dirname(list.files(source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]]))
  copy_root <- file.path(copy, substring(store_root, nchar(source_dir) + 2L))
  operations <- vapply(record$steps, function(s) s$operation, "")
  states <- vapply(record$steps, function(s) s$state, "")
  prepared <- harness$kitPrepareReplayProject(copy_root, states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))], api)
  files <- c(list.files(file.path(source_dir, "import_sources"), full.names = TRUE), file.path(source_dir, "design_matrix.tab"))
  sums <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
  inputs <- character()
  for (input in record$inputs) {
    hit <- files[sums == (input$sha256 %||% "")]
    if (length(hit)) inputs[[input$role]] <- hit[[1L]]
  }
  replay <- suppressWarnings(suppressMessages(harness$kitReplayJourney(record, prepared, store_root, inputs, fixture$paths,
    withr::local_tempdir(.local_envir = parent.frame()), api)))
  list(record = record, history = history, replay = replay)
}

for (class in c("lipid", "metab")) {
  binding <- kit172Binding(class)
  prefix <- binding$prefix

  # ---- E14: no control found in A2 (Q6.2 (c), no usable controls); today A2 is dropped -------------------------------
  test_that(paste("KIT-242 E14:", class, "legacy RUV-III-C with no control found in A2 keeps A2 uncorrected, with a classed warning"), {
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(14L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    run <- kit242Conditions(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = list(A1 = 1L, A2 = 1L),
      ctrl = list(A1 = sprintf("A1_F%03d", 1:8), A2 = c("absent_1", "absent_2"))), "omics_ruv_not_applied")
    expect_identical(names(methods::slot(run$value, binding$data_slot)), c("A1", "A2"))
    expect_length(run$kept, 1L)
    if (length(run$kept)) {
      expect_identical(run$kept[[1L]]$assay, "A2")
      expect_identical(run$kept[[1L]]$reason, "no_usable_controls")
    }
    a2 <- kit242Matrix(run$value, binding, "A2")
    expect_false(is.null(a2))
    if (!is.null(a2)) expect_equal(a2, kit242Matrix(object, binding, "A2"), tolerance = 0)
    expect_false(kit242Equal(kit242Matrix(object, binding, "A1"), kit242Matrix(run$value, binding, "A1")))
  })

  # ---- E15: no control in any assay; today the input object is returned with no condition ---------------------------
  test_that(paste("KIT-242 E15:", class, "legacy RUV-III-C with no control in any assay keeps both, each with a classed warning"), {
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(15L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    run <- kit242Conditions(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = list(A1 = 1L, A2 = 1L),
      ctrl = list(A1 = c("absent_1", "absent_2"), A2 = c("absent_3", "absent_4"))), "omics_ruv_not_applied")
    expect_identical(vapply(run$kept, function(w) w$assay, ""), c("A1", "A2"))
    for (assay in c("A1", "A2")) expect_equal(kit242Matrix(run$value, binding, assay), kit242Matrix(object, binding, assay), tolerance = 0)
  })

  # ---- E16: an invalid k stops (an argument; every option); today the assay is dropped with a logger line -------------
  test_that(paste("KIT-242 E16:", class, "legacy RUV-III-C with k = -1 for A2 stops with omics_ruv_failed (invalid_k)"), {
    values <- lapply(kit242Draw(16L, n_is = c(A1 = 0L, A2 = 0L)), log2)
    object <- kit242Tibbles(kit242Object(values, class), binding)
    controls <- list(A1 = sprintf("A1_F%03d", 1:8), A2 = sprintf("A2_F%03d", 1:8))
    err <- expect_error(suppressWarnings(suppressMessages(ruvIII_C_Varying(object, ruv_grouping_variable = "group",
      ruv_number_k = list(A1 = 1L, A2 = -1L), ctrl = controls))), class = "omics_ruv_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$reason, "invalid_k")
  })

  # ---- E18: legacy Automatic, A2 with B_6 empty (item 1a i to iii); today A2 is dropped (G0 C4) -----------------------
  test_that(paste("KIT-242 E18:", class, "legacy Automatic RUV corrects A2 without B_6, keeps B_6 empty there and in A1 keeps it"), {
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(kit242Object(kit242RuvEmptySampleValues(18L), class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)), error = function(e) e)
    expect_false(inherits(run, "error"))
    expect_identical(owner$current_state, paste0(prefix, "_ruv_corrected"))
    after <- owner$getState(owner$current_state)
    a2 <- kit242Matrix(after, binding, "A2")
    expect_false(is.null(a2))
    expect_true(all(is.finite(kit242Matrix(after, binding, "A1")[, "B_6"])))
    if (is.null(a2) || inherits(run, "error")) return(invisible())
    expect_true(all(is.na(a2[, "B_6"])))
    reference <- kit172Quiet(kit242LegacyRuvReference(binding, kit242OneAssayWithout(normalised, binding, "A2", "B_6"),
      kit172LegacyRuv("automatic")))
    ref <- kit242Matrix(reference$object, binding, "A2")
    others <- setdiff(kit242Runs(), "B_6")
    expect_identical(sort(rownames(a2)), sort(rownames(ref)))
    expect_equal(a2[rownames(ref), others], ref[, others], tolerance = 0)
    expect_identical(as.integer(run$value$best_k$A2), as.integer(reference$k$A2))
    expect_identical(run$value$ruv_samples_excluded$A2, "B_6")
  })

  # ---- E19: legacy Manual on the same object; today omics_ruv_insufficient_controls stops the run ---------------------
  test_that(paste("KIT-242 E19:", class, "legacy Manual RUV corrects A2 without B_6 and keeps B_6 empty there"), {
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(kit242Object(kit242RuvEmptySampleValues(19L), class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit242ManualRuv(), config)), error = function(e) e)
    expect_false(inherits(run, "error"))
    if (inherits(run, "error")) return(invisible())
    a2 <- kit242Matrix(owner$getState(paste0(prefix, "_ruv_corrected")), binding, "A2")
    expect_false(is.null(a2))
    if (is.null(a2)) return(invisible())
    expect_true(all(is.na(a2[, "B_6"])))
    reference <- kit172Quiet(kit242LegacyRuvReference(binding, kit242OneAssayWithout(normalised, binding, "A2", "B_6"), kit242ManualRuv()))
    ref <- kit242Matrix(reference$object, binding, "A2")
    others <- setdiff(kit242Runs(), "B_6")
    expect_equal(a2[rownames(ref), others], ref[, others], tolerance = 0)
  })

  # ---- E20: the declared design (items iv to vi); today A2 cannot be given factors ----------------------------------
  test_that(paste("KIT-242 E20:", class, "the declared design computes A2's factors on the 11 samples and leaves B_6 out of A2's DA"), {
    skip_if_not_installed("ruv")
    object <- kit242Tibbles(kit242Object(lapply(kit242RuvEmptySampleValues(20L), log2), class), binding)
    run <- tryCatch(kit172Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user")), error = function(e) e)
    expect_false(inherits(run, "error"))
    if (inherits(run, "error")) return(invisible())
    others <- setdiff(kit242Runs(), "B_6")
    w <- run$factors$A2$W
    expect_identical(rownames(w), others)
    expect_identical(run$factors$A2$samples_excluded, "B_6")
    reference <- kit172Quiet(runOmicsRuvFactors(kit242OneAssayWithout(object, binding, "A2", "B_6"), "groups_all_features", 1L, "user"))
    expect_equal(w, reference$factors$A2$W, tolerance = 0)
    corrected <- kit242Matrix(run$corrected, binding, "A2")
    expect_true(all(is.na(corrected[, "B_6"])))
    expect_equal(corrected[, others], kit242Matrix(reference$corrected, binding, "A2")[rownames(corrected), others], tolerance = 0)
    design <- as.data.frame(kit172Design())
    expr <- kit242TableMatrix(object, binding, "A2")
    da <- applyOmicsRuvFactorsToDaInput(expr, design, "~ 0 + group", "Run", run$factors$A2, "A2", "group")
    ref_da <- applyOmicsRuvFactorsToDaInput(expr[, others], design, "~ 0 + group", "Run", reference$factors$A2, "A2", "group")
    expect_identical(colnames(da$expr_matrix), others)
    expect_equal(da$expr_matrix, ref_da$expr_matrix, tolerance = 0)
    expect_identical(as.character(da$design_matrix$Run), others)
    kept <- as.character(ref_da$design_matrix$Run) %in% others
    expect_equal(da$design_matrix$W1, ref_da$design_matrix$W1[kept], tolerance = 0)
    expect_identical(da$formula_string, ref_da$formula_string)
  })

  # ---- E22: Q6.1, a feature with exactly one value; today the search fails at the control alignment ----------------
  test_that(paste("KIT-242 E22:", class, "a feature with one value is not a control and the search and correction of A2 run"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(22L, n_is = c(A1 = 0L, A2 = 0L))
    v$A2["A2_F001", -1L] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)), error = function(e) e)
    expect_false(inherits(run, "error"))
    if (inherits(run, "error")) return(invisible())
    after <- owner$getState(paste0(prefix, "_ruv_corrected"))
    a2 <- kit242Matrix(after, binding, "A2")
    expect_false(is.null(a2))
    if (!is.null(a2)) expect_false(kit242Equal(kit242Matrix(normalised, binding, "A2"), a2))
    expect_identical(run$value$ruv_status$A2, "corrected")
    controls <- run$optimisation$A2$control_genes_index
    expect_false(isTRUE(unname(controls["A2_F001"])))
  })

  # ---- E23: Q6.1, a feature with no value; today NIPALS refuses and A2 is dropped -------------------------------------
  test_that(paste("KIT-242 E23:", class, "a feature with no value is left out of the search and correction and kept empty"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(23L, n_is = c(A1 = 0L, A2 = 0L))
    v$A2["A2_F001", ] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config, method = "none")
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)), error = function(e) e)
    expect_false(inherits(run, "error"))
    if (inherits(run, "error")) return(invisible())
    a2 <- kit242Matrix(owner$getState(paste0(prefix, "_ruv_corrected")), binding, "A2")
    expect_false(is.null(a2))
    if (is.null(a2)) return(invisible())
    expect_true("A2_F001" %in% rownames(a2))
    if ("A2_F001" %in% rownames(a2)) expect_true(all(is.na(a2["A2_F001", ])))
    expect_false(kit242Equal(kit242Matrix(normalised, binding, "A2"), a2))
    expect_identical(run$value$ruv_features_excluded$A2, "A2_F001")
  })

  # ---- E24: Q6.2 (c), no usable controls (C12 a-decl's object): A2 kept without RUV, recorded and told -----------------
  test_that(paste("KIT-242 E24:", class, "an assay with no usable controls keeps its values without RUV, recorded per assay"), {
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(kit242Object(kit242RuvIncompleteValues(24L), class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)), error = function(e) e)
    expect_false(inherits(run, "error"))
    if (inherits(run, "error")) return(invisible())
    after <- owner$getState(paste0(prefix, "_ruv_corrected"))
    a2 <- kit242Matrix(after, binding, "A2")
    expect_false(is.null(a2))
    if (!is.null(a2)) expect_equal(a2, kit242Matrix(normalised, binding, "A2"), tolerance = 0)
    expect_false(kit242Equal(kit242Matrix(normalised, binding, "A1"), kit242Matrix(after, binding, "A1")))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_ruv_corrected"))
    expect_true(isTRUE(outcomes$ruv_applied))
    expect_identical(outcomes$ruv_applied_per_assay, list(A1 = TRUE, A2 = FALSE))
    expect_true(outcomes$ruv_not_applied_reason_per_assay$A2 %in% c("no_usable_controls", "no_factor_supported"))
    expect_true(is.na(outcomes$k_per_assay$A2))
    expect_match(run$value$ruv_not_applied_text$A2, "RUV-III was not applied to assay A2", fixed = TRUE)
  })

  for (managed in c(TRUE, FALSE)) {
    test_that(paste("KIT-242 E24:", class, if (managed) "managed" else "in process", "tab names A2 as not corrected and warns"), {
      skip_if_not(workflowOwnerCapable())
      tab <- kit242Tab(class, kit242Object(kit242RuvIncompleteValues(24L), class), managed,
        list(apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "automatic", log_offset = 0))
      expect_true(isTRUE(tab$ruv_complete))
      expect_true(any(grepl("RUV-III was not applied to assay A2", tab$log, fixed = TRUE)), info = paste(tab$log, collapse = " | "))
      types <- vapply(tab$notes, function(note) as.character(note$type), "")
      messages <- vapply(tab$notes, function(note) as.character(note$message), "")
      expect_true(any(types == "warning" & grepl("RUV-III was not applied to assay A2", messages, fixed = TRUE)))
    })
  }

  # ---- E25: Q6.2 (c), a failed computation after k and controls are set; today A2 is dropped (`:700`) -----------------
  test_that(paste("KIT-242 E25:", class, "a RUV-III-C error in A2 stops with omics_ruv_failed naming A2 and the step"), {
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(25L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    real <- RUVIIIC::RUVIII_C_Varying
    local_mocked_bindings(RUVIII_C_Varying = function(k, Y, M, toCorrect, potentialControls) {
      if (any(startsWith(colnames(Y), "A2_"))) stop("planted RUV-III-C failure")
      real(k = k, Y = Y, M = M, toCorrect = toCorrect, potentialControls = potentialControls)
    }, .package = "MultiScholaR")
    err <- expect_error(suppressWarnings(suppressMessages(ruvIII_C_Varying(object, ruv_grouping_variable = "group",
      ruv_number_k = list(A1 = 1L, A2 = 1L), ctrl = list(A1 = sprintf("A1_F%03d", 1:8), A2 = sprintf("A2_F%03d", 1:8))))),
      class = "omics_ruv_failed")
    expect_identical(err$assay, "A2")
    expect_identical(err$reason, "correction_error")
  })

  # ---- E25e (VAL-290): no elbow at any percentage: kept without RUV, with the notice; today both assays are dropped ----
  test_that(paste("KIT-242 E25e:", class, "when no number of factors is supported by the data, every assay keeps its values, recorded"), {
    skip_if_not(workflowOwnerCapable())
    local_mocked_bindings(findBestKElbow = function(...) NA_integer_)
    fixture <- kit242Owner(kit242Object(kit242Draw(125L, n_is = c(A1 = 0L, A2 = 0L)), class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)), error = function(e) e)
    expect_false(inherits(run, "error"))
    if (inherits(run, "error")) return(invisible())
    after <- owner$getState(paste0(prefix, "_ruv_corrected"))
    for (assay in c("A1", "A2")) expect_equal(kit242Matrix(after, binding, assay), kit242Matrix(normalised, binding, assay), tolerance = 0)
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_ruv_corrected"))
    expect_false(isTRUE(outcomes$ruv_applied))
    expect_identical(outcomes$ruv_applied_per_assay, list(A1 = FALSE, A2 = FALSE))
    expect_identical(outcomes$ruv_not_applied_reason_per_assay, list(A1 = "no_factor_supported", A2 = "no_factor_supported"))
    expect_match(run$value$ruv_not_applied_text$A1, "no number of unwanted factors was supported by the data", fixed = TRUE)
  })

  test_that(paste("KIT-242 E25e:", class, "in process tab says RUV-III was not applied"), {
    skip_if_not(workflowOwnerCapable())
    local_mocked_bindings(findBestKElbow = function(...) NA_integer_)
    tab <- kit242Tab(class, kit242Object(kit242Draw(125L, n_is = c(A1 = 0L, A2 = 0L)), class), FALSE,
      list(apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "automatic", log_offset = 0))
    expect_true(any(startsWith(tab$log, "RUV-III not applied")), info = paste(tab$log, collapse = " | "))
    expect_false(any(tab$log == "RUV-III correction applied"))
  })

  # ---- ruling (7) (KIT-190's leftover, parked row 48 addendum (3)): the search succeeds, the correction fails ---------
  test_that(paste("KIT-242 ruling 7:", class, "when the correction fails after the search, the tab neither completes RUV nor shows k as removed"), {
    skip_if_not(workflowOwnerCapable())
    local_mocked_bindings(RUVIII_C_Varying = function(k, Y, M, toCorrect, potentialControls) stop("planted RUV-III-C failure"),
      .package = "MultiScholaR")
    tab <- kit242Tab(class, kit242Object(kit242Draw(127L, n_is = c(A1 = 0L, A2 = 0L)), class), FALSE,
      list(apply_itsd = FALSE, norm_method = "cyclicloess", ruv_mode = "automatic", log_offset = 0))
    expect_false(isTRUE(tab$ruv_complete))
    expect_false(any(tab$log == "RUV-III correction applied"))
    for (assay in c("A1", "A2")) {
      summary <- omicNormRuvSummary(tab$ruv_results[[assay]])
      expect_false(startsWith(summary, "Best k") || startsWith(summary, "Factors removed"), info = paste(assay, summary))
    }
  })

  # ---- NA3: the RUV record never names an assay the RUV state lacks, nor "applied" for one left as it was ----------
  test_that(paste("KIT-242 NA3:", class, "every assay the legacy RUV record names is in the RUV state, FALSE where uncorrected"), {
    skip_if_not(workflowOwnerCapable())
    fixture <- kit242Owner(kit242Object(kit242RuvIncompleteValues(12L), class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    try(suppressWarnings(suppressMessages(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config))),
      silent = TRUE)
    # r1.5 (K242-RC4-N06): under Q6.2 (c) the step completes.
    expect_identical(owner$current_state, paste0(prefix, "_ruv_corrected"))
    if (!identical(owner$current_state, paste0(prefix, "_ruv_corrected"))) return(invisible())
    before <- owner$getState(paste0(prefix, "_normalized"))
    after <- owner$getState(paste0(prefix, "_ruv_corrected"))
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_ruv_corrected"))
    named <- names(outcomes$k_per_assay)
    expect_setequal(named, c("A1", "A2"))
    for (assay in named) {
      corrected <- kit242Matrix(after, binding, assay)
      expect_false(is.null(corrected), info = paste(assay, "named by the record, present in the RUV state"))
      applied <- outcomes$ruv_applied_per_assay[[assay]]
      expect_true(is.logical(applied) && length(applied) == 1L, info = paste(assay, "per-assay outcome present"))
      if (!is.null(corrected) && kit242Equal(kit242Matrix(before, binding, assay), corrected)) {
        expect_false(isTRUE(applied), info = paste(assay, "unchanged, so not recorded as RUV-corrected"))
      }
    }
  })

  # ---- E11c: older records (OD-4 (c)): the RUV marker and its earlier rule --------------------------------------------
  test_that(paste("KIT-242 E11c:", class, "the RUV marker: current 'by_cause', earlier 'drop' reproduces the earlier result, and replays"), {
    marker <- workflowReplayRuleMarkers$ruv_fallback
    expect_identical(marker$option, "MultiScholaR.ruv_fallback")
    expect_identical(marker$current, "by_cause")
    expect_identical(marker$earlier, "drop")
    expect_identical(workflowOperationEntry(paste0(prefix, "_ruv"))$declaration$fixed$ruv_fallback, "by_cause")
    expect_identical(omicsRuvOperation(binding, "declared")$fixed$ruv_fallback, "by_cause")
    # (ii) under the earlier token E14's call gives the earlier result: A2 dropped, no classed warning.
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(14L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    earlier <- withr::with_options(kit242RuvRule("drop"), kit242Conditions(ruvIII_C_Varying(object, ruv_grouping_variable = "group",
      ruv_number_k = list(A1 = 1L, A2 = 1L), ctrl = list(A1 = sprintf("A1_F%03d", 1:8), A2 = c("absent_1", "absent_2"))),
      "omics_ruv_not_applied"))
    expect_length(earlier$kept, 0L)
    expect_identical(names(methods::slot(earlier$value, binding$data_slot)), "A1")
    # (iii) a record made under each rule replays with no mismatch (the earlier one without the marker).
    skip_if_not(workflowOwnerCapable())
    harness <- new.env(parent = globalenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    for (rule in c("drop", "by_cause")) {
      root <- withr::local_tempdir()
      fixture <- kit172Quiet(kit172Owner(root, kit242Object(kit242RuvIncompleteValues(24L), class), class))
      owner <- fixture$owner
      config <- fixture$workflow$config_list
      withr::with_options(kit242RuvRule(if (identical(rule, "drop")) "drop" else NULL), {
        kit242RuvChain(binding, owner, config)
        kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config))
      })
      replayed <- kit242ReplayRecord(fixture, root, harness, api)
      step <- Filter(function(s) identical(s$state, paste0(prefix, "_ruv_corrected")), replayed$record$steps)[[1L]]
      has_marker <- any(vapply(step$parameters, function(p) identical(p$name, "ruv_fallback"), logical(1)))
      expect_identical(has_marker, identical(rule, "by_cause"), info = rule)
      expect_identical(replayed$replay$states, replayed$history, info = rule)
      expect_identical(as.character(replayed$replay$mismatches), character(), info = rule)
    }
  })

  # ---- G4: legacy Automatic where every assay corrects: identical under the earlier and the current rule ------------
  # ---- G8: legacy Manual, no sparse feature: the same ---------------------------------------------------------------
  for (mode in c("automatic", "manual")) {
    test_that(paste("KIT-242", if (identical(mode, "automatic")) "G4:" else "G8:", class, mode,
      "legacy RUV where every assay corrects is identical under the earlier and the current rule"), {
      skip_if_not(workflowOwnerCapable())
      object <- kit242Object(kit242Draw(if (identical(mode, "automatic")) 44L else 48L, n_is = c(A1 = 0L, A2 = 0L)), class)
      parameters <- if (identical(mode, "automatic")) kit172LegacyRuv("automatic") else kit242ManualRuv()
      journey <- function(rule) {
        fixture <- kit242Owner(object, class, env = parent.frame())
        owner <- fixture$owner
        config <- fixture$workflow$config_list
        withr::with_options(kit242RuvRule(rule), {
          kit242RuvChain(binding, owner, config)
          kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, parameters, config))
        })
        state <- owner$getState(paste0(prefix, "_ruv_corrected"))
        list(tables = methods::slot(state, binding$data_slot), design = state@design_matrix,
          outcomes = kit242Outcomes(owner, paste0(prefix, "_ruv_corrected")))
      }
      current <- journey(NULL)
      earlier <- journey("drop")
      # RUV corrected every assay (the guard's precondition).
      expect_identical(unlist(current$outcomes$ruv_status_per_assay), c(A1 = "corrected", A2 = "corrected"))
      expect_identical(current$tables, earlier$tables)
      expect_identical(current$design, earlier$design)
      added <- c("ruv_applied_per_assay", "ruv_not_applied_reason_per_assay")
      expect_identical(current$outcomes[setdiff(names(current$outcomes), added)], earlier$outcomes)
    })
  }

  # ---- G5: the declared design on an object with no empty sample: unchanged -------------------------------------------
  test_that(paste("KIT-242 G5:", class, "the declared design with no empty sample is identical under the earlier and the current rule"), {
    skip_if_not_installed("ruv")
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(45L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    current <- kit172Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user"))
    earlier <- withr::with_options(kit242RuvRule("drop"), kit172Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user")))
    expect_identical(current$factors, earlier$factors)
    expect_identical(methods::slot(current$corrected, binding$data_slot), methods::slot(earlier$corrected, binding$data_slot))
  })

  # ---- E25m (the coordinator's ruling 2 after P0): the Manual search finds no control set for an assay it cannot read --
  test_that(paste("KIT-242 E25m:", class, "a Manual search with no control set for A2 (duplicated ids) stops, naming the cause"), {
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(51L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    tables <- methods::slot(object, binding$data_slot)
    tables$A2 <- rbind(tables$A2, tables$A2[c(3L, 7L), , drop = FALSE])
    methods::slot(object, binding$data_slot, check = FALSE) <- tables
    parameters <- kit242ManualRuv()
    results <- kit172Quiet(binding$optimise_ruv(object, "manual", parameters[setdiff(names(parameters), "ruv_mode")]))
    # On the code before Part B the function does not exist: the row fails by value there, not by an error.
    searchOutcomes <- get0("omicsRuvSearchOutcomes", envir = asNamespace("MultiScholaR"), mode = "function")
    expect_true(is.function(searchOutcomes), info = "omicsRuvSearchOutcomes is absent: the search's outcome is not decided by cause")
    if (!is.function(searchOutcomes)) return(invisible())
    outcome <- searchOutcomes(results, "manual", object)
    expect_identical(outcome$A2$action, "stop")
    expect_identical(outcome$A2$reason, "search_error")
    expect_match(outcome$A2$text, "found no control features for assay A2", fixed = TRUE)
    expect_match(outcome$A2$text, "2 duplicated feature identifiers (A2_F003, A2_F007)", fixed = TRUE)
    expect_match(outcome$A2$text, "Run 'Resolve Duplicates' in the Quality Control tab", fixed = TRUE)
    expect_identical(outcome$A1$action, "correct")
  })

  # ---- G9: a feature with values that RUV-III-C returns all non-finite is still removed (`:713-715`), as today --------
  test_that(paste("KIT-242 G9:", class, "a feature with values returned all non-finite by RUV-III-C is removed, as today"), {
    object <- kit242Tibbles(kit242Object(lapply(kit242Draw(49L, n_is = c(A1 = 0L, A2 = 0L)), log2), class), binding)
    real <- RUVIIIC::RUVIII_C_Varying
    local_mocked_bindings(RUVIII_C_Varying = function(k, Y, M, toCorrect, potentialControls) {
      out <- real(k = k, Y = Y, M = M, toCorrect = toCorrect, potentialControls = potentialControls)
      if ("A2_F040" %in% colnames(out)) out[, "A2_F040"] <- NA_real_
      out
    }, .package = "MultiScholaR")
    out <- suppressWarnings(suppressMessages(ruvIII_C_Varying(object, ruv_grouping_variable = "group",
      ruv_number_k = list(A1 = 1L, A2 = 1L), ctrl = list(A1 = sprintf("A1_F%03d", 1:8), A2 = sprintf("A2_F%03d", 1:8)))))
    a2 <- kit242Matrix(out, binding, "A2")
    expect_false(is.null(a2))
    if (!is.null(a2)) expect_false("A2_F040" %in% rownames(a2))
  })
}

# ==== VAL-373 (Q-C): the one bounded correction after the implementation audit ========================================
# E17 through the DIA normalisation declaration (audit I04): the sample with no value that leaves the design is recorded
# and told. E27 (audit I02, VAL-373 (Q-A) option (b)): an assay in which a group is left with fewer than two samples
# holding values keeps its values without RUV, recorded and shown as `no_usable_controls` with the groups named; the
# other assays are corrected, on the Automatic and the Manual search; KIT-275's VE shape (group B with one sample in
# every assay) keeps every assay without RUV.

test_that("KIT-242 E17 (VAL-373 Q-C): the DIA normalisation step records and tells a sample removed for having no value", {
  y <- kit242ProteinValues(17L)
  y[, "B_6"] <- NA_real_
  object <- kit242ProteinObject(y)
  store <- list(checkpoint = function(x, name, from = NULL) x)
  parameters <- list(norm_method = "cyclicloess", ruv_mode = "skip", auto_percentage_min = 1, auto_percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 10,
    ruv_k = 1L, grouping = "group")
  run <- kit242Conditions(diaNormalizationOperation()$kernel(object, parameters, store), "omics_normalisation_sample_removed")
  expect_length(run$kept, 0L) # captured inside the step, not escaping as a warning
  expect_identical(run$value$samples_removed, "B_6")
  expect_false("B_6" %in% as.character(run$value$normalized@design_matrix$Run))
  # The record: the commit's outcomes name the sample.
  local_mocked_bindings(readDIAStateMetadata = function(manager, reference = NULL, ...) list(data = list(products = list())))
  commit <- diaNormalizationCommit(list(current_state = "initial"), list(norm_method = "cyclicloess", ruv_mode = "skip"),
    list(value = run$value, reference = NULL))
  expect_identical(commit$config$outcomes$samples_removed, "B_6")
  # The tab: a warning notification names the sample.
  notes <- character()
  finishProtNormNormalization(TRUE, list(ruv_mode = "skip", norm_method = "cyclicloess"), list(),
    list(samples_removed = run$value$samples_removed), buildProtNormCompletionNotification,
    function(text, type = "default", ...) notes <<- c(notes, paste(type, text)), function(...) NULL)
  expect_true(any(grepl("^warning .*B_6.*removed from the design", notes)), info = paste(notes, collapse = " | "))
})

for (class in c("lipid", "metab")) {
  binding <- kit172Binding(class)
  prefix <- binding$prefix

  for (mode in c("automatic", "manual")) {
    test_that(paste("KIT-242 E27 (VAL-373 Q-A):", class, mode, "A2 with group B left with one sample keeps its values without RUV, A1 corrected"), {
      skip_if_not(workflowOwnerCapable())
      v <- kit242Draw(27L, n_is = c(A1 = 0L, A2 = 0L))
      v$A2[, c("B_2", "B_3", "B_4", "B_5", "B_6")] <- NA_real_
      fixture <- kit242Owner(kit242Object(v, class), class)
      owner <- fixture$owner
      config <- fixture$workflow$config_list
      kit242RuvChain(binding, owner, config)
      normalised <- owner$getState(paste0(prefix, "_normalized"))
      parameters <- if (identical(mode, "automatic")) kit172LegacyRuv("automatic") else kit242ManualRuv()
      run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, parameters, config)), error = function(e) e)
      expect_false(inherits(run, "error"), info = if (inherits(run, "error")) conditionMessage(run) else "")
      if (inherits(run, "error")) return(invisible())
      after <- owner$getState(paste0(prefix, "_ruv_corrected"))
      expect_equal(kit242Matrix(after, binding, "A2"), kit242Matrix(normalised, binding, "A2"), tolerance = 0)
      expect_false(kit242Equal(kit242Matrix(normalised, binding, "A1"), kit242Matrix(after, binding, "A1")))
      outcomes <- kit242Outcomes(owner, paste0(prefix, "_ruv_corrected"))
      expect_identical(outcomes$ruv_applied_per_assay, list(A1 = TRUE, A2 = FALSE))
      expect_identical(outcomes$ruv_not_applied_reason_per_assay$A2, "no_usable_controls")
      expect_identical(outcomes$ruv_not_applied_groups_per_assay, list(A2 = "B"))
      expect_match(run$value$ruv_not_applied_text$A2, "the groups B have fewer than two samples with values", fixed = TRUE)
    })
  }

  test_that(paste("KIT-242 E27 (VAL-373 Q-A):", class, "KIT-275's VE shape: group B with one sample in every assay keeps every assay without RUV"), {
    skip_if_not(workflowOwnerCapable())
    v <- kit242Draw(127L + 100L, n_is = c(A1 = 0L, A2 = 0L))
    for (a in c("A1", "A2")) v[[a]][, c("B_2", "B_3", "B_4", "B_5", "B_6")] <- NA_real_
    fixture <- kit242Owner(kit242Object(v, class), class)
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    kit242RuvChain(binding, owner, config)
    normalised <- owner$getState(paste0(prefix, "_normalized"))
    run <- tryCatch(kit172Quiet(runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)), error = function(e) e)
    expect_false(inherits(run, "error"), info = if (inherits(run, "error")) conditionMessage(run) else "")
    if (inherits(run, "error")) return(invisible())
    after <- owner$getState(paste0(prefix, "_ruv_corrected"))
    for (assay in c("A1", "A2")) expect_equal(kit242Matrix(after, binding, assay), kit242Matrix(normalised, binding, assay), tolerance = 0)
    outcomes <- kit242Outcomes(owner, paste0(prefix, "_ruv_corrected"))
    expect_false(isTRUE(outcomes$ruv_applied))
    expect_identical(outcomes$ruv_not_applied_reason_per_assay, list(A1 = "no_usable_controls", A2 = "no_usable_controls"))
    expect_identical(outcomes$ruv_not_applied_groups_per_assay, list(A1 = "B", A2 = "B"))
  })
}
