# KIT-076: lipid and metab DA as one declaration per contrast on a saved per-assay fit: the kernel's halves give the
# all-contrast kernel's numbers, the manifest is keyed by contrast, the module's list is rebuilt from it, a failed
# contrast fails alone, and a v1 record still resolves its operation.
library(testthat)

omicsDaContrastsFixture <- function(omic) {
  samples <- paste0("S", seq_len(12L))
  design <- data.frame(Run = samples, group = rep(c("A", "B", "C"), each = 4L), batch = rep(c("b1", "b2"), 6L), stringsAsFactors = FALSE)
  assay <- function(seed, id, annotation) {
    withr::with_seed(seed, {
      values <- matrix(stats::rnorm(60L * 12L, 10, 1), nrow = 60L, dimnames = list(NULL, samples))
      values[seq_len(10L), design$group == "B"] <- values[seq_len(10L), design$group == "B"] + 2
    })
    table <- data.frame(paste0("F", seed, "_", seq_len(60L)), paste("Name", seq_len(60L)), stringsAsFactors = FALSE)
    names(table) <- c(id, annotation)
    cbind(table, as.data.frame(values, check.names = FALSE))
  }
  if (identical(omic, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = list(pos = assay(1L, "lipid_id", "lipid_name"), neg = assay(2L, "lipid_id", "lipid_name")),
      design_matrix = design, sample_id = "Run", group_id = "group", lipid_id_column = "lipid_id", annotation_id_column = "lipid_name")
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = list(pos = assay(1L, "metabolite_id", "metabolite_name"), neg = assay(2L, "metabolite_id", "metabolite_name")),
      design_matrix = design, sample_id = "Run", group_id = "group", metabolite_id_column = "metabolite_id", annotation_id_column = "metabolite_name")
  }
}

omicsDaContrastsTable <- function(contrasts = c("groupB-groupA", "groupC-groupA", "groupC-groupB")) {
  data.frame(contrasts = contrasts, friendly_names = gsub("group", "", sub("-", "_vs_", contrasts)), stringsAsFactors = FALSE)
}

omicsDaContrastsRun <- function(omic, object = omicsDaContrastsFixture(omic), contrasts = omicsDaContrastsTable(), lfc = 0, formula = "~ 0 + group", env = parent.frame()) {
  workflow <- omicsTestOwner(object, "da_input", env = env)
  manager <- workflow$state_manager
  wrapper <- if (identical(omic, "lipid")) runLipidReferenceDaAnalysis else runMetabReferenceDaAnalysis
  run <- suppressWarnings(suppressMessages(wrapper(manager, getWorkflowStateReference(manager),
    list(formula_string = formula, da_q_val_thresh = 0.05, treat_lfc_cutoff = lfc), contrasts)))
  list(manager = manager, run = run)
}

omicsDaKernel <- function(omic) if (identical(omic, "lipid")) runLipidsDA else runMetabolitesDA

for (omic in c("lipid", "metab")) {
  long_table <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"

  test_that(paste("KIT-076:", omic, "per-contrast DA on saved fits equals the all-contrast kernel"), {
    for (lfc in c(0, 0.5)) {
      object <- omicsDaContrastsFixture(omic)
      direct <- suppressWarnings(suppressMessages(omicsDaKernel(omic)(theObject = object, contrasts_tbl = omicsDaContrastsTable(),
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = lfc)))
      result <- omicsDaContrastsRun(omic, object, lfc = lfc)$run
      expect_identical(result$results[[long_table]], direct[[long_table]])
      expect_identical(result$results$per_assay_results, direct$per_assay_results)
      expect_identical(result$results$significant_counts, direct$significant_counts)
      expect_identical(result$results$qvalue_warnings, direct$qvalue_warnings)
      for (assay in c("pos", "neg")) {
        expect_identical(result$results$contrasts_results[[assay]]$results, direct$contrasts_results[[assay]]$results)
        # The coefficient view Glimma reads: the same names and row order, the same values.
        expect_identical(colnames(result$results$contrasts_results[[assay]]$fit.eb$coefficients), omicsDaContrastsTable()$contrasts)
        expect_identical(rownames(result$results$contrasts_results[[assay]]$fit.eb$coefficients), rownames(direct$contrasts_results[[assay]]$fit.eb$coefficients))
        expect_equal(unname(result$results$contrasts_results[[assay]]$fit.eb$coefficients), unname(direct$contrasts_results[[assay]]$fit.eb$coefficients), tolerance = 0)
      }
      expect_length(result$results$failed_contrasts, 0L)
    }
  })

  test_that(paste("KIT-076:", omic, "manifest v2 is keyed by contrast and shares one saved fit per assay"), {
    out <- omicsDaContrastsRun(omic)
    manifest <- out$run$manifest
    context <- out$manager$.context
    expect_identical(manifest$schema, paste0(omic, ".da_results.v2"))
    expect_identical(names(manifest$results), omicsDaContrastsTable()$contrasts)
    expect_identical(manifest$assays, c("pos", "neg"))
    for (assay in c("pos", "neg")) {
      fits <- vapply(manifest$results, function(entry) entry$assays[[assay]]$fit$entity_id, "")
      expect_length(unique(fits), 1L)
    }
    found <- findSavedWorkflowEntities(context$store, list(kind = "da_fit"), max_matches = 10L)
    expect_length(found$matches, 2L)
    # Each contrast's saved fit.eb holds that contrast only.
    value <- withWorkflowStateScope(out$manager, function() hydrateSavedWorkflowEntity(context$store, manifest$results[[1L]]$products$da_contrast)$value)
    expect_identical(colnames(value$assays$pos$fit.eb$coefficients), "groupB-groupA")
    expect_s4_class(value$assays$pos$fit.eb, "MArrayLM")
    # The head holds the manifest, and the record labels each entry by contrast and assay.
    head_manifest <- withWorkflowStateScope(out$manager, function() resolveWorkflowStateValue(context, context$snapshot$da_results))
    expect_identical(head_manifest$schema, manifest$schema)
    record <- buildWorkflowAnalysisRecord(out$manager)
    da <- record$analyses$differential_abundance
    expect_identical(da$operation, paste0(omic, "_da_contrast"))
    expect_identical(vapply(da$results, function(r) r$contrast, ""), rep(omicsDaContrastsTable()$contrasts, each = 2L))
    expect_identical(vapply(da$results, function(r) r$assay, ""), rep(c("pos", "neg"), 3L))
    expect_identical(vapply(da$results, function(r) r$rows, 1L), rep(60L, 6L))
  })

  test_that(paste("KIT-076:", omic, "a failing contrast fails alone and is recorded"), {
    object <- omicsDaContrastsFixture(omic)
    binding <- omicsDeclarationBindings()[[if (identical(omic, "lipid")) "lipidomics" else "metabolomics"]]
    # A contrast the preflight cannot see: it fails in the contrast half only.
    halves <- paste0("apply", if (identical(omic, "lipid")) "Lipid" else "Metab", "DaContrasts")
    original <- get(halves, envir = asNamespace("MultiScholaR"))
    failing <- function(fit, contrast_strings, ...) {
      if (identical(contrast_strings, "groupC-groupA")) stop("limma could not fit this contrast")
      original(fit, contrast_strings, ...)
    }
    do.call(testthat::local_mocked_bindings, c(stats::setNames(list(failing), halves), list(.package = "MultiScholaR")))
    direct <- suppressWarnings(suppressMessages(omicsDaKernel(omic)(theObject = object,
      contrasts_tbl = omicsDaContrastsTable(c("groupB-groupA", "groupC-groupB")), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)))
    out <- omicsDaContrastsRun(omic, object)
    results <- out$run$results
    expect_identical(names(results$failed_contrasts), "groupC-groupA")
    expect_identical(results$failed_contrasts[["groupC-groupA"]], list(pos = "limma could not fit this contrast", neg = "limma could not fit this contrast"))
    expect_identical(unique(results[[long_table]]$comparison), c("groupB-groupA", "groupC-groupB"))
    # The other contrasts carry exactly what a run without the failing contrast gives.
    expect_identical(results$per_assay_results, direct$per_assay_results)
    manifest <- out$run$manifest
    expect_length(manifest$results[["groupC-groupA"]]$assays, 0L)
    record <- buildWorkflowAnalysisRecord(out$manager)
    failures <- record$analyses$differential_abundance$failures
    expect_identical(vapply(failures, function(f) f$contrast, ""), rep("groupC-groupA", 2L))
    expect_identical(vapply(failures, function(f) f$error, ""), rep("limma could not fit this contrast", 2L))
    expect_identical(validateWorkflowAnalysisRecord(record)$valid, TRUE)
  })

  test_that(paste("KIT-076:", omic, "an undefined contrast term fails alone; when every contrast fails the run stops and the head is unchanged"), {
    object <- omicsDaContrastsFixture(omic)
    # The lipid kernel's preflight refuses it before the contrast half; the metab kernel has none, so the term check does.
    out <- omicsDaContrastsRun(omic, object, contrasts = omicsDaContrastsTable(c("groupB-groupA", "groupZ-groupA")))
    expect_identical(names(out$run$results$failed_contrasts), "groupZ-groupA")
    expect_match(out$run$results$failed_contrasts[["groupZ-groupA"]]$pos, "groupZ")
    expect_identical(unique(out$run$results[[long_table]]$comparison), "groupB-groupA")
    expect_identical(buildWorkflowAnalysisRecord(out$manager)$analyses$differential_abundance$failures[[1L]]$contrast, "groupZ-groupA")
    workflow <- omicsTestOwner(omicsDaContrastsFixture(omic), "da_input")
    manager <- workflow$state_manager
    wrapper <- if (identical(omic, "lipid")) runLipidReferenceDaAnalysis else runMetabReferenceDaAnalysis
    expect_error(suppressWarnings(suppressMessages(wrapper(manager, getWorkflowStateReference(manager),
      list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), omicsDaContrastsTable("groupZ-groupA")))), "groupZ")
    expect_null(manager$.context$snapshot$da_results)
  })

  test_that(paste("KIT-076:", omic, "a v1 manifest still reads in the record and a v1 record's operation still replays per contrast"), {
    out <- omicsDaContrastsRun(omic, contrasts = omicsDaContrastsTable("groupB-groupA"))
    expect_identical(workflowOutputDeclaration(paste0(omic, ".da_results.v1"))$operation, paste0(omic, "_da_output"))
    expect_false(is.null(workflowOperationEntry(paste0(omic, "_da_output"))$declaration))
    context <- out$manager$.context
    # A v1 manifest on the head: the record keeps its reading, assay names under `contrast`.
    v1 <- list(schema = paste0(omic, ".da_results.v1"), source = out$run$manifest$source, contrasts = omicsDaContrastsTable("groupB-groupA"),
      results = list(pos = list(assay = "pos", tables = list(), rows = 60L)), parameters = out$run$manifest$parameters, model = list(model = "limma"))
    withWorkflowStateScope(out$manager, function() {
      reference <- saveWorkflowStateValue(context, v1, "da-result-references", codec = "rds")
      snapshot <- context$snapshot
      snapshot$da_results <- reference
      publishWorkflowStateHead(context, snapshot)
    })
    da <- buildWorkflowAnalysisRecord(out$manager)$analyses$differential_abundance
    expect_identical(da$operation, paste0(omic, "_da_output"))
    expect_identical(da$results[[1L]]$contrast, "pos")
    expect_null(da$failures)
  })
}
