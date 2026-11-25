# KIT-249 gate G0 (ticket r1.2, section 8): observe, at the base commit, what each route of runTestsContrasts puts in
# fdr_qvalue when qvalue() fails. No package code: branch kit-249-g0 carries only this script and is never merged.
# Nothing here judges a route: each row records what the route returned; the comparison with the ticket's expectations
# is made by reading the tables.
#
# Usage (run-tests scripts lane, from the package root, after the package commit is installed):
#   Rscript tools/kit249/g0.R <out dir> [<c shift>]
# <c shift>: added to every c_j of the null proteins (0 by default). Section 8's instrument rule: if qvalue() does not
# error on B_vs_A, G0 is re-sent with the shift raised by 0.10, at most three times; no expectation changes.
#
# The constructed comparison (section 8): groups A, B, C of 4 samples (S1..S12, `~ 0 + group`), proteins P01..P30, no
# missing value, no random numbers. mu_j = 20 + 0.25 (j - 1); sigma_j = 0.20 (1 + 0.05 ((j - 1) mod 5));
# e_A = (-1.5, -0.5, 0.5, 1.5), e_B = (0.5, 1.5, -1.5, -0.5), e_C = (1.5, -1.5, -0.5, 0.5);
# A = mu + sigma e_A; B = mu + beta + sigma e_B with beta_j = 2 (P01..P06) and c_j sigma_j, c_j = 0.30 + 0.05 (j - 7)
# (P07..P30, true nulls written as a fixed realisation); C = mu + gamma + sigma e_C with gamma_j = 2 (P01..P06), else 0.
#
# Routes: G0-a runTestsContrasts topTable (treat_lfc_cutoff NA, trend and robust TRUE); G0-b the same with
# treat_lfc_cutoff 0 (topTreat, the application's default); G0-c S4 differentialAbundanceAnalysis on a
# ProteinQuantitativeData (workflow type LFQ), one comparison per call as the application's kernel; G0-d the
# application's kernel in process (runDIAReferenceDaAnalysis on a DIA reference owner, KIT-200 T6's owner pattern), its
# manifest and its analysis record; G0-e protein_deAnalysisWrapperFunction (daAnalysisWrapperFunction is not run:
# KIT-218 has not merged at the base); G0-f the lipidomics and metabolomics exported S4 differentialAbundanceAnalysis
# (the list method) with the matrix as one assay. Controls: G0-c1 adjustProtDaResultTable on G0-b's B_vs_A table;
# G0-c2 applyLipidDaContrasts on the matrix.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
out <- args[[1L]]
shift <- if (length(args) == 2L) as.numeric(args[[2L]]) else 0
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)

routes <- quote({

PLANTED <- sprintf("P%02d", 1:6)
NOTIFICATION <- "Some contrasts used the existing BH fallback after q-value estimation failed."
FULL <- c(B_vs_A = "B_vs_A=groupB-groupA", C_vs_A = "C_vs_A=groupC-groupA")
oneLine <- function(x) gsub("[\t\r\n]+", " ", paste(as.character(x), collapse = " "))

g0Matrix <- function(shift = 0) {
  j <- 1:30
  mu <- 20 + 0.25 * (j - 1)
  sigma <- 0.20 * (1 + 0.05 * ((j - 1) %% 5))
  beta <- ifelse(j <= 6, 2, (0.30 + 0.05 * (j - 7) + shift) * sigma)
  gamma <- ifelse(j <= 6, 2, 0)
  e <- list(A = c(-1.5, -0.5, 0.5, 1.5), B = c(0.5, 1.5, -1.5, -0.5), C = c(1.5, -1.5, -0.5, 0.5))
  x <- cbind(mu + outer(sigma, e$A), mu + beta + outer(sigma, e$B), mu + gamma + outer(sigma, e$C))
  dimnames(x) <- list(sprintf("P%02d", j), paste0("S", 1:12))
  x
}

g0Design <- function() {
  runs <- paste0("S", 1:12)
  data.frame(Run = runs, group = rep(c("A", "B", "C"), each = 4L), replicates = runs, row.names = runs,
    stringsAsFactors = FALSE)
}

g0Object <- function(x, args = list(globalParameters = list(workflow_type = "LFQ"))) {
  design <- g0Design()
  rownames(design) <- NULL
  ProteinQuantitativeData(protein_quant_table = data.frame(Protein.Ids = rownames(x), x, row.names = NULL,
    check.names = FALSE, stringsAsFactors = FALSE), design_matrix = design, sample_id = "Run", group_id = "group",
    technical_replicate_id = "replicates", protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(Protein.Ids = rownames(x), gene_names = paste0("GENE_", rownames(x)),
      stringsAsFactors = FALSE), args = args)
}

quiet <- function(expr) {
  value <- NULL
  utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
  value
}

# The kernel's own value, captured from inside a route (the wrappers and the lipid/metab S4 method keep no
# qvalue_warnings of their own).
captureRunTests <- function(code_fn, extra = list()) {
  seen <- new.env(parent = emptyenv())
  run_tests <- get("runTestsContrasts", envir = asNamespace("MultiScholaR"))
  bindings <- c(list(runTestsContrasts = function(...) {
    value <- run_tests(...)
    seen$value <- value
    value
  }), extra)
  runner <- function() {
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "MultiScholaR", .env = environment())))
    code_fn()
  }
  result <- tryCatch(runner(), error = function(e) e)
  list(result = result, kernel = seen$value)
}

# One table's row of g0-routes.tsv, and its proteins' rows.
ROWS <- list()
PROTEINS <- list()
record <- function(route, comparison, table, warnings = NULL, extra = list(), id = NULL, bh_column = "fdr_value_bh_adjustment") {
  row <- list(route = route, comparison = comparison, error = "")
  if (inherits(table, "error") || is.null(table)) {
    row$error <- if (is.null(table)) "no table" else oneLine(conditionMessage(table))
    ROWS[[length(ROWS) + 1L]] <<- as.data.frame(c(row, extra), stringsAsFactors = FALSE)
    return(invisible())
  }
  table <- as.data.frame(table)
  ids <- if (!is.null(id) && id %in% names(table)) as.character(table[[id]]) else rownames(table)
  p <- table$raw_pvalue
  q <- table$fdr_qvalue
  bh <- table[[bh_column]]
  valid <- !is.na(p) & is.finite(p)
  planted <- ids %in% PLANTED
  tried <- function(expr) tryCatch(list(ok = TRUE, value = expr), error = function(e) list(ok = FALSE, value = conditionMessage(e)))
  storey <- tried(qvalue::qvalue(p[valid]))
  boot <- tried(qvalue::qvalue(p[valid], pi0.method = "bootstrap"))
  half <- tried(qvalue::qvalue(p[valid], lambda = 0.5)$pi0)
  row <- c(row, list(rows = nrow(table), valid_p = sum(valid), max_p = sprintf("%.6f", max(p[valid])),
    max_p_planted = sprintf("%.3g", max(p[valid & planted])), min_p_null = sprintf("%.6f", min(p[valid & !planted])),
    qvalue_errors = !storey$ok, qvalue_message = if (storey$ok) "" else oneLine(storey$value),
    bootstrap_errors = !boot$ok, bootstrap_message = if (boot$ok) "" else oneLine(boot$value),
    pi0_lambda_0.5 = if (half$ok) sprintf("%.6f", half$value) else paste("error:", oneLine(half$value)),
    nonNA_fdr_qvalue = sum(!is.na(q)), identical_q_bh = identical(q, bh),
    q_equals_p_adjust_BH = isTRUE(all.equal(q[valid], stats::p.adjust(p[valid], "BH"), tolerance = 0)),
    q_equals_storey = storey$ok && isTRUE(all.equal(q[valid], storey$value$qvalues, tolerance = 0)),
    qvalue_warnings = if (is.null(warnings)) "NA (not returned by this route)" else
      if (length(warnings)) oneLine(paste(names(warnings), collapse = ";")) else "<none>",
    sig_planted_of_6 = sum(!is.na(q) & q < 0.05 & planted), sig_null_of_24 = sum(!is.na(q) & q < 0.05 & !planted),
    bh_sig_planted_of_6 = sum(!is.na(bh) & bh < 0.05 & planted), bh_sig_null_of_24 = sum(!is.na(bh) & bh < 0.05 & !planted)))
  ROWS[[length(ROWS) + 1L]] <<- as.data.frame(c(row, extra), stringsAsFactors = FALSE)
  PROTEINS[[length(PROTEINS) + 1L]] <<- data.frame(route = route, comparison = comparison, protein = ids,
    logFC = sprintf("%.6f", table$logFC), raw_pvalue = sprintf("%.6g", p), fdr_qvalue = sprintf("%.6g", q),
    bh = sprintf("%.6g", bh), stringsAsFactors = FALSE)
}

warningsFor <- function(warnings, full) warnings[intersect(names(warnings), full)]

runAll <- function(shift) {
  x <- g0Matrix(shift)
  design <- g0Design()

  # G0-a and G0-b: runTestsContrasts directly
  for (route in c("G0-a", "G0-b")) {
    lfc <- if (identical(route, "G0-a")) NA else 0
    value <- tryCatch(quiet(runTestsContrasts(x, unname(FULL), design, "~ 0 + group", treat_lfc_cutoff = lfc,
      eBayes_trend = TRUE, eBayes_robust = TRUE)), error = function(e) e)
    for (cmp in names(FULL)) {
      if (inherits(value, "error")) { record(route, cmp, value); next }
      record(route, cmp, value$results[[FULL[[cmp]]]], warningsFor(value$qvalue_warnings, FULL[[cmp]]))
    }
    if (identical(route, "G0-b") && !inherits(value, "error")) {
      # G0-c1 (control): adjustProtDaResultTable on G0-b's B_vs_A table, from its raw p-values
      tbl <- value$results[[FULL[["B_vs_A"]]]]
      tbl$P.Value <- tbl$raw_pvalue
      tbl$raw_pvalue <- NULL
      tbl$fdr_qvalue <- NULL
      tbl$fdr_value_bh_adjustment <- NULL
      adjusted <- tryCatch(adjustProtDaResultTable(tbl, "raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment",
        messageFn = function(...) invisible(NULL)), error = function(e) e)
      if (inherits(adjusted, "error")) record("G0-c1", "B_vs_A", adjusted) else
        record("G0-c1", "B_vs_A", adjusted$table, extra = list(note = paste("qvalue_failed", adjusted$qvalue_failed)))
    }
  }

  # G0-c: the S4 method, one comparison per call (as the application's kernel)
  object <- g0Object(x)
  for (cmp in names(FULL)) {
    value <- tryCatch(quiet(differentialAbundanceAnalysis(object,
      contrasts_tbl = data.frame(contrasts = FULL[[cmp]], stringsAsFactors = FALSE), formula_string = "~ 0 + group",
      da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
      args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")), error = function(e) e)
    if (inherits(value, "error")) { record("G0-c", cmp, value); next }
    counts <- tryCatch(oneLine(paste(value$num_sig_da_molecules$status, value$num_sig_da_molecules$counts, sep = ": ")),
      error = function(e) paste("error:", conditionMessage(e)))
    record("G0-c", cmp, value$contrasts_results_table, value$qvalue_warnings %||% list(), id = "Protein.Ids",
      extra = list(note = paste("count table:", counts)))
  }

  # G0-d: the application's kernel in process on a DIA reference owner, its manifest and its analysis record
  capable <- workflowOwnerCapable()
  if (!capable) {
    for (cmp in names(FULL)) record("G0-d", cmp, simpleError("workflowOwnerCapable() is FALSE"))
  } else {
    for (publish in c(TRUE, FALSE)) {
      root <- normalizePath(tempfile("kit249-g0-owner-"), mustWork = FALSE)
      dir.create(root)
      paths <- list(source_dir = file.path(root, "source"), protein_qc_dir = file.path(root, "protein_qc"),
        time_dir = file.path(root, "time"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
      for (path in paths) dir.create(path, recursive = TRUE, showWarnings = FALSE)
      manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit249-dia", "DIA"), WorkflowState$new(audit_enabled = FALSE))
      manager$.context$reference_states <- TRUE
      manager$setWorkflowType("DIA")
      manager$.context$experiment_paths <- paths
      setDIAStepProvenance(manager, "design_builder")
      # KIT-200 T6's owner object (args empty), held as the correlation filter's state (T6 (d)'s reference owner)
      manager$saveState("correlation_filtered", g0Object(x, list()), list(), "fixture")
      workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list(),
        tab_status = list()))
      da <- new.env(parent = emptyenv())
      contrasts <- data.frame(contrasts = names(FULL), full_format = unname(FULL), friendly_names = names(FULL),
        stringsAsFactors = FALSE)
      workflow$contrasts_tbl <- contrasts
      da$contrasts_available <- contrasts$contrasts
      da$current_s4_object <- getDIAStateReference(manager)
      settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
      manifest <- tryCatch(quiet(runDIAReferenceDaAnalysis(settings, da, workflow, paths, contrasts, NULL,
        publish = publish)), error = function(e) e)
      label <- if (publish) "G0-d" else "G0-d (publish FALSE)"
      if (inherits(manifest, "error")) {
        for (cmp in names(FULL)) record(label, cmp, manifest)
        closeDIAStateAdapter(manager)
        next
      }
      flagged <- unlist(lapply(manifest$results, function(result) result$qvalue_warnings), recursive = FALSE)
      shown <- if (length(flagged)) NOTIFICATION else "<no warning notification>"
      fdr <- tryCatch({
        rec <- buildWorkflowAnalysisRecord(manager)
        entries <- rec$analyses$differential_abundance$results
        stats::setNames(vapply(entries, function(e) e$fdr_method %||% "NULL", ""),
          vapply(entries, function(e) e$contrast %||% "", ""))
      }, error = function(e) c(error = conditionMessage(e)))
      for (name in names(manifest$results)) {
        tables <- withDIAStateScope(manager, function() lapply(c(results = "contrasts_results_table",
          long = "da_proteins_long"), function(field) as.data.frame(materializeDIATable(manager$.context,
            manifest$results[[name]]$tables[[field]]))))
        long <- tables$long
        record(label, name, tables$results, manifest$results[[name]]$qvalue_warnings %||% list(), id = "Protein.Ids",
          extra = list(note = paste0("manifest qvalue_warnings: ", if (length(manifest$results[[name]]$qvalue_warnings))
            oneLine(names(manifest$results[[name]]$qvalue_warnings)) else "<none>",
            "; record fdr_method: ", if (name %in% names(fdr)) fdr[[name]] else oneLine(paste(names(fdr), fdr, sep = "=")),
            "; notification: ", shown, "; da_proteins_long non-NA fdr_qvalue: ",
            if ("fdr_qvalue" %in% names(long)) sum(!is.na(long$fdr_qvalue)) else "no column",
            " of ", nrow(long))))
      }
      closeDIAStateAdapter(manager)
      break
    }
  }

  # G0-e: protein_deAnalysisWrapperFunction (plots before the DA step replaced, as KIT-200 T9)
  wrapper <- captureRunTests(function() quiet(protein_deAnalysisWrapperFunction(object,
    contrasts_tbl = tibble::tibble(contrasts = unname(FULL)), formula_string = "~ 0 + group", de_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE, args_row_id = "Protein.Ids")),
    extra = list(plotRle = function(...) ggplot2::ggplot(), plotPca = function(...) ggplot2::ggplot(),
      plotNumOfValuesNoLog = function(...) ggplot2::ggplot()))
  wrapper_note <- paste("wrapper:", if (inherits(wrapper$result, "error")) paste("stopped:",
    oneLine(conditionMessage(wrapper$result))) else "completed")
  for (cmp in names(FULL)) {
    if (is.null(wrapper$kernel)) { record("G0-e", cmp, wrapper$result, extra = list(note = wrapper_note)); next }
    record("G0-e", cmp, wrapper$kernel$results[[FULL[[cmp]]]], warningsFor(wrapper$kernel$qvalue_warnings, FULL[[cmp]]),
      extra = list(note = wrapper_note))
  }

  # G0-f: lipidomics and metabolomics exported S4 differentialAbundanceAnalysis (list method), one assay
  assay_design <- g0Design()
  rownames(assay_design) <- NULL
  assay <- data.frame(Name = rownames(x), x, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(assay) <- NULL
  objects <- list(
    lipidomics = createLipidomicsAssayData(lipid_data = list(Assay1 = assay), design_matrix = assay_design,
      lipid_id_column = "Name", sample_id = "Run", group_id = "group"),
    metabolomics = createMetaboliteAssayData(metabolite_data = list(Assay1 = assay), design_matrix = assay_design,
      metabolite_id_column = "Name", sample_id = "Run", group_id = "group"))
  plain <- c(B_vs_A = "groupB-groupA", C_vs_A = "groupC-groupA")
  for (omic in names(objects)) {
    run <- captureRunTests(function() quiet(differentialAbundanceAnalysis(list(Assay1 = objects[[omic]]),
      contrasts_tbl = data.frame(contrasts = unname(plain), stringsAsFactors = FALSE), formula_string = "~ 0 + group",
      treat_lfc_cutoff = 0)))
    label <- paste("G0-f", omic)
    note <- paste("method:", if (inherits(run$result, "error")) paste("stopped:", oneLine(conditionMessage(run$result))) else "completed")
    for (cmp in names(plain)) {
      table <- if (!inherits(run$result, "error")) {
        tables <- run$result$Assay1@contrasts_results_table
        tables[[plain[[cmp]]]] %||% tables[[which(names(plain) == cmp)]]
      } else run$result
      record(label, cmp, table, if (is.null(run$kernel)) NULL else warningsFor(run$kernel$qvalue_warnings, plain[[cmp]]),
        id = "Name", extra = list(note = note))
    }
  }

  # G0-c2 (control): applyLipidDaContrasts on the matrix
  fit <- suppressWarnings(limma::lmFit(x, stats::model.matrix(~ 0 + group, design)))
  lipid <- tryCatch(quiet(applyLipidDaContrasts(fit, unname(plain), treat_lfc_cutoff = 0)), error = function(e) e)
  for (cmp in names(plain)) {
    if (inherits(lipid, "error")) { record("G0-c2", cmp, lipid); next }
    record("G0-c2", cmp, lipid$results[[plain[[cmp]]]], warningsFor(lipid$qvalue_warnings, plain[[cmp]]),
      bh_column = "fdr_value_bh")
  }

  list(routes = ROWS, proteins = PROTEINS)
}

})

env <- new.env(parent = asNamespace("MultiScholaR"))
eval(routes, env)
git <- function(...) suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
provenance <- c(sprintf("head\t%s", git("rev-parse", "HEAD")[[1L]]), sprintf("R\t%s", R.version.string),
  sprintf("qvalue\t%s", as.character(utils::packageVersion("qvalue"))),
  sprintf("limma\t%s", as.character(utils::packageVersion("limma"))),
  sprintf("MultiScholaR\t%s (installed, %s)", as.character(utils::packageVersion("MultiScholaR")),
    find.package("MultiScholaR")), sprintf("c_shift\t%s", format(shift)),
  sprintf("started\t%s", format(Sys.time(), tz = "UTC")))
result <- env$runAll(shift)
bind <- function(rows) {
  columns <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(row) {
    for (column in setdiff(columns, names(row))) row[[column]] <- ""
    row[columns]
  }))
}
utils::write.table(bind(result$routes), file.path(out, "g0-routes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA")
if (length(result$proteins)) utils::write.table(do.call(rbind, result$proteins), file.path(out, "g0-proteins.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
utils::write.table(data.frame(protein = rownames(env$g0Matrix(shift)), env$g0Matrix(shift), check.names = FALSE),
  file.path(out, "g0-matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(provenance, sprintf("finished\t%s", format(Sys.time(), tz = "UTC"))), file.path(out, "provenance.txt"))
cat(readLines(file.path(out, "g0-routes.tsv")), sep = "\n")
