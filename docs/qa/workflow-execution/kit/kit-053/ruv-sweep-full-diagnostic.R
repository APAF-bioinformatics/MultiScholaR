# KIT-053: diagnostic on the halted full-corpus recorder store (the route on its shipped defaults: norm_method = none).
# The RUV sweep found no valid percentage; this script hydrates the correlation_filtered state, shows the per-sample offset
# by group and the control counts before and after cyclic loess. Usage: KIT053_INSTALLED=1 R_LIBS=<library> Rscript <this> (after
# oracle-record.R full has run with norm_method = none and halted at protein_ruv_covariates).
suppressMessages(library(MultiScholaR))
repo <- normalizePath("."); ns <- asNamespace("MultiScholaR")
for (f in ls(ns)) if (!exists(f, envir = globalenv(), inherits = TRUE) && is.function(get(f, ns))) assign(f, get(f, ns), envir = globalenv())
harness <- new.env(parent = globalenv()); sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), harness)
bindings <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
store <- harness$kitFindStore(file.path(bindings$harness$root, "oracles", "limpa_full", "project")); store_root <- dirname(store)
identity <- readSavedStoreIdentity(file.path(store_root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
withDIAStateScope(manager, function() {
  obj <- manager$getState("correlation_filtered")
  design <- obj@design_matrix; samples <- as.character(design[[obj@sample_id]])
  Y <- as.matrix(tibble::column_to_rownames(as.data.frame(obj@protein_quant_table), obj@protein_id_column))[, samples]
  g <- design$group[match(samples, design[[obj@sample_id]])]
  cm <- colMeans(Y); cat("column means by group:\n"); print(tapply(cm, g, function(x) round(c(mean = mean(x), sd = sd(x)), 3)))
  gm <- t(apply(Y, 1, function(y) tapply(y, g, mean))); cat("sign of (treated - control) across proteins:", table(sign(gm[, "treated"] - gm[, "control"])), "\n")
  t0 <- Sys.time()
  norm <- suppressMessages(normaliseBetweenSamples(obj, normalisation_method = "cyclicloess"))
  cat("cyclicloess", round(as.numeric(Sys.time() - t0), 1), "s\n")
  Yn <- as.matrix(tibble::column_to_rownames(as.data.frame(norm@protein_quant_table), norm@protein_id_column))[, samples]
  ps <- apply(Yn, 1, function(y) summary(stats::aov(y ~ factor(g)))[[1]][["Pr(>F)"]][1]); bh <- p.adjust(ps, "BH")
  cat("after loess: BH q>0.05:", sum(bh > 0.05), " median p", signif(median(ps), 3), "\n")
  for (p in c(1, 4, 10, 20, 30)) { ctl <- getNegCtrlProtAnova(norm, ruv_grouping_variable = "group", percentage_as_neg_ctrl = p, ruv_qval_cutoff = 0.05, ruv_fdr_method = "BH"); cat("p", p, "controls", sum(ctl, na.rm = TRUE), "\n") }
})
