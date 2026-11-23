# KIT-172: fixture F172 and the declared-operation journeys of the plain-log2 tests (test-kit172-*.R) and of the
# ticket's scripts (docs/qa/workflow-execution/kit/kit-172/). Fixture F172 (ticket section 6.4): 520 features of
# which 5 are internal-standard rows (`IS_`), two groups of 6 samples, feature levels log-uniform so that the values
# span 0.01 to 1e6, log2-normal noise (SD 0.25) in every cell, 2% of the observed cells exactly 1 (at most 2 per
# feature and group), 10% of the cells missing at random (at most 3 of 6 per feature and group). Both classes are
# built from the same draw. Every journey runs through the declared operations on a session owner, the default
# (partitioned) variant where the step has one, as the application runs them.

kit172Design <- function() {
  runs <- c(sprintf("A_%d", 1:6), sprintf("B_%d", 1:6))
  data.frame(Run = runs, group = rep(c("A", "B"), each = 6L), batch = "B1", replicates = c(1:6, 1:6),
    tech_rep_group = runs, stringsAsFactors = FALSE)
}

kit172Contrasts <- function() {
  data.frame(contrasts = "groupB-groupA", friendly_names = "B_vs_A", full_format = "B_vs_A=groupB-groupA",
    stringsAsFactors = FALSE)
}

#' Cells (row, column) taken in random order, at most `per_group` per feature and group, `count` in all, among `allowed`.
kit172PickCells <- function(allowed, groups, per_group, count) {
  cells <- which(allowed, arr.ind = TRUE)
  cells <- cells[sample.int(nrow(cells)), , drop = FALSE]
  key <- paste(cells[, 1L], groups[cells[, 2L]])
  rank <- stats::ave(seq_along(key), key, FUN = seq_along)
  cells <- cells[rank <= per_group, , drop = FALSE]
  cells[seq_len(min(count, nrow(cells))), , drop = FALSE]
}

#' F172's values: a numeric matrix (features by the design's runs) and the planted truth.
#' @param planted T1b: one unwanted factor W (-1, -1, -1, +1, +1, +1 within each group, loadings N(0, 0.5^2) on log2)
#'   and true effects (10% of the features at +-1 log2 in group B).
#' @param missing Whether 10% of the cells are missing.
kit172Values <- function(seed, planted = FALSE, missing = TRUE, n = 520L, n_itsd = 5L) {
  set.seed(seed)
  design <- kit172Design()
  groups <- design$group
  n_features <- n - n_itsd
  level <- stats::runif(n_features, log2(0.01), log2(1e6))
  logs <- matrix(level, n_features, nrow(design)) + matrix(stats::rnorm(n_features * nrow(design), 0, 0.25), n_features)
  effects <- rep(0, n_features)
  loadings <- rep(0, n_features)
  if (planted) {
    changed <- sample.int(n_features, round(0.1 * n_features))
    effects[changed] <- sample(c(-1, 1), length(changed), replace = TRUE)
    logs[, groups == "B"] <- logs[, groups == "B"] + effects
    loadings <- stats::rnorm(n_features, 0, 0.5)
    logs <- logs + outer(loadings, rep(c(-1, -1, -1, 1, 1, 1), 2L))
  }
  values <- 2^logs
  if (missing) {
    gone <- kit172PickCells(matrix(TRUE, n_features, nrow(design)), groups, 3L, round(0.1 * length(values)))
    values[gone] <- NA_real_
  }
  ones <- kit172PickCells(!is.na(values), groups, 2L, round(0.02 * sum(!is.na(values))))
  values[ones] <- 1
  standards <- matrix(2^(13 + stats::rnorm(n_itsd * nrow(design), 0, 0.1)), n_itsd)
  out <- rbind(values, standards)
  dimnames(out) <- list(c(sprintf("F%03d", seq_len(n_features)), sprintf("IS_%d", seq_len(n_itsd))), design$Run)
  list(values = out, effects = stats::setNames(c(effects, rep(0, n_itsd)), rownames(out)),
    loadings = stats::setNames(c(loadings, rep(0, n_itsd)), rownames(out)))
}

kit172Binding <- function(class) if (identical(class, "lipid")) lipidomicsDeclarationBinding else metabolomicsDeclarationBinding

#' The assay table of a values matrix for one class (one assay, `LCMS_Pos`).
kit172Table <- function(values, class) {
  ids <- rownames(values)
  annotation <- ifelse(startsWith(ids, "IS_"), "ISTD", "PC")
  table <- if (identical(class, "lipid")) data.frame(LipidName = ids, LipidClass = annotation, stringsAsFactors = FALSE) else
    data.frame(Metabolite = ids, Ontology = annotation, stringsAsFactors = FALSE)
  cbind(table, as.data.frame(unname(values), stringsAsFactors = FALSE) |> stats::setNames(colnames(values)))
}

#' The S4 object of a values matrix (`values` may be a named list of matrices, one per assay).
kit172Object <- function(values, class) {
  if (!is.list(values)) values <- list(LCMS_Pos = values)
  assays <- lapply(values, kit172Table, class = class)
  design <- kit172Design()
  args <- list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
  if (identical(class, "lipid")) {
    methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = design, lipid_id_column = "LipidName",
      annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  } else {
    methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = design, metabolite_id_column = "Metabolite",
      annotation_id_column = "Ontology", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
      database_identifier_type = "Unknown", internal_standard_regex = "^IS_", args = args)
  }
}

#' A session owner holding `object` as its design state (the declarations fixtures' owner).
kit172Owner <- function(root, object, class) {
  fixture <- if (identical(class, "lipid")) lipidDeclarationsOwner(root, object) else metabDeclarationsOwner(root, object)
  fixture$workflow$contrasts_tbl <- kit172Contrasts()
  fixture
}

#' The application's default log2 offset (the tab's configuration): 0 from KIT-172 on, 1 before.
kit172DefaultOffset <- function(class) {
  config <- if (identical(class, "lipid")) omicGuiLipidNorm() else omicGuiMetabNorm()
  config$panels$log2$offset
}

kit172Quiet <- function(expr) suppressWarnings(suppressMessages(expr))

#' One journey through the declared operations: intensity filter (percentile 1, proportion 0.5; optional), QC
#' finalise, ITSD (median, the `IS_` standards; optional), log2, between-sample normalisation, RUV (optional), DA.
#' @return The owner, the fixture and every step's result.
kit172Journey <- function(root, object, class, itsd = TRUE, method = "cyclicloess", filter = TRUE,
    offset = kit172DefaultOffset(class), ruv = NULL, da = TRUE) {
  fixture <- kit172Owner(root, object, class)
  owner <- fixture$owner
  config <- fixture$workflow$config_list
  binding <- kit172Binding(class)
  prefix <- binding$prefix
  out <- list(fixture = fixture, owner = owner, binding = binding, prefix = prefix)
  kit172Quiet({
    if (filter) out$filter <- runWorkflowOmicsReferenceIntensityFilter(binding, owner, 1, 0.5, config)
    commitWorkflowOmicsStoreStep(owner, paste0(prefix, "_qc_finalize"), paste0(prefix, "_qc_complete"), config, "QC complete")
    if (itsd) out$itsd <- runWorkflowOmicsReferenceItsdNormalisation(binding, owner, "median", NULL, config)
    out$log2 <- runWorkflowOmicsReferenceLog2(binding, owner, offset, config)
    out$between <- runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, method, config)
    if (!is.null(ruv)) out$ruv <- runWorkflowOmicsReferenceRuv(binding, owner, ruv, config)
    if (da) out$da <- runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner),
      list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit172Contrasts())
  })
  out
}

#' A state's assay as a matrix (rows by id, the design's runs as columns).
kit172StateMatrix <- function(journey, state, assay = "LCMS_Pos") {
  object <- journey$owner$getState(state)
  table <- methods::slot(object, journey$binding$data_slot)[[assay]]
  id <- methods::slot(object, journey$binding$id_slot)
  runs <- kit172Design()$Run
  matrix(as.matrix(as.data.frame(table)[, runs, drop = FALSE]), ncol = length(runs), dimnames = list(table[[id]], runs))
}

#' The DA long table of a journey, with the id column named `feature`.
kit172DaTable <- function(journey) {
  long <- as.data.frame(journey$da$results[[journey$binding$long_table]])
  id <- intersect(c("lipid_id", "metabolite_id"), names(long))[[1L]]
  long$feature <- long[[id]]
  long
}

#' The legacy RUV-III parameters at the panel's defaults (mode `automatic` or `manual`).
kit172LegacyRuv <- function(mode = "automatic") {
  list(ruv_mode = mode, ruv_grouping_variable = "group", percentage_min = 1, percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2L,
    manual_percentage = 5)
}

#' F172 (seed 172001) mapped onto 0.01 to 50 (log-affine), the cells that were exactly 1 set to 1 again, and five
#' features with no value in group A.
kit172ConcentrationValues <- function() {
  values <- kit172Values(172001L)$values
  features <- !startsWith(rownames(values), "IS_")
  ones <- values == 1 & !is.na(values)
  logs <- log2(values[features, ])
  range <- range(logs, na.rm = TRUE)
  values[features, ] <- 2^(log2(0.01) + (logs - range[[1L]]) / diff(range) * (log2(50) - log2(0.01)))
  values[ones] <- 1
  values[sprintf("F%03d", 10:14), kit172Design()$group == "A"] <- NA_real_
  values
}
