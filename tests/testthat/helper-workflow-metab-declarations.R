# KIT-032: a deterministic metabolomics fixture (the KIT-031 lipid fixture's twin) large enough for every declared step (the manual RUV-III mode, the
# per-group correlation pairs, limma over two assays): 24 samples (2 groups x 6 biological x 2 technical replicates,
# 3 batches), two assays of 40 metabolites with internal-standard rows and duplicate ids, dense intensities on the linear
# scale. Used by test-workflow-metab-declarations.R.

metabDeclarationsDesign <- function() {
  groups <- rep(c("WT", "KO"), each = 12L)
  bio <- rep(rep(seq_len(6L), each = 2L), times = 2L)
  tech <- rep(seq_len(2L), times = 12L)
  samples <- sprintf("%s_%d_%d", groups, bio, tech)
  data.frame(Run = samples, group = groups, batch = rep(c("B1", "B2", "B3"), length.out = 24L),
    replicates = bio, tech_rep_group = paste(groups, bio, sep = "_"), stringsAsFactors = FALSE)
}

metabDeclarationsAssay <- function(seed, prefix, design = metabDeclarationsDesign(), n = 40L) {
  set.seed(seed)
  ids <- sprintf("%s_M%02d", prefix, seq_len(n))
  classes <- rep(c("PC", "PE", "TG", "SM", "CE"), length.out = n)
  base <- 2^stats::runif(n, 12, 20)
  effect <- ifelse(seq_len(n) %% 5L == 0L, 2.5, 1)
  values <- vapply(seq_len(nrow(design)), function(j) {
    scale <- 2^stats::rnorm(1, 0, 0.2)
    batch <- 2^(as.integer(factor(design$batch[[j]])) * 0.15)
    noise <- 2^stats::rnorm(n, 0, 0.25)
    base * ifelse(design$group[[j]] == "KO", effect, 1) * scale * batch * noise
  }, numeric(n))
  colnames(values) <- design$Run
  standards <- matrix(2^stats::rnorm(4L * nrow(design), 16, 0.1), nrow = 4L, dimnames = list(NULL, design$Run))
  table <- data.frame(Metabolite = c(ids, paste0("IS_", prefix, "_", seq_len(4L)), ids[c(3L, 7L)]),
    Ontology = c(classes, rep("ISTD", 4L), classes[c(3L, 7L)]), stringsAsFactors = FALSE)
  cbind(table, as.data.frame(rbind(values, standards, values[c(3L, 7L), , drop = FALSE] * 0.6), check.names = FALSE))
}

metabDeclarationsObject <- function() {
  design <- metabDeclarationsDesign()
  methods::new("MetaboliteAssayData",
    metabolite_data = list(LCMS_Pos = metabDeclarationsAssay(31L, "pos", design), GCMS = metabDeclarationsAssay(32L, "gc", design)),
    design_matrix = design, metabolite_id_column = "Metabolite", annotation_id_column = "Ontology",
    sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group",
    database_identifier_type = "Unknown", internal_standard_regex = "^IS_",
    args = list(deAnalysisParameters = list(formula_string = "~ 0 + group", de_q_val_thresh = 0.05, treat_lfc_cutoff = 0)))
}

metabDeclarationsContrasts <- function() {
  data.frame(contrasts = "groupKO-groupWT", friendly_names = "KO_vs_WT", full_format = "KO_vs_WT=groupKO-groupWT",
    stringsAsFactors = FALSE)
}

#' A metabolomics session owner in a temporary project with the import's source copies and the design file on disk, the
#' design state committed with its provenance and inputs (what the design step does on the owner).
metabDeclarationsOwner <- function(root, object = metabDeclarationsObject(), config = list(qc = list(mode = "declarations"))) {
  source_dir <- file.path(root, "source"); dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(source_dir, "import_sources"), showWarnings = FALSE)
  for (assay in names(object@metabolite_data)) {
    utils::write.table(object@metabolite_data[[assay]], file.path(source_dir, "import_sources", paste0(assay, "_export.txt")),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  utils::write.table(object@design_matrix, file.path(source_dir, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
  paths <- list(source_dir = source_dir, metabolite_qc_dir = file.path(root, "qc"), da_output_dir = file.path(root, "da"),
    publication_graphs_dir = file.path(root, "pub"))
  for (p in paths[-1]) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  session <- omicsOwnerSession()
  workflow <- list2env(list(state_update_trigger = NULL, config_list = config, contrasts_tbl = metabDeclarationsContrasts(),
    tab_status = list(quality_control = "pending", normalization = "pending")))
  workflow$state_manager <- createDIASessionOwner(workflow, paths, session, omic_type = "metabolomics")
  shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, "METABOLOMICS", "msdial"))
  owner <- workflow$state_manager
  recordMetabDesignProvenance(owner, object, source_dir = source_dir)
  owner$saveState("metab_raw_data_s4", object, config, "Initial MetaboliteAssayData S4 object created after design matrix")
  list(owner = owner, workflow = workflow, session = session, paths = paths, source_dir = source_dir, object = object)
}

metabDeclarationsRuvParameters <- function() {
  list(ruv_mode = "manual", ruv_grouping_variable = "group", percentage_min = 1, percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 1L, manual_percentage = 20)
}

#' The whole metabolomics journey on the owner through the wrappers the modules call (the GUI's order and state names).
metabDeclarationsJourney <- function(fixture, ruv = metabDeclarationsRuvParameters()) {
  owner <- fixture$owner; config <- fixture$workflow$config_list
  out <- list()
  out$intensity <- runMetabReferenceIntensityFilter(owner, 1, 0.5, config, hydrate = FALSE)
  out$duplicates <- runMetabReferenceDuplicateResolution(owner, config, hydrate = FALSE)
  commitWorkflowOmicsStoreStep(owner, "metab_qc_finalize", "metab_qc_complete", config, "QC processing complete - ready for normalization")
  out$itsd <- runMetabReferenceItsdNormalisation(owner, "median", NULL, config, hydrate = FALSE)
  out$log2 <- runMetabReferenceLog2(owner, 1, config, hydrate = FALSE)
  out$between <- runMetabReferenceBetweenSampleNormalisation(owner, "cyclicloess", config, hydrate = FALSE)
  out$ruv <- runMetabReferenceRuv(owner, ruv, config, hydrate = FALSE)
  out$correlation <- runMetabReferenceCorrelation(owner, "group", 0.5, config, hydrate = FALSE)
  out$da <- runMetabReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
    list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$workflow$contrasts_tbl)
  out
}
