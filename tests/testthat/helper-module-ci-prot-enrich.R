module_ci_prot_enrich_design <- function(samples = c("S1", "S2", "S3", "S4")) {
  design <- data.frame(
    Run = samples,
    group = c("A", "A", "B", "B"),
    batch = c("Batch1", "Batch2", "Batch1", "Batch2"),
    stringsAsFactors = FALSE
  )
  rownames(design) <- design$Run
  design
}

module_ci_prot_enrich_da_table <- function(identifier_case = "uniprot",
                                           result_case = "standard") {
  base <- data.frame(
    uniprot_acc = c("P04637", "P10600", "Q07817", "P55957", "O14757", "P38398"),
    gene_name = c("TP53", "APP", "BCL2L1", "BID", "CHEK2", "BRCA1"),
    ensembl_id = c("ENSG00000141510", "ENSG00000142192", "ENSG00000171552", "ENSG00000015475", "ENSG00000183765", "ENSG00000012048"),
    log2FC = c(2.1, -1.8, 0.7, 1.4, -2.2, 0.2),
    raw_pvalue = c(0.001, 0.002, 0.04, 0.03, 0.004, 0.5),
    fdr_qvalue = c(0.005, 0.006, 0.049, 0.04, 0.01, 0.8),
    stringsAsFactors = FALSE
  )

  base$log2norm.S1.A <- c(10, 11, 12, 13, 14, 15)
  base$log2norm.S2.A <- base$log2norm.S1.A + 0.2
  base$log2norm.S3.B <- base$log2norm.S1.A + c(2, -2, 0.5, 1.5, -2, 0)
  base$log2norm.S4.B <- base$log2norm.S3.B + 0.2

  base <- switch(identifier_case,
    uniprot = base,
    gene_symbol = transform(base, uniprot_acc = gene_name),
    ensembl_like = transform(base, uniprot_acc = ensembl_id),
    duplicate_ids = {
      base$uniprot_acc[2] <- base$uniprot_acc[1]
      base
    },
    missing_annotation = transform(base, gene_name = c("TP53", NA, "", "BID", NA, "BRCA1")),
    mixed_case = transform(base, uniprot_acc = tolower(uniprot_acc), gene_name = c("Tp53", "app", "Bcl2l1", "bid", "Chek2", "brca1")),
    no_mappable_ids = transform(base, uniprot_acc = paste0("NO_MAP_", seq_len(nrow(base))), gene_name = paste0("NO_MAP_GENE_", seq_len(nrow(base)))),
    empty = base[0, , drop = FALSE],
    stop(sprintf("Unknown proteomics enrichment identifier fixture: %s", identifier_case), call. = FALSE)
  )

  switch(result_case,
    standard = base,
    no_hits = transform(base, fdr_qvalue = rep(0.75, nrow(base))),
    all_hits = transform(base, fdr_qvalue = seq(0.001, 0.006, length.out = nrow(base))),
    threshold_boundary = transform(base, fdr_qvalue = c(0, 0.001, 0.05, 0.1, 1, 0.049), log2FC = c(0, 1, -1, 2, -2, 0.5)),
    base
  )
}

module_ci_prot_enrich_contrasts <- function(raw = "groupB-groupA",
                                            friendly = "B_vs_A") {
  data.frame(
    contrasts = raw,
    full_format = paste0(friendly, "=", raw),
    friendly_names = friendly,
    stringsAsFactors = FALSE
  )
}

module_ci_prot_enrich_object <- function(protein_id_column = "uniprot_acc",
                                         da_table = module_ci_prot_enrich_da_table(),
                                         design = module_ci_prot_enrich_design()) {
  id_values <- if (protein_id_column %in% names(da_table)) {
    da_table[[protein_id_column]]
  } else {
    da_table$uniprot_acc
  }
  protein_quant_table <- data.frame(
    feature_id = id_values,
    S1 = seq_len(nrow(da_table)) + 10,
    S2 = seq_len(nrow(da_table)) + 11,
    S3 = seq_len(nrow(da_table)) + 12,
    S4 = seq_len(nrow(da_table)) + 13,
    stringsAsFactors = FALSE
  )
  names(protein_quant_table)[[1]] <- protein_id_column
  protein_id_table <- data.frame(
    feature_id = id_values,
    gene_names = da_table$gene_name,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(protein_id_table)[[1]] <- protein_id_column

  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = protein_quant_table,
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    protein_id_column = protein_id_column,
    protein_id_table = protein_id_table,
    args = list(
      globalParameters = list(workflow_type = "DIA", report_template = "DIA_report"),
      enrichmentAnalysisUI = list()
    )
  )
}

module_ci_prot_enrich_da_results <- function(identifier_case = "uniprot",
                                             result_case = "standard",
                                             contrast = "groupB-groupA") {
  methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = contrast),
    da_data = stats::setNames(
      list(module_ci_prot_enrich_da_table(identifier_case, result_case)),
      contrast
    ),
    design_matrix = module_ci_prot_enrich_design()
  )
}

module_ci_prot_enrich_annotations <- function(identifier_case = "uniprot",
                                              organism = "human") {
  da_table <- module_ci_prot_enrich_da_table(identifier_case)
  organism_name <- switch(organism,
    human = "Homo sapiens (Human)",
    mouse = "Mus musculus (Mouse)",
    custom = "Synthetic test organism",
    "Homo sapiens (Human)"
  )
  taxon_id <- switch(organism,
    human = "9606",
    mouse = "10090",
    custom = "999999",
    "9606"
  )

  data.frame(
    Entry = toupper(da_table$uniprot_acc),
    gene_names = toupper(da_table$gene_name),
    Organism = organism_name,
    OX = paste0("ID=", taxon_id),
    go_id_go_biological_process = "GO:0006915; GO:0007049",
    go_term_go_biological_process = "apoptotic process; cell cycle",
    stringsAsFactors = FALSE
  )
}

module_ci_prot_enrich_method <- function(backend = "gprofiler2",
                                         taxon = "9606") {
  if (identical(backend, "gprofiler2")) {
    return(list(
      method = "gprofiler2",
      supported = TRUE,
      species_id = "hsapiens",
      species_name = "Homo sapiens"
    ))
  }

  list(
    method = "clusterprofiler",
    supported = FALSE,
    species_id = NULL,
    species_name = paste("Taxon ID", taxon)
  )
}

module_ci_prot_enrich_gprofiler_table <- function(direction = c("positive", "negative"),
                                                  result_case = "standard",
                                                  contrast = "B_vs_A",
                                                  taxon = "9606") {
  direction <- match.arg(direction)
  table <- buildProtEnrichDeterministicGprofilerResult(direction)$result
  table$Description <- table$term_name
  table$pvalue <- table$p_value
  table$p.adjust <- table$p_value
  table$qvalue <- table$p_value
  table$Count <- table$term_size
  table$data_source <- table$source
  table$selected_contrast <- contrast
  table$organism_taxid <- taxon
  table$up_cutoff <- 1
  table$down_cutoff <- 0
  table$q_cutoff <- 0.05
  if (identical(result_case, "empty")) {
    table <- table[0, , drop = FALSE]
  }
  if (identical(result_case, "high_q")) {
    table$pvalue <- table$p.adjust <- table$qvalue <- 0.2
  }
  table
}

module_ci_prot_enrich_cluster_table <- function(direction = c("up", "down"),
                                                result_case = "standard",
                                                contrast = "B_vs_A",
                                                taxon = "999999") {
  direction <- match.arg(direction)
  table <- buildProtEnrichDeterministicClusterProfilerResult(direction)@result
  table$selected_contrast <- contrast
  table$organism_taxid <- taxon
  table$up_cutoff <- 1
  table$down_cutoff <- 0
  table$q_cutoff <- 0.05
  if (identical(result_case, "empty")) {
    table <- table[0, , drop = FALSE]
  }
  if (identical(result_case, "high_q")) {
    table$pvalue <- table$p.adjust <- table$qvalue <- 0.2
  }
  table
}

module_ci_prot_enrich_results <- function(backend = "gprofiler2",
                                          result_case = "standard",
                                          contrast = "groupB-groupA") {
  enrichment_results <- createEnrichmentResults(tibble::tibble(contrasts = contrast))

  if (identical(backend, "gprofiler2")) {
    enrichment_results@enrichment_data <- stats::setNames(
      list(list(
        up = list(result = module_ci_prot_enrich_gprofiler_table("positive", result_case)),
        down = list(result = module_ci_prot_enrich_gprofiler_table("negative", result_case))
      )),
      contrast
    )
  } else {
    enrichment_results@enrichment_data <- stats::setNames(
      list(list(
        up = methods::new("ProtEnrichDeterministicResult", result = module_ci_prot_enrich_cluster_table("up", result_case)),
        down = methods::new("ProtEnrichDeterministicResult", result = module_ci_prot_enrich_cluster_table("down", result_case))
      )),
      contrast
    )
  }

  enrichment_results@enrichment_plotly <- stats::setNames(
    list(list(
      up = buildProtEnrichDeterministicPlot(backend, "up"),
      down = buildProtEnrichDeterministicPlot(backend, "down")
    )),
    contrast
  )
  enrichment_results@enrichment_plots <- stats::setNames(list(list(up = NULL, down = NULL)), contrast)
  enrichment_results@enrichment_summaries <- stats::setNames(
    list(list(
      up = list(backend = backend, direction = "up"),
      down = list(backend = backend, direction = "down")
    )),
    contrast
  )

  enrichment_results
}

module_ci_prot_enrich_paths <- function() {
  root <- tempfile("module-ci-prot-enrich-")
  paths <- list(
    results_dir = file.path(root, "results"),
    results_summary_dir = file.path(root, "results_summary"),
    publication_graphs_dir = file.path(root, "publication_graphs"),
    time_dir = file.path(root, "time"),
    qc_dir = file.path(root, "qc"),
    da_output_dir = file.path(root, "da_output"),
    pathway_dir = file.path(root, "pathway_enrichment"),
    source_dir = file.path(root, "source"),
    feature_qc_dir = file.path(root, "feature_qc")
  )
  lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)
  paths
}

module_ci_prot_enrich_write_pathway_tables <- function(paths,
                                                       backend = "gprofiler2",
                                                       contrast = "B_vs_A") {
  if (identical(backend, "gprofiler2")) {
    up <- module_ci_prot_enrich_gprofiler_table("positive", contrast = contrast, taxon = "9606")
    down <- module_ci_prot_enrich_gprofiler_table("negative", contrast = contrast, taxon = "9606")
  } else {
    up <- module_ci_prot_enrich_cluster_table("up", contrast = contrast, taxon = "999999")
    down <- module_ci_prot_enrich_cluster_table("down", contrast = contrast, taxon = "999999")
  }

  readr::write_tsv(up, file.path(paths$pathway_dir, paste0(contrast, "_up_enrichment_results.tsv")))
  readr::write_tsv(down, file.path(paths$pathway_dir, paste0(contrast, "_down_enrichment_results.tsv")))
  invisible(paths$pathway_dir)
}

module_ci_prot_enrich_archive <- function(backend = "gprofiler2",
                                          result_case = "standard") {
  archive <- tempfile(fileext = ".zip")
  temp_dir <- tempfile("module-ci-prot-enrich-archive-")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  method_info <- module_ci_prot_enrich_method(
    backend,
    taxon = if (identical(backend, "gprofiler2")) "9606" else "999999"
  )

  if (identical(backend, "gprofiler2")) {
    gprofiler_results <- module_ci_prot_enrich_gprofiler_table("positive", result_case)
    cluster_results <- NULL
  } else {
    gprofiler_results <- NULL
    cluster_results <- module_ci_prot_enrich_cluster_table("up", result_case)
  }

  writeProtEnrichResultsDownloadArchive(
    file = archive,
    selectedContrast = "B_vs_A",
    methodInfo = method_info,
    organismTaxid = if (identical(backend, "gprofiler2")) "9606" else "999999",
    upCutoff = 1,
    downCutoff = 0,
    qCutoff = 0.05,
    gprofilerResults = gprofiler_results,
    clusterprofilerResults = cluster_results,
    stringdbResults = NULL,
    tempDir = temp_dir
  )

  list(
    archive = archive,
    listing = utils::unzip(archive, list = TRUE),
    temp_dir = temp_dir
  )
}

module_ci_prot_enrich_datatable_capture <- function(data, ...) {
  structure(list(data = data, args = list(...)), class = "module_ci_datatable")
}

module_ci_prot_enrich_format_round_capture <- function(table, columns, digits) {
  attr(table, "rounded_columns") <- columns
  attr(table, "rounded_digits") <- digits
  table
}

module_ci_prot_enrich_workflow <- function(object = module_ci_prot_enrich_object(),
                                           da_results = list(`groupB-groupA` = module_ci_prot_enrich_da_table())) {
  manager <- new.env(parent = emptyenv())
  manager$current_state <- "normalized"
  manager$getState <- function(state_name) {
    if (identical(state_name, "normalized")) {
      return(object)
    }
    NULL
  }

  workflow <- new.env(parent = emptyenv())
  workflow$state_manager <- manager
  workflow$da_analysis_results_list <- da_results
  workflow$tab_status <- list(enrichment_analysis = "pending")
  workflow$taxon_id <- "9606"
  workflow$mixed_species_analysis <- list(enabled = FALSE)
  workflow
}
