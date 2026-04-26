# testthat for Proteomics Volcano Plot
# Phase 4 of Proteomics GUI Test Strategy
# fidelity-coverage-compare: shared

if (!methods::isClass("MockProtDaVolcanoObject")) {
  methods::setClass(
    "MockProtDaVolcanoObject",
    slots = c(
      protein_quant_table = "data.frame",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character"
    )
  )
}

build_prot_volcano_long <- function(comparison = "T_vs_C", ids = c("P1", "P2", "P3", "P4")) {
  data.frame(
    Protein.Ids = ids,
    comparison = rep(comparison, length(ids)),
    log2FC = c(2.5, -3.0, 0.5, -0.2)[seq_along(ids)],
    fdr_qvalue = c(0, 0.005, 0.5, 0.8)[seq_along(ids)],
    raw_pvalue = c(0.0001, 0.0005, 0.1, 0.2)[seq_along(ids)],
    custom.display = paste0("display-", seq_along(ids)),
    stringsAsFactors = FALSE
  )
}

expect_prot_glimma_widget <- function(result) {
  expect_s3_class(result, "htmlwidget")
  expect_s3_class(result, "glimmaXY")
}

localProtVolcanoBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

test_that("prepareDataForVolcanoPlot covers protein and metabolite inputs", {
  protein_ready <- prepareDataForVolcanoPlot(
    input_table = data.frame(
      Protein.Ids = c("P1:isoform", "P2", "P3"),
      q.mod = c(0.01, 0.02, 0.5),
      log2FC = c(1.5, -2, 0.2),
      stringsAsFactors = FALSE
    ),
    protein_id_column = Protein.Ids,
    uniprot_table = data.frame(
      uniprot_acc_first = c("P1", "P2", "P3"),
      gene_name = c("GENE1 full", "GENE2 alternate", "GENE3"),
      stringsAsFactors = FALSE
    ),
    uniprot_protein_id_column = uniprot_acc_first,
    gene_name_column = gene_name,
    number_of_genes = 2,
    fdr_threshold = 0.05,
    fdr_column = q.mod,
    log2FC_column = log2FC
  )
  expect_equal(protein_ready$colour, c("red", "blue", "grey"))
  expect_equal(as.character(protein_ready$label), c("Significant Increase", "Significant Decrease", "Not significant"))
  expect_equal(protein_ready$gene_name_significant[1:2], c("GENE1", "GENE2"))

  metabolite_ready <- prepareDataForVolcanoPlot(
    input_table = data.frame(
      Name = c("M1", "M2", "M3"),
      q.mod = c(0.01, 0.02, 0.5),
      log2FC = c(1.5, -2, 0.2),
      stringsAsFactors = FALSE
    ),
    protein_id_column = Name,
    uniprot_table = data.frame(),
    number_of_genes = 1,
    fdr_threshold = 0.05,
    fdr_column = q.mod,
    log2FC_column = log2FC
  )
  expect_equal(metabolite_ready$colour, c("red", "blue", "grey"))
  expect_equal(metabolite_ready$gene_name_significant[1], "M1")
})

test_that("writeInteractiveVolcanoPlotProteomics covers static and widget writers", {
  captured <- new.env(parent = emptyenv())
  captured$static <- list()
  captured$widgets <- list()

  localProtVolcanoBinding(
    environment(writeInteractiveVolcanoPlotProteomics),
    "getGlimmaVolcanoProteomics",
    function(volcano_plot_tab, contrast_name, counts_tbl = NULL, groups = NULL, output_dir = NULL, ...) {
      captured$static[[length(captured$static) + 1L]] <- list(
        contrast_name = contrast_name,
        row_count = nrow(volcano_plot_tab),
        counts_dim = dim(counts_tbl),
        groups = groups,
        output_dir = output_dir
      )
      invisible(TRUE)
    }
  )
  localProtVolcanoBinding(
    environment(writeInteractiveVolcanoPlotProteomicsWidget),
    "getGlimmaVolcanoProteomicsWidget",
    function(fit.eb, coef, volcano_plot_tab, ...) {
      captured$widgets[[length(captured$widgets) + 1L]] <- list(
        coef = coef,
        labels = volcano_plot_tab$label,
        genes = volcano_plot_tab$gene_name
      )
      list(coef = coef, rows = nrow(volcano_plot_tab))
    }
  )

  da_long <- data.frame(
    Protein.Ids = c("P1-1", "P2", "P3"),
    uniprot_acc = c("P1", "P2", "P3"),
    comparison = c("A_vs_B", "A_vs_B", "B_vs_C"),
    fdr_qvalue = c(0.01, 0.2, 0.03),
    raw_pvalue = c(0.001, 0.05, 0.01),
    log2FC = c(1.4, -0.2, -1.5),
    stringsAsFactors = FALSE
  )
  uniprot_tbl <- data.frame(
    Entry = c("P1", "P2", "P3"),
    gene_names = c("GENE1 full", "GENE2", "GENE3"),
    gene_names_column = c("GENE1 full", "", NA_character_),
    stringsAsFactors = FALSE
  )

  writeInteractiveVolcanoPlotProteomics(
    da_proteins_long = da_long,
    uniprot_tbl = uniprot_tbl,
    publication_graphs_dir = tempdir(),
    args_row_id = "missing_id",
    counts_tbl = matrix(seq_len(6), nrow = 3, dimnames = list(c("P1-1", "P2", "P3"), c("S1", "S2"))),
    groups = c(S1 = "A", S2 = "B")
  )
  expect_equal(vapply(captured$static, `[[`, character(1), "contrast_name"), c("A_vs_B", "B_vs_C"))
  expect_true(dir.exists(captured$static[[1]]$output_dir))
  expect_equal(captured$static[[1]]$counts_dim, c(2L, 2L))

  widget_result <- writeInteractiveVolcanoPlotProteomicsWidget(
    da_proteins_long = da_long,
    uniprot_tbl = uniprot_tbl,
    fit.eb = list(coefficients = matrix(1:4, ncol = 2)),
    args_row_id = "uniprot_acc",
    gene_names_column = "gene_names_column"
  )
  expect_length(widget_result, 2)
  expect_equal(vapply(captured$widgets, `[[`, integer(1), "coef"), 1:2)
  expect_true("Sig., logFC >= 1" %in% captured$widgets[[1]]$labels)
})

test_that("generateProtDAVolcanoPlotGlimma works with captured snapshot", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp08_volcano_input.rds")
  
  if (file.exists(cp_file)) {
    data <- tryCatch(
      readRDS(cp_file),
      error = function(e) {
        skip(sprintf("Snapshot cp08 is unreadable: %s", conditionMessage(e)))
      }
    )
    
    # Run the function
    # Note: we need to handle potential NULLs in captured data if any
    p <- do.call(generateProtDAVolcanoPlotGlimma, data)
    
    # Assertions
    expect_s3_class(p, "htmlwidget")
    expect_s3_class(p, "glimmaXY")
  } else {
    skip("Snapshot cp08 not found")
  }
})

test_that("generateProtDAVolcanoPlotGlimma preserves guard exits and fallback contrast behavior", {
  expect_null(generateProtDAVolcanoPlotGlimma(NULL, selected_contrast = "T_vs_C"))
  expect_null(generateProtDAVolcanoPlotGlimma(list(), selected_contrast = "T_vs_C"))
  expect_null(generateProtDAVolcanoPlotGlimma(
    da_results_list = list(da_proteins_long = build_prot_volcano_long()),
    selected_contrast = NULL
  ))

  invalid_rows <- build_prot_volcano_long(ids = c("", NA, "P3", "P4"))
  invalid_rows$log2FC <- c(NA_real_, NA_real_, NA_real_, NA_real_)
  expect_null(generateProtDAVolcanoPlotGlimma(
    da_results_list = list(da_proteins_long = invalid_rows),
    selected_contrast = "T_vs_C"
  ))

  fallback_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(da_proteins_long = build_prot_volcano_long("available_contrast")),
    selected_contrast = "missing_contrast",
    display_columns = "custom.display"
  )
  expect_prot_glimma_widget(fallback_result)
})

test_that("generateProtDAVolcanoPlotGlimma preserves fuzzy, prefix, and UniProt annotation paths", {
  fuzzy_rows <- build_prot_volcano_long(
    comparison = "T-vs-C",
    ids = c("sp|P11111-1|ALPHA", "tr|P22222|BETA", "lcl|P33333", "P44444")
  )
  fuzzy_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(da_proteins_long = fuzzy_rows),
    selected_contrast = "T_vs_C",
    uniprot_tbl = data.frame(
      Entry = c("P11111", "P22222", "P33333"),
      gene_names = c("GENE1 alternate", "GENE2;ALT", "GENE3:ALT"),
      stringsAsFactors = FALSE
    ),
    display_columns = c("best_uniprot_acc", "custom.display")
  )
  expect_prot_glimma_widget(fuzzy_result)

  prefix_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(da_proteins_long = build_prot_volcano_long("Dose=High_vs_Low")),
    selected_contrast = "Dose",
    uniprot_tbl = data.frame(
      Accession = "not-a-supported-id-column",
      Symbol = "not-a-supported-gene-column",
      stringsAsFactors = FALSE
    )
  )
  expect_prot_glimma_widget(prefix_result)
})

test_that("generateProtDAVolcanoPlotGlimma preserves S4 count synchronization paths", {
  protein_object <- methods::new(
    "MockProtDaVolcanoObject",
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2", "P3", "P4"),
      "Sample A" = c(11, 21, 31, 41),
      "Sample B" = c(12, 22, 32, 42),
      check.names = FALSE
    ),
    design_matrix = data.frame(
      Run = c("sample a", "sample b"),
      group = c("Control", "Treatment"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group"
  )

  synchronized_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(
      da_proteins_long = build_prot_volcano_long(),
      theObject = protein_object
    ),
    selected_contrast = "T_vs_C",
    args_row_id = "Protein.Ids"
  )
  expect_prot_glimma_widget(synchronized_result)

  missing_group_object <- protein_object
  missing_group_object@group_id <- "missing_group"
  missing_group_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(
      da_proteins_long = build_prot_volcano_long(),
      theObject = missing_group_object
    ),
    selected_contrast = "T_vs_C",
    args_row_id = "Protein.Ids"
  )
  expect_prot_glimma_widget(missing_group_result)

  unmatched_samples_object <- protein_object
  unmatched_samples_object@design_matrix <- data.frame(
    Run = c("Other A", "Other B"),
    group = c("Control", "Treatment"),
    stringsAsFactors = FALSE
  )
  unmatched_samples_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(
      da_proteins_long = build_prot_volcano_long(),
      theObject = unmatched_samples_object
    ),
    selected_contrast = "T_vs_C",
    args_row_id = "Protein.Ids"
  )
  expect_prot_glimma_widget(unmatched_samples_result)

  missing_count_id_object <- protein_object
  missing_count_id_object@protein_quant_table <- data.frame(
    OtherId = c("P1", "P2", "P3", "P4"),
    SampleA = c(11, 21, 31, 41),
    SampleB = c(12, 22, 32, 42),
    stringsAsFactors = FALSE
  )
  missing_count_id_result <- generateProtDAVolcanoPlotGlimma(
    da_results_list = list(
      da_proteins_long = build_prot_volcano_long(),
      theObject = missing_count_id_object
    ),
    selected_contrast = "T_vs_C",
    args_row_id = "Protein.Ids"
  )
  expect_prot_glimma_widget(missing_count_id_result)
})

test_that("generateProtDAVolcanoPlotGlimma handles mock data", {
  # Mock CP08 input
  mock_cp08 <- list(
    da_results_list = list(
      da_proteins_long = data.frame(
        Protein.Ids = c("P1", "P2", "P3", "P4"),
        comparison = rep("T_vs_C", 4),
        log2FC = c(2.5, -3.0, 0.5, -0.2),
        fdr_qvalue = c(0.001, 0.005, 0.5, 0.8),
        raw_pvalue = c(0.0001, 0.0005, 0.1, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    selected_contrast = "T_vs_C",
    da_q_val_thresh = 0.05,
    args_row_id = "Protein.Ids"
  )
  
  # Run function
  result <- do.call(generateProtDAVolcanoPlotGlimma, mock_cp08)
  
  # Assertions
  expect_s3_class(result, "htmlwidget")
  expect_s3_class(result, "glimmaXY")
})

test_that("generateProtDAVolcanoStatic handles mock data", {
  # Mock inputs
  da_results_list <- list(
    da_proteins_long = data.frame(
      Protein.Ids = paste0("P", 1:100),
      comparison = rep("T_vs_C", 100),
      log2FC = c(2.2, -2.4, rep(0.2, 98)),
      fdr_qvalue = c(0.01, 0.02, rep(0.5, 98)),
      raw_pvalue = c(0.001, 0.002, rep(0.2, 98)),
      stringsAsFactors = FALSE
    )
  )
  
  # Run function
  p <- generateProtDAVolcanoStatic(
    da_results_list = da_results_list,
    selected_contrast = "T_vs_C",
    da_q_val_thresh = 0.05,
    lfc_threshold = 1.0
  )
  
  # Assertions
  expect_s3_class(p, "ggplot")

  fallback_id_plot <- generateProtDAVolcanoStatic(
    da_results_list = list(
      da_proteins_long = data.frame(
        feature_id = c("F1", "F2"),
        comparison = "T_vs_C",
        log2FC = c(0.2, -0.3),
        fdr_qvalue = c(0.8, 0.9),
        raw_pvalue = c(0.8, 0.9),
        stringsAsFactors = FALSE
      )
    ),
    selected_contrast = "T_vs_C"
  )
  expect_s3_class(fallback_id_plot, "ggplot")
})

test_that("generateProtDAVolcanoStatic preserves guard, contrast parsing, and no-label paths", {
  expect_null(generateProtDAVolcanoStatic(NULL, selected_contrast = "T_vs_C"))
  expect_null(generateProtDAVolcanoStatic(list(), selected_contrast = "T_vs_C"))
  expect_null(generateProtDAVolcanoStatic(
    list(da_proteins_long = build_prot_volcano_long()),
    selected_contrast = NULL
  ))

  parsed_plot <- generateProtDAVolcanoStatic(
    list(
      da_proteins_long = data.frame(
        uniprot_acc = c("P1", "P2", "P3"),
        gene_name = c("GENE1", "", NA_character_),
        comparison = "Dose",
        log2FC = c(1.5, -1.5, 0.2),
        fdr_qvalue = c(0.01, 0.02, 0.5),
        raw_pvalue = c(0.001, 0.002, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    selected_contrast = "Dose=High_vs_Low",
    show_labels = FALSE
  )
  expect_s3_class(parsed_plot, "ggplot")
  expect_identical(parsed_plot$labels$title, "Volcano Plot: Dose")

  localProtVolcanoBinding(
    parent.env(environment(generateProtDAVolcanoStatic)),
    "fixed",
    stringr::fixed
  )

  expect_null(generateProtDAVolcanoStatic(
    list(da_proteins_long = build_prot_volcano_long("available")),
    selected_contrast = "missing"
  ))

  fuzzy_plot <- generateProtDAVolcanoStatic(
    list(
      da_proteins_long = data.frame(
        uniprot_acc = c("P1", "P2", "P3"),
        gene_name = c("GENE1", "", NA_character_),
        comparison = "DoseA",
        log2FC = c(1.5, -1.5, 0.2),
        fdr_qvalue = c(0.01, 0.02, 0.5),
        raw_pvalue = c(0.001, 0.002, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    selected_contrast = "Dose",
    show_labels = FALSE
  )
  expect_s3_class(fuzzy_plot, "ggplot")
  expect_identical(fuzzy_plot$labels$title, "Volcano Plot: DoseA")
})

# APAF Bioinformatics | test-prot-08-volcano.R | Approved | 2026-03-13
