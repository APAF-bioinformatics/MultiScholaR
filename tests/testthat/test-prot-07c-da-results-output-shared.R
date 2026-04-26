# fidelity-coverage-compare: shared
library(testthat)

getSharedProtDaResultsOutput <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("outputDaResultsAllContrasts", envir = package_ns, inherits = FALSE)) {
    return(get("outputDaResultsAllContrasts", envir = package_ns, inherits = FALSE))
  }

  outputDaResultsAllContrasts
}

localSharedProtDaOutputNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

newSharedProtDaOutputCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$vroom_files <- character()
  captured$excel_files <- character()
  captured$plot_files <- character()
  captured$volcano_titles <- character()
  captured
}

withSharedProtDaOutputMocks <- function(captured, code) {
  mock_frame <- parent.frame()
  package_ns <- asNamespace("MultiScholaR")

  localSharedProtDaOutputNamespaceBinding(
    package_ns,
    "plotOneVolcanoNoVerticalLines",
    function(data, title, ...) {
      captured$volcano_titles <- c(captured$volcano_titles, title)
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point()
    },
    mock_frame
  )

  testthat::local_mocked_bindings(
    vroom_write = function(x, file, ...) {
      dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
      utils::write.table(x, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
      captured$vroom_files <- c(captured$vroom_files, basename(file))
      invisible(x)
    },
    .package = "vroom",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    write_xlsx = function(x, path, ...) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      writeLines("mock xlsx", path)
      captured$excel_files <- c(captured$excel_files, basename(path))
      invisible(path)
    },
    .package = "writexl",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    ggsave = function(filename, plot = ggplot2::last_plot(), ...) {
      dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      writeLines("mock plot", filename)
      captured$plot_files <- c(captured$plot_files, basename(filename))
      invisible(filename)
    },
    .package = "ggplot2",
    .env = mock_frame
  )

  force(code)
}

captureSharedProtDaOutputMessages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  list(value = value, messages = messages)
}

makeSharedProtDaObject <- function(args = list()) {
  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = data.frame(
      uniprot_acc = c("P1", "P2", "P3"),
      S1 = c(10, 11, 12),
      S2 = c(20, 21, 22),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "uniprot_acc",
    protein_id_table = data.frame(
      uniprot_acc = c("P1", "P2", "P3"),
      gene_names = c("GENE1", "GENE2", "GENE3"),
      stringsAsFactors = FALSE
    ),
    args = args
  )
}

makeSharedDaProteinRows <- function(include_gene_symbol = FALSE) {
  rows <- data.frame(
    uniprot_acc = c("P1-1", "P2-2", "P3-3"),
    log2FC = c(1.6, -1.2, 0.2),
    fdr_qvalue = c(0.01, 0.2, 0.03),
    raw_pvalue = c(0.001, 0.04, 0.02),
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_gene_symbol)) {
    rows$Gene.Symbol <- c("ALT1:description", "", NA_character_)
  }

  rows
}

makeSharedNumSigTable <- function(significant = TRUE) {
  if (isTRUE(significant)) {
    return(data.frame(
      status = c("Up", "Down", "Not significant"),
      counts = c(2, 1, 0),
      comparison = c("A_vs_B", "A_vs_B", "A_vs_B"),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    status = "Not significant",
    counts = 0,
    comparison = "A_vs_B",
    stringsAsFactors = FALSE
  )
}

makeSharedProtDaSummaryInput <- function() {
  data.frame(
    logFC = c(1.5, -0.6, 0.25, -1.2),
    fdr_qvalue = c(0.01, 0.02, 0.4, 0.001),
    raw_pvalue = c(0.001, 0.01, 0.2, 0.0005),
    stringsAsFactors = FALSE,
    row.names = paste0("P", 1:4)
  )
}

test_that("proteomics DA summary helpers count and format significant data", {
  summary_input <- makeSharedProtDaSummaryInput()

  counts <- countStatDaGenes(
    summary_input,
    lfc_thresh = 0,
    q_val_thresh = 0.05,
    log_fc_column = logFC,
    q_value_column = fdr_qvalue
  )
  expect_equal(counts$status, c("Not significant", "Significant and Up", "Significant and Down"))
  expect_equal(counts$counts, c(1, 1, 2))

  helper_result <- countStatDaGenesHelper(
    da_table = list("A_vs_B=log_intensity" = summary_input),
    description = "limma",
    facet_column = analysis_type,
    comparison_column = "comparison",
    expression_column = "expression",
    log_fc_column = logFC,
    q_value_column = fdr_qvalue
  )
  expect_equal(unique(helper_result$analysis_type), "limma")
  expect_equal(unique(helper_result$comparison), "A_vs_B")
  expect_equal(unique(helper_result$expression), "log_intensity")

  table_result <- printCountDaGenesTable(
    list(list("A_vs_B=log_intensity" = summary_input)),
    list("limma"),
    formula_string = NA_character_
  )
  expect_s3_class(table_result$plot, "ggplot")
  expect_equal(nrow(table_result$table), 3)

  significant <- getSignificantData(
    list(list("A_vs_B=log_intensity" = summary_input)),
    list("limma"),
    row_id = uniprot_acc,
    p_value_column = raw_pvalue,
    q_value_column = fdr_qvalue,
    log_fc_column = logFC,
    comparison_column = "comparison",
    expression_column = "expression",
    facet_column = analysis_type,
    q_val_thresh = 0.05
  )
  expect_equal(significant$uniprot_acc, rownames(summary_input))
  expect_equal(unique(significant$comparison), "A_vs_B")
  expect_equal(unique(significant$expression), "log_intensity")
  expect_s3_class(significant$colour, "factor")
  expect_true(all(c("black", "orange", "blue", "purple") %in% levels(significant$colour)))
})

test_that("proteomics DA long-format helpers preserve grouping and result extraction", {
  design <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  norm_counts <- data.frame(
    S1 = c(10, 20),
    S2 = c(11, 21),
    S3 = c(30, 40),
    S4 = c(31, 41),
    row.names = c("P1", "P2"),
    check.names = FALSE
  )
  raw_counts <- norm_counts * 10
  lfc_qval_tbl <- data.frame(
    uniprot_acc = c("P1", "P2"),
    lqm = c(2, 1),
    colour = c("purple", "black"),
    analysis_type = "limma",
    comparison = "A_vs_B",
    log_intensity = "groupA-groupB",
    fdr_qvalue = c(0.01, 0.2),
    raw_pvalue = c(0.001, 0.05),
    log2FC = c(1.5, -0.2),
    stringsAsFactors = FALSE
  )

  long_result <- createDaResultsLongFormat(
    lfc_qval_tbl = lfc_qval_tbl,
    norm_counts_input_tbl = norm_counts,
    raw_counts_input_tbl = raw_counts,
    row_id = "uniprot_acc",
    sample_id = "Run",
    group_id = "group",
    group_pattern = "S",
    design_matrix_norm = design,
    design_matrix_raw = design,
    expression_column = log_intensity,
    protein_id_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      gene_names = c("GENE1", "GENE2"),
      stringsAsFactors = FALSE
    )
  )

  expect_equal(long_result$numerator, c("A", "A"))
  expect_equal(long_result$denominator, c("B", "B"))
  expect_true(all(c("log2norm.S1.A", "raw.S3.B", "gene_names") %in% names(long_result)))

  grouping <- getTypeOfGrouping(design, group_id = "group", sample_id = "Run")
  expect_equal(grouping$A, c("S1", "S2"))
  expect_equal(grouping$B, c("S3", "S4"))

  extracted <- extractResults(list(
    contrast_a = list(results = data.frame(value = 1)),
    contrast_b = list(results = data.frame(value = 2))
  ))
  expect_named(extracted, c("contrast_a", "contrast_b"))
  expect_equal(extracted$contrast_b$value, 2)
})

test_that("protein DA output method writes annotated results and summary graphics", {
  output_fn <- getSharedProtDaResultsOutput()
  captured <- newSharedProtDaOutputCapture()
  output_dir <- file.path(tempdir(), "prot-da-output")
  graphs_dir <- file.path(tempdir(), "prot-da-graphs")
  unlink(c(output_dir, graphs_dir), recursive = TRUE, force = TRUE)

  da_results <- list(
    "A vs B" = list(
      da_proteins_long = makeSharedDaProteinRows(),
      num_sig_da_molecules_first_go = list(table = makeSharedNumSigTable(significant = TRUE))
    )
  )
  uniprot_tbl <- data.frame(
    Entry = c("P1", "P2", "P3"),
    gene_names = c("GENE1 full", "", NA_character_),
    stringsAsFactors = FALSE
  )
  prot_obj <- makeSharedProtDaObject(
    args = list(outputDaResultsAllContrasts = list(
      da_q_val_thresh = 0.05,
      fdr_column = "fdr_qvalue",
      log2fc_column = "log2FC"
    ))
  )

  result <- withSharedProtDaOutputMocks(captured, {
    captureSharedProtDaOutputMessages(
      output_fn(
        theObject = prot_obj,
        da_results_list_all_contrasts = da_results,
        uniprot_tbl = uniprot_tbl,
        da_output_dir = output_dir,
        publication_graphs_dir = graphs_dir,
        file_prefix = "da_proteins",
        args_row_id = "uniprot_acc",
        gene_names_column = "gene_names",
        uniprot_id_column = "Entry"
      )
    )
  })

  expect_true(result$value)
  expect_true("da_proteins_A_vs_B_long_annot.tsv" %in% captured$vroom_files)
  expect_true("da_proteins_A_vs_B_long_annot.xlsx" %in% captured$excel_files)
  expect_true("da_proteins_num_sig_da_molecules.tab" %in% captured$vroom_files)
  expect_true("da_proteins_num_sig_da_molecules.xlsx" %in% captured$excel_files)
  expect_true("A vs B" %in% captured$volcano_titles)
  expect_true(all(c("A_vs_B.png", "A_vs_B.pdf", "da_proteins_num_sig_da_molecules.png") %in% captured$plot_files))
  expect_true(any(grepl("Created combined PDF", result$messages, fixed = TRUE)))
  expect_true(any(grepl("Wrote NumSigDE table", result$messages, fixed = TRUE)))
})

test_that("protein DA output method covers fallback annotations and no-signal summaries", {
  output_fn <- getSharedProtDaResultsOutput()
  captured <- newSharedProtDaOutputCapture()
  output_dir <- file.path(tempdir(), "prot-da-output-fallback")
  graphs_dir <- file.path(tempdir(), "prot-da-graphs-fallback")
  unlink(c(output_dir, graphs_dir), recursive = TRUE, force = TRUE)

  da_results <- list(
    "C:D" = list(
      da_proteins_long = makeSharedDaProteinRows(include_gene_symbol = TRUE),
      num_sig_da_molecules_first_go = list(table = makeSharedNumSigTable(significant = FALSE))
    ),
    "empty" = list(da_proteins_long = NULL)
  )

  prot_obj <- makeSharedProtDaObject(
    args = list(outputDaResultsAllContrasts = list(
      da_q_val_thresh = 0.05,
      fdr_column = "fdr_qvalue",
      log2fc_column = "log2FC"
    ))
  )
  result <- withSharedProtDaOutputMocks(captured, {
    captureSharedProtDaOutputMessages(
      output_fn(
        theObject = prot_obj,
        da_results_list_all_contrasts = da_results,
        uniprot_tbl = data.frame(Entry = character(), stringsAsFactors = FALSE),
        da_output_dir = output_dir,
        publication_graphs_dir = graphs_dir,
        file_prefix = "fallback",
        args_row_id = "uniprot_acc",
        gene_names_column = "Gene.Symbol",
        uniprot_id_column = "Entry"
      )
    )
  })

  expect_true(result$value)
  expect_true("fallback_C_D_long_annot.tsv" %in% captured$vroom_files)
  expect_true("fallback_C_D_long_annot.xlsx" %in% captured$excel_files)
  expect_true("fallback_num_sig_da_molecules.tab" %in% captured$vroom_files)
  expect_true("C:D" %in% captured$volcano_titles)
  expect_true(any(grepl("No significant DE molecules found", result$messages, fixed = TRUE)))
})

test_that("protein DA output method handles contrast lists without output tables", {
  output_fn <- getSharedProtDaResultsOutput()
  captured <- newSharedProtDaOutputCapture()
  output_dir <- file.path(tempdir(), "prot-da-output-empty")
  graphs_dir <- file.path(tempdir(), "prot-da-graphs-empty")
  unlink(c(output_dir, graphs_dir), recursive = TRUE, force = TRUE)

  da_results <- list(skipped = list(da_proteins_long = NULL))
  result <- withSharedProtDaOutputMocks(captured, {
    captureSharedProtDaOutputMessages(
      output_fn(
        theObject = makeSharedProtDaObject(),
        da_results_list_all_contrasts = da_results,
        uniprot_tbl = data.frame(),
        da_output_dir = output_dir,
        publication_graphs_dir = graphs_dir,
        file_prefix = "empty",
        args_row_id = "uniprot_acc",
        gene_names_column = "missing_gene_names",
        uniprot_id_column = "Entry"
      )
    )
  })

  expect_true(result$value)
  expect_length(captured$vroom_files, 0L)
  expect_length(captured$volcano_titles, 0L)
  expect_true(any(grepl("No NumSigDE data found", result$messages, fixed = TRUE)))
})
