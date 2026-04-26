# fidelity-coverage-compare: shared
library(testthat)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

test_that("general volcano plotting helpers preserve current ggplot assembly", {
  volcano_data <- data.frame(
    logFC = c(2, -1.5, 0.2),
    lqm = c(4, 3, 0.5),
    label = c("Up", "Down", "Neutral"),
    colour = c("red", "blue", "grey"),
    analysis_type = c("proteomics", "proteomics", "proteomics"),
    comparison = c("A_vs_B", "A_vs_B", "A_vs_B"),
    stringsAsFactors = FALSE
  )

  expect_s3_class(plotOneVolcano(volcano_data, input_title = "Volcano"), "ggplot")
  expect_s3_class(plotOneVolcanoNoVerticalLines(volcano_data, input_title = "Volcano"), "ggplot")
  expect_s3_class(
    plotVolcano(
      selected_data = transform(
        volcano_data,
        colour = factor(colour, levels = c("red", "blue", "grey"))
      ),
      formula_string = "analysis_type ~ comparison"
    ),
    "ggplot"
  )
})

test_that("general volcano helpers preserve protein labeling and Glimma export routing", {
  volcano_input <- data.frame(
    Protein.Ids = c("P1-1", "P2", "P3"),
    best_uniprot_acc = c("P1", "P2", "P3"),
    gene_name = c("GENE1", "GENE2", ""),
    PROTEIN_NAMES = c("Protein 1", "Protein 2", "Protein 3"),
    fdr_qvalue = c(0.01, 0.02, 0.5),
    log2FC = c(2.1, -1.8, 0.1),
    raw_pvalue = c(0.001, 0.002, 0.2),
    stringsAsFactors = FALSE
  )

  uniprot_table <- data.frame(
    Entry = c("P1", "P2", "P3"),
    gene_name = c("GENE1", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )

  labelled_plot <- printOneVolcanoPlotWithProteinLabel(
    input_table = volcano_input,
    uniprot_table = uniprot_table,
    protein_id_column = Protein.Ids,
    uniprot_protein_id_column = Entry,
    gene_name_column = gene_name,
    number_of_genes = 2,
    include_protein_label = TRUE
  )
  expect_s3_class(labelled_plot, "ggplot")

  skip_if_not_installed("Glimma")
  loadNamespace("Glimma")

  glimma_namespace <- asNamespace("Glimma")
  output_dir <- tempfile("glimma-prot-", tmpdir = getwd())
  dir.create(output_dir)
  withr::defer(unlink(output_dir, recursive = TRUE, force = TRUE))

  localNamespaceBinding(
    glimma_namespace,
    "glimmaXY",
    function(..., html) {
      writeLines("<html><head></head><body>mock glimma</body></html>", html)
      structure(list(html = html), class = "glimmaXY")
    }
  )

  getGlimmaVolcanoProteomics(
    volcano_plot_tab = volcano_input,
    uniprot_column = best_uniprot_acc,
    gene_name_column = gene_name,
    display_columns = "PROTEIN_NAMES",
    counts_tbl = NULL,
    groups = c("A", "B"),
    output_dir = output_dir,
    contrast_name = "A_vs_B"
  )

  expect_true(file.exists(file.path(output_dir, "A_vs_B.html")))
})
