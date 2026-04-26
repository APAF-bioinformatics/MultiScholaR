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

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

if (!methods::isClass("mockEnrichmentResultsCarrier")) {
  methods::setClass(
    "mockEnrichmentResultsCarrier",
    slots = c(
      enrichment_data = "list",
      enrichment_plotly = "list",
      enrichment_summaries = "list"
    )
  )
}

localFakeClusterProfilerGeneral <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  old_options <- options(
    multischolar.fake_clusterProfiler.enricher = NULL,
    multischolar.fake_clusterProfiler.GSEA = NULL,
    multischolar.fake_clusterProfiler.enricher_result = NULL,
    multischolar.fake_clusterProfiler.gsea_result = NULL
  )

  fake_lib <- file.path(tempdir(), "multischolar-fake-clusterprofiler-general-lib")
  fake_pkg <- file.path(tempdir(), "multischolar-fake-clusterprofiler-general-src")
  installed_pkg <- file.path(fake_lib, "clusterProfiler")

  if ("clusterProfiler" %in% loadedNamespaces()) {
    try(unloadNamespace("clusterProfiler"), silent = TRUE)
  }

  if (!dir.exists(installed_pkg)) {
    unlink(fake_pkg, recursive = TRUE)
    dir.create(file.path(fake_pkg, "R"), recursive = TRUE, showWarnings = FALSE)
    dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

    writeLines(
      c(
        "Package: clusterProfiler",
        "Version: 0.0.1",
        "Title: Fake clusterProfiler test double",
        "Description: Minimal namespace used by MultiScholaR enrichment coverage tests.",
        "License: MIT",
        "Encoding: UTF-8",
        "Depends: methods"
      ),
      file.path(fake_pkg, "DESCRIPTION")
    )
    writeLines(
      c(
        "export(enricher,GSEA)",
        "exportClasses(enrichResult,gseaResult)"
      ),
      file.path(fake_pkg, "NAMESPACE")
    )
    writeLines(
      c(
        "methods::setClass('enrichResult', slots = c(result = 'data.frame'))",
        "methods::setClass('gseaResult', slots = c(result = 'data.frame'))",
        "",
        ".call_handler <- function(handler, args) {",
        "  formal_names <- names(formals(handler))",
        "  if (is.null(formal_names) || '...' %in% formal_names) {",
        "    return(do.call(handler, args))",
        "  }",
        "  do.call(handler, args[intersect(names(args), formal_names)])",
        "}",
        "",
        "enricher <- function(gene, universe = NULL, TERM2GENE = NULL, TERM2NAME = NULL, ...) {",
        "  handler <- getOption('multischolar.fake_clusterProfiler.enricher')",
        "  args <- list(gene = gene, universe = universe, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)",
        "  if (is.function(handler)) {",
        "    return(.call_handler(handler, args))",
        "  }",
        "  result <- getOption('multischolar.fake_clusterProfiler.enricher_result')",
        "  if (is.null(result)) {",
        "    result <- data.frame()",
        "  }",
        "  methods::new('enrichResult', result = result)",
        "}",
        "",
        "GSEA <- function(geneList, TERM2GENE = NULL, TERM2NAME = NULL, ...) {",
        "  handler <- getOption('multischolar.fake_clusterProfiler.GSEA')",
        "  args <- list(geneList = geneList, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, ...)",
        "  if (is.function(handler)) {",
        "    return(.call_handler(handler, args))",
        "  }",
        "  result <- getOption('multischolar.fake_clusterProfiler.gsea_result')",
        "  if (is.null(result)) {",
        "    result <- data.frame()",
        "  }",
        "  methods::new('gseaResult', result = result)",
        "}"
      ),
      file.path(fake_pkg, "R", "clusterProfiler.R")
    )
    utils::install.packages(fake_pkg, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
  }

  .libPaths(c(fake_lib, old_lib_paths))

  withr::defer({
    options(old_options)
    if ("clusterProfiler" %in% loadedNamespaces()) {
      try(unloadNamespace("clusterProfiler"), silent = TRUE)
    }
    .libPaths(old_lib_paths)
  }, envir = env)

  invisible(fake_lib)
}

test_that("general enrichment helper conversions preserve lookup and grouping behavior", {
  annotation_tbl <- data.frame(
    annotation_name = c("Term A", "Term B", "Term A"),
    annotation_id = c("GO:1", "GO:2", "GO:1"),
    protein_id = c("P1", "P2", "P3"),
    group = c("Positive", "Positive", "Negative"),
    gene_name = c("GeneA alpha", "GeneB beta", "GeneC gamma"),
    stringsAsFactors = FALSE
  )

  dict <- buildAnnotationIdToAnnotationNameDictionary(
    annotation_tbl,
    annotation_name,
    annotation_id
  )

  expect_equal(convertIdToAnnotation("GO:1", dict), "Term A")
  expect_error(convertIdToAnnotation("GO:missing", dict), "subscript out of bounds")
  expect_equal(convertProteinAccToGeneSymbol(c("P1", "P9"), c(P1 = "GeneA")), "GeneA/NA")

  gene_sets <- buildOneProteinToAnnotationList(annotation_tbl, annotation_id, protein_id)
  expect_equal(gene_sets[["GO:1"]], c("P1", "P3"))

  grouped_tables <- listifyTableByColumn(annotation_tbl, group)
  expect_equal(sort(names(grouped_tables)), c("Negative", "Positive"))
  expect_equal(nrow(grouped_tables$Positive), 2L)

  lookup <- getUniprotAccToGeneSymbolDictionary(
    annotation_tbl,
    protein_id,
    gene_name,
    protein_id
  )
  expect_equal(lookup[["P1"]], "GeneA")
})

test_that("GO and gene-set helper pipelines preserve filtered enrichment structure", {
  one_go_enrichment <- makeFunctionWithOverrides(
    oneGoEnrichment,
    list(
      enricher = function(...) {
        data.frame(
          ID = "GO:1",
          Description = "placeholder",
          GeneRatio = "2/3",
          BgRatio = "3/3",
          pvalue = 0.01,
          p.adjust = 0.02,
          qvalue = 0.02,
          geneID = "P1/P2",
          Count = 2,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  run_gsea <- makeFunctionWithOverrides(
    runGsea,
    list(
      geneIds = function(x) x$gene_sets,
      GSEA = function(geneList, TERM2GENE) {
        list(gene_list = geneList, term_to_gene = TERM2GENE)
      }
    )
  )
  run_enricher_override <- makeFunctionWithOverrides(
    runEnricher,
    list(
      geneIds = function(x) x$gene_sets,
      enricher = function(gene, TERM2GENE, ...) {
        data.frame(
          ID = "SetA",
          Count = length(unique(TERM2GENE$gene)),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  go_annot <- data.frame(
    uniprot_acc = c("P1", "P2", "P3", "P1", "P2"),
    go_id = c("GO:1", "GO:1", "GO:1", "GO:2", "GO:2"),
    go_type = c("BP", "BP", "BP", "BP", "BP"),
    stringsAsFactors = FALSE
  )
  background <- data.frame(uniprot_acc = c("P1", "P2", "P3"), stringsAsFactors = FALSE)
  output <- one_go_enrichment(
    go_annot = go_annot,
    background_list = background,
    go_aspect = "BP",
    query_list = c("P1", "P2", "P3"),
    id_to_annotation_dictionary = c(`GO:1` = "Term A", `GO:2` = "Term B"),
    annotation_id = go_id,
    protein_id = uniprot_acc,
    aspect_column = go_type,
    p_val_thresh = 0.05,
    min_gene_set_size = 2,
    max_gene_set_size = 4,
    get_cluster_profiler_object = TRUE
  )

  expect_equal(output$output_table$term[[1]], "Term A")
  expect_identical(output$cluster_profiler_object$ID[[1]], "GO:1")

  gsea_input <- list(
    contrast_a = c(P1 = 1.2, P2 = 0.8, P3 = 0.3),
    contrast_b = c("P1", "P2", "P3")
  )
  gene_sets <- list(
    hallmark = list(gene_sets = list(SetA = c("P1", "P2"), SetB = c("P2", "P3", "P4")))
  )

  gsea_result <- run_gsea(
    "hallmark",
    "contrast_a",
    gsea_input,
    gene_sets,
    min_set_size = 2,
    max_set_size = 3
  )
  expect_equal(sort(unique(gsea_result$term_to_gene$term)), c("SetA", "SetB"))
  expect_equal(sort(unique(gsea_result$term_to_gene$gene)), c("P1", "P2", "P3"))

  enricher_result <- run_enricher_override(
    "hallmark",
    "contrast_b",
    gsea_input,
    gene_sets,
    min_set_size = 2,
    max_set_size = 3
  )
  expect_equal(enricher_result$Count[[1]], 3L)
})

test_that("plotting and duplicate helpers preserve expected enrichment summaries", {
  input_tbl <- data.frame(
    term = c("Path A", "Path A", "Path B"),
    qvalue = c(0.01, 0.02, 0.03),
    gene_set = c("Positive genes", "Negative genes", "Both genes"),
    stringsAsFactors = FALSE
  )

  cleaned <- cleanDuplicatesEnrichment(input_tbl)
  expect_equal(nrow(cleaned), 2L)
  expect_true("neg_log_10_fdr" %in% names(cleaned))

  plot_obj <- plotEnrichmentBarplot(input_tbl)
  expect_s3_class(plot_obj, "ggplot")

  expect_equal(
    get_param_change_message(
      "foldChange",
      data.frame(
        listname = "color.params",
        present = "foldChange",
        original = "foldChange",
        row.names = "foldChange",
        stringsAsFactors = FALSE
      )
    ),
    "Use 'color.params = list(foldChange = your_value)' instead of 'foldChange'.\n The foldChange parameter will be removed in the next version."
  )

  p <- list(data = data.frame(name = c("cat", "gene1", "gene2"), stringsAsFactors = FALSE))
  alpha_updated <- node_add_alpha(
    p,
    hilight_category = "cat",
    hilight_gene = "gene2",
    alpha_nohilight = 0.2,
    alpha_hilight = 0.9
  )
  expect_equal(alpha_updated$data$alpha, c(0.9, 0.2, 0.9))

  expect_equal(get_enrichplot_color(2), c("#e06663", "#327eba"))
  expect_equal(get_enrichplot_color(3), c("#e06663", "white", "#327eba"))
  expect_error(get_enrichplot_color(4), "'n' should be 2 or 3")

  scale_two <- set_enrichplot_color(colors = c("#111111", "#999999"))
  scale_three <- set_enrichplot_color(colors = c("#111111", "#555555", "#999999"))
  expect_s3_class(scale_two, "Scale")
  expect_s3_class(scale_three, "Scale")
  expect_equal(get_ggrepel_segsize(default = 0.4), 0.4)

  graph_tbl <- list2df(list(A = c("P1", "P2"), B = "P3"))
  expect_equal(nrow(graph_tbl), 3L)
  list_to_graph <- makeFunctionWithOverrides(
    list2graph,
    list(graph.data.frame = igraph::graph_from_data_frame)
  )
  expect_true(igraph::gorder(list_to_graph(list(A = c("P1", "P2"), B = "P3"))) >= 3L)
})

test_that("UniProt and contrast aggregation helpers preserve caching and alias behavior", {
  package_ns <- asNamespace("MultiScholaR")
  cache_root <- tempfile("general-enrichment-cache-")
  dir.create(cache_root, recursive = TRUE)
  withr::defer(unlink(cache_root, recursive = TRUE, force = TRUE))

  localNamespaceBindings(
    package_ns,
    list(
      batchQueryEvidenceGeneId = function(...) {
        data.frame(
          Entry = c("P1", "P2"),
          Gene.Ontology.IDs = c("GO:1", "GO:2"),
          Protein.existence = c("1", "2"),
          Protein.names = c("Prot 1", "Prot 2"),
          stringsAsFactors = FALSE
        )
      },
      uniprotGoIdToTerm = function(.data, ...) .data,
      enrichProteinsPathwaysHelper = function(da_analysis_results, ...) {
        data.frame(
          ID = "GO:1",
          directionality = "positive",
          stringsAsFactors = FALSE
        )
      }
    )
  )

  protein_ids <- data.frame(uniprot_acc = c("P1", "P2"), stringsAsFactors = FALSE)
  cache_file <- file.path(cache_root, "uniprot-cache.rds")

  downloaded <- download_uniprot_data(
    protein_ids = protein_ids,
    cache_file = cache_file,
    uniprot_handle = list()
  )
  expect_true(file.exists(cache_file))
  expect_true(all(c("Protein_existence", "Protein_names") %in% names(downloaded)))

  localNamespaceBinding(
    package_ns,
    "batchQueryEvidenceGeneId",
    function(...) stop("batchQueryEvidenceGeneId should not be called when cache is complete")
  )
  cached_only <- download_uniprot_data(
    protein_ids = data.frame(uniprot_acc = "P1", stringsAsFactors = FALSE),
    cache_file = cache_file,
    uniprot_handle = list()
  )
  expect_equal(nrow(cached_only), 2L)

  enrich_tbl <- enrichProteinsPathways(
    da_analysis_results_list = NULL,
    de_analysis_results_list = list(
      contrast_a = list(),
      contrast_b = list()
    ),
    taxon_id = "9606"
  )
  expect_equal(sort(unique(enrich_tbl$comparison)), c("contrast_a", "contrast_b"))
})

test_that("UniProt GO parsing, gprofiler execution, and result accessors preserve behavior", {
  parsed <- uniprotGoIdToTermSimple(
    data.frame(
      UNIPROTKB = c("P1", "P2"),
      GO.IDs = c("GO:1; GO:2", "GO:3"),
      Gene.Names = c("GeneA", "GeneB"),
      stringsAsFactors = FALSE
    ),
    go_id_column = GO.IDs,
    sep = "; ",
    goterms = c(`GO:1` = "Term 1", `GO:2` = "Term 2", `GO:3` = "Term 3"),
    gotypes = c(`GO:1` = "BP", `GO:2` = "MF", `GO:3` = "CC")
  )

  expect_equal(nrow(parsed), 3L)
  expect_equal(sort(unique(parsed$go_type)), c("Biological Process", "Cellular Compartment", "Molecular Function"))

  gost_calls <- new.env(parent = emptyenv())
  gost_calls$count <- 0L
  gost_calls$queries <- list()
  local_mocked_bindings(
    gost = function(query, custom_bg, ...) {
      gost_calls$count <- gost_calls$count + 1L
      gost_calls$queries[[gost_calls$count]] <- list(query = query, custom_bg = custom_bg)
      if (gost_calls$count == 1L) {
        stop("transient failure")
      }
      list(result = data.frame(term_id = "GO:1", stringsAsFactors = FALSE))
    },
    .package = "gprofiler2"
  )

  success <- suppressWarnings(
    perform_enrichment(
      data_subset = data.frame(protein_id = c("P1", NA, "", "P2"), stringsAsFactors = FALSE),
      species = "hsapiens",
      threshold = 0.05,
      sources = c("GO:BP"),
      domain_scope = "annotated",
      custom_bg = c("B1", NA, "B2"),
      exclude_iea = FALSE,
      max_retries = 2,
      wait_time = 0,
      protein_id_column = "protein_id"
    )
  )

  expect_equal(gost_calls$count, 2L)
  expect_equal(gost_calls$queries[[2]]$query, c("P1", "P2"))
  expect_equal(gost_calls$queries[[2]]$custom_bg, c("B1", "B2"))
  expect_true(!is.null(success$result))
  expect_null(
    suppressWarnings(
      perform_enrichment(
        data_subset = data.frame(protein_id = character(), stringsAsFactors = FALSE),
        species = "hsapiens",
        threshold = 0.05,
        sources = c("GO:BP"),
        domain_scope = "annotated",
        custom_bg = character(),
        exclude_iea = FALSE,
        max_retries = 1,
        wait_time = 0,
        protein_id_column = "protein_id"
      )
    )
  )

  pathway_dir <- tempfile("general-enrichment-plots-")
  dir.create(pathway_dir, recursive = TRUE)
  withr::defer(unlink(pathway_dir, recursive = TRUE, force = TRUE))

  local_mocked_bindings(
    ggplotly = function(plot, tooltip = "text") list(kind = "plotly", tooltip = tooltip),
    .package = "plotly"
  )

  enrichment_result <- list(
    meta = list(query_metadata = list(user_threshold = 0.05)),
    result = data.frame(
      source = c("GO:BP", "KEGG", "REAC"),
      p_value = c(0.001, 0.01, 0.02),
      term_name = c("Term 1", "Term 2", "Term 3"),
      term_id = c("GO:1", "KEGG:1", "REAC:1"),
      intersection_size = c(2, 3, 4),
      term_size = c(5, 6, 7),
      parents = I(list(c("PARENT1"), c("PARENT2"), c("PARENT3"))),
      stringsAsFactors = FALSE
    )
  )

  plots <- generate_enrichment_plots(
    enrichment_result = enrichment_result,
    contrast = "contrast_a",
    direction = "up",
    pathway_dir = pathway_dir
  )
  expect_s3_class(plots$static, "ggplot")
  expect_equal(plots$interactive$kind, "plotly")
  expect_true(file.exists(file.path(pathway_dir, "contrast_a_up_enrichment_results.tsv")))

  empty_plots <- generate_enrichment_plots(
    enrichment_result = list(meta = list(query_metadata = list(user_threshold = 0.05)), result = data.frame()),
    contrast = "contrast_a",
    direction = "down",
    pathway_dir = pathway_dir
  )
  expect_null(empty_plots$static)
  expect_null(empty_plots$interactive)

  summary_df <- summarize_enrichment(enrichment_result)
  expect_equal(summary_df$total[[1]], 3L)
  expect_equal(summary_df$GO_BP[[1]], 1L)
  expect_equal(summary_df$KEGG[[1]], 1L)
  expect_equal(summary_df$REAC[[1]], 1L)

  mock_results <- methods::new(
    "mockEnrichmentResultsCarrier",
    enrichment_data = list(contrast_a = list(up = data.frame(term = "A"), down = data.frame(term = "B"))),
    enrichment_plotly = list(contrast_a = list(up = list(kind = "plotly-up"), down = list(kind = "plotly-down"))),
    enrichment_summaries = list(
      contrast_a = list(
        up = data.frame(total = 3, GO_BP = 1, GO_CC = 0, GO_MF = 1, KEGG = 1, REAC = 0),
        down = data.frame(total = 1, GO_BP = 0, GO_CC = 1, GO_MF = 0, KEGG = 0, REAC = 0)
      )
    )
  )

  expect_equal(getEnrichmentResult(mock_results, "contrast_a", "up")$term[[1]], "A")
  expect_equal(getEnrichmentPlotly(mock_results, "contrast_a", "down")$kind, "plotly-down")
  expect_equal(getEnrichmentSummary(mock_results)$up_total[[1]], 3)
  expect_error(getEnrichmentResult(mock_results, "missing", "up"), "Contrast not found")
  expect_error(getEnrichmentPlotly(mock_results, "contrast_a", "sideways"), "Direction must be 'up' or 'down'")
})

if (!methods::isClass("mockClusterProfilerEnrichResultShared")) {
  methods::setClass(
    "mockClusterProfilerEnrichResultShared",
    slots = c(result = "data.frame")
  )
}

buildSharedDaResultsForEnrichment <- function(contrast_names, contrast_tables) {
  methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrast = contrast_names),
    da_data = stats::setNames(contrast_tables, contrast_names),
    design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE)
  )
}

test_that("processEnrichments validates exclude_iea before doing any work", {
  da_results <- buildSharedDaResultsForEnrichment(
    contrast_names = "internal_contrast",
    contrast_tables = list(data.frame(
      `Protein.IDs` = "P1:iso1",
      fdr_qvalue = 0.01,
      log2FC = 1.5,
      stringsAsFactors = FALSE
    ))
  )

  expect_error(
    processEnrichments(
      da_results = da_results,
      taxon_id = "9606",
      exclude_iea = NULL,
      pathway_dir = tempdir(),
      protein_id_column = "Protein.IDs"
    ),
    "exclude_iea must be explicitly set to TRUE or FALSE"
  )

  expect_error(
    processEnrichments(
      da_results = da_results,
      taxon_id = "9606",
      exclude_iea = "nope",
      pathway_dir = tempdir(),
      protein_id_column = "Protein.IDs"
    ),
    "exclude_iea must be a logical value"
  )
})

test_that("processEnrichments supported-organism branch preserves contrast remapping and artifact routing", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$perform <- list()
  captured$saved_html <- character()
  captured$saved_png <- character()

  localNamespaceBindings(
    package_ns,
    list(
      perform_enrichment = function(data_subset,
                                    species,
                                    threshold,
                                    sources,
                                    domain_scope,
                                    custom_bg,
                                    exclude_iea,
                                    protein_id_column,
                                    correction_method) {
        captured$perform[[length(captured$perform) + 1L]] <- list(
          ids = data_subset[[protein_id_column]],
          species = species,
          threshold = threshold,
          custom_bg = custom_bg,
          exclude_iea = exclude_iea,
          correction_method = correction_method
        )

        list(
          result = data.frame(
            source = c("GO:BP", "KEGG"),
            stringsAsFactors = FALSE
          )
        )
      },
      generate_enrichment_plots = function(enrichment_result, contrast, direction, pathway_dir) {
        list(
          static = ggplot2::ggplot(
            data.frame(x = 1, y = seq_len(nrow(enrichment_result$result))),
            ggplot2::aes(x, y)
          ) + ggplot2::geom_point(),
          interactive = list(kind = paste0(contrast, "::", direction), path = pathway_dir)
        )
      },
      summarize_enrichment = function(enrichment_result) {
        data.frame(
          total = nrow(enrichment_result$result),
          GO_BP = sum(enrichment_result$result$source == "GO:BP"),
          GO_CC = 0,
          GO_MF = 0,
          KEGG = sum(enrichment_result$result$source == "KEGG"),
          REAC = sum(enrichment_result$result$source == "REAC"),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  local_mocked_bindings(
    saveWidget = function(widget, file, selfcontained = TRUE) {
      captured$saved_html <- c(captured$saved_html, file)
      invisible(file)
    },
    .package = "htmlwidgets"
  )
  local_mocked_bindings(
    ggsave = function(filename, plot, width, height, dpi, bg = NULL, ...) {
      captured$saved_png <- c(captured$saved_png, filename)
      invisible(filename)
    },
    .package = "ggplot2"
  )

  pathway_dir <- tempfile("general-enrichment-supported-")
  dir.create(pathway_dir, recursive = TRUE)
  withr::defer(unlink(pathway_dir, recursive = TRUE, force = TRUE))

  da_results <- buildSharedDaResultsForEnrichment(
    contrast_names = c("internal.alpha", "internal.beta"),
    contrast_tables = list(
      data.frame(
        `Protein.IDs` = c("P1:iso1", "P2:iso2", "P3:iso3"),
        fdr_qvalue = c(0.01, 0.02, 0.20),
        log2FC = c(2.5, -1.4, 0.1),
        stringsAsFactors = FALSE
      ),
      data.frame(
        `Protein.IDs` = c("P4:iso1", "P5:iso2"),
        fdr_qvalue = c(0.01, 0.02),
        log2FC = c(1.2, 1.8),
        stringsAsFactors = FALSE
      )
    )
  )

  enrichment_results <- processEnrichments(
    da_results = da_results,
    taxon_id = "9606",
    q_cutoff = 0.05,
    pathway_dir = pathway_dir,
    exclude_iea = FALSE,
    protein_id_column = "Protein.IDs",
    contrast_names = c("short_alpha", "short_beta"),
    correction_method = "bonferroni"
  )

  expect_equal(names(enrichment_results@enrichment_data), c("internal.alpha", "internal.beta"))
  expect_equal(names(enrichment_results@enrichment_plotly), c("internal.alpha", "internal.beta"))
  expect_equal(captured$perform[[1]]$ids, "P1")
  expect_equal(captured$perform[[2]]$ids, "P2")
  expect_equal(captured$perform[[1]]$custom_bg, c("P1:iso1", "P2:iso2", "P3:iso3"))
  expect_equal(captured$perform[[1]]$species, "hsapiens")
  expect_false(captured$perform[[1]]$exclude_iea)
  expect_equal(captured$perform[[1]]$correction_method, "bonferroni")
  expect_equal(getEnrichmentResult(enrichment_results, "internal.alpha", "up")$result$source[[1]], "GO:BP")
  expect_equal(getEnrichmentPlotly(enrichment_results, "internal.alpha", "down")$kind, "internal.alpha::down")
  expect_equal(sort(getEnrichmentSummary(enrichment_results)$contrast), c("internal.alpha", "internal.beta"))
  expect_equal(getEnrichmentSummary(enrichment_results)$up_total, c(2, 2))
  expect_true(any(grepl("short_alpha_up_enrichment_plot.html$", captured$saved_html)))
  expect_true(any(grepl("short_alpha_down_enrichment_plot.png$", captured$saved_png)))
  expect_true(any(grepl("short_beta_up_enrichment_plot.html$", captured$saved_html)))
  expect_false(any(grepl("short_beta_down_enrichment_plot.html$", captured$saved_html)))
})

test_that("processEnrichments unsupported-organism branch preserves GO annotation enrichment and plot assembly", {
  localFakeClusterProfilerGeneral()
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$save_widget <- character()
  captured$save_png <- character()
  captured$plotly <- character()

  local_mocked_bindings(
    enricher = function(gene, universe, TERM2GENE, TERM2NAME, pvalueCutoff, pAdjustMethod, ...) {
      methods::new(
        "mockClusterProfilerEnrichResultShared",
        result = data.frame(
          ID = if (gene[[1]] == "P1") "GO:1001" else "GO:2001",
          Description = if (gene[[1]] == "P1") "Cell cycle term" else "Signal term",
          qvalue = 0.01,
          Count = length(gene),
          GeneRatio = sprintf("%d/%d", length(gene), length(gene)),
          BgRatio = sprintf("%d/%d", length(gene), length(universe)),
          stringsAsFactors = FALSE
        )
      )
    },
    .package = "clusterProfiler"
  )
  local_mocked_bindings(
    ggplotly = function(plot, tooltip = "text") {
      captured$plotly <- c(captured$plotly, tooltip)
      list(kind = "plotly", tooltip = tooltip)
    },
    .package = "plotly"
  )
  local_mocked_bindings(
    saveWidget = function(widget, file, selfcontained = TRUE) {
      captured$save_widget <- c(captured$save_widget, file)
      invisible(file)
    },
    .package = "htmlwidgets"
  )
  local_mocked_bindings(
    ggsave = function(filename, plot, width, height, dpi, bg = NULL, ...) {
      captured$save_png <- c(captured$save_png, filename)
      invisible(filename)
    },
    .package = "ggplot2"
  )

  pathway_dir <- tempfile("general-enrichment-unsupported-")
  dir.create(pathway_dir, recursive = TRUE)
  withr::defer(unlink(pathway_dir, recursive = TRUE, force = TRUE))

  da_results <- buildSharedDaResultsForEnrichment(
    contrast_names = "internal.gamma",
    contrast_tables = list(data.frame(
      `Protein.IDs` = c("P1", "P2", "P3"),
      fdr_qvalue = c(0.01, 0.02, 0.5),
      log2FC = c(2.3, -1.8, 0.1),
      stringsAsFactors = FALSE
    ))
  )
  go_annotations <- data.frame(
    Entry = c("P1", "P2"),
    go_id_go_biological_process = c("GO:1001", "GO:1001"),
    go_term_go_biological_process = c("Cell cycle term", "Cell cycle term"),
    go_id_go_molecular_function = c("GO:2001", "GO:2001"),
    go_term_go_molecular_function = c("Signal term", "Signal term"),
    go_id_go_cellular_compartment = c("GO:3001", "GO:3001"),
    go_term_go_cellular_compartment = c("Compartment term", "Compartment term"),
    stringsAsFactors = FALSE
  )

  enrichment_results <- processEnrichments(
    da_results = da_results,
    taxon_id = "999999",
    q_cutoff = 0.05,
    pathway_dir = pathway_dir,
    go_annotations = go_annotations,
    exclude_iea = FALSE,
    protein_id_column = "Protein.IDs",
    contrast_names = "short_gamma"
  )

  expect_equal(names(enrichment_results@enrichment_data), "short_gamma")
  expect_s3_class(enrichment_results@enrichment_plots$short_gamma$up, "ggplot")
  expect_equal(getEnrichmentPlotly(enrichment_results, "short_gamma", "down")$kind, "plotly")
  expect_s4_class(enrichment_results, "EnrichmentResults")
  expect_true(any(grepl("short_gamma_up_enrichment_plot.html$", captured$save_widget)))
  expect_true(any(grepl("short_gamma_down_enrichment_plot.png$", captured$save_png)))
  expect_equal(captured$plotly, c("text", "text"))
})

test_that("Revigo enrichment helpers preserve filtering, clustering, and file outputs", {
  package_ns <- asNamespace("MultiScholaR")
  temp_dir <- tempfile("general-enrichment-helpers-")
  dir.create(temp_dir, recursive = TRUE)
  withr::defer(unlink(temp_dir, recursive = TRUE, force = TRUE))

  stringi_env <- environment(queryRevigo)
  while (!exists("stri_replace_all_fixed", envir = stringi_env, inherits = FALSE) &&
         environmentIsLocked(stringi_env) &&
         !identical(parent.env(stringi_env), emptyenv())) {
    stringi_env <- parent.env(stringi_env)
  }

  localNamespaceBinding(
    stringi_env,
    "stri_replace_all_fixed",
    function(str, pattern, replacement) {
      gsub(pattern, replacement, str, fixed = TRUE)
    }
  )
  localNamespaceBinding(
    stringi_env,
    "discard",
    function(.x, .p) {
      predicate <- purrr::as_mapper(.p)
      .x[!vapply(.x, predicate, logical(1))]
    }
  )
  localNamespaceBinding(
    stringi_env,
    "html_nodes",
    function(x, css) {
      rvest::html_nodes(x, css)
    }
  )
  localNamespaceBinding(
    stringi_env,
    "html_table",
    function(x, ...) {
      rvest::html_table(x, ...)
    }
  )

  post_capture <- new.env(parent = emptyenv())
  local_mocked_bindings(
    POST = function(url, body, encode) {
      post_capture$url <- url
      post_capture$body <- body
      post_capture$encode <- encode
      list(
        text = paste0(
          "<html><body><table>",
          "<tr><th>Term ID</th><th>Name</th><th>Eliminated</th></tr>",
          "<tr><td>GO:1</td><td>Term 1</td><td>False</td></tr>",
          "<tr><td>GO:2</td><td>Term 2</td><td>True</td></tr>",
          "</table></body></html>"
        )
      )
    },
    content = function(x, encoding = NULL, ...) x$text,
    .package = "httr"
  )

  revigo_file <- file.path(temp_dir, "revigo.html")
  revigo_tbl <- queryRevigo(
    input_list = c("GO:1", "GO:2"),
    cutoff = 0.4,
    speciesTaxon = 9606,
    temp_file = revigo_file
  )

  expect_identical(post_capture$encode, "form")
  expect_match(post_capture$body$goList, "GO:1\nGO:2", fixed = TRUE)
  expect_identical(revigo_tbl$`Term ID`, c("GO:1", "GO:2"))
  expect_false(file.exists(revigo_file))

  populated_file <- file.path(temp_dir, "enrichment.tsv")
  empty_file <- file.path(temp_dir, "empty.tsv")
  readr::write_tsv(
    tibble::tibble(
      ID = "GO:1",
      names_of_genes_list = "positive_list",
      min_gene_set_size = 2L,
      max_gene_set_size = 20L,
      pvalue = 0.01,
      qvalue = 0.02,
      p.adjust = 0.02,
      term = "Term 1"
    ),
    populated_file
  )
  readr::write_tsv(
    tibble::tibble(
      ID = character(),
      names_of_genes_list = character(),
      min_gene_set_size = integer(),
      max_gene_set_size = integer()
    ),
    empty_file
  )

  enriched_from_files <- readEnrichmentResultFiles(
    table_of_files = tibble::tibble(
      file_name = c(populated_file, empty_file),
      comparison = c("contrast_a", "contrast_b"),
      Analysis_Type = c("Core", "Core")
    ),
    file_names_column = file_name,
    go_type = "KEGG"
  )

  expect_identical(enriched_from_files$annotation_id, "GO:1")
  expect_identical(enriched_from_files$gene_set, "positive_list")
  expect_identical(enriched_from_files$go_type, "KEGG")
  expect_identical(enriched_from_files$comparison, "contrast_a")

  localNamespaceBinding(
    package_ns,
    "queryRevigo",
    function(input_list, cutoff = 0.7, speciesTaxon = 9606, temp_file = NA) {
      data.frame(
        `Term ID` = input_list,
        Name = paste("Term", seq_along(input_list)),
        Eliminated = c("False", rep("True", max(length(input_list) - 1L, 0L))),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  filtered <- filterResultsWithRevigo(
    enriched_results_tbl = data.frame(
      Analysis_Type = c("Core", "Core"),
      comparison = c("contrast_a", "contrast_a"),
      gene_set = c("positive_list", "positive_list"),
      go_type = c("BP", "BP"),
      annotation_id = c("GO:1", "GO:2"),
      qvalue = c(0.01, 0.02),
      stringsAsFactors = FALSE
    ),
    added_columns = "Analysis_Type",
    revigo_cutoff = 0.6
  )
  expect_identical(filtered$annotation_id, "GO:1")

  scholar_filtered <- filterResultsWithRevigoScholar(
    enriched_results_tbl = data.frame(
      Analysis_Type = c("Core", "Core"),
      go_type = c("BP", "BP"),
      ID = c("GO:1", "GO:2"),
      score = c(1, 2),
      stringsAsFactors = FALSE
    ),
    added_columns = "Analysis_Type"
  )
  expect_identical(scholar_filtered$go_id, "GO:1")

  localNamespaceBinding(
    package_ns,
    "queryRevigo",
    function(input_list, cutoff = 0.7, speciesTaxon = 9606, temp_file = NA) {
      data.frame()
    }
  )
  expect_warning(
    filterResultsWithRevigo(
      enriched_results_tbl = data.frame(
        Analysis_Type = "Core",
        comparison = "contrast_a",
        gene_set = "positive_list",
        go_type = "BP",
        annotation_id = "GO:1",
        stringsAsFactors = FALSE
      ),
      added_columns = "Analysis_Type"
    ),
    "Revigo summarizatio"
  )

  heatmap_tbl <- data.frame(
    Analysis_Type = c("Core", "Core", "Core", "Aux"),
    comparison = c("A", "A", "B", "B"),
    annotation_id = c("GO:1", "GO:1", "GO:2", "GO:3"),
    term = c("Term 1", "Term 1", "Term 2", "Term 3"),
    p.adjust = c(0.01, 0.03, 0.04, 0.05),
    pvalue = c(0.01, 0.03, 0.04, 0.05),
    qvalue = c(0.01, 0.03, 0.04, 0.05),
    gene_set = c("positive_list", "negative_list", "shared", "positive_list"),
    go_type = c("BP", "BP", "BP", "MF"),
    min_set_size = c(2L, 2L, 2L, 2L),
    max_set_size = c(20L, 20L, 20L, 20L),
    gene_symbol = c(strrep("GENE", 2000), "A", "B", "C"),
    stringsAsFactors = FALSE
  )

  clustered <- clusterPathways(
    input_table = heatmap_tbl,
    added_columns = "Analysis_Type",
    remove_duplicted_entries = "merge"
  )
  expect_true("shared" %in% clustered$gene_set)

  heatmap_plot <- getEnrichmentHeatmap(
    input_table = clustered,
    x_axis = comparison,
    input_go_type = "BP",
    input_plot_title = "BP terms",
    xaxis_levels = c("A", "B")
  )
  expect_s3_class(heatmap_plot, "ggplot")

  heatmap_list <- drawListOfFunctionalEnrichmentHeatmaps(
    enriched_results_tbl = heatmap_tbl,
    added_columns = "Analysis_Type",
    set_size_min = 2L,
    set_size_max = 20L,
    x_axis = Analysis_Type,
    analysis_column = Analysis_Type,
    remove_duplicted_entries = "merge",
    xaxis_levels = c("A", "B")
  )
  expect_equal(sort(names(heatmap_list)), c("BP", "MF"))
  expect_true(all(vapply(heatmap_list, inherits, logical(1), what = "ggplot")))

  written_xlsx <- character()
  local_mocked_bindings(
    write_xlsx = function(x, path) {
      written_xlsx <<- c(written_xlsx, path)
      writeLines("xlsx", path)
      invisible(path)
    },
    .package = "writexl"
  )

  saveFilteredFunctionalEnrichmentTable(
    enriched_results_tbl = heatmap_tbl,
    set_size_min = 2L,
    set_size_max = 20L,
    results_dir = temp_dir,
    file_name = "functional-enrichment",
    list_of_columns_to_trim = "gene_symbol"
  )

  expect_true(file.exists(file.path(temp_dir, "functional-enrichment.tab")))
  expect_true(file.exists(file.path(temp_dir, "functional-enrichment_unfiltered.tab")))
  expect_equal(length(written_xlsx), 2L)
  expect_true(all(file.exists(written_xlsx)))
})

test_that("GO enrichment wrapper preserves merged enrichment output and gene symbol expansion", {
  localFakeClusterProfilerGeneral()
  options(
    multischolar.fake_clusterProfiler.enricher = function(gene, universe, TERM2GENE, TERM2NAME, ...) {
      methods::new(
        "enrichResult",
        result = data.frame(
          ID = TERM2NAME$go_id[[1]],
          Description = TERM2NAME$go_term[[1]],
          GeneRatio = sprintf("%d/%d", min(2L, length(gene)), length(gene)),
          BgRatio = sprintf("%d/%d", length(unique(TERM2GENE$uniprot_acc)), length(universe)),
          pvalue = 0.01,
          p.adjust = 0.02,
          qvalue = 0.02,
          geneID = paste(utils::head(intersect(gene, TERM2GENE$uniprot_acc), 2L), collapse = "/"),
          Count = length(unique(TERM2GENE$uniprot_acc)),
          stringsAsFactors = FALSE
        )
      )
    }
  )

  go_annotations <- data.frame(
    uniprot_acc = c("P1", "P2", "P3", "P1", "P2", "P3", "P1", "P2", "P3"),
    go_id = c("GO:BP1", "GO:BP1", "GO:BP1", "GO:MF1", "GO:MF1", "GO:MF1", "GO:CC1", "GO:CC1", "GO:CC1"),
    go_term = c("Cell cycle", "Cell cycle", "Cell cycle", "Binding", "Binding", "Binding", "Organelle", "Organelle", "Organelle"),
    go_type = c(
      "Biological Process",
      "Biological Process",
      "Biological Process",
      "Molecular Function",
      "Molecular Function",
      "Molecular Function",
      "Cellular Compartment",
      "Cellular Compartment",
      "Cellular Compartment"
    ),
    stringsAsFactors = FALSE
  )
  uniprot_data <- data.frame(
    Entry = c("P1", "P2", "P3"),
    gene_names = c("GeneA alpha", "GeneB beta", "GeneC gamma"),
    stringsAsFactors = FALSE
  )

  enrichment_results <- runOneGoEnrichmentInOutFunction(
    significant_proteins = c("P1", "P2"),
    background_proteins = c("P1", "P2", "P3"),
    go_annotations = go_annotations,
    uniprot_data = uniprot_data,
    p_val_thresh = 0.05,
    min_gene_set_size = 2L,
    max_gene_set_size = 4L,
    min_sig_gene_set_size = 1L
  )

  expect_equal(sort(enrichment_results$ID), c("GO:BP1", "GO:CC1", "GO:MF1"))
  expect_true(all(c("GeneA", "GeneB") %in% unlist(strsplit(enrichment_results$gene_names[[1]], ";"))))
  expect_true(all(enrichment_results$p.adjust <= 0.02))
  expect_equal(nrow(enrichment_results), 3L)
})

test_that("general enrichment fallback helpers preserve utility and no-Revigo branches", {
  expect_equal(parseNumList("1,2;3:4"), c(1L, 2L, 3L, 4L))
  expect_equal(parseNumList("5"), 5L)

  go_annot <- data.frame(
    uniprot_acc = c("P1", "P2", "P1", "P2"),
    go_id = c("GO:1", "GO:1", "GO:2", "GO:2"),
    go_type = c("BP", "BP", "MF", "MF"),
    stringsAsFactors = FALSE
  )
  background <- data.frame(uniprot_acc = c("P1", "P2"), stringsAsFactors = FALSE)

  one_go_default <- makeFunctionWithOverrides(
    oneGoEnrichment,
    list(
      enricher = function(...) {
        data.frame(
          ID = "GO:1",
          Description = "placeholder",
          GeneRatio = "2/2",
          BgRatio = "2/2",
          pvalue = 0.01,
          p.adjust = 0.02,
          qvalue = 0.02,
          geneID = "P1/P2",
          Count = 2,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  default_output <- one_go_default(
    go_annot = go_annot,
    background_list = background,
    go_aspect = NA,
    query_list = c("P1", "P2"),
    id_to_annotation_dictionary = c(`GO:1` = "Term A", `GO:2` = "Term B"),
    annotation_id = go_id,
    protein_id = uniprot_acc,
    aspect_column = go_type,
    p_val_thresh = 0.05,
    min_gene_set_size = 2,
    max_gene_set_size = 4,
    get_cluster_profiler_object = FALSE
  )
  expect_equal(default_output$term[[1]], "Term A")
  expect_false("cluster_profiler_object" %in% names(default_output))

  one_go_null <- makeFunctionWithOverrides(
    oneGoEnrichment,
    list(enricher = function(...) NULL)
  )
  expect_true(is.na(
    one_go_null(
      go_annot = go_annot,
      background_list = background,
      go_aspect = NA,
      query_list = c("P1", "P2"),
      id_to_annotation_dictionary = c(`GO:1` = "Term A"),
      annotation_id = go_id,
      protein_id = uniprot_acc,
      aspect_column = go_type,
      p_val_thresh = 0.05,
      min_gene_set_size = 2,
      max_gene_set_size = 4
    )
  ))

  localFakeClusterProfilerGeneral()
  options(
    multischolar.fake_clusterProfiler.enricher = function(...) {
      methods::new("enrichResult", result = data.frame())
    }
  )
  expect_null(
    runOneGoEnrichmentInOutFunction(
      significant_proteins = c("P1", "P2"),
      background_proteins = c("P1", "P2"),
      go_annotations = data.frame(
        uniprot_acc = c("P1", "P2"),
        go_id = c("GO:1", "GO:1"),
        go_term = c("Term A", "Term A"),
        go_type = c("Biological Process", "Biological Process"),
        stringsAsFactors = FALSE
      ),
      uniprot_data = data.frame(
        Entry = c("P1", "P2"),
        gene_names = c("GeneA", "GeneB"),
        stringsAsFactors = FALSE
      ),
      p_val_thresh = 0.05,
      min_gene_set_size = 2L,
      max_gene_set_size = 4L,
      min_sig_gene_set_size = 3L
    )
  )

  package_ns <- asNamespace("MultiScholaR")
  localNamespaceBinding(
    package_ns,
    "queryRevigo",
    function(input_list, cutoff = 0.7, speciesTaxon = 9606, temp_file = NA) {
      data.frame(
        `Term ID` = input_list,
        Name = paste("Term", seq_along(input_list)),
        Eliminated = "False",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )

  enrichment_tbl <- data.frame(
    Analysis_Type = "Core",
    comparison = "contrast_a",
    gene_set = "positive_list",
    go_type = "BP",
    annotation_id = "GO:1",
    stringsAsFactors = FALSE
  )
  expect_equal(
    filterResultsWithRevigo(
      enriched_results_tbl = enrichment_tbl,
      added_columns = "Analysis_Type",
      is_run_revigo = FALSE
    ),
    enrichment_tbl
  )

  scholar_tbl <- data.frame(
    Analysis_Type = "Core",
    go_type = "BP",
    ID = "GO:1",
    stringsAsFactors = FALSE
  )
  expect_equal(
    filterResultsWithRevigoScholar(
      enriched_results_tbl = scholar_tbl,
      added_columns = "Analysis_Type",
      is_run_revigo = FALSE
    ),
    scholar_tbl
  )
})

test_that("general enrichment scholar plotting helpers preserve export and edge branches", {
  heatmap_tbl <- data.frame(
    Analysis_Type = c("Core", "Core", "Core", "Aux"),
    comparison = c("A", "A", "B", "B"),
    annotation_id = c("GO:1", "GO:1", "GO:2", "GO:3"),
    term = c("Term 1", "Term 1", "Term 2", "Term 3"),
    p.adjust = c(0.01, 0.03, 0.04, 0.05),
    pvalue = c(0.01, 0.03, 0.04, 0.05),
    qvalue = c(0.01, 0.03, 0.04, 0.05),
    gene_set = c("positive_list", "negative_list", "shared", "positive_list"),
    go_type = c("BP", "BP", "BP", "MF"),
    min_set_size = c(2L, 2L, 2L, 2L),
    max_set_size = c(20L, 20L, 20L, 20L),
    gene_symbol = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    evaluateBestMinMaxGeneSetSize(heatmap_tbl, "Analysis_Type"),
    "ggplot"
  )

  expect_error(
    getEnrichmentHeatmap(
      input_table = heatmap_tbl,
      x_axis = comparison,
      input_go_type = "BP",
      input_plot_title = "BP terms",
      xaxis_levels = c("A", "missing-level")
    )
  )

  expect_s3_class(
    getEnrichmentHeatmap(
      input_table = heatmap_tbl,
      x_axis = comparison,
      input_go_type = NA,
      input_plot_title = "All terms",
      xaxis_levels = NA
    ),
    "ggplot"
  )

  scholar_fun <- makeFunctionWithOverrides(
    drawListOfFunctionalEnrichmentHeatmapsScholar,
    list(
      clusterPathways = function(input_table, added_columns, remove_duplicted_entries = TRUE) {
        input_table
      },
      getEnrichmentHeatmap = function(...) {
        ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
          ggplot2::geom_point()
      }
    )
  )

  scholar_heatmaps <- scholar_fun(
    enriched_results_tbl = heatmap_tbl,
    added_columns = "Analysis_Type",
    x_axis = comparison,
    analysis_column = Analysis_Type,
    facet_by_column = "Analysis_Type",
    remove_duplicted_entries = "delete",
    xaxis_levels = c("A", "B")
  )
  expect_equal(sort(names(scholar_heatmaps)), c("BP", "MF"))
  expect_true(all(vapply(scholar_heatmaps, inherits, logical(1), what = "ggplot")))

  saved_calls <- new.env(parent = emptyenv())
  saved_calls$items <- list()
  save_heatmaps <- makeFunctionWithOverrides(
    saveListOfFunctionalEnrichmentHeatmaps,
    list(
      savePlot = function(plot, base_dir, plot_name, width, height, limitsize) {
        saved_calls$items[[length(saved_calls$items) + 1L]] <<- list(
          base_dir = base_dir,
          plot_name = plot_name,
          width = width,
          height = height,
          limitsize = limitsize
        )
        invisible(plot_name)
      }
    )
  )

  plot_list <- list(
    BP = ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(),
    MF = ggplot2::ggplot(data.frame(x = 1, y = 2), ggplot2::aes(x, y)) + ggplot2::geom_point()
  )

  save_heatmaps(
    list_of_heatmaps = plot_list,
    results_dir = tempdir(),
    file_name = "heatmap",
    plot_width = 7,
    plot_height = 8
  )
  expect_length(saved_calls$items, 2L)
  expect_error(
    save_heatmaps(
      list_of_heatmaps = plot_list,
      results_dir = tempdir(),
      file_name = "heatmap",
      plot_width = c(7, 8, 9),
      plot_height = 8
    ),
    "Length of plot_width"
  )

  expect_s3_class(
    enrichedPathwayBarPlot(
      input_table = heatmap_tbl,
      input_go_type = NA,
      remove_duplicted_entries = TRUE,
      added_columns = "Analysis_Type"
    ),
    "ggplot"
  )
  expect_error(
    enrichedPathwayBarPlot(
      input_table = heatmap_tbl,
      input_go_type = "missing",
      remove_duplicted_entries = TRUE,
      added_columns = "Analysis_Type"
    ),
    "is not in the input_table"
  )

  saved_files <- character()
  go_bar_fun <- makeFunctionWithOverrides(
    enrichedGoTermBarPlot,
    list(
      filtered_enrich_revigo = heatmap_tbl,
      ggsave = function(plot, filename, width, height) {
        saved_files <<- c(saved_files, filename)
        invisible(filename)
      }
    )
  )
  go_bar_fun(
    input_table = heatmap_tbl,
    output_dir = tempdir(),
    analysis_type = "GO",
    file_suffix = c("png", "pdf"),
    width = 6,
    height = 4
  )
  expect_true(length(saved_files) >= 4L)
})
