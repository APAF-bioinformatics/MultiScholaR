library(testthat)
library(shiny)

test_that("MCI-009.1 DA session reload accepts workflow variants and refuses invalid sessions", {
  valid_variants <- list(
    DIA = list(workflow_type = "DIA", use_limpa = FALSE, report_template = "DIANN_report.rmd"),
    DIA_limpa = list(workflow_type = "DIA", use_limpa = TRUE, report_template = "DIANN_limpa_report.rmd"),
    TMT = list(workflow_type = "TMT", use_limpa = FALSE, report_template = "TMT_report.rmd"),
    LFQ_MaxQuant = list(workflow_type = "LFQ_MaxQuant", use_limpa = FALSE, report_template = "LFQ_report.rmd"),
    LFQ_FragPipe = list(workflow_type = "LFQ_FragPipe", use_limpa = FALSE, report_template = "LFQ_report.rmd")
  )

  module_ci_prot_da_with_globals(code = {
    for (variant_name in names(valid_variants)) {
      variant <- valid_variants[[variant_name]]
      paths <- module_ci_prot_da_paths()
      payload <- do.call(module_ci_prot_da_session_payload, variant)
      module_ci_prot_da_write_session(paths$source_dir, payload)

      da_data <- shiny::reactiveValues(
        current_s4_object = NULL,
        contrasts_available = NULL,
        formula_from_s4 = NULL
      )
      workflow_data <- module_ci_prot_da_workflow(payload$current_s4_object)
      captured <- module_ci_prot_da_shiny_mocks()

      testServer(
        module_ci_prot_da_load_module(da_data, workflow_data, paths),
        {
          session$setInputs(load_filtered_session = 1)
          session$flushReact()

          expect_s4_class(da_data$current_s4_object, "ProteinQuantitativeData")
          expect_identical(
            da_data$current_s4_object@args$globalParameters$workflow_type,
            variant$workflow_type
          )
          expect_identical(
            da_data$current_s4_object@args$globalParameters$use_limpa,
            variant$use_limpa
          )
          expect_identical(
            da_data$current_s4_object@args$globalParameters$report_template,
            variant$report_template
          )
          expect_identical(da_data$contrasts_available, "groupB-groupA")
          expect_identical(da_data$formula_from_s4, "~ 0 + group")
          expect_identical(captured$text_updates$formula_string, "~ 0 + group")
          expect_identical(workflow_data$state_manager$current_state, "correlation_filtered")
          expect_identical(workflow_data$tab_status$normalization, "complete")
          expect_identical(workflow_data$tab_status$differential_expression, "pending")
          expect_s3_class(workflow_data$state_update_trigger, "POSIXct")
          expect_equal(nrow(workflow_data$uniprot_dat_cln), 6L)
          expect_true(workflow_data$accession_cleanup_results$cleanup_applied)
          expect_identical(workflow_data$config_list$globalParameters$workflow_type, variant$workflow_type)
          expect_true(any(vapply(captured$notifications, function(notification) {
            identical(notification$type, "message")
          }, logical(1))))
        }
      )
    }

    invalid_cases <- list(
      missing_session = function(paths) {
        dir.create(paths$source_dir, recursive = TRUE)
        "No exported session data found"
      },
      stale_source_dir = function(paths) {
        paths$source_dir <- file.path(paths$source_dir, "missing")
        "Could not find the source directory"
      },
      malformed_rds = function(paths) {
        dir.create(paths$source_dir, recursive = TRUE)
        saveRDS(list(not_a_session = TRUE), file.path(paths$source_dir, "filtered_session_data_latest.rds"))
        "ProteinQuantitativeData"
      }
    )

    for (case_name in names(invalid_cases)) {
      paths <- module_ci_prot_da_paths()
      expected_error <- invalid_cases[[case_name]](paths)
      da_data <- shiny::reactiveValues(current_s4_object = NULL, contrasts_available = NULL)
      workflow_data <- module_ci_prot_da_workflow()
      captured <- module_ci_prot_da_shiny_mocks()

      testServer(
        module_ci_prot_da_load_module(da_data, workflow_data, paths),
        {
          session$setInputs(load_filtered_session = 1)
          session$flushReact()

          expect_null(da_data$current_s4_object)
          expect_null(da_data$contrasts_available)
          expect_identical(workflow_data$tab_status$differential_expression, "locked")
          error_messages <- vapply(captured$notifications, `[[`, character(1), "message")
          expect_true(any(grepl(expected_error, error_messages, fixed = TRUE)), info = case_name)
        }
      )
    }
  })
})

test_that("MCI-009.2 DA formula and contrast validation covers valid, invalid, duplicate, and reversed cases", {
  two_group_object <- module_ci_prot_da_object()
  two_group <- validateProtDaModelAndContrasts(
    two_group_object,
    "~ 0 + group",
    module_ci_prot_da_contrasts("two_group")
  )
  expect_identical(two_group$contrast_strings, "B_vs_A=groupB-groupA")
  expect_true(all(c("groupA", "groupB") %in% colnames(two_group$design_model)))

  reversed <- validateProtDaModelAndContrasts(
    two_group_object,
    "~ 0 + group",
    module_ci_prot_da_contrasts("reversed")
  )
  expect_identical(reversed$contrast_strings, "A_vs_B=groupA-groupB")

  raw_only <- validateProtDaModelAndContrasts(
    two_group_object,
    "~ 0 + group",
    module_ci_prot_da_contrasts("raw_only")
  )
  expect_identical(raw_only$contrast_strings, "B_vs_A=groupB-groupA")

  batch_adjusted <- validateProtDaModelAndContrasts(
    two_group_object,
    "~ 0 + group + batch",
    module_ci_prot_da_contrasts("batch_adjusted")
  )
  expect_true(any(grepl("^batch", colnames(batch_adjusted$design_model))))

  multi_design <- module_ci_prot_da_design(
    samples = c("A1", "A2", "B1", "B2", "C1", "C2"),
    groups = c("A", "A", "B", "B", "C", "C"),
    batches = c("Batch1", "Batch2", "Batch1", "Batch2", "Batch1", "Batch2")
  )
  multi_group <- validateProtDaModelAndContrasts(
    module_ci_prot_da_object(design = multi_design),
    "~ 0 + group",
    module_ci_prot_da_contrasts("multi_group")
  )
  expect_identical(
    multi_group$contrast_strings,
    c("B_vs_A=groupB-groupA", "C_vs_A=groupC-groupA")
  )

  expect_error(
    validateProtDaModelAndContrasts(two_group_object, "~ 0 + missing_term", module_ci_prot_da_contrasts("two_group")),
    "Invalid DA formula",
    fixed = TRUE
  )
  expect_error(
    validateProtDaModelAndContrasts(two_group_object, "~ 0 + group", module_ci_prot_da_contrasts("empty")),
    "empty contrast",
    fixed = TRUE
  )
  expect_error(
    validateProtDaModelAndContrasts(two_group_object, "~ 0 + group", module_ci_prot_da_contrasts("duplicate")),
    "duplicate contrast definitions",
    fixed = TRUE
  )
  expect_error(
    validateProtDaModelAndContrasts(two_group_object, "~ 0 + group", module_ci_prot_da_contrasts("invalid_term")),
    "absent from the model matrix",
    fixed = TRUE
  )

  invalid_handler_cases <- list(
    invalid_formula = list(formula = "~ 0 + missing_term", contrasts = module_ci_prot_da_contrasts("two_group")),
    duplicate_contrast = list(formula = "~ 0 + group", contrasts = module_ci_prot_da_contrasts("duplicate")),
    empty_contrast = list(formula = "~ 0 + group", contrasts = module_ci_prot_da_contrasts("empty")),
    invalid_contrast_term = list(formula = "~ 0 + group", contrasts = module_ci_prot_da_contrasts("invalid_term"))
  )

  for (case_name in names(invalid_handler_cases)) {
    case <- invalid_handler_cases[[case_name]]
    object <- module_ci_prot_da_object()
    da_data <- shiny::reactiveValues(
      da_results_list = NULL,
      contrasts_available = case$contrasts$contrasts,
      analysis_complete = FALSE,
      current_s4_object = object
    )
    workflow_data <- module_ci_prot_da_workflow(object)
    paths <- module_ci_prot_da_paths()
    stats_called <- FALSE
    captured <- module_ci_prot_da_shiny_mocks()

    module_ci_prot_da_with_globals(list(contrasts_tbl = case$contrasts), {
      module_ci_prot_da_with_namespace_bindings(
        list(
          differentialAbundanceAnalysis = function(...) {
            stats_called <<- TRUE
            stop("statistical code should not run", call. = FALSE)
          },
          outputDaResultsAllContrasts = function(...) {
            stop("export should not run", call. = FALSE)
          },
          .capture_checkpoint = function(value, checkpoint_id, label) invisible(value)
        ),
        {
          testServer(
            module_ci_prot_da_run_module(da_data, workflow_data, paths),
            {
              session$setInputs(
                formula_string = case$formula,
                da_q_val_thresh = 0.05,
                treat_lfc_cutoff = 0
              )
              session$setInputs(run_da_analysis = 1)
              session$flushReact()

              expect_false(stats_called, info = case_name)
              expect_false(isTRUE(da_data$analysis_complete), info = case_name)
              expect_null(da_data$da_results_list, info = case_name)
              expect_true(any(vapply(captured$notifications, function(notification) {
                identical(notification$type, "error")
              }, logical(1))), info = case_name)
            }
          )
        }
      )
    })
  }
})

test_that("MCI-009.3 DA statistical output matrix preserves schemas across q-value edge cases", {
  object <- module_ci_prot_da_object()
  design <- object@design_matrix
  data <- as.matrix(object@protein_quant_table[, setdiff(names(object@protein_quant_table), object@protein_id_column)])
  rownames(data) <- object@protein_quant_table[[object@protein_id_column]]

  current_variant <- "standard"
  top_table_for_variant <- function(variant) {
    fixture <- module_ci_prot_da_result_table(variant)
    output <- data.frame(
      logFC = fixture$log2FC,
      AveExpr = seq_len(nrow(fixture)),
      t = seq_len(nrow(fixture)) / 2,
      P.Value = fixture$raw_pvalue,
      adj.P.Val = fixture$fdr_value_bh_adjustment,
      B = seq_len(nrow(fixture)) / 10,
      row.names = fixture$uniprot_acc
    )
    output
  }

  module_ci_prot_da_with_namespace_bindings(
    list(
      makeContrasts = function(contrasts, levels) {
        matrix(1, nrow = length(levels), ncol = length(contrasts), dimnames = list(levels, contrasts))
      },
      lmFit = function(object, design, block = NULL, correlation = NULL, ...) {
        list(object = object, design = design, block = block, correlation = correlation)
      },
      contrasts.fit = function(fit, contrasts) {
        fit$contrasts <- contrasts
        fit
      },
      eBayes = function(fit, trend = FALSE, robust = FALSE, ...) {
        fit$trend <- trend
        fit$robust <- robust
        fit
      },
      topTable = function(fit, coef, n = Inf, ...) top_table_for_variant(current_variant),
      qvalue = function(p) list(q = pmin(1, as.numeric(p) * 1.5))
    ),
    {
      variants <- c("standard", "no_significant", "all_significant", "tied_pvalues", "missing_values", "small_n")
      for (variant in variants) {
        current_variant <- variant
        result <- runTestsContrasts(
          data = data,
          contrast_strings = "B_vs_A=groupB-groupA",
          design_matrix = design,
          formula_string = "~ 0 + group",
          treat_lfc_cutoff = NA
        )
        table <- result$results[[1L]]

        expect_true(all(c("raw_pvalue", "fdr_qvalue", "fdr_value_bh_adjustment") %in% names(table)), info = variant)
        expect_equal(nrow(table), nrow(module_ci_prot_da_result_table(variant)), info = variant)
        expect_false(any(is.nan(table$fdr_qvalue)), info = variant)
        if (identical(variant, "small_n")) {
          expect_equal(table$fdr_qvalue, p.adjust(table$raw_pvalue, method = "BH"))
        }
        if (identical(variant, "missing_values")) {
          expect_true(any(is.na(table$fdr_qvalue)))
          expect_true(all(is.na(table$fdr_qvalue[!is.finite(table$raw_pvalue) | is.na(table$raw_pvalue)])))
        }

        da_long <- module_ci_prot_da_result_table(variant)
        module_ci_prot_da_assert_result_schema(da_long)
        counts <- countStatDaGenes(
          da_long,
          lfc_thresh = 0,
          q_val_thresh = 0.05,
          log_fc_column = log2FC,
          q_value_column = fdr_qvalue
        )
        expect_setequal(counts$status, c("Not significant", "Significant and Up", "Significant and Down"))
      }
    }
  )
})

test_that("MCI-009.4 DA render matrix covers volcano, heatmap, table, thresholds, and empty states", {
  captured <- module_ci_prot_da_render_mocks()

  table_output <- new.env(parent = emptyenv())
  table_data <- module_ci_prot_da_render_state(module_ci_prot_da_result_table("standard"))
  da_server_table_render_handler(
    input = list(
      table_contrast = "B_vs_A",
      table_significance = "up",
      da_q_val_thresh = 0.05,
      treat_lfc_cutoff = 1,
      table_max_rows = 1
    ),
    output = table_output,
    session = list(),
    da_data = table_data
  )
  expect_identical(table_output$da_results_table$data$uniprot_acc, "P1")
  expect_match(table_output$da_summary_stats, "Significant (q < 0.050): 4", fixed = TRUE)
  expect_match(table_output$da_summary_stats, "Up-regulated: 2", fixed = TRUE)
  expect_match(table_output$da_summary_stats, "Down-regulated: 1", fixed = TRUE)

  threshold_output <- new.env(parent = emptyenv())
  da_server_table_render_handler(
    input = list(
      table_contrast = "B_vs_A",
      table_significance = "significant",
      da_q_val_thresh = 0.01,
      treat_lfc_cutoff = 1,
      table_max_rows = 25
    ),
    output = threshold_output,
    session = list(),
    da_data = table_data
  )
  expect_identical(nrow(threshold_output$da_results_table$data), 0L)
  expect_match(threshold_output$da_summary_stats, "Significant (q < 0.010): 0", fixed = TRUE)

  empty_output <- new.env(parent = emptyenv())
  empty_data <- module_ci_prot_da_render_state(module_ci_prot_da_result_table("no_significant"))
  da_server_table_render_handler(
    input = list(
      table_contrast = "missing",
      table_significance = "significant",
      da_q_val_thresh = 0.05,
      treat_lfc_cutoff = 1,
      table_max_rows = 25
    ),
    output = empty_output,
    session = list(),
    da_data = empty_data
  )
  expect_identical(empty_output$da_results_table$data$Message, "No results available for selected contrast")
  expect_identical(empty_output$da_summary_stats, "No results available for selected contrast")

  heatmap_output <- new.env(parent = emptyenv())
  heatmap_data <- module_ci_prot_da_render_state(module_ci_prot_da_result_table("standard"))
  heatmap_args <- NULL
  module_ci_prot_da_with_namespace_bindings(
    list(
      .capture_checkpoint = function(payload, checkpoint_id, label) invisible(payload),
      generateProtDAHeatmap = function(...) {
        heatmap_args <<- list(...)
        list(
          plot = "heatmap-plot",
          row_clusters = c(P1 = 1, P2 = 2),
          col_clusters = c(S1 = 1, S2 = 1)
        )
      },
      save_heatmap_products = function(...) invisible(NULL)
    ),
    {
      da_server_heatmap_render_handler(
        input = list(
          heatmap_contrast = "B_vs_A",
          heatmap_top_n = 3,
          heatmap_cluster_method = "average",
          heatmap_distance_method = "euclidean",
          heatmap_clustering = "both",
          heatmap_scaling = "row",
          heatmap_color_scheme = "Viridis",
          heatmap_show_labels = TRUE,
          da_q_val_thresh = 0.05,
          heatmap_tree_cut_method = "k_clusters",
          heatmap_n_clusters = 2,
          heatmap_cut_height = 0.8,
          heatmap_min_cluster_size = 2,
          save_heatmap = NULL
        ),
        output = heatmap_output,
        session = list(),
        ns = function(id) id,
        da_data = heatmap_data,
        experiment_paths = list(publication_graphs_dir = tempdir())
      )
    }
  )
  expect_identical(heatmap_output$heatmap_plot, "heatmap-plot")
  expect_identical(heatmap_args$top_n_genes, 3)
  expect_identical(heatmap_args$selected_contrast, "B_vs_A")
  expect_match(heatmap_output$cluster_summary, "Total Clusters: 2", fixed = TRUE)

  empty_heatmap_output <- new.env(parent = emptyenv())
  module_ci_prot_da_with_namespace_bindings(
    list(
      .capture_checkpoint = function(payload, checkpoint_id, label) invisible(payload),
      generateProtDAHeatmap = function(...) NULL
    ),
    {
      expect_no_error(da_server_heatmap_render_handler(
        input = list(
          heatmap_contrast = "B_vs_A",
          heatmap_top_n = 20,
          heatmap_cluster_method = "average",
          heatmap_distance_method = "euclidean",
          heatmap_clustering = "row",
          heatmap_scaling = "row",
          heatmap_color_scheme = "Viridis",
          heatmap_show_labels = TRUE,
          da_q_val_thresh = 0.0001,
          heatmap_tree_cut_method = "k_clusters",
          heatmap_n_clusters = 2,
          heatmap_cut_height = 0.8,
          heatmap_min_cluster_size = 2,
          save_heatmap = NULL
        ),
        output = empty_heatmap_output,
        session = list(),
        ns = function(id) id,
        da_data = heatmap_data,
        experiment_paths = list(publication_graphs_dir = tempdir())
      ))
    }
  )

  volcano_output <- new.env(parent = emptyenv())
  volcano_data <- module_ci_prot_da_render_state(module_ci_prot_da_result_table("standard"))
  volcano_args <- list()
  module_ci_prot_da_with_namespace_bindings(
    list(
      .capture_checkpoint = function(payload, checkpoint_id, label) invisible(payload),
      generateProtDAVolcanoPlotGlimma = function(...) {
        volcano_args$interactive <<- list(...)
        "interactive-volcano"
      },
      generateProtDAVolcanoStatic = function(...) {
        volcano_args$static <<- list(...)
        "static-volcano"
      }
    ),
    {
      da_server_volcano_render_handler(
        input = list(
          volcano_contrast = "B_vs_A",
          da_q_val_thresh = 0.05,
          treat_lfc_cutoff = 1,
          volcano_show_labels = TRUE,
          volcano_label_top_n = 5
        ),
        output = volcano_output,
        session = list(),
        ns = function(id) id,
        da_data = volcano_data,
        experiment_paths = list(da_output_dir = tempdir())
      )
    }
  )
  expect_identical(volcano_output$volcano_glimma, "interactive-volcano")
  expect_identical(volcano_output$volcano_plot_static, "static-volcano")
  expect_identical(volcano_args$interactive$selected_contrast, "B_vs_A")
  expect_identical(volcano_args$static$lfc_threshold, 1)
})

test_that("MCI-009.5 DA run and export matrix preserves raw/friendly labels, intensity columns, and report metadata", {
  object <- module_ci_prot_da_object()
  contrasts_tbl <- rbind(
    module_ci_prot_da_contrasts("two_group"),
    module_ci_prot_da_contrasts("reversed")
  )
  paths <- module_ci_prot_da_paths()
  dir.create(paths$da_output_dir, recursive = TRUE)
  dir.create(paths$publication_graphs_dir, recursive = TRUE)

  da_data <- shiny::reactiveValues(
    da_results_list = NULL,
    contrasts_available = contrasts_tbl$contrasts,
    analysis_complete = FALSE,
    current_s4_object = object
  )
  workflow_data <- module_ci_prot_da_workflow(object)
  captured <- module_ci_prot_da_shiny_mocks()
  calls <- new.env(parent = emptyenv())
  calls$analysis <- list()
  calls$output <- list()

  analysis_stub <- function(theObject, contrasts_tbl, formula_string, da_q_val_thresh,
                            treat_lfc_cutoff, qvalue_column, raw_pvalue_column) {
    contrast_string <- contrasts_tbl$contrasts[[1L]]
    friendly <- sub("=.*$", "", contrast_string)
    raw <- sub("^[^=]+=", "", contrast_string)
    calls$analysis[[length(calls$analysis) + 1L]] <- list(
      full_format = contrast_string,
      friendly = friendly,
      raw = raw,
      formula_string = formula_string,
      da_q_val_thresh = da_q_val_thresh,
      treat_lfc_cutoff = treat_lfc_cutoff
    )
    list(
      theObject = theObject,
      da_proteins_long = module_ci_prot_da_result_table("standard", comparison = friendly),
      qvalue_warnings = NULL
    )
  }

  output_stub <- function(theObject, da_results_list_all_contrasts, uniprot_tbl,
                          da_output_dir, publication_graphs_dir, file_prefix,
                          args_row_id, gene_names_column, uniprot_id_column) {
    calls$output[[length(calls$output) + 1L]] <- list(
      theObject = theObject,
      da_results_list_all_contrasts = da_results_list_all_contrasts,
      da_output_dir = da_output_dir,
      publication_graphs_dir = publication_graphs_dir,
      file_prefix = file_prefix,
      args_row_id = args_row_id,
      gene_names_column = gene_names_column,
      uniprot_id_column = uniprot_id_column
    )
    TRUE
  }

  module_ci_prot_da_with_globals(list(contrasts_tbl = contrasts_tbl), {
    module_ci_prot_da_with_namespace_bindings(
      list(
        differentialAbundanceAnalysis = analysis_stub,
        outputDaResultsAllContrasts = output_stub,
        .capture_checkpoint = function(value, checkpoint_id, label) invisible(value)
      ),
      {
        testServer(
          module_ci_prot_da_run_module(da_data, workflow_data, paths),
          {
            session$setInputs(
              formula_string = "~ 0 + group",
              da_q_val_thresh = 0.05,
              treat_lfc_cutoff = 0.5
            )
            session$setInputs(run_da_analysis = 1)
            session$flushReact()

            expect_true(da_data$analysis_complete)
            expect_identical(
              vapply(calls$analysis, `[[`, character(1), "full_format"),
              contrasts_tbl$full_format
            )
            expect_identical(names(da_data$da_results_list$individual_contrasts), contrasts_tbl$contrasts)
            expect_identical(
              unique(da_data$da_results_list$da_proteins_long$comparison),
              contrasts_tbl$friendly_names
            )
            module_ci_prot_da_assert_result_schema(da_data$da_results_list$da_proteins_long)
            expect_identical(names(captured$select_updates$volcano_contrast), contrasts_tbl$friendly_names)
            expect_identical(names(captured$select_updates$heatmap_contrast), contrasts_tbl$friendly_names)
            expect_identical(names(captured$select_updates$table_contrast), contrasts_tbl$friendly_names)
            expect_equal(workflow_data$da_ui_params$q_value_threshold, 0.05)
            expect_equal(workflow_data$da_ui_params$log_fold_change_cutoff, 0.5)
            expect_true(workflow_data$da_ui_params$treat_enabled)
            expect_identical(workflow_data$tab_status$enrichment_analysis, "pending")
            expect_length(calls$output, 1L)
            expect_identical(calls$output[[1L]]$file_prefix, "da_proteins")
            expect_identical(names(calls$output[[1L]]$da_results_list_all_contrasts), contrasts_tbl$contrasts)
            expect_true("da_q_val_thresh" %in% names(calls$output[[1L]]$theObject@args$outputDaResultsAllContrasts))
          }
        )
      }
    )
  })

  export_captured <- new.env(parent = emptyenv())
  export_captured$vroom <- list()
  export_captured$xlsx <- list()
  module_ci_prot_da_with_namespace_bindings(
    list(
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) default_value,
      plotOneVolcanoNoVerticalLines = function(...) {
        ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point()
      }
    ),
    {
      testthat::local_mocked_bindings(
        vroom_write = function(x, path, ...) {
          export_captured$vroom[[basename(path)]] <<- x
          writeLines("mock tsv", path)
          invisible(path)
        },
        .package = "vroom"
      )
      testthat::local_mocked_bindings(
        write_xlsx = function(x, path, ...) {
          export_captured$xlsx[[basename(path)]] <<- x
          writeLines("mock xlsx", path)
          invisible(path)
        },
        .package = "writexl"
      )
      testthat::local_mocked_bindings(
        ggsave = function(filename, plot, ...) {
          writeLines("mock plot", filename)
          invisible(filename)
        },
        .package = "ggplot2"
      )

      export_object <- module_ci_prot_da_object()
      export_object@args$outputDaResultsAllContrasts <- list(
        da_q_val_thresh = 0.05,
        fdr_column = "fdr_qvalue",
        log2fc_column = "log2FC"
      )
      export_ok <- outputDaResultsAllContrasts(
        theObject = export_object,
        da_results_list_all_contrasts = list(
          "groupB-groupA" = module_ci_prot_da_result_list("standard", "B_vs_A"),
          "groupA-groupB" = module_ci_prot_da_result_list("no_significant", "A_vs_B")
        ),
        uniprot_tbl = data.frame(
          Entry = paste0("P", 1:6),
          gene_names = paste0("GENE", 1:6),
          stringsAsFactors = FALSE
        ),
        da_output_dir = paths$da_output_dir,
        publication_graphs_dir = paths$publication_graphs_dir,
        file_prefix = "da_proteins",
        args_row_id = "uniprot_acc",
        gene_names_column = "gene_names",
        uniprot_id_column = "Entry"
      )
      expect_true(export_ok)
    }
  )
  expect_true("da_proteins_groupB-groupA_long_annot.tsv" %in% names(export_captured$vroom))
  expect_true("da_proteins_groupA-groupB_long_annot.tsv" %in% names(export_captured$vroom))
  expect_true("da_proteins_groupB-groupA_long_annot.xlsx" %in% names(export_captured$xlsx))
  exported <- export_captured$vroom[["da_proteins_groupB-groupA_long_annot.tsv"]]
  module_ci_prot_da_assert_result_schema(exported)
  expect_true("gene_name" %in% names(exported))
  expect_true(any(grepl("^log2norm\\.", names(exported))))
})

test_that("MCI-009.6 DA output schema remains consumable by enrichment and summary/report consumers", {
  da_output_dir <- tempfile("module-ci-prot-da-output-")
  dir.create(da_output_dir, recursive = TRUE)
  da_table <- module_ci_prot_da_result_table("standard", "B_vs_A")
  readr::write_tsv(da_table, file.path(da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv"))
  writeLines("mock xlsx", file.path(da_output_dir, "da_proteins_groupB-groupA_long_annot.xlsx"))

  contrasts_tbl <- data.frame(
    contrasts = "groupB-groupA",
    full_format = "B_vs_A=groupB-groupA",
    friendly_names = "B_vs_A",
    stringsAsFactors = FALSE
  )
  design <- module_ci_prot_da_design()

  enrichment_ready <- createDAResultsForEnrichment(
    contrasts_tbl = contrasts_tbl,
    design_matrix = design,
    da_output_dir = da_output_dir
  )
  expect_s4_class(enrichment_ready, "da_results_for_enrichment")
  expect_identical(names(enrichment_ready@da_data), "groupB-groupA")
  expect_s3_class(enrichment_ready@da_data[["groupB-groupA"]], "spec_tbl_df")
  expect_true(all(c("uniprot_acc", "comparison", "fdr_qvalue", "log2FC") %in% names(enrichment_ready@da_data[["groupB-groupA"]])))

  contrast_choices <- buildProtEnrichContrastChoices(
    daResultsList = enrichment_ready@da_data,
    contrastsTbl = contrasts_tbl
  )
  expect_identical(contrast_choices$contrastsAvailable, "B_vs_A")
  resolved <- resolveProtEnrichSelectedDaResults(
    selectedContrast = "B_vs_A",
    daResultsData = enrichment_ready@da_data,
    contrastsTbl = contrasts_tbl
  )
  expect_s3_class(resolved$selectedDaResults, "spec_tbl_df")
  expect_identical(resolved$rawContrastName, "groupB-groupA")
  expect_equal(nrow(resolved$selectedDaResults), nrow(da_table))

  summary_da_files <- list.files(
    path = da_output_dir,
    pattern = "(da|de)_.+_long_annot\\.xlsx$",
    full.names = FALSE
  )
  expect_identical(summary_da_files, "da_proteins_groupB-groupA_long_annot.xlsx")
  summary_index_description <- basename(summary_da_files) |>
    stringr::str_remove("^da_") |>
    stringr::str_remove("^de_") |>
    stringr::str_remove("_long_annot\\.xlsx$")
  expect_identical(summary_index_description, "proteins_groupB-groupA")
})
