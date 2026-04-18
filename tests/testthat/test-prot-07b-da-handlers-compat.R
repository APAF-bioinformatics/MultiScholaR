library(testthat)
library(shiny)

make_mock_protein_da_object <- function(args = list()) {
  ProteinQuantitativeData(
    protein_quant_table = data.frame(
      uniprot_acc = c("P1", "P2"),
      S1 = c(10, 11),
      S2 = c(20, 21),
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
      uniprot_acc = c("P1", "P2"),
      gene_names = c("GENE1", "GENE2"),
      stringsAsFactors = FALSE
    ),
    args = args
  )
}

with_global_contrasts_tbl <- function(value, code) {
  had_existing <- exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
  if (had_existing) {
    old_value <- get("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)
  }

  if (is.null(value)) {
    if (had_existing) {
      rm("contrasts_tbl", envir = .GlobalEnv)
    }
  } else {
    assign("contrasts_tbl", value, envir = .GlobalEnv)
  }

  on.exit({
    if (had_existing) {
      assign("contrasts_tbl", old_value, envir = .GlobalEnv)
    } else if (exists("contrasts_tbl", envir = .GlobalEnv, inherits = FALSE)) {
      rm("contrasts_tbl", envir = .GlobalEnv)
    }
  }, add = TRUE)

  force(code)
}

with_function_overrides <- function(env, replacements, code) {
  names_vec <- names(replacements)
  old_values <- lapply(names_vec, function(name) get(name, envir = env, inherits = FALSE))
  names(old_values) <- names_vec
  locked <- vapply(names_vec, function(name) bindingIsLocked(name, env), logical(1))

  for (name in names_vec) {
    if (bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    assign(name, replacements[[name]], envir = env)
  }

  on.exit({
    for (name in names_vec) {
      if (bindingIsLocked(name, env)) {
        unlockBinding(name, env)
      }
      assign(name, old_values[[name]], envir = env)
      if (locked[[name]]) {
        lockBinding(name, env)
      }
    }
  }, add = TRUE)

  force(code)
}

test_that("writeInteractiveVolcanoPlotProteomicsMain accepts legacy de_* aliases", {
  obj <- make_mock_protein_da_object()
  legacy_results <- list(
    theObject = obj,
    de_proteins_long = data.frame(
      uniprot_acc = c("P1", "P2"),
      comparison = c("A_vs_B", "A_vs_B"),
      log2FC = c(1.5, -1.2),
      fdr_qvalue = c(0.01, 0.02),
      raw_pvalue = c(0.001, 0.002),
      stringsAsFactors = FALSE
    ),
    contrasts_results = list(fit.eb = "fit-placeholder")
  )

  captured <- NULL
  wrapper_env <- environment(writeInteractiveVolcanoPlotProteomicsMain)
  with_function_overrides(
    wrapper_env,
    list(
      writeInteractiveVolcanoPlotProteomics = function(...) {
        captured <<- list(...)
        captured
      },
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) {
        default_value
      },
      updateParamInObject = function(theObject, param_name_string) {
        theObject
      }
    ),
    {
      result <- writeInteractiveVolcanoPlotProteomicsMain(
        da_analysis_results_list = NULL,
        de_analysis_results_list = legacy_results,
        theObject = obj,
        uniprot_tbl = NULL,
        args_row_id = "uniprot_acc",
        de_q_val_thresh = 0.01
      )

      expect_type(result, "list")
      expect_identical(result[[1]], legacy_results$de_proteins_long)
      expect_identical(result$fit.eb, "fit-placeholder")
      expect_equal(result$da_q_val_thresh, 0.01)
      expect_equal(rownames(result$counts_tbl), c("P1", "P2"))
      expect_equal(unname(result$groups), c("A", "B"))
    }
  )
})

test_that("writeInteractiveVolcanoPlotProteomicsMain errors clearly for missing inputs", {
  obj <- make_mock_protein_da_object()
  wrapper_env <- environment(writeInteractiveVolcanoPlotProteomicsMain)
  with_function_overrides(
    wrapper_env,
    list(
      checkParamsObjectFunctionSimplify = function(theObject, param_name_string, default_value = NULL) {
        default_value
      },
      updateParamInObject = function(theObject, param_name_string) {
        theObject
      }
    ),
    {
      expect_error(
        writeInteractiveVolcanoPlotProteomicsMain(
          da_analysis_results_list = NULL,
          de_analysis_results_list = NULL,
          theObject = obj,
          uniprot_tbl = NULL
        ),
        "A DA/DE analysis results list must be supplied."
      )

      expect_error(
        writeInteractiveVolcanoPlotProteomicsMain(
          da_analysis_results_list = list(theObject = obj),
          theObject = obj,
          uniprot_tbl = NULL
        ),
        "Results list does not contain the expected DA/DE volcano inputs."
      )
    }
  )
})

test_that("da_server_init_handlers auto-generates prefixed contrasts from workflow state", {
  with_global_contrasts_tbl(NULL, {
    selected_tab <- shiny::reactiveVal(NULL)
    obj <- make_mock_protein_da_object(
      args = list(
        deAnalysisParameters = list(
          formula_string = "~ 0 + group"
        )
      )
    )
    workflow_data <- shiny::reactiveValues(
      state_manager = list(
        current_state = "normalized",
        getState = function(state) obj
      ),
      state_update_trigger = NULL
    )

    init_test_module <- function(id, workflow_data, selected_tab) {
      shiny::moduleServer(id, function(input, output, session) {
        da_data <- shiny::reactiveValues(
          da_results_list = NULL,
          contrasts_available = NULL,
          analysis_complete = FALSE,
          current_s4_object = NULL,
          formula_from_s4 = NULL,
          current_row_clusters = NULL,
          current_col_clusters = NULL
        )

        da_server_init_handlers(input, output, session, da_data, workflow_data, selected_tab)
        da_data
      })
    }

    testServer(
      init_test_module,
      args = list(workflow_data = workflow_data, selected_tab = selected_tab),
      {
        session$setInputs(formula_string = "~ 0 + group")
        selected_tab("da")
        session$flushReact()

        expect_s4_class(da_data$current_s4_object, "ProteinQuantitativeData")
        expect_equal(da_data$formula_from_s4, "~ 0 + group")
        expect_equal(as.vector(da_data$contrasts_available), "groupA-groupB")
      }
    )
  })
})

test_that("da_server_init_handlers normalizes user contrasts to group-prefixed format", {
  with_global_contrasts_tbl(
    data.frame(comparison = "A-B", stringsAsFactors = FALSE),
    {
      selected_tab <- shiny::reactiveVal(NULL)
      obj <- make_mock_protein_da_object(
        args = list(
          deAnalysisParameters = list(
            formula_string = "~ 0 + group"
          )
        )
      )
      workflow_data <- shiny::reactiveValues(
        state_manager = list(
          current_state = "normalized",
          getState = function(state) obj
        ),
        state_update_trigger = NULL
      )

      init_test_module <- function(id, workflow_data, selected_tab) {
        shiny::moduleServer(id, function(input, output, session) {
          da_data <- shiny::reactiveValues(
            da_results_list = NULL,
            contrasts_available = NULL,
            analysis_complete = FALSE,
            current_s4_object = NULL,
            formula_from_s4 = NULL,
            current_row_clusters = NULL,
            current_col_clusters = NULL
          )

          da_server_init_handlers(input, output, session, da_data, workflow_data, selected_tab)
          da_data
        })
      }

      testServer(
        init_test_module,
        args = list(workflow_data = workflow_data, selected_tab = selected_tab),
        {
          session$setInputs(formula_string = "~ 0 + group")
          selected_tab("da")
          session$flushReact()

          expect_equal(as.vector(da_data$contrasts_available), "groupA-groupB")
        }
      )
    }
  )
})
