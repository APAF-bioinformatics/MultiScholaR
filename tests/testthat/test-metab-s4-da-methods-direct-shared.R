# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

newMetabS4DaObject <- function(group_values = c("1A", "1A", "2B", "2B"), args = list()) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(data.frame(
      Name = c("M1", "M2"),
      S1 = c(10, 20),
      S2 = c(12, 18),
      S3 = c(30, 5),
      S4 = c(32, 7),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )),
    metabolite_id_column = "Name",
    annotation_id_column = "Name",
    database_identifier_type = "Name",
    internal_standard_regex = "",
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = group_values,
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = args
  )
}

test_that("metabolite S4 DA methods preserve list validation and helper orchestration", {
  list_method <- makeFunctionWithOverrides(
    methods::selectMethod("differentialAbundanceAnalysis", "list"),
    list(
      differentialAbundanceAnalysisHelper = function(obj, ...) {
        paste("ok", obj@group_id)
      }
    )
  )

  expect_error(list_method(list("bad-object")), "MetaboliteAssayData", fixed = TRUE)

  list_result <- list_method(list(Assay1 = newMetabS4DaObject(group_values = c("A", "A", "B", "B"))))
  expect_identical(list_result$Assay1, "ok group")
})

test_that("metabolite S4 DA helper preserves numeric-group remapping and namespace fallback", {
  helper_calls <- new.env(parent = emptyenv())

  helper_method <- makeFunctionWithOverrides(
    methods::selectMethod("differentialAbundanceAnalysisHelper", "MetaboliteAssayData"),
    list(
      checkParamsObjectFunctionSimplify = function(theObject, name, default) {
        value <- theObject@args[[name]]
        if (is.null(value)) default else value
      },
      updateParamInObject = function(theObject, name) theObject,
      runTestsContrasts = function(data_matrix,
                                   contrast_strings,
                                   design_matrix,
                                   formula_string,
                                   treat_lfc_cutoff,
                                   eBayes_trend,
                                   eBayes_robust) {
        helper_calls$rownames <- rownames(data_matrix)
        helper_calls$contrast_strings <- contrast_strings
        helper_calls$design_groups <- design_matrix$group
        helper_calls$formula_string <- formula_string
        helper_calls$treat_lfc_cutoff <- treat_lfc_cutoff
        helper_calls$eBayes_trend <- eBayes_trend
        helper_calls$eBayes_robust <- eBayes_robust

        list(
          fit.eb = structure(list(marker = "fit"), class = "mock_fit"),
          results = list(
            "grp_2B-grp_1A" = data.frame(
              comparison = c("grp_2B-grp_1A", "grp_2B-grp_1A"),
              logFC = c(1.5, -1.2),
              raw_pvalue = c(0.01, 0.02),
              fdr_qvalue = c(0.02, 0.03),
              row.names = c("M1", "M2"),
              check.names = FALSE
            )
          )
        )
      },
      rownames_to_column = tibble::rownames_to_column,
      message = function(...) invisible(NULL),
      print = function(...) invisible(NULL)
    )
  )

  result_object <- helper_method(newMetabS4DaObject(args = list(
    contrasts_tbl = data.frame(contrasts = "2B-1A", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group",
    treat_lfc_cutoff = "0.5",
    eBayes_trend = "FALSE",
    eBayes_robust = "TRUE",
    args_group_pattern = "(\\d+)"
  )))

  expect_s4_class(result_object, "MetabolomicsDifferentialAbundanceResults")
  expect_identical(helper_calls$rownames, c("M1", "M2"))
  expect_identical(helper_calls$contrast_strings, "grp_2B-grp_1A")
  expect_identical(helper_calls$design_groups, c("grp_1A", "grp_1A", "grp_2B", "grp_2B"))
  expect_identical(helper_calls$formula_string, "~ 0 + group")
  expect_identical(helper_calls$treat_lfc_cutoff, 0.5)
  expect_false(helper_calls$eBayes_trend)
  expect_true(helper_calls$eBayes_robust)
  expect_identical(result_object@contrasts_results_table[[1]]$comparison, c("2B-1A", "2B-1A"))
  expect_identical(result_object@theObject@design_matrix$group, c("grp_1A", "grp_1A", "grp_2B", "grp_2B"))
})
