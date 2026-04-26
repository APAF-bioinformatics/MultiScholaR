# fidelity-coverage-compare: shared
library(testthat)
library(ggplot2)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

getCategoricalColourPalette <- get(
  "getCategoricalColourPalette",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
getOneContinousPalette <- get(
  "getOneContinousPalette",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
getContinousColourRules <- get(
  "getContinousColourRules",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
getCategoricalAndContinuousColourRules <- get(
  "getCategoricalAndContinuousColourRules",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
apafTheme <- get(
  "apafTheme",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
get_color_palette <- get(
  "get_color_palette",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
printPValuesDistribution <- get(
  "printPValuesDistribution",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
gg_save_logging <- get(
  "gg_save_logging",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
summarizeQCPlot <- get(
  "summarizeQCPlot",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

if (!methods::isClass("MockQCFigure")) {
  methods::setClass(
    "MockQCFigure",
    slots = c(
      rle_plots = "list",
      pca_plots = "list",
      density_plots = "list",
      pearson_plots = "list"
    )
  )
}

test_that("general plotting palette helpers preserve colour rule generation", {
  metadata_tbl <- data.frame(
    Run = c("S1", "S2", "S3"),
    Group = c("A", "B", NA),
    Score = c(1, 2, NA),
    Machine = c("M1", "M1", "M2"),
    stringsAsFactors = FALSE
  )

  categorical_palette <- getCategoricalColourPalette()
  expect_length(categorical_palette, 74)

  continuous_palette <- getOneContinousPalette(
    metadata_tbl = metadata_tbl,
    column_name = "Score",
    palette_name = "Blues",
    na_colour = "grey90"
  )
  expect_identical(unname(names(continuous_palette)[1:2]), c("1", "2"))
  expect_true(any(is.na(names(continuous_palette))))
  expect_identical(tail(names(continuous_palette), 1), "NA")
  expect_identical(tail(unname(continuous_palette), 1), "grey90")

  continuous_rules <- getContinousColourRules(
    metadata_tbl = metadata_tbl,
    metadata_column_labels = c("Score Label"),
    metadata_column_selected = c("Score"),
    continous_scale_columns = c("Score"),
    na_colour = "grey80"
  )
  expect_named(continuous_rules, "Score Label")
  expect_identical(unname(names(continuous_rules[[1]])[1:2]), c("1", "2"))
  expect_identical(tail(names(continuous_rules[[1]]), 1), "NA")
  expect_identical(tail(unname(continuous_rules[[1]]), 1), "grey80")

  combined_rules <- makeFunctionWithOverrides(
    getCategoricalAndContinuousColourRules,
    list(
      getCategoricalColourRules = function(metadata_tbl,
                                           metadata_column_labels,
                                           metadata_column_selected,
                                           categorical_columns,
                                           ms_machine_column,
                                           columns_to_exclude,
                                           na_colour = "white") {
        list(
          "Group Label" = c(A = "#111111", B = "#222222", "NA" = na_colour),
          "Machine Label" = c(M1 = "#333333", M2 = "#444444")
        )
      }
    )
  )(
    metadata_tbl = metadata_tbl,
    metadata_column_labels = c("Group Label", "Score Label", "Machine Label"),
    metadata_column_selected = c("Group", "Score", "Machine"),
    categorical_columns = c("Group", "Machine"),
    continous_scale_columns = c("Score"),
    ms_machine_column = "Machine",
    sample_id_column = Run,
    columns_to_exclude = "Machine",
    na_colour = "grey70"
  )

  expect_true("Group Label" %in% names(combined_rules))
  expect_true("Score Label" %in% names(combined_rules))
  expect_false("Machine Label" %in% names(combined_rules))

  theme_obj <- apafTheme()
  expect_s3_class(theme_obj, "theme")

  generated_palette <- get_color_palette(4, "#00AA55")
  expect_length(generated_palette, 4)
  expect_identical(generated_palette[[1]], "#00AA55")
})

test_that("general plotting reporting helpers preserve histogram, save logging, and QC summary behavior", {
  selected_data <- data.frame(
    raw_pvalue = c(0.001, 0.02, 0.2, 0.5),
    comparison = c("A_vs_B", "A_vs_B", "A_vs_C", "A_vs_C"),
    is_ruv_applied = c("yes", "yes", "no", "no"),
    stringsAsFactors = FALSE
  )

  faceted_plot <- printPValuesDistribution(
    selected_data = selected_data,
    p_value_column = raw_pvalue,
    formula_string = "is_ruv_applied ~ comparison"
  )
  expect_s3_class(faceted_plot, "ggplot")
  expect_identical(faceted_plot$labels$x, "P-value")

  plain_plot <- printPValuesDistribution(
    selected_data = selected_data,
    p_value_column = raw_pvalue,
    formula_string = NA_character_
  )
  expect_s3_class(plain_plot, "ggplot")

  file_name_part <- tempfile("plot-output-")
  withr::defer({
    unlink(paste0(file_name_part, ".png"), force = TRUE)
    unlink(paste0(file_name_part, ".pdf"), force = TRUE)
  })

  gg_save_logging(
    input_plot = plain_plot,
    file_name_part = file_name_part,
    plots_format = c(".png", ".pdf"),
    width = 4,
    height = 5
  )

  expect_true(file.exists(paste0(file_name_part, ".png")))
  expect_true(file.exists(paste0(file_name_part, ".pdf")))

  qc_figure <- methods::new(
    "MockQCFigure",
    rle_plots = list(rle = ggplot(selected_data, aes(raw_pvalue, raw_pvalue))),
    pca_plots = list(pca = ggplot(selected_data, aes(raw_pvalue, raw_pvalue))),
    density_plots = list(density = ggplot(selected_data, aes(raw_pvalue, raw_pvalue))),
    pearson_plots = list(pearson = ggplot(selected_data, aes(raw_pvalue, raw_pvalue)))
  )

  summary_output <- capture.output(summarizeQCPlot(qc_figure))
  expect_true(any(grepl("^RLE Plots:", summary_output)))
  expect_true(any(grepl("^PCA Plots:", summary_output)))
  expect_true(any(grepl("^Density Plots:", summary_output)))
  expect_true(any(grepl("^Pearson Correlation Plots:", summary_output)))
})
