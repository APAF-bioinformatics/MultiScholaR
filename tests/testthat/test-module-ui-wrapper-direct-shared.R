# fidelity-coverage-compare: shared
library(testthat)

renderDirectUi <- function(ui) {
  htmltools::renderTags(ui)$html
}

expectDirectUiContains <- function(ui, patterns) {
  html <- renderDirectUi(ui)
  expect_true(nzchar(html))

  for (pattern in patterns) {
    expect_match(html, pattern, fixed = TRUE)
  }

  invisible(html)
}

test_that("direct module UI wrappers render their current controls", {
  expectDirectUiContains(
    mod_metab_da_ui("metab-da-direct"),
    c(
      "Differential Abundance Analysis",
      "metab-da-direct-formula_string",
      "metab-da-direct-run_da_analysis",
      "metab-da-direct-heatmap_plot"
    )
  )

  expectDirectUiContains(
    mod_prot_enrich_ui("prot-enrich-direct"),
    c(
      "Enrichment Analysis Settings",
      "prot-enrich-direct-selected_contrast",
      "prot-enrich-direct-run_enrichment_analysis",
      "prot-enrich-direct-gprofiler_results_table"
    )
  )

  expectDirectUiContains(
    mod_prot_qc_ui("prot-qc-direct"),
    c(
      "Quality Control &amp; Filtering",
      "prot-qc-direct-dynamic_qc_tabs"
    )
  )

  expectDirectUiContains(
    mod_prot_summary_ui("prot-summary-direct"),
    c(
      "Session Summary &amp; Report Generation",
      "prot-summary-direct-save_workflow_args",
      "prot-summary-direct-enable_github"
    )
  )

  expectDirectUiContains(
    mod_lipid_design_ui("lipid-design-direct"),
    c(
      "Design Matrix Builder",
      "lipid-design-direct-show_import_modal",
      "lipid-design-direct-contrasts_preview"
    )
  )

  expectDirectUiContains(
    mod_metab_qc_duplicates_ui("metab-dup-direct"),
    c(
      "Duplicate Feature Resolution",
      "metab-dup-direct-detect_duplicates",
      "metab-dup-direct-filter_plot"
    )
  )

  expectDirectUiContains(
    mod_metab_qc_itsd_ui("metab-itsd-direct"),
    c(
      "Internal Standard QC",
      "metab-itsd-direct-analyze_is",
      "metab-itsd-direct-is_viz_tabs"
    )
  )

  expectDirectUiContains(
    mod_metab_summary_ui("metab-summary-direct"),
    c(
      "Session Summary &amp; Report Generation",
      "metab-summary-direct-generate_report",
      "metab-summary-direct-session_summary"
    )
  )

  expectDirectUiContains(
    mod_prot_design_ui("prot-design-direct"),
    c(
      "Design Matrix Builder",
      "prot-design-direct-show_import_modal",
      "prot-design-direct-contrasts_preview"
    )
  )

  expectDirectUiContains(
    mod_metab_design_ui("metab-design-direct"),
    c(
      "Design Matrix Builder",
      "metab-design-direct-show_import_modal",
      "metab-design-direct-assays_preview"
    )
  )

  expectDirectUiContains(
    mod_lipid_summary_ui("lipid-summary-direct"),
    c(
      "Session Summary &amp; Report Generation",
      "lipid-summary-direct-push_to_github",
      "lipid-summary-direct-session_summary"
    )
  )
})
