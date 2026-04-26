# fidelity-coverage-compare: shared
library(testthat)

renderSharedUi <- function(ui) {
  htmltools::renderTags(ui)$html
}

expectSharedUiContains <- function(ui, patterns) {
  html <- renderSharedUi(ui)
  expect_true(nzchar(html))
  for (pattern in patterns) {
    expect_match(html, pattern, fixed = TRUE)
  }
  invisible(html)
}

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

makeWorkflowUiStub <- function(label) {
  force(label)
  function(id) {
    shiny::div(class = "workflow-ui-stub", `data-module-id` = id, label)
  }
}

test_that("application and standalone workflow UIs render their core controls", {
  expectSharedUiContains(
    app_ui(NULL),
    c("MultiScholaR", "dark_mode_toggle", "dynamic_menu")
  )

  expectSharedUiContains(
    mod_metab_norm_ui("metab-norm"),
    c(
      "Normalization Options",
      "RUV-III Batch Correction",
      "Run Normalization Pipeline"
    )
  )

  expectSharedUiContains(
    mod_prot_norm_ui("prot-norm"),
    c(
      "Normalization Options",
      "RUV-III Batch Correction",
      "Run Normalization &amp; RUV"
    )
  )

  expectSharedUiContains(
    mod_prot_da_ui("prot-da"),
    c(
      "DA Analysis Settings",
      "Volcano Plot",
      "DA Results Table"
    )
  )

  expectSharedUiContains(
    mod_lipid_da_ui("lipid-da"),
    c(
      "Differential Abundance Analysis",
      "Volcano Plot",
      "Download All Results"
    )
  )

  expectSharedUiContains(
    mod_metab_import_ui("metab-import"),
    c(
      "Metabolomics Data Import",
      "Step 1: Select Vendor Format",
      "sample_cols_pattern"
    )
  )

  expectSharedUiContains(
    mod_prot_import_ui("prot-import"),
    c(
      "Setup &amp; Import Data",
      "Step 1: Import Searched Data File",
      "format_override"
    )
  )

  expectSharedUiContains(
    mod_prot_design_builder_ui("prot-builder"),
    c(
      "Design Matrix Builder",
      "Rename Samples",
      "Assign Technical Replicates"
    )
  )

  expectSharedUiContains(
    mod_lipid_design_builder_ui("lipid-builder"),
    c(
      "Design Matrix Builder",
      "Information &amp; Assignments",
      "Defined Contrasts"
    )
  )
})

test_that("top-level workflow UIs render orchestrator tabs with stubbed module bodies", {
  package_ns <- asNamespace("MultiScholaR")

  localNamespaceBindings(
    package_ns,
    list(
      mod_prot_import_ui = makeWorkflowUiStub("proteomics-import"),
      mod_prot_design_ui = makeWorkflowUiStub("proteomics-design"),
      mod_prot_qc_ui = makeWorkflowUiStub("proteomics-qc"),
      mod_prot_norm_ui = makeWorkflowUiStub("proteomics-norm"),
      mod_prot_da_ui = makeWorkflowUiStub("proteomics-da"),
      mod_prot_enrich_ui = makeWorkflowUiStub("proteomics-enrich"),
      mod_prot_summary_ui = makeWorkflowUiStub("proteomics-summary"),
      mod_lipid_import_ui = makeWorkflowUiStub("lipidomics-import"),
      mod_lipid_design_ui = makeWorkflowUiStub("lipidomics-design"),
      mod_lipid_qc_ui = makeWorkflowUiStub("lipidomics-qc"),
      mod_lipid_norm_ui = makeWorkflowUiStub("lipidomics-norm"),
      mod_lipid_da_ui = makeWorkflowUiStub("lipidomics-da"),
      mod_lipid_summary_ui = makeWorkflowUiStub("lipidomics-summary"),
      mod_metab_import_ui = makeWorkflowUiStub("metabolomics-import"),
      mod_metab_design_ui = makeWorkflowUiStub("metabolomics-design"),
      mod_metab_qc_ui = makeWorkflowUiStub("metabolomics-qc"),
      mod_metab_norm_ui = makeWorkflowUiStub("metabolomics-norm"),
      mod_metab_da_ui = makeWorkflowUiStub("metabolomics-da"),
      mod_metab_summary_ui = makeWorkflowUiStub("metabolomics-summary")
    )
  )

  expectSharedUiContains(
    mod_proteomics_ui("proteomics"),
    c(
      "Setup &amp; Import",
      "Differential Expression",
      "Session Summary &amp; Report",
      "proteomics-enrich"
    )
  )

  expectSharedUiContains(
    mod_lipidomics_ui("lipidomics"),
    c(
      "Lipidomics Workflow",
      "Differential Analysis",
      "Summary &amp; Export",
      "lipidomics-summary"
    )
  )

  expectSharedUiContains(
    mod_metabolomics_ui("metabolomics"),
    c(
      "Metabolomics Workflow",
      "Differential Analysis",
      "Summary &amp; Export",
      "metabolomics-summary"
    )
  )
})
