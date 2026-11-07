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

# KIT-102 W4: the shared shell (`mod_omic_workflow_ui`) and the import frame (`mod_omic_import_ui`).

# Every tab module UI stubbed; the stub shows the module id and, for trx, the profile it was given.
localWorkflowTabUiStubs <- function(.local_envir = parent.frame()) {
  stub <- function(label) {
    force(label)
    function(id, profile = NULL) shiny::div(class = "stub", `data-module-id` = id, `data-profile` = profile$omic_type, label)
  }
  names <- c(paste0("mod_prot_", c("import", "design", "qc", "norm", "da", "enrich", "summary"), "_ui"),
    paste0("mod_lipid_", c("import", "design", "qc", "norm", "da", "summary"), "_ui"),
    paste0("mod_metab_", c("import", "design", "qc", "norm", "da", "summary"), "_ui"))
  localNamespaceBindings(asNamespace("MultiScholaR"), stats::setNames(lapply(names, stub), names), .local_envir = .local_envir)
}

# The tabset's random ids (`tab-<n>-<i>`, `data-tabsetid`) normalised.
normalisedWorkflowUi <- function(ui) {
  html <- renderSharedUi(ui)
  gsub('data-tabsetid="[0-9]+"', 'data-tabsetid="N"', gsub("tab-[0-9]+-", "tab-N-", html))
}

test_that("KIT-102 W4: the four workflow UIs are byte-identical to the pre-KIT-102 shells (stubbed tabs)", {
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
  localWorkflowTabUiStubs()
  # sha256 of each shell's normalised HTML, recorded on HEAD 60334e4b (before the shared shell) with the same stubs.
  expected <- c(
    proteomics = "edc57bef8434a45bf0d4a4da439482cc323fa36119663861039c7c4e4077c6c5",
    lipidomics = "69f51d7ed330568a7050ee7ab546e84b08d413a535deca16a60f6b46fbfe4de6",
    metabolomics = "276b077dc30ecf2998b3b57221fc6b520d2647efd19a4d4cc02e7f2f070d8473",
    transcriptomics = "0fd18964c6acf4be24677f6e0626064c67b214993450b05b99e4506ef0ebb96e")
  for (omic in names(expected)) {
    ui <- get(paste0("mod_", omic, "_ui"), envir = asNamespace("MultiScholaR"))
    html <- normalisedWorkflowUi(ui(paste0(omic, "_workflow")))
    expect_identical(digest::digest(html, algo = "sha256", serialize = FALSE), expected[[omic]], info = html)
  }
})

test_that("KIT-102 W4: the shared shell renders each omic's tabs from its configuration", {
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
  localWorkflowTabUiStubs()
  for (omic in c("proteomics", "lipidomics", "metabolomics", "transcriptomics")) {
    wrapper <- get(paste0("mod_", omic, "_ui"), envir = asNamespace("MultiScholaR"))
    expect_identical(normalisedWorkflowUi(mod_omic_workflow_ui("w", omicGuiConfig(omic))), normalisedWorkflowUi(wrapper("w")),
      info = omic)
  }
  # The values come from the configuration, not from the shell.
  config <- omicGuiConfig("lipidomics")
  config$shell$tabs[[3L]]$title <- " Checks"
  config$shell$tabs[[3L]]$testid <- "lipid-tab-checks"
  html <- renderSharedUi(mod_omic_workflow_ui("w", config))
  expect_match(html, "Checks", fixed = TRUE)
  expect_match(html, "data-testid=\"lipid-tab-checks\"", fixed = TRUE)
  expect_false(grepl("Quality Control", html, fixed = TRUE))
  # An invalid configuration is refused before anything renders.
  config$shell$tabset_id <- NULL
  expect_error(mod_omic_workflow_ui("w", config), class = "omic_gui_config_invalid")
})

test_that("KIT-102 W4: the component registry resolves every configured tab and refuses an unknown id", {
  for (omic in c("proteomics", "lipidomics", "metabolomics", "transcriptomics")) {
    config <- omicGuiConfig(omic)
    for (id in c(vapply(config$shell$tabs, `[[`, character(1), "component"), config$shell$extras)) {
      component <- omicGuiComponent(id)
      expect_true(is.character(component$ui) && length(component$ui) == 1L, info = id)
      expect_true(is.character(component$server) && length(component$server) == 1L, info = id)
      expect_true(exists(component$ui, envir = asNamespace("MultiScholaR"), mode = "function"), info = id)
      expect_true(exists(component$server, envir = asNamespace("MultiScholaR"), mode = "function"), info = id)
    }
  }
  expect_identical(omicGuiComponent("trx_norm")$ui, "mod_lipid_norm_ui")
  expect_error(omicGuiComponent("lipid_enrich"), class = "omic_gui_config_invalid")
})

test_that("KIT-102 W4: the import frame renders the title and the configured source panel, with no DOM id of its own", {
  seen <- new.env(parent = emptyenv())
  # The panels are phase C's (KIT-102 W5, W6): their UI functions are injected through the panel registry.
  panels <- omicImportSourcePanel
  stubs <- list(
    feature_table = function(ns, config) {
      seen$feature_table <- config$shell$omic_type
      shiny::tagList(shiny::fluidRow(shiny::column(6, shiny::textInput(ns("a"), "A"))), shiny::hr(), shiny::div("footer"))
    },
    proteomics = function(ns, config) shiny::div(id = ns("prot_panel"), "proteomics panel"))
  local_mocked_bindings(omicImportSourcePanel = function(source_panel) {
    utils::modifyList(panels(source_panel), list(ui = stubs[[source_panel]]))
  }, .package = "MultiScholaR")
  for (omic in c("lipidomics", "metabolomics", "transcriptomics")) {
    config <- omicGuiConfig(omic)
    html <- renderSharedUi(mod_omic_import_ui("w-import", config))
    expected <- renderSharedUi(shiny::tagList(shinyjs::useShinyjs(), shiny::fluidRow(shiny::column(12, shiny::wellPanel(
      shiny::h3(config$import$title), shiny::fluidRow(shiny::column(6, shiny::textInput("w-import-a", "A"))), shiny::hr(),
      shiny::div("footer"))))))
    expect_identical(html, expected, info = omic)
    expect_identical(seen$feature_table, omic)
  }
  html <- renderSharedUi(mod_omic_import_ui("w-setup_import", omicGuiConfig("proteomics")))
  expect_match(html, "<h3>Setup &amp; Import Data</h3>", fixed = TRUE)
  expect_match(html, "id=\"w-setup_import-prot_panel\"", fixed = TRUE)
  ids <- regmatches(html, gregexpr(" id=\"[^\"]+\"", html))[[1L]]
  expect_identical(ids, " id=\"w-setup_import-prot_panel\"")
})
