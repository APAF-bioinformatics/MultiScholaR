# KIT-077: lipid and metab figures as store entities: every file is written from a saved, rooted rendered asset, a hit writes
# the stored bytes without rendering, the rebuildable types keep their plot data, a reopened owner serves its figures, and
# the two copy defects stay repaired.
library(testthat)

omicsFigureAssets <- function(manager) {
  found <- findSavedWorkflowEntities(manager$.context$store, list(kind = "rendered_asset"), max_matches = 500L)
  vapply(found$matches, function(reference) reference$entity_id, "")
}

omicsFigureProducts <- function(manager) {
  found <- findSavedWorkflowEntities(manager$.context$store, list(kind = "scientific_product"), max_matches = 500L)
  Filter(function(value) is.list(value) && !is.null(value$figure_type), lapply(found$matches, function(reference) {
    hydrateSavedWorkflowEntity(manager$.context$store, reference)$value
  }))
}

# A reopened owner on a session, as the app attaches one.
omicsFigureReopen <- function(manager, omic_type) {
  directory <- manager$.context$store$directory
  identity <- readSavedStoreIdentity(file.path(directory, "store.json"))
  closeDIAStateAdapter(manager)
  reopened <- newDIAStateAdapter(newWorkflowEntityStore(dirname(directory), identity$project_id, identity$workflow_id), reopen = TRUE)
  attachDIAStateSession(reopened, omicsOwnerSession(), new.env(parent = emptyenv()), omic_type = omic_type)
}

omicsFigureBytes <- function(path) readBin(path, "raw", n = file.size(path))

test_that("KIT-077: a miss renders, saves and roots; a hit writes the stored bytes without rendering", {
  manager <- omicsTestOwner(omicsOwnerLipidObject(), "s")$state_manager
  root <- withr::local_tempdir()
  path <- file.path(root, "figure.png")
  renders <- 0L
  render <- function(target) {
    renders <<- renders + 1L
    ggplot2::ggsave(target, ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, x)) + ggplot2::geom_point(), width = 2, height = 2, dpi = 72)
    TRUE
  }
  write <- function(content = "a") writeWorkflowOmicsFigure(path, render, "lipidomics", content = content, figure = "figure", width = 2, height = 2, dpi = 72)
  withWorkflowStateScope(manager, function() {
    expect_identical(write(), path)
    expect_identical(renders, 1L)
    assets <- omicsFigureAssets(manager)
    expect_length(assets, 1L)
    expect_true(assets %in% names(manager$.context$snapshot$products))
    stored <- omicsFigureBytes(path)
    unlink(path)
    expect_identical(write(), path)
    expect_identical(renders, 1L)
    expect_identical(omicsFigureBytes(path), stored)
    # Different content is a different figure.
    expect_identical(write("b"), path)
    expect_identical(renders, 2L)
  })
})

test_that("KIT-077: nothing to draw saves nothing; a figure needs an omics owner of its own omic to be saved", {
  workflow <- omicsTestOwner(omicsOwnerLipidObject(), "s")
  manager <- workflow$state_manager
  root <- withr::local_tempdir()
  path <- file.path(root, "figure.png")
  draw <- function(target) {
    grDevices::png(target); graphics::plot(1); grDevices::dev.off()
    TRUE
  }
  withWorkflowStateScope(manager, function() {
    expect_null(writeWorkflowOmicsFigure(path, function(target) NULL, "lipidomics", content = 1, figure = "none"))
    expect_length(omicsFigureAssets(manager), 0L)
    # The other omic's figure under this owner is rendered to its file only.
    expect_identical(writeWorkflowOmicsFigure(path, draw, "metabolomics", content = 1, figure = "other"), path)
    expect_length(omicsFigureAssets(manager), 0L)
  })
  # No owner in reach: as before.
  unlink(path)
  expect_identical(writeWorkflowOmicsFigure(path, draw, "lipidomics", content = 1, figure = "plain"), path)
  expect_true(file.exists(path))
  expect_null(workflowOmicsFigureOwner("lipidomics"))
})

test_that("KIT-077: the session finds its owner by omic, and not once it is closed", {
  skip_if_not(workflowOwnerCapable())
  source_dir <- withr::local_tempdir()
  session <- omicsOwnerSession()
  workflow <- new.env(parent = emptyenv())
  manager <- createDIASessionOwner(workflow, list(source_dir = source_dir), session, omic_type = "lipidomics")
  shiny::withReactiveDomain(session, {
    expect_identical(workflowOmicsFigureOwner("lipidomics"), manager)
    expect_null(workflowOmicsFigureOwner("metabolomics"))
  })
  closeDIAStateAdapter(manager)
  shiny::withReactiveDomain(session, expect_null(workflowOmicsFigureOwner("lipidomics")))
})

for (omic in c("lipid", "metab")) {
  omic_type <- if (identical(omic, "lipid")) "lipidomics" else "metabolomics"

  test_that(paste("KIT-077:", omic, "QC figures are saved with their plot data, and a reopened owner writes them without plotting"), {
    object <- omicsDaContrastsFixture(omic)
    manager <- omicsTestOwner(object, "s")$state_manager
    qc_dir <- withr::local_tempdir()
    paths <- if (identical(omic, "lipid")) list(lipid_qc_dir = qc_dir) else list(metabolite_qc_dir = qc_dir)
    generate <- if (identical(omic, "lipid")) generateLipidQcPlots else generateMetabQcPlots
    run <- function(owner) withWorkflowStateScope(owner, function() suppressWarnings(suppressMessages(
      generate(object, paths, stage = "post_norm", grouping_variable = "group"))))
    saved <- run(manager)
    expect_length(saved, 8L) # four figures for each of the two assays
    expect_true(all(file.exists(unlist(saved))))
    expect_length(omicsFigureAssets(manager), 8L)
    expect_true(all(omicsFigureAssets(manager) %in% names(manager$.context$snapshot$products)))
    # Plot data for the four rebuildable types, equal to the plot's own data.
    products <- withWorkflowStateScope(manager, function() omicsFigureProducts(manager))
    expect_setequal(vapply(products, `[[`, "", "figure_type"), c("pca", "rle", "density", "correlation"))
    expect_true(all(vapply(products, `[[`, "", "schema") == "omics_qc_plot_product_v1"))
    pca <- Filter(function(p) identical(p$figure_type, "pca") && identical(p$assay, "pos"), products)[[1L]]
    direct <- suppressWarnings(suppressMessages(plotPca(object, grouping_variable = "group", label_column = "", shape_variable = NULL, title = "", font_size = 8)))
    expect_equal(pca$data, as.data.frame(direct$pos$data))
    stored <- lapply(saved, omicsFigureBytes)
    unlink(unlist(saved))
    # Reopened: every file comes back byte for byte with the plotting functions refusing to run.
    reopened <- omicsFigureReopen(manager, omic_type)
    withr::defer(closeDIAStateAdapter(reopened))
    refuse <- function(...) stop("replotted")
    testthat::local_mocked_bindings(plotPca = refuse, plotRle = refuse, plotDensity = refuse, plotPearson = refuse, .package = "MultiScholaR")
    again <- run(reopened)
    expect_identical(lapply(again, omicsFigureBytes), stored)
    expect_length(omicsFigureAssets(reopened), 8L)
  })

  test_that(paste("KIT-077:", omic, "DA figures are saved, and a second export writes the stored bytes"), {
    skip_if_not_installed("ComplexHeatmap")
    object <- omicsDaContrastsFixture(omic)
    manager <- omicsTestOwner(object, "s")$state_manager
    results <- suppressWarnings(suppressMessages(omicsDaKernel(omic)(theObject = object,
      contrasts_tbl = omicsDaContrastsTable(c("groupB-groupA", "groupB-groupC")), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)))
    root <- withr::local_tempdir()
    writer <- if (identical(omic, "lipid")) outputLipidDaResultsAllContrasts else outputMetabDaResultsAllContrasts
    run <- function() withWorkflowStateScope(manager, function() suppressWarnings(suppressMessages(writer(da_results_list = results,
      da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "pub"), da_q_val_thresh = 0.05, lfc_threshold = 1))))
    expect_true(run())
    figures <- list.files(file.path(root, "pub"), recursive = TRUE, full.names = TRUE, pattern = "\\.(png|pdf)$")
    # Two assays x two contrasts x (volcano, heatmap) x (png, pdf), two combined PDFs, the NumSig bar plot.
    expect_length(figures, 19L)
    assets <- omicsFigureAssets(manager)
    expect_length(assets, 19L)
    expect_true(all(file.size(figures) > 0))
    stored <- lapply(figures, omicsFigureBytes)
    unlink(figures)
    # Identical PDF bytes (which carry a timestamp) and volcano PNGs (whose labels are placed at random) can only be the stored ones.
    expect_true(run())
    expect_identical(lapply(figures, omicsFigureBytes), stored)
    expect_setequal(omicsFigureAssets(manager), assets)
    products <- withWorkflowStateScope(manager, function() omicsFigureProducts(manager))
    expect_true(all(vapply(products, `[[`, "", "schema") == "omics_figure_v1"))
    expect_true(all(vapply(products, function(p) is.null(p$data), logical(1))))
  })
}

test_that("KIT-077: filtering QC plots and their grid are saved; a grid keeps its bytes only", {
  manager <- omicsTestOwner(omicsOwnerLipidObject(), "s")$state_manager
  root <- withr::local_tempdir()
  plots <- lapply(1:2, function(i) ggplot2::ggplot(data.frame(x = 1:3, y = 1:3 * i), ggplot2::aes(x, y)) + ggplot2::geom_point())
  names(plots) <- c("first", "second")
  withWorkflowStateScope(manager, function() {
    suppressMessages(saveLipidFilteringPlots(plots, "step", actual_save_dir = root, return_grid = TRUE))
    expect_setequal(list.files(root), c("step_first.png", "step_second.png", "step_combined_plots.png"))
    expect_length(omicsFigureAssets(manager), 3L)
    products <- omicsFigureProducts(manager)
    expect_setequal(vapply(products, `[[`, "", "figure_type"), c("filtering", "filtering_grid"))
    expect_true(all(vapply(products, function(p) is.null(p$data), logical(1))))
  })
})

test_that("KIT-077: the composite's png and pdf are saved from one render, and its rds is written on a hit too", {
  manager <- omicsTestOwner(omicsOwnerLipidObject(), "s")$state_manager
  root <- withr::local_tempdir()
  plot <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  input <- file.path(root, "input.txt"); writeLines("per-assay figure", input)
  renders <- 0L
  render <- function(...) { renders <<- renders + 1L; savePlot(...) }
  write <- function() withWorkflowStateScope(manager, function() writeWorkflowOmicsCompositeFigure(plot, root, "lipidomics_composite_QC_figure", width = 3, height = 3,
    plot_files = input, labels = list(rows = list(), columns = "a"), omic_type = "lipidomics", savePlotFn = render))
  files <- file.path(root, paste0("lipidomics_composite_QC_figure.", c("png", "pdf", "rds")))
  write()
  expect_identical(renders, 1L)
  expect_true(all(file.exists(files)))
  expect_length(omicsFigureAssets(manager), 2L)
  stored <- omicsFigureBytes(files[[2L]])
  unlink(files)
  write()
  expect_identical(renders, 1L)
  expect_true(all(file.exists(files)))
  expect_identical(omicsFigureBytes(files[[2L]]), stored)
})

test_that("KIT-077: the metab composite reaches QC_figures under the name the report reads", {
  dirs <- module_ci_metab_summary_paths()
  config <- module_ci_metab_summary_config("combined", ruv_mode = "automatic")
  module_ci_metab_summary_write_publication_inputs(dirs, config, include_optional = TRUE)
  workflow <- module_ci_metab_summary_workflow(config)
  failed <- module_ci_metab_summary_with_project_dirs(dirs, {
    utils::capture.output(result <- suppressWarnings(suppressMessages(copyToResultsSummary("metabolomics", "kit077",
      contrasts_tbl = workflow$contrasts_tbl, design_matrix = workflow$design_matrix, force = TRUE))))
    result
  })
  expect_true(file.exists(file.path(dirs$metabolomics$results_summary_dir, "QC_figures", "metabolomics_composite_QC_figure.png")))
  expect_false(any(vapply(failed, function(f) grepl("composite_QC_figure", paste(unlist(f), collapse = " ")), logical(1))))
})

test_that("KIT-077: the lipid composite copies with no spurious failure", {
  dirs <- module_ci_lipid_summary_paths()
  config <- module_ci_lipid_summary_config("pos_neg", ruv_mode = "automatic")
  module_ci_lipid_summary_write_publication_inputs(dirs, config, include_optional = TRUE)
  workflow <- module_ci_lipid_summary_workflow(config)
  failed <- module_ci_lipid_summary_with_project_dirs(dirs, {
    utils::capture.output(result <- suppressWarnings(suppressMessages(copyToResultsSummary("lipidomics", "kit077",
      contrasts_tbl = workflow$contrasts_tbl, design_matrix = workflow$design_matrix, force = TRUE))))
    result
  })
  expect_true(file.exists(file.path(dirs$lipidomics$results_summary_dir, "QC_figures", "lipidomics_composite_QC_figure.png")))
  expect_true(file.exists(file.path(dirs$lipidomics$results_summary_dir, "QC_figures", "lipidomics_composite_QC_figure.pdf")))
  expect_false(any(vapply(failed, function(f) grepl("composite_QC_figure", paste(unlist(f), collapse = " ")), logical(1))))
})
