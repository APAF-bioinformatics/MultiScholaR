# fidelity-coverage-compare: shared
library(testthat)
library(MultiScholaR)

test_that("generateMetabDAVolcanoPlotGlimma handles mock data", {
  # Mock data similar to repro script
  assay_data <- data.frame(
    metabolite_id = paste0("M", 1:50),
    metabolite_name = paste0("Metab", 1:50),
    S1 = rnorm(50), S2 = rnorm(50), S3 = rnorm(50), S4 = rnorm(50),
    check.names = FALSE
  )
  design <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B")
  )
  theObject <- new("MetaboliteAssayData",
    metabolite_data = list(Assay1 = assay_data),
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "metabolite_name"
  )
  da_metabolites_long <- data.frame(
    metabolite_id = paste0("M", 1:50),
    metabolite_name = paste0("Metab", 1:50),
    assay = "Assay1",
    comparison = "groupB - groupA",
    friendly_name = "B_vs_A",
    logFC = rnorm(50),
    raw_pvalue = runif(50),
    fdr_qvalue = runif(50),
    significant = "Not sig."
  )
  
  # Dummy fit
  library(limma)
  dummy_mat <- matrix(rnorm(200), 50, 4)
  rownames(dummy_mat) <- paste0("M", 1:50)
  fit <- lmFit(dummy_mat, model.matrix(~0+group, data=design))
  fit2 <- contrasts.fit(fit, makeContrasts(groupB-groupA, levels=model.matrix(~0+group, data=design)))
  fit2 <- eBayes(fit2)
  
  da_results_list <- list(
    theObject = theObject,
    da_metabolites_long = da_metabolites_long,
    contrasts_results = list(Assay1 = list(fit.eb = fit2))
  )
  
  p <- generateMetabDAVolcanoPlotGlimma(da_results_list, "B_vs_A", "Assay1")
  expect_s3_class(p, "glimmaXY")
})

test_that("generateLipidDAVolcanoPlotGlimma handles mock data", {
  # Mock data similar to metab
  assay_data <- data.frame(
    lipid_id = paste0("L", 1:50),
    lipid_name = paste0("Lipid", 1:50),
    S1 = rnorm(50), S2 = rnorm(50), S3 = rnorm(50), S4 = rnorm(50),
    check.names = FALSE
  )
  design <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B")
  )
  theObject <- new("LipidomicsAssayData",
    lipid_data = list(Assay1 = assay_data),
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    lipid_id_column = "lipid_id",
    annotation_id_column = "lipid_name"
  )
  da_lipids_long <- data.frame(
    lipid_id = paste0("L", 1:50),
    lipid_name = paste0("Lipid", 1:50),
    assay = "Assay1",
    comparison = "groupB - groupA",
    friendly_name = "B_vs_A",
    logFC = rnorm(50),
    raw_pvalue = runif(50),
    fdr_qvalue = runif(50),
    significant = "Not sig."
  )
  
  # Dummy fit
  library(limma)
  dummy_mat <- matrix(rnorm(200), 50, 4)
  rownames(dummy_mat) <- paste0("L", 1:50)
  fit <- lmFit(dummy_mat, model.matrix(~0+group, data=design))
  fit2 <- contrasts.fit(fit, makeContrasts(groupB-groupA, levels=model.matrix(~0+group, data=design)))
  fit2 <- eBayes(fit2)
  
  da_results_list <- list(
    theObject = theObject,
    da_lipids_long = da_lipids_long,
    contrasts_results = list(Assay1 = list(fit.eb = fit2))
  )
  
  p <- generateLipidDAVolcanoPlotGlimma(da_results_list, "B_vs_A", "Assay1")
  expect_s3_class(p, "glimmaXY")
})
