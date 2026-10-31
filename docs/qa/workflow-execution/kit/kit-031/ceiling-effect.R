# The intensity filter's retained rows on the e2e lipid lane (10 lipids x 6 samples, LCMS_Pos + LCMS_Neg) and on the
# KIT-031 declarations fixture, before (worktree) and after (tree) the ceiling() removal; usage: Rscript ceiling-effect.R <tree>
tree <- commandArgs(TRUE)[[1]]
suppressMessages(pkgload::load_all(tree, quiet = TRUE))
lane <- file.path(tree, "tests/testdata/e2e/lipid_canonical")
read_lane <- function(f) { d <- utils::read.delim(file.path(lane, f), check.names = FALSE, stringsAsFactors = FALSE); d }
design <- data.frame(Run = c("WT_1","WT_2","WT_3","KO_1","KO_2","KO_3"), group = rep(c("WT","KO"), each = 3), tech_rep_group = paste0(rep(c("WT","KO"), each = 3), "_", 1:3), stringsAsFactors = FALSE)
object <- methods::new("LipidomicsAssayData", lipid_data = list(LCMS_Pos = read_lane("lipidsearch_lcms_pos.txt"), LCMS_Neg = read_lane("lipidsearch_lcms_neg.txt")),
  design_matrix = design, lipid_id_column = "LipidName", annotation_id_column = "LipidClass", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group", args = list())
for (pct in c(1, 5, 10, 25)) for (prop in c(0.5)) {
  f <- suppressMessages(lipidIntensityFiltering(object, lipids_intensity_cutoff_percentile = pct, lipids_proportion_of_samples_below_cutoff = prop))
  cat(sprintf("e2e lane percentile %2d proportion %.1f: retained %s\n", pct, prop, paste(vapply(f@lipid_data, nrow, 1L), collapse = "/")))
}
source(file.path(tree, "tests/testthat/helper-workflow-lipid-declarations.R"))
fx <- lipidDeclarationsObject()
for (pct in c(1, 5)) { f <- suppressMessages(lipidIntensityFiltering(fx, lipids_intensity_cutoff_percentile = pct, lipids_proportion_of_samples_below_cutoff = 0.5))
  cat(sprintf("KIT-031 fixture percentile %2d: retained %s\n", pct, paste(vapply(f@lipid_data, nrow, 1L), collapse = "/"))) }
