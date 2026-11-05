# KIT-087: PCA with missing values, ours against mixOmics, on the matrices nipals-compare.R extracted from stored states.
# Usage (R_LIBS = the kit_087 library): pca-compare.R call <input rds> <old|new> <features: 0 for all> | compare <input rds> <features>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); mode <- args[[1L]]; input <- readRDS(args[[2L]])
features <- as.integer(args[[if (mode == "call") 4L else 3L]]); Y <- input$Y
if (features > 0L && features < ncol(Y)) Y <- Y[, seq_len(features), drop = FALSE]
stem <- sub("\\.rds$", sprintf(".pca%d", features), args[[2L]])
fields <- function(p) list(variates = p$variates$X, proportion = p$prop_expl_var$X, sdev = p$sdev, loadings = p$loadings$X, var_tot = p$var.tot)
if (mode == "call") {
  which_one <- args[[3L]]; gc(reset = TRUE)
  elapsed <- system.time(result <- suppressWarnings(suppressMessages(if (which_one == "old") mixOmics::pca(Y, ncomp = 2L) else MultiScholaR:::pcaWithMissing(Y, ncomp = 2L))))[["elapsed"]]
  result <- fields(result); saveRDS(result, paste0(stem, ".", which_one, ".result"), compress = FALSE)
  cat(sprintf("PCACALL %s %s %s features %d elapsed %.1f s digest %s\n", input$label, input$assay, which_one, ncol(Y), elapsed, digest::digest(result)))
} else {
  old <- readRDS(paste0(stem, ".old.result")); new <- readRDS(paste0(stem, ".new.result"))
  cat(sprintf("PCACOMPARE %s %s samples %d features %d missing %.4f identical %s\n", input$label, input$assay, nrow(Y), ncol(Y), mean(is.na(Y)), identical(old, new)))
}
