# KIT-087: imputation on the first <features> features of a stored matrix, for assays where mixOmics at full size passed its 3,600 s cap.
# Usage (R_LIBS = the kit_087 library): nipals-subset.R <input rds> <features>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); input <- readRDS(args[[1L]]); features <- as.integer(args[[2L]])
Y <- input$Y[, seq_len(features), drop = FALSE]
quietly <- function(x) suppressWarnings(suppressMessages(x))
old_seconds <- system.time(old <- quietly(mixOmics::impute.nipals(Y, ncomp = input$ncomp)))[["elapsed"]]
new_seconds <- system.time(new <- quietly(MultiScholaR:::imputeNipalsMatrix(Y, ncomp = input$ncomp)))[["elapsed"]]
cat(sprintf("SUBSET %s %s samples %d features %d ncomp %d missing %.4f identical %s old %.1f s new %.1f s\n", input$label, input$assay, nrow(Y), ncol(Y), input$ncomp, mean(is.na(Y)), identical(old, new), old_seconds, new_seconds))
