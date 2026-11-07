# KIT-101 revision 2.3, F-19 follow-up: does a between-sample normalisation chosen on the DIA_LIMPA route reach its DA?
# The route's protein normalisation step calls normaliseBetweenSamples on the ProteinQuantitativeData
# (computeDIAReferenceProteinNormalisation); the DA's dpcDE branch fits limpa_dpc_quant_results$quantified_elist.
# This applies the method the DIA route uses (cyclicloess) to the pre-RUV state dumped by f19_pipeline.R and compares
# the protein table, the EList E and the DA calls with the unnormalised state. Diagnosis only.
# Usage: R_LIBS=<kit_100 library> F19_OUT=<f19_pipeline.R output dir> Rscript .../f19_norm_check.R

suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
OUT <- Sys.getenv("F19_OUT", file.path(tempdir(), "f19"))
limpa <- readRDS(file.path(OUT, "limpa_route_states.rds"))
pre <- limpa$states[["correlation_filtered"]]
contrast <- as.character(limpa$contrasts_tbl[[1L]][[1L]])
planted <- sprintf("P%05d", 1:40)
invisible(utils::capture.output(normed <- suppressMessages(normaliseBetweenSamples(pre, normalisation_method = "cyclicloess"))))
samples <- setdiff(names(pre@protein_quant_table), pre@protein_id_column)
cat(sprintf("protein table changed by cyclicloess: max abs diff %.3f\n",
    max(abs(as.matrix(normed@protein_quant_table[samples]) - as.matrix(pre@protein_quant_table[samples])), na.rm = TRUE)))
cat(sprintf("EList E changed: %s\n", !identical(normed@args$limpa_dpc_quant_results$quantified_elist$E,
    pre@args$limpa_dpc_quant_results$quantified_elist$E)))
calls <- function(object) {
    single <- data.frame(contrasts = contrast, stringsAsFactors = FALSE)
    arguments <- list(contrasts_tbl = single, formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)", args_row_id = object@protein_id_column,
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper"))
        if (is.null(object@args[[section]])) object@args[[section]] <- arguments
    invisible(utils::capture.output(r <- suppressMessages(suppressWarnings(differentialAbundanceAnalysis(object, contrasts_tbl = single,
        formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, qvalue_column = "fdr_qvalue",
        raw_pvalue_column = "raw_pvalue")))))
    t <- r$contrasts_results_table
    sig <- regmatches(t[[object@protein_id_column]], regexpr("P[0-9]{5}", t[[object@protein_id_column]]))[t$fdr_qvalue < 0.05]
    sprintf("model %s, sig %d, planted %d, FP %d", r$da_model$model, length(sig), sum(sig %in% planted), sum(!sig %in% planted))
}
cat("DA, unnormalised state:          ", calls(pre), "\n")
cat("DA, cyclicloess-normalised state:", calls(normed), "\n")
