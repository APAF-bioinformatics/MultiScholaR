# KIT-125, classic proteomics: null splits of the healthy volunteers of Zenodo 10.5281/zenodo.17607826 (a FragPipe
# reprocessing of the COMBAT blood atlas, PXD023175; CC BY 4.0) through the code paths of KIT-124's simulator S1
# (docs/qa/workflow-execution/kit/kit-124/sim/s1_classic.R). Usage and outputs: k125_common.R.
#
# Data, in place of S1's `simulate`:
#   file    <data dir>/combined_protein.tsv, the "<sample> MaxLFQ Intensity" columns of the samples named HC_*
#           (30 healthy volunteers, one plasma sample each: the biological unit; no technical replicates, one plate)
#   values  0 -> missing; log2. No further between-sample normalisation (MaxLFQ is normalised by FragPipe; S1 has none)
#   pool    the HC samples with at least half the median number of quantified proteins among the 30 (a rule against
#           failed injections, fixed before any run and blind to the splits)
#   split   12 units, 6 and 6 (k125_common.R); false groups "control" (A) and "treated" (B)
#   filter  S1's: a protein is kept when at least half of each false group's samples are observed
# Arms (S1's functions, copied; the numbers are S1's):
#   1  no RUV: runTestsContrasts as the default path calls it ("~ 0 + group", TREAT cutoff 0, trend, robust)
#   3  today's design: findBestNegCtrlPercentage (the GUI's defaults) -> ruvIII_C_Varying -> the same test
#   7  recommended: all proteins as controls, the false groups as replicate sets, k = E =
#      sva::num.sv(complete proteins, ~ group, "be", B = 20, seed = stream 0), computeRuvCovariates on the proteins
#      observed in every sample, W in the model of the same test on all proteins. E = 0: no RUV.
# The package kernels of arm 7 are on branch kit-126 and not merged; this is the simulator's own arm-7 path on the
# pinned tree, so that the result is comparable with S1's.

local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
    source(file.path(dirname(file), "k125_common.R"))
})
args <- k125Args()
suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
try(RUVIIIC::omp_set_num_threads(1L), silent = TRUE)
quietly <- function(expr) invisible(suppressWarnings(suppressMessages(utils::capture.output(expr, type = "output"))))

# ---- S1's pinned test settings (s1_classic.R:158-171) --------------------------------------------------------------
test_settings <- list(formula_string = "~ 0 + group", treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
    q_threshold = 0.05, contrast = "grouptreated-groupcontrol")
local({
    lines <- readLines(file.path("inst", "config", "config.ini"))
    from <- which(trimws(lines) == "[deAnalysisParameters]")
    stopifnot(length(from) == 1L)
    block <- lines[seq(from + 1L, length(lines))]
    value <- function(key) trimws(sub("^[^=]*=", "", grep(paste0("^\\s*", key, "\\s*="), block, value = TRUE)[[1L]]))
    ok <- identical(value("formula_string"), test_settings$formula_string) &&
        as.double(value("treat_lfc_cutoff")) == test_settings$treat_lfc_cutoff &&
        identical(tolower(value("eBayes_trend")), "true") && identical(tolower(value("eBayes_robust")), "true") &&
        as.double(value("de_q_val_thresh")) == test_settings$q_threshold
    if (!ok) stop("inst/config/config.ini [deAnalysisParameters] no longer holds the settings S1 pins")
})

# ---- the public table ----------------------------------------------------------------------------------------------
input <- file.path(args$data_dir, "combined_protein.tsv")
raw <- utils::read.delim(input, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
columns <- grep("^HC_.* MaxLFQ Intensity$", names(raw), value = TRUE)
all_x <- as.matrix(raw[, columns]); storage.mode(all_x) <- "double"
all_x[!is.finite(all_x) | all_x <= 0] <- NA_real_
all_x <- log2(all_x)
dimnames(all_x) <- list(raw[["Protein ID"]], sub(" MaxLFQ Intensity$", "", columns))
stopifnot(!anyDuplicated(rownames(all_x)), ncol(all_x) == 30L)
quantified <- colSums(is.finite(all_x))
pool <- colnames(all_x)[quantified >= stats::median(quantified) / 2]

# ---- S1's functions (s1_classic.R:297-334, 383-428), unchanged but for the data they close over -------------------
make_object <- function(x, design, technical) {
    quant <- as.data.frame(c(list(Protein.Group = rownames(x)), as.data.frame(x)), check.names = FALSE, stringsAsFactors = FALSE)
    object_design <- design; rownames(object_design) <- NULL
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(x), stringsAsFactors = FALSE), design_matrix = object_design,
        sample_id = "Run", group_id = "group", technical_replicate_id = technical,
        args = list(globalParameters = list(workflow_type = "LFQ")))
}
as_matrix <- function(object) {
    q <- object@protein_quant_table
    m <- as.matrix(q[, setdiff(names(q), object@protein_id_column), drop = FALSE]); rownames(m) <- q[[object@protein_id_column]]; m
}
refuse <- function(message) stop(structure(class = c("s1_refusal", "error", "condition"), list(message = message, call = NULL)))
is_refusal <- function(e) inherits(e, c("s1_refusal", "prot_ruv_insufficient_controls", "prot_ruv_empty_result")) ||
    grepl("No valid percentage found", conditionMessage(e), fixed = TRUE)
da <- function(matrix, design, covariates = NULL, terms = character()) {
    if (!is.null(covariates)) {
        design <- cbind(design, as.data.frame(covariates[rownames(design), , drop = FALSE]))
        terms <- c(terms, colnames(covariates))
    }
    formula <- paste(c(test_settings$formula_string, terms), collapse = " + ")
    value <- NULL
    quietly(value <- runTestsContrasts(matrix, contrast_strings = test_settings$contrast, design_matrix = design,
        formula_string = formula, weights = NA, treat_lfc_cutoff = test_settings$treat_lfc_cutoff,
        eBayes_trend = test_settings$eBayes_trend, eBayes_robust = test_settings$eBayes_robust,
        group_id = "group", technical_replicate_id = "tech_rep_group"))
    tab <- value$results[[1L]]
    list(p = tab$raw_pvalue, q = tab$fdr_qvalue)
}

prepare <- function(units, seed) {
    samples <- c(units$A, units$B)
    groups <- rep(c("control", "treated"), each = 6L)
    x <- all_x[, samples, drop = FALSE]
    keep <- rowSums(is.finite(x[, groups == "control"])) >= 3 & rowSums(is.finite(x[, groups == "treated"])) >= 3
    x <- x[keep, , drop = FALSE]
    runs <- sprintf("run_%03d", seq_along(samples)); colnames(x) <- runs
    design <- data.frame(Run = runs, group = factor(groups, levels = c("control", "treated")), replicates = seq_along(runs),
        tech_rep_group = paste0("u", seq_along(runs)), no_tech = paste0("i", seq_along(runs)), stringsAsFactors = FALSE)
    rownames(design) <- runs
    complete <- rowSums(!is.finite(x)) == 0
    xc <- x[complete, , drop = FALSE]
    n_samples <- ncol(x); n_groups <- 2L
    all_controls <- stats::setNames(rep(TRUE, nrow(x)), rownames(x))
    stream <- function(arm) seed * 100L + arm
    E <- function() {
        if (nrow(xc) < 2L) refuse("no proteins observed in every sample: E is undefined")
        value <- NULL
        quietly(value <- sva::num.sv(xc, stats::model.matrix(~ group, design), method = "be", B = 20, seed = stream(0L)))
        as.integer(value)
    }
    no_ruv <- list(same_as_no_ruv = TRUE)
    in_model <- function(controls, k, technical) {
        if (k == 0L) return(no_ruv)
        if (n_samples - n_groups - k < 1L) refuse(sprintf("no residual degree of freedom for %d covariates with %d samples and %d groups", k, n_samples, n_groups))
        usable <- unname(controls[complete])
        if (sum(usable) <= k) refuse(sprintf("k = %d is not below the %d controls observed in every sample", k, sum(usable)))
        blocks <- length(unique(as.character(design$group)))
        if (k > n_samples - blocks) refuse(sprintf("k = %d exceeds samples minus replicate sets (%d)", k, n_samples - blocks))
        value <- NULL
        quietly(value <- computeRuvCovariates(make_object(xc, design, technical), k = k, controls = usable, ruv_grouping_variable = "group"))
        c(da(x, design, value$W), list(k = k))
    }
    list(arms = list(
        `1` = function() c(da(x, design), list(k = 0L)),
        `3` = function() {
            value <- NULL
            quietly(value <- findBestNegCtrlPercentage(normalised_protein_matrix_obj = make_object(x, design, "tech_rep_group"),
                percentage_range = seq(1, 20, by = 1), separation_metric = "max_difference", k_penalty_weight = 0.5,
                adaptive_k_penalty = TRUE, ruv_grouping_variable = "group", verbose = FALSE))
            controls <- as.logical(value$best_control_genes_index); controls[is.na(controls)] <- FALSE
            k <- validateProtNormRuvK(value$best_k, "ruv_k")
            corrected <- NULL
            quietly(corrected <- ruvIII_C_Varying(make_object(x, design, "tech_rep_group"), ruv_grouping_variable = "group",
                ruv_number_k = k, ctrl = stats::setNames(controls, rownames(x))))
            c(da(as_matrix(corrected), design), list(k = k))
        },
        `7` = function() in_model(all_controls, E(), "no_tech")))
}

k125Run("classic", args, pool = pool, batch = NULL, types = "random",
    designs = c(`1` = "no RUV", `3` = "today's design", `7` = "recommended"), prepare = prepare, is_refusal = is_refusal,
    inputs = input, packages = c("MultiScholaR", "limma", "ruv", "RUVIIIC", "sva", "qvalue", "pkgload"),
    notes = sprintf("pool rule: HC samples with at least %.1f quantified proteins (half the median of the 30); quantified proteins per pooled sample: %d to %d",
        stats::median(quantified) / 2, min(quantified[pool]), max(quantified[pool])))
