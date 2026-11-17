# KIT-125, DIA_LIMPA: null splits of the untreated single HeLa cells of PRIDE PXD024043 (Brunner et al. 2022; CC0)
# through the code paths of KIT-124's simulator S2 (kit-124/sim/s2_limpa.R). Usage and outputs: k125_common.R.
#
# Data, in place of S2's `simulateData`:
#   file    <data dir>/DIANN1.8_SingleCells_CellCycle.zip, member 20210919_DIANN_SingleCellOutput.tsv (the DIA-NN 1.8
#           main report). Read once; the precursor matrix of the pool is kept in <data dir>/prepared/ (outside the
#           repository) and reused.
#   rows    Q.Value, PG.Q.Value and, where the columns exist, Global.Q.Value and Global.PG.Q.Value at most 0.01 (the
#           shipped thresholds); Precursor.Normalised above 0, log2; protein = Protein.Group, precursor = Precursor.Id
#   pool    the runs whose name carries the label UB (the paper's untreated cells, one cell per run: the biological
#           unit) with at least 500 precursors (the shipped per-sample cutoff); batch = the acquisition date in the
#           run name (the paper: two independent cell culture batches)
#   split   12 units, 6 and 6; types random and batch (k125_common.R); false groups "A" and "B"
#   filter  precursors observed in at least 2 runs of each false group (the shipped peptide filter's counts:
#           min_reps_per_group = 2, min_groups = 2)
# Per split, S2's `quantify`: proteinMissingValueImputationLimpa on the split's 12 runs (engine fastlimpadpc, the
# shipped engine; dpc_slope 0.8, chunk 1000), no protein normalisation (shipped: none). Then the arms (S2's
# functions, copied; the numbers are S2's):
#   1  no RUV: the route's test (limpa::dpcDE, convertDpcDAToStandardFormat, q by the package's rule)
#   3  today's correction mode: findBestNegCtrlPercentageCapped -> computeRuvCovariates -> test on the corrected EList
#   4  today's covariates mode (the route's shipped default): the same sweep and factors, W1..Wk in the model
#   7  recommended: all proteins as controls, the false groups as replicate sets, k = E = sva::num.sv(quantified
#      matrix, ~ group, "be", B = 20, seed = 100 * seed + 99), covariates. E = 0: no RUV.
# The package kernel of arm 7 is being written on branch kit-129-k-dia; this is the simulator's own arm-7 path on
# the pinned tree.

local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
    source(file.path(dirname(file), "k125_common.R"))
})
args <- k125Args()
suppressMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
for (package in c("limpa", "fastlimpadpc", "limma", "ruv", "sva", "qvalue")) loadNamespace(package)

FORMULA <- "~ 0 + group"
CONTRAST <- "B_vs_A=groupB-groupA"
REFUSALS <- "No valid percentage found|needs more than [0-9]+ control proteins|value is invalid"
MEMBER <- "20210919_DIANN_SingleCellOutput.tsv"

# ---- the route's settings, read from the package and pinned (s2_limpa.R:190-207) ----------------------------------
shipped <- readWorkflowShippedDefaults()
RUV_PARAMETERS <- resolveDIARuvCovariatesParameters(list(), grouping = "group", ruv_mode = "covariates")
QUANT <- list(engine = resolveLimpaDpcEngine(shipped$limpaDpcQuant$engine),
    dpc_slope = as.numeric(shipped$limpaDpcQuant$dpc_slope), chunk = as.numeric(shipped$limpaDpcQuant$chunk))
stopifnot(
    identical(RUV_PARAMETERS[c("percentage_min", "percentage_max", "num_components_to_impute",
        "ruv_grouping_variable", "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight",
        "adaptive_k_penalty")], list(percentage_min = 1, percentage_max = 30, num_components_to_impute = 5,
        ruv_grouping_variable = "group", ruv_qval_cutoff = 0.05, ruv_fdr_method = "BH",
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE)),
    identical(QUANT, list(engine = "fastlimpadpc", dpc_slope = 0.8, chunk = 1000)),
    identical(trimws(shipped$deAnalysisParameters$formula_string), FORMULA),
    as.numeric(shipped$deAnalysisParameters$treat_lfc_cutoff) == 0,
    as.numeric(shipped$deAnalysisParameters$de_q_val_thresh) == 0.05,
    tolower(shipped$deAnalysisParameters$eBayes_trend) == "true",
    tolower(shipped$deAnalysisParameters$eBayes_robust) == "true")
quietly <- function(expr) {
    suppressWarnings(suppressMessages(utils::capture.output(value <- expr, type = "output")))
    invisible(NULL)
}

# ---- the public report: the pool's precursor matrix, built once -----------------------------------------------------
input <- file.path(args$data_dir, "DIANN1.8_SingleCells_CellCycle.zip")
prepared_file <- file.path(args$data_dir, "prepared", "ub_precursors.rds")
if (!file.exists(prepared_file)) {
    dir.create(dirname(prepared_file), showWarnings = FALSE)
    extracted <- file.path(args$data_dir, "extracted")
    if (!file.exists(file.path(extracted, MEMBER))) utils::unzip(input, files = MEMBER, exdir = extracted, junkpaths = TRUE, unzip = "unzip")
    header <- names(data.table::fread(file.path(extracted, MEMBER), nrows = 0L))
    q_columns <- intersect(c("Q.Value", "PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value"), header)
    report <- data.table::fread(file.path(extracted, MEMBER), select = c("Run", "Protein.Group", "Stripped.Sequence", "Precursor.Id",
        "Precursor.Normalised", q_columns), showProgress = FALSE)
    runs_in_report <- length(unique(report$Run))
    for (column in q_columns) report <- report[report[[column]] <= 0.01]
    report <- report[is.finite(Precursor.Normalised) & Precursor.Normalised > 0 & grepl("_UB_", Run, fixed = TRUE)]
    stopifnot(!anyDuplicated(report[, c("Run", "Precursor.Id")]))
    wide <- data.table::dcast(report, Protein.Group + Stripped.Sequence + Precursor.Id ~ Run, value.var = "Precursor.Normalised")
    y <- log2(as.matrix(wide[, -(1:3)]))
    saveRDS(list(y = y, protein = wide$Protein.Group, sequence = wide$Stripped.Sequence, precursor = wide$Precursor.Id,
        q_columns = q_columns, runs_in_report = runs_in_report, input_sha256 = unname(tools::sha256sum(input))), prepared_file)
    rm(report, wide, y)
}
prepared <- readRDS(prepared_file)
stopifnot(identical(prepared$input_sha256, unname(tools::sha256sum(input))))
per_run <- colSums(is.finite(prepared$y))
pool <- colnames(prepared$y)[per_run >= 500L]
batch <- regmatches(pool, regexpr("20[0-9]{6}", pool))
stopifnot(length(batch) == length(pool), length(unique(batch)) == 2L)

# ---- S2's functions (s2_limpa.R:311-361, 481-496), unchanged but for the data they close over ------------------------
quantify <- function(y, protein, sequence, precursor, design) {
    dimnames(y) <- list(paste(protein, precursor, sep = "%"), design$Run)
    object <- methods::new("PeptideQuantitativeData")
    object@peptide_data <- data.frame(Protein.Ids = protein, Stripped.Sequence = precursor, stringsAsFactors = FALSE)
    object@peptide_matrix <- y
    object@protein_id_column <- "Protein.Ids"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@is_logged_data <- TRUE
    object@design_matrix <- design
    object@sample_id <- "Run"
    object@group_id <- "group"
    object@technical_replicate_id <- "tech_rep_group"
    object@args <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))
    counts <- tibble::tibble(Protein.Ids = unique(protein),
        peptide_count = as.integer(tapply(sequence, factor(protein, levels = unique(protein)), function(s) length(unique(s)))),
        peptidoform_count = as.integer(table(factor(protein, levels = unique(protein)))))
    quantified <- NULL
    quietly(quantified <- proteinMissingValueImputationLimpa(object, dpc_slope = QUANT$dpc_slope,
        verbose = FALSE, chunk = QUANT$chunk, engine = QUANT$engine, peptide_summary = counts))
    quantified
}
routeTest <- function(elist, design, covariates = NULL) {
    formula_used <- FORMULA
    if (!is.null(covariates)) {
        design <- cbind(design[setdiff(names(design), colnames(covariates))],
            as.data.frame(covariates)[rownames(design), , drop = FALSE])
        formula_used <- paste(FORMULA, "+", paste(colnames(covariates), collapse = " + "))
    }
    model <- stats::model.matrix(stats::as.formula(formula_used), design)
    samples <- rownames(model)
    replicate <- resolveProtTechnicalReplicateBlock(design[samples, , drop = FALSE], "group", "tech_rep_group")
    block <- if (isTRUE(replicate$technical)) replicate$block else NULL
    fit <- NULL
    result <- NULL
    quietly(fit <- limpa::dpcDE(elist[, samples], model, plot = FALSE, block = block))
    quietly(result <- convertDpcDAToStandardFormat(dpc_fit = fit, contrast_strings = CONTRAST,
        design_matrix = model, eBayes_trend = TRUE, eBayes_robust = TRUE))
    table <- result$results[[1L]]
    list(p = table$raw_pvalue, q = table$fdr_qvalue)
}

prepare <- function(units, seed) {
    samples <- c(units$A, units$B)
    y <- prepared$y[, samples, drop = FALSE]
    seen <- is.finite(y)
    keep <- rowSums(seen[, 1:6]) >= 2L & rowSums(seen[, 7:12]) >= 2L
    runs <- sprintf("run_%03d", 1:12)
    design <- data.frame(Run = runs, group = rep(c("A", "B"), each = 6L), replicates = rep(1:6, 2L), stringsAsFactors = FALSE)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    quantified <- quantify(y[keep, , drop = FALSE], prepared$protein[keep], prepared$sequence[keep], prepared$precursor[keep], design)
    elist <- quantified@args$limpa_dpc_quant_results$quantified_elist
    stopifnot(!anyNA(elist$E), identical(colnames(elist$E), runs),
        identical(rownames(elist$E), quantified@protein_quant_table[["Protein.Ids"]]))
    rownames(design) <- design$Run
    k_e <- NULL
    quietly(k_e <- sva::num.sv(elist$E, stats::model.matrix(~ group, design), method = "be", B = 20, seed = 100L * seed + 99L))
    k_e <- as.integer(k_e)
    sweep <- NULL
    today <- function() {
        if (is.null(sweep)) {
            set.seed(100L * seed + 3L)                      # the sweep runs on arm 3's stream, as in S2
            sweep <<- tryCatch({
                value <- NULL
                quietly(value <- findBestNegCtrlPercentageCapped(quantified,
                    percentage_range = seq(RUV_PARAMETERS$percentage_min, RUV_PARAMETERS$percentage_max, by = 1),
                    num_components_to_impute = RUV_PARAMETERS$num_components_to_impute,
                    ruv_grouping_variable = RUV_PARAMETERS$ruv_grouping_variable,
                    ruv_qval_cutoff = RUV_PARAMETERS$ruv_qval_cutoff, ruv_fdr_method = RUV_PARAMETERS$ruv_fdr_method,
                    separation_metric = RUV_PARAMETERS$separation_metric,
                    k_penalty_weight = RUV_PARAMETERS$k_penalty_weight,
                    adaptive_k_penalty = isTRUE(RUV_PARAMETERS$adaptive_k_penalty), verbose = FALSE))
                value
            }, error = function(e) e)
        }
        if (inherits(sweep, "condition")) stop(sweep)
        sweep
    }
    factorsOf <- function(k, controls) {
        value <- NULL
        quietly(value <- computeRuvCovariates(quantified, k, controls, RUV_PARAMETERS$ruv_grouping_variable))
        value
    }
    inModel <- function(k, controls) {
        factors <- factorsOf(k, controls)
        c(routeTest(elist, design, factors$W), list(k = factors$k))
    }
    subtracted <- function(k, controls) {
        factors <- factorsOf(k, controls)
        corrected <- elist
        corrected$E <- factors$adjusted[rownames(elist$E), colnames(elist$E), drop = FALSE]
        c(routeTest(corrected, design), list(k = factors$k))
    }
    list(arms = list(
        `1` = function() c(routeTest(elist, design), list(k = 0L)),
        `3` = function() subtracted(today()$best_k, today()$best_control_genes_index),
        `4` = function() inModel(today()$best_k, today()$best_control_genes_index),
        `7` = function() if (k_e > 0L) inModel(k_e, rep(TRUE, nrow(elist$E))) else list(same_as_no_ruv = TRUE)))
}

k125Run("dia_limpa", args, pool = pool, batch = batch, types = c("random", "batch"),
    designs = c(`1` = "no RUV", `3` = "today's correction mode", `4` = "today's covariates mode (default)", `7` = "recommended"),
    prepare = prepare, is_refusal = function(e) grepl(REFUSALS, conditionMessage(e)), inputs = input,
    packages = c("MultiScholaR", "limpa", "fastlimpadpc", "limma", "ruv", "sva", "qvalue", "statmod", "pkgload"),
    notes = c(sprintf("report: %d runs; q-value columns filtered at 0.01: %s; UB runs in the matrix: %d; precursors per UB run: %d to %d (median %.0f); pool rule: at least 500",
        prepared$runs_in_report, paste(prepared$q_columns, collapse = ", "), ncol(prepared$y), min(per_run), max(per_run), stats::median(per_run)),
        sprintf("DPC engine: %s (RAYON_NUM_THREADS set to 1 by the package inside the method)", QUANT$engine)))
