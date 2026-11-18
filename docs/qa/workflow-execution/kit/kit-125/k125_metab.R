# KIT-125, metabolomics: null splits of the adult volunteers of MetaboLights MTBLS404 (Thevenot et al., "sacurine",
# human urine, LC-HRMS negative mode; EMBL-EBI terms of use) through the code paths of KIT-124's simulator S4
# (kit-124/sim/s4_lipid_metab.R), on the metabolomics class. Usage and outputs: k125_common.R.
#
# Data, in place of S4's `generateData`:
#   files   <data dir>/m_MTBLS404_sacurine_v2_maf.tsv (120 identified metabolites; one column per MS assay),
#           s_MTBLS404.txt (material type per sample), a_MTBLS404_LC-MS_sacurine.txt (MS assay name, batch)
#   pool    the MS assays of material type "biological material" (184 volunteers, one urine sample each: the
#           biological unit; pooled QC injections and blanks are not used); batch = Factor Value[Batch] (ne1, ne2)
#   values  as deposited (the study's processed intensities); 0 or missing -> missing; log2; one assay ("neg").
#           No normalisation and no standards: S4's package arms take the log2 table as it is
#   split   12 units, 6 and 6; types random and batch (k125_common.R); false groups "A" and "B"; ids S01..S12
#   filter  S4's presence rule: a feature is kept when it is observed at least once in each false group
# Arms (S4's functions, copied for one assay; the numbers are S4's):
#   1  no correction: the default path's kernel (limma, trend and robust, no TREAT; q by the package's rule)
#   3  today's design: runPerAssayRuvOptimization (the tab's defaults) -> ruvIII_C_Varying (the groups as replicate
#      sets) -> runMetabolitesDA on the corrected object
#   7  recommended: all features as controls, the false groups as replicate sets, k = E = sva::num.sv(features
#      complete in every sample, ~ group, "be", B = 20, seed * 100 + 99); W from ruv::RUVIII on the complete
#      features, the controls centred (S4's ruvW), in the model of the kernel on the uncorrected values. E = 0: arm 1.
# S4 simulated the lipid class and proved at start-up that the lipid and metabolite kernels give identical output;
# this script calls the metabolite functions and repeats the kernel part of that proof on its first split's data.
# The package kernels of arm 7 are on branch kit-127 and not merged; this is the simulator's own arm-7 path on the
# pinned tree.

if (!identical(Sys.getenv("OMP_NUM_THREADS"), "1")) stop("run with OMP_NUM_THREADS=1 (RUVIIIC is threaded)")
local({
    file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
    source(file.path(dirname(file), "k125_common.R"))
})
args <- k125Args()
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
logger::log_threshold(logger::FATAL)
logger::log_threshold(logger::FATAL, namespace = "MultiScholaR")
for (package in c("limma", "qvalue", "ruv", "RUVIIIC", "sva", "statmod")) loadNamespace(package)

ASSAY <- "neg"
ID_COL <- "database_identifier"
CONTRAST <- "groupB-groupA"
RUV_PARAMS <- list(ruv_grouping_variable = "group", percentage_min = 1, percentage_max = 20,
    separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, manual_k = 2,
    manual_percentage = 5)

# ---- the public table ----------------------------------------------------------------------------------------------
inputs <- file.path(args$data_dir, c("m_MTBLS404_sacurine_v2_maf.tsv", "s_MTBLS404.txt", "a_MTBLS404_LC-MS_sacurine.txt"))
maf <- utils::read.delim(inputs[[1L]], check.names = FALSE, quote = "", stringsAsFactors = FALSE)
samples <- utils::read.delim(inputs[[2L]], check.names = FALSE, stringsAsFactors = FALSE)
assays <- utils::read.delim(inputs[[3L]], check.names = FALSE, stringsAsFactors = FALSE)
assays <- assays[assays[["MS Assay Name"]] %in% names(maf), ]
material <- samples[["Factor Value[Material type]"]][match(assays[["Sample Name"]], samples[["Sample Name"]])]
subjects <- assays[material == "biological material", ]
stopifnot(nrow(subjects) == 184L, !anyDuplicated(subjects[["Sample Name"]]), !anyDuplicated(subjects[["MS Assay Name"]]))
pool <- subjects[["MS Assay Name"]]
batch <- subjects[["Factor Value[Batch]"]]
all_y <- as.matrix(maf[, pool]); storage.mode(all_y) <- "double"
all_y[!is.finite(all_y) | all_y <= 0] <- NA_real_
all_y <- log2(all_y)
rownames(all_y) <- sprintf("%s_f%03d", ASSAY, seq_len(nrow(all_y)))

# ---- S4's helpers (s4_lipid_metab.R:189-218, 349-390, 446-468), for the metabolite class and one assay ---------------
quiet <- function(expr) {
    con <- file(nullfile(), "w")
    sink(con)
    sink(con, type = "message")
    on.exit({
        sink(type = "message")
        sink()
        close(con)
    })
    suppressWarnings(suppressMessages(expr))
}
refusal <- function(reason) stop(structure(class = c("s4_refusal", "error", "condition"), list(message = reason, call = NULL)))
is_refusal <- function(e) inherits(e, c("s4_refusal", "omics_ruv_insufficient_controls"))
assayTable <- function(m) {
    table <- tibble::as_tibble(as.data.frame(m, check.names = FALSE))
    table <- tibble::add_column(table, !!ID_COL := rownames(m), .before = 1L)
    tibble::add_column(table, name = rownames(m), .after = 1L)
}
metabObject <- function(m, design) {
    methods::new("MetaboliteAssayData", metabolite_data = stats::setNames(list(assayTable(m)), ASSAY), metabolite_id_column = ID_COL,
        annotation_id_column = "name", design_matrix = design[, c("Sample_ID", "group")], sample_id = "Sample_ID",
        group_id = "group", args = list())
}
assayMatrix <- function(table, samples) {
    m <- as.matrix(as.data.frame(table[, samples], check.names = FALSE))
    rownames(m) <- table[[ID_COL]]
    m
}
kernelTable <- function(result) {
    t <- result$results[[CONTRAST]]
    data.frame(feature = rownames(t), p = t$P.Value, q = t$fdr_qvalue, stringsAsFactors = FALSE)
}
daKernel <- function(m, design, formula, kernel = runTestsContrastsMetabDA) {
    if (anyNA(m)) {
        x <- stats::model.matrix(stats::as.formula(formula), design)
        seen <- is.finite(m)
        pattern <- drop(seen %*% 2^(seq_len(ncol(m)) - 1))
        estimable <- vapply(split(seq_len(nrow(m)), pattern), function(i) qr(x[seen[i[[1L]], ], , drop = FALSE])$rank == ncol(x), logical(1L))
        m <- m[estimable[as.character(pattern)], , drop = FALSE]
    }
    kernelTable(quiet(kernel(data = m, contrast_strings = CONTRAST, design_matrix = design,
        formula_string = formula, sample_id_col = "Sample_ID", treat_lfc_cutoff = NA, eBayes_trend = TRUE,
        eBayes_robust = TRUE)))
}
groupM <- function(design) {
    m <- ruv::replicate.matrix(stats::setNames(design$group, design$Sample_ID))
    rownames(m) <- design$Sample_ID
    m
}
# S4's ruvW: W of RUV-III on complete columns, the controls centred, as computeRuvCovariates forms it.
ruvW <- function(y, m, ctl, k, residual_df) {
    if (is.na(k)) refusal("k could not be estimated")
    if (k > sum(ctl)) refusal(sprintf("k = %d above %d complete controls", k, sum(ctl)))
    if (k > nrow(y) - ncol(m)) refusal(sprintf("k = %d above the rank %d the replicates give", k, nrow(y) - ncol(m)))
    if (residual_df - k < 1L) refusal(sprintf("no residual degrees of freedom with k = %d", k))
    ac <- ruv::RUVIII(Y = y, M = m[rownames(y), , drop = FALSE], ctl = ctl, k = k, return.info = TRUE)$fullalpha[seq_len(k), ctl, drop = FALSE]
    w <- sweep(y[, ctl, drop = FALSE], 2L, colMeans(y[, ctl, drop = FALSE])) %*% t(ac) %*% solve(ac %*% t(ac))
    rownames(w) <- rownames(y)
    w <- scale(w, center = TRUE, scale = FALSE)
    colnames(w) <- paste0("W", seq_len(k))
    w
}

kernels_checked <- FALSE
prepare <- function(units, seed) {
    ids <- sprintf("S%02d", 1:12)
    design <- data.frame(Sample_ID = ids, group = rep(c("A", "B"), each = 6L), stringsAsFactors = FALSE)
    y <- all_y[, c(units$A, units$B), drop = FALSE]
    colnames(y) <- ids
    seen <- is.finite(y)
    y <- y[rowSums(seen[, 1:6]) >= 1L & rowSums(seen[, 7:12]) >= 1L, , drop = FALSE]
    complete <- rowSums(!is.finite(y)) == 0L
    if (!kernels_checked) {
        # the lipid and metabolite kernels give the same table on this data (S4's start-up check, the kernel part)
        stopifnot(identical(daKernel(y, design, "~ 0 + group"), daKernel(y, design, "~ 0 + group", runTestsContrastsLipidDA)))
        kernels_checked <<- TRUE
    }
    list(arms = list(
        `1` = function() c(as.list(daKernel(y, design, "~ 0 + group")[c("p", "q")]), list(k = 0L)),
        `3` = function() {
            object <- metabObject(y, design)
            sweep <- quiet(runPerAssayRuvOptimization(object, "automatic", RUV_PARAMS))
            if (!isTRUE(sweep[[ASSAY]]$success)) refusal(sweep[[ASSAY]]$error %||% "the sweep found no valid percentage")
            best_k <- extractBestKPerAssay(sweep)
            corrected <- quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = best_k, ctrl = extractCtrlPerAssay(sweep)))
            if (identical(corrected@metabolite_data[[ASSAY]], object@metabolite_data[[ASSAY]])) stop("ruvIII_C_Varying left the assay uncorrected")
            da <- quiet(runMetabolitesDA(corrected, data.frame(contrasts = CONTRAST, friendly_names = CONTRAST, stringsAsFactors = FALSE),
                "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_metabolites_long
            rows <- da[da$assay == ASSAY, ]
            if (!nrow(rows)) stop("runMetabolitesDA returned no rows")
            list(p = rows$P.Value, q = rows$fdr_qvalue, k = as.integer(best_k[[ASSAY]]))
        },
        `7` = function() {
            dat <- y[complete, , drop = FALSE]
            E <- if (nrow(dat) < 20L) NA_integer_ else tryCatch(as.integer(quiet(sva::num.sv(dat, stats::model.matrix(~ factor(design$group)),
                method = "be", B = 20, seed = seed * 100L + 99L))), error = function(e) NA_integer_)
            if (!is.na(E) && E == 0L) return(list(same_as_no_ruv = TRUE))
            yc <- t(dat)
            w <- ruvW(yc, groupM(design), rep(TRUE, ncol(yc)), E, residual_df = nrow(design) - 2L)
            table <- daKernel(y, cbind(design, as.data.frame(w[design$Sample_ID, , drop = FALSE])),
                paste("~ 0 + group +", paste(colnames(w), collapse = " + ")))
            list(p = table$p, q = table$q, k = E)
        }))
}

k125Run("metab", args, pool = pool, batch = batch, types = c("random", "batch"),
    designs = c(`1` = "no RUV", `3` = "today's design", `7` = "recommended"), prepare = prepare, is_refusal = is_refusal,
    inputs = inputs, packages = c("MultiScholaR", "limma", "ruv", "RUVIIIC", "sva", "qvalue", "pkgload"),
    notes = sprintf("features in the table: %d; missing share among the pool's values: %.4f", nrow(all_y), mean(is.na(all_y))))
