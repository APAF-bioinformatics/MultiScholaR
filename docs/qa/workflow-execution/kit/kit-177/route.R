# KIT-177 shared route runner (sourced by effect.R, replay-check.R and t1b.R; never by the package).
#
# The label-free and TMT default path as it ships after KIT-180 (VAL-26), run in process as declared operations, as the
# tabs run them: the design state (made here, either transform, or by the package's own design step), Accession Cleanup
# (no FASTA table: recorded as skipped), Protein Intensity Filter, Duplicate Removal, Protein Replicate Filter (each
# compute half, then the shared protein finishing step without publication), normalisation with RUV (the tab's compute
# and commit halves), the sample correlation filter, and the DA (`runDIAReferenceDaAnalysis`, publication deferred).
# Source into an environment whose parent is the MultiScholaR namespace (or the global environment after load_all).

k177Config <- function() {
    path <- system.file("config", "config.ini", package = "MultiScholaR")
    utils::capture.output(config <- suppressMessages(readConfigFile(file = path, create_workers = FALSE)))
    config
}

#' A fresh on-disk owner of one workflow type, with experiment paths under its root.
k177NewManager <- function(root, workflow_type) {
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, paste0("k177-", basename(root)), workflow_type),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(workflow_type)
    paths <- list(source_dir = file.path(root, "src"), protein_qc_dir = file.path(root, "pqc"),
        da_output_dir = file.path(root, "da"), publication_graphs_dir = file.path(root, "pub"),
        time_dir = file.path(root, "time"), results_dir = file.path(root, "res"))
    for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
    manager$.context$experiment_paths <- paths
    manager
}

k177Close <- function(manager) closeDIAStateAdapter(manager)

#' The design step's transform applied to a wide table, in the script (KIT-177 T5's construction).
#' `old`: log2(x + 1) on every cell (the base's step). `new`: the rule of section 3 (0 and NaN missing and counted,
#' then log2). Returns `list(wide, fixed)`.
k177Transform <- function(wide, id_column, transform) {
    runs <- setdiff(names(wide), id_column)
    if (identical(transform, "old")) {
        for (run in runs) wide[[run]] <- log2(wide[[run]] + 1)
        return(list(wide = wide, fixed = list(pivot = "protein x run", transform = "log2(x + 1)")))
    }
    zeros <- 0L
    nans <- 0L
    for (run in runs) {
        x <- as.numeric(wide[[run]])
        zero <- !is.na(x) & x == 0
        nan <- is.nan(x)
        zeros <- zeros + sum(zero)
        nans <- nans + sum(nan)
        x[zero | nan] <- NA_real_
        wide[[run]] <- log2(x)
    }
    list(wide = wide, fixed = list(pivot = "protein x run", transform = "log2(x); 0 missing",
        nonpositive_set_to_na = zeros, nonfinite_set_to_na = nans))
}

#' Commit `protein_s4_initial` made in the script from a wide intensity table and a design (Run, group, replicates).
k177DesignState <- function(manager, wide, id_column, design, transform, config) {
    made <- k177Transform(wide, id_column, transform)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    s4 <- ProteinQuantitativeData(protein_quant_table = made$wide, protein_id_column = id_column,
        protein_id_table = dplyr::distinct(made$wide, !!rlang::sym(id_column)), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group", args = config)
    setDIAStepProvenance(manager, "design_builder", fixed = made$fixed)
    manager$saveState("protein_s4_initial", s4, config, "Initial Protein S4 object created in the KIT-177 script")
    invisible(made$fixed)
}

#' The normalisation tab's declared parameters (the proteomics tab's rendered defaults, `omicGuiConfig("proteomics")`).
k177NormParameters <- function(norm_method = "cyclicloess", ruv_mode = "automatic", ruv_percentage = 5, ruv_k = NULL) {
    list(norm_method = norm_method, ruv_mode = ruv_mode, auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE,
        ruv_percentage = ruv_percentage, ruv_k = ruv_k, grouping = "group")
}

#' The path from the design state: QC (the tabs' defaults unless given), normalisation, correlation filter, DA.
#' @param qc `list(min_reps_per_group, min_groups, percentile)` for the intensity filter, or `intensity = FALSE` to
#'   leave the filter out (T2).
#' @param norm `k177NormParameters()` output.
#' @param contrast `list(name, full)` e.g. `list(name = "treated.vs.control", full = "treated.vs.control=grouptreated-groupcontrol")`.
#' @return `list(states, da, k, percentage, control_count, ruv_applied, columns)`.
k177RunPath <- function(manager, norm, contrast, qc = list(min_reps_per_group = 2, min_groups = 2, percentile = 1),
    intensity = TRUE, correlation = 0.5, treat = 0) {
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    paths <- manager$.context$experiment_paths
    withDIAStateScope(manager, function() {
        finish <- function(value) finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
        finish(computeDIAReferenceAccessionCleanup(manager, getDIAStateReference(manager),
            list(delimiter = ";", aggregation_method = "mean")))
        if (isTRUE(intensity)) finish(computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager),
            list(strict_mode = FALSE, min_reps_per_group = qc$min_reps_per_group, min_groups = qc$min_groups,
                proteins_intensity_cutoff_percentile = qc$percentile)))
        finish(computeDIAReferenceProteinDuplicateRemoval(manager, getDIAStateReference(manager),
            list(aggregation_method = "mean")))
        completeDIAProteinAction(workflow, computeDIAReferenceProteinReplicateFilter(manager, getDIAStateReference(manager),
            list(grouping_variable = "group", parallel_cores = 4L), variant = "rows"), publish = FALSE)
        run <- computeDIAReferenceNormalization(manager, getDIAStateReference(manager), norm)
        commit <- diaNormalizationCommit(manager, norm, run)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        saved <- run$value
        correlated <- computeDIAReferenceCorrelation(manager, saved$final, "group", correlation)
        commit <- diaCorrelationCommit(manager, correlated$final, correlation, FALSE)
        if (!is.null(commit)) commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
        contrasts <- data.frame(contrasts = contrast$name, full_format = contrast$full, friendly_names = contrast$name,
            stringsAsFactors = FALSE)
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- getDIAStateReference(manager)
        da_data$contrasts_available <- contrasts$contrasts
        workflow$contrasts_tbl <- contrasts
        manifest <- runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = treat), da_data, workflow, paths, contrasts = contrasts, publish = FALSE)
        outputs <- readDIAContrastOutput(manifest$results[[1L]], c("da_proteins_long", "contrasts_results_table"))
        controls <- saved$controls
        list(states = as.character(unlist(manager$.context$snapshot$state_history)),
            da = as.data.frame(outputs$da_proteins_long), results = as.data.frame(outputs$contrasts_results_table %||% data.frame()),
            k = saved$best_k %||% NA, percentage = saved$percentage %||% NA,
            control_count = if (is.logical(controls)) sum(controls, na.rm = TRUE) else length(controls[!is.na(controls)]),
            ruv_applied = saved$ruv_result$ruv_applied %||% NA)
    })
}

#' Ticket 9.1's generator, fixed before the run. `set.seed(seed)` with R 4.6.1's default kinds; draws in this order:
#' z (2,000 standard normals); the changed set `sample(2000, 200)` (its first 100 +1 log2, the rest -1 log2, applied to
#' group B's mean); then column by column (runs 1 to 5 group A, 6 to 10 group B) e1 ~ N(0, 0.15^2) for the 2,000
#' proteins, then e2 ~ N(0, 3^2) for the 2,000 proteins; value = mean x exp(e1) + e2; a value at or below 0 is missing.
k177Simulate91 <- function(seed, m) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    n <- 2000L
    z <- stats::rnorm(n)
    mu <- m * 2^(1.5 * z)
    changed <- sample(n, 200L)
    sign <- integer(n)
    sign[changed[1:100]] <- 1L
    sign[changed[101:200]] <- -1L
    mu_b <- mu * 2^sign
    values <- matrix(NA_real_, n, 10L)
    for (j in seq_len(10L)) {
        mean_j <- if (j <= 5L) mu else mu_b
        e1 <- stats::rnorm(n, 0, 0.15)
        e2 <- stats::rnorm(n, 0, 3)
        values[, j] <- mean_j * exp(e1) + e2
    }
    values[!is.na(values) & values <= 0] <- NA_real_
    runs <- sprintf("run%02d", seq_len(10L))
    colnames(values) <- runs
    wide <- tibble::as_tibble(cbind(data.frame(Protein.Ids = sprintf("P%05d", seq_len(n)), stringsAsFactors = FALSE),
        as.data.frame(values)))
    cuts <- stats::quantile(mu, c(1 / 3, 2 / 3), names = FALSE)
    list(wide = wide, sign = sign, mu = mu, third = as.integer(cut(mu, c(-Inf, cuts, Inf), labels = FALSE)),
        design = data.frame(Run = runs, group = rep(c("control", "treated"), each = 5L), replicates = rep(1:5, 2L),
            stringsAsFactors = FALSE))
}

#' Ticket 9.1's measures of one route run against the planted truth (BH on the reported p over the tested proteins).
k177Measure91 <- function(table, sim) {
    index <- as.integer(sub("^P", "", table$id))
    tested <- !is.na(table$p)
    q <- rep(NA_real_, nrow(table))
    q[tested] <- stats::p.adjust(table$p[tested], "BH")
    called <- tested & q < 0.05
    sign <- sim$sign[index]
    third <- sim$third[index]
    changed <- tested & sign != 0L
    share <- function(a, b) if (b > 0) a / b else NA_real_
    out <- list(tested = sum(tested), changed_tested = sum(changed),
        bias = if (any(changed)) mean(sign[changed] * table$logFC[changed]) else NA_real_,
        bias_up = if (any(changed & sign > 0L)) mean(table$logFC[changed & sign > 0L]) else NA_real_,
        bias_down = if (any(changed & sign < 0L)) mean(-table$logFC[changed & sign < 0L]) else NA_real_,
        calls = sum(called), fdp = share(sum(called & sign == 0L), sum(called)),
        tpr = sum(called & sign != 0L) / sum(sim$sign != 0L))
    for (t in 1:3) {
        out[[paste0("fdp_third", t)]] <- share(sum(called & sign == 0L & third == t), sum(called & third == t))
        out[[paste0("tpr_third", t)]] <- share(sum(called & sign != 0L & third == t), sum(sim$sign != 0L & sim$third == t))
        out[[paste0("bias_third", t)]] <- if (any(changed & third == t)) mean(sign[changed & third == t] *
            table$logFC[changed & third == t]) else NA_real_
    }
    out$calls_ids <- table$id[called]
    out
}

#' The DA table's columns the measurements need: protein id, logFC, t (when the kernel reports it) and the raw p.
k177DaColumns <- function(result, id_column) {
    pick <- function(table, names) {
        hit <- names[names %in% names(table)]
        if (length(hit)) table[[hit[[1L]]]] else NULL
    }
    da <- result$da
    ids <- as.character(pick(da, c(id_column, "uniprot_acc", "Protein.Ids")))
    logfc <- as.numeric(pick(da, c("log2FC", "logFC")))
    p <- as.numeric(pick(da, c("raw_pvalue", "P.Value")))
    q <- pick(da, c("fdr_qvalue"))
    t <- pick(da, c("t", "t.mod"))
    if (is.null(t) && nrow(result$results)) {
        rt <- result$results
        rid <- as.character(pick(rt, c(id_column, "uniprot_acc", "Protein.Ids")))
        tv <- pick(rt, c("t", "t.mod"))
        if (!is.null(tv) && length(rid)) t <- as.numeric(tv)[match(ids, rid)]
    }
    data.frame(id = ids, logFC = logfc, t = if (is.null(t)) NA_real_ else as.numeric(t), p = p,
        q = if (is.null(q)) NA_real_ else as.numeric(q), stringsAsFactors = FALSE)
}

#' Fixture F177 (ticket section 5; the same construction as tests/testthat/test-prot-design-log2.R `k177F177`), with
#' T1b's additions when `planted = TRUE`: one unwanted factor balanced across the groups (w = +1, -1 alternating within
#' each group) with per-protein loadings N(0, 0.5^2) on the log2 scale, and 10 % of proteins changed by +1 or -1 log2 in
#' group B (half each). Draw order: F177's draws, then the loadings, then the changed set and its signs.
k177F177 <- function(seed, n = 400L, range = c(0.5, 1e7), planted = FALSE) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    level <- exp(stats::runif(n, log(range[[1L]]), log(range[[2L]])))
    values <- level * 2^matrix(stats::rnorm(n * 10L, 0, 0.25), n, 10L)
    for (i in seq_len(n)) for (g in list(1:5, 6:10)) {
        missing <- sample(g, sample(0:2, 1L, prob = c(0.6, 0.3, 0.1)))
        values[i, missing] <- NA_real_
        observed <- setdiff(g, missing)
        ones <- observed[stats::runif(length(observed)) < 0.02]
        values[i, utils::head(ones, 2L)] <- 1
    }
    sign <- integer(n)
    if (planted) {
        loading <- stats::rnorm(n, 0, 0.5)
        w <- rep(c(1, -1, 1, -1, 1), 2L)
        changed <- sample(n, as.integer(n / 10))
        sign[changed] <- sample(rep(c(1L, -1L), length.out = length(changed)))
        values <- values * 2^(outer(loading, w))
        values[, 6:10] <- values[, 6:10] * 2^sign
    }
    runs <- sprintf("S%02d", 1:10)
    list(ids = sprintf("P%05d", seq_len(n)), runs = runs, values = values, sign = sign,
        design = data.frame(Run = runs, group = rep(c("control", "treated"), each = 5L), replicates = rep(1:5, 2L),
            stringsAsFactors = FALSE))
}

#' F177's wide table for one shape at scale s (the identifier column `Protein.Ids`).
k177F177Wide <- function(sim, s = 1) {
    values <- sim$values * s
    colnames(values) <- sim$runs
    tibble::as_tibble(cbind(data.frame(Protein.Ids = sim$ids, stringsAsFactors = FALSE), as.data.frame(values)))
}
