# KIT-177 R0, ticket section 9.1: the effect of the offset on level-dependent synthetic data against a planted truth,
# both transforms, before the fix (descriptive, no pass rule). Also the five e2e fixtures through the same path.
#
# Usage, from the package root on a run-tests runner (the commit under test installed; its R/ is the base's):
#   Rscript docs/qa/workflow-execution/kit/kit-177/effect.R [<out dir>] [<ceiling minutes>]
#   <out dir>          default docs/qa/workflow-execution/kit/kit-177/effect
#   <ceiling minutes>  default 187 (ticket section 11, R0: the measurement's ceiling)
#
# Generator and measures: route.R `k177Simulate91`, `k177Measure91` (fixed in the ticket, 9.1). m in {10, 40, 150, 1000,
# 10000}; 20 seeds 20261123 to 20261142. Per seed, m and configuration, the design state is made in this script both
# ways (log2(x + 1), as the base's step makes it; log2(x) with 0 missing, the head's rule) and run through the
# label-free default path (route.R `k177RunPath`): cyclicloess with RUV automatic (configuration 1, the shipped
# defaults) or RUV skip (configuration 2), correlation filter 0.5, DA ~ 0 + group, TREAT lfc 0; calls at BH q < 0.05 on
# the reported p. The runs go to four worker processes (PSOCK, fresh R sessions on the installed package). The first seed
# is timed; if the projected total exceeds the ceiling the script stops there and writes the projection (ticket 11: a
# projection above the ceiling stops for the owner).
# Written to <out dir>: effect-runs.tsv (one row per route run), effect-pairs.tsv (per seed, m and configuration: calls
# that differ between the transforms, largest change in logFC), effect.tsv (mean, min and max over seeds per m, transform
# and configuration), effect-fixtures.tsv, effect-projection.tsv.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-177/effect"
ceiling_minutes <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 187
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
route_file <- normalizePath("docs/qa/workflow-execution/kit/kit-177/route.R")
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
started <- Sys.time()
cat("KIT-177 9.1 at", commit, "ceiling", ceiling_minutes, "minutes\n")
scales <- c(10, 40, 150, 1000, 10000)
seeds <- 20261123:20261142
configs <- c(automatic = "automatic", skip = "skip")

setupWorker <- function(route_file) {
    suppressPackageStartupMessages(loadNamespace("MultiScholaR"))
    kenv <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(route_file, kenv)
    kenv$config <- kenv$k177Config()
    assign("kenv", kenv, envir = globalenv())
    invisible(TRUE)
}

# One task: one seed, one m, one configuration, both transforms.
runTask <- function(task) {
    kenv <- get("kenv", envir = globalenv())
    sim <- kenv$k177Simulate91(task$seed, task$m)
    rows <- list()
    calls <- list()
    tables <- list()
    for (transform in c("old", "new")) {
        root <- tempfile("k177-91-")
        t0 <- Sys.time()
        manager <- NULL
        result <- tryCatch({
            manager <- kenv$k177NewManager(root, "LFQ")
            kenv$k177DesignState(manager, sim$wide, "Protein.Ids", sim$design, transform, kenv$config)
            utils::capture.output(run <- suppressMessages(suppressWarnings(kenv$k177RunPath(manager,
                kenv$k177NormParameters(ruv_mode = task$config),
                list(name = "treated.vs.control", full = "treated.vs.control=grouptreated-groupcontrol")))))
            run
        }, error = function(e) e, finally = {
            if (!is.null(manager)) try(kenv$k177Close(manager), silent = TRUE)
        })
        unlink(root, recursive = TRUE)
        seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        base <- data.frame(m = task$m, seed = task$seed, transform = transform, config = task$config,
            seconds = round(seconds, 2), stringsAsFactors = FALSE)
        if (inherits(result, "error")) {
            base$status <- paste("error", class(result)[[1L]], gsub("[\t\r\n]+", " ", substr(conditionMessage(result), 1L, 200L)))
            rows[[transform]] <- base
            next
        }
        table <- kenv$k177DaColumns(result, "Protein.Ids")
        measures <- kenv$k177Measure91(table, sim)
        calls[[transform]] <- measures$calls_ids
        tables[[transform]] <- table
        measures$calls_ids <- NULL
        base$status <- "ok"
        rows[[transform]] <- cbind(base, as.data.frame(measures), k = result$k, percentage = result$percentage,
            control_count = result$control_count, t_reported = any(!is.na(table$t)))
    }
    pair <- data.frame(m = task$m, seed = task$seed, config = task$config, calls_differ = NA_integer_,
        max_abs_logfc_change = NA_real_, stringsAsFactors = FALSE)
    if (length(tables) == 2L) {
        pair$calls_differ <- length(union(setdiff(calls$old, calls$new), setdiff(calls$new, calls$old)))
        both <- merge(tables$old[, c("id", "logFC")], tables$new[, c("id", "logFC")], by = "id", suffixes = c(".old", ".new"))
        pair$max_abs_logfc_change <- if (nrow(both)) max(abs(both$logFC.old - both$logFC.new), na.rm = TRUE) else NA_real_
    }
    list(rows = rows, pair = pair)
}

bindRows <- function(parts) {
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) return(NULL)
    columns <- unique(unlist(lapply(parts, names)))
    do.call(rbind, lapply(parts, function(part) {
        for (column in setdiff(columns, names(part))) part[[column]] <- NA
        part[, columns]
    }))
}
results <- list()
writeTables <- function() {
    runs <- bindRows(unlist(lapply(results, `[[`, "rows"), recursive = FALSE))
    pairs <- bindRows(lapply(results, `[[`, "pair"))
    if (!is.null(runs)) utils::write.table(runs, file.path(out_dir, "effect-runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(pairs)) utils::write.table(pairs, file.path(out_dir, "effect-pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(list(runs = runs, pairs = pairs))
}

cluster <- parallel::makePSOCKcluster(4L)
invisible(parallel::clusterCall(cluster, setupWorker, route_file))
tasks <- function(seed) unlist(lapply(scales, function(m) lapply(configs, function(config) list(seed = seed, m = m,
    config = config))), recursive = FALSE)

# The first seed, timed.
first_started <- Sys.time()
results <- c(results, parallel::parLapplyLB(cluster, tasks(seeds[[1L]]), runTask))
writeTables()
first_minutes <- as.numeric(difftime(Sys.time(), first_started, units = "mins"))
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
projected <- elapsed + first_minutes * (length(seeds) - 1L) + 5
projection <- data.frame(first_seed_minutes = round(first_minutes, 2), elapsed_minutes = round(elapsed, 2),
    projected_total_minutes = round(projected, 1), ceiling_minutes = ceiling_minutes, proceed = projected <= ceiling_minutes)
utils::write.table(projection, file.path(out_dir, "effect-projection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(projection, row.names = FALSE)
if (projected > ceiling_minutes) {
    cat("KIT-177 9.1: projection above the ceiling after the first seed; stopped for the owner\n")
    parallel::stopCluster(cluster)
    quit(status = 0L)
}
for (seed in seeds[-1L]) {
    results <- c(results, parallel::parLapplyLB(cluster, tasks(seed), runTask))
    writeTables()
    cat("seed", seed, "done at", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "minutes\n")
}
tables <- writeTables()

# Summary: mean, min and max over seeds per m, transform and configuration.
runs <- tables$runs
ok <- runs[runs$status == "ok", , drop = FALSE]
metrics <- intersect(c("tested", "changed_tested", "bias", "bias_up", "bias_down", "calls", "fdp", "tpr", "fdp_third1",
    "tpr_third1", "bias_third1", "fdp_third2", "tpr_third2", "bias_third2", "fdp_third3", "tpr_third3", "bias_third3", "k",
    "percentage", "control_count"), names(ok))
summary <- list()
for (m in scales) for (transform in c("old", "new")) for (config in configs) {
    part <- ok[ok$m == m & ok$transform == transform & ok$config == config, , drop = FALSE]
    for (metric in metrics) {
        x <- suppressWarnings(as.numeric(part[[metric]]))
        summary[[length(summary) + 1L]] <- data.frame(m = m, transform = transform, config = config, metric = metric,
            seeds = sum(!is.na(x)), mean = if (any(!is.na(x))) signif(mean(x, na.rm = TRUE), 6L) else NA,
            min = if (any(!is.na(x))) signif(min(x, na.rm = TRUE), 6L) else NA,
            max = if (any(!is.na(x))) signif(max(x, na.rm = TRUE), 6L) else NA, stringsAsFactors = FALSE)
    }
    errors <- sum(runs$m == m & runs$transform == transform & runs$config == config & runs$status != "ok")
    summary[[length(summary) + 1L]] <- data.frame(m = m, transform = transform, config = config, metric = "runs_failed",
        seeds = length(seeds), mean = errors, min = NA, max = NA, stringsAsFactors = FALSE)
}
pairs <- tables$pairs
for (m in scales) for (config in configs) for (metric in c("calls_differ", "max_abs_logfc_change")) {
    x <- pairs[[metric]][pairs$m == m & pairs$config == config]
    summary[[length(summary) + 1L]] <- data.frame(m = m, transform = "old_vs_new", config = config, metric = metric,
        seeds = sum(!is.na(x)), mean = if (any(!is.na(x))) signif(mean(x, na.rm = TRUE), 6L) else NA,
        min = if (any(!is.na(x))) signif(min(x, na.rm = TRUE), 6L) else NA,
        max = if (any(!is.na(x))) signif(max(x, na.rm = TRUE), 6L) else NA, stringsAsFactors = FALSE)
}
utils::write.table(do.call(rbind, summary), file.path(out_dir, "effect.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The five e2e fixtures through the same path (both configurations, both transforms): calls changed, largest change.
fixtures <- list(
    prot_lfq = list(file = "proteinGroups.txt", reader = "importMaxQuantData", type = "LFQ"),
    prot_tmt = list(file = "pd_tmt_peptides.tsv", reader = "importProteomeDiscovererTMTData", type = "TMT"),
    prot_lfq_fragpipe = list(file = "combined_protein.tsv", reader = "importFragPipeData", type = "LFQ"),
    prot_tmt_pd = list(file = "pd_tmt_peptides.tsv", reader = "importProteomeDiscovererTMTData", type = "TMT"),
    prot_lfq_maxquant = list(file = "proteinGroups.txt", reader = "importMaxQuantData", type = "LFQ"))
fixtureTask <- function(task) {
    kenv <- get("kenv", envir = globalenv())
    fixture <- task$fixture
    spec <- task$spec
    directory <- file.path(task$root, "tests", "testdata", "e2e", fixture)
    utils::capture.output(imported <- suppressMessages(get(spec$reader, envir = asNamespace("MultiScholaR"))(file.path(directory, spec$file))))
    mapping <- imported$column_mapping
    long <- as.data.frame(imported$data)
    wide <- tidyr::pivot_wider(long, id_cols = dplyr::all_of(mapping$protein_col), names_from = dplyr::all_of(mapping$run_col),
        values_from = dplyr::all_of(mapping$quantity_col))
    names(wide)[[1L]] <- "Protein.Ids"
    runs <- setdiff(names(wide), "Protein.Ids")
    values <- unlist(wide[runs])
    sample <- sub("^.*((WT|KO)_[0-9]+)$", "\\1", runs)
    design <- utils::read.delim(file.path(directory, "design.tsv"), stringsAsFactors = FALSE)
    design <- data.frame(Run = runs, group = design$group[match(sample, design$sample)],
        replicates = as.integer(sub("^.*_([0-9]+)$", "\\1", sample)), stringsAsFactors = FALSE)
    rows <- list()
    tables <- list()
    for (transform in c("old", "new")) {
        root <- tempfile("k177-fx-")
        manager <- NULL
        result <- tryCatch({
            manager <- kenv$k177NewManager(root, spec$type)
            kenv$k177DesignState(manager, wide, "Protein.Ids", design, transform, kenv$config)
            utils::capture.output(run <- suppressMessages(suppressWarnings(kenv$k177RunPath(manager,
                kenv$k177NormParameters(ruv_mode = task$config), list(name = "KO.vs.WT", full = "KO.vs.WT=groupKO-groupWT")))))
            run
        }, error = function(e) e, finally = {
            if (!is.null(manager)) try(kenv$k177Close(manager), silent = TRUE)
        })
        unlink(root, recursive = TRUE)
        if (inherits(result, "error")) {
            rows[[transform]] <- paste("error", class(result)[[1L]], gsub("[\t\r\n]+", " ", substr(conditionMessage(result), 1L, 160L)))
            next
        }
        table <- kenv$k177DaColumns(result, "Protein.Ids")
        tested <- !is.na(table$p)
        q <- rep(NA_real_, nrow(table))
        q[tested] <- stats::p.adjust(table$p[tested], "BH")
        table$called <- tested & q < 0.05
        tables[[transform]] <- table
        rows[[transform]] <- "ok"
    }
    out <- data.frame(fixture = fixture, config = task$config, proteins = nrow(wide),
        min_positive_quantity = min(values[!is.na(values) & values > 0]), zeros = sum(!is.na(values) & values == 0),
        status_old = rows$old, status_new = rows$new, tested_old = NA, tested_new = NA, calls_old = NA, calls_new = NA,
        calls_changed = NA, max_abs_logfc_change = NA, stringsAsFactors = FALSE)
    if (length(tables) == 2L) {
        out$tested_old <- sum(!is.na(tables$old$p))
        out$tested_new <- sum(!is.na(tables$new$p))
        out$calls_old <- sum(tables$old$called)
        out$calls_new <- sum(tables$new$called)
        a <- tables$old$id[tables$old$called]
        b <- tables$new$id[tables$new$called]
        out$calls_changed <- length(union(setdiff(a, b), setdiff(b, a)))
        both <- merge(tables$old[, c("id", "logFC")], tables$new[, c("id", "logFC")], by = "id", suffixes = c(".old", ".new"))
        out$max_abs_logfc_change <- if (nrow(both)) signif(max(abs(both$logFC.old - both$logFC.new), na.rm = TRUE), 6L) else NA
    }
    out
}
fixture_tasks <- unlist(lapply(names(fixtures), function(name) lapply(configs, function(config) list(fixture = name,
    spec = fixtures[[name]], config = config, root = normalizePath(".")))), recursive = FALSE)
parallel::clusterExport(cluster, "fixtureTask")
fixture_rows <- parallel::parLapplyLB(cluster, fixture_tasks, function(task) tryCatch(fixtureTask(task),
    error = function(e) data.frame(fixture = task$fixture, config = task$config, status_old = paste("error", conditionMessage(e)),
        stringsAsFactors = FALSE)))
utils::write.table(bindRows(fixture_rows), file.path(out_dir, "effect-fixtures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
parallel::stopCluster(cluster)
cat("KIT-177 9.1 done in", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "minutes\n")
