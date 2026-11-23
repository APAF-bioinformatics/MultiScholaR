# KIT-177 T1b (ticket section 5, measured; the rule fixed in the ticket): unit dependence of the later steps on plain
# log2, on the fix commit.
#
# Usage, from the package root on a run-tests runner (the fix commit installed):
#   Rscript docs/qa/workflow-execution/kit/kit-177/t1b.R [<out dir>]
#
# F177 with one planted unwanted factor and planted effects (route.R `k177F177(planted = TRUE)`), seeds 177101 to
# 177120, TMT and LFQ shapes; the design state by the package's own step (`buildProtDesignStateCheckpoint`) at scale s,
# then the label-free default path (route.R `k177RunPath`) in five configurations: `scale` (RUV skipped), legacy RUV
# `automatic`, legacy RUV `manual` (5 %, k 2), and the new designs `groups_all_features` and `first_pass_nulls` (gate on,
# k rule `num_sv_be`), each with cyclic loess except `scale`. Factors 2^-10, 2^10 (from the design step, intensity
# filter included) and 0.001, 1000 (from after the filter, on the protein set the filter kept at s = 1), each against
# s = 1 on the same protein set. Per configuration, shape, seed and factor: the largest change in logFC and in t, the
# proteins whose call changes at the route's rule (the reported fdr_qvalue < 0.05), the chosen k and control share.
# A seed is `changed` when, at any factor, a call, k or the control share differs from s = 1.
# Written: t1b-seeds.tsv (one row per configuration, shape, seed and factor), t1b.tsv (per configuration and shape: seeds
# changed of 20, seeds with a call change, largest change in t, the rule's outcome).

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-177/t1b"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
route_file <- normalizePath("docs/qa/workflow-execution/kit/kit-177/route.R")
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
started <- Sys.time()
cat("KIT-177 T1b at", commit, "\n")
seeds <- 177101:177120
configs <- list(
    scale = list(norm = "scale", ruv = "skip"),
    legacy_automatic = list(norm = "cyclicloess", ruv = "automatic"),
    legacy_manual = list(norm = "cyclicloess", ruv = "manual", percentage = 5, k = 2),
    groups_all_features = list(norm = "cyclicloess", ruv = "automatic", design = "groups_all_features"),
    first_pass_nulls = list(norm = "cyclicloess", ruv = "automatic", design = "first_pass_nulls"))

setupWorker <- function(route_file) {
    suppressPackageStartupMessages(loadNamespace("MultiScholaR"))
    options(MultiScholaR.ruv_design_panel = TRUE)
    kenv <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(route_file, kenv)
    assign("kenv", kenv, envir = globalenv())
    invisible(TRUE)
}

# One task: one seed, shape and configuration, all factors.
runTask <- function(task) {
    kenv <- get("kenv", envir = globalenv())
    ns <- asNamespace("MultiScholaR")
    sim <- kenv$k177F177(task$seed, planted = TRUE)
    config <- task$config
    norm <- kenv$k177NormParameters(norm_method = config$norm, ruv_mode = config$ruv, ruv_percentage = config$percentage %||% 5,
        ruv_k = config$k)
    if (!is.null(config$design)) norm <- c(norm, list(ruv_design = config$design, ruv_k_rule = "num_sv_be"))
    runAt <- function(s, intensity = TRUE, keep = NULL) {
        wide <- kenv$k177F177Wide(sim, s)
        if (!is.null(keep)) wide <- wide[wide$Protein.Ids %in% keep, , drop = FALSE]
        long <- tidyr::pivot_longer(wide, -"Protein.Ids", names_to = "Run", values_to = "Intensity")
        root <- tempfile("k177-t1b-")
        manager <- kenv$k177NewManager(root, task$shape)
        on.exit({ try(kenv$k177Close(manager), silent = TRUE); unlink(root, recursive = TRUE) }, add = TRUE)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        workflow$data_cln <- as.data.frame(long)
        workflow$column_mapping <- list(protein_col = "Protein.Ids", run_col = "Run", quantity_col = "Intensity")
        workflow$design_matrix <- sim$design
        workflow$config_list <- list(globalParameters = list(workflow_type = task$shape))
        get("withDIAStateScope", ns)(manager, function() get("buildProtDesignStateCheckpoint", ns)(workflow, task$shape, "Save"))
        utils::capture.output(run <- suppressMessages(suppressWarnings(kenv$k177RunPath(manager, norm,
            list(name = "treated.vs.control", full = "treated.vs.control=grouptreated-groupcontrol"), intensity = intensity))))
        kept <- if (isTRUE(intensity)) get("withDIAStateScope", ns)(manager, function() {
            object <- manager$getState("protein_intensity_filtered")
            as.character(object@protein_quant_table$Protein.Ids)
        }) else NULL
        list(da = kenv$k177DaColumns(run, "Protein.Ids"), k = run$k, percentage = run$percentage, kept = kept)
    }
    compare <- function(reference, other, s) {
        both <- merge(reference$da, other$da, by = "id", suffixes = c(".1", ".s"))
        calls1 <- both$id[!is.na(both$q.1) & both$q.1 < 0.05]
        callss <- both$id[!is.na(both$q.s) & both$q.s < 0.05]
        data.frame(config = task$name, shape = task$shape, seed = task$seed, factor = s, proteins = nrow(both),
            max_abs_logfc_change = max(abs(both$logFC.s - both$logFC.1), na.rm = TRUE),
            max_abs_t_change = max(abs(both$t.s - both$t.1), na.rm = TRUE),
            calls_changed = length(union(setdiff(calls1, callss), setdiff(callss, calls1))) +
                abs(nrow(reference$da) - nrow(other$da)),
            k_1 = reference$k %||% NA, k_s = other$k %||% NA, share_1 = reference$percentage %||% NA,
            share_s = other$percentage %||% NA, status = "ok", stringsAsFactors = FALSE)
    }
    out <- list()
    one <- tryCatch(runAt(1), error = function(e) e)
    if (inherits(one, "error")) return(data.frame(config = task$name, shape = task$shape, seed = task$seed, factor = NA,
        status = paste("error at s = 1:", gsub("[\t\r\n]+", " ", substr(conditionMessage(one), 1L, 200L))), stringsAsFactors = FALSE))
    one_after <- NULL
    for (s in c(2^-10, 2^10, 0.001, 1000)) {
        exact <- s %in% c(2^-10, 2^10)
        row <- tryCatch({
            if (!exact && is.null(one_after)) one_after <- runAt(1, intensity = FALSE, keep = one$kept)
            other <- if (exact) runAt(s) else runAt(s, intensity = FALSE, keep = one$kept)
            compare(if (exact) one else one_after, other, s)
        }, error = function(e) data.frame(config = task$name, shape = task$shape, seed = task$seed, factor = s,
            status = paste("error:", gsub("[\t\r\n]+", " ", substr(conditionMessage(e), 1L, 200L))), stringsAsFactors = FALSE))
        out[[length(out) + 1L]] <- row
    }
    bindRows(out)
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

cluster <- parallel::makePSOCKcluster(4L)
invisible(parallel::clusterCall(cluster, setupWorker, route_file))
parallel::clusterExport(cluster, "bindRows")
tasks <- list()
for (seed in seeds) for (shape in c("TMT", "LFQ")) for (name in names(configs))
    tasks[[length(tasks) + 1L]] <- list(seed = seed, shape = shape, name = name, config = configs[[name]])
rows <- parallel::parLapplyLB(cluster, tasks, runTask)
parallel::stopCluster(cluster)
seeds_table <- bindRows(rows)
utils::write.table(seeds_table, file.path(out_dir, "t1b-seeds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The rule (ticket T1b; KIT-172 r1.2 for `scale`, KIT-172 D8 for the legacy RUV rows): own ticket if a call, the chosen k
# or the control share changes in more than 1 seed of 20; otherwise a stated limitation (scale: the warning of D6) with
# the largest change in t. The new designs' rows are reported.
summary <- list()
for (name in names(configs)) for (shape in c("TMT", "LFQ")) {
    part <- seeds_table[seeds_table$config == name & seeds_table$shape == shape, , drop = FALSE]
    ok <- part[part$status == "ok", , drop = FALSE]
    changed_row <- ok$calls_changed > 0 | !mapply(identical, ok$k_1, ok$k_s) | !mapply(identical, ok$share_1, ok$share_s)
    seeds_changed <- length(unique(ok$seed[changed_row]))
    seeds_call <- length(unique(ok$seed[ok$calls_changed > 0]))
    rule <- if (name %in% c("groups_all_features", "first_pass_nulls")) "reported" else
        if (seeds_changed > 1L) "own ticket" else if (identical(name, "scale")) "warning of D6" else "stated limitation"
    summary[[length(summary) + 1L]] <- data.frame(config = name, shape = shape, seeds = length(unique(part$seed)),
        rows_failed = sum(part$status != "ok"), seeds_changed = seeds_changed, seeds_with_call_change = seeds_call,
        max_abs_t_change = if (nrow(ok)) signif(max(ok$max_abs_t_change, na.rm = TRUE), 6L) else NA,
        max_abs_logfc_change = if (nrow(ok)) signif(max(ok$max_abs_logfc_change, na.rm = TRUE), 6L) else NA,
        rule_outcome = rule, stringsAsFactors = FALSE)
}
utils::write.table(do.call(rbind, summary), file.path(out_dir, "t1b.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("KIT-177 T1b done in", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "minutes\n")
