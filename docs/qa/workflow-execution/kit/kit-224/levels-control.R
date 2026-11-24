# KIT-224 S1 (ticket section 10; acceptance criterion A3): the legacy lipidomics RUV route, base against head.
# Committed with the tests; runs in R2's scripts jobs on the fix commit, one job per part.
#
# Usage, from the package root of the head commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-224/levels-control.R <base commit> <out dir> <part>
#   <part>  `synthetic` (criteria (i) to (v), reports (a) and (b)) or `units` (criterion (vi))
#
# Mechanism (KIT-214 S1's): `git worktree add <tmp> <base>`; each side's tree is loaded once to compile, then the work
# is split into chunks, each run by a fresh Rscript process (this script, `--child`) whose working directory is that
# tree's root, with pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE); four processes at a time. Each child
# writes <out>/raw/<part>-<side>-<chunk>.rds; the comparison is made here, in base R.
#
# Part `synthetic`. Inputs: 50 seeded synthetic lipid objects r = 1..50 (seed 224000000 + r): two assays (LCMS_Pos,
# LCMS_Neg) of 300 features, KIT-172's design (groups A and B, six runs each), levels N(20, 2^2) on log2, one
# sample-wise unwanted factor uncorrelated with the groups (loadings N(0, 0.5^2)), noise SD per feature
# sqrt(4 * 0.25^2 / chi^2_4) (KIT-185's law); r even: 5 % of the cells missing, drawn among the lowest 20 % of the
# values (at most two per feature and group); r > 25: 20 % of the features carry a planted effect of +-1 log2 in group B.
# Each object is run under two variance laws, stated here before the run: `k185` (the SD above) and `level` (the SD
# above times 2^(max(0, (20 - level) / 4)): doubling every 4 log2 below level 20), from the same random draws.
# Route (the declared operations on a session owner, as the application runs them): intensity filter (percentile 1,
# proportion 0.5), QC finalise, log2 (offset 0), cyclic loess, legacy RUV Automatic (the panel's defaults,
# kit172LegacyRuv("automatic")), correlation filter (Pearson pairs within `group`, threshold 0.5), DA (~ 0 + group,
# B - A, q < 0.05) at fold-change cutoff 0 and 0.5. On the head, law `k185` also under
# MultiScholaR.ruv_feature_levels = "shifted" (side `head_shifted`). The declared path: runOmicsRuvFactors on each
# object's log2 values (law `k185`), designs `groups_all_features` and `supplied_controls` (the first 100 features of
# each assay, which carry no planted effect only when r <= 25; the design does not need them null), k = 1, rule `user`,
# and runLipidsDA with its factors.
# Criteria (compare.tsv; criterion, object, law, setting, value, pass):
#   (i)   at lipid_ruv_corrected, head - base is constant per feature on observed cells: range < 1e-9;
#   (ii)  every logFC equal to 1e-10 where the correlation filter keeps the same samples on both sides;
#   (iii) head: each feature's mean over its corrected samples equals its lipid_normalized mean there, to 1e-9;
#   (iv)  head_shifted equals base bit for bit: lipid_ruv_corrected matrices and the DA tables (identical());
#   (v)   the declared path's corrected objects and DA tables identical on base and head (sha256 of serialize()).
# Reported, not judged: changed-calls.tsv (per object, law and cutoff: calls on each side, calls changed, samples the
# correlation filter removed on each side, k and control share on each side); rates.tsv ((b): per law and side, null
# objects: share with any call at q < 0.05 and mean number of calls; effect objects: false discovery proportion and power;
# each with a Monte Carlo 95 % interval, and the paired head - base difference); priors.tsv (per object, law and assay:
# limma's prior degrees of freedom d0 on each side and the median ratio head / base of the prior variances, recomputed
# with limma::lmFit and eBayes(trend = TRUE, robust = TRUE) on the matrix the DA read, ~ 0 + group, contrast B - A).
#
# Part `units` (criterion (vi), stop rule 3): KIT-172 T1b's legacy_manual journey (kit-172/t1b-t11.R's `run`, ITSD on,
# cyclic loess, offset 0, gate on), complete data, lipid seeds 172101 to 172120 and metab seed 172101, units 2^-10, 2^10
# (with the intensity filter) and 0.001, 1000 (after it, on the features kept at 1), each against 1: calls (sign of logFC
# where q < 0.05), max |delta t|, max |delta logFC| per row and side (units.tsv). Head: every call identical and
# max |delta t| < 1e-6 in all 84 rows (compare.tsv); base recorded.
#
# Exit status non-zero when a child fails or any criterion row does not pass.

args <- commandArgs(trailingOnly = TRUE)

S1_RUNS <- c(sprintf("A_%d", 1:6), sprintf("B_%d", 1:6))
S1_CONTRASTS <- data.frame(contrasts = "groupB-groupA", friendly_names = "B_vs_A", full_format = "B_vs_A=groupB-groupA",
    stringsAsFactors = FALSE)

# ---- Child: one side, one chunk ----------------------------------------------------------------------------------------
childMain <- function(tree, side, part, chunk, chunks, out) {
    setwd(tree)
    suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
    options(MultiScholaR.ruv_design_panel = TRUE)
    quiet <- function(expr) {
        value <- NULL
        utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
        value
    }
    binding <- lipidomicsDeclarationBinding
    stateMats <- function(owner, state) {
        object <- owner$getState(state)
        id <- methods::slot(object, binding$id_slot)
        lapply(methods::slot(object, binding$data_slot), function(table) {
            table <- as.data.frame(table)
            runs <- intersect(S1_RUNS, names(table))
            m <- as.matrix(table[, runs, drop = FALSE])
            rownames(m) <- as.character(table[[id]])
            m
        })
    }
    daTable <- function(da) {
        long <- as.data.frame(da$results[[binding$long_table]])
        long <- long[intersect(c("lipid_id", "assay", "comparison", "logFC", "t", "P.Value", "fdr_qvalue", "AveExpr",
            "significant"), names(long))]
        long[order(long$assay, long$comparison, long$lipid_id), , drop = FALSE]
    }
    results <- list()

    if (identical(part, "synthetic")) {
        values <- function(r, law) {
            set.seed(224000000L + r)
            groups <- rep(c("A", "B"), each = 6L)
            f <- stats::residuals(stats::lm(stats::rnorm(12L) ~ factor(groups)))
            f <- f / sqrt(mean(f^2))
            mats <- list()
            truth <- list()
            for (assay in c("LCMS_Pos", "LCMS_Neg")) {
                p <- 300L
                level <- stats::rnorm(p, 20, 2)
                sd <- sqrt(4 * 0.25^2 / stats::rchisq(p, 4))
                if (identical(law, "level")) sd <- sd * 2^(pmax(0, (20 - level) / 4))
                effect <- rep(0, p)
                if (r > 25L) {
                    changed <- sample.int(p, round(0.2 * p))
                    effect[changed] <- sample(c(-1, 1), length(changed), replace = TRUE)
                }
                loadings <- stats::rnorm(p, 0, 0.5)
                y <- level + outer(effect, as.numeric(groups == "B")) + outer(loadings, f) + matrix(stats::rnorm(p * 12L), p) * sd
                ids <- sprintf("%s_F%03d", assay, seq_len(p))
                dimnames(y) <- list(ids, S1_RUNS)
                if (r %% 2L == 0L) {
                    low <- which(y <= stats::quantile(y, 0.2), arr.ind = TRUE)
                    low <- low[sample.int(nrow(low)), , drop = FALSE]
                    key <- paste(low[, 1L], groups[low[, 2L]])
                    low <- low[stats::ave(seq_along(key), key, FUN = seq_along) <= 2L, , drop = FALSE]
                    y[low[seq_len(min(nrow(low), round(0.05 * length(y)))), , drop = FALSE]] <- NA_real_
                }
                mats[[assay]] <- y
                truth[[assay]] <- stats::setNames(effect, ids)
            }
            list(log2 = mats, truth = truth)
        }
        priors <- function(mats) {
            lapply(mats, function(m) {
                groups <- factor(sub("_.*$", "", colnames(m)), levels = c("A", "B"))
                design <- stats::model.matrix(~ 0 + groups)
                colnames(design) <- levels(groups)
                fit <- limma::lmFit(m, design)
                fit <- limma::contrasts.fit(fit, limma::makeContrasts(B - A, levels = design))
                fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
                list(s2_prior = stats::setNames(rep_len(fit$s2.prior, nrow(m)), rownames(m)),
                    df_prior = stats::median(fit$df.prior))
            })
        }
        journey <- function(input) {
            root <- tempfile("kit224-s1-")
            fixture <- kit172Owner(root, kit172Object(input, "lipid"), "lipid")
            owner <- fixture$owner
            on.exit({
                closeDIAStateAdapter(owner)
                unlink(root, recursive = TRUE)
            }, add = TRUE)
            config <- fixture$workflow$config_list
            run <- quiet({
                runWorkflowOmicsReferenceIntensityFilter(binding, owner, 1, 0.5, config)
                commitWorkflowOmicsStoreStep(owner, "lipid_qc_finalize", "lipid_qc_complete", config, "QC complete")
                runWorkflowOmicsReferenceLog2(binding, owner, 0, config)
                runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, owner, "cyclicloess", config)
                ruv <- runWorkflowOmicsReferenceRuv(binding, owner, kit172LegacyRuv("automatic"), config)
                normalized <- stateMats(owner, "lipid_normalized")
                corrected <- stateMats(owner, "lipid_ruv_corrected")
                runWorkflowOmicsReferenceCorrelation(binding, owner, "group", 0.5, config)
                filtered <- stateMats(owner, owner$current_state)
                da <- lapply(c(cutoff_0 = 0, cutoff_0.5 = 0.5), function(cutoff) {
                    daTable(runWorkflowOmicsReferenceDaAnalysis(binding, owner, getWorkflowStateReference(owner),
                        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = cutoff), S1_CONTRASTS))
                })
                one <- function(x) if (is.null(x)) NA_character_ else paste(names(x), unlist(x), sep = "=", collapse = ";")
                list(normalized = normalized, corrected = corrected, kept = lapply(filtered, colnames), da = da,
                    priors = priors(filtered), k = one(ruv$value$best_k), share = one(ruv$value$percentage))
            })
            run
        }
        declared <- function(input) {
            object <- kit172Object(input, "lipid")
            supplied <- lapply(input, function(m) rownames(m)[1:100])
            one <- function(design) {
                run <- if (identical(design, "supplied_controls")) {
                    quiet(runOmicsRuvFactors(object, design, 1L, "user", ruv_control_ids = supplied))
                } else quiet(runOmicsRuvFactors(object, design, 1L, "user"))
                da <- quiet(runLipidsDA(object, S1_CONTRASTS[c("contrasts", "friendly_names")], "~ 0 + group", 0.05, 0, TRUE, TRUE,
                    ruv_factors = run$factors))$da_lipids_long
                sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)
                c(corrected = sha(run$corrected), da = sha(as.data.frame(da)))
            }
            list(groups_all_features = one("groups_all_features"), supplied_controls = one("supplied_controls"))
        }
        tasks <- seq_len(50L)
        tasks <- tasks[(tasks - 1L) %% chunks + 1L == chunk]
        for (r in tasks) {
            for (law in c("k185", "level")) {
                v <- values(r, law)
                variants <- if (identical(side, "head") && identical(law, "k185")) c("head", "head_shifted") else side
                for (variant in variants) {
                    key <- paste(r, law, variant, sep = "|")
                    results[[key]] <- tryCatch({
                        options <- if (identical(variant, "head_shifted")) list(MultiScholaR.ruv_feature_levels = "shifted") else list()
                        x <- withr::with_options(options, journey(lapply(v$log2, function(m) 2^m)))
                        x$truth <- v$truth
                        x$r <- r
                        x$law <- law
                        x$side <- variant
                        if (identical(law, "k185") && !identical(variant, "head_shifted")) x$declared <- declared(v$log2)
                        x
                    }, error = function(e) list(r = r, law = law, side = variant, error = conditionMessage(e)))
                    cat(sprintf("KIT-224 S1 child %s: object %d law %s %s\n", side, r, law,
                        if (is.null(results[[key]]$error)) "done" else paste("error:", results[[key]]$error)))
                }
            }
        }
    } else {
        factors <- c(`2^-10` = 2^-10, `2^10` = 2^10, `0.001` = 0.001, `1000` = 1000)
        calls <- function(da) ifelse(!is.na(da$fdr_qvalue) & da$fdr_qvalue < 0.05, sign(da$logFC), 0)
        run <- function(values, class, filter = TRUE) {
            root <- tempfile("kit224-s1u-")
            journey <- kit172Journey(root, kit172Object(values, class), class, itsd = TRUE, method = "cyclicloess",
                filter = filter, offset = 0, ruv = kit172LegacyRuv("manual"))
            on.exit({
                closeDIAStateAdapter(journey$owner)
                unlink(root, recursive = TRUE)
            }, add = TRUE)
            list(kept = rownames(kit172StateMatrix(journey, paste0(journey$prefix, if (filter) "_intensity_filtered" else "_qc_complete"))),
                da = kit172DaTable(journey))
        }
        tasks <- c(lapply(172101:172120, function(seed) list(seed = seed, class = "lipid")), list(list(seed = 172101L, class = "metab")))
        tasks <- tasks[(seq_along(tasks) - 1L) %% chunks + 1L == chunk]
        for (task in tasks) {
            key <- paste(task$seed, task$class, side, sep = "|")
            results[[key]] <- tryCatch({
                values <- kit172Values(task$seed, planted = TRUE, missing = FALSE)$values
                reference <- kit172Quiet(run(values, task$class))
                rows <- lapply(names(factors), function(name) {
                    s <- factors[[name]]
                    power <- name %in% c("2^-10", "2^10")
                    input <- if (power) values * s else values[reference$kept, , drop = FALSE] * s
                    got <- kit172Quiet(run(input, task$class, filter = power))
                    k <- function(d) paste(d$assay, d$comparison, d$feature)
                    at <- match(k(reference$da), k(got$da))
                    data.frame(seed = task$seed, class = task$class, factor = name, side = side,
                        rows = nrow(reference$da), matched = sum(!is.na(at)),
                        call_changes = sum(calls(got$da)[at] != calls(reference$da), na.rm = TRUE),
                        max_dt = max(abs(got$da$t[at] - reference$da$t), na.rm = TRUE),
                        max_dlogfc = max(abs(got$da$logFC[at] - reference$da$logFC), na.rm = TRUE),
                        calls_1 = sum(calls(reference$da) != 0), stringsAsFactors = FALSE)
                })
                do.call(rbind, rows)
            }, error = function(e) data.frame(seed = task$seed, class = task$class, factor = NA_character_, side = side,
                error = conditionMessage(e), stringsAsFactors = FALSE))
            cat(sprintf("KIT-224 S1 child %s: seed %d %s done\n", side, task$seed, task$class))
        }
    }
    saveRDS(results, file.path(out, "raw", sprintf("%s-%s-%d.rds", part, side, chunk)))
    quit(status = 0L)
}

if (length(args) >= 1L && identical(args[[1L]], "--child")) {
    childMain(args[[2L]], args[[3L]], args[[4L]], as.integer(args[[5L]]), as.integer(args[[6L]]), args[[7L]])
}

# ---- Parent: both sides, then the comparison ----------------------------------------------------------------------------
stopifnot(length(args) == 3L, grepl("^[0-9a-f]{40}$", args[[1L]]), args[[3L]] %in% c("synthetic", "units"))
base_commit <- args[[1L]]
part <- args[[3L]]
dir.create(file.path(args[[2L]], "raw"), recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[2L]])
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-224 S1", part, ": head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")
base_root <- file.path(tempdir(), "kit224-s1-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_root), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_root, "R")))
trees <- c(base = base_root, head = head_root)
started <- Sys.time()
# Compile each tree once before the parallel children (they would race on the same object files).
for (side in names(trees)) {
    code <- system2(rscript, c("-e", shQuote(sprintf("suppressPackageStartupMessages(pkgload::load_all('%s', quiet = TRUE))",
        trees[[side]]))), stdout = file.path(out, paste0(part, "-compile-", side, ".txt")),
        stderr = file.path(out, paste0(part, "-compile-", side, ".txt")))
    cat("KIT-224 S1 compile", side, "exit", code, "\n")
}
chunks <- 4L
jobs <- expand.grid(chunk = seq_len(chunks), side = names(trees), stringsAsFactors = FALSE)
status <- unlist(parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    side <- jobs$side[[i]]
    chunk <- jobs$chunk[[i]]
    log <- file.path(out, sprintf("%s-%s-%d-console.txt", part, side, chunk))
    system2(rscript, c(shQuote(script), "--child", shQuote(trees[[side]]), side, part, chunk, chunks, shQuote(out)),
        stdout = log, stderr = log)
}, mc.cores = 4L, mc.preschedule = FALSE))
sides <- data.frame(part = part, side = jobs$side, chunk = jobs$chunk, exit = status,
    commit = ifelse(jobs$side == "base", base_commit, head_commit), stringsAsFactors = FALSE)
utils::write.table(sides, file.path(out, paste0(part, "-sides.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("KIT-224 S1 %s: children done in %.1f minutes, exits %s\n", part,
    as.numeric(difftime(Sys.time(), started, units = "mins")), paste(status, collapse = ",")))

raw <- list()
for (side in names(trees)) for (chunk in seq_len(chunks)) {
    path <- file.path(out, "raw", sprintf("%s-%s-%d.rds", part, side, chunk))
    if (file.exists(path)) raw <- c(raw, readRDS(path))
}
write_tsv <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
row <- function(criterion, object, law, setting, value, pass) {
    data.frame(criterion = criterion, object = as.character(object), law = law, setting = setting,
        value = as.character(value), pass = isTRUE(pass), stringsAsFactors = FALSE)
}
compare <- list()

if (identical(part, "units")) {
    units <- do.call(rbind, lapply(raw, function(x) {
        for (c in setdiff(c("seed", "class", "factor", "side", "rows", "matched", "call_changes", "max_dt", "max_dlogfc",
            "calls_1", "error"), names(x))) x[[c]] <- NA
        x[c("seed", "class", "factor", "side", "rows", "matched", "call_changes", "max_dt", "max_dlogfc", "calls_1", "error")]
    }))
    units <- units[order(units$side, units$class, units$seed, units$factor), , drop = FALSE]
    write_tsv(units, "units.tsv")
    head_rows <- units[units$side == "head", , drop = FALSE]
    for (i in seq_len(nrow(head_rows))) {
        u <- head_rows[i, ]
        ok <- is.na(u$error) && isTRUE(u$call_changes == 0L) && isTRUE(u$max_dt < 1e-6) && isTRUE(u$matched == u$rows)
        compare[[length(compare) + 1L]] <- row("vi_head_unit_invariant", paste(u$seed, u$class), "", u$factor,
            sprintf("calls changed %s; max |dt| %s; max |dlogFC| %s%s", u$call_changes, signif(u$max_dt, 4), signif(u$max_dlogfc, 4),
                if (is.na(u$error)) "" else paste("; error:", u$error)), ok)
    }
    compare[[length(compare) + 1L]] <- row("vi_rows", "all", "", "head rows", nrow(head_rows), nrow(head_rows) == 84L)
    base_rows <- units[units$side == "base", , drop = FALSE]
    cat(sprintf("KIT-224 S1 units: base rows with call changes %d of %d (largest %s); head rows with call changes %d of %d\n",
        sum(base_rows$call_changes > 0, na.rm = TRUE), nrow(base_rows), max(base_rows$call_changes, na.rm = TRUE),
        sum(head_rows$call_changes > 0, na.rm = TRUE), nrow(head_rows)))
} else {
    getResult <- function(r, law, side) raw[[paste(r, law, side, sep = "|")]]
    errors <- Filter(function(x) !is.null(x$error), raw)
    for (x in errors) compare[[length(compare) + 1L]] <- row("child_error", x$r, x$law, x$side, x$error, FALSE)
    calls <- function(da) ifelse(!is.na(da$fdr_qvalue) & da$fdr_qvalue < 0.05, sign(da$logFC), 0)
    key <- function(d) paste(d$assay, d$comparison, d$lipid_id)
    changed <- list()
    per_object <- list()
    prior_rows <- list()
    for (r in seq_len(50L)) for (law in c("k185", "level")) {
        b <- getResult(r, law, "base")
        h <- getResult(r, law, "head")
        if (is.null(b) || is.null(h) || !is.null(b$error) || !is.null(h$error)) {
            compare[[length(compare) + 1L]] <- row("missing_result", r, law, "", "a side has no result", FALSE)
            next
        }
        for (assay in names(h$corrected)) {
            hc <- h$corrected[[assay]]
            bc <- b$corrected[[assay]]
            same_shape <- identical(dimnames(hc), dimnames(bc))
            d <- if (same_shape) hc - bc else NULL
            spread <- if (same_shape) max(apply(d, 1L, function(x) { x <- x[is.finite(x)]; if (length(x)) diff(range(x)) else 0 })) else NA
            compare[[length(compare) + 1L]] <- row("i_constant_per_feature", r, law, assay,
                if (same_shape) sprintf("%.3g", spread) else "dimnames differ", same_shape && spread < 1e-9)
            before <- h$normalized[[assay]][rownames(hc), colnames(hc), drop = FALSE]
            both <- is.finite(hc) & is.finite(before)
            mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
            level <- max(abs(mean_in(hc) - mean_in(before)), na.rm = TRUE)
            compare[[length(compare) + 1L]] <- row("iii_head_levels_kept", r, law, assay, sprintf("%.3g", level), level < 1e-9)
            bp <- b$priors[[assay]]
            hp <- h$priors[[assay]]
            prior_rows[[length(prior_rows) + 1L]] <- data.frame(object = r, law = law, assay = assay, d0_base = bp$df_prior,
                d0_head = hp$df_prior, median_s2prior_ratio = if (identical(names(hp$s2_prior), names(bp$s2_prior)))
                    stats::median(hp$s2_prior / bp$s2_prior) else NA_real_, stringsAsFactors = FALSE)
        }
        same_samples <- identical(b$kept, h$kept)
        removed <- function(x) paste(setdiff(S1_RUNS, unique(unlist(x$kept))), collapse = ",")
        for (cutoff in names(h$da)) {
            hd <- h$da[[cutoff]]
            bd <- b$da[[cutoff]]
            at <- match(key(hd), key(bd))
            if (same_samples) {
                dl <- max(abs(hd$logFC - bd$logFC[at]), na.rm = TRUE)
                compare[[length(compare) + 1L]] <- row("ii_logfc_equal", r, law, cutoff, sprintf("%.3g", dl),
                    !anyNA(at) && nrow(hd) == nrow(bd) && dl < 1e-10)
            }
            changed[[length(changed) + 1L]] <- data.frame(object = r, law = law, cutoff = cutoff, missing = r %% 2L == 0L,
                effect = r > 25L, calls_base = sum(calls(bd) != 0), calls_head = sum(calls(hd) != 0),
                calls_changed = sum(calls(hd) != calls(bd)[at], na.rm = TRUE), same_samples = same_samples,
                removed_base = removed(b), removed_head = removed(h), k_base = b$k, k_head = h$k, share_base = b$share,
                share_head = h$share, stringsAsFactors = FALSE)
            if (identical(cutoff, "cutoff_0")) {
                for (side_name in c("base", "head")) {
                    x <- if (identical(side_name, "base")) b else h
                    da <- x$da[[cutoff]]
                    truth <- unlist(unname(x$truth))
                    truth <- truth[da$lipid_id]
                    called <- calls(da) != 0
                    true_called <- called & !is.na(truth) & truth != 0 & sign(da$logFC) == sign(truth)
                    per_object[[length(per_object) + 1L]] <- data.frame(object = r, law = law, side = side_name,
                        effect = r > 25L, any_call = any(called), calls = sum(called),
                        fdp = if (sum(called)) sum(called & (is.na(truth) | truth == 0)) / sum(called) else 0,
                        power = if (sum(truth != 0, na.rm = TRUE)) sum(true_called) / sum(truth != 0, na.rm = TRUE) else NA_real_,
                        stringsAsFactors = FALSE)
                }
            }
        }
        if (identical(law, "k185")) {
            hs <- getResult(r, law, "head_shifted")
            ok <- !is.null(hs) && is.null(hs$error) && identical(hs$corrected, b$corrected) && identical(hs$da, b$da)
            compare[[length(compare) + 1L]] <- row("iv_earlier_rule_equals_base", r, law, "corrected and DA",
                if (is.null(hs) || !is.null(hs$error)) "no result" else ok, ok)
            for (design in names(h$declared)) {
                compare[[length(compare) + 1L]] <- row("v_declared_identical", r, law, design,
                    paste(h$declared[[design]], collapse = ","), identical(h$declared[[design]], b$declared[[design]]))
            }
        }
    }
    write_tsv(do.call(rbind, changed), "changed-calls.tsv")
    write_tsv(do.call(rbind, prior_rows), "priors.tsv")
    objects <- do.call(rbind, per_object)
    write_tsv(objects, "per-object-rates.tsv")
    interval <- function(x) {
        x <- x[!is.na(x)]
        m <- mean(x)
        half <- if (length(x) > 1L) 1.96 * stats::sd(x) / sqrt(length(x)) else NA_real_
        c(n = length(x), mean = m, lower = m - half, upper = m + half)
    }
    rates <- list()
    for (law in c("k185", "level")) {
        for (effect in c(FALSE, TRUE)) {
            sel <- objects[objects$law == law & objects$effect == effect, , drop = FALSE]
            measures <- if (effect) c("fdp", "power", "calls") else c("any_call", "calls")
            for (measure in measures) {
                bv <- as.numeric(sel[sel$side == "base", measure][order(sel$object[sel$side == "base"])])
                hv <- as.numeric(sel[sel$side == "head", measure][order(sel$object[sel$side == "head"])])
                for (what in c("base", "head", "head_minus_base")) {
                    x <- switch(what, base = bv, head = hv, head_minus_base = if (length(hv) == length(bv)) hv - bv else NA_real_)
                    s <- interval(x)
                    rates[[length(rates) + 1L]] <- data.frame(law = law, objects = if (effect) "effect" else "null", measure = measure,
                        side = what, n = s[["n"]], mean = s[["mean"]], lower_95 = s[["lower"]], upper_95 = s[["upper"]],
                        stringsAsFactors = FALSE)
                }
            }
        }
    }
    write_tsv(do.call(rbind, rates), "rates.tsv")
    cat(sprintf("KIT-224 S1 synthetic: %d object-law-cutoff rows; %d with a call changed\n", length(changed),
        sum(vapply(changed, function(x) x$calls_changed > 0, logical(1)))))
}

result <- do.call(rbind, compare)
write_tsv(result, if (identical(part, "units")) "compare-units.tsv" else "compare.tsv")
cat(sprintf("KIT-224 S1 %s: %d criterion rows, %d not passed\n", part, nrow(result), sum(!result$pass)))
if (any(!result$pass)) print(utils::head(result[!result$pass, ], 40L), right = FALSE)
quit(status = if (any(status != 0L) || any(!result$pass)) 1L else 0L)
