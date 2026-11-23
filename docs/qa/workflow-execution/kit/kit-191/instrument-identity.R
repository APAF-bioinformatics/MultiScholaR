# KIT-191 U3 (ticket 10): the product's permutation test against the comparator KIT-183's simulation measured. The
# functions `fisherZ`, `pcgseT` and `comparators` are parsed from the committed docs/qa/workflow-execution/kit/kit-183/
# simulate-k183.R (their deparsed SHA-256 recorded) and run, with K$permutations = 1000, on three fixtures (n = 6, 12,
# 20 biological samples; G = 300 features; 30 sets of 5 to 40) with perm_seed = the product's seed. Expected: p-values
# identical except where a count differs on a reordering that leaves the factor unchanged (counted by the product,
# decided by rounding in the script) or where | |t*| - |t| | <= 1e-9; each such reordering listed with its set (expected:
# in the 6-sample fixture only); t within 1e-10.
#
# Usage (from the installed head's package root): Rscript instrument-identity.R <out dir>
# Writes <out dir>/functions.tsv, identity.tsv, draws.tsv, verdict.txt.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
out_dir <- arguments[[1L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(library(MultiScholaR))
ns <- asNamespace("MultiScholaR")

script <- "docs/qa/workflow-execution/kit/kit-183/simulate-k183.R"
exprs <- parse(script, keep.source = FALSE)
wanted <- c("fisherZ", "pcgseT", "comparators")
sim <- new.env(parent = globalenv())
for (e in exprs) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]]) && as.character(e[[2L]]) %in% wanted) {
        eval(e, sim)
    }
}
stopifnot(all(wanted %in% ls(sim)))
assign("K", list(permutations = 1000L, mofa_min_size = 10L), envir = sim)
functions <- data.frame(name = wanted, sha256 = vapply(wanted, function(name)
    digest::digest(paste(deparse(get(name, sim)), collapse = "\n"), algo = "sha256", serialize = FALSE), character(1)))
utils::write.table(functions, file.path(out_dir, "functions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

seed <- 20261121L
identity_rows <- list()
draw_rows <- list()
for (n in c(6L, 12L, 20L)) {
    g <- 300L
    set.seed(20261191L + n)
    f <- stats::rnorm(n)
    y <- matrix(stats::rnorm(g * n), g, n)
    y[1:20, ] <- y[1:20, ] + 0.8 * matrix(f, 20L, n, byrow = TRUE)
    ids <- sprintf("C%05d", seq_len(g))
    rownames(y) <- ids
    sizes <- sample(5:40, 30L, replace = TRUE)
    sets <- lapply(seq_len(30L), function(i) if (i <= 3L) ids[sample(1:20, min(sizes[[i]], 20L))] else sample(ids, sizes[[i]]))
    names(sets) <- sprintf("S%02d", seq_len(30L))
    # The script's comparator, as the simulation called it (weights_view NULL: MOFA2's default left out).
    script_out <- sim$comparators(list(y = y), sets, f, seed, NULL)
    # The product.
    term2gene <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE))
    product <- ns$integrationEnrichmentPermutation(list(y = y, f = f, counts = list()), term2gene,
        data.frame(term = names(sets), name = names(sets)), list(min_size = 1L, max_size = 1e6, permutations = 1000L,
            time_limit = Inf), seed)
    terms <- product$terms[match(names(sets), product$terms$ID), , drop = FALSE]
    # Each reordering's |t*| against |t|, by the script's arithmetic and the product's.
    index <- lapply(sets, function(s) match(s, ids))
    perms <- withr::with_seed(seed, vapply(seq_len(1000L), function(i) sample(n), integer(n)))
    fs <- scale(matrix(f[perms], n, 1000L))
    zp <- sim$fisherZ((t(scale(t(y))) %*% fs) / (n - 1))
    t_script <- sim$pcgseT(zp, index)
    t_obs_script <- script_out$pcgse_perm_t
    members <- lapply(index, sort)
    t_product <- ns$integrationEnrichmentSetT(zp, members)$t
    t_obs_product <- ns$integrationEnrichmentSetT(ns$integrationEnrichmentFisherZ(as.numeric(stats::cor(t(y), f))), members)$t[, 1L]
    unchanged <- vapply(seq_len(1000L), function(b) identical(f[perms[, b]], f), logical(1))
    script_ge <- abs(t_script) >= abs(t_obs_script)
    product_ge <- abs(t_product) >= abs(t_obs_product)
    product_ge[, unchanged] <- TRUE
    differ <- which(script_ge != product_ge, arr.ind = TRUE)
    if (nrow(differ)) {
        draw_rows[[length(draw_rows) + 1L]] <- data.frame(samples = n, set = names(sets)[differ[, 1L]], draw = differ[, 2L],
            unchanged = unchanged[differ[, 2L]], script_counts = script_ge[differ], product_counts = product_ge[differ],
            abs_t_star_minus_abs_t = abs(t_product[differ]) - abs(t_obs_product[differ[, 1L]]))
    }
    identity_rows[[length(identity_rows) + 1L]] <- data.frame(samples = n, set = names(sets), size = lengths(sets),
        p_script = script_out$pcgse_perm, p_product = terms$pvalue, p_identical = script_out$pcgse_perm == terms$pvalue,
        t_script = t_obs_script, t_product = t_obs_product, t_difference = abs(t_obs_script - t_obs_product),
        unchanged_draws = sum(unchanged), product_identity_draws = product$permutation$identity_draws)
}
identity <- do.call(rbind, identity_rows)
draws <- if (length(draw_rows)) do.call(rbind, draw_rows) else data.frame(samples = integer(), set = character(),
    draw = integer(), unchanged = logical(), script_counts = logical(), product_counts = logical(),
    abs_t_star_minus_abs_t = numeric())
utils::write.table(format(identity, digits = 17), file.path(out_dir, "identity.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)
utils::write.table(format(draws, digits = 17), file.path(out_dir, "draws.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
explained <- draws$unchanged | abs(draws$abs_t_star_minus_abs_t) <= 1e-9
checks <- c(
    t_within_1e10 = all(identity$t_difference <= 1e-10),
    every_difference_explained = all(explained),
    p_identical_where_no_listed_draw = all(identity$p_identical | paste(identity$samples, identity$set) %in%
        paste(draws$samples, draws$set)),
    unchanged_counted_by_product = all(identity$unchanged_draws == identity$product_identity_draws),
    unchanged_only_at_6 = all(identity$unchanged_draws[identity$samples > 6L] == 0L))
writeLines(c(sprintf("%s\t%s", names(checks), checks), sprintf("listed draws\t%d", nrow(draws)),
    sprintf("sets with p not identical\t%d", sum(!identity$p_identical))), file.path(out_dir, "verdict.txt"))
print(checks)
cat("listed draws:", nrow(draws), "; sets with p not identical:", sum(!identity$p_identical), "\n")
stopifnot(all(checks))
