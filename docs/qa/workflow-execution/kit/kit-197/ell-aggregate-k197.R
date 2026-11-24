# KIT-197, after the implementation audit (K197-I-R1-01; the owner's ruling of 2026-11-24): the ratio rho_O of the
# added set's unchanged proteins that hold a value at or below the filter's intensity threshold in a compared group
# (ell = 1) and of the rest (ell = 0), each aggregated over the declared R7 strata with the strata's tested unchanged
# proteins as weights, against the unchanged proteins kept under both rules; pooled over the ten deciding settings
# (vd) and the same ten pairs with random missingness (mcar). POST HOC: written after the results were seen; reported,
# not deciding. No expected value is asserted here other than the self-check against committed figures.
#
# Usage, from the package root (no package loaded; base R and stats only):
#   Rscript docs/qa/workflow-execution/kit/kit-197/ell-aggregate-k197.R <kit-197 dir> <sa|sb> <out.tsv>
#
# Reads only the committed rows (<kit-197 dir>/<range>/rows/: A.tsv.gz to D.tsv.gz, s2metrics.tsv), the committed
# pooled table for the self-check (<range>/results/pooled.tsv for sb, descriptive.tsv for sa) and the committed
# scoring script score-k197.R, whose definitions (atoms, estimands, resampling) are evaluated unchanged up to its
# "main" section, so the bootstrap is the declared one: set.seed(197000001, Mersenne-Twister, Inversion, Rejection),
# one B x R index matrix per cell A to D in that order, shared by every mask and statistic; type-7 percentiles;
# delete-one-replicate jackknife. Self-check first: the pooled rho_O estimate, 2.5 % and 97.5 % limits and jackknife
# SE (vd and mcar) are re-derived and the script stops if any differs from the committed figure by more than a
# relative 1e-9. Writes only <out.tsv>.

args <- commandArgs(TRUE)
if (length(args) != 3L || !(args[[2L]] %in% c("sa", "sb")))
    stop("usage: ell-aggregate-k197.R <kit-197 dir> <sa|sb> <out.tsv>", call. = FALSE)
kit_dir <- args[[1L]]
range_name <- args[[2L]]
out_file <- args[[3L]]
SELF_CHECK_TOLERANCE <- 1e-9

# The scoring script's definitions, unchanged, in their own environment (everything before `args <- commandArgs`).
sc <- new.env(parent = globalenv())
exprs <- parse(file.path(kit_dir, "score-k197.R"), keep.source = FALSE)
main_at <- which(vapply(exprs, \(e) identical(e, quote(args <- commandArgs(TRUE))), logical(1L)))
if (length(main_at) != 1L) stop("score-k197.R: its main section was not found", call. = FALSE)
for (e in exprs[seq_len(main_at - 1L)]) eval(e, sc)

# The added set's unchanged proteins by ell, aggregated over the declared strata (sums of the strata's counts are the
# tested-null-weighted mean of the strata's rates), against O's rate.
ellStats <- function(S) {
    g <- function(name) S[, name, drop = TRUE]
    a_o <- g("p05_kept_both") / g("n0_kept_both")
    by_ell <- function(l, prefix) Reduce(`+`, lapply(sc$STRATA, \(st) g(paste0(prefix, "added_", st, "_", l))))
    n1 <- by_ell(1L, "sn0_")
    n0 <- by_ell(0L, "sn0_")
    cbind(rho_O_ell1 = (by_ell(1L, "sp05_") / n1) / a_o, rho_O_ell0 = (by_ell(0L, "sp05_") / n0) / a_o,
        share_ell1 = n1 / (n0 + n1), n0_ell1_per_r = n1 / g("one"), n0_ell0_per_r = n0 / g("one"))
}
score_stats <- sc$computeStats
sc$computeStats <- function(S) cbind(score_stats(S)[, "rho_O", drop = FALSE], ellStats(S))

rows_dir <- file.path(kit_dir, range_name, "rows")
input <- sc$readInput(rows_dir)
metrics <- input$metrics
replicates <- lapply(stats::setNames(sc$CELLS, sc$CELLS), \(cellname)
    sort(unique(metrics$r[metrics$cell == cellname])))
masks <- c(sc$VD, sc$MCAR)
atoms <- lapply(stats::setNames(sc$CELLS, sc$CELLS), \(cellname) lapply(stats::setNames(masks, masks), \(mask)
    sc$atomMatrix(input$rows[input$rows$cell == cellname & input$rows$mask == mask, , drop = FALSE],
        metrics[metrics$cell == cellname & metrics$mask == mask, , drop = FALSE], replicates[[cellname]])))
weights <- sc$bootstrapWeights(lapply(replicates, length))
poolParts <- function(masks_of) lapply(stats::setNames(names(sc$POOL), names(sc$POOL)), \(cellname)
    Reduce(`+`, lapply(masks_of(sc$POOL[[cellname]]), \(mask) atoms[[cellname]][[mask]])))
pools <- list(vd = sc$resampleStats(poolParts(identity), weights),
    mcar = sc$resampleStats(poolParts(\(m) sub("^vd", "mcar", m)), weights))
table <- rbind(sc$summariseStats(pools$vd, "pool", "ten", "vd"), sc$summariseStats(pools$mcar, "pool", "ten", "mcar"))
table$upper95 <- c(sc$quant(pools$vd$boot, 0.95), sc$quant(pools$mcar$boot, 0.95))

# Self-check against the committed pooled rho_O.
ref_file <- file.path(kit_dir, range_name, "results", if (range_name == "sb") "pooled.tsv" else "descriptive.tsv")
ref <- utils::read.delim(ref_file, stringsAsFactors = FALSE)
for (arm in c("vd", "mcar")) for (column in c("estimate", "lo95", "hi95", "jack_se")) {
    mine <- table[[column]][table$mask == arm & table$stat == "rho_O"]
    theirs <- as.numeric(ref[[column]][ref$scope == "pool" & ref$mask == arm & ref$stat == "rho_O"])
    if (length(mine) != 1L || length(theirs) != 1L || !isTRUE(abs(mine - theirs) <= SELF_CHECK_TOLERANCE * abs(theirs)))
        stop(sprintf("self-check failed: %s %s %.17g against committed %s", arm, column, mine,
            paste(theirs, collapse = ",")), call. = FALSE)
}

table$range <- range_name
table$note <- "post hoc (after the results were seen); reported, not deciding"
sc$writeTsv(table[c("range", "scope", "cell", "mask", "stat", "R", "estimate", "lo95", "hi95", "upper95", "boot_se",
    "jack_se", "boot_nonfinite", "note")], out_file)
cat(sprintf("KIT-197 ell aggregate (%s): self-check passed (pooled rho_O re-derived); written %s\n", range_name,
    out_file))
