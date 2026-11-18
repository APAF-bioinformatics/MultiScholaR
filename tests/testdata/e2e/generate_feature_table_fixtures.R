# E2E Fixture Generation Script: the feature-table lanes (KIT-145)
#
# DEVELOPER TOOL - not run in CI.
#
# Writes enlarged synthetic fixtures for the lanes metab_lc, metab_gc and metab_combined: 16 metabolites per assay file
# and six samples, so that every assay reaches the between-sample normalisation with more than the 10 features cyclic
# loess needs after the lanes' intensity filter (decision M17; the committed hand-made fixtures have 3 per assay, on
# which the fit is decided by rounding). All values are synthetic and come from the seeds in `METAB_FIXTURE_RECIPES`; nothing is read.
#
# Until the single re-record the files are written to `pending-kit-145/<fixture_dir>/` beside the committed lane
# directories, which stay as they are (their entry-point oracles and the manifest's sha256 are frozen on them). The
# re-record moves each pending directory over its lane directory and records `METAB_FIXTURE_RECIPES[[lane]]` as that
# lane's `fixture_recipe` in manifest.json; `generateMetabCustomFixture(lane, lane_dir)` has the signature of the
# writers in generate_fixtures.R, so it can then be listed there under the import tool `custom`.
#
# Needs withr. Run from the repository root:
#   Rscript tests/testdata/e2e/generate_feature_table_fixtures.R
#
# Idempotent: running twice produces byte-identical output files.

METAB_FIXTURE_PENDING_DIR <- "pending-kit-145"

#' The recipes, in the shape of a manifest lane's `fixture_recipe`. Each file has its own seed (the lane's seed plus
#' the file's position), a feature-name prefix and its annotation columns.
#' @noRd
METAB_FIXTURE_RECIPES <- list(
    metab_lc = list(
        generator = "tests/testdata/e2e/generate_feature_table_fixtures.R", profile = "metab_custom", seed = 20261118L,
        biological_per_group = 3L, features = 16L,
        files = list(
            lcms_pos_features.tsv = list(prefix = "LC_Pos_Feature", kind = "lc", adduct = "[M+H]+"),
            lcms_neg_features.tsv = list(prefix = "LC_Neg_Feature", kind = "lc", adduct = "[M-H]-"))),
    metab_gc = list(
        generator = "tests/testdata/e2e/generate_feature_table_fixtures.R", profile = "metab_custom", seed = 20261128L,
        biological_per_group = 3L, features = 16L,
        files = list(
            gcms_features.tsv = list(prefix = "GC_Feature", kind = "gc"))),
    metab_combined = list(
        generator = "tests/testdata/e2e/generate_feature_table_fixtures.R", profile = "metab_custom", seed = 20261138L,
        biological_per_group = 3L, features = 16L,
        files = list(
            combined_lcms_features.tsv = list(prefix = "Combined_LC_Feature", kind = "lc", platform = "LCMS_Pos"),
            combined_gcms_features.tsv = list(prefix = "Combined_GC_Feature", kind = "gc", platform = "GCMS"))))

#' One 'custom' format feature table: annotation columns, then WT_1..n and KO_1..n intensities.
#'
#' Abundances are log-normal over about two orders of magnitude (so a curve over abundance has a range to be fitted
#' on), each sample carries its own small loading effect, and three features are about twice as abundant in KO and
#' three about half (the others unchanged). No missing values, no duplicated names and no internal-standard names
#' (the lanes' `^ITSD|_d[0-9]+$` pattern finds none, as on the earlier fixtures).
#' @param spec A file entry of the recipe (`prefix`, `kind`, and `adduct` or `platform`).
#' @param features Number of features.
#' @param per_group Samples per group.
#' @param seed Seed of this file.
#' @return A data frame, one row per feature.
#' @noRd
metabCustomFeatureTable <- function(spec, features, per_group, seed) {
    withr::with_seed(seed, {
        samples <- c(paste0("WT_", seq_len(per_group)), paste0("KO_", seq_len(per_group)))
        level <- stats::runif(features, 11, 17)
        effect <- sample(rep(c(1, -1, 0), c(3L, 3L, features - 6L)))
        loading <- stats::rnorm(length(samples), 0, 0.15)
        noise <- matrix(stats::rnorm(features * length(samples), 0, 0.08), nrow = features)
        is_ko <- rep(c(0, 1), each = per_group)
        log2_values <- outer(level, rep(1, length(samples))) + outer(effect, is_ko) + outer(rep(1, features), loading) + noise
        annotation <- data.frame(Feature.Name = sprintf("%s_%03d", spec$prefix, seq_len(features)), stringsAsFactors = FALSE)
        if (identical(spec$kind, "lc")) annotation[["m.z"]] <- sprintf("%.4f", stats::runif(features, 100, 900))
        annotation[["Retention.Time"]] <- sprintf("%.2f", sort(stats::runif(features, 1, 15)))
        if (identical(spec$kind, "gc")) annotation[["Similarity"]] <- as.integer(round(stats::runif(features, 700, 950)))
        if (!is.null(spec$adduct)) annotation[["Adduct"]] <- spec$adduct
        if (!is.null(spec$platform)) annotation[["Platform"]] <- spec$platform
        values <- as.data.frame(round(2^log2_values, 1))
        names(values) <- samples
        cbind(annotation, values)
    })
}

#' Write one metabolomics 'custom' lane: its feature tables and `design.tsv`.
#' @param lane A manifest lane, with `fixture_recipe` as in `METAB_FIXTURE_RECIPES`.
#' @param lane_dir Directory to write into (created if missing).
#' @return `lane_dir`, invisibly.
#' @noRd
generateMetabCustomFixture <- function(lane, lane_dir) {
    recipe <- lane$fixture_recipe
    per_group <- as.integer(recipe$biological_per_group)
    dir.create(lane_dir, recursive = TRUE, showWarnings = FALSE)
    write <- function(table, name) {
        utils::write.table(table, file = file.path(lane_dir, name), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
    }
    for (index in seq_along(recipe$files)) {
        write(metabCustomFeatureTable(recipe$files[[index]], as.integer(recipe$features), per_group, as.integer(recipe$seed) + index),
            names(recipe$files)[[index]])
    }
    write(data.frame(sample = c(paste0("WT_", seq_len(per_group)), paste0("KO_", seq_len(per_group))),
        group = rep(c("WT", "KO"), each = per_group), stringsAsFactors = FALSE), "design.tsv")
    cat("  wrote", paste(sort(list.files(lane_dir)), collapse = ", "), "\n")
    invisible(lane_dir)
}

#' Write the pending fixtures of every recipe under `<root>/pending-kit-145/`.
#' @param root The e2e fixture root (tests/testdata/e2e).
#' @noRd
generateMetabPendingFixtures <- function(root) {
    if (!requireNamespace("withr", quietly = TRUE)) stop("Package 'withr' is required.")
    for (lane_id in names(METAB_FIXTURE_RECIPES)) {
        cat("[", lane_id, "]\n", sep = "")
        generateMetabCustomFixture(list(lane_id = lane_id, fixture_recipe = METAB_FIXTURE_RECIPES[[lane_id]]),
            file.path(root, METAB_FIXTURE_PENDING_DIR, lane_id))
    }
    invisible(NULL)
}

# Run when executed as a script (not when sourced for its functions).
local({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) && identical(basename(sub("^--file=", "", file_arg[[1L]])), "generate_feature_table_fixtures.R")) {
        generateMetabPendingFixtures(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
    }
})
