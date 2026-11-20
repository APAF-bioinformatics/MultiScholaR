# KIT-137: tier 1 of the synthetic multi-omics corpus for integration (the specification is
# docs/qa/workflow-execution/kit/kit-060/AUDIT-integration-touchpoints.md, section 6).
#
# Four omics measured on the same twenty biological samples (groups WT and KO, ten each):
# proteomics (one view, one run per biological sample), transcriptomics (one view), lipidomics
# (views LCMS_Pos and LCMS_Neg) and metabolomics (views LCMS_Pos and GCMS); the last three carry
# two technical replicates per biological sample, as their generators lay them out.
#
# The generators are called, never edited (tools/profiling/ is harness-locked):
#   proteomics       kitProteinRecipe("lfq", missingness = "mnar"), kitProteinMatrix
#   lipidomics       kitOmicsRecipe("lipid_screen"), kitOmicsDesign, kitOmicsAssayMatrix
#   metabolomics     kitOmicsRecipe("metab_screen"), kitOmicsDesign, kitOmicsAssayMatrix
#   transcriptomics  kitOmicsRecipe("trx_screen"), kitOmicsDesign, kitTrxCounts
# Their own group, differential and batch effects are set to zero through recipe fields (the
# lipid, metabolite and transcript batches would otherwise share one assignment across three
# omics). They supply base abundances, biological and technical noise, sample scale or library
# size, and the missing values: intensity-dependent for proteomics, the lowest-abundance 20% of
# features missing in 15% of their runs for lipidomics and metabolomics, none for counts.
# Everything known is added by this composer on the log2 scale, from one seed:
#   Z1     the group (KO = 1); every view; 15% of features
#   Z2     continuous; every view; 10% of features
#   Z3     continuous; proteomics and transcriptomics only; 10% of features
#   batch  each omic's own assignment of the biological samples, balanced over the groups and
#          different per omic; 30% of features, a N(0, 0.5) shift per feature and batch
# Z2 and Z3 are drawn standard normal per biological sample, orthogonalised in the sample against
# the intercept, Z1 and each other, and scaled to unit variance. Loadings are a random sign times
# U(0.8, 1.6) for Z1 and U(0.5, 1.0) for Z2 and Z3. Technical replicates share their biological
# sample's scores. Mass-spectrometry views: log2 intensity plus the effects, then run medians
# equalised. Transcriptomics: counts scaled by 2^effect and rounded, then log2 CPM
# (log2((count + 0.5) / (library + 1) * 1e6)). The biological sample of recipe$dropped is removed
# from that omic's files after generation (decision D4: it stays missing, not dropped elsewhere).
#
# Layout: each <omic>/ folder is one omic's integration input in the form of the contract
# (R/utils_workflow_integration_input.R; KIT-147: contract v2, Parquet; readIntegrationInput()
# accepts it), the tables under the contract's fixed names, written by the package's own Parquet
# writer (writeIntegrationInputParquet: format 2.6, snappy, row groups of 10,000 rows), so the
# package must be installed or loaded:
#   <omic>/matrix_<view>.parquet    features in rows; feature_id, then one double column per run
#                               (the omic's own run names, in a seeded order); log2 values, 4
#                               decimals; NA kept
#   <omic>/features_<view>.parquet  feature_id, display_name (text)
#   <omic>/samples.parquet      Run, group, replicates, tech_rep_group, batch (text); one row per
#                               matrix column, in column order. tech_rep_group is the cross-omic key
#                               in the omics' own format, <group>_<replicate> unpadded (KO_5), and the
#                               manifest names it as the biological-sample column, as the omics'
#                               writers do (KIT-062)
#   <omic>/integration_input.json  the contract's manifest: matrix_role "normalised", no
#                               correction, the Parquet layout, every table's sha256. Written here,
#                               not through writeIntegrationInput(), so that the file holds nothing
#                               that varies by machine, package version or time: no version, and
#                               written_at fixed to the corpus's date
#   truth_scores.tsv            biological_sample (the tech_rep_group values), group, Z1, Z2, Z3
#                               (6 decimals)
#   truth_loadings.tsv          omic, view, feature_id, Z1, Z2, Z3 loadings and batch_B1..B3
#                               effects (NA where the omic has no such batch)
#   truth.json                  views, active omics per factor, batch per omic and sample, the
#                               dropped sample, run, feature and missing-cell counts, every recipe
#   manifest.json               bytes and sha256 of every other file
# Views per omic, named as the omics' writers name them (INT-4: <omic>_<assay> for an omic with
# assays, the omic's name for a single-view omic): proteomics/proteomics,
# transcriptomics/transcriptomics, lipidomics/lipidomics_LCMS_Pos, lipidomics/lipidomics_LCMS_Neg,
# metabolomics/metabolomics_LCMS_Pos, metabolomics/metabolomics_GCMS.
#
# Tier 2 (the same design run through the omics journeys from vendor files) is not built here.
#
# Usage, from the repository root, with MultiScholaR installed (its contract writes the Parquet
# files; the committed copy is written on a runner, never on the laptop):
#   Rscript tests/testdata/integration/compose_integration_corpus.R [directory]
# (default directory tests/testdata/integration/corpus, the committed copy).

INTEGRATION_CORPUS_SEED <- 20261117L

# The route each omic's integration input names (KIT-060 audit 4.2) and its values' scale.
INTEGRATION_CORPUS_WORKFLOWS <- c(proteomics = "LFQ", transcriptomics = "TRANSCRIPTOMICS",
    lipidomics = "LIPIDOMICS", metabolomics = "METABOLOMICS")
INTEGRATION_CORPUS_SCALES <- c(proteomics = "log2", transcriptomics = "log2_cpm",
    lipidomics = "log2", metabolomics = "log2")

integrationCorpusRecipe <- function(seed = INTEGRATION_CORPUS_SEED) {
    all_omics <- c("proteomics", "transcriptomics", "lipidomics", "metabolomics")
    list(
        corpus = "KIT-137 tier 1",
        synthetic = TRUE,
        seed = as.integer(seed),
        biological_per_group = 10L,
        groups = c("WT", "KO"),
        effect_group = "KO",
        omics = list(
            proteomics = list(generator = "lfq", features = 400L, seed_offset = 10000L,
                batches = c("B1", "B2")),
            transcriptomics = list(generator = "trx_screen", features = 400L, seed_offset = 40000L,
                run_prefix = "trx_", batches = c("B1", "B2", "B3")),
            lipidomics = list(generator = "lipid_screen", features = 150L, seed_offset = 20000L,
                run_prefix = "lip_", feature_prefix = "lipid", base_mean = 16,
                missing_features = 0.2, batches = c("B1", "B2", "B3")),
            metabolomics = list(generator = "metab_screen", features = 150L, seed_offset = 30000L,
                run_prefix = "met_", feature_prefix = "metab", base_mean = 15,
                missing_features = 0.2, batches = c("B1", "B2"))
        ),
        factors = list(
            Z1 = list(active = all_omics, share = 0.15, abs_loading = c(0.8, 1.6)),
            Z2 = list(active = all_omics, share = 0.10, abs_loading = c(0.5, 1.0)),
            Z3 = list(active = c("proteomics", "transcriptomics"), share = 0.10,
                abs_loading = c(0.5, 1.0))
        ),
        batch = list(share = 0.3, sd = 0.5),
        dropped = list(omic = "metabolomics", biological_sample = "KO_7"),
        digits = 4L
    )
}

# --- Generators ---

integrationGenerators <- function(repo) {
    gen <- new.env(parent = globalenv())
    for (file in c("kit_omics_corpus.R", "kit_protein_corpus.R")) {
        sys.source(file.path(repo, "tools", "profiling", file), envir = gen)
    }
    gen
}

# The biological-sample key in the omics' format (`paste(group, replicates, sep = "_")`,
# R/func_omics_design_shared.R); kept internally as `biological_sample`, written as tech_rep_group.
integrationBiologicalKey <- function(group, replicate) {
    paste(group, as.integer(replicate), sep = "_")
}

integrationSamples <- function(runs, group, replicate) {
    key <- integrationBiologicalKey(group, replicate)
    data.frame(Run = runs, group = group, replicates = as.integer(replicate),
        tech_rep_group = key, biological_sample = key, stringsAsFactors = FALSE)
}

integrationBaseProteomics <- function(gen, recipe, spec) {
    n <- length(recipe$groups) * recipe$biological_per_group
    r <- gen$kitProteinRecipe("lfq", samples = n, proteins = spec$features, missingness = "mnar")
    r$seed <- recipe$seed + spec$seed_offset
    r$differential_fraction <- 0
    m <- gen$kitProteinMatrix(r)
    replicate <- as.integer(stats::ave(seq_len(n), m$group, FUN = seq_along))
    values <- log2(m$values)
    values[m$missing] <- NA_real_
    # Run names as the LFQ import cleans MaxQuant's sample names (wt_001).
    dimnames(values) <- list(sprintf("P%05d", seq_len(spec$features)),
        tolower(sprintf("%s_%03d", m$group, replicate)))
    list(generator = r, scale = "log2", views = list(proteomics = values),
        display = list(proteomics = sprintf("GENE%d", seq_len(spec$features))),
        samples = integrationSamples(colnames(values), m$group, replicate))
}

integrationBaseMs <- function(gen, recipe, spec, omic) {
    r <- gen$kitOmicsRecipe(spec$generator, recipe$biological_per_group, spec$features,
        standards = 0L, missing_features = spec$missing_features)
    r$seed <- recipe$seed + spec$seed_offset
    r$differential_fraction <- 0
    r$batch_effect_sd <- 0
    design <- gen$kitOmicsDesign(r)
    runs <- paste0(spec$run_prefix, design$sample)
    views <- lapply(seq_along(r$assays), \(i) {
        m <- gen$kitOmicsAssayMatrix(r, design, r$seed + i, base_mean = spec$base_mean)
        values <- log2(m$values)
        dimnames(values) <- list(
            sprintf("%s_%s_%03d", spec$feature_prefix, tolower(r$assays[[i]]), seq_len(nrow(values))),
            runs)
        values
    })
    # INT-4: the writers' view name, <omic>_<assay> (integrationInputViewName in the package).
    names(views) <- paste(omic, gsub("[^A-Za-z0-9_-]", "_", r$assays), sep = "_")
    list(generator = r, scale = "log2", views = views, display = lapply(views, rownames),
        samples = integrationSamples(runs, design$group, design$replicate))
}

integrationBaseTrx <- function(gen, recipe, spec) {
    r <- gen$kitOmicsRecipe("trx_screen", recipe$biological_per_group, spec$features,
        low_count_fraction = 0)
    r$seed <- recipe$seed + spec$seed_offset
    r$differential_fraction <- 0
    r$batch_fraction <- 0
    r$zero_genes <- 0L
    r$constant_genes <- 0L
    design <- gen$kitOmicsDesign(r)
    # The generator's groups (RPMI, Sera) are this corpus's (WT, KO), in that order.
    group <- unname(stats::setNames(recipe$groups, r$groups)[design$group])
    runs <- paste0(spec$run_prefix, group, sub("^[^_]+", "", design$sample))
    counts <- gen$kitTrxCounts(r, design)$counts
    dimnames(counts) <- list(sprintf("TRXG%05d", seq_len(spec$features)), runs)
    list(generator = r, scale = "counts", views = list(transcriptomics = counts),
        display = list(transcriptomics = sprintf("SYM%05d", seq_len(spec$features))),
        samples = integrationSamples(runs, group, design$replicate))
}

# --- Known effects ---

integrationOrthogonal <- function(x, basis) {
    r <- qr.resid(qr(basis), x)
    r / stats::sd(r)
}

integrationFactorScores <- function(recipe) {
    per <- recipe$biological_per_group
    group <- rep(recipe$groups, each = per)
    bio <- integrationBiologicalKey(group, rep(seq_len(per), times = length(recipe$groups)))
    z1 <- as.numeric(group == recipe$effect_group)
    draws <- withr::with_seed(recipe$seed, matrix(stats::rnorm(2L * length(bio)), ncol = 2L))
    z2 <- integrationOrthogonal(draws[, 1L], cbind(1, z1))
    z3 <- integrationOrthogonal(draws[, 2L], cbind(1, z1, z2))
    data.frame(biological_sample = bio, group = group, Z1 = z1, Z2 = z2, Z3 = z3,
        stringsAsFactors = FALSE)
}

integrationBatches <- function(scores, batches, seed) {
    withr::with_seed(seed, {
        out <- stats::setNames(character(nrow(scores)), scores$biological_sample)
        for (g in unique(scores$group)) {
            i <- which(scores$group == g)
            out[i] <- sample(rep(batches, length.out = length(i)))
        }
        out
    })
}

integrationViewLoadings <- function(recipe, p, omic, batches, seed) {
    withr::with_seed(seed, {
        factors <- vapply(recipe$factors, \(f) {
            l <- numeric(p)
            if (omic %in% f$active) {
                k <- round(p * f$share)
                l[sample.int(p, k)] <- sample(c(-1, 1), k, replace = TRUE) *
                    stats::runif(k, f$abs_loading[[1L]], f$abs_loading[[2L]])
            }
            l
        }, numeric(p))
        batch <- matrix(0, p, length(batches), dimnames = list(NULL, batches))
        affected <- sample.int(p, round(p * recipe$batch$share))
        batch[affected, ] <- stats::rnorm(length(affected) * length(batches), 0, recipe$batch$sd)
        list(factors = factors, batch = batch)
    })
}

integrationApplyEffects <- function(values, scale, samples, scores, loadings) {
    s <- as.matrix(scores[match(samples$biological_sample, scores$biological_sample),
        colnames(loadings$factors)])
    effect <- loadings$factors %*% t(s) + loadings$batch[, samples$batch, drop = FALSE]
    if (identical(scale, "counts")) {
        counts <- round(values * 2^effect)
        return(log2(sweep(counts + 0.5, 2L, colSums(counts) + 1, "/") * 1e6))
    }
    x <- values + effect
    run_median <- apply(x, 2L, stats::median, na.rm = TRUE)
    sweep(x, 2L, run_median - stats::median(run_median), "-")
}

integrationBuildOmic <- function(gen, recipe, omic, index, scores) {
    spec <- recipe$omics[[omic]]
    base <- switch(spec$generator,
        lfq = integrationBaseProteomics(gen, recipe, spec),
        trx_screen = integrationBaseTrx(gen, recipe, spec),
        integrationBaseMs(gen, recipe, spec, omic))
    batch_of <- integrationBatches(scores, spec$batches, recipe$seed + 100L + index)
    samples <- base$samples
    samples$batch <- unname(batch_of[samples$biological_sample])
    loadings <- lapply(seq_along(base$views), \(v) integrationViewLoadings(recipe,
        nrow(base$views[[v]]), omic, spec$batches, recipe$seed + 10L * index + v))
    names(loadings) <- names(base$views)
    views <- Map(\(values, l) integrationApplyEffects(values, base$scale, samples, scores, l),
        base$views, loadings)
    order <- withr::with_seed(recipe$seed + 200L + index, sample.int(nrow(samples)))
    dropped <- omic == recipe$dropped$omic &
        samples$biological_sample[order] == recipe$dropped$biological_sample
    keep <- order[!dropped]
    list(generator = base$generator, samples = samples[keep, , drop = FALSE],
        views = lapply(views, \(x) x[, keep, drop = FALSE]), display = base$display,
        loadings = loadings, batch_of = batch_of)
}

# --- Writers ---

integrationWriteTable <- function(table, path) {
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    path
}

# KIT-147: the package's contract (constants and Parquet writer), from its namespace (installed or
# loaded with pkgload).
integrationContract <- function(name) get(name, envir = asNamespace("MultiScholaR"))

integrationWriteParquet <- function(table, path) {
    integrationContract("writeIntegrationInputParquet")(table, path)
    path
}

integrationWriteOmic <- function(directory, omic, built, digits) {
    out <- file.path(directory, omic)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    columns <- c("Run", "group", "replicates", "tech_rep_group", "batch")
    text <- integrationContract("integrationInputTextTable")
    paths <- integrationWriteParquet(text(built$samples[, columns]), file.path(out, "samples.parquet"))
    for (view in names(built$views)) {
        values <- built$views[[view]]
        matrix_table <- cbind(data.frame(feature_id = rownames(values), stringsAsFactors = FALSE),
            as.data.frame(round(values, digits)))
        matrix_table <- integrationContract("integrationInputMatrixTable")(matrix_table, colnames(values))
        feature_table <- data.frame(feature_id = rownames(values),
            display_name = built$display[[view]], stringsAsFactors = FALSE)
        paths <- c(paths,
            integrationWriteParquet(matrix_table, file.path(out, sprintf("matrix_%s.parquet", view))),
            integrationWriteParquet(text(feature_table), file.path(out, sprintf("features_%s.parquet", view))))
    }
    c(paths, integrationWriteInputManifest(out, omic, built, paths))
}

# The omic's contract v2 manifest (fields of KIT-060 audit 4.2) over the tables just written.
integrationWriteInputManifest <- function(out, omic, built, paths) {
    layout <- integrationContract("INTEGRATION_INPUT_PARQUET")
    views <- lapply(names(built$views), \(view) list(name = view,
        matrix_file = sprintf("matrix_%s.parquet", view), features_file = sprintf("features_%s.parquet", view),
        rows = nrow(built$views[[view]]), columns = ncol(built$views[[view]]),
        row_groups = as.integer(ceiling(nrow(built$views[[view]]) / layout$row_group_rows))))
    manifest <- list(
        schema = integrationContract("INTEGRATION_INPUT_SCHEMA"),
        omic_type = omic,
        workflow = INTEGRATION_CORPUS_WORKFLOWS[[omic]],
        matrix_role = "normalised",
        correction = NULL,
        scale = INTEGRATION_CORPUS_SCALES[[omic]],
        source_state = list(name = "synthetic", synthetic = TRUE,
            note = "KIT-137 tier 1: written by the composer, not by an omic's journey"),
        sample_id_column = "Run",
        group_column = "group",
        technical_replicate_column = "tech_rep_group",
        biological_sample_column = "tech_rep_group",
        feature_id_column = "feature_id",
        display_name_column = "display_name",
        views = views,
        samples_file = "samples.parquet",
        parquet = layout,
        files = stats::setNames(lapply(paths, \(path) digest::digest(file = path, algo = "sha256")),
            basename(paths)),
        software = list(package = "MultiScholaR",
            composer = "tests/testdata/integration/compose_integration_corpus.R"),
        written_at = "2026-11-17T00:00:00Z",
        synthetic = TRUE)
    path <- file.path(out, "integration_input.json")
    jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
    path
}

integrationLoadingTable <- function(built, batch_levels) {
    rows <- do.call(c, lapply(names(built), \(omic) lapply(names(built[[omic]]$views), \(view) {
        l <- built[[omic]]$loadings[[view]]
        batch <- matrix(NA_real_, nrow(l$factors), length(batch_levels),
            dimnames = list(NULL, paste0("batch_", batch_levels)))
        batch[, paste0("batch_", colnames(l$batch))] <- l$batch
        cbind(data.frame(omic = omic, view = view, feature_id = rownames(built[[omic]]$views[[view]]),
            stringsAsFactors = FALSE), as.data.frame(l$factors), as.data.frame(batch))
    })))
    do.call(rbind, rows)
}

integrationWriteTruth <- function(directory, recipe, scores, built) {
    batch_levels <- sort(unique(unlist(lapply(recipe$omics, `[[`, "batches"))))
    rounded <- scores
    rounded[c("Z2", "Z3")] <- lapply(scores[c("Z2", "Z3")], round, 6L)
    loadings <- integrationLoadingTable(built, batch_levels)
    loadings[-(1:3)] <- lapply(loadings[-(1:3)], round, 6L)
    truth <- list(
        corpus = recipe$corpus, synthetic = TRUE, seed = recipe$seed,
        biological_samples = scores$biological_sample,
        views = lapply(built, \(b) names(b$views)),
        factors = list(
            Z1 = c(recipe$factors$Z1, list(scores = "group indicator, KO = 1")),
            Z2 = c(recipe$factors$Z2, list(scores = "N(0, 1), orthogonal in the sample to Z1")),
            Z3 = c(recipe$factors$Z3, list(scores = "N(0, 1), orthogonal in the sample to Z1, Z2"))),
        batches = lapply(built, \(b) as.list(b$batch_of)),
        dropped = recipe$dropped,
        runs = lapply(built, \(b) nrow(b$samples)),
        features = lapply(built, \(b) lapply(b$views, nrow)),
        missing_cells = lapply(built, \(b) lapply(b$views, \(x) sum(is.na(x)))),
        generators = lapply(built, `[[`, "generator"),
        recipe = recipe)
    truth_path <- file.path(directory, "truth.json")
    jsonlite::write_json(truth, truth_path, auto_unbox = TRUE, pretty = TRUE, digits = NA,
        null = "null")
    c(integrationWriteTable(rounded, file.path(directory, "truth_scores.tsv")),
        integrationWriteTable(loadings, file.path(directory, "truth_loadings.tsv")),
        truth_path)
}

integrationWriteManifest <- function(directory, recipe, files) {
    relative <- substring(files, nchar(directory) + 2L)
    manifest <- list(corpus = recipe$corpus, synthetic = TRUE, seed = recipe$seed,
        composer = "tests/testdata/integration/compose_integration_corpus.R",
        files = stats::setNames(lapply(files, \(path) list(bytes = file.size(path),
            sha256 = digest::digest(file = path, algo = "sha256"))), relative))
    jsonlite::write_json(manifest, file.path(directory, "manifest.json"), auto_unbox = TRUE,
        pretty = TRUE, digits = NA)
    manifest
}

# Compose the tier-1 integration corpus into `directory`; returns the manifest (a list).
composeIntegrationCorpus <- function(directory, seed = INTEGRATION_CORPUS_SEED, repo = ".") {
    recipe <- integrationCorpusRecipe(seed)
    gen <- integrationGenerators(repo)
    scores <- integrationFactorScores(recipe)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    omics <- names(recipe$omics)
    built <- lapply(seq_along(omics), \(i) integrationBuildOmic(gen, recipe, omics[[i]], i, scores))
    names(built) <- omics
    files <- unlist(lapply(omics, \(o) integrationWriteOmic(directory, o, built[[o]], recipe$digits)))
    files <- c(files, integrationWriteTruth(directory, recipe, scores, built))
    integrationWriteManifest(directory, recipe, files)
}

# --- The truth check, without a factor model ---

integrationAverageReplicates <- function(m, key) {
    levels <- unique(key)
    out <- vapply(levels, \(k) rowMeans(m[, key == k, drop = FALSE], na.rm = TRUE), numeric(nrow(m)))
    out[is.nan(out)] <- NA_real_
    out
}

# Per feature, the F-test p-value of adding z to the null design, on the observed samples.
integrationFeaturePvalues <- function(y, null, z) {
    full <- cbind(null, z)
    apply(y, 1L, \(v) {
        ok <- !is.na(v)
        q0 <- qr(null[ok, , drop = FALSE])
        q1 <- qr(full[ok, , drop = FALSE])
        df <- sum(ok) - q1$rank
        if (df < 3L || q1$rank == q0$rank) return(NA_real_)
        rss0 <- sum(qr.resid(q0, v[ok])^2)
        rss1 <- sum(qr.resid(q1, v[ok])^2)
        stats::pf((rss0 - rss1) / (rss1 / df), 1, df, lower.tail = FALSE)
    })
}

# Each feature adjusted for the group and the omic's batch, scaled to unit variance, missing set to
# the feature's mean (0).
integrationAdjustedScaled <- function(y, z1, batch) {
    x <- stats::model.matrix(~ z1 + batch, data.frame(z1 = z1, batch = batch))
    out <- t(apply(y, 1L, \(v) {
        ok <- !is.na(v)
        r <- numeric(length(v))
        r[ok] <- qr.resid(qr(x[ok, , drop = FALSE]), v[ok])
        spread <- if (sum(ok) > 1L) stats::sd(r[ok]) else 0
        # A feature observed too rarely to leave a residual carries no information here.
        if (spread > 0) r / spread else numeric(length(v))
    }))
    dimnames(out) <- dimnames(y)
    out
}

# Recover the known factors without a factor model. Returns the absolute correlation of Z2 with
# the first principal component of the concatenated adjusted, scaled views (biological samples
# present in every omic), and per view and factor the share of features associated with the
# factor at `p_threshold` (Z1 against the batch alone; Z2 and Z3 against Z1 and the batch).
checkIntegrationCorpusTruth <- function(directory, p_threshold = 0.001) {
    truth <- jsonlite::read_json(file.path(directory, "truth.json"), simplifyVector = TRUE)
    scores <- utils::read.delim(file.path(directory, "truth_scores.tsv"), stringsAsFactors = FALSE)
    factors <- names(truth$factors)
    adjusted <- list()
    association <- list()
    for (omic in names(truth$views)) {
        samples <- as.data.frame(arrow::read_parquet(file.path(directory, omic, "samples.parquet")))
        for (view in truth$views[[omic]]) {
            table <- as.data.frame(arrow::read_parquet(file.path(directory, omic,
                sprintf("matrix_%s.parquet", view))))
            m <- as.matrix(table[-1L])
            rownames(m) <- table$feature_id
            y <- integrationAverageReplicates(m, samples$tech_rep_group[match(colnames(m), samples$Run)])
            s <- scores[match(colnames(y), scores$biological_sample), ]
            batch <- samples$batch[match(colnames(y), samples$tech_rep_group)]
            batch_x <- stats::model.matrix(~batch, data.frame(batch = batch))
            key <- paste(omic, view, sep = "/")
            association[[key]] <- do.call(rbind, lapply(factors, \(k) {
                null <- if (k == "Z1") batch_x else cbind(batch_x, Z1 = s$Z1)
                p <- integrationFeaturePvalues(y, null, s[[k]])
                data.frame(view_key = key, factor = k, active = omic %in% truth$factors[[k]]$active,
                    features = nrow(y), fraction = mean(p < p_threshold, na.rm = TRUE),
                    stringsAsFactors = FALSE)
            }))
            adjusted[[key]] <- integrationAdjustedScaled(y, s$Z1, batch)
        }
    }
    shared <- Reduce(intersect, lapply(adjusted, colnames))
    stacked <- do.call(rbind, lapply(adjusted, \(x) x[, shared, drop = FALSE]))
    pc1 <- stats::prcomp(t(stacked), center = TRUE, scale. = FALSE, rank. = 1L)$x[, 1L]
    z2 <- scores$Z2[match(shared, scores$biological_sample)]
    association <- do.call(rbind, association)
    rownames(association) <- NULL
    list(z2_pc1_abs_cor = abs(stats::cor(pc1, z2)), shared_samples = length(shared),
        association = association)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    directory <- if (length(args)) args[[1L]] else file.path("tests", "testdata", "integration", "corpus")
    manifest <- composeIntegrationCorpus(directory)
    for (path in names(manifest$files)) {
        cat(sprintf("%-40s %8d bytes  %s\n", path, as.integer(manifest$files[[path]]$bytes),
            substr(manifest$files[[path]]$sha256, 1L, 12L)))
    }
    check <- checkIntegrationCorpusTruth(directory)
    cat(sprintf("Z2 against PC1: |r| = %.4f on %d shared biological samples\n",
        check$z2_pc1_abs_cor, check$shared_samples))
    print(check$association, row.names = FALSE)
}
