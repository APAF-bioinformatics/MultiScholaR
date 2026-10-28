# KIT-068: reproducible lipidomics and metabolomics corpora in the vendor layouts the GUI imports (the
# kit_protein_corpus.R pattern: a recipe, withr::with_seed, a manifest with per-file sha256).
# - lipid_screen: two LipidSearch exports (LCMS_Pos, LCMS_Neg; tab-separated .txt as the e2e fixture) with the
#   columns importLipidSearchData reads (LipidName, LipidClass, FattyAcid, IonType, BaseRt, CalcMz, Grade) and the
#   sample columns; internal-standard rows whose LipidClass is "ISTD" and LipidName starts with "IS_" (the
#   lipidsearch is_pattern `^IS_|ISTD|Internal|_IS` is applied to the annotation column by the ITSD
#   normalisation); duplicate LipidName rows under a second adduct (the duplicate step has work to do).
# - metab_screen: two MS-DIAL exports (LCMS_Pos, GCMS; tab-separated .txt) with the columns
#   importMetabMSDIALData reads (Alignment ID, Average Rt(min), Average Mz, Metabolite name, Adduct type) and the
#   sample columns; internal standards named "IS_<name>_d<k>" (the msdial is_pattern `^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:`).
#   Alignment IDs are unique, as MS-DIAL writes them, so the duplicate step resolves nothing on this corpus (recorded).
# Shared design: two groups (WT, KO), twelve biological samples each measured twice (technical-replicate pairs),
# three batches balanced over groups, a declared differential fraction and a batch effect so that RUV-III has
# something to find. No missing values (the between-sample normalisation runs dense). design.tsv carries
# sample, group, batch, tech_rep_group (the biological sample) and replicates (the biological replicate number
# within the group, the column the kit driver's technical-replicate assignment reads, as the DIA corpora name it).
# Sample names are already clean (letters, digits, underscores), so the import runs with sanitize_names off.
# Usage: Rscript tools/profiling/kit_omics_corpus.R <lipid_screen|metab_screen> <directory> [biological samples per group] [features per assay]

kitOmicsRecipe <- function(kind, biological_per_group = 6L, features = NULL) {
    stopifnot(kind %in% c("lipid_screen", "metab_screen"))
    seed <- switch(kind, lipid_screen = 20261029L, metab_screen = 20261030L)
    features <- as.integer(features %||% switch(kind, lipid_screen = 200L, metab_screen = 300L))
    list(kind = kind, synthetic = TRUE, seed = seed, biological_per_group = as.integer(biological_per_group),
        technical_replicates = 2L, groups = c("WT", "KO"), batches = c("B1", "B2", "B3"),
        features_per_assay = features, standards_per_assay = 6L, duplicate_fraction = if (kind == "lipid_screen") 0.05 else 0,
        differential_fraction = 0.15, batch_effect_sd = 0.35, sample_scale_sd = 0.3, technical_noise_sd = 0.12,
        residual_sd = 0.3, assays = switch(kind, lipid_screen = c("LCMS_Pos", "LCMS_Neg"), metab_screen = c("LCMS_Pos", "GCMS")))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# The design: biological samples alternate groups; technical replicates are consecutive; batches cycle over the
# biological samples so each batch holds both groups.
kitOmicsDesign <- function(recipe) {
    biological <- expand.grid(replicate = seq_len(recipe$biological_per_group), group = recipe$groups, stringsAsFactors = FALSE)
    biological <- biological[order(biological$replicate, biological$group), , drop = FALSE]
    biological$tech_rep_group <- sprintf("%s_%02d", biological$group, biological$replicate)
    biological$batch <- recipe$batches[((seq_len(nrow(biological)) - 1L) %% length(recipe$batches)) + 1L]
    design <- biological[rep(seq_len(nrow(biological)), each = recipe$technical_replicates), , drop = FALSE]
    design$sample <- sprintf("%s_r%d", design$tech_rep_group, rep(seq_len(recipe$technical_replicates), times = nrow(biological)))
    rownames(design) <- NULL
    design[, c("sample", "group", "batch", "tech_rep_group", "replicate")]
}

# One assay's intensity matrix (features x samples, raw scale) and its standards (standards x samples).
kitOmicsAssayMatrix <- function(recipe, design, assay_seed, base_mean) {
    withr::with_seed(assay_seed, {
        p <- recipe$features_per_assay; n <- nrow(design)
        base <- stats::rnorm(p, mean = base_mean, sd = 1.8)
        differential <- seq_len(p) <= round(p * recipe$differential_fraction)
        effect <- ifelse(differential, stats::rnorm(p, 0, 1.1), 0)
        batch_effect <- stats::rnorm(length(recipe$batches), 0, recipe$batch_effect_sd); names(batch_effect) <- recipe$batches
        # The batch effect is feature-specific in sign (a loading), as a real batch effect is.
        batch_loading <- stats::rnorm(p, 1, 0.3)
        biological_noise <- matrix(stats::rnorm(p * length(unique(design$tech_rep_group)), 0, recipe$residual_sd), nrow = p)
        colnames(biological_noise) <- unique(design$tech_rep_group)
        sample_scale <- stats::rnorm(n, 0, recipe$sample_scale_sd)
        log2 <- outer(base, rep(0, n), "+") + outer(effect, as.numeric(design$group == "KO")) +
            outer(batch_loading, batch_effect[design$batch]) + biological_noise[, design$tech_rep_group, drop = FALSE] +
            matrix(sample_scale, nrow = p, ncol = n, byrow = TRUE) + matrix(stats::rnorm(p * n, 0, recipe$technical_noise_sd), nrow = p)
        standards_base <- stats::rnorm(recipe$standards_per_assay, base_mean + 1, 0.5)
        standards <- outer(standards_base, rep(0, n), "+") + matrix(sample_scale, nrow = recipe$standards_per_assay, ncol = n, byrow = TRUE) +
            matrix(stats::rnorm(recipe$standards_per_assay * n, 0, 0.05), nrow = recipe$standards_per_assay)
        list(values = round(2^log2, 1), standards = round(2^standards, 1), differential = differential, sample_scale = sample_scale)
    })
}

kitWriteOmicsManifest <- function(directory, recipe, files, rows, design, notes) {
    manifest <- list(profile = recipe$kind, synthetic = TRUE, samples = nrow(design), rows = as.integer(rows),
        biological_samples = length(unique(design$tech_rep_group)), technical_replicates = recipe$technical_replicates,
        batches = length(recipe$batches), assays = recipe$assays, recipe = recipe, expected_retained_samples = nrow(design),
        notes = notes, files = stats::setNames(lapply(files, function(path) list(bytes = file.size(path),
            sha256 = digest::digest(file = path, algo = "sha256"))), basename(files)))
    jsonlite::write_json(manifest, file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    manifest
}

writeKitLipidCorpus <- function(directory, recipe = kitOmicsRecipe("lipid_screen")) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    design <- kitOmicsDesign(recipe)
    classes <- list(LCMS_Pos = c("PC", "PE", "TG", "SM", "CE"), LCMS_Neg = c("PI", "PE", "PS", "FA", "PG"))
    adducts <- list(LCMS_Pos = c("[M+H]+", "[M+NH4]+"), LCMS_Neg = c("[M-H]-", "[M+HCOO]-"))
    files <- character(); rows <- 0L
    for (i in seq_along(recipe$assays)) {
        assay <- recipe$assays[[i]]
        m <- kitOmicsAssayMatrix(recipe, design, recipe$seed + i, base_mean = 16)
        p <- recipe$features_per_assay
        withr::with_seed(recipe$seed + 10L * i, {
            class <- sample(classes[[assay]], p, replace = TRUE)
            carbons <- sample(30:40, p, replace = TRUE); doubles <- sample(0:6, p, replace = TRUE)
            name <- sprintf("%s %d:%d", class, carbons, doubles)
            # Distinct lipid species: the same class and total composition may repeat; a repeat is a duplicate for the
            # duplicate step only when the LipidName repeats, which the fraction below controls explicitly.
            name <- make.unique(name, sep = "_")
            n_dup <- round(p * recipe$duplicate_fraction)
            dup_rows <- if (n_dup > 0L) sample(seq_len(p), n_dup) else integer()
            fatty <- sprintf("%d:%d/%d:%d", carbons %/% 2, doubles %/% 2, carbons - carbons %/% 2, doubles - doubles %/% 2)
            rt <- round(stats::runif(p, 1, 18), 2); mz <- round(stats::runif(p, 400, 1000), 3)
        })
        table <- data.frame(LipidName = name, LipidClass = class, FattyAcid = fatty, IonType = adducts[[assay]][[1L]],
            BaseRt = rt, CalcMz = mz, Grade = "A", stringsAsFactors = FALSE, check.names = FALSE)
        values <- as.data.frame(m$values); names(values) <- design$sample
        table <- cbind(table, values)
        if (length(dup_rows)) {
            # The duplicate carries the second adduct and a lower intensity, so the resolution keeps the first row.
            duplicates <- table[dup_rows, , drop = FALSE]
            duplicates$IonType <- adducts[[assay]][[2L]]
            duplicates$CalcMz <- round(duplicates$CalcMz + 17.027, 3)
            duplicates[design$sample] <- round(duplicates[design$sample] * 0.6, 1)
            table <- rbind(table, duplicates)
        }
        standards <- data.frame(LipidName = sprintf("IS_%s_d7", sprintf("%s %d:%d", classes[[assay]][seq_len(recipe$standards_per_assay) %% length(classes[[assay]]) + 1L],
            34L + seq_len(recipe$standards_per_assay), 1L)), LipidClass = "ISTD", FattyAcid = "", IonType = adducts[[assay]][[1L]],
            BaseRt = round(seq(2, 16, length.out = recipe$standards_per_assay), 2), CalcMz = round(seq(500, 900, length.out = recipe$standards_per_assay), 3),
            Grade = "A", stringsAsFactors = FALSE, check.names = FALSE)
        standard_values <- as.data.frame(m$standards); names(standard_values) <- design$sample
        standards <- cbind(standards, standard_values)
        table <- rbind(table, standards)
        table <- withr::with_seed(recipe$seed + 100L * i, table[sample.int(nrow(table)), , drop = FALSE])
        rownames(table) <- NULL
        path <- file.path(directory, sprintf("lipidsearch_%s.txt", tolower(assay)))
        utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        files <- c(files, path); rows <- rows + nrow(table)
    }
    utils::write.table(design[, c("sample", "group", "batch", "tech_rep_group")] |> transform(replicates = design$replicate),
        file.path(directory, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    kitWriteOmicsManifest(directory, recipe, c(files, file.path(directory, "design.tsv")), rows, design,
        notes = sprintf("%d lipids per assay in two LipidSearch exports (%s), %d duplicate LipidName rows under a second adduct at 0.6 x intensity, %d internal-standard rows per assay (LipidClass ISTD, LipidName IS_*), dense (no missing values); %d biological samples x %d technical replicates in %d batches; %d%% of lipids differential (KO vs WT)",
            recipe$features_per_assay, paste(recipe$assays, collapse = ", "), round(recipe$features_per_assay * recipe$duplicate_fraction),
            recipe$standards_per_assay, length(unique(design$tech_rep_group)), recipe$technical_replicates, length(recipe$batches),
            round(100 * recipe$differential_fraction)))
}

writeKitMetabCorpus <- function(directory, recipe = kitOmicsRecipe("metab_screen")) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    design <- kitOmicsDesign(recipe)
    files <- character(); rows <- 0L
    for (i in seq_along(recipe$assays)) {
        assay <- recipe$assays[[i]]
        m <- kitOmicsAssayMatrix(recipe, design, recipe$seed + i, base_mean = 15)
        p <- recipe$features_per_assay
        withr::with_seed(recipe$seed + 10L * i, {
            rt <- round(stats::runif(p, 0.5, 25), 2); mz <- round(stats::runif(p, 80, 900), 4)
            named <- stats::runif(p) < 0.6
        })
        name <- ifelse(named, sprintf("Metabolite_%s_%03d", assay, seq_len(p)), "Unknown")
        table <- if (assay == "GCMS") data.frame(`Alignment ID` = seq_len(p) - 1L, `Average Rt(min)` = rt, `Metabolite name` = name,
            stringsAsFactors = FALSE, check.names = FALSE) else data.frame(`Alignment ID` = seq_len(p) - 1L, `Average Rt(min)` = rt,
            `Average Mz` = mz, `Metabolite name` = name, `Adduct type` = "[M+H]+", stringsAsFactors = FALSE, check.names = FALSE)
        values <- as.data.frame(m$values); names(values) <- design$sample
        table <- cbind(table, values)
        standards <- table[seq_len(recipe$standards_per_assay), , drop = FALSE]
        standards$`Alignment ID` <- p + seq_len(recipe$standards_per_assay) - 1L
        standards$`Metabolite name` <- sprintf("IS_Standard_%s_%d_d%d", assay, seq_len(recipe$standards_per_assay), 3L + seq_len(recipe$standards_per_assay))
        standard_values <- as.data.frame(m$standards); names(standard_values) <- design$sample
        standards[design$sample] <- standard_values
        table <- rbind(table, standards)
        path <- file.path(directory, sprintf("msdial_%s.txt", tolower(assay)))
        utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        files <- c(files, path); rows <- rows + nrow(table)
    }
    utils::write.table(design[, c("sample", "group", "batch", "tech_rep_group")] |> transform(replicates = design$replicate),
        file.path(directory, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    kitWriteOmicsManifest(directory, recipe, c(files, file.path(directory, "design.tsv")), rows, design,
        notes = sprintf("%d features per assay in two MS-DIAL exports (%s; the GC export without m/z and adduct columns), %d internal-standard rows per assay (Metabolite name IS_*_d<k>), unique Alignment IDs (the duplicate step resolves nothing), dense (no missing values); %d biological samples x %d technical replicates in %d batches; %d%% of features differential (KO vs WT)",
            recipe$features_per_assay, paste(recipe$assays, collapse = ", "), recipe$standards_per_assay,
            length(unique(design$tech_rep_group)), recipe$technical_replicates, length(recipe$batches), round(100 * recipe$differential_fraction)))
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    if (length(args) < 2L) stop("usage: kit_omics_corpus.R <lipid_screen|metab_screen> <directory> [biological samples per group] [features per assay]")
    kind <- args[[1L]]; directory <- args[[2L]]
    recipe <- kitOmicsRecipe(kind, if (length(args) > 2L) as.integer(args[[3L]]) else 6L, if (length(args) > 3L) as.integer(args[[4L]]) else NULL)
    manifest <- switch(kind, lipid_screen = writeKitLipidCorpus(directory, recipe), metab_screen = writeKitMetabCorpus(directory, recipe))
    cat(kind, directory, "samples", manifest$samples, "rows", manifest$rows, "\n")
    for (name in names(manifest$files)) cat("  ", name, manifest$files[[name]]$bytes, "bytes", substr(manifest$files[[name]]$sha256, 1, 12), "\n")
}
