# KIT-021: reproducible protein-level corpora in the real tool output formats (no production activation).
# - LFQ: a MaxQuant `proteinGroups.txt` (dotted headers as read by importMaxQuantData: Protein.IDs,
#   Majority.protein.IDs, Gene.names, LFQ.intensity.<sample>, Potential.contaminant, Reverse).
# - TMT: a Proteome Discoverer protein export (Accession, Description, Abundance: F<plex>: <channel>, <sample>),
#   TMT10 plexes, as read by importProteomeDiscovererTMTData (runs named <channel>_<sample>).
# Each corpus directory gets design.tsv (sample, group, batch), proteins.fasta and manifest.json with sha256 sums.
# The design names the runs as the import sanitises them (janitor::make_clean_names, the GUI's default): lowercase,
# non-alphanumerics to "_", an "x" before a leading digit (a TMT run "126_WT_001" is "x126_wt_001").
# Usage: Rscript tools/profiling/kit_protein_corpus.R <lfq|tmt|lfq_fixture|tmt_fixture> <directory> [samples] [proteins]

kitProteinRecipe <- function(kind, samples = 300L, proteins = 8000L) {
    seed <- switch(kind, lfq = 20261020L, tmt = 20261021L, stop("kind: lfq | tmt"))
    list(kind = kind, synthetic = TRUE, seed = seed, samples = as.integer(samples), proteins = as.integer(proteins),
        groups = c("WT", "KO"), differential_fraction = 0.1, missing_fraction = if (kind == "lfq") 0.08 else 0.02,
        flagged_rows = if (kind == "lfq") 50L else 0L, channels = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131"))
}

kitProteinMatrix <- function(recipe) {
    withr::with_seed(recipe$seed, {
        n <- recipe$samples; p <- recipe$proteins
        group <- rep(recipe$groups, length.out = n)
        batch <- paste0("B", ((seq_len(n) - 1L) %/% 10L) + 1L)
        base <- stats::rnorm(p, mean = 20, sd = 2.2)
        differential <- seq_len(p) <= round(p * recipe$differential_fraction)
        effect <- ifelse(differential, stats::rnorm(p, 0, 1.2), 0)
        batch_effect <- stats::rnorm(length(unique(batch)), 0, 0.15); names(batch_effect) <- unique(batch)
        sample_scale <- stats::rnorm(n, 0, 0.25)
        log2 <- outer(base, rep(0, n), "+") + outer(effect, as.numeric(group == "KO")) +
            matrix(batch_effect[batch], nrow = p, ncol = n, byrow = TRUE) + matrix(sample_scale, nrow = p, ncol = n, byrow = TRUE) +
            matrix(stats::rnorm(p * n, 0, 0.4), nrow = p, ncol = n)
        values <- round(2^log2, 1)
        missing <- matrix(stats::runif(p * n) < recipe$missing_fraction, nrow = p, ncol = n)
        list(values = values, missing = missing, group = group, batch = batch, differential = differential)
    })
}

kitWriteFasta <- function(path, ids, genes) {
    lines <- character(0)
    for (i in seq_along(ids)) {
        seq <- paste(rep(c("MPEPAKLVRAAAGSTDEFHIKLMNPQRSTVWY"), 2L), collapse = "")
        lines <- c(lines, sprintf(">sp|%s|KIT_%s Synthetic protein %d OS=Homo sapiens OX=9606 GN=%s PE=1 SV=1", ids[i], ids[i], i, genes[i]), seq)
    }
    writeLines(lines, path)
}

kitWriteManifest <- function(directory, recipe, files, rows, samples, notes = NULL) {
    manifest <- list(profile = recipe$kind, synthetic = TRUE, samples = as.integer(samples), rows = as.integer(rows),
        recipe = recipe[setdiff(names(recipe), "channels")], expected_retained_samples = as.integer(samples), notes = notes,
        files = stats::setNames(lapply(files, function(path) list(bytes = file.size(path),
            sha256 = digest::digest(file = path, algo = "sha256"))), basename(files)))
    jsonlite::write_json(manifest, file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    manifest
}

writeKitLfqCorpus <- function(directory, recipe = kitProteinRecipe("lfq")) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    m <- kitProteinMatrix(recipe)
    p <- recipe$proteins; n <- recipe$samples
    samples <- sprintf("%s_%03d", m$group, stats::ave(seq_len(n), m$group, FUN = seq_along))
    ids <- sprintf("P%05d", seq_len(p)); genes <- sprintf("GENE%d", seq_len(p))
    values <- m$values; values[m$missing] <- 0   # MaxQuant writes 0 for a missing LFQ intensity
    table <- data.frame(Protein.IDs = ids, Majority.protein.IDs = ids, Gene.names = genes, stringsAsFactors = FALSE, check.names = FALSE)
    intensities <- as.data.frame(values); names(intensities) <- paste0("LFQ.intensity.", samples)
    table <- cbind(table, intensities)
    table$Potential.contaminant <- ""; table$Reverse <- ""
    # Flagged rows the reader drops (contaminants and decoys), appended as MaxQuant lists them.
    if (recipe$flagged_rows > 0L) {
        flagged <- table[seq_len(recipe$flagged_rows), ]
        flagged$Protein.IDs <- flagged$Majority.protein.IDs <- paste0("CON__", flagged$Protein.IDs)
        flagged$Potential.contaminant <- "+"
        decoys <- table[seq_len(recipe$flagged_rows) + recipe$flagged_rows, ]
        decoys$Protein.IDs <- decoys$Majority.protein.IDs <- paste0("REV__", decoys$Protein.IDs)
        decoys$Reverse <- "+"
        table <- rbind(table, flagged, decoys)
    }
    report <- file.path(directory, "proteinGroups.txt")
    utils::write.table(table, report, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    design <- data.frame(sample = janitor::make_clean_names(samples), group = m$group, batch = m$batch, stringsAsFactors = FALSE)
    utils::write.table(design, file.path(directory, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    kitWriteFasta(file.path(directory, "proteins.fasta"), ids, genes)
    kitWriteManifest(directory, recipe, file.path(directory, c("proteinGroups.txt", "design.tsv", "proteins.fasta")), nrow(table), n,
        notes = sprintf("%d protein groups + %d contaminant and %d decoy rows (dropped by the reader); %.0f%% missing written as 0; %d%% of proteins differential (KO vs WT)",
            p, recipe$flagged_rows, recipe$flagged_rows, 100 * recipe$missing_fraction, round(100 * recipe$differential_fraction)))
}

writeKitTmtCorpus <- function(directory, recipe = kitProteinRecipe("tmt")) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    stopifnot(recipe$samples %% length(recipe$channels) == 0L)
    m <- kitProteinMatrix(recipe)
    p <- recipe$proteins; n <- recipe$samples
    plex <- ((seq_len(n) - 1L) %/% length(recipe$channels)) + 1L
    channel <- recipe$channels[((seq_len(n) - 1L) %% length(recipe$channels)) + 1L]
    labels <- sprintf("%s_%03d", m$group, stats::ave(seq_len(n), m$group, FUN = seq_along))
    ids <- sprintf("P%05d", seq_len(p)); genes <- sprintf("GENE%d", seq_len(p))
    values <- m$values; values[m$missing] <- NA   # Proteome Discoverer leaves a missing abundance empty
    table <- data.frame(Accession = ids, Description = sprintf("Synthetic protein %d OS=Homo sapiens GN=%s", seq_len(p), genes), stringsAsFactors = FALSE, check.names = FALSE)
    abundances <- as.data.frame(values)
    names(abundances) <- sprintf("Abundance: F%d: %s, %s", plex, channel, labels)
    table <- cbind(table, abundances)
    report <- file.path(directory, "pd_proteins.tsv")
    utils::write.table(table, report, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    # The reader names each run <channel>_<sample label>; the design uses those names and the plex as the batch.
    design <- data.frame(sample = janitor::make_clean_names(paste0(channel, "_", labels)), group = m$group, batch = paste0("F", plex), stringsAsFactors = FALSE)
    utils::write.table(design, file.path(directory, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    kitWriteFasta(file.path(directory, "proteins.fasta"), ids, genes)
    kitWriteManifest(directory, recipe, file.path(directory, c("pd_proteins.tsv", "design.tsv", "proteins.fasta")), nrow(table), n,
        notes = sprintf("%d proteins across %d TMT10 plexes (%d channels each); %.0f%% missing written empty; %d%% of proteins differential (KO vs WT); runs are <channel>_<sample>",
            p, max(plex), length(recipe$channels), 100 * recipe$missing_fraction, round(100 * recipe$differential_fraction)))
}

# The e2e fixtures as bound corpora (the screen lifecycle members): copied with a manifest.
writeKitFixtureCorpus <- function(directory, kind, repo = normalizePath(".")) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    source <- file.path(repo, "tests/testdata/e2e", if (kind == "lfq") "prot_lfq" else "prot_tmt")
    report <- if (kind == "lfq") "proteinGroups.txt" else "pd_tmt_peptides.tsv"
    for (name in c(report, "proteins.fasta")) file.copy(file.path(source, name), file.path(directory, name), overwrite = TRUE)
    design <- utils::read.delim(file.path(source, "design.tsv"), stringsAsFactors = FALSE)
    if (kind == "tmt") {
        # The fixture's design names samples WT_1; the reader names runs <channel>_<sample>.
        header <- names(utils::read.delim(file.path(source, report), check.names = FALSE, nrows = 1L))
        runs <- sub("Abundance: F[0-9]+: ([0-9]+[A-Z]?), (.+)", "\\1_\\2", grep("^Abundance: ", header, value = TRUE))
        design$sample <- vapply(design$sample, function(s) runs[endsWith(runs, paste0("_", s))][[1L]], "")
    }
    design$sample <- janitor::make_clean_names(design$sample)
    utils::write.table(design, file.path(directory, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    recipe <- list(kind = paste0(kind, "_fixture"), synthetic = TRUE, seed = 20260428L, samples = nrow(design), proteins = 5L, source = source)
    kitWriteManifest(directory, recipe, file.path(directory, c(report, "design.tsv", "proteins.fasta")), 5L, nrow(design),
        notes = "the e2e fixture, copied; the screen lifecycle corpus for the protein route")
}

if (sys.nframe() == 0L) {
    args <- commandArgs(TRUE)
    if (length(args) < 2L) stop("usage: kit_protein_corpus.R <lfq|tmt|lfq_fixture|tmt_fixture> <directory> [samples] [proteins]")
    kind <- args[[1L]]; directory <- args[[2L]]
    manifest <- switch(kind,
        lfq = writeKitLfqCorpus(directory, kitProteinRecipe("lfq", if (length(args) > 2L) as.integer(args[[3L]]) else 300L, if (length(args) > 3L) as.integer(args[[4L]]) else 8000L)),
        tmt = writeKitTmtCorpus(directory, kitProteinRecipe("tmt", if (length(args) > 2L) as.integer(args[[3L]]) else 300L, if (length(args) > 3L) as.integer(args[[4L]]) else 8000L)),
        lfq_fixture = writeKitFixtureCorpus(directory, "lfq"), tmt_fixture = writeKitFixtureCorpus(directory, "tmt"),
        stop("kind: lfq | tmt | lfq_fixture | tmt_fixture"))
    cat(kind, directory, "samples", manifest$samples, "rows", manifest$rows, "\n")
    for (name in names(manifest$files)) cat("  ", name, manifest$files[[name]]$bytes, "bytes", substr(manifest$files[[name]]$sha256, 1, 12), "\n")
}
