# E2E Fixture Generation Script
#
# DEVELOPER TOOL — not run in CI.
#
# Writes every generated e2e fixture from `manifest.json` and then checks every committed fixture file against the
# sha256 the manifest records for it (KIT-101 F-06: the generator reproduces the committed fixtures byte for byte, or
# the manifest says which are hand-made and why). It stops with an error if any file differs.
#
# - Proteomics lanes (KIT-101 revision 2.2, F-07): prot_dia, prot_dia_limpa (DIA-NN report.tsv), prot_tmt (Proteome
#   Discoverer TMT), prot_lfq (MaxQuant proteinGroups.txt) and prot_lfq_fragpipe (FragPipe combined_protein.tsv). One
#   synthetic protein universe (400 accessions, tryptic peptides, a matching FASTA) shared by the five lanes; each
#   lane's protein abundances come from the KIT-021 protein corpus recipe (tools/profiling/kit_protein_corpus.R:
#   kitProteinRecipe/kitProteinMatrix, 10% of proteins differential KO vs WT) at the lane's size and seed
#   (`fixture_recipe` in the manifest). The DIA-NN lanes expand each protein into precursors with every column the
#   import, the Q-value filter and the precursor rollup read.
# - trx_screen: the KIT-040 corpus recipe (tools/profiling/kit_omics_corpus.R: kitOmicsRecipe("trx_screen", 3, 300),
#   withr::with_seed(20261106)), written by the corpus tool itself.
# - Hand-made lanes (`fixture_recipe$generator == "none"`: lipid_canonical, metab_lc, metab_gc, metab_combined): not
#   written; their committed files are checked only. The manifest records why.
# - Retained fixtures (`retained_fixtures` in the manifest): the earlier six-sample proteomics fixtures in prot_lfq/ and
#   prot_tmt/, still read by tests/testthat/test-workflow-native-session.R and by the KIT-021 lifecycle members
#   (tools/profiling/kit_protein_corpus.R writeKitFixtureCorpus), are rewritten by their original writers below; the
#   hand-made seed_* files are checked only.
#
# Needs jsonlite, withr and digest. Run from the repository root:
#   Rscript tests/testdata/e2e/generate_fixtures.R
#
# Idempotent: running twice produces byte-identical output files.

set.seed(20260428L)

# ---------------------------------------------------------------------------
# Bootstrap: resolve paths without testthat context
# ---------------------------------------------------------------------------

.script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0L) {
        dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
        # Running via source() — use sys.frame to find caller file
        src <- tryCatch(
            normalizePath(sys.frame(1L)$ofile)
            , error = function(e) NULL
        )
        if (!is.null(src)) {
            dirname(src)
        } else {
            # Last resort: assume CWD is project root
            file.path(getwd(), "tests", "testdata", "e2e")
        }
    }
}

FIXTURE_ROOT <- .script_dir()
cat("[generate_fixtures] fixture root:", FIXTURE_ROOT, "\n")

# ---------------------------------------------------------------------------
# Helper: read the manifest directly (no testthat dependency)
# ---------------------------------------------------------------------------

.read_manifest <- function() {
    manifest_path <- file.path(FIXTURE_ROOT, "manifest.json")
    if (!file.exists(manifest_path)) {
        stop("manifest.json not found at: ", manifest_path)
    }
    jsonlite::read_json(manifest_path, simplifyVector = FALSE)
}

# ---------------------------------------------------------------------------
# Helper: write TSV with consistent format (no row.names, tab separator)
# ---------------------------------------------------------------------------

.write_tsv <- function(df, path) {
    utils::write.table(
        df
        , file = path
        , sep = "\t"
        , row.names = FALSE
        , col.names = TRUE
        , quote = FALSE
        , na = ""
    )
    invisible(path)
}

# ---------------------------------------------------------------------------
# Proteomics lanes (KIT-101 revision 2.2)
# ---------------------------------------------------------------------------

#' The KIT-021 protein corpus tool (kitProteinRecipe, kitProteinMatrix), read from this source tree.
#' @return Environment holding the tool's functions.
#' @noRd
kitProteinTools <- function() {
    tool <- file.path(FIXTURE_ROOT, "..", "..", "..", "tools", "profiling", "kit_protein_corpus.R")
    if (!file.exists(tool)) stop("the KIT-021 protein corpus tool is not in this source tree: ", tool)
    env <- new.env(parent = globalenv())
    sys.source(tool, env)
    env
}

#' The protein universe shared by the five proteomics lanes: accessions, gene and entry names, the peptides a
#' DIA-NN report observes (1 to 4 per protein, unique, tryptic) and full sequences (the observed peptides plus two
#' unobserved ones). No random draw outside `seed`.
#' @param proteins Number of proteins.
#' @param seed Seed of the universe.
#' @return List of ids, genes, entries, names, observed (list of peptides) and sequences.
#' @noRd
protFixtureUniverse <- function(proteins, seed) {
    withr::with_seed(seed, {
        residues <- strsplit("ACDEFGHILMNQSTVWY", "", fixed = TRUE)[[1L]]
        observed <- sample(1:4, proteins, replace = TRUE, prob = c(0.15, 0.3, 0.3, 0.25))
        per_protein <- observed + 2L
        candidates <- vapply(seq_len(2L * sum(per_protein)), \(i) {
            body <- paste(sample(residues, sample(6:13, 1L), replace = TRUE), collapse = "")
            paste0(body, sample(c("K", "R"), 1L))
        }, character(1))
        peptides <- unique(candidates)[seq_len(sum(per_protein))]
        by_protein <- split(peptides, rep(seq_len(proteins), per_protein))
        list(
            ids = sprintf("P%05d", seq_len(proteins)),
            genes = sprintf("GENE%d", seq_len(proteins)),
            entries = sprintf("SYN%d_HUMAN", seq_len(proteins)),
            names = sprintf("Synthetic protein %d", seq_len(proteins)),
            observed = lapply(seq_len(proteins), \(i) by_protein[[i]][seq_len(observed[[i]])]),
            sequences = vapply(by_protein, \(p) paste0("M", paste(p, collapse = "")), character(1), USE.NAMES = FALSE)
        )
    })
}

#' The lane's protein abundance matrix from the KIT-021 recipe, with its columns in design order (WT_1..WT_n,
#' KO_1..KO_n).
#' @param lane Manifest lane (its `fixture_recipe` gives kit_recipe, seed, biological_per_group, proteins).
#' @param kit The KIT-021 tool environment.
#' @return List of design, values (linear, proteins x samples), missing (logical), differential, recipe.
#' @noRd
protFixtureMatrix <- function(lane, kit) {
    shape <- lane$fixture_recipe
    per_group <- as.integer(shape$biological_per_group)
    recipe <- kit$kitProteinRecipe(shape$kit_recipe, samples = 2L * per_group, proteins = as.integer(shape$proteins))
    recipe$seed <- as.integer(shape$seed)
    if (!is.null(shape$flagged_rows)) recipe$flagged_rows <- as.integer(shape$flagged_rows)
    m <- kit$kitProteinMatrix(recipe)
    label <- sprintf("%s_%d", m$group, stats::ave(seq_along(m$group), m$group, FUN = seq_along))
    groups <- unlist(lane$groups)
    design <- data.frame(
        sample = sprintf("%s_%d", rep(groups, each = per_group), rep(seq_len(per_group), length(groups))),
        group = rep(groups, each = per_group),
        stringsAsFactors = FALSE
    )
    at <- match(design$sample, label)
    if (anyNA(at)) stop("lane groups differ from the KIT-021 recipe groups: ", lane$lane_id)
    list(design = design, values = m$values[, at, drop = FALSE], missing = m$missing[, at, drop = FALSE],
        differential = m$differential, recipe = recipe)
}

#' Write the universe as a UniProt-style FASTA.
#' @noRd
writeProtFixtureFasta <- function(path, universe) {
    headers <- sprintf(">sp|%s|%s %s OS=Homo sapiens OX=9606 GN=%s PE=1 SV=1",
        universe$ids, universe$entries, universe$names, universe$genes)
    writeLines(as.vector(rbind(headers, universe$sequences)), path, useBytes = TRUE)
    invisible(path)
}

#' Write the lane's design.tsv and FASTA, the files every proteomics lane has.
#' @noRd
writeProtFixtureCommon <- function(lane, lane_dir, universe, design) {
    .write_tsv(design, file.path(lane_dir, "design.tsv"))
    writeProtFixtureFasta(file.path(lane_dir, lane$fasta_file), universe)
    cat("  wrote design.tsv,", lane$fasta_file, "\n")
}

#' A proteomics lane's inputs: the universe, the abundance matrix and the KIT-021 tool.
#' @noRd
protFixtureInputs <- function(lane) {
    kit <- kitProteinTools()
    shape <- lane$fixture_recipe
    list(
        universe = protFixtureUniverse(as.integer(shape$proteins), as.integer(shape$universe_seed)),
        matrix = protFixtureMatrix(lane, kit)
    )
}

#' DIA-NN main report (long format): one row per precursor and run, the columns the import (Protein.Group,
#' Stripped.Sequence, Run, Precursor.Normalised), the Q-value filter (Q.Value, Global.Q.Value, Global.PG.Q.Value,
#' Proteotypic) and the precursor rollup (Protein.Ids, Modified.Sequence, Precursor.Id, Precursor.Charge,
#' Precursor.Quantity, Precursor.Normalised) read. A missing precursor has no row, as DIA-NN writes it. A few
#' precursors fail the run or global q-value, a few peptides are not proteotypic and four unchanged proteins fail the
#' global protein-group q-value, so the confidence filter removes rows.
#' Precursor.Quantity is the raw quantity and keeps the recipe's per-sample loading (kitProteinMatrix sample_scale).
#' Precursor.Normalised removes that loading as DIA-NN's cross-run normalisation would (KIT-101 F-19): each run is
#' centred on its median log2 Precursor.Quantity, to the median of the run medians, with a small residual run error
#' (run_scale, SD 0.1 log2).
#' @noRd
generateProtDiannFixture <- function(lane, lane_dir) {
    inputs <- protFixtureInputs(lane)
    universe <- inputs$universe
    lm <- inputs$matrix
    design <- lm$design
    report <- withr::with_seed(as.integer(lane$fixture_recipe$seed) + 1L, {
        n_obs <- lengths(universe$observed)
        peptide <- data.frame(protein = rep(seq_along(n_obs), n_obs), sequence = unlist(universe$observed),
            stringsAsFactors = FALSE)
        peptide$proteotypic <- as.integer(stats::runif(nrow(peptide)) >= 0.03)
        peptide$offset <- stats::rnorm(nrow(peptide), 0, 0.8)
        third <- which(stats::runif(nrow(peptide)) < 0.25)
        precursor <- peptide[c(seq_len(nrow(peptide)), third), , drop = FALSE]
        precursor$charge <- c(rep(2L, nrow(peptide)), rep(3L, length(third)))
        precursor <- precursor[order(precursor$protein, precursor$sequence, precursor$charge), , drop = FALSE]
        precursor$offset <- precursor$offset + stats::rnorm(nrow(precursor), 0, 0.3)
        precursor$global_q <- signif(stats::runif(nrow(precursor), 1e-4, 5e-3), 3)
        precursor$global_q[stats::runif(nrow(precursor)) < 0.01] <- 0.02
        protein_global_pg_q <- signif(stats::runif(length(n_obs), 1e-4, 2e-3), 3)
        protein_global_pg_q[sample(which(!lm$differential), 4L)] <- 0.03
        modified <- gsub("C", "C(UniMod:4)", precursor$sequence, fixed = TRUE)
        run_scale <- 2^stats::rnorm(nrow(design), 0, 0.1)
        rows <- lapply(seq_len(nrow(design)), \(s) {
            log2_quantity <- log2(lm$values[precursor$protein, s]) + precursor$offset +
                stats::rnorm(nrow(precursor), 0, 0.15)
            q_value <- signif(stats::runif(nrow(precursor), 1e-4, 8e-3), 3)
            q_value[stats::runif(nrow(precursor)) < 0.02] <- 0.03
            pg_q_value <- signif(stats::runif(length(n_obs), 1e-4, 5e-3), 3)
            absent <- lm$missing[precursor$protein, s] |
                stats::runif(nrow(precursor)) < 0.03 + 0.1 * stats::plogis(17 - log2_quantity)
            quantity <- round(2^log2_quantity, 1)
            data.frame(
                Run = design$sample[[s]],
                Protein.Group = universe$ids[precursor$protein],
                Protein.Ids = universe$ids[precursor$protein],
                Protein.Names = universe$entries[precursor$protein],
                Genes = universe$genes[precursor$protein],
                Modified.Sequence = modified,
                Stripped.Sequence = precursor$sequence,
                Precursor.Id = paste0(modified, precursor$charge),
                Precursor.Charge = precursor$charge,
                Q.Value = q_value,
                Global.Q.Value = precursor$global_q,
                PG.Q.Value = pg_q_value[precursor$protein],
                Global.PG.Q.Value = protein_global_pg_q[precursor$protein],
                Proteotypic = precursor$proteotypic,
                Precursor.Quantity = quantity,
                Precursor.Normalised = NA_real_,
                stringsAsFactors = FALSE
            )[!absent, , drop = FALSE]
        })
        report <- do.call(rbind, rows)
        # F-19: the per-run loading removed (median centring of each run), the residual run error kept.
        run_median <- tapply(log2(report$Precursor.Quantity), report$Run, stats::median)
        shift <- stats::median(run_median) - run_median[report$Run]
        report$Precursor.Normalised <- round(report$Precursor.Quantity * 2^shift * run_scale[match(report$Run, design$sample)], 1)
        report
    })
    .write_tsv(report, file.path(lane_dir, lane$seed_file))
    cat("  wrote", lane$seed_file, "(", nrow(report), "rows )\n")
    writeProtFixtureCommon(lane, lane_dir, universe, design)
}

#' Proteome Discoverer TMT export: one row per protein (Annotated.Sequence, Master.Protein.Accessions) and one TMT10
#' plex, `Abundance: F1: <channel>, <sample>`; a missing abundance is left empty, as Proteome Discoverer writes it.
#' @noRd
generateProtPdTmtFixture <- function(lane, lane_dir) {
    inputs <- protFixtureInputs(lane)
    universe <- inputs$universe
    lm <- inputs$matrix
    values <- lm$values
    values[lm$missing] <- NA
    channels <- lm$recipe$channels[seq_len(nrow(lm$design))]
    table <- data.frame(
        Annotated.Sequence = vapply(universe$observed, `[[`, character(1), 1L),
        Master.Protein.Accessions = universe$ids,
        stringsAsFactors = FALSE
    )
    abundances <- as.data.frame(values)
    names(abundances) <- sprintf("Abundance: F1: %s, %s", channels, lm$design$sample)
    table <- cbind(table, abundances)
    .write_tsv(table, file.path(lane_dir, lane$seed_file))
    cat("  wrote", lane$seed_file, "(", nrow(table), "rows )\n")
    writeProtFixtureCommon(lane, lane_dir, universe, lm$design)
}

#' MaxQuant proteinGroups.txt in the dotted form the reader handles, a missing LFQ intensity written as 0, and
#' contaminant (CON__) and decoy (REV__) rows interleaved, one of each leading the table, as the KIT-021 corpus writer
#' lays them out (the reader drops them).
#' @noRd
generateProtMaxQuantFixture <- function(lane, lane_dir) {
    inputs <- protFixtureInputs(lane)
    universe <- inputs$universe
    lm <- inputs$matrix
    values <- lm$values
    values[lm$missing] <- 0
    table <- data.frame(Protein.IDs = universe$ids, Majority.protein.IDs = universe$ids,
        Protein.names = universe$names, Gene.names = universe$genes, stringsAsFactors = FALSE)
    intensities <- as.data.frame(values)
    names(intensities) <- paste0("LFQ.intensity.", lm$design$sample)
    table <- cbind(table, intensities)
    table$Potential.contaminant <- ""
    table$Reverse <- ""
    flagged_rows <- lm$recipe$flagged_rows
    flagged <- table[seq_len(flagged_rows), , drop = FALSE]
    flagged$Protein.IDs <- flagged$Majority.protein.IDs <- paste0("CON__", flagged$Protein.IDs)
    flagged$Potential.contaminant <- "+"
    decoys <- table[seq_len(flagged_rows) + flagged_rows, , drop = FALSE]
    decoys$Protein.IDs <- decoys$Majority.protein.IDs <- paste0("REV__", decoys$Protein.IDs)
    decoys$Reverse <- "+"
    table <- rbind(table, flagged, decoys)
    table <- withr::with_seed(lm$recipe$seed + 1L, table[sample.int(nrow(table)), , drop = FALSE])
    lead <- c(which(table$Potential.contaminant == "+")[[1L]], which(table$Reverse == "+")[[1L]])
    table <- table[c(lead, setdiff(seq_len(nrow(table)), lead)), , drop = FALSE]
    .write_tsv(table, file.path(lane_dir, lane$seed_file))
    cat("  wrote", lane$seed_file, "(", nrow(table), "rows )\n")
    writeProtFixtureCommon(lane, lane_dir, universe, lm$design)
}

#' FragPipe combined_protein.tsv: protein annotation columns and `<sample> MaxLFQ Intensity`, a missing intensity
#' written as 0.
#' @noRd
generateProtFragPipeFixture <- function(lane, lane_dir) {
    inputs <- protFixtureInputs(lane)
    universe <- inputs$universe
    lm <- inputs$matrix
    values <- lm$values
    values[lm$missing] <- 0
    coverage <- round(100 * vapply(universe$observed, \(p) sum(nchar(p)), numeric(1)) / nchar(universe$sequences), 1)
    table <- data.frame(
        Protein = sprintf("sp|%s|%s", universe$ids, universe$entries),
        `Protein ID` = universe$ids,
        `Entry Name` = universe$entries,
        Gene = universe$genes,
        Length = nchar(universe$sequences),
        Organism = "Homo sapiens",
        `Protein Description` = universe$names,
        `Protein Existence` = "1:Experimental evidence at protein level",
        Coverage = coverage,
        `Protein Probability` = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    intensities <- as.data.frame(values)
    names(intensities) <- paste0(lm$design$sample, " MaxLFQ Intensity")
    table <- cbind(table, intensities)
    .write_tsv(table, file.path(lane_dir, lane$seed_file))
    cat("  wrote", lane$seed_file, "(", nrow(table), "rows )\n")
    writeProtFixtureCommon(lane, lane_dir, universe, lm$design)
}

# ---------------------------------------------------------------------------
# Retained six-sample proteomics fixtures (prot_lfq/, prot_tmt/): their original writers, unchanged, so the files
# the native-session tests and the KIT-021 lifecycle members read are reproduced byte for byte.
# ---------------------------------------------------------------------------

.write_proteomics_fasta <- function(path) {
    lines <- c(
        ">sp|P00001|E2E_PROT1 Protein Alpha OS=Homo sapiens OX=9606 GN=GENE1 PE=1 SV=1",
        "MPEPAKLVRAAA",
        ">sp|P00002|E2E_PROT2 Protein Beta OS=Homo sapiens OX=9606 GN=GENE2 PE=1 SV=1",
        "MPEPBKLVRAAA",
        ">sp|P00003|E2E_PROT3 Protein Gamma OS=Homo sapiens OX=9606 GN=GENE3 PE=1 SV=1",
        "MPEPCKLVRAAA",
        ">sp|P00004|E2E_PROT4 Protein Delta OS=Homo sapiens OX=9606 GN=GENE4 PE=1 SV=1",
        "MPEPDKLVRAAA",
        ">sp|P00005|E2E_PROT5 Protein Epsilon OS=Homo sapiens OX=9606 GN=GENE5 PE=1 SV=1",
        "MPEPEKLVRAAA",
        ">sp|P00006|E2E_PROT6 Protein Zeta OS=Homo sapiens OX=9606 GN=GENE6 PE=1 SV=1",
        "MPEPFKLVRAAA"
    )
    writeLines(lines, path, useBytes = TRUE)
    invisible(path)
}

.make_design <- function(lane) {
    groups <- unlist(lane$groups)
    n_per_group <- lane$sample_count %/% length(groups)
    samples <- unlist(lapply(groups, \(g) paste0(g, "_", seq_len(n_per_group))))
    group_labels <- rep(groups, each = n_per_group)
    data.frame(sample = samples, group = group_labels, stringsAsFactors = FALSE)
}

# Two "DA-significant" features: feature 1 is up in KO (~2x), feature 2 is down in KO (~0.5x).
.base_intensities <- function(n_features, n_wt, n_ko, base_mean = 1000.0, noise_sd = 50.0) {
    set.seed(20260428L)  # reset per call for reproducibility
    wt_mat <- matrix(
        rnorm(n_features * n_wt, mean = base_mean, sd = noise_sd)
        , nrow = n_features
        , ncol = n_wt
    )
    ko_mat <- matrix(
        rnorm(n_features * n_ko, mean = base_mean, sd = noise_sd)
        , nrow = n_features
        , ncol = n_ko
    )
    # Inject 2x fold change for feature 1 (KO up)
    if (n_features >= 1L) ko_mat[1L, ] <- ko_mat[1L, ] * 2.0
    # Inject 0.5x fold change for feature 2 (KO down)
    if (n_features >= 2L) ko_mat[2L, ] <- ko_mat[2L, ] * 0.5
    list(wt = round(wt_mat, 1L), ko = round(ko_mat, 1L))
}

.generate_retained_prot_tmt <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    proteins <- data.frame(
        accession = c("P00001", "P00002", "P00003", "P00004", "P00005")
        , gene    = c("GENE1",  "GENE2",  "GENE3",  "GENE4",  "GENE5")
        , desc    = c("Protein Alpha", "Protein Beta", "Protein Gamma", "Protein Delta", "Protein Epsilon")
        , stringsAsFactors = FALSE
    )
    n_proteins <- nrow(proteins)
    intensities <- .base_intensities(n_proteins, n_wt, n_ko)

    # TMT channel names (TMT10: 126, 127N, 127C, 128N, 128C, 129N for 6 samples)
    tmt_channels <- c("126", "127N", "127C", "128N", "128C", "129N")
    all_samples <- c(wt_samples, ko_samples)
    abundance_cols <- sprintf(
        "Abundance: F%d: %s, %s",
        seq_along(all_samples),
        tmt_channels[seq_along(all_samples)],
        all_samples
    )

    # One peptide per protein (simplified — TMT uses protein-level abundances here)
    peptides <- paste0("PEP", seq_len(n_proteins), "K")

    df <- data.frame(
        Annotated.Sequence       = peptides
        , Master.Protein.Accessions = proteins$accession
        , stringsAsFactors = FALSE
    )

    for (i in seq_along(all_samples)) {
        is_wt <- i <= n_wt
        grp_mat <- if (is_wt) intensities$wt else intensities$ko
        col_idx <- if (is_wt) i else i - n_wt
        df[[abundance_cols[i]]] <- grp_mat[, col_idx]
    }

    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(df, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(df), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")

    if (!is.null(lane$fasta_file)) {
        .write_proteomics_fasta(file.path(lane_dir, lane$fasta_file))
        cat("  wrote", lane$fasta_file, "\n")
    }
}

.generate_retained_prot_lfq <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    proteins <- data.frame(
        Protein.IDs = c("P00001", "P00002", "P00003", "P00004", "P00005")
        , Gene.names = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
        , stringsAsFactors = FALSE
    )
    n_proteins <- nrow(proteins)
    intensities <- .base_intensities(n_proteins, n_wt, n_ko)

    all_samples <- c(wt_samples, ko_samples)

    if (identical(lane$import_tool, "fragpipe")) {
        entry_names <- c("ALPHA_HUMAN", "BETA_HUMAN", "GAMMA_HUMAN", "DELTA_HUMAN", "EPSIL_HUMAN")
        descriptions <- c("Protein Alpha", "Protein Beta", "Protein Gamma", "Protein Delta", "Protein Epsilon")
        df <- data.frame(
            Protein = paste0("sp|", proteins$Protein.IDs, "|", entry_names)
            , `Protein ID` = proteins$Protein.IDs
            , `Entry Name` = entry_names
            , Gene = proteins$Gene.names
            , Length = c(350L, 280L, 420L, 195L, 510L)
            , Organism = "Homo sapiens"
            , `Protein Description` = descriptions
            , `Protein Existence` = 1L
            , Coverage = c(42.5, 38.2, 35.1, 51.3, 29.8)
            , check.names = FALSE
            , stringsAsFactors = FALSE
        )
        lfq_cols <- paste0(all_samples, " MaxLFQ Intensity")
    } else {
        df <- proteins
        lfq_cols <- paste0("LFQ.intensity.", all_samples)
    }

    for (i in seq_along(all_samples)) {
        is_wt <- i <= n_wt
        grp_mat <- if (is_wt) intensities$wt else intensities$ko
        col_idx <- if (is_wt) i else i - n_wt
        df[[lfq_cols[i]]] <- grp_mat[, col_idx]
    }

    if (!identical(lane$import_tool, "fragpipe")) {
        df$Potential.contaminant <- ""
        df$Reverse <- ""
    }

    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(df, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(df), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")

    if (!is.null(lane$fasta_file)) {
        .write_proteomics_fasta(file.path(lane_dir, lane$fasta_file))
        cat("  wrote", lane$fasta_file, "\n")
    }
}

# The retained fixtures' original lane definitions (manifest entries before revision 2.2).
.RETAINED_WRITERS <- list(
    list(fixture_dir = "prot_tmt", import_tool = "pd_tmt", seed_file = "pd_tmt_peptides.tsv", writer = .generate_retained_prot_tmt),
    list(fixture_dir = "prot_lfq", import_tool = "maxquant", seed_file = "proteinGroups.txt", writer = .generate_retained_prot_lfq),
    list(fixture_dir = "prot_lfq", import_tool = "fragpipe", seed_file = "seed_combined_protein.tsv", writer = .generate_retained_prot_lfq)
)

generateRetainedFixtures <- function() {
    for (spec in .RETAINED_WRITERS) {
        cat("[retained ", spec$fixture_dir, "/", spec$seed_file, "]\n", sep = "")
        lane <- list(import_tool = spec$import_tool, seed_file = spec$seed_file, fasta_file = "proteins.fasta",
            sample_count = 6L, groups = list("WT", "KO"))
        spec$writer(lane, file.path(FIXTURE_ROOT, spec$fixture_dir))
    }
}

# ---------------------------------------------------------------------------
# Transcriptomics
# ---------------------------------------------------------------------------

# KIT-101: the transcriptomics screen lane, the KIT-040 trx_screen corpus at a
# smaller shape (3 biological samples per group, 2 technical replicates, 300
# genes), written by the corpus tool itself so the recipe is not copied.
.generate_trx_screen <- function(lane, lane_dir) {
    corpus_tool <- file.path(FIXTURE_ROOT, "..", "..", "..", "tools", "profiling", "kit_omics_corpus.R")
    if (!file.exists(corpus_tool)) stop("the KIT-040 corpus tool is not in this source tree: ", corpus_tool)
    corpus_env <- new.env(parent = globalenv())
    sys.source(corpus_tool, corpus_env)
    shape <- lane$fixture_recipe
    recipe <- corpus_env$kitOmicsRecipe("trx_screen", as.integer(shape$biological_per_group), as.integer(shape$features))
    corpus_env$writeKitTrxCorpus(lane_dir, recipe)
    cat("  wrote", paste(sort(list.files(lane_dir)), collapse = ", "), "\n")
}

# ---------------------------------------------------------------------------
# Dispatch table: import_tool -> generator function. Lanes whose fixture_recipe$generator is "none" are hand-made
# and are not written.
# ---------------------------------------------------------------------------

.GENERATORS <- list(
    diann       = generateProtDiannFixture
    , pd_tmt    = generateProtPdTmtFixture
    , maxquant  = generateProtMaxQuantFixture
    , fragpipe  = generateProtFragPipeFixture
    , trx_counts = .generate_trx_screen
)

# ---------------------------------------------------------------------------
# Check: every file the manifest records has its recorded sha256.
# ---------------------------------------------------------------------------

#' Compare each recorded file's sha256 with the file on disk.
#' @param recorded Named list: path relative to the fixture root -> sha256.
#' @return Character vector of mismatches (empty when all match).
#' @noRd
checkFixtureSha256 <- function(recorded) {
    problems <- character()
    for (path in names(recorded)) {
        full <- file.path(FIXTURE_ROOT, path)
        actual <- if (file.exists(full)) digest::digest(file = full, algo = "sha256") else "<missing>"
        if (!identical(actual, recorded[[path]])) problems <- c(problems, sprintf("%s: recorded %s, now %s", path, recorded[[path]], actual))
    }
    problems
}

recordedFixtureSha256 <- function(manifest) {
    recorded <- list()
    for (lane in manifest$lanes) {
        for (name in names(lane$files_sha256)) recorded[[file.path(lane$fixture_dir, name)]] <- lane$files_sha256[[name]]
    }
    retained <- manifest$retained_fixtures
    c(recorded, retained$regenerated_sha256, retained$hand_made_sha256)
}

# ---------------------------------------------------------------------------
# Main generation loop
# ---------------------------------------------------------------------------

generate_all_fixtures <- function() {
    for (pkg in c("jsonlite", "withr", "digest")) {
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
    }

    cat("\n=== E2E Fixture Generation ===\n")

    manifest <- .read_manifest()
    cat("Loaded manifest with", length(manifest$lanes), "lanes.\n\n")

    errors <- character()

    for (lane in manifest$lanes) {
        lane_id  <- lane$lane_id
        lane_dir <- file.path(FIXTURE_ROOT, lane$fixture_dir)

        cat("[", lane_id, "]\n", sep = "")

        if (identical(lane$fixture_recipe$generator, "none")) {
            cat("  hand-made: not written (checked below)\n")
            next
        }

        if (!dir.exists(lane_dir)) {
            dir.create(lane_dir, recursive = TRUE)
            cat("  created directory:", lane$fixture_dir, "\n")
        }

        generator <- .GENERATORS[[lane$import_tool]]
        if (is.null(generator)) {
            msg <- paste0("No generator for import_tool '", lane$import_tool, "' in lane '", lane_id, "'")
            cat("  ERROR:", msg, "\n")
            errors <- c(errors, msg)
            next
        }

        tryCatch(
            generator(lane, lane_dir)
            , error = function(e) {
                msg <- paste0("Generator failed for lane '", lane_id, "': ", e$message)
                cat("  ERROR:", msg, "\n")
                errors <<- c(errors, msg)
            }
        )
    }

    tryCatch(generateRetainedFixtures(), error = function(e) {
        msg <- paste0("Retained fixture writer failed: ", e$message)
        cat("  ERROR:", msg, "\n")
        errors <<- c(errors, msg)
    })

    cat("\n=== Checking every recorded fixture sha256 ===\n")
    mismatches <- checkFixtureSha256(recordedFixtureSha256(manifest))
    for (m in mismatches) cat("  MISMATCH", m, "\n")
    errors <- c(errors, mismatches)
    if (!length(mismatches)) cat("  all", length(recordedFixtureSha256(manifest)), "recorded files match\n")

    if (length(errors) > 0L) {
        stop("Fixture generation completed with ", length(errors), " error(s).")
    }

    cat("\n=== Fixture generation complete ===\n")
    invisible(NULL)
}

# Run immediately when sourced or executed as a script
generate_all_fixtures()
