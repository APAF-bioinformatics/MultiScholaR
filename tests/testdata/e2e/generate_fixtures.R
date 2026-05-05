# E2E Fixture Generation Script
#
# DEVELOPER TOOL — not run in CI.
#
# Regenerates all 8 lane fixture seed files and design.tsv from scratch using
# deterministic synthetic data. Run this script whenever fixture formats change
# or new lanes are added.
#
# Usage:
#   source("tests/testdata/e2e/generate_fixtures.R")
#   # or from repo root:
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
    raw <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
    stats::setNames(raw$lanes, vapply(raw$lanes, `[[`, character(1L), "lane_id"))
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

# ---------------------------------------------------------------------------
# Helper: build a design.tsv data.frame for a lane
# ---------------------------------------------------------------------------

.make_design <- function(lane) {
    groups <- unlist(lane$groups)
    n_per_group <- lane$sample_count %/% length(groups)
    samples <- unlist(lapply(groups, \(g) paste0(g, "_", seq_len(n_per_group))))
    group_labels <- rep(groups, each = n_per_group)
    data.frame(sample = samples, group = group_labels, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Synthetic signal helpers
#
# Two "DA-significant" features per lane: feature 1 is up in KO (~2x),
# feature 2 is down in KO (~0.5x). Remaining features are near-null.
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Lane generators — one function per import_tool format
# ---------------------------------------------------------------------------

.generate_prot_diann <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    proteins <- data.frame(
        Protein.Group  = c("P00001", "P00002", "P00003", "P00004", "P00005")
        , Protein.Ids  = c("P00001", "P00002", "P00003", "P00004", "P00005")
        , Protein.Names = c("Protein Alpha", "Protein Beta", "Protein Gamma", "Protein Delta", "Protein Epsilon")
        , Genes        = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
        , stringsAsFactors = FALSE
    )
    n_proteins <- nrow(proteins)

    # Two precursors per protein (suffix K charge 2, R charge 3)
    peptide_suffixes <- list(c("AK", 2L), c("VR", 3L))
    intensities <- .base_intensities(n_proteins, n_wt, n_ko)

    rows <- vector("list", n_proteins * length(peptide_suffixes) * (n_wt + n_ko))
    row_idx <- 0L

    for (s_idx in seq_len(n_wt + n_ko)) {
        is_wt <- s_idx <= n_wt
        sample_name <- if (is_wt) wt_samples[s_idx] else ko_samples[s_idx - n_wt]
        grp_mat <- if (is_wt) intensities$wt else intensities$ko
        col_idx <- if (is_wt) s_idx else s_idx - n_wt

        for (p_idx in seq_len(n_proteins)) {
            prot <- proteins[p_idx, ]
            for (pep in peptide_suffixes) {
                row_idx <- row_idx + 1L
                pep_seq <- paste0("PEP", p_idx, pep[[1L]])
                pep_id  <- paste0("PEP", p_idx, pep[[1L]], pep[[2L]])
                rows[[row_idx]] <- data.frame(
                    Protein.Group      = prot$Protein.Group
                    , Protein.Ids      = prot$Protein.Ids
                    , Protein.Names    = prot$Protein.Names
                    , Genes            = prot$Genes
                    , Precursor.Id     = pep_id
                    , Modified.Sequence = pep_seq
                    , Stripped.Sequence = pep_seq
                    , Precursor.Charge = pep[[2L]]
                    , Precursor.Quantity = grp_mat[p_idx, col_idx]
                    , Q.Value          = 0.001
                    , PG.Q.Value       = 0.001
                    , Run              = sample_name
                    , stringsAsFactors = FALSE
                )
            }
        }
    }

    report <- do.call(rbind, rows)
    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(report, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(report), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")

    if (!is.null(lane$fasta_file)) {
        .write_proteomics_fasta(file.path(lane_dir, lane$fasta_file))
        cat("  wrote", lane$fasta_file, "\n")
    }
}

.generate_prot_tmt <- function(lane, lane_dir) {
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
    abundance_cols <- paste0("Abundance.", tmt_channels[seq_along(all_samples)])

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

.generate_prot_lfq <- function(lane, lane_dir) {
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
    lfq_cols <- paste0("LFQ.intensity.", all_samples)

    df <- proteins
    for (i in seq_along(all_samples)) {
        is_wt <- i <= n_wt
        grp_mat <- if (is_wt) intensities$wt else intensities$ko
        col_idx <- if (is_wt) i else i - n_wt
        df[[lfq_cols[i]]] <- grp_mat[, col_idx]
    }
    df$Potential.contaminant <- ""
    df$Reverse <- ""

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

.generate_metab_lc <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    assay_names <- unlist(lane$assays)  # e.g. c("LCMS_Pos", "LCMS_Neg")
    n_features_per_assay <- 3L

    # mz values per platform
    mz_pos <- c(180.0634, 132.0423, 204.0899)
    mz_neg <- c(175.0552, 191.0556, 269.0454)
    rt_vals <- c(2.45, 3.12, 1.87)
    adducts_pos <- c("[M+H]+", "[M-H]-", "[M+H]+")
    adducts_neg <- c("[M-H]-", "[M+H]+", "[M-H]-")

    rows <- list()
    feat_counter <- 0L

    for (assay in assay_names) {
        is_pos <- grepl("Pos", assay, ignore.case = TRUE)
        intensities <- .base_intensities(n_features_per_assay, n_wt, n_ko, base_mean = 5000.0, noise_sd = 250.0)
        mz_vals   <- if (is_pos) mz_pos else mz_neg
        adducts   <- if (is_pos) adducts_pos else adducts_neg

        for (f_idx in seq_len(n_features_per_assay)) {
            feat_counter <- feat_counter + 1L
            row <- data.frame(
                Feature.Name   = sprintf("Feature_%03d", feat_counter)
                , m.z          = mz_vals[f_idx]
                , Retention.Time = rt_vals[f_idx]
                , Adduct       = adducts[f_idx]
                , stringsAsFactors = FALSE
            )
            for (i in seq_len(n_wt)) {
                row[[wt_samples[i]]] <- intensities$wt[f_idx, i]
            }
            for (i in seq_len(n_ko)) {
                row[[ko_samples[i]]] <- intensities$ko[f_idx, i]
            }
            rows[[feat_counter]] <- row
        }
    }

    df <- do.call(rbind, rows)

    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(df, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(df), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")
}

.generate_metab_gc <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    n_features <- 3L
    intensities <- .base_intensities(n_features, n_wt, n_ko, base_mean = 10000.0, noise_sd = 300.0)
    rt_vals     <- c(5.23, 8.67, 12.34)
    similarity  <- c(850L, 920L, 780L)

    rows <- vector("list", n_features)
    for (f_idx in seq_len(n_features)) {
        row <- data.frame(
            Feature.Name    = sprintf("GC_Feature_%03d", f_idx)
            , Retention.Time = rt_vals[f_idx]
            , Similarity    = similarity[f_idx]
            , stringsAsFactors = FALSE
        )
        for (i in seq_len(n_wt)) {
            row[[wt_samples[i]]] <- intensities$wt[f_idx, i]
        }
        for (i in seq_len(n_ko)) {
            row[[ko_samples[i]]] <- intensities$ko[f_idx, i]
        }
        rows[[f_idx]] <- row
    }

    df <- do.call(rbind, rows)

    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(df, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(df), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")
}

.generate_metab_combined <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    assay_names <- unlist(lane$assays)  # e.g. c("LCMS_Pos", "LCMS_Neg", "GCMS")

    platform_meta <- list(
        LCMS_Pos = list(mz = 180.0634, rt = 2.45, prefix = "Feature")
        , LCMS_Neg = list(mz = 132.0423, rt = 3.12, prefix = "Feature")
        , GCMS     = list(mz = NA_real_, rt = 5.23, prefix = "GC_Feature")
    )

    rows <- list()
    feat_counter <- 0L

    for (assay in assay_names) {
        meta <- platform_meta[[assay]]
        feat_counter <- feat_counter + 1L
        intensities <- .base_intensities(1L, n_wt, n_ko, base_mean = 7000.0, noise_sd = 300.0)

        row <- data.frame(
            Feature.Name    = sprintf("%s_%03d", meta$prefix, feat_counter)
            , Platform      = assay
            , m.z           = meta$mz
            , Retention.Time = meta$rt
            , stringsAsFactors = FALSE
        )
        for (i in seq_len(n_wt)) {
            row[[wt_samples[i]]] <- intensities$wt[1L, i]
        }
        for (i in seq_len(n_ko)) {
            row[[ko_samples[i]]] <- intensities$ko[1L, i]
        }
        rows[[feat_counter]] <- row
    }

    df <- do.call(rbind, rows)
    # Inject DA signal: row 1 up 2x, row 2 down 0.5x
    if (nrow(df) >= 2L) {
        ko_cols <- paste0(ko_samples)
        df[1L, ko_cols] <- df[1L, ko_cols] * 2.0
        df[2L, ko_cols] <- df[2L, ko_cols] * 0.5
    }

    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(df, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(df), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")
}

.generate_lipid_canonical <- function(lane, lane_dir) {
    design <- .make_design(lane)
    wt_samples <- design$sample[design$group == "WT"]
    ko_samples <- design$sample[design$group == "KO"]
    n_wt <- length(wt_samples)
    n_ko <- length(ko_samples)

    lipids <- data.frame(
        LipidMolec  = c("PC(36:2)", "TG(54:3)", "CE(18:1)")
        , Class     = c("PC", "TG", "CE")
        , FattyAcid = c("18:1/18:1", "18:1/18:1/18:1", "18:1")
        , Grade     = c("A", "A", "B")
        , stringsAsFactors = FALSE
    )
    n_lipids <- nrow(lipids)
    intensities <- .base_intensities(n_lipids, n_wt, n_ko, base_mean = 40000.0, noise_sd = 2000.0)

    df <- lipids[, c("LipidMolec", "Class", "FattyAcid")]
    for (i in seq_len(n_wt)) {
        df[[paste0(wt_samples[i], ".MeanArea")]] <- intensities$wt[, i]
    }
    for (i in seq_len(n_ko)) {
        df[[paste0(ko_samples[i], ".MeanArea")]] <- intensities$ko[, i]
    }
    df$Grade <- lipids$Grade

    seed_path <- file.path(lane_dir, lane$seed_file)
    .write_tsv(df, seed_path)
    cat("  wrote", basename(seed_path), "(", nrow(df), "rows )\n")

    design_path <- file.path(lane_dir, "design.tsv")
    .write_tsv(design, design_path)
    cat("  wrote design.tsv\n")
}

# ---------------------------------------------------------------------------
# Dispatch table: import_tool -> generator function
# ---------------------------------------------------------------------------

.GENERATORS <- list(
    diann       = .generate_prot_diann
    , pd_tmt    = .generate_prot_tmt
    , maxquant  = .generate_prot_lfq
    , msdial    = function(lane, lane_dir) {
        # Dispatch on lane_id for msdial variants
        if (lane$lane_id == "metab_gc") {
            .generate_metab_gc(lane, lane_dir)
        } else if (lane$lane_id == "metab_combined") {
            .generate_metab_combined(lane, lane_dir)
        } else {
            .generate_metab_lc(lane, lane_dir)
        }
    }
    , lipidsearch = .generate_lipid_canonical
)

# ---------------------------------------------------------------------------
# Main generation loop
# ---------------------------------------------------------------------------

generate_all_fixtures <- function() {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")
    }

    cat("\n=== E2E Fixture Generation ===\n")
    cat("Seed: 20260428\n\n")

    manifest <- .read_manifest()
    cat("Loaded manifest with", length(manifest), "lanes.\n\n")

    errors <- character()

    for (lane in manifest) {
        lane_id  <- lane$lane_id
        lane_dir <- file.path(FIXTURE_ROOT, lane$fixture_dir)

        cat("[", lane_id, "]\n", sep = "")

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

    if (length(errors) > 0L) {
        cat("\n=== GENERATION ERRORS ===\n")
        for (err in errors) cat("  -", err, "\n")
        stop("Fixture generation completed with ", length(errors), " error(s).")
    }

    cat("\n=== Validating generated fixtures ===\n")

    # Validate using integrity helper (requires testthat context when sourced from tests/)
    valid <- tryCatch({
        # Attempt to call validate_e2e_fixtures if available
        if (exists("validate_e2e_fixtures", mode = "function")) {
            result <- validate_e2e_fixtures()
            if (!result$valid) {
                cat("Validation errors:\n")
                for (err in result$errors) cat("  -", err, "\n")
                FALSE
            } else {
                cat("All fixtures validated successfully.\n")
                TRUE
            }
        } else {
            cat("validate_e2e_fixtures() not in scope — skipping integrity check.\n")
            cat("Run devtools::test() or source helper-e2e-fixture-integrity.R to validate.\n")
            TRUE
        }
    }, error = function(e) {
        cat("Validation skipped (no testthat context):", e$message, "\n")
        TRUE
    })

    if (valid) {
        cat("\n=== Fixture generation complete ===\n")
        cat("All", length(manifest), "lanes generated from set.seed(20260428).\n")
    }

    invisible(NULL)
}

# Run immediately when sourced or executed as a script
generate_all_fixtures()
