# Reproducible synthetic input generation for PERF-107; no production activation.

#' Define fixed synthetic training and held-out cohort profiles.
#' @param profile Either calibration or heldout.
#' @param small Use the bounded eight-sample correctness fixture.
#' @return Frozen profile with distinct sample and feature identities.
#' @noRd
perf107SyntheticRecipe <- function(profile, small = FALSE) {
    stopifnot(profile %in% c("calibration", "heldout"), is.logical(small), length(small) == 1L)
    held <- identical(profile, "heldout")
    seed <- if (held) 20261009L else 20261008L
    samples <- if (small) 8L else 300L
    universe <- if (small) 400L else 60000L
    counts <- withr::with_seed(seed, pmin(universe,
        pmax(500L, as.integer(round(exp(stats::rnorm(samples, log(16000), if (held) 0.85 else 0.65)))))))
    counts[1:8] <- if (small) c(40L, 99L, 100L, 101L, 170L, 220L, 300L, 400L) else
        c(500L, 2999L, 3000L, 3001L, 5000L, 12000L, 40000L, 60000L)
    list(profile = profile, synthetic = TRUE, seed = seed, samples = samples,
        universe = universe, counts = counts, cutoff = if (small) 100L else 3000L,
        proteins = if (small) 20L else 2500L,
        feature_offset = if (held) 1000000L else 0L,
        annotation_tail_bytes = if (held) 2048L else 512L)
}

#' Construct the finite catalogue using the existing fixture sequence encoder.
#' @param recipe Frozen synthetic recipe.
#' @return Feature catalogue; not a full assay.
#' @noRd
perf107SyntheticCatalogue <- function(recipe) {
    i <- seq_len(recipe$universe)
    # Non-uniform groups: high-occupancy groups coexist with many smaller groups.
    protein <- pmin(recipe$proteins, 1L + as.integer(recipe$proteins * (i / (recipe$universe + 1))^1.7))
    id <- sprintf("SYN_%s_P%05d", recipe$profile, protein)
    sequence <- workflowCheckpointSequences(i + recipe$feature_offset)
    width <- ifelse(protein %% 97L == 0L, recipe$annotation_tail_bytes,
        ifelse(protein %% 13L == 0L, 128L, 24L))
    data.frame(feature = i, Protein.Group = id,
        Protein.Ids = ifelse(protein %% 11L == 0L, paste0(id, ";", id, "_ISO2"), id),
        Stripped.Sequence = sequence,
        Protein.Description = paste0("synthetic_", strrep("a", width)),
        stringsAsFactors = FALSE)
}

#' Generate one complete sample without retaining other samples.
#' @param recipe,catalogue Frozen inputs.
#' @param sample Integer sample index.
#' @return DIA-NN-compatible precursor rows, bounded by twice the catalogue size.
#' @noRd
perf107SyntheticSample <- function(recipe, catalogue, sample) {
    stopifnot(sample %in% seq_len(recipe$samples))
    withr::with_seed(recipe$seed + sample, {
        selected <- sample.int(recipe$universe, recipe$counts[[sample]], replace = FALSE)
        feature <- catalogue[selected, , drop = FALSE]
        run <- sprintf("%s_run_%03d", recipe$profile, sample)
        # Two charge states for a minority of peptides; identity remains unique.
        at <- c(seq_len(nrow(feature)), which(selected %% 17L == 0L))
        form <- c(rep(2L, nrow(feature)), rep(3L, sum(selected %% 17L == 0L)))
        feature <- feature[at, , drop = FALSE]
        log_quantity <- 20 + (feature$feature %% 37L) / 8 +
            (sample %% 3L == 0L) * (feature$feature %% 5L == 0L) +
            stats::rnorm(length(at), sd = 0.45)
        quantity <- as.double(2^log_quantity)
        zero <- stats::runif(length(at)) < 0.02
        quantity[zero] <- 0
        missing <- stats::runif(length(at)) < 0.04 + 0.10 * stats::plogis(21 - log_quantity)
        missing <- missing | feature$feature %% 997L == 0L |
            (sample %% 3L == 0L & feature$feature %% 29L == 0L)
        quantity[missing] <- NA_real_
        # Normalized values deliberately differ, including at raw-missing sites.
        normal <- as.double(2^log_quantity * (1 + sample / 1000))
        modified <- ifelse(form == 3L, paste0("M[Oxidation]", substring(feature$Stripped.Sequence, 2)),
            feature$Stripped.Sequence)
        result <- data.frame(Run = run, Protein.Group = feature$Protein.Group,
            Protein.Ids = feature$Protein.Ids, Stripped.Sequence = feature$Stripped.Sequence,
            Modified.Sequence = modified, Precursor.Charge = form,
            Precursor.Id = paste0(modified, "_", form),
            Precursor.Quantity = quantity, Precursor.Normalised = normal,
            Q.Value = 0.001, Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001,
            PG.Q.Value = 0.001, Proteotypic = 1L,
            Protein.Description = feature$Protein.Description, stringsAsFactors = FALSE)
        rownames(result) <- NULL
        result
    })
}

#' Write and independently scan a cohort in existing bounded reader batches.
#' @param directory Fresh owned output directory.
#' @param recipe Frozen generation recipe.
#' @param fail_after_sample Optional small failure-instrumentation seam.
#' @return Manifest with input hashes and independently observed sample counts.
#' @noRd
writePerf107SyntheticCohort <- function(directory, recipe, fail_after_sample = NULL) {
    if (dir.exists(directory) || !dir.create(directory, recursive = TRUE)) stop("output already exists")
    started <- proc.time()[["elapsed"]]
    completed <- FALSE
    on.exit(jsonlite::write_json(list(status = if (completed) "complete" else "failed",
        elapsed_seconds = proc.time()[["elapsed"]] - started, synthetic = TRUE),
        file.path(directory, "generation-status.json"), auto_unbox = TRUE, pretty = TRUE), add = TRUE)
    catalogue <- perf107SyntheticCatalogue(recipe)
    report <- file.path(directory, "report.tsv")
    expected <- vector("list", recipe$samples)
    for (sample in seq_len(recipe$samples)) {
        block <- perf107SyntheticSample(recipe, catalogue, sample)
        expected[[sample]] <- data.frame(Run = block$Run[[1L]], rows = nrow(block),
            distinct_peptides = recipe$counts[[sample]], missing = sum(is.na(block$Precursor.Quantity)),
            zero = sum(block$Precursor.Quantity == 0, na.rm = TRUE))
        # Keep write batches bounded even when a complete sample is larger.
        for (start in seq.int(1L, nrow(block), by = 10000L)) {
            at <- seq.int(start, min(nrow(block), start + 9999L))
            vroom::vroom_write(block[at, , drop = FALSE], report, delim = "\t",
                append = file.exists(report), col_names = !file.exists(report),
                progress = FALSE, num_threads = 8L)
        }
        if (!is.null(fail_after_sample) && sample == fail_after_sample) stop("injected generation failure")
    }
    expected <- do.call(rbind, expected)
    design <- data.frame(sample = expected$Run,
        group = rep(c("control", "treated", "contrast"), length.out = recipe$samples),
        replicates = ceiling(seq_len(recipe$samples) / 9),
        batch = paste0("batch_", 1L + (seq_len(recipe$samples) - 1L) %/% 30L))
    vroom::vroom_write(design, file.path(directory, "design.tsv"), delim = "\t", progress = FALSE)
    vroom::vroom_write(expected, file.path(directory, "expected-sample-counts.tsv"), delim = "\t", progress = FALSE)
    fasta <- file(file.path(directory, "proteins.fasta"), "w")
    tryCatch({
        for (protein in unique(catalogue$Protein.Group)) {
            subset <- catalogue[catalogue$Protein.Group == protein, , drop = FALSE]
            sequence <- paste0(subset$Stripped.Sequence, collapse = "")
            for (id in unique(unlist(strsplit(subset$Protein.Ids, ";", fixed = TRUE)))) {
                writeLines(c(paste0(">", id, " synthetic PERF-107 protein"), sequence), fasta)
            }
        }
    }, finally = close(fasta))
    # Independent file scan: no regeneration used as the observed result.
    observed <- stats::setNames(vector("list", recipe$samples), expected$Run)
    seen <- new.env(parent = emptyenv())
    current <- NULL
    callback <- readr::SideEffectChunkCallback$new(function(data, position) {
        stopifnot(nrow(data) <= 10000L, all(is.finite(data$Precursor.Quantity) | is.na(data$Precursor.Quantity)),
            all(data$Precursor.Quantity >= 0 | is.na(data$Precursor.Quantity)))
        for (run in unique(data$Run)) {
            part <- data[data$Run == run, , drop = FALSE]
            if (!identical(run, current)) {
                if (!is.null(current)) observed[[current]]$distinct_peptides <<- length(ls(seen))
                stopifnot(is.null(observed[[run]]))
                current <<- run
                seen <<- new.env(parent = emptyenv())
                observed[[run]] <<- list(rows = 0L, missing = 0L, zero = 0L)
            }
            keys <- paste(part$Protein.Group, part$Stripped.Sequence, sep = "|")
            for (key in unique(keys)) assign(key, TRUE, envir = seen)
            observed[[run]]$rows <<- observed[[run]]$rows + nrow(part)
            observed[[run]]$missing <<- observed[[run]]$missing + sum(is.na(part$Precursor.Quantity))
            observed[[run]]$zero <<- observed[[run]]$zero + sum(part$Precursor.Quantity == 0, na.rm = TRUE)
        }
    })
    readr::read_tsv_chunked(report, callback, chunk_size = 10000L,
        col_types = readr::cols(.default = readr::col_character(),
            Precursor.Quantity = readr::col_double()), progress = FALSE)
    observed[[current]]$distinct_peptides <- length(ls(seen))
    for (i in seq_len(nrow(expected))) for (field in c("rows", "missing", "zero", "distinct_peptides")) {
        if (!identical(as.integer(observed[[expected$Run[[i]]]][[field]]),
            as.integer(expected[[field]][[i]]))) stop(paste("readback differs",
                expected$Run[[i]], field, observed[[expected$Run[[i]]]][[field]], expected[[field]][[i]]))
    }
    files <- file.path(directory, c("report.tsv", "design.tsv", "proteins.fasta", "expected-sample-counts.tsv"))
    result <- list(synthetic = TRUE, profile = recipe$profile, recipe = recipe,
        scope = "Synthetic report inputs; not measured real-study applicability or post-QC fingerprints",
        rows = sum(expected$rows), samples = nrow(expected), verified_bounded_scan = TRUE,
        expected_filter_distinct_population = sum(expected$distinct_peptides[expected$distinct_peptides >= recipe$cutoff]),
        expected_retained_samples = sum(expected$distinct_peptides >= recipe$cutoff),
        expectation_boundary = "distinct protein/stripped-peptide counts in report; rebind after any upstream filtering",
        missing = sum(expected$missing), zeros = sum(expected$zero),
        files = stats::setNames(lapply(files, function(path) list(path = path,
            bytes = file.size(path), sha256 = digest::digest(file = path, algo = "sha256"))), basename(files)))
    jsonlite::write_json(result, file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
    completed <- TRUE
    result
}

#' Generate the two frozen 300-sample synthetic qualification cohorts.
#' @param run_dir,source Existing guarded runner output and source snapshot.
#' @return Both verified cohort manifests.
#' @noRd
runPerf107SyntheticCorpus <- function(run_dir, source) {
    sys.source(file.path(source, "tests/testthat/helper-workflow-checkpoint-fixtures.R"),
        environment(runPerf107SyntheticCorpus))
    profiles <- c("calibration", "heldout")
    result <- stats::setNames(lapply(profiles, function(profile) {
        writePerf107SyntheticCohort(file.path(run_dir, profile), perf107SyntheticRecipe(profile))
    }), profiles)
    saveRDS(result, file.path(run_dir, "synthetic-corpus.rds"))
    result
}
