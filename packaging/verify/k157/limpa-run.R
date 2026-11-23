# KIT-157 (branch run-k157-rc only): DIA_LIMPA quantification on the fastlimpadpc engine, saved for the fixed comparison
# rule (packaging/verify/compare-results.R) and timed.
# Usage: <bundle>/bin/Rscript --vanilla limpa-run.R <out.rds> <library first on the path, or ->
# With a library given, fastlimpadpc loads from it (the previous pin); with "-", from the bundle's library (this pin).
# Cases: the KIT-129 run-path fixture of test-kit129-run-path.R item 5 (seed 129, 120 proteins, 5 per group); three data
# sets of simulator S2's core size (limpa::simProteinDataSet, 4,000 precursors of 1,000 proteins, 6 per group, two
# groups, one unwanted factor); one larger data set for the time (20,000 precursors, 10 per group). Simulation and
# quantification as test-kit129-run-path.R (kit129LimpaSimulate, kit129LimpaQuantify, copied), i.e. the package's own
# proteinMissingValueImputationLimpa(engine = "fastlimpadpc") as the route calls it. Each case's quantified object is
# saved as a plain list of its slots. Prints K157RUN lines only.
arguments <- commandArgs(trailingOnly = TRUE)
out <- arguments[[1L]]
if (length(arguments) > 1L && !identical(arguments[[2L]], "-")) .libPaths(c(arguments[[2L]], .libPaths()))
suppressPackageStartupMessages(library(MultiScholaR))

simulate <- function(seed, proteins = 200L, per_group = 6L, groups = 2L, factors = 1L, strength = 3, units = 0L) {
    set.seed(seed)
    sim <- limpa::simProteinDataSet(n.peptides = proteins * 4L, n.groups = groups, samples.per.group = per_group,
        peptides.per.protein = 4L, mu.range = c(2, 10), prop.de = 0, dpc.slope = 0.7, prop.missing = 0.25)
    y <- sim$other$E.complete
    samples <- ncol(y)
    ids <- unique(sim$genes$Protein)
    protein_of <- match(sim$genes$Protein, ids)
    group_index <- sim$targets$Group
    contrast <- c(-0.5, 0.5, 0)[group_index]
    n_effects <- round(0.2 * proteins)
    delta <- numeric(proteins)
    delta[sample.int(proteins, n_effects)] <- rep(c(1, -1), c(ceiling(0.5 * n_effects), n_effects - ceiling(0.5 * n_effects)))
    term <- outer(delta, contrast)
    if (factors > 0L) {
        scores <- matrix(stats::rnorm(samples * factors), samples)
        loadings <- sqrt(strength * 0.09) * matrix(stats::rnorm(proteins * factors), proteins)
        term <- term + loadings %*% t(scores)
    }
    within_group <- stats::ave(group_index, group_index, FUN = seq_along)
    replicates <- within_group
    if (units > 0L) {
        replicates <- ceiling(within_group / (per_group / units))
        unit <- (group_index - 1L) * units + replicates
        term <- term + matrix(stats::rnorm(proteins * max(unit), 0, 0.3), proteins)[, unit]
    }
    y <- y + term[protein_of, , drop = FALSE]
    probability <- stats::plogis(sim$dpc[[1L]] + y * sim$dpc[[2L]], lower.tail = FALSE)
    y[stats::rbinom(length(y), prob = probability, size = 1L) == 1L] <- NA
    order <- sample.int(nrow(y))
    runs <- sprintf("run_%03d", seq_len(samples))
    design <- data.frame(Run = runs, group = c("A", "B", "C")[group_index], replicates = replicates, stringsAsFactors = FALSE)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    list(y = y[order, , drop = FALSE], protein = sim$genes$Protein[order], precursor = rownames(y)[order], design = design)
}

quantify <- function(sim, engine) {
    y <- sim$y
    dimnames(y) <- list(paste(sim$protein, sim$precursor, sep = "%"), sim$design$Run)
    object <- methods::new("PeptideQuantitativeData")
    object@peptide_data <- data.frame(Protein.Ids = sim$protein, Stripped.Sequence = sim$precursor, stringsAsFactors = FALSE)
    object@peptide_matrix <- y
    object@protein_id_column <- "Protein.Ids"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@is_logged_data <- TRUE
    object@design_matrix <- sim$design
    object@sample_id <- "Run"
    object@group_id <- "group"
    object@technical_replicate_id <- "tech_rep_group"
    object@args <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))
    counts <- tibble::tibble(Protein.Ids = unique(sim$protein), peptide_count = 4L, peptidoform_count = 4L)
    quantified <- NULL
    utils::capture.output(quantified <- suppressWarnings(suppressMessages(proteinMissingValueImputationLimpa(object,
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = engine, peptide_summary = counts))))
    quantified
}

plain <- function(x) {
    if (isS4(x) && !is.data.frame(x)) {
        slots <- methods::slotNames(x)
        return(stats::setNames(lapply(slots, function(name) plain(methods::slot(x, name))), slots))
    }
    if (is.environment(x)) return("<environment>")
    if (is.function(x)) return("<function>")
    if (is.list(x) && !is.data.frame(x)) return(lapply(x, plain))
    x
}

cases <- list(
    item5 = list(seed = 129L, proteins = 120L, per_group = 5L),
    s2_1 = list(seed = 1L, proteins = 1000L, per_group = 6L),
    s2_2 = list(seed = 2L, proteins = 1000L, per_group = 6L),
    s2_3 = list(seed = 3L, proteins = 1000L, per_group = 6L),
    large = list(seed = 157L, proteins = 5000L, per_group = 10L))
results <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    sim <- simulate(case$seed, proteins = case$proteins, per_group = case$per_group)
    seconds <- system.time(quantified <- quantify(sim, "fastlimpadpc"))[["elapsed"]]
    results[[name]] <- plain(quantified)
    cat(sprintf("K157RUN case %s: %d precursors x %d runs; %s; %.2f s\n", name, nrow(sim$y), ncol(sim$y),
        paste(class(quantified), collapse = "/"), seconds))
}
dll <- getLoadedDLLs()[["fastlimpadpc"]][["path"]]
cat(sprintf("K157RUN fastlimpadpc from %s; RemoteSha %s; library md5 %s; RAYON_NUM_THREADS %s\n",
    dirname(find.package("fastlimpadpc")), utils::packageDescription("fastlimpadpc")$RemoteSha %||% "none",
    unname(tools::md5sum(dll)), Sys.getenv("RAYON_NUM_THREADS")))
saveRDS(results, out)
