# KIT-183 simulation (ticket section 6.2; declaration docs/qa/workflow-execution/kit/kit-183/declaration.md).
# Run on the MOFA2 evidence lane (branch run-mofa-evidence-k183, .github/workflows/k183-evidence.yml, which runs this
# file from the package commit under test; the lane's .github/mofa-smoke/k183.R is a byte copy), with MultiScholaR
# installed at that commit, MOFA2 and the bundle's Python.
#
# Usage: Rscript simulate-k183.R <out dir> <task> [<cell n> <first replicate> <last replicate>]
#   task  replicates  the replicates of one cell, each: generate, write the inputs, import, train (shipped defaults,
#                     training_seeds = 1), analyse F_main, F_weak, F_none and F_ind with every test and comparator;
#                     one line per replicate in status.tsv, aggregated rows in the other tables
#         seeds       F3: one replicate of the cell trained with 1 and with 5 seeds; factors and weights compared
#         memory      8.7: TestValues (pages of 10,000 rows) and camera on 20,000 features x 500 samples, 2,000 sets
#         fail        a replicate with a declared invalid argument (n = 7): recorded as did not complete
# Every table is appended to, one file per table, in <out dir>.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 2L)
k183_out <- arguments[[1L]]
k183_task <- arguments[[2L]]
k183_arguments <- arguments
dir.create(k183_out, recursive = TRUE, showWarnings = FALSE)
k183_out <- normalizePath(k183_out)

local(envir = new.env(parent = asNamespace("MultiScholaR")), {
# Prefixed names outside: the namespace and its imports would otherwise be searched first.
out <- k183_out
task <- k183_task
arguments <- k183_arguments

# --- Declared constants (declaration.md) ---

K <- list(
    seed_base = c("12" = 20261183, "20" = 20261183 + 100000, "6" = 20261183 + 200000),
    ind_seed = 20261187, perm_seed = 20261185, set_seed = 20261188,
    views = list(
        transcriptomics = list(omic = "transcriptomics", view = "transcriptomics", features = 4000L, mu = c(5, 2.5),
            a = 0.25, b = 1.0, m = 0, s = 2, z2 = TRUE, censor = NA_real_, blocks_per = 3L, sizes = c(10L, 30L, 100L),
            block_sizes = c(10L, 30L, 100L), c2_sizes = c(10L, 30L), c4_sizes = c(10L, 30L, 100L),
            p_sizes = c(10L, 30L, 100L), id_type = "gene", database = "string"),
        proteomics = list(omic = "proteomics", view = "proteomics", features = 1800L, mu = c(20, 2),
            a = 0.15, b = 0.6, m = 14, s = 2, z2 = TRUE, censor = 17, blocks_per = 2L, sizes = c(10L, 30L, 100L),
            block_sizes = c(10L, 30L, 100L), c2_sizes = c(10L, 15L), c4_sizes = c(10L, 30L), p_sizes = c(10L, 30L),
            id_type = "uniprot", database = "string"),
        metabolomics = list(omic = "metabolomics", view = "metabolomics_LCMS", features = 250L, mu = c(16, 2),
            a = 0.15, b = 0.6, m = 14, s = 2, z2 = FALSE, censor = 13, blocks_per = 1L, sizes = c(5L, 10L, 20L),
            block_sizes = c(5L, 10L), c2_sizes = c(5L, 10L), c4_sizes = c(5L, 10L, 20L), p_sizes = 5L, id_type = "kegg",
            database = "kegg")),
    # Declaration, dated section 2026-11-24 (after L1 run 45219821547, before S): proteomics C2 sizes 10, 15 (was 10,
    # 30), C4 and P sizes 10, 30 (was 10, 30, 100); metabolomics P size 5 (was 5, 10, 20): those formed fewer than 90 %
    # of their sets in L1.
    rhos = c(0.1, 0.3, 0.6), z1_share = 0.10, z2_share = 0.05, z1_loading = c(0.5, 1.5), z2_loading = c(0.3, 1.0),
    multi_identifiers = 31L, metab_identifiers = 204L, adduct_correlation = 0.9, censor_slope = 1.5,
    c0_per_size = 10L, c2_per_size = 5L, c3_per_size = 5L, c3_sizes = c(5L, 10L), c4_per_size = 5L, p_per_size = 5L,
    min_out = 5L, permutations = 1000L, mofa_min_size = 10L, alpha = 0.05,
    tests = c("camera", "gsea", "camera_pr", "camera_floor", "pcgse_cor", "pcgse_perm", "mofa_default"))

append <- function(table, name) {
    if (is.null(table) || !nrow(table)) return(invisible(NULL))
    path <- file.path(out, paste0(name, ".tsv"))
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(path),
        append = file.exists(path))
}

# --- The generator (6.2) ---

# The generator is data: `generator` holds the constants and the data sets (K below); a later ticket passes another.
simulateReplicate <- function(n, r, generator = K) {
    K <- generator
    if (n %% 2L != 0L || n < 4L) stop("the number of biological samples must be even and at least 4: ", n)
    seed <- K$seed_base[[as.character(n)]] + r
    set.seed(seed)
    samples <- sprintf("S%02d", seq_len(n))
    group <- rep(c("A", "B"), each = n / 2L)
    z1 <- ifelse(group == "A", -1, 1) + stats::rnorm(n, 0, 0.25)
    z2 <- stats::rnorm(n)
    z2 <- stats::residuals(stats::lm(z2 ~ z1))
    z2 <- z2 / stats::sd(z2)
    views <- lapply(K$views, function(spec) simulateView(spec, n, z1, z2, K))
    # n = 12: the transcriptomics data set lacks one biological sample (the last), kept as missing (keep_missing).
    present <- lapply(views, function(v) samples)
    if (n == 12L) present$transcriptomics <- samples[-n]
    list(n = n, r = r, seed = seed, samples = samples, group = group, z1 = z1, z2 = z2, views = views, present = present)
}

simulateView <- function(spec, n, z1, z2, generator = K) {
    K <- generator
    g <- spec$features
    mu <- stats::rnorm(g, spec$mu[[1L]], spec$mu[[2L]])
    sigma <- spec$a + spec$b * exp(-(mu - spec$m) / spec$s)
    # Identifiers: one per feature, except in metabolomics, where 31 of 204 identifiers carry two or three features
    # (alternately 2 and 3): the first 16 of them independent features, the last 15 adducts of one compound.
    if (spec$omic == "metabolomics") {
        multi_sizes <- rep_len(c(2L, 3L), K$multi_identifiers)
        unit_sizes <- c(multi_sizes, rep(1L, g - sum(multi_sizes)))
        unit <- rep(seq_along(unit_sizes), unit_sizes)
        kind <- rep(c(rep("independent", 16L), rep("adduct", 15L), rep("single", length(unit_sizes) - 31L)), unit_sizes)
        identifier <- sprintf("C%05d", 10000L + unit)
        feature_id <- sprintf("m%04d", seq_len(g))
    } else {
        unit <- seq_len(g)
        kind <- rep("single", g)
        identifier <- if (spec$omic == "proteomics") sprintf("P%05d", seq_len(g)) else sprintf("G%05d", seq_len(g))
        feature_id <- identifier
    }
    # Loadings, in units of each feature's sigma: Z1 on 10% of the loading units (features; an adduct group is one
    # unit sharing one loading), half positive; Z2 (transcriptomics, proteomics) on 5% of the other units.
    units <- unique(unit)
    l1_unit <- stats::setNames(numeric(length(units)), units)
    l2_unit <- l1_unit
    active1 <- sample(units, round(K$z1_share * length(units)))
    signs1 <- rep_len(c(1, -1), length(active1))
    l1_unit[as.character(active1)] <- signs1 * stats::runif(length(active1), K$z1_loading[[1L]], K$z1_loading[[2L]])
    if (spec$z2) {
        rest <- setdiff(units, active1)
        active2 <- sample(rest, round(K$z2_share * length(units)))
        l2_unit[as.character(active2)] <- sample(c(-1, 1), length(active2), TRUE) *
            stats::runif(length(active2), K$z2_loading[[1L]], K$z2_loading[[2L]])
    }
    l1 <- unname(l1_unit[as.character(unit)])
    l2 <- unname(l2_unit[as.character(unit)])
    # Independent multi-feature identifiers: each feature its own loading draw (10% active).
    independent <- which(kind == "independent")
    if (length(independent)) {
        on <- stats::runif(length(independent)) < K$z1_share
        l1[independent] <- ifelse(on, sample(c(-1, 1), length(independent), TRUE) *
            stats::runif(length(independent), K$z1_loading[[1L]], K$z1_loading[[2L]]), 0)
    }
    e <- matrix(stats::rnorm(g * n), g, n)
    for (u in unique(unit[kind == "adduct"])) {
        rows <- which(unit == u)
        shared <- stats::rnorm(n)
        e[rows, ] <- sqrt(K$adduct_correlation) * matrix(shared, length(rows), n, byrow = TRUE) +
            sqrt(1 - K$adduct_correlation) * e[rows, , drop = FALSE]
    }
    values <- function(e) mu + sigma * (outer(l1, z1) + outer(l2, z2) + e)
    y <- values(e)
    # Left-censoring on the value (proteomics, metabolomics).
    mask <- matrix(FALSE, g, n)
    if (!is.na(spec$censor)) mask <- matrix(stats::runif(g * n), g, n) < 1 / (1 + exp(K$censor_slope * (y - spec$censor)))
    complete <- rowSums(mask) == 0L
    # Null blocks among complete features with zero loading on Z1 and Z2 (single-feature identifiers), disjoint, in
    # random order; each member's noise becomes sqrt(rho) u_b + sqrt(1 - rho) e (variance unchanged). The mask was
    # drawn on the values before the block latent (declaration: marginal distribution unchanged).
    zero <- l1 == 0 & l2 == 0 & kind == "single"
    pool <- sample(which(complete & zero))
    blocks <- list()
    taken <- 0L
    for (rho in K$rhos) for (size in spec$block_sizes) for (k in seq_len(spec$blocks_per)) {
        if (taken + size > length(pool)) {
            blocks[[length(blocks) + 1L]] <- list(rho = rho, size = size, members = integer(), formed = FALSE)
            next
        }
        members <- pool[taken + seq_len(size)]
        taken <- taken + size
        u <- stats::rnorm(n)
        e[members, ] <- sqrt(rho) * matrix(u, size, n, byrow = TRUE) + sqrt(1 - rho) * e[members, , drop = FALSE]
        blocks[[length(blocks) + 1L]] <- list(rho = rho, size = size, members = members, formed = TRUE)
    }
    y <- values(e)
    y[mask] <- NA_real_
    in_block <- seq_len(g) %in% unlist(lapply(blocks, `[[`, "members"))
    free <- which(complete & zero & !in_block)
    c2_pool <- free[mu[free] <= stats::quantile(mu[free], 0.2, names = FALSE)]
    list(spec = spec, y = y, mu = mu, sigma = sigma, l1 = l1, l2 = l2, kind = kind, unit = unit, identifier = identifier,
        feature_id = feature_id, complete = complete, zero = zero, blocks = blocks, in_block = in_block, free = free,
        c2_pool = c2_pool)
}

# --- Inputs, import and training through the package (the shipped path) ---

writeReplicateInputs <- function(sim, root) {
    integration_dir <- file.path(root, "integration")
    for (name in names(sim$views)) {
        v <- sim$views[[name]]
        spec <- v$spec
        keep <- sim$samples %in% sim$present[[name]]
        samples <- data.frame(Run = sim$samples[keep], group = sim$group[keep], tech_rep_group = sim$samples[keep],
            stringsAsFactors = FALSE)
        matrix <- data.frame(feature_id = v$feature_id, v$y[, keep, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
        names(matrix) <- c("feature_id", sim$samples[keep])
        features <- data.frame(feature_id = v$feature_id, display_name = v$feature_id, stringsAsFactors = FALSE)
        if (spec$omic == "metabolomics") features$kegg <- v$identifier
        manifest <- list(omic_type = spec$omic, workflow = toupper(spec$omic), matrix_role = "normalised",
            scale = if (spec$omic == "transcriptomics") "log2_cpm" else "log2", correction = NULL,
            source_state = list(name = "kit183-simulation", synthetic = TRUE), sample_id_column = "Run",
            group_column = "group", technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
        dir <- file.path(integration_dir, spec$omic)
        dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)
        writeIntegrationInput(dir, manifest, samples, stats::setNames(list(matrix), spec$view),
            stats::setNames(list(features), spec$view))
    }
    integration_dir
}

# The training step is data: `training` holds the model parameters the request sets (KIT-183: training_seeds = 1).
trainReplicate <- function(sim, root, training = list(training_seeds = 1L)) {
    integration_dir <- writeReplicateInputs(sim, root)
    manager <- newIntegrationStateOwner(file.path(root, "store"))
    scan <- integrationScanProject(integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    warnings <- character()
    run <- withCallingHandlers(runIntegrationModel(manager, training),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    list(manager = manager, run = run, warnings = warnings)
}

# --- The comparators (D7, VAL-53 (5); reported only) ---

fisherZ <- function(r) atanh(pmax(pmin(r, 1 - 1e-12), -1 + 1e-12))

# PCGSE's t of a set's mean difference against the other features (`.pcgse_ttest`, MOFA2 1.22.0): mean in the set
# minus mean outside, over the pooled standard deviation, times sqrt(vif/m1 + 1/m2). `stat` is features x columns
# (one column per permutation); the result is sets x columns. Computed from the sums and sums of squares in and out of
# each set (the same quantities as the per-set loop of .pcgse_ttest).
pcgseT <- function(stat, index, vif = NULL) {
    stat <- as.matrix(stat)
    g <- nrow(stat)
    member <- matrix(0, length(index), g)
    for (i in seq_along(index)) member[i, index[[i]]] <- 1
    m1 <- rowSums(member)
    m2 <- g - m1
    sum_in <- member %*% stat
    sq_in <- member %*% stat^2
    sum_all <- matrix(colSums(stat), length(index), ncol(stat), byrow = TRUE)
    sq_all <- matrix(colSums(stat^2), length(index), ncol(stat), byrow = TRUE)
    mean_a <- sum_in / m1
    mean_b <- (sum_all - sum_in) / m2
    ss_a <- sq_in - m1 * mean_a^2
    ss_b <- (sq_all - sq_in) - m2 * mean_b^2
    pooled <- sqrt((ss_a + ss_b) / (m1 + m2 - 2))
    inflation <- if (is.null(vif)) rep(1, length(index)) else vif
    (mean_a - mean_b) / (pooled * sqrt(inflation / m1 + 1 / m2))
}

comparators <- function(values, sets, f, perm_seed, weights_view) {
    y <- values$y
    universe <- rownames(y)
    index <- lapply(sets, function(s) match(s, universe))
    r <- as.numeric(stats::cor(t(y), f))
    z <- fisherZ(r)
    n <- ncol(y)
    # (i) correlation-adjusted: VIF from the marginal correlation matrix of the set's features, t on n - 2 df.
    vif <- vapply(index, function(i) {
        m <- length(i)
        cm <- stats::cor(t(y[i, , drop = FALSE]))
        1 + (m - 1) * (sum(cm) - m) / (m * (m - 1))
    }, numeric(1))
    t_cor <- pcgseT(z, index, vif)[, 1L]
    p_cor <- 2 * stats::pt(-abs(t_cor), df = n - 2)
    # (ii) permutation of the factor's scores across S_d (the model not refitted), the uncorrected t's |t| against its
    # permutation distribution: p = (1 + #{|t*| >= |t|}) / (1 + permutations).
    t_obs <- pcgseT(z, index)[, 1L]
    ys <- t(scale(t(y)))
    perms <- withr::with_seed(perm_seed, vapply(seq_len(K$permutations), function(i) sample(n), integer(n)))
    fs <- matrix(f[perms], n, K$permutations)
    fs <- scale(fs)
    rp <- (ys %*% fs) / (n - 1)
    zp <- fisherZ(rp)
    t_perm <- pcgseT(zp, index)
    p_perm <- (1 + rowSums(abs(t_perm) >= abs(t_obs))) / (1 + K$permutations)
    out <- data.frame(set = names(sets), pcgse_cor = p_cor, pcgse_perm = p_perm, stringsAsFactors = FALSE)
    out$pcgse_cor_t <- t_cor
    out$pcgse_perm_t <- t_obs
    # (iii) MOFA2's default: the uncorrected parametric t on the absolute weights of the model's features of the data
    # set (its universe), sets with at least 10 of them (min.size = 10).
    out$mofa_default <- NA_real_
    out$mofa_default_t <- NA_real_
    if (!is.null(weights_view)) {
        stat <- abs(weights_view$weight)
        model_index <- lapply(sets, function(s) which(weights_view$identifier %in% s))
        ok <- lengths(model_index) >= K$mofa_min_size
        if (any(ok)) {
            t_mofa <- pcgseT(stat, model_index[ok])[, 1L]
            m1 <- lengths(model_index[ok])
            out$mofa_default_t[ok] <- t_mofa
            out$mofa_default[ok] <- 2 * stats::pt(-abs(t_mofa), df = length(stat) - 2)
        }
    }
    out
}

# --- Sets of the classes (6.2), drawn per replicate, data set and factor ---
# The classes are data: CLASSES is a list of functions, each drawing its sets through `ctx$add` and counting what it
# formed through `ctx$count`, in the list's order (the order fixes the random stream). A later ticket (KIT-184) adds
# classes by passing another list; KIT-183's list is the one below.

drawPerSize <- function(ctx, class, sizes, per_size, pool, draw) {
    for (size in sizes) {
        made <- 0L
        for (k in seq_len(per_size)) if (length(pool) >= size) {
            ctx$add(sprintf("%s_%d_%d", class, size, k), draw(size), class, size)
            made <- made + 1L
        }
        ctx$count(class, size, per_size, made)
    }
}

CLASSES <- list(
    # C0: uniform from U_cam.
    C0 = function(ctx) drawPerSize(ctx, "C0", ctx$v$spec$sizes, K$c0_per_size, ctx$universe,
        function(size) sample(ctx$universe, size)),
    # C1: each formed null block; C1-out its members the model did not select; C4b those it selected.
    C1 = function(ctx) {
        v <- ctx$v
        for (b in seq_along(v$blocks)) {
            block <- v$blocks[[b]]
            ctx$count("C1", block$size, 1L, as.integer(block$formed))
            if (!block$formed) next
            ctx$add(sprintf("C1_%g_%d_%d", block$rho, block$size, b), ctx$id[block$members], "C1", block$size, block$rho)
            used <- v$feature_id[block$members] %in% ctx$selected
            if (v$spec$omic != "metabolomics" && sum(!used) >= K$min_out) {
                ctx$add(sprintf("C1out_%g_%d_%d", block$rho, block$size, b), ctx$id[block$members[!used]], "C1out",
                    block$size, block$rho)
            }
            if (sum(used) >= K$min_out) {
                ctx$add(sprintf("C4b_%g_%d_%d", block$rho, block$size, b), ctx$id[block$members[used]], "C4b", block$size,
                    block$rho)
            }
        }
    },
    # C2: the lowest-abundance fifth of the complete zero-loading features in no block.
    C2 = function(ctx) drawPerSize(ctx, "C2", ctx$v$spec$c2_sizes, K$c2_per_size, ctx$v$c2_pool,
        function(size) ctx$id[sample(ctx$v$c2_pool, size)]),
    # C3 (metabolomics): identifiers carrying two or three zero-loading features, present in U_cam.
    C3 = function(ctx) {
        v <- ctx$v
        if (v$spec$omic != "metabolomics") return(invisible(NULL))
        multi <- unique(v$unit[v$kind != "single"])
        ok <- vapply(multi, function(u) {
            rows <- which(v$unit == u)
            all(v$l1[rows] == 0 & v$l2[rows] == 0) && ctx$id[rows[[1L]]] %in% ctx$universe
        }, logical(1))
        pool <- unique(ctx$id[v$unit %in% multi[ok]])
        drawPerSize(ctx, "C3", K$c3_sizes, K$c3_per_size, pool, function(size) sample(pool, size))
    },
    # C4: zero-loading complete features in no block that the model selected.
    C4 = function(ctx) {
        pool <- ctx$v$free[ctx$v$feature_id[ctx$v$free] %in% ctx$selected]
        drawPerSize(ctx, "C4", ctx$v$spec$c4_sizes, K$c4_per_size, pool, function(size) ctx$id[sample(pool, size)])
    },
    # P-up, P-down: half planted (positive or negative Z1 loading, complete, single-feature identifiers), half from
    # complete zero-loading features in no block. Power on F_main only; drawn on every factor for one code path.
    P = function(ctx) {
        v <- ctx$v
        positive <- which(v$complete & v$l1 > 0 & v$kind == "single" & ctx$in_u)
        negative <- which(v$complete & v$l1 < 0 & v$kind == "single" & ctx$in_u)
        for (direction in c("up", "down")) {
            planted <- if (direction == "up") positive else negative
            for (size in v$spec$p_sizes) {
                made <- 0L
                half <- size %/% 2L
                for (k in seq_len(K$p_per_size)) if (length(planted) >= half && length(v$free) >= size - half) {
                    ctx$add(sprintf("P%s_%d_%d", direction, size, k), ctx$id[c(sample(planted, half),
                        sample(v$free, size - half))], paste0("P", direction), size, direction = direction)
                    made <- made + 1L
                }
                ctx$count(paste0("P", direction), size, K$p_per_size, made)
            }
        }
    })

drawSets <- function(v, values, selected, classes = CLASSES) {
    sets <- list()
    meta <- list()
    formed <- list()
    ctx <- new.env()
    ctx$v <- v
    ctx$universe <- rownames(values$y)
    ctx$id <- v$identifier
    ctx$in_u <- v$identifier %in% ctx$universe
    ctx$selected <- selected
    ctx$add <- function(name, members, class, size, rho = NA_real_, direction = NA_character_) {
        sets[[name]] <<- unique(members)
        meta[[name]] <<- data.frame(set = name, class = class, size = length(unique(members)), nominal = size, rho = rho,
            planted = direction, stringsAsFactors = FALSE)
    }
    ctx$count <- function(class, size, wanted, made) formed[[length(formed) + 1L]] <<- data.frame(class = class,
        size = size, wanted = wanted, formed = made, stringsAsFactors = FALSE)
    for (draw in classes) draw(ctx)
    list(sets = sets, meta = do.call(rbind, meta), formed = do.call(rbind, formed))
}

# --- The arms (tests and comparators) ---
# The arms are data: ARMS is a list of functions of one context (the data set's values, ranked list, sets and
# settings), each returning rows `set`, `test`, `pvalue`, `direction`. A later ticket (KIT-184) adds arms by passing
# another list; KIT-183's list is the one below.

armRows <- function(terms, test, direction = terms$Direction) {
    if (is.null(terms) || !nrow(terms)) return(NULL)
    data.frame(set = terms$set, test = test, pvalue = terms$pvalue, direction = direction, stringsAsFactors = FALSE)
}

entryPoint <- function(ctx, test, settings = ctx$settings) {
    terms <- integrationEnrichmentPathwayTest(ctx$ranked, ctx$term2gene, ctx$term2name, 20261121L, settings, test,
        if (test != "gsea") ctx$values)$terms
    terms$set <- terms$ID
    terms
}

# The whole path for metabolomics on a learned factor: runIntegrationEnrichment with a synthetic KEGG table.
wholePath <- function(ctx, test) {
    sets <- ctx$sets
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    mapping <- list(database = "KEGG", release = "KEGG simulation", sha256 = strrep("0", 64L), url = "simulation",
        downloaded_at = "2026-11-23T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
    run <- runIntegrationEnrichment(ctx$manager, list(factor = ctx$factor, view = ctx$view, database = "kegg",
        id_column = "kegg", test = test), mapping = mapping)
    terms <- as.data.frame(run$run$value$terms)
    terms$set <- ctx$term2name$term[match(terms$Description, ctx$term2name$name)]
    terms
}

ARMS <- list(
    camera = function(ctx) armRows(if (ctx$full_path) wholePath(ctx, "camera") else entryPoint(ctx, "camera"), "camera"),
    gsea = function(ctx) {
        terms <- if (ctx$full_path) wholePath(ctx, "gsea") else entryPoint(ctx, "gsea")
        armRows(terms, "gsea", ifelse(terms$NES > 0, "Up", "Down"))
    },
    camera_pr = function(ctx) armRows(entryPoint(ctx, "camera_pr"), "camera_pr"),
    # The sensitivity arm: allow.neg.cor = FALSE.
    camera_floor = function(ctx) armRows(entryPoint(ctx, "camera", c(ctx$settings, list(allow_neg_cor = FALSE))),
        "camera_floor"),
    # The C3 arm: camera with the collapse of the pre-ranked test (largest absolute correlation among complete features).
    camera_rcollapse = function(ctx) {
        c3 <- names(ctx$sets)[startsWith(names(ctx$sets), "C3_")]
        if (!length(c3)) return(NULL)
        v <- ctx$v
        rk <- ctx$ranking$ranking
        cand <- which(v$complete & v$identifier %in% rownames(ctx$values$y))
        corr <- abs(rk$correlation[match(v$feature_id[cand], rk$feature_id)])
        ord <- cand[order(v$identifier[cand], -corr, cand)]
        keep <- sort(ord[!duplicated(v$identifier[ord])])
        y_r <- v$y[keep, match(ctx$values$samples, ctx$sim$samples), drop = FALSE]
        rownames(y_r) <- v$identifier[keep]
        terms <- integrationEnrichmentCamera(list(y = y_r, design = ctx$values$design, counts = list()),
            ctx$term2gene[ctx$term2gene$term %in% c3, , drop = FALSE], ctx$term2name[ctx$term2name$term %in% c3, ,
                drop = FALSE], ctx$settings, "camera")$terms
        terms$set <- terms$ID
        armRows(terms, "camera_rcollapse")
    },
    # The comparators (D7, reported only).
    comparators = function(ctx) {
        weights_view <- if (ctx$kind != "F_ind") {
            rows <- ctx$weights[ctx$weights$view == ctx$view, , drop = FALSE]
            data.frame(identifier = ctx$v$identifier[match(rows$feature_id, ctx$v$feature_id)], weight = rows[[ctx$factor]],
                stringsAsFactors = FALSE)
        }
        comp <- comparators(ctx$values, ctx$sets, ctx$values$f, K$perm_seed + ctx$r, weights_view)
        do.call(rbind, lapply(c("pcgse_cor", "pcgse_perm", "mofa_default"), function(test) {
            stat <- comp[[paste0(test, "_t")]]
            data.frame(set = comp$set, test = test, pvalue = comp[[test]],
                direction = ifelse(is.na(stat), NA_character_, ifelse(stat > 0, "Up", "Down")), stringsAsFactors = FALSE)
        }))
    })

# --- One data set and one factor ---

analyse <- function(sim, v, fit, factor, kind, value_aug, object, export, settings, r, classes = CLASSES,
    arms = ARMS) {
    spec <- v$spec
    manager <- fit$manager
    view <- spec$view
    p <- list(database = spec$database, omic = spec$omic, id_type = spec$id_type, id_column = if (spec$omic == "metabolomics") "kegg")
    weights <- fit$weights
    selected <- weights$feature_id[weights$view == view]
    ranking <- integrationEnrichmentRank(object, manager, export$source, value_aug, export$parameters, factor, view,
        "pearson", 5L)
    mapped <- integrationEnrichmentMappedList(ranking, p, object@features[[view]])
    values <- integrationEnrichmentTestValues(object, manager, export$source, value_aug, export$parameters, ranking,
        mapped, view)
    drawn <- drawSets(v, values, selected, classes)
    sets <- drawn$sets
    term2gene <- data.frame(term = rep(names(sets), lengths(sets)), gene = unlist(sets, use.names = FALSE),
        stringsAsFactors = FALSE)
    term2name <- data.frame(term = names(sets), name = names(sets), stringsAsFactors = FALSE)
    ranked <- mapped$table[c("identifier", "value")]
    ctx <- list(sim = sim, v = v, kind = kind, factor = factor, view = view, manager = manager, ranking = ranking,
        values = values, ranked = ranked, sets = sets, term2gene = term2gene, term2name = term2name, settings = settings,
        r = r, weights = weights, full_path = spec$omic == "metabolomics" && kind != "F_ind")
    long <- do.call(rbind, lapply(arms, function(arm) arm(ctx)))
    long <- long[!is.na(long$pvalue), , drop = FALSE]
    meta <- drawn$meta
    long <- merge(long, meta, by = "set")
    long$padj <- stats::ave(long$pvalue, long$test, FUN = function(p) stats::p.adjust(p, "BH"))
    # camera's diagnostic: within-set variance of z against camera's pooled variance.
    fit <- limma::eBayes(limma::lmFit(values$y, values$design), trend = TRUE)
    z <- limma::zscoreT(fit$t[, "factor"], df = fit$df.total, approx = TRUE, method = "hill")
    names(z) <- rownames(values$y)
    g <- length(z)
    diag <- do.call(rbind, lapply(names(sets), function(s) {
        inside <- names(z) %in% sets[[s]]
        m <- sum(inside)
        delta <- mean(z[inside]) - mean(z[!inside])
        pooled <- ((g - 1) * stats::var(z) - delta^2 * m * (g - m) / g) / (g - 2)
        data.frame(set = s, ratio = stats::var(z[inside]) / pooled, stringsAsFactors = FALSE)
    }))
    diag <- merge(diag, meta[c("set", "class", "rho")], by = "set")
    list(long = long, formed = drawn$formed, diag = diag, counts = values$counts)
}

# --- One replicate ---

orientationSign <- function(weights, factor) {
    # The declared rule, coded here independently of the kernel: the weight of largest absolute value over every data
    # set is positive (the first such feature in the model's data set and feature order on a tie).
    w <- weights[[factor]]
    i <- which.max(abs(w))
    if (w[[i]] < 0) -1L else 1L
}

runReplicate <- function(n, r, cell_label, design = DESIGN) {
    started <- Sys.time()
    sim <- simulateReplicate(n, r, design$generator)
    root <- tempfile("k183-")
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    t0 <- Sys.time()
    fit <- trainReplicate(sim, root, design$training)
    on.exit(try(closeWorkflowStateOwner(fit$manager), silent = TRUE), add = TRUE, after = FALSE)
    train_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    manager <- fit$manager
    value <- fit$run$value
    fit$weights <- integrationModelWeights(manager, value)
    scores <- value$factors
    at <- match(scores$biological_sample, sim$samples)
    factors <- setdiff(names(scores), "biological_sample")
    oriented <- vapply(factors, function(k) orientationSign(fit$weights, k) * scores[[k]], numeric(nrow(scores)))
    c1 <- abs(stats::cor(oriented, sim$z1[at]))[, 1L]
    c2 <- abs(stats::cor(oriented, sim$z2[at]))[, 1L]
    f_main <- factors[[which.max(c1)]]
    others <- setdiff(factors, f_main)
    f_weak <- others[[which.max(c2[others])]]
    remaining <- setdiff(others, f_weak)
    f_none <- remaining[[which.min(pmax(c1[remaining], c2[remaining]))]]
    # F_ind: one N(0, 1) value per biological sample, drawn independently of the data.
    value_aug <- value
    value_aug$factors$F_ind <- withr::with_seed(K$ind_seed + r, stats::rnorm(nrow(scores)))
    export <- manager$.context$snapshot$exports$integration_mofa_model
    object <- integrationImportedObject(manager, fields = c("features", "design_matrix", "inputs"))
    settings <- integrationGseaSettings(p_adjust_method = "BH", min_size = 5L, max_size = 500L, p_cutoff = 1, tie_jitter = 0)
    kinds <- c(F_main = f_main, F_weak = f_weak, F_none = f_none, F_ind = "F_ind")
    variance <- value$variance_explained
    factor_rows <- do.call(rbind, lapply(names(kinds), function(kind) {
        k <- kinds[[kind]]
        do.call(rbind, lapply(names(design$generator$views), function(name) {
            view <- design$generator$views[[name]]$view
            rows <- variance[variance$view == view, , drop = FALSE]
            data.frame(kind = kind, factor = k, view = view,
                cor_z1 = if (kind == "F_ind") abs(stats::cor(value_aug$factors$F_ind, sim$z1[at])) else c1[[k]],
                cor_z2 = if (kind == "F_ind") abs(stats::cor(value_aug$factors$F_ind, sim$z2[at])) else c2[[k]],
                r2_percent = if (kind == "F_ind") NA_real_ else rows$r2_percent[rows$factor == k],
                largest_r2_percent = max(rows$r2_percent), stringsAsFactors = FALSE)
        }))
    }))
    # The orientation and Z1's sign for the scorer's expected direction (computed here, not by the kernel).
    raw_cor_main <- stats::cor(scores[[f_main]], sim$z1[at])
    sign_main <- orientationSign(fit$weights, f_main)
    rows <- list()
    formed <- list()
    diag <- list()
    counts <- list()
    analysis_seconds <- 0
    for (name in names(design$generator$views)) {
        v <- sim$views[[name]]
        for (kind in names(kinds)) {
            t1 <- Sys.time()
            a <- analyse(sim, v, fit, kinds[[kind]], kind, value_aug, object, export, settings, r, design$classes,
                design$arms)
            analysis_seconds <- analysis_seconds + as.numeric(difftime(Sys.time(), t1, units = "secs"))
            a$long$view <- v$spec$view
            a$long$kind <- kind
            rows[[length(rows) + 1L]] <- a$long
            a$formed$view <- v$spec$view
            a$formed$kind <- kind
            formed[[length(formed) + 1L]] <- a$formed
            a$diag$view <- v$spec$view
            a$diag$kind <- kind
            diag[[length(diag) + 1L]] <- a$diag
            counts[[length(counts) + 1L]] <- data.frame(view = v$spec$view, kind = kind, a$counts, stringsAsFactors = FALSE)
        }
    }
    long <- do.call(rbind, rows)
    long$sig <- long$pvalue <= K$alpha
    long$bh <- long$padj <= K$alpha
    # Error-rate rows: per data set, factor, test, class (and rho): sets and sets with p <= 0.05.
    key <- c("view", "kind", "test", "class", "rho")
    long$rho_key <- ifelse(is.na(long$rho), "", as.character(long$rho))
    err <- stats::aggregate(cbind(sets = 1L, sig = as.integer(sig)) ~ view + kind + test + class + rho_key, data = long,
        FUN = sum)
    # Planted sets: called at BH <= 0.05, and the direction reported.
    planted <- long[long$class %in% c("Pup", "Pdown"), , drop = FALSE]
    # Column names only inside the formula: the data hold a column `planted` (the planted direction), which would mask
    # the data frame of that name (L1 runs 45148992728 and 45151737782).
    power <- if (nrow(planted)) stats::aggregate(cbind(sets = 1L, called = as.integer(bh),
        called_up = as.integer(bh & direction %in% "Up"),
        called_down = as.integer(bh & direction %in% "Down")) ~ view + kind + test + class + nominal,
        data = planted, FUN = sum)
    # Family-wise: any BH <= 0.05 among the C0 and C2 sets of the call (BH recomputed over them).
    nulls <- long[long$class %in% c("C0", "C2"), , drop = FALSE]
    fw <- do.call(rbind, lapply(split(nulls, list(nulls$view, nulls$kind, nulls$test), drop = TRUE), function(d)
        data.frame(view = d$view[[1L]], kind = d$kind[[1L]], test = d$test[[1L]], any_bh = any(stats::p.adjust(d$pvalue,
            "BH") <= K$alpha), stringsAsFactors = FALSE)))
    diag <- do.call(rbind, diag)
    diag$rho_key <- ifelse(is.na(diag$rho), "", as.character(diag$rho))
    diag <- stats::aggregate(ratio ~ view + kind + class + rho_key, data = diag, FUN = mean)
    # Censoring losses.
    losses <- do.call(rbind, lapply(names(design$generator$views), function(name) {
        v <- sim$views[[name]]
        loaded <- v$l1 != 0 | v$l2 != 0
        data.frame(view = v$spec$view, zero_lost = mean(!v$complete[!loaded]), loaded_lost = mean(!v$complete[loaded]),
            missing_share = mean(is.na(v$y)), stringsAsFactors = FALSE)
    }))
    stamp <- function(d) if (is.null(d) || !nrow(d)) d else cbind(data.frame(cell = cell_label, n = n, r = r), d)
    append(stamp(err), "errors")
    append(stamp(power), "power")
    append(stamp(fw), "familywise")
    append(stamp(diag), "diagnostic")
    append(stamp(do.call(rbind, formed)), "formed")
    append(stamp(do.call(rbind, counts)), "counts")
    append(stamp(losses), "losses")
    append(stamp(factor_rows), "factors")
    append(stamp(data.frame(f_main = f_main, sign_main = sign_main, raw_cor_main_z1 = raw_cor_main,
        warnings = gsub("[\t\r\n]+", " ", paste(unique(fit$warnings), collapse = " | ")), stringsAsFactors = FALSE)), "orientation")
    data.frame(train_seconds = train_seconds, analysis_seconds = analysis_seconds,
        total_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
}

safeReplicate <- function(n, r, cell_label, design = DESIGN) {
    started <- Sys.time()
    # The error's call stack (the innermost calls) is kept with its message, so a failure can be located.
    stack <- character()
    outcome <- tryCatch(withCallingHandlers(runReplicate(n, r, cell_label, design), error = function(e) {
        calls <- vapply(sys.calls(), function(call) paste(deparse(call, nlines = 1L), collapse = ""), character(1))
        stack <<- utils::tail(calls, 12L)
    }), error = function(e) e)
    ok <- !inherits(outcome, "error")
    if (!ok) outcome$message <- paste(conditionMessage(outcome), "| calls:", paste(substr(stack, 1L, 160L), collapse = " <- "))
    append(data.frame(cell = cell_label, n = n, r = r, status = if (ok) "completed" else "did_not_complete",
        message = if (ok) "" else gsub("[\t\r\n]+", " ", conditionMessage(outcome)),
        train_seconds = if (ok) outcome$train_seconds else NA_real_,
        analysis_seconds = if (ok) outcome$analysis_seconds else NA_real_,
        total_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), stringsAsFactors = FALSE), "status")
    invisible(ok)
}

# --- KIT-183's design: the generator, the training step, the classes and the arms, as data ---

DESIGN <- list(generator = K, training = list(training_seeds = 1L), classes = CLASSES, arms = ARMS)

# --- Tasks ---

if (task == "replicates") {
    n <- as.integer(arguments[[3L]])
    first <- as.integer(arguments[[4L]])
    last <- as.integer(arguments[[5L]])
    # No replicate starts after the job's deadline (K183_DEADLINE, epoch seconds); those are recorded as not started.
    deadline <- as.numeric(Sys.getenv("K183_DEADLINE", unset = "Inf"))
    for (r in seq.int(first, last)) {
        if (as.numeric(Sys.time()) > deadline) {
            append(data.frame(cell = paste0("n", n), n = n, r = r, status = "not_started", message = "job deadline",
                train_seconds = NA_real_, analysis_seconds = NA_real_, total_seconds = NA_real_), "status")
            next
        }
        safeReplicate(n, r, paste0("n", n))
    }
} else if (task == "fail") {
    safeReplicate(7L, 1L, "fail")
} else if (task == "seeds") {
    n <- as.integer(arguments[[3L]])
    r <- as.integer(arguments[[4L]])
    sim <- simulateReplicate(n, r)
    one <- trainReplicate(sim, tempfile("k183-s1-"), list(training_seeds = 1L))
    five <- trainReplicate(sim, tempfile("k183-s5-"), list(training_seeds = 5L))
    w1 <- integrationModelWeights(one$manager, one$run$value)
    w5 <- integrationModelWeights(five$manager, five$run$value)
    append(data.frame(n = n, r = r, factors_identical = identical(one$run$value$factors, five$run$value$factors),
        weights_identical = identical(w1, w5), factors_max_abs_difference = max(abs(as.matrix(one$run$value$factors[-1L]) -
            as.matrix(five$run$value$factors[-1L]))), stringsAsFactors = FALSE), "seeds")
    closeWorkflowStateOwner(one$manager)
    closeWorkflowStateOwner(five$manager)
} else if (task == "memory") {
    g <- 20000L
    n <- 500L
    set.seed(20261189)
    f <- stats::rnorm(n)
    samples <- sprintf("S%03d", seq_len(n))
    ids <- sprintf("G%05d", seq_len(g))
    table <- data.frame(feature_id = ids, matrix(stats::rnorm(g * n, 10, 1), g, n), check.names = FALSE)
    names(table) <- c("feature_id", samples)
    object <- methods::new("IntegrationViewsData", views = list(big = table[0L, , drop = FALSE]),
        features = list(big = data.frame(feature_id = ids, display_name = ids)),
        design_matrix = data.frame(omic = "transcriptomics", run = samples, group = "A", biological_sample = samples),
        inputs = list(transcriptomics = list(omic = "transcriptomics", views = "big")))
    value <- list(factors = data.frame(biological_sample = samples, Factor1 = f))
    ranking <- list(factor = "Factor1", min_samples = 5L, direction = list(sign = 1L))
    mapped <- list(mapping = data.frame(feature_id = ids, identifier = ids, status = "mapped"))
    pages <- function(view, columns, kernel) {
        for (start in seq.int(1L, g, by = 10000L)) kernel(table[start:min(g, start + 9999L), columns], start - 1L)
    }
    rm_peak <- function() {
        line <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
        as.numeric(gsub("[^0-9]", "", line)) / 1024
    }
    before <- rm_peak()
    t0 <- Sys.time()
    values <- integrationEnrichmentTestValues(object, NULL, NULL, value, list(), ranking, mapped, "big", pages = pages)
    sets <- withr::with_seed(20261190L, lapply(seq_len(2000L), function(i) sample(ids, sample(5:500, 1L))))
    term2gene <- data.frame(term = rep(sprintf("S%04d", seq_along(sets)), lengths(sets)), gene = unlist(sets))
    term2name <- data.frame(term = unique(term2gene$term), name = unique(term2gene$term))
    settings <- integrationGseaSettings(p_adjust_method = "BH", min_size = 5L, max_size = 500L, p_cutoff = 1, tie_jitter = 0)
    camera <- integrationEnrichmentPathwayTest(data.frame(identifier = ids, value = 0), term2gene, term2name, 1L, settings,
        "camera", values)
    seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    append(data.frame(features = g, samples = n, sets = length(sets), tested = nrow(camera$terms), seconds = seconds,
        peak_rss_mib_before = before, peak_rss_mib_after = rm_peak(), stringsAsFactors = FALSE), "memory")
} else {
    stop("unknown task: ", task)
}
})
