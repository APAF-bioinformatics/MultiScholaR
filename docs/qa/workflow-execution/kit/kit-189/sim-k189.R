# KIT-189 (VAL-75, VAL-86): enrichment significance on the proteomics clusterProfiler route, q-value against
# Benjamini-Hochberg. Evidence code under docs/qa: no package code. Runs on a GitHub-hosted runner (run-tests
# `scripts`) with MultiScholaR installed from the bound commit, from the package root.
#
# Usage:
#   Rscript docs/qa/workflow-execution/kit/kit-189/sim-k189.R <part> <cell> <first> <last> <out dir> [<replicate list>]
#   part    a label for the rows (s: counted run S; d: determinism job D; l1: the blind lifecycle; e7k, e7c, e7u: E7)
#   cell    a cell index 1 to 24 (KIT-189 5.2), a comma list of indices, or "all"
#   first, last  replicate range (S: 1 to R; L1: 9001 to 9005; E7: 9011 to 9050)
#   out dir  where the job's tables go; the rows go to <out dir>/rows (part l1: to k189-blind/rows, never uploaded)
#   replicate list  optional file, one replicate number per line: only those replicates of the range run (a
#            continuation job, KIT-189 10 "Attempts")
# Environment: K189_WORKERS (default 4), K189_JOB (the job's name in provenance-<job>.json; default <part>-<cell>).
#
# Each replicate draws its list (every draw before the route is called) under
# set.seed(189000000 + 100000 c + r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection"),
# calls processEnrichments() once per propagation value ("none", "is_a_part_of_regulates") on the same list, and
# writes its 14 rows (2 values x 7 arms) to rows/part-<cell>-<r>.tsv.tmp, renamed to rows/part-<cell>-<r>.tsv once
# complete. At the end the parts of the replicates this run ran are assembled into rows/<part>-<cell>-<first>-<last>.tsv.gz
# and removed. Part l1 also runs the checks E0 to E4 and E6 (KIT-189 7.4) and writes only completion, timing,
# mechanics counts, check results and the SHA-256 of the rows it wrote (blind lifecycle, VAL-11).
# The design, arms, truth and estimands are KIT-189 revision 1.2, sections 5 and 7, and its declaration.

args <- commandArgs(TRUE)
if (!length(args) %in% c(5L, 6L)) {
    stop("usage: Rscript sim-k189.R <part> <cell> <first> <last> <out dir> [<replicate list>]")
}
part <- args[[1L]]
cell_arg <- args[[2L]]
first_rep <- as.integer(args[[3L]])
last_rep <- as.integer(args[[4L]])
out_dir <- args[[5L]]
replicate_list <- if (length(args) == 6L) as.integer(readLines(args[[6L]])) else NULL
stopifnot(grepl("^[a-z0-9]+$", part), !is.na(first_rep), !is.na(last_rep), first_rep <= last_rep,
    is.null(replicate_list) || !anyNA(replicate_list))
lifecycle <- identical(part, "l1")
workers <- as.integer(Sys.getenv("K189_WORKERS", "4"))
job <- Sys.getenv("K189_JOB", paste(part, gsub(",", "_", cell_arg), sep = "-"))
started_at <- Sys.time()
elapsed <- function() proc.time()[["elapsed"]]
t_start <- elapsed()
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rows_dir <- if (lifecycle) file.path("k189-blind", "rows") else file.path(out_dir, "rows")
dir.create(rows_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(MultiScholaR))
h <- new.env(parent = asNamespace("MultiScholaR"))
sys.source(file.path("tests", "testthat", "helper-go-propagation.R"), envir = h)
ns <- asNamespace("MultiScholaR")
build <- get("buildUniprotGoGeneSets", ns)
sortIds <- get("sortEnrichmentGeneIds", ns)
calculateQvalue <- get("calculate_qvalue", asNamespace("enrichit"))
# The kernel's plot files are not written (no estimand depends on them), as in KIT-173's effect script.
restore_widget <- h$goPropReplaceBinding("MultiScholaR", "saveSelfContainedWidget", function(...) invisible(NULL))
restore_ggsave <- h$goPropReplaceBinding("ggplot2", "ggsave", function(...) invisible(NULL))
t_loaded <- elapsed()

seedKit <- function(seed) set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
sortC <- function(x) sort(x, method = "radix")
sha256File <- function(path) digest::digest(file = path, algo = "sha256")
clean <- function(x) substr(gsub("[\"\t\r\n]+", " ", x), 1L, 300L)

# ---- cells (KIT-189 5.2) -------------------------------------------------------------------------------------------
CELLS <- data.frame(
    idx = 1:24,
    bg = c(rep(c("B1", "B2", "B3"), each = 3L), rep(c("B1", "B2", "B3"), each = 5L)),
    scenario = c(rep("null", 9L), rep(c("small", "mid", "large", "mid", "mid"), 3L)),
    L = c(rep(c(25L, 126L, 400L), 3L), rep(c(126L, 126L, 126L, 25L, 400L), 3L)),
    R = c(rep(2000L, 9L), rep(c(400L, 400L, 400L, 200L, 200L), 3L)),
    stringsAsFactors = FALSE)
STRATA <- list(small = c(10L, 30L), mid = c(31L, 100L), large = c(101L, 500L))
VALUES <- c("none", "is_a_part_of_regulates")
ARMS <- c("Q", "BH", "BHf", "Qf", "BY", "Holm", "raw")
CUT <- 0.05
cells_run <- if (identical(cell_arg, "all")) CELLS$idx else as.integer(strsplit(cell_arg, ",", fixed = TRUE)[[1L]])
stopifnot(!anyNA(cells_run), all(cells_run %in% CELLS$idx))

# ---- backgrounds (KIT-189 4) ---------------------------------------------------------------------------------------
b2Entries <- function(raw) {
    entries_all <- sortC(unique(raw$Entry))
    seedKit(189000001)
    sortC(entries_all[sample.int(length(entries_all), 2000L)])
}
loadBackground <- function(code) {
    t0 <- elapsed()
    raw <- h$goPropReadTsv(if (identical(code, "B3")) "uniprot-yeast-559292.tsv.gz" else "uniprot-ecoli-83333.tsv.gz")
    b2 <- NULL
    if (identical(code, "B2")) {
        b2 <- b2Entries(raw)
        committed <- file.path("docs", "qa", "workflow-execution", "kit", "kit-189", "b2-entries.tsv")
        if (file.exists(committed)) {
            stopifnot(identical(utils::read.delim(committed, colClasses = "character")$entry, b2))
        }
        raw <- raw[raw$Entry %in% b2, , drop = FALSE]
    }
    annotations <- h$goPropConvert(raw)
    entries <- sortC(unique(raw$Entry))
    sets <- lapply(stats::setNames(VALUES, VALUES), function(value) build(annotations, value))
    universes <- lapply(sets, function(s) sortC(intersect(entries, unique(s$term2gene[[2L]]))))
    U <- sortC(intersect(universes[[1L]], universes[[2L]]))
    prop <- sets$is_a_part_of_regulates
    bp_terms <- unique(prop$all_terms$TERM[prop$all_terms$ONTOLOGY == "BP"])
    t2g <- prop$term2gene[prop$term2gene[[1L]] %in% bp_terms & prop$term2gene[[2L]] %in% U, , drop = FALSE]
    members <- lapply(split(t2g[[2L]], t2g[[1L]]), function(x) sortC(unique(x)))
    sizes <- lengths(members)
    candidates <- lapply(STRATA, function(r) sortC(names(sizes)[sizes >= r[[1L]] & sizes <= r[[2L]]]))
    list(code = code, raw = raw, b2 = b2, annotations = annotations, entries = entries, sets = sets,
        universes = universes, U = U, members = members, candidates = candidates, seconds = elapsed() - t0)
}
bgs_needed <- unique(CELLS$bg[CELLS$idx %in% cells_run])
if (lifecycle) bgs_needed <- c("B1", "B2", "B3")
BG <- lapply(stats::setNames(bgs_needed, bgs_needed), loadBackground)
if ("B2" %in% names(BG) && lifecycle) {
    utils::write.table(data.frame(entry = BG$B2$b2), file.path(out_dir, "b2-entries.tsv"), sep = "\t", quote = FALSE,
        row.names = FALSE)
}

# ---- one replicate's list (KIT-189 5.2) ----------------------------------------------------------------------------
drawReplicate <- function(cell, r) {
    seed <- 189000000 + 100000 * cell$idx + r
    seedKit(seed)
    bg <- BG[[cell$bg]]
    U <- bg$U
    L <- cell$L
    if (identical(cell$scenario, "null")) {
        return(list(seed = seed, list = sortC(U[sample.int(length(U), L)]), term = NA_character_, T = character(),
            s = 0L))
    }
    candidates <- bg$candidates[[cell$scenario]]
    term <- candidates[[sample.int(length(candidates), 1L)]]
    members <- bg$members[[term]]
    s <- min(L %/% 2L, length(members) %/% 2L)
    inside <- members[sample.int(length(members), s)]
    rest <- U[!U %in% members]
    outside <- rest[sample.int(length(rest), L - s)]
    list(seed = seed, list = sortC(c(inside, outside)), term = term, T = members, s = s)
}

# ---- the route (KIT-189 5.2; KIT-173 effect-measurement.R:74-85) ---------------------------------------------------
runRoute <- function(bg, list_ids, value) {
    table <- data.frame(Protein.IDs = bg$entries, fdr_qvalue = ifelse(bg$entries %in% list_ids, 0.01, 0.5),
        log2FC = ifelse(bg$entries %in% list_ids, 1, 0), stringsAsFactors = FALSE)
    da <- methods::new("da_results_for_enrichment", contrasts = tibble::tibble(contrast = "effect"),
        da_data = list(effect = table), design_matrix = data.frame(sample = "S1", stringsAsFactors = FALSE))
    pathway_dir <- tempfile("k189-pw-")
    dir.create(pathway_dir)
    on.exit(unlink(pathway_dir, recursive = TRUE), add = TRUE)
    route_warnings <- character()
    results <- NULL
    t0 <- elapsed()
    out <- tryCatch({
        utils::capture.output(results <- withCallingHandlers(suppressMessages(processEnrichments(da, "999999",
            q_cutoff = CUT, pathway_dir = pathway_dir, go_annotations = bg$annotations, exclude_iea = FALSE,
            protein_id_column = "Protein.IDs", contrast_names = "effect", go_propagation = value)),
            warning = function(w) {
                route_warnings <<- c(route_warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }))
        list(status = "ok", up = results@enrichment_data$effect$up, message = "")
    }, error = function(e) list(status = "crashed", up = NULL, message = clean(conditionMessage(e))))
    out$seconds <- elapsed() - t0
    swallowed <- grep("Error processing up-regulated genes", route_warnings, value = TRUE, fixed = TRUE)
    if (identical(out$status, "ok") && length(swallowed)) {
        out$status <- "crashed"
        out$message <- clean(swallowed[[1L]])
    } else if (identical(out$status, "ok") && (is.null(out$up) || !nrow(out$up@result))) {
        out$status <- "empty"
    }
    out
}

# ---- the q-value chain repeated (enrichit 0.2.0 R/utilities.R:14-73; KIT-189 5.5 (g)) ----------------------------
qvalueChain <- function(p) {
    q <- rep(NA_real_, length(p))
    ok <- is.finite(p) & !is.na(p) & p >= 0 & p <= 1
    pp <- p[ok]
    if (!length(pp)) return(list(stage = "none_valid", q = q, pi0 = NA_real_, messages = ""))
    if (all(pp == 1)) {
        q[ok] <- 1
        return(list(stage = "all_one", q = q, pi0 = NA_real_, messages = ""))
    }
    if (all(pp == 0)) {
        q[ok] <- 0
        return(list(stage = "all_zero", q = q, pi0 = NA_real_, messages = ""))
    }
    stages <- list(
        function() qvalue::qvalue(pp),
        function() qvalue::qvalue(pp, lambda = 0.05, pi0.method = "bootstrap"),
        function() qvalue::qvalue(pp, lambda = seq(0, 0.9, 0.05), pi0.method = "bootstrap"),
        function() qvalue::qvalue(pp, lambda = seq(0, 0.95, 0.05), pi0.method = "bootstrap"))
    messages <- character()
    for (i in seq_along(stages)) {
        obj <- suppressWarnings(tryCatch(stages[[i]](), error = function(e) e))
        if (inherits(obj, "qvalue")) {
            q[ok] <- obj$qvalues
            return(list(stage = as.character(i), q = q, pi0 = obj$pi0, messages = clean(paste(messages, collapse = " | "))))
        }
        messages <- c(messages, paste0(i, ": ", if (inherits(obj, "error")) conditionMessage(obj) else "no qvalue object"))
    }
    list(stage = "na", q = q, pi0 = NA_real_, messages = clean(paste(messages, collapse = " | ")))
}
ratioOf <- function(q, adj) {
    keep <- !is.na(q) & !is.na(adj) & adj > 0
    if (!any(keep)) return(c(value = NA_real_, spread = NA_real_))
    r <- q[keep] / adj[keep]
    c(value = r[[1L]], spread = max(abs(r - r[[1L]])) / r[[1L]])
}
relDiff <- function(a, b) if (!length(a)) 0 else max(abs(a - b) / pmax(abs(b), 1e-300))

# ---- the truth (KIT-189 5.2): exact integer comparisons, ER as a number for the record -----------------------------
truthOf <- function(v, vt, u, t, L, s) {
    # doubles: the products reach about 4e12, exact in double precision (below 2^53), beyond the integer range
    v <- as.numeric(v)
    vt <- as.numeric(vt)
    u <- as.numeric(u)
    t <- as.numeric(t)
    L <- as.numeric(L)
    s <- as.numeric(s)
    if (t == 0) {
        expected0 <- L * v / u
        er <- ifelse(v > 0, expected0 / expected0, NA_real_)
        return(list(er = er, nonnull = rep(FALSE, length(v)), material = rep(FALSE, length(v))))
    }
    lhs <- s * vt * (u - t) * u + (L - s) * (v - vt) * t * u
    rhs <- L * v * t * (u - t)
    expected <- s * vt / t + (L - s) * (v - vt) / (u - t)
    er <- ifelse(v > 0, expected / (L * v / u), NA_real_)
    list(er = er, nonnull = v > 0 & lhs > rhs, material = v > 0 & lhs > rhs & 2 * lhs >= 3 * rhs)
}

# ---- one call's 7 rows (KIT-189 5.3, 5.5, 7.4 E2) ------------------------------------------------------------------
emptyRows <- function(base, status, message) {
    data.frame(base, arm = ARMS, status = status, message = message, m = NA_integer_, m_kept = 0L,
        calls = if (identical(status, "empty")) 0L else NA_integer_,
        calls_null = if (identical(status, "empty")) 0L else NA_integer_,
        calls_nonnull = if (identical(status, "empty")) 0L else NA_integer_,
        calls_below15 = if (identical(status, "empty")) 0L else NA_integer_,
        calls_material = if (identical(status, "empty")) 0L else NA_integer_,
        planted_tested = FALSE, planted_called = FALSE, n_p_ge_05 = 0L, n_p_ge_095 = 0L, n_p_eq_1 = 0L,
        q_all_na = FALSE, pi0_hat = NA_real_, pi0_spread = NA_real_, chain_stage = NA_character_,
        chain_pi0 = NA_real_, chain_messages = NA_character_, chain_equal = NA, median_er_null = NA_real_,
        er_dev = NA_real_, e2_ok = NA, e2_fail = NA_character_, stringsAsFactors = FALSE)
}
analyseCall <- function(rc, bg, draw, L, base) {
    base$seconds <- round(rc$seconds, 2)
    if (!identical(rc$status, "ok")) {
        rows <- emptyRows(base[setdiff(names(base), "seconds")], rc$status, rc$message)
        rows$seconds <- base$seconds
        return(rows)
    }
    up <- rc$up
    res <- up@result
    gs <- up@geneSets
    ids <- as.character(res$ID)
    p <- as.numeric(res$pvalue)
    adj <- as.numeric(res$p.adjust)
    q <- as.numeric(res$qvalue)
    m <- length(gs)
    mk <- nrow(res)
    ext <- c(p, rep(1, m - mk))
    bhf <- stats::p.adjust(ext, "BH")[seq_len(mk)]
    qf_full <- suppressWarnings(calculateQvalue(ext))
    qf <- qf_full[seq_len(mk)]
    called <- list(Q = !is.na(q) & q < CUT, BH = adj < CUT, BHf = bhf < CUT, Qf = !is.na(qf) & qf < CUT,
        BY = stats::p.adjust(p, "BY") < CUT, Holm = stats::p.adjust(p, "holm") < CUT, raw = p < CUT)
    # truth per tested term
    sets_k <- gs[ids]
    len <- lengths(sets_k)
    genes <- unlist(sets_k, use.names = FALSE)
    grp <- rep.int(seq_along(ids), len)
    v <- tabulate(grp[genes %in% bg$U], nbins = mk)
    vt <- tabulate(grp[genes %in% draw$T], nbins = mk)
    truth <- truthOf(v, vt, length(bg$U), length(draw$T), L, draw$s)
    null_term <- !truth$nonnull
    # mechanism
    q_all_na <- all(is.na(q))
    chain <- qvalueChain(p)
    chain_f <- qvalueChain(ext)
    rq <- ratioOf(q, adj)
    rqf <- ratioOf(qf, bhf)
    # E2
    universe_n <- length(up@universe)
    n_list <- length(intersect(up@gene, up@universe))
    ph <- stats::phyper(res$Count - 1, len, universe_n - len, n_list, lower.tail = FALSE)
    gene_ratio <- strsplit(as.character(res$GeneRatio), "/", fixed = TRUE)
    bg_ratio <- strsplit(as.character(res$BgRatio), "/", fixed = TRUE)
    margins_same <- all(as.numeric(vapply(gene_ratio, `[`, "", 2L)) == n_list) &&
        all(as.numeric(vapply(bg_ratio, `[`, "", 1L)) == len) && all(as.numeric(vapply(bg_ratio, `[`, "", 2L)) == universe_n)
    checks <- c(
        margins = margins_same,
        bh = relDiff(adj, stats::p.adjust(p, "BH")) <= 1e-12,
        pi0 = is.na(rq[["spread"]]) || rq[["spread"]] <= 1e-12,
        chain = identical(as.numeric(chain$q), q),
        chain_f = identical(as.numeric(chain_f$q), as.numeric(qf_full)),
        count = all(res$Count > 0) && mk == sum(res$Count > 0),
        phyper = relDiff(p, ph) <= 1e-12,
        nest_bhf = all(!called$BHf | called$BH),
        nest_bh = q_all_na || all(!called$BH | called$Q))
    e2_fail <- paste(names(checks)[!checks], collapse = ";")
    planted_at <- if (is.na(draw$term)) integer() else which(ids == draw$term)
    rows <- lapply(ARMS, function(arm) {
        hit <- called[[arm]]
        is_q <- arm %in% c("Q", "Qf")
        ch <- if (identical(arm, "Q")) chain else chain_f
        rr <- if (identical(arm, "Q")) rq else rqf
        data.frame(base[setdiff(names(base), "seconds")], arm = arm, status = "ok", message = "", m = m, m_kept = mk,
            calls = sum(hit), calls_null = sum(hit & null_term), calls_nonnull = sum(hit & truth$nonnull),
            calls_below15 = sum(hit & !truth$material), calls_material = sum(hit & truth$material),
            planted_tested = length(planted_at) > 0L, planted_called = length(planted_at) > 0L && any(hit[planted_at]),
            n_p_ge_05 = sum(p >= 0.5), n_p_ge_095 = sum(p >= 0.95), n_p_eq_1 = sum(p == 1), q_all_na = q_all_na,
            pi0_hat = if (is_q) unname(rr[["value"]]) else NA_real_,
            pi0_spread = if (is_q) unname(rr[["spread"]]) else NA_real_,
            chain_stage = if (is_q) ch$stage else NA_character_, chain_pi0 = if (is_q) ch$pi0 else NA_real_,
            chain_messages = if (is_q) ch$messages else NA_character_,
            chain_equal = if (identical(arm, "Q")) checks[["chain"]] else if (identical(arm, "Qf")) checks[["chain_f"]] else NA,
            median_er_null = if (length(draw$T) && any(null_term & !is.na(truth$er))) stats::median(truth$er[null_term], na.rm = TRUE) else NA_real_,
            er_dev = if (!length(draw$T) && any(!is.na(truth$er))) max(abs(truth$er - 1), na.rm = TRUE) else NA_real_,
            e2_ok = all(checks), e2_fail = e2_fail, seconds = base$seconds, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# ---- one replicate: draw, two calls, 14 rows, the part written then renamed (KIT-189 7.1, R01) ---------------------
cellRow <- function(idx) CELLS[CELLS$idx == idx, , drop = FALSE]
partPath <- function(idx, r) file.path(rows_dir, sprintf("part-%d-%d.tsv", idx, r))
runReplicate <- function(task) {
    t0 <- elapsed()
    cell <- as.list(cellRow(task$cell))
    tryCatch({
        draw <- drawReplicate(cell, task$r)
        bg <- BG[[cell$bg]]
        rows <- do.call(rbind, lapply(VALUES, function(value) {
            rc <- runRoute(bg, draw$list, value)
            base <- data.frame(part = part, cell = cell$idx, background = cell$bg, scenario = cell$scenario, L = cell$L,
                replicate = task$r, seed = sprintf("%.0f", draw$seed),
                list_sha = digest::digest(draw$list, algo = "sha256"), propagation = value,
                planted_term = draw$term, planted_size = length(draw$T), planted_s = draw$s, seconds = 0,
                stringsAsFactors = FALSE)
            analyseCall(rc, bg, draw, cell$L, base)
        }))
        target <- partPath(cell$idx, task$r)
        tmp <- paste0(target, ".tmp")
        utils::write.table(rows, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        stopifnot(file.rename(tmp, target))
        list(cell = cell$idx, r = task$r, ok = TRUE, message = "", wall = elapsed() - t0,
            statuses = table(factor(rows$status[rows$arm == "Q"], levels = c("ok", "empty", "crashed"))),
            e2_fail = sum(rows$e2_ok[rows$arm == "Q"] %in% FALSE), er_dev = suppressWarnings(max(rows$er_dev, na.rm = TRUE)),
            call_seconds = stats::setNames(rows$seconds[rows$arm == "Q"], rows$propagation[rows$arm == "Q"]))
    }, error = function(e) list(cell = cell$idx, r = task$r, ok = FALSE, message = clean(conditionMessage(e)),
        wall = elapsed() - t0))
}

# ---- run ------------------------------------------------------------------------------------------------------------
reps <- first_rep:last_rep
if (!is.null(replicate_list)) reps <- reps[reps %in% replicate_list]
tasks <- unlist(lapply(cells_run, function(idx) lapply(reps, function(r) list(cell = idx, r = r))), recursive = FALSE)
writeLines(sprintf("job %s: part %s, cells %s, replicates %d to %d (%d to run), %d workers; set-up %.1f s",
    job, part, paste(cells_run, collapse = ","), first_rep, last_rep, length(reps), workers, elapsed() - t_start))
t_par <- elapsed()
summaries <- if (workers > 1L && length(tasks) > 1L) {
    parallel::mclapply(tasks, runReplicate, mc.cores = workers, mc.preschedule = TRUE)
} else lapply(tasks, runReplicate)
wall_par <- elapsed() - t_par
okOf <- function(x) is.list(x) && isTRUE(x$ok)
n_ok <- sum(vapply(summaries, okOf, logical(1)))
writeLines(sprintf("job %s: %d of %d replicates written in %.1f s", job, n_ok, length(tasks), wall_par))

# assemble the parts of the replicates this run ran, per cell
assembled <- list()
for (idx in cells_run) {
    parts <- partPath(idx, reps)
    parts <- parts[file.exists(parts)]
    if (!length(parts)) next
    target <- file.path(rows_dir, sprintf("%s-%d-%d-%d.tsv.gz", part, idx, first_rep, last_rep))
    header <- readLines(parts[[1L]], n = 1L)
    body <- unlist(lapply(parts, function(f) readLines(f)[-1L]), use.names = FALSE)
    con <- gzfile(target, "w")
    writeLines(c(header, body), con)
    close(con)
    file.remove(parts)
    assembled[[as.character(idx)]] <- data.frame(cell = idx, file = basename(target), replicates = length(parts),
        lines = length(body), sha256 = sha256File(target), stringsAsFactors = FALSE)
}
assembled <- do.call(rbind, assembled)
failures <- do.call(rbind, lapply(summaries, function(x) if (okOf(x)) NULL else data.frame(
    cell = if (is.list(x)) x$cell else NA, replicate = if (is.list(x)) x$r else NA,
    message = if (is.list(x)) x$message else clean(paste(as.character(x), collapse = " ")), stringsAsFactors = FALSE)))

# ---- part l1: the blind lifecycle record and the checks E0 to E4, E6 (KIT-189 7.4) ---------------------------------
checks <- NULL
if (lifecycle) {
    oks <- Filter(okOf, summaries)
    per_rep <- do.call(rbind, lapply(oks, function(x) data.frame(cell = x$cell, replicate = x$r, wall = x$wall,
        ok = x$statuses[["ok"]], empty = x$statuses[["empty"]], crashed = x$statuses[["crashed"]],
        e2_fail = x$e2_fail, er_dev = x$er_dev, none_seconds = x$call_seconds[["none"]],
        prop_seconds = x$call_seconds[["is_a_part_of_regulates"]], stringsAsFactors = FALSE)))
    cell_table <- do.call(rbind, lapply(cells_run, function(idx) {
        x <- per_rep[per_rep$cell == idx, , drop = FALSE]
        a <- if (is.null(assembled)) NULL else assembled[assembled$cell == idx, , drop = FALSE]
        data.frame(cellRow(idx)[, c("idx", "bg", "scenario", "L")], replicates_written = nrow(x),
            calls_ok = sum(x$ok), calls_empty = sum(x$empty), calls_crashed = sum(x$crashed), e2_failed_calls = sum(x$e2_fail),
            rep_wall_median = stats::median(x$wall), none_seconds_median = stats::median(x$none_seconds),
            prop_seconds_median = stats::median(x$prop_seconds),
            rows_sha256 = if (is.null(a) || !nrow(a)) NA_character_ else a$sha256, stringsAsFactors = FALSE)
    }))
    utils::write.table(cell_table, file.path(out_dir, "lifecycle-cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    timing <- data.frame(startup_seconds = t_loaded - t_start,
        setup_B1 = BG$B1$seconds, setup_B2 = BG$B2$seconds, setup_B3 = BG$B3$seconds,
        replicates = length(tasks), written = n_ok, parallel_wall = wall_par, sum_rep_wall = sum(per_rep$wall),
        speedup = sum(per_rep$wall) / wall_par, workers = workers,
        cell2_rep_wall_median = stats::median(per_rep$wall[per_rep$cell == 2L]), stringsAsFactors = FALSE)
    utils::write.table(timing, file.path(out_dir, "lifecycle-timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

    checkRow <- function(id, pass, detail) data.frame(id = id, pass = pass, detail = clean(detail), stringsAsFactors = FALSE)
    # E1: the installed code against the fetched sources (KIT-189 2.1)
    e1 <- tryCatch({
        work <- tempfile("k189-e1-")
        dir.create(work)
        archives <- list(
            enrichit = c("https://cran.r-project.org/src/contrib/Archive/enrichit/enrichit_0.2.0.tar.gz",
                "fb2192cc3beca97e08618ac01e28a1ecc4adb93a2b3a5078da60413eb9234107"),
            qvalue = c("https://bioconductor.org/packages/3.23/bioc/src/contrib/qvalue_2.44.0.tar.gz",
                "e87a8b32fa7fa9021956212ab631a7c55fe410b724c0a54d99e1afc9f8fbbd0e"),
            clusterProfiler = c(
                "https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/clusterProfiler/clusterProfiler_4.20.0.tar.gz",
                "0450496ea91145af91eb038ff1e9f4e1b5dd3769b0de6fd53f6a01ce8d1aee78"))
        for (pkg in names(archives)) {
            dest <- file.path(work, basename(archives[[pkg]][[1L]]))
            utils::download.file(archives[[pkg]][[1L]], dest, mode = "wb", quiet = TRUE)
            got <- sha256File(dest)
            if (!identical(got, archives[[pkg]][[2L]])) stop(sprintf("E1 archive hash differs for %s: %s", pkg, got))
            utils::untar(dest, exdir = work)
        }
        file_hashes <- c(
            "enrichit/R/ora_gson.R" = "0404dd407870d287b8d47e8735f1f13e7b01a6e36a0f2dc0ba178b6eb8cb0d06",
            "enrichit/R/utilities.R" = "9a7ac4a035976ec06bcdd85f96341897e405c3121fae3e94d57493720fab4877",
            "enrichit/src/enrichit.cpp" = "1b478d9ae1f2d5cad836013f3d07dbd6fc141f23a1c5b90ec22adfd690ad1d83",
            "qvalue/R/qvalue.R" = "8831b60897542b6603f9922fce883357202143872200efbec60dcfdf78d9f655",
            "qvalue/R/pi0est.R" = "57b456ff997d44a044de565c11287eb9fe9b634952904ac43127f271d211115d")
        file_ok <- vapply(names(file_hashes), function(f) identical(sha256File(file.path(work, f)), file_hashes[[f]]),
            logical(1))
        targets <- list(c("enrichit", "R/ora_gson.R", "ora_gson"), c("enrichit", "R/utilities.R", "calculate_qvalue"),
            c("qvalue", "R/qvalue.R", "qvalue"), c("qvalue", "R/pi0est.R", "pi0est"), c("qvalue", "R/lfdr.R", "lfdr"),
            c("clusterProfiler", "R/enricher.R", "enricher"))
        compared <- vapply(targets, function(tg) {
            exprs <- parse(file.path(work, tg[[1L]], tg[[2L]]), keep.source = FALSE)
            hit <- Filter(function(e) is.call(e) && (identical(e[[1L]], as.name("<-")) || identical(e[[1L]], as.name("="))) &&
                identical(e[[2L]], as.name(tg[[3L]])) && is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function")),
                as.list(exprs))
            if (length(hit) != 1L) stop(sprintf("E1: %d definitions of %s in %s", length(hit), tg[[3L]], tg[[2L]]))
            # the function expression evaluated in a new environment (no bindings; base as its parent); no body is run
            from_source <- eval(hit[[1L]][[3L]], new.env(parent = baseenv()))
            installed <- get(tg[[3L]], envir = asNamespace(tg[[1L]]))
            same <- identical(deparse(formals(from_source)), deparse(formals(installed))) &&
                identical(deparse(body(from_source)), deparse(body(installed)))
            same
        }, logical(1))
        names(compared) <- vapply(targets, function(tg) paste0(tg[[1L]], "::", tg[[3L]]), character(1))
        checkRow("E1", all(compared), paste0("archives hash ok; files ", paste(names(file_ok), file_ok, sep = "=", collapse = ", "),
            "; deparse identical: ", paste(names(compared), compared, sep = "=", collapse = ", ")))
    }, error = function(e) checkRow("E1", NA, paste("instrument:", conditionMessage(e))))

    # E2: over every L1 call (the per-call flags written with the rows)
    l1_rows <- do.call(rbind, lapply(assembled$file, function(f) utils::read.delim(file.path(rows_dir, f),
        colClasses = "character", na.strings = "NA", quote = "")))
    calls_q <- l1_rows[l1_rows$arm == "Q" & l1_rows$status == "ok", , drop = FALSE]
    fail_kinds <- unlist(strsplit(calls_q$e2_fail[nzchar(calls_q$e2_fail) & !is.na(calls_q$e2_fail)], ";", fixed = TRUE))
    e2 <- checkRow("E2", nrow(calls_q) > 0L && all(calls_q$e2_ok == "TRUE"),
        sprintf("%d ok calls checked (of %d calls); failing calls %d; by check: %s", nrow(calls_q),
            sum(l1_rows$arm == "Q"), sum(calls_q$e2_ok != "TRUE"),
            if (length(fail_kinds)) paste(names(table(fail_kinds)), table(fail_kinds), sep = "=", collapse = ", ") else "none"))
    # pairing: both propagation values received the same list in every replicate
    pairs <- tapply(l1_rows$list_sha, paste(l1_rows$cell, l1_rows$replicate), function(x) length(unique(x)))
    pairing <- checkRow("pairing", all(pairs == 1L), sprintf("%d replicates; one list per replicate in %d", length(pairs),
        sum(pairs == 1L)))

    # E0: KIT-173's null draw 5 (effect-measurement.R drawList, arm 4, draw 5) under "none"
    e0 <- tryCatch({
        b1 <- BG$B1
        entries173 <- sort(b1$raw$Entry)
        prop <- b1$sets$is_a_part_of_regulates
        planted173 <- lapply(c("GO:0032790", "GO:0006396", "GO:0044282"), function(term)
            sort(unique(prop$term2gene[[2L]][prop$term2gene[[1L]] == term])))
        null_size <- as.integer(stats::median(2L * (lengths(planted173) %/% 2L)))
        seedKit(20261123L + 100L * 4L + 5L)
        list173 <- sort(sample(entries173, null_size))
        rc <- runRoute(b1, list173, "none")
        draw173 <- list(seed = NA, list = list173, term = NA_character_, T = character(), s = 0L)
        base <- data.frame(part = "e0", cell = 0L, background = "B1", scenario = "null", L = null_size, replicate = 5L,
            seed = "20261528", list_sha = "", propagation = "none", planted_term = NA_character_, planted_size = 0L,
            planted_s = 0L, seconds = 0, stringsAsFactors = FALSE)
        r0 <- analyseCall(rc, b1, draw173, null_size, base)
        q0 <- r0[r0$arm == "Q", ]
        b0 <- r0[r0$arm == "BH", ]
        checkRow("E0", identical(rc$status, "ok") && q0$m_kept == 156L && q0$calls == 120L && b0$calls == 0L,
            sprintf("list of %d; status %s; tested %s, q-value calls %s, BH calls %s (expected 156, 120, 0); pi0_hat %s; stage %s (chain pi0 %s; earlier stages: %s)",
                null_size, rc$status, q0$m_kept, q0$calls, b0$calls, format(q0$pi0_hat, digits = 6), q0$chain_stage,
                format(q0$chain_pi0, digits = 6), q0$chain_messages))
    }, error = function(e) checkRow("E0", NA, paste("instrument:", conditionMessage(e))))

    # E3: route parity against a direct clusterProfiler::enricher call, one list per background and value
    e3 <- tryCatch({
        details <- character()
        same_all <- TRUE
        for (code in c("B1", "B2", "B3")) {
            idx <- c(B1 = 2L, B2 = 5L, B3 = 8L)[[code]]
            draw <- drawReplicate(as.list(cellRow(idx)), first_rep)
            bg <- BG[[code]]
            for (value in VALUES) {
                rc <- runRoute(bg, draw$list, value)
                direct <- suppressWarnings(suppressMessages(clusterProfiler::enricher(gene = draw$list,
                    universe = bg$entries, TERM2GENE = bg$sets[[value]]$term2gene, TERM2NAME = bg$sets[[value]]$term2name,
                    pvalueCutoff = CUT, pAdjustMethod = "BH")))
                a <- rc$up@result
                b <- direct@result
                cols <- setdiff(names(a), "geneID")
                gene_sets_same <- identical(names(a), names(b)) && identical(rownames(a), rownames(b)) &&
                    all(mapply(function(x, y) setequal(strsplit(x, "/", fixed = TRUE)[[1L]], strsplit(y, "/", fixed = TRUE)[[1L]]),
                        as.character(a$geneID), as.character(b$geneID)))
                same <- identical(rc$status, "ok") && identical(a[, cols], b[, cols]) && isTRUE(gene_sets_same)
                same_all <- same_all && same
                details <- c(details, sprintf("%s/%s %s (%d terms)", code, value, if (same) "identical" else "DIFFERENT", nrow(a)))
            }
        }
        checkRow("E3", same_all, paste(details, collapse = "; "))
    }, error = function(e) checkRow("E3", NA, paste("instrument:", conditionMessage(e))))

    # E4: A-BHf and A-Qf on a constructed vector against the hand values
    e4 <- tryCatch({
        p4 <- c(0.001, 0.01, 0.2)
        ext4 <- c(p4, rep(1, 10L - 3L))
        bhf4 <- stats::p.adjust(ext4, "BH")[1:3]
        qf4 <- suppressWarnings(calculateQvalue(ext4))[1:3]
        hand <- c(0.01, 0.05, 2 / 3)
        chain4 <- qvalueChain(ext4)
        ok4 <- isTRUE(all.equal(bhf4, hand, tolerance = 1e-14)) && isTRUE(all.equal(qf4, hand, tolerance = 1e-14)) &&
            identical(chain4$stage, "1")
        checkRow("E4", ok4, sprintf("BHf %s; Qf %s; hand %s (Qf = pi0 x BHf with pi0 = 1: every pi0(lambda) at lambda >= 0.5 exceeds 1); Qf stage %s, pi0 %s",
            paste(format(bhf4, digits = 16), collapse = ","), paste(format(qf4, digits = 16), collapse = ","),
            paste(format(hand, digits = 16), collapse = ","), chain4$stage, format(chain4$pi0, digits = 16)))
    }, error = function(e) checkRow("E4", NA, paste("instrument:", conditionMessage(e))))

    # E6: the truth function
    e6 <- tryCatch({
        planted6 <- truthOf(v = c(20, 10), vt = c(20, 0), u = 100, t = 20, L = 10, s = 5L)
        null6 <- truthOf(v = c(10, 37, 500), vt = c(0, 0, 0), u = 4000, t = 0, L = 126, s = 0L)
        null_dev <- suppressWarnings(max(as.numeric(l1_rows$er_dev[l1_rows$scenario == "null"]), na.rm = TRUE))
        ok6 <- abs(planted6$er[[1L]] - 2.5) <= 1e-12 && isTRUE(planted6$nonnull[[1L]]) && isTRUE(planted6$material[[1L]]) &&
            abs(planted6$er[[2L]] - 0.625) <= 1e-12 && !planted6$nonnull[[2L]] && all(null6$er == 1) && !any(null6$nonnull) &&
            is.finite(null_dev) && null_dev <= 1e-12
        checkRow("E6", ok6, sprintf("planted case ER of T %s (hand 2.5), of a disjoint term %s (hand 0.625); constructed null ER %s; largest |ER - 1| over the L1 null calls %s",
            format(planted6$er[[1L]], digits = 16), format(planted6$er[[2L]], digits = 16),
            paste(null6$er, collapse = ","), format(null_dev, digits = 6)))
    }, error = function(e) checkRow("E6", NA, paste("instrument:", conditionMessage(e))))

    checks <- rbind(e0, e1, e2, pairing, e3, e4, e6)
    utils::write.table(checks, file.path(out_dir, "checks-l1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(sprintf("check %s: %s", checks$id, checks$pass))
}

# ---- provenance ----------------------------------------------------------------------------------------------------
tables <- c(list.files(out_dir, pattern = "\\.(tsv|gz)$", full.names = TRUE, recursive = TRUE),
    if (lifecycle) file.path(rows_dir, assembled$file) else character())
versions <- vapply(c("MultiScholaR", "clusterProfiler", "enrichit", "qvalue", "GO.db", "AnnotationDbi"),
    function(p) as.character(utils::packageVersion(p)), character(1))
cpu <- tryCatch(sub("^model name\\s*:\\s*", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[[1L]]),
    error = function(e) NA_character_)
provenance <- list(job = job, part = part, cells = cells_run, first = first_rep, last = last_rep,
    replicate_list = if (is.null(replicate_list)) NULL else replicate_list, replicates_run = length(reps),
    seeds = "189000000 + 100000 c + r (Mersenne-Twister, Inversion, Rejection); B2 189000001",
    commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_),
    r_version = R.version.string, versions = as.list(versions), blas = extSoftVersion()[["BLAS"]],
    lapack = La_library(), runner_image = paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion")), cpu = cpu,
    workers = workers, start = format(started_at, "%Y-%m-%dT%H:%M:%S%z"), end = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    replicates_written = n_ok, replicate_failures = if (is.null(failures)) 0L else nrow(failures),
    exit_status = if (is.null(failures)) 0L else 1L,
    script_sha256 = sha256File(file.path("docs", "qa", "workflow-execution", "kit", "kit-189", "sim-k189.R")),
    tables = stats::setNames(as.list(vapply(tables, sha256File, character(1))), basename(tables)))
jsonlite::write_json(provenance, file.path(out_dir, sprintf("provenance-%s.json", job)), auto_unbox = TRUE,
    pretty = TRUE, digits = NA)
if (!is.null(failures)) {
    utils::write.table(failures, file.path(out_dir, sprintf("failures-%s.tsv", job)), sep = "\t", quote = FALSE,
        row.names = FALSE)
}
restore_widget()
restore_ggsave()
if (!is.null(failures)) {
    writeLines(sprintf("job %s: %d replicates failed in the script (see failures-%s.tsv)", job, nrow(failures), job))
    quit(status = 1L)
}
