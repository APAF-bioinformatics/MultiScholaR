# KIT-196 script S1 (ticket section 9 "S1", acceptance criteria A2 (i) to (iv), and D6): KIT-188's lipidomics and
# metabolomics journeys NO and OR on data sets 401 to 403, base against head. Evidence code under docs/qa; no package code.
#
# Usage, from the package root of the commit under test on a run-tests runner (full git history, the KIT-188 branch fetched):
#   Rscript docs/qa/workflow-execution/kit/kit-196/da-sparse-control.R <part L|M> <base commit> [--base-only]
#       [--base-records=<dir>] [<out dir>]
#   --base-records=<dir>  the base side is not run again: its records are read from <dir>/s1-<part>-base.tsv, written by
#       this script's base side on the same <base commit> (R2 reuses R1's; the base side's code is unchanged)
#   part L: cells L12 (deciding) and L64 (deciding, ba); part M: M12 (deciding) and M168 (deciding, ba)
#   <out dir>  default docs/qa/workflow-execution/kit/kit-196/s1 (made absolute)
# Child process (started by the parent):
#   Rscript da-sparse-control.R --unit <side> <tree> <commit> <cell> <data set> <job> <rds>
#
# Mechanism (ticket 9, S1, revision 1.2):
# - Base is a detached git worktree of <base commit>; head is this tree. Per side, the parent prepares the tree as the
#   run-tests jobs do (src/build_dia_rsum.R when the DIA manifest is absent) and compiles its native code once
#   (pkgbuild::compile_dll, in its own process), before any child starts.
# - Each unit (one cell, data set and job, with every condition of KIT-188's job: sim-k188.R conditionsOf) runs in a fresh
#   Rscript child, four at a time (parallel::mclapply in the parent, which loads no package and owns no session). The child
#   loads its side's tree (pkgload::load_all(tree, export_all = TRUE, helpers = FALSE, compile = FALSE)), checks the
#   SHA-256 of sim-k188-69fd01bc.R, evaluates its definitions (everything before its "child process" section), and asserts
#   that MultiScholaR's namespace is its side's tree. Every row carries the side's commit.
# - Per condition, the NO journey then the OR journey, exactly as sim-k188.R's runConditionL makes them (journeyL, the
#   defaults routeDefaults(omic), T_ids from the first journey's kept features). The DA step's call is wrapped in the global
#   environment (journeyL calls it there) to record, without changing it, the DA input of every assay (the head state's
#   assay matrices as the per-contrast kernel builds them: binding$da_assay_input) and the DA result or its error.
# Records per journey and assay: status; the full stop's message; the manifest's errors per contrast and assay; the count
# of features without values in the DA input (no value in any fitted sample); whether the positive residual df are
# unequal; SHA-256 of the journey's da table and of its long table; for every assay with an error, the fit half and
# contrast half re-run directly (binding$da_fit, binding$da_apply with TRUE, TRUE) in this script's own tryCatch, with the
# condition's message and call; on head, the comparison of the assay's tested rows with the same halves on the DA input
# without the features without values (A2 (iii)), and D6 (zero-df features, descriptive, each computation in its own
# tryCatch).
# Written to <out>: s1-<part>-<side>.tsv, s1-units-<part>-<side>.tsv, oracle-stops-<part>.tsv (A2 (i)),
# partial-failures-<part>.tsv, and with head: compare-<part>.tsv (journey, field, base, head, equal, criterion),
# d6-zero-df-<part>.tsv, d6-summary-<part>.tsv; summary-<part>.txt. Synthetic data only.
# Exit status non-zero when a side's preparation fails or a child fails to return its records.

S1_COPY <- "sim-k188-69fd01bc.R"
S1_COPY_SHA <- "c0fbc4ab34925dbbbfee90169dee3e7af60f24b3d29c1ba94aa83175b52a1c22"
S1_K188_COMMIT <- "69fd01bc4172331daf902642a29517bbf42995c9"
S1_K188_ROWS <- "docs/qa/workflow-execution/kit/kit-188/l1/%s/rows/l-%s-401-403-l1.tsv.gz"
S1_PARTS <- list(
    L = list(dir = "l", cells = list(L12 = "deciding", L64 = c("deciding", "ba"))),
    M = list(dir = "m", cells = list(M12 = "deciding", M168 = c("deciding", "ba"))))
S1_REPS <- 401:403
S1_WORKERS <- 4L
S1_COLUMNS <- c("side", "commit", "cell", "job", "scenario", "mask", "rep", "arm", "status", "class", "step", "seconds",
    "full_stop", "stop_message", "da_sha", "long_sha", "miss_share", "assay", "n_T", "n_features", "n_samples",
    "n_without_values", "without_values", "unequal_df", "err_contrasts", "err_message", "direct_message", "direct_fn",
    "direct_call", "reduced_has_na", "cmp_status", "cmp_max_rel", "cmp_ids_ok", "untested_ok", "untested_counts")
S1_D6_COLUMNS <- c("side", "commit", "cell", "job", "scenario", "mask", "rep", "arm", "assay", "contrast", "error",
    "n_features", "n_without_values", "n_zero_df", "ref_n", "alt_n", "ref_unequal_df", "alt_unequal_df",
    "ref_downweighted", "alt_downweighted", "ref_s2prior_median", "alt_s2prior_median", "ref_dfprior_min",
    "ref_dfprior_median", "ref_dfprior_max", "alt_dfprior_min", "alt_dfprior_median", "alt_dfprior_max", "n_both",
    "p_absdiff_max", "p_absdiff_median", "calls_ref", "calls_alt", "calls_only_ref", "calls_only_alt", "zero_df_called_ref")
S1_KEY <- c("cell", "job", "scenario", "mask", "rep", "arm")

s1_argv <- commandArgs(trailingOnly = TRUE)
s1_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
s1_dir <- dirname(s1_file)
s1_rscript <- file.path(R.home("bin"), "Rscript")

s1Sha <- function(path) digest::digest(file = path, algo = "sha256")
s1Flat <- function(x, n = 500L) substr(gsub("[\r\n\t]+", " ", paste(as.character(x), collapse = " ")), 1L, n)
s1Row <- function(values, columns) {
    row <- stats::setNames(as.list(rep(NA, length(columns))), columns)
    given <- intersect(names(values), columns)
    row[given] <- lapply(values[given], function(v) if (is.null(v) || !length(v)) NA else v[[1L]])
    as.data.frame(row, stringsAsFactors = FALSE)
}
s1Bind <- function(rows, columns) {
    if (!length(rows)) return(as.data.frame(stats::setNames(replicate(length(columns), character(), simplify = FALSE), columns)))
    do.call(rbind, rows)
}
# A table's digest from its column values (attributes removed), so that no store or class attribute enters it.
s1TableSha <- function(table) {
    if (is.null(table)) return("")
    table <- as.data.frame(table)
    columns <- lapply(table, function(x) {
        if (is.factor(x)) x <- as.character(x)
        attributes(x) <- NULL
        x
    })
    digest::digest(list(names = names(table), rows = nrow(table), columns = columns), algo = "sha256")
}

# ---- child process: one unit ------------------------------------------------------------------------------------------
s1Child <- function(a) {
    stopifnot(length(a) == 8L)
    side <- a[[2L]]
    tree <- normalizePath(a[[3L]])
    commit <- a[[4L]]
    cell <- a[[5L]]
    r <- as.integer(a[[6L]])
    job <- a[[7L]]
    rds <- a[[8L]]
    copy <- file.path(s1_dir, S1_COPY)
    started <- proc.time()[[3L]]
    setwd(tree)
    suppressPackageStartupMessages(suppressMessages(pkgload::load_all(tree, export_all = TRUE, helpers = FALSE,
        compile = FALSE, quiet = TRUE)))
    if (!identical(s1Sha(copy), S1_COPY_SHA)) stop("the copy of sim-k188.R does not have its recorded SHA-256")
    lines <- readLines(copy)
    end <- grep("^# ---- child process", lines)[[1L]] - 1L
    for (expression in parse(text = lines[seq_len(end)], keep.source = FALSE)) eval(expression, envir = globalenv())
    ns_path <- normalizePath(getNamespaceInfo("MultiScholaR", "path"))
    if (!identical(ns_path, tree)) stop(sprintf("MultiScholaR is loaded from %s, not from the side's tree %s", ns_path, tree))
    s1InstallDaCapture()
    ctx <- list(side = side, commit = commit, cell = cell, job = job, rep = r)
    d <- routeDefaults(CELLS[[cell]]$omic)
    conditions <- conditionsOf("L", cell, job)
    rows <- list()
    d6 <- list()
    for (s in unique(conditions$scenario)) {
        sim <- generate(cell, s, r)
        for (m in conditions$mask[conditions$scenario == s]) {
            out <- s1Condition(sim, m, job, d, ctx)
            rows <- c(rows, out$rows)
            d6 <- c(d6, out$d6)
        }
    }
    list(rows = s1Bind(rows, S1_COLUMNS), d6 = s1Bind(d6, S1_D6_COLUMNS), seconds = proc.time()[[3L]] - started,
        ns_path = ns_path, commit = commit)
}

# The DA step's wrapper in the global environment (journeyL's lookup finds it before the attached namespace): the DA
# input of every assay before the step, then the namespace's function, unchanged; its result or the error that escaped it.
s1InstallDaCapture <- function() {
    capture <- new.env()
    assign("s1_capture", capture, envir = globalenv())
    original <- get("runWorkflowOmicsReferenceDaAnalysis", envir = asNamespace("MultiScholaR"))
    wrapper <- function(binding, manager, source, settings, contrasts) {
        capture$inputs <- tryCatch(s1DaInputs(binding, manager), error = function(e) list())
        capture$input_error <- NULL
        tryCatch({
            da <- original(binding, manager, source, settings, contrasts)
            capture$da <- da
            da
        }, error = function(e) {
            capture$stop <- conditionMessage(e)
            stop(e)
        })
    }
    assign("runWorkflowOmicsReferenceDaAnalysis", wrapper, envir = globalenv())
}

#' Per assay, the matrix the per-contrast kernel fits (omicsDaContrastOperation: binding$da_assay_input on the state's
#' assay table and design) with the arguments of its fit half.
s1DaInputs <- function(binding, manager) {
    object <- manager$getState(manager$current_state)
    design <- object@design_matrix
    tables <- methods::slot(object, binding$data_slot)
    inputs <- lapply(stats::setNames(names(tables), names(tables)), function(assay) {
        input <- suppressMessages(binding$da_assay_input(tables[[assay]], assay, design, object@sample_id,
            methods::slot(object, binding$id_slot), object@annotation_id_column))
        if (is.null(input)) return(NULL)
        list(m = input$expr_matrix, design = design, sample_id = object@sample_id, group_id = object@group_id,
            technical_replicate_id = object@technical_replicate_id)
    })
    Filter(Negate(is.null), inputs)
}
s1Fit <- function(binding, input, m) {
    suppressWarnings(suppressMessages(binding$da_fit(m, input$design, "~ 0 + group", input$sample_id,
        group_id = input$group_id, technical_replicate_id = input$technical_replicate_id)))
}
s1Apply <- function(binding, fit, contrast) suppressWarnings(suppressMessages(binding$da_apply(fit, contrast, NA, TRUE, TRUE)))
s1UnequalDf <- function(df) {
    positive <- df[!is.na(df) & df > 0]
    length(positive) > 0L && !identical(min(positive), max(positive))
}

#' One journey of journeyL, with the capture reset before it.
s1RunJourney <- function(run) {
    capture <- get("s1_capture", envir = globalenv())
    rm(list = ls(capture, all.names = TRUE), envir = capture)
    j <- run()
    list(j = j, inputs = capture$inputs, da = capture$da, stop = capture$stop)
}

#' One condition: the NO journey, then the OR journey on T_ids (sim-k188.R runConditionL, arms NO and OR).
s1Condition <- function(sim, mask, job, d, ctx) {
    values <- maskedValues(sim, mask)
    miss <- mean(is.na(values[sim$feat_ids, ]))
    requests <- REQUESTS_L(d)
    tmp <- function(arm) tempfile(sprintf("k188-%s-%s-", sim$cell, arm))
    key <- c(ctx, list(scenario = sim$scenario, mask = mask))
    no <- s1RunJourney(function() journeyL(sim, values, "NO", requests$NO, d, tmp("NO")))
    T_ids <- if (!is.null(no$j$kept)) setdiff(no$j$kept, sim$is_ids)
    out <- s1Records(no, "NO", T_ids, sim, miss, key)
    if (!is.null(T_ids)) {
        oracle <- maskedValues(sim, mask, oracle = TRUE)
        or <- s1RunJourney(function() journeyL(sim, oracle, "OR", list(gate = FALSE, ruv = NULL), d, tmp("OR"), filter = FALSE,
            itsd = FALSE, method = "none", correlation = FALSE, ids = T_ids))
        out_or <- s1Records(or, "OR", T_ids, sim, miss, key)
        out <- list(rows = c(out$rows, out_or$rows), d6 = c(out$d6, out_or$d6))
    }
    out
}

#' The records of one journey: one row per assay of its DA input (one row with assay NA when no DA input was reached).
s1Records <- function(run, arm, T_ids, sim, miss, key) {
    j <- run$j
    binding <- bindingOf(sim$p$omic)
    contrasts <- contrastsOf(sim)$contrasts
    id_column <- if (identical(sim$p$omic, "lipid")) "lipid_id" else "metabolite_id"
    long <- if (!is.null(run$da)) as.data.frame(run$da$results[[binding$long_table]]) else NULL
    errors <- if (!is.null(run$da)) lapply(run$da$manifest$results, `[[`, "errors") else list()
    full_stop <- identical(j$status, "crashed") && identical(j$step, "da")
    journey <- c(key, list(arm = arm, status = j$status, class = s1Flat(j$class %||% "", 300L), step = j$step %||% NA,
        seconds = round(j$seconds %||% NA, 1), full_stop = full_stop, stop_message = s1Flat(run$stop %||% ""),
        da_sha = s1TableSha(j$da), long_sha = s1TableSha(long), miss_share = miss))
    inputs <- run$inputs %||% list()
    if (!length(inputs)) return(list(rows = list(s1Row(journey, S1_COLUMNS)), d6 = list()))
    rows <- list()
    d6 <- list()
    for (assay in names(inputs)) {
        input <- inputs[[assay]]
        m <- input$m
        common <- intersect(colnames(m), as.character(input$design[[input$sample_id]]))
        observed <- rowSums(!is.na(m[, common, drop = FALSE]))
        without <- rownames(m)[observed == 0L]
        fit <- tryCatch(s1Fit(binding, input, m), error = function(e) e)
        row <- c(journey, list(assay = assay, n_T = if (is.null(T_ids)) NA_integer_ else sum(sim$assay_of[T_ids] == assay),
            n_features = nrow(m), n_samples = length(common), n_without_values = length(without),
            without_values = s1Flat(without), unequal_df = if (inherits(fit, "error")) NA else s1UnequalDf(fit$df.residual)))
        # the manifest's errors of this assay (a full stop: every contrast, the joined message)
        if (full_stop) {
            row$err_contrasts <- length(contrasts)
            row$err_message <- s1Flat(run$stop)
        } else {
            mine <- Filter(Negate(is.null), lapply(errors, function(e) e[[assay]]))
            row$err_contrasts <- length(mine)
            row$err_message <- s1Flat(unique(unlist(mine)))
        }
        # the direct re-run of the halves for an assay with an error (P02)
        if (row$err_contrasts > 0L) {
            if (inherits(fit, "error")) {
                row$direct_message <- s1Flat(paste("fit half:", conditionMessage(fit)))
            } else {
                first <- if (full_stop) contrasts[[1L]] else names(Filter(function(e) !is.null(e[[assay]]), errors))[[1L]]
                e <- tryCatch({
                    s1Apply(binding, fit, first)
                    NULL
                }, error = function(e) e)
                row$direct_message <- if (is.null(e)) "<completed>" else s1Flat(conditionMessage(e))
                call <- if (is.null(e)) NULL else conditionCall(e)
                row$direct_fn <- if (is.null(call)) "" else s1Flat(deparse(call[[1L]]), 200L)
                row$direct_call <- if (is.null(call)) "" else s1Flat(deparse(call), 300L)
            }
        }
        # head: the tested rows against the halves on the input without the features without values (A2 (iii)); D6
        if (identical(key$side, "head") && !is.null(long) && row$err_contrasts == 0L && !inherits(fit, "error")) {
            row <- c(row, s1Compare(binding, input, m, common, without, long, assay, contrasts, id_column))
            d6 <- c(d6, s1D6(binding, input, m, fit, without, contrasts, c(key, list(arm = arm, assay = assay))))
        }
        rows[[length(rows) + 1L]] <- s1Row(row, S1_COLUMNS)
    }
    list(rows = rows, d6 = d6)
}

#' A2 (iii): per contrast, the journey's rows with AveExpr present against the contrast half on the reduced input
#' (logFC, t, P.Value, fdr_qvalue), and its rows with AveExpr missing against the features without values.
s1Compare <- function(binding, input, m, common, without, long, assay, contrasts, id_column) {
    reduced <- m[!rownames(m) %in% without, , drop = FALSE]
    out <- list(reduced_has_na = anyNA(reduced[, common, drop = FALSE]), cmp_status = "not_run", cmp_max_rel = NA_real_,
        cmp_ids_ok = NA, untested_ok = NA, untested_counts = "")
    fit <- tryCatch(s1Fit(binding, input, reduced), error = function(e) e)
    if (inherits(fit, "error")) {
        out$cmp_status <- s1Flat(paste("reduced fit:", conditionMessage(fit)), 200L)
        return(out)
    }
    identical_all <- TRUE
    rel_max <- 0
    ids_ok <- TRUE
    untested_ok <- TRUE
    counts <- character()
    for (contrast in contrasts) {
        table <- long[long$assay == assay & long$comparison == contrast, , drop = FALSE]
        untested <- as.character(table[[id_column]][is.na(table$AveExpr)])
        counts <- c(counts, as.character(length(untested)))
        untested_ok <- untested_ok && length(untested) == length(without) && setequal(untested, without)
        tested <- table[!is.na(table$AveExpr), , drop = FALSE]
        ref <- tryCatch(s1Apply(binding, fit, contrast)$results[[contrast]], error = function(e) e)
        if (inherits(ref, "error")) {
            out$cmp_status <- s1Flat(paste("reduced contrast half:", conditionMessage(ref)), 200L)
            return(out)
        }
        index <- match(rownames(ref), as.character(tested[[id_column]]))
        ids_ok <- ids_ok && !anyNA(index) && nrow(tested) == nrow(ref)
        for (column in c("logFC", "t", "P.Value", "fdr_qvalue")) {
            head_values <- if (anyNA(index)) rep(NA_real_, nrow(ref)) else as.numeric(tested[[column]][index])
            ref_values <- as.numeric(ref[[column]])
            identical_all <- identical_all && identical(head_values, ref_values)
            both <- !is.na(head_values) & !is.na(ref_values)
            if (any(is.na(head_values) != is.na(ref_values))) rel_max <- Inf
            if (any(both)) {
                rel <- abs(head_values[both] - ref_values[both]) / pmax(abs(ref_values[both]), .Machine$double.xmin)
                rel_max <- max(rel_max, rel)
            }
        }
    }
    out$cmp_ids_ok <- ids_ok
    out$untested_ok <- untested_ok
    out$untested_counts <- paste(counts, collapse = ",")
    out$cmp_max_rel <- rel_max
    out$cmp_status <- if (identical_all && ids_ok) "identical" else if (ids_ok && !out$reduced_has_na && rel_max <= 1e-12) "rounding" else
        "differs"
    out
}

#' D6 (VAL-116 (7)): limma's robust step and zero-df features, descriptive. Reference: the contrast half on the journey's
#' DA input (TRUE, TRUE; the features without values removed by the fix). Alternative: the same halves on that input
#' without its zero-df features (df.residual 0 and a value). Each computation in its own tryCatch.
s1D6 <- function(binding, input, m, fit, without, contrasts, key) {
    zero <- tryCatch(rownames(fit$coefficients)[!is.na(fit$df.residual) & fit$df.residual == 0 & !is.na(fit$Amean)],
        error = function(e) e)
    if (inherits(zero, "error")) {
        return(list(s1Row(c(key, list(error = s1Flat(paste("zero-df:", conditionMessage(zero)), 200L))), S1_D6_COLUMNS)))
    }
    alt_fit <- if (length(zero)) tryCatch(s1Fit(binding, input, m[!rownames(m) %in% zero, , drop = FALSE]), error = function(e) e) else fit
    prior <- function(eb, drop) {
        df_prior <- as.numeric(eb$df.prior)
        keep <- if (length(df_prior) > 1L) !rownames(eb$coefficients) %in% drop else TRUE
        list(n = nrow(eb$coefficients), unequal = s1UnequalDf(eb$df.residual),
            down = if (length(df_prior) > 1L) sum(df_prior[keep] < max(df_prior)) else 0L,
            s2 = stats::median(as.numeric(eb$s2.prior)), min = min(df_prior), med = stats::median(df_prior), max = max(df_prior))
    }
    lapply(contrasts, function(contrast) {
        row <- c(key, list(contrast = contrast, error = "", n_features = nrow(m), n_without_values = length(without),
            n_zero_df = length(zero)))
        computed <- tryCatch({
            if (inherits(alt_fit, "error")) stop("alternative fit: ", conditionMessage(alt_fit))
            ref <- s1Apply(binding, fit, contrast)
            alt <- if (length(zero)) s1Apply(binding, alt_fit, contrast) else ref
            pr <- prior(ref$fit.eb, zero)
            pa <- prior(alt$fit.eb, character())
            tr <- ref$results[[contrast]]
            ta <- alt$results[[contrast]]
            both <- intersect(rownames(tr)[is.finite(tr$P.Value)], rownames(ta)[is.finite(ta$P.Value)])
            dp <- abs(tr[both, "P.Value"] - ta[both, "P.Value"])
            called_ref <- !is.na(tr[both, "fdr_qvalue"]) & tr[both, "fdr_qvalue"] < 0.05
            called_alt <- !is.na(ta[both, "fdr_qvalue"]) & ta[both, "fdr_qvalue"] < 0.05
            zero_rows <- intersect(zero, rownames(tr))
            list(ref_n = pr$n, alt_n = pa$n, ref_unequal_df = pr$unequal, alt_unequal_df = pa$unequal,
                ref_downweighted = pr$down, alt_downweighted = pa$down, ref_s2prior_median = pr$s2, alt_s2prior_median = pa$s2,
                ref_dfprior_min = pr$min, ref_dfprior_median = pr$med, ref_dfprior_max = pr$max, alt_dfprior_min = pa$min,
                alt_dfprior_median = pa$med, alt_dfprior_max = pa$max, n_both = length(both),
                p_absdiff_max = if (length(both)) max(dp) else NA_real_, p_absdiff_median = if (length(both)) stats::median(dp) else NA_real_,
                calls_ref = sum(called_ref), calls_alt = sum(called_alt), calls_only_ref = sum(called_ref & !called_alt),
                calls_only_alt = sum(called_alt & !called_ref),
                zero_df_called_ref = sum(!is.na(tr[zero_rows, "fdr_qvalue"]) & tr[zero_rows, "fdr_qvalue"] < 0.05))
        }, error = function(e) list(error = s1Flat(conditionMessage(e), 200L)))
        s1Row(utils::modifyList(row, computed), S1_D6_COLUMNS)
    })
}

# ---- parent ---------------------------------------------------------------------------------------------------------------
s1Prepare <- function(tree) {
    log <- tempfile("kit196-prepare-", fileext = ".txt")
    if (file.exists(file.path(tree, "src", "build_dia_rsum.R")) &&
        !file.exists(file.path(tree, "inst", "duckdb", "linux_amd64", "manifest.rds"))) {
        status <- system(sprintf("cd %s && %s --vanilla src/build_dia_rsum.R > %s 2>&1", shQuote(tree), shQuote(s1_rscript), shQuote(log)))
        if (!identical(as.integer(status), 0L)) return(sprintf("build_dia_rsum.R exit %s: %s", status, s1Flat(utils::tail(readLines(log), 5L))))
    }
    status <- system2(s1_rscript, c("-e", shQuote(sprintf("pkgbuild::compile_dll('%s', force = TRUE, quiet = TRUE)", tree))),
        stdout = log, stderr = log)
    if (!identical(as.integer(status), 0L)) return(sprintf("compile_dll exit %s: %s", status, s1Flat(utils::tail(readLines(log), 5L))))
    ""
}

s1RunSide <- function(part, side, tree, commit, out) {
    units <- do.call(rbind, lapply(names(S1_PARTS[[part]]$cells), function(cell) {
        do.call(rbind, lapply(S1_PARTS[[part]]$cells[[cell]], function(job) data.frame(cell = cell, job = job, rep = S1_REPS,
            stringsAsFactors = FALSE)))
    }))
    weight <- c(L12 = 1, L64 = 2.5, M12 = 0.5, M168 = 3)[units$cell] * c(deciding = 1, ba = 0.3)[units$job]
    units <- units[order(-weight), , drop = FALSE]
    logs <- file.path(tempdir(), sprintf("kit196-s1-%s-%s", part, side))
    dir.create(logs, showWarnings = FALSE)
    runChild <- function(i) {
        u <- units[i, ]
        tag <- sprintf("%s-%s-%d", u$cell, u$job, u$rep)
        rds <- file.path(logs, paste0(tag, ".rds"))
        log <- file.path(logs, paste0(tag, ".log"))
        t0 <- Sys.time()
        code <- system2(s1_rscript, c(shQuote(s1_file), "--unit", side, shQuote(tree), commit, u$cell, u$rep, u$job, shQuote(rds)),
            stdout = log, stderr = log)
        res <- if (file.exists(rds)) readRDS(rds) else list(error = sprintf("child exit %s, no result", code))
        res$exit <- code
        res$wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        if (!is.null(res$error) || !identical(as.integer(code), 0L)) res$log_tail <- s1Flat(utils::tail(readLines(log, warn = FALSE), 25L), 4000L)
        res
    }
    results <- parallel::mclapply(seq_len(nrow(units)), runChild, mc.cores = S1_WORKERS, mc.preschedule = FALSE)
    rows <- list()
    d6 <- list()
    timing <- list()
    failed <- character()
    for (i in seq_len(nrow(units))) {
        res <- results[[i]]
        if (!is.list(res)) res <- list(error = s1Flat(paste("worker failed:", paste(as.character(res), collapse = " "))))
        u <- units[i, ]
        error <- res$error %||% ""
        if (nzchar(error) || is.null(res$rows)) failed <- c(failed, sprintf("%s %s %d: %s | %s", u$cell, u$job, u$rep, error, res$log_tail %||% ""))
        if (!is.null(res$rows)) rows[[length(rows) + 1L]] <- res$rows
        if (!is.null(res$d6) && nrow(res$d6)) d6[[length(d6) + 1L]] <- res$d6
        timing[[i]] <- data.frame(side = side, commit = commit, cell = u$cell, job = u$job, rep = u$rep, exit = res$exit %||% NA,
            wall_seconds = round(res$wall %||% NA, 1), unit_seconds = round(res$seconds %||% NA, 1), ns_path = res$ns_path %||% "",
            error = s1Flat(error, 300L), stringsAsFactors = FALSE)
    }
    rows <- if (length(rows)) do.call(rbind, rows) else s1Bind(list(), S1_COLUMNS)
    d6 <- if (length(d6)) do.call(rbind, d6) else s1Bind(list(), S1_D6_COLUMNS)
    timing <- do.call(rbind, timing)
    s1Write(rows, file.path(out, sprintf("s1-%s-%s.tsv", part, side)))
    s1Write(timing, file.path(out, sprintf("s1-units-%s-%s.tsv", part, side)))
    list(rows = rows, d6 = d6, timing = timing, failed = failed)
}

s1Write <- function(table, path) utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
s1Journeys <- function(rows) rows[!duplicated(rows[S1_KEY]), , drop = FALSE]
s1KeyOf <- function(x) do.call(paste, c(lapply(x[S1_KEY], as.character), sep = "|"))

#' A2 (i): base's OR full stops, and n_T and miss_share of every OR journey, against KIT-188's rows (deciding and ba jobs).
s1OracleStops <- function(part, base_rows, out) {
    k188 <- do.call(rbind, lapply(names(S1_PARTS[[part]]$cells), function(cell) {
        path <- tempfile(fileext = ".tsv.gz")
        status <- system2("git", c("show", sprintf("%s:%s", S1_K188_COMMIT, sprintf(S1_K188_ROWS, S1_PARTS[[part]]$dir, cell))), stdout = path)
        if (!identical(as.integer(status), 0L) || !file.exists(path) || file.size(path) == 0) {
            stop(sprintf("KIT-188's rows of %s could not be read at %s (git show exit %s)", cell, S1_K188_COMMIT, status))
        }
        x <- utils::read.delim(gzfile(path), colClasses = "character", na.strings = "NA", quote = "", comment.char = "")
        x[x$arm == "OR" & x$job %in% S1_PARTS[[part]]$cells[[cell]], c("cell", "job", "scenario", "mask", "rep", "assay", "status",
            "class", "n_T", "miss_share"), drop = FALSE]
    }))
    k188 <- k188[!duplicated(k188[c("cell", "job", "scenario", "mask", "rep", "assay")]), , drop = FALSE]
    mine <- base_rows[base_rows$arm == "OR", , drop = FALSE]
    key <- function(x) paste(x$cell, x$job, x$scenario, x$mask, x$rep, x$assay, sep = "|")
    keys <- union(key(k188), key(mine))
    k <- k188[match(keys, key(k188)), , drop = FALSE]
    s <- mine[match(keys, key(mine)), , drop = FALSE]
    s1_miss <- ifelse(is.na(s$miss_share), NA_character_, as.character(signif(as.numeric(s$miss_share), 7L)))
    miss_equal <- !is.na(k$miss_share) & !is.na(s1_miss) & (k$miss_share == s1_miss |
        abs(as.numeric(k$miss_share) - as.numeric(s1_miss)) <= 1e-12 * pmax(abs(as.numeric(k$miss_share)), 1e-300))
    table <- data.frame(key = keys, k188_status = k$status, s1_status = s$status, k188_class = s1Flat_v(k$class),
        s1_class = s1Flat_v(s$class), stop_equal = !is.na(k$status) & !is.na(s$status) & ((k$status == "crashed") == (s$status == "crashed")),
        k188_n_T = k$n_T, s1_n_T = s$n_T, n_T_equal = !is.na(k$n_T) & !is.na(s$n_T) & as.character(k$n_T) == as.character(s$n_T),
        k188_miss_share = k$miss_share, s1_miss_share = s1_miss, miss_equal = miss_equal,
        s1_n_without_values = s$n_without_values, s1_err_contrasts = s$err_contrasts, stringsAsFactors = FALSE)
    s1Write(table, file.path(out, sprintf("oracle-stops-%s.tsv", part)))
    table
}
s1Flat_v <- function(x) vapply(x, function(v) if (is.na(v)) NA_character_ else s1Flat(v, 200L), character(1L), USE.NAMES = FALSE)

#' A2 (ii) to (iv), base against head, per journey (and per assay where the criterion is per assay).
s1Compare2 <- function(base, head) {
    rows <- list()
    add <- function(journey, field, b, h, equal, criterion) {
        rows[[length(rows) + 1L]] <<- data.frame(journey = journey, field = field, base = s1Flat(b, 300L), head = s1Flat(h, 300L),
            equal = isTRUE(equal), criterion = criterion, stringsAsFactors = FALSE)
    }
    bj <- s1Journeys(base)
    hj <- s1Journeys(head)
    for (k in union(s1KeyOf(bj), s1KeyOf(hj))) {
        b <- base[s1KeyOf(base) == k, , drop = FALSE]
        h <- head[s1KeyOf(head) == k, , drop = FALSE]
        if (!nrow(b) || !nrow(h)) {
            add(k, "present", nrow(b), nrow(h), FALSE, "journey on both sides")
            next
        }
        b_errors <- sum(as.integer(b$err_contrasts), na.rm = TRUE)
        h_errors <- sum(as.integer(h$err_contrasts), na.rm = TRUE)
        b_completes <- identical(b$status[[1L]], "ok") && b_errors == 0L
        b_failed <- isTRUE(as.logical(b$full_stop[[1L]])) || b_errors > 0L
        if (b_completes) {
            add(k, "da_sha", b$da_sha[[1L]], h$da_sha[[1L]], identical(b$da_sha[[1L]], h$da_sha[[1L]]), "A2 (ii): identical")
            add(k, "long_sha", b$long_sha[[1L]], h$long_sha[[1L]], identical(b$long_sha[[1L]], h$long_sha[[1L]]), "A2 (ii): identical")
        } else if (b_failed) {
            add(k, "completes", sprintf("%s; errors %d", b$status[[1L]], b_errors), sprintf("%s; errors %d", h$status[[1L]], h_errors),
                identical(h$status[[1L]], "ok") && h_errors == 0L, "A2 (iii): completes on head")
            for (i in seq_len(nrow(h))) {
                a <- h$assay[[i]]
                bi <- b[b$assay %in% a, , drop = FALSE]
                ok <- h$cmp_status[[i]] %in% c("identical", "rounding") &&
                    (identical(h$cmp_status[[i]], "identical") || identical(as.logical(h$reduced_has_na[[i]]), FALSE))
                add(paste(k, a), "tested_rows", sprintf("without values %s", if (nrow(bi)) bi$n_without_values[[1L]] else NA),
                    sprintf("%s; max relative %s; reduced input has a missing value %s", h$cmp_status[[i]], h$cmp_max_rel[[i]], h$reduced_has_na[[i]]),
                    ok, "A2 (iii): identical where the reduced input holds a missing value, else relative 1e-12")
                add(paste(k, a), "untested_rows", h$without_values[[i]], sprintf("%s (per contrast %s)", h$untested_ok[[i]], h$untested_counts[[i]]),
                    identical(as.logical(h$untested_ok[[i]]), TRUE), "A2 (iii): rows with AveExpr missing are the features without values")
            }
        } else {
            add(k, "status", b$status[[1L]], h$status[[1L]], identical(b$status[[1L]], h$status[[1L]]), "no DA reached on base: same status")
        }
        # (iv) head: no failure with either message
        h_msg <- paste(c(h$err_message, h$stop_message), collapse = " ")
        add(k, "head_messages", "", s1Flat(h_msg, 200L),
            !grepl("prior.weights contain NA values|covariate contains NA or infinite values", h_msg), "A2 (iv): no head failure with either message")
    }
    do.call(rbind, rows)
}

#' A2 (iv), base: every assay with an error has a feature without values, and the direct re-run gives the message of its
#' df case.
s1BaseDiagnosis <- function(base) {
    e <- base[!is.na(base$err_contrasts) & as.integer(base$err_contrasts) > 0L, , drop = FALSE]
    if (!nrow(e)) return(data.frame(journey = character(), assay = character(), n_without_values = integer(), unequal_df = logical(),
        direct_message = character(), expected = character(), holds = logical(), stringsAsFactors = FALSE))
    expected <- ifelse(as.logical(e$unequal_df), "prior.weights contain NA values", "covariate contains NA or infinite values")
    data.frame(journey = s1KeyOf(e), assay = e$assay, full_stop = e$full_stop, n_without_values = as.integer(e$n_without_values),
        unequal_df = as.logical(e$unequal_df), direct_message = e$direct_message, direct_fn = e$direct_fn, expected = expected,
        holds = as.integer(e$n_without_values) >= 1L & e$direct_message == expected, stringsAsFactors = FALSE)
}

s1D6Summary <- function(d6) {
    if (!nrow(d6)) return(d6)
    measures <- c("ref_downweighted", "alt_downweighted", "ref_s2prior_median", "alt_s2prior_median", "ref_dfprior_min",
        "ref_dfprior_median", "alt_dfprior_min", "alt_dfprior_median", "p_absdiff_max", "p_absdiff_median", "calls_ref", "calls_alt",
        "calls_only_ref", "calls_only_alt", "zero_df_called_ref", "n_zero_df")
    groups <- unique(d6[c("cell", "job", "mask")])
    do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
        g <- d6[d6$cell == groups$cell[[i]] & d6$job == groups$job[[i]] & d6$mask == groups$mask[[i]], , drop = FALSE]
        ok <- g[is.na(g$error) | !nzchar(g$error), , drop = FALSE]
        values <- unlist(lapply(measures, function(name) {
            x <- suppressWarnings(as.numeric(ok[[name]]))
            x <- x[!is.na(x)]
            stats::setNames(c(if (length(x)) stats::median(x) else NA, if (length(x)) max(x) else NA), paste0(name, c("_median", "_max")))
        }))
        cbind(groups[i, , drop = FALSE], data.frame(rows = nrow(g), errors = nrow(g) - nrow(ok), as.list(values), check.names = FALSE))
    }))
}

s1Main <- function(a) {
    flags <- a[startsWith(a, "--")]
    positional <- a[!startsWith(a, "--")]
    stopifnot(length(positional) >= 2L, positional[[1L]] %in% names(S1_PARTS), grepl("^[0-9a-f]{40}$", positional[[2L]]))
    part <- positional[[1L]]
    base_commit <- positional[[2L]]
    base_only <- "--base-only" %in% flags
    out <- if (length(positional) >= 3L) positional[[3L]] else "docs/qa/workflow-execution/kit/kit-196/s1"
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(out)
    head_root <- normalizePath(".")
    head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    started <- Sys.time()
    say <- function(...) {
        line <- paste0("KIT-196 S1 ", part, ": ", sprintf(...))
        cat(line, "\n", sep = "")
        cat(line, "\n", sep = "", file = file.path(out, sprintf("summary-%s.txt", part)), append = TRUE)
    }
    unlink(file.path(out, sprintf("summary-%s.txt", part)))
    say("head %s, base %s, base only %s, copy SHA-256 %s", head_commit, base_commit, base_only, s1Sha(file.path(s1_dir, S1_COPY)))
    records <- sub("^--base-records=", "", flags[startsWith(flags, "--base-records=")])
    base_tree <- file.path(tempdir(), "kit196-s1-base")
    if (!length(records)) {
        invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
        if (!dir.exists(file.path(base_tree, "R"))) stop("the base worktree could not be made")
        base_tree <- normalizePath(base_tree)
    }
    say("git diff --name-only 98df468b %s -- R/: %s", base_commit,
        paste(system2("git", c("diff", "--name-only", "98df468b", base_commit, "--", "R/"), stdout = TRUE), collapse = " "))
    sides <- list(base = list(tree = base_tree, commit = base_commit))
    if (!base_only) sides$head <- list(tree = head_root, commit = head_commit)
    results <- list()
    failed <- character()
    if (length(records)) {
        path <- file.path(records[[1L]], sprintf("s1-%s-base.tsv", part))
        base_rows <- utils::read.delim(path, quote = "", comment.char = "", na.strings = "NA", stringsAsFactors = FALSE)
        if (!all(base_rows$commit == base_commit)) stop("the base records were not made on ", base_commit)
        results$base <- list(rows = base_rows, failed = character(), timing = data.frame())
        sides$base <- NULL
        say("base side read from %s (SHA-256 %s): %d journeys", path, s1Sha(path), nrow(s1Journeys(base_rows)))
    }
    for (side in names(sides)) {
        t0 <- Sys.time()
        prepared <- s1Prepare(sides[[side]]$tree)
        if (nzchar(prepared)) {
            failed <- c(failed, paste(side, "preparation:", prepared))
            say("%s: preparation failed: %s", side, prepared)
            next
        }
        results[[side]] <- s1RunSide(part, side, sides[[side]]$tree, sides[[side]]$commit, out)
        failed <- c(failed, if (length(results[[side]]$failed)) paste(side, results[[side]]$failed))
        r <- results[[side]]$rows
        j <- s1Journeys(r)
        say("%s: %d journeys (%d units, %d failed to return) in %.1f min; full stops %d; journeys with an assay error %d",
            side, nrow(j), nrow(results[[side]]$timing), length(results[[side]]$failed),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(as.logical(j$full_stop), na.rm = TRUE),
            length(unique(s1KeyOf(r[!is.na(r$err_contrasts) & as.integer(r$err_contrasts) > 0L, , drop = FALSE]))))
    }
    for (f in failed) say("FAILED %s", s1Flat(f, 2000L))
    base <- results$base$rows
    if (!is.null(base) && nrow(base)) {
        oracle <- tryCatch(s1OracleStops(part, base, out), error = function(e) {
            failed <<- c(failed, paste("A2 (i):", conditionMessage(e)))
            say("A2 (i) could not be computed: %s", conditionMessage(e))
            NULL
        })
        if (!is.null(oracle)) {
            stops_k <- sort(unique(sub("\\|[^|]*$", "", oracle$key[oracle$k188_status %in% "crashed"])))
            stops_s <- sort(unique(sub("\\|[^|]*$", "", oracle$key[oracle$s1_status %in% "crashed"])))
            say("A2 (i): OR full stops KIT-188 %d, S1 base %d, list equal %s; n_T equal %d of %d; miss_share equal %d of %d; verdict %s",
                length(stops_k), length(stops_s), identical(stops_k, stops_s), sum(oracle$n_T_equal), nrow(oracle), sum(oracle$miss_equal),
                nrow(oracle), if (all(oracle$stop_equal) && all(oracle$n_T_equal) && all(oracle$miss_equal)) "HOLDS" else "FAILS")
            for (s in stops_s) say("A2 (i) S1 OR full stop: %s", s)
        }
        partial <- base[!as.logical(base$full_stop) & !is.na(base$err_contrasts) & as.integer(base$err_contrasts) > 0L, , drop = FALSE]
        s1Write(partial, file.path(out, sprintf("partial-failures-%s.tsv", part)))
        say("base partial failures (an assay's error in a journey that completed): %d journeys, %d assays",
            length(unique(s1KeyOf(partial))), nrow(partial))
        diagnosis <- s1BaseDiagnosis(base)
        s1Write(diagnosis, file.path(out, sprintf("base-diagnosis-%s.tsv", part)))
        say("A2 (iv) base: %d assays with an error, %d with a feature without values (counts %s), direct message as its df case %d; verdict %s",
            nrow(diagnosis), sum(diagnosis$n_without_values >= 1L), paste(sort(diagnosis$n_without_values), collapse = ","),
            sum(diagnosis$direct_message == diagnosis$expected), if (nrow(diagnosis) && all(diagnosis$holds)) "HOLDS" else if (!nrow(diagnosis)) "no errors" else "FAILS")
        say("direct messages on base: %s", paste(names(table(diagnosis$direct_message)), table(diagnosis$direct_message), sep = " x", collapse = "; "))
    }
    if (!base_only && !is.null(results$head) && !is.null(base)) {
        compare <- s1Compare2(base, results$head$rows)
        s1Write(compare, file.path(out, sprintf("compare-%s.tsv", part)))
        for (criterion in unique(compare$criterion)) {
            x <- compare[compare$criterion == criterion, , drop = FALSE]
            say("%s: %d rows, %d hold", criterion, nrow(x), sum(x$equal))
        }
        s1Write(results$head$d6, file.path(out, sprintf("d6-zero-df-%s.tsv", part)))
        s1Write(s1D6Summary(results$head$d6), file.path(out, sprintf("d6-summary-%s.tsv", part)))
        say("D6: %d rows, %d with an error (descriptive only)", nrow(results$head$d6),
            sum(!is.na(results$head$d6$error) & nzchar(results$head$d6$error)))
    }
    if (!length(records)) invisible(system2("git", c("worktree", "remove", "--force", shQuote(base_tree)), stdout = TRUE, stderr = TRUE))
    say("R %s, limma %s, BLAS %s; wall %.1f min; exit %d", R.version.string, as.character(utils::packageVersion("limma")),
        extSoftVersion()[["BLAS"]], as.numeric(difftime(Sys.time(), started, units = "mins")), as.integer(length(failed) > 0L))
    quit(save = "no", status = if (length(failed)) 1L else 0L)
}

if (length(s1_argv) >= 1L && identical(s1_argv[[1L]], "--unit")) {
    result <- tryCatch(s1Child(s1_argv), error = function(e) list(error = s1Flat(paste(class(e)[[1L]], conditionMessage(e)), 2000L)))
    saveRDS(result, s1_argv[[8L]])
    quit(save = "no", status = if (is.null(result$error)) 0L else 1L)
}
s1Main(s1_argv)
