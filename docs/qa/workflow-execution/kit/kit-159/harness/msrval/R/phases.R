# The phases that run cells (KIT-159 section 7.8, steps 5 to 9): the ticket's prepare.R, the cell
# table, the cell runner by arm kind (7.5) with the kernel identity gate, per-unit metrics, the
# run-level intervals and classification (4.4 to 4.7), blind scoring in the lifecycle, and the run
# folder (5.7). `val-run.R` wires these to the declaration check, retrieval and the evidence check.

#' Source the ticket's pinned R code into one isolated environment (7.8 step 5)
#'
#' @param code_files Paths of the declaration's `code` files; files under a `tests/` folder and
#'   files that are not `.R` are skipped.
#' @return The environment.
#' @export
valTicketEnv <- function(code_files) {
    env <- new.env(parent = globalenv())
    keep <- grepl("\\.[Rr]$", code_files) & !grepl("(^|/)tests/", code_files)
    for (f in code_files[keep]) sys.source(f, envir = env)
    env
}

#' Run the ticket's preparation: `prepareInputs`, `buildTruth`, `featureKey` (7.8 step 5)
#'
#' The key universe is built as 4.1 fixes it: every feature of the declared input role
#' (`prepareInputs()$features`, a data frame `assay, feature_id`), mapped by the ticket's
#' `featureKey(feature_id, arm)` with `arm = NA`; only features whose status is `ok` enter it, each
#' with its own assay. `universe.tsv` (`assay, key`) is written into `work_dir` (it is
#' feature-level and stays in the job); its SHA-256, its row count and the mapping's counts by
#' status and reason go into provenance (`inputs.universe`). Every prepared file is hashed.
#'
#' @param env From [valTicketEnv()].
#' @param data_dir The fetched inputs.
#' @param work_dir The preparation's private folder.
#' @return A list: `units` (per unit: `files`, `design`, `contrasts`, ...), `truth` (named list),
#'   `universe` (`assay`, `key`), `universe_counts` (`status`, `reason`, `n`), `universe_sha256`,
#'   `features`, `prepared` (path, sha256), `result_tables`.
#' @export
valPrepare <- function(env, data_dir, work_dir) {
    for (fn in c("prepareInputs", "buildTruth", "featureKey")) {
        if (!is.function(env[[fn]])) {
            valAbort("the ticket's code defines no {fn}()", "val_spec_invalid")
        }
    }
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    prep <- env$prepareInputs(data_dir, work_dir)
    if (!is.list(prep[["units"]]) || is.null(names(prep[["units"]]))) {
        valAbort("prepareInputs() must return a named list `units`", "val_spec_invalid")
    }
    feats <- prep[["features"]]
    if (!is.data.frame(feats) || !all(c("assay", "feature_id") %in% names(feats)) ||
            nrow(feats) == 0L) {
        valAbort(paste("prepareInputs() must return `features`, data.frame(assay, feature_id):",
                       "every feature of the declared input role (4.1)"), "val_spec_invalid")
    }
    truth <- env$buildTruth(data_dir)
    if (is.data.frame(truth)) truth <- list(primary = truth)
    for (id in names(truth)) valCheckTruth(truth[[id]], id)
    uni <- valUniverse(env, feats)
    universe_path <- file.path(work_dir, "universe.tsv")
    utils::write.table(uni$universe, universe_path, sep = "\t", quote = FALSE, row.names = FALSE)
    # R5: a unit with no files gives NULL here, and file.exists(NULL) is an error
    files <- as.character(unique(unlist(lapply(prep[["units"]], function(u) unlist(u[["files"]])))))
    files <- files[file.exists(files)]
    prepared <- data.frame(path = files, sha256 = vapply(files, valSha256, character(1)),
                           stringsAsFactors = FALSE, row.names = NULL)
    list(units = prep[["units"]], truth = truth, universe = uni$universe,
         universe_counts = uni$counts, universe_sha256 = valSha256(universe_path),
         features = feats, prepared = prepared,
         result_tables = prep[["result_tables"]] %||% list())
}

#' A truth table's required columns and unique keys
#' @noRd
valCheckTruth <- function(t1, id) {
    if (!is.data.frame(t1) || !all(c("key", "class", "expected", "truth_status") %in% names(t1))) {
        valAbort("truth table {id} lacks key, class, expected or truth_status", "val_truth_invalid")
    }
    if (anyDuplicated(t1$key)) {
        valAbort("truth table {id} holds a key more than once", "val_truth_invalid")
    }
    invisible(TRUE)
}

KEY_STATUSES <- c("ok", "unmapped", "ambiguous")

#' The ticket's `featureKey(feature_id, arm)` under its contract (4.1)
#'
#' @param env The ticket's code environment.
#' @param feature_id Feature identifiers.
#' @param arm The arm id (`NA` when the universe is built).
#' @return `data.frame(key, status, reason)`, one row per feature.
#' @noRd
valFeatureKeys <- function(env, feature_id, arm) {
    fk <- env$featureKey(as.character(feature_id), arm)
    if (!is.data.frame(fk) || !all(c("key", "status") %in% names(fk)) ||
            nrow(fk) != length(feature_id)) {
        valAbort(paste("featureKey(feature_id, arm) must return data.frame(key, status, reason),",
                       "one row per feature (4.1)"), "val_spec_invalid")
    }
    status <- as.character(fk$status)
    if (!all(status %in% KEY_STATUSES)) {
        valAbort("featureKey() gave a status outside ok, unmapped and ambiguous",
                 "val_spec_invalid")
    }
    reason <- if (is.null(fk[["reason"]])) rep(NA_character_, nrow(fk)) else
        as.character(fk$reason)
    data.frame(key = as.character(fk$key), status = status, reason = reason,
               stringsAsFactors = FALSE)
}

#' Counts of a key mapping by status and reason
#' @noRd
valKeyCounts <- function(fk) {
    if (nrow(fk) == 0L) {
        return(data.frame(status = character(), reason = character(), n = integer(),
                          stringsAsFactors = FALSE))
    }
    reason <- ifelse(is.na(fk$reason), "", fk$reason)
    tab <- stats::aggregate(list(n = rep(1L, nrow(fk))),
                            list(status = fk$status, reason = reason), sum)
    tab <- tab[order(tab$status, tab$reason, method = "radix"), , drop = FALSE]
    rownames(tab) <- NULL
    tab
}

#' The key universe of the declared input role's features (4.1)
#' @noRd
valUniverse <- function(env, features) {
    fk <- valFeatureKeys(env, features$feature_id, NA_character_)
    ok <- fk$status == "ok" & !is.na(fk$key) & nzchar(fk$key)
    u <- unique(data.frame(assay = as.character(features$assay)[ok], key = fk$key[ok],
                           stringsAsFactors = FALSE))
    u <- u[order(u$assay, u$key, method = "radix"), , drop = FALSE]
    rownames(u) <- NULL
    list(universe = u, counts = valKeyCounts(fk))
}

#' Rows of an arm's keyed table whose (assay, key) is `ok` and outside the universe
#' @noRd
valNotInUniverse <- function(tab, universe) {
    u_assay <- universe$assay %||% rep("all", nrow(universe))
    u <- paste(u_assay, universe$key, sep = "\r")
    sum(tab$key_status == "ok" & !is.na(tab$key) &
            !(paste(tab$assay, tab$key, sep = "\r") %in% u))
}

#' A declared data set and arm by id
#' @noRd
valDeclDataset <- function(decl, id) {
    ds <- Filter(function(d) identical(d[["id"]], id), decl[["datasets"]])
    if (length(ds) != 1L) valAbort("data set {id} is not declared once", "val_spec_invalid")
    ds[[1L]]
}

#' @noRd
valDeclArm <- function(decl, id) {
    a <- Filter(function(x) identical(x[["id"]], id), decl[["arms"]])
    if (length(a) != 1L) valAbort("arm {id} is not declared once", "val_spec_invalid")
    a[[1L]]
}

#' The cells of a phase: units x arms (7.6, 7.8)
#'
#' One `run` per (unit run, arm), so rows naming one run are executed once (A39). The seed is the
#' data set's seed plus the run's position among the unit table's runs, the same for every arm of a
#' unit (paired arms see the same random stream). In phase `lifecycle` the units are the declared
#' `lifecycle_units`, plus `not_under_test.arms` on `not_under_test.units` (C60).
#'
#' @param decl The declaration (a list).
#' @param dataset Data set id.
#' @param units The unit table (data frame or path).
#' @param phase The phase.
#' @param unit_ids Restrict to these units, or `NULL`.
#' @param arms Restrict to these arm ids, or `NULL` for all.
#' @return Data frame: `cell`, `unit`, `run`, `arm`, `kind`, `route`, `seed`.
#' @export
valPhaseCells <- function(decl, dataset, units, phase, unit_ids = NULL, arms = NULL) {
    ds <- valDeclDataset(decl, dataset)
    units <- valUnits(units)
    all_runs <- unique(units$run)
    arm_ids <- vapply(decl[["arms"]], function(a) a[["id"]], character(1))
    arms <- arms %||% arm_ids
    if (!all(arms %in% arm_ids)) valAbort("an arm is not declared", "val_spec_invalid")
    if (!is.null(unit_ids)) units <- units[units$unit %in% unit_ids, , drop = FALSE]
    grid <- expand.grid(i = seq_len(nrow(units)), a = seq_along(arms), KEEP.OUT.ATTRS = FALSE)
    u <- units[grid$i, , drop = FALSE]
    arm <- arms[grid$a]
    if (identical(phase, "lifecycle") && length(ds[["lifecycle_units"]]) > 0L) {
        life <- unlist(ds[["lifecycle_units"]])
        nut <- ds[["not_under_test"]] %||% list()
        keep <- u$unit %in% life |
            (u$unit %in% unlist(nut[["units"]]) & arm %in% unlist(nut[["arms"]]))
        u <- u[keep, , drop = FALSE]
        arm <- arm[keep]
    }
    kind <- vapply(arm, function(a) valDeclArm(decl, a)[["kind"]], character(1))
    route <- vapply(arm, function(a) valDeclArm(decl, a)[["route"]] %||% ds[["route"]] %||% "none",
                    character(1))
    seed <- as.integer(ds[["seed"]]) + match(u$run, all_runs) - 1L
    data.frame(cell = paste(u$unit, arm, sep = ":"), unit = u$unit,
               run = paste(u$run, arm, sep = "@"),
               arm = arm, kind = unname(kind), route = unname(route), seed = seed,
               stringsAsFactors = FALSE, row.names = NULL)
}

#' A cell's private work folder
#' @noRd
valCellDir <- function(work_dir, run) file.path(work_dir, gsub("[^A-Za-z0-9._-]", "_", run))

#' Write a DA table in the standard columns (7.3)
#' @noRd
valWriteDa <- function(da, path) {
    cols <- c("feature_id", "assay", "contrast", "log2fc", "p", "s_own", "call_statistic",
              "called_own", "mean_expression")
    if (is.null(da[["assay"]])) da$assay <- "all"
    for (col in setdiff(cols, names(da))) da[[col]] <- NA
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(da[, cols], path, sep = "\t", quote = FALSE, row.names = FALSE)
    path
}

#' Signal a refusal from a caught condition
#' @noRd
valRefuseFrom <- function(e, step) {
    valArmRefused(step, class(e)[[1L]],
                  digest::digest(valMaskMessage(conditionMessage(e)), algo = "sha256",
                                 serialize = FALSE))
}

#' Run one cell by its arm's kind (7.5; the chunk's `runCell`)
#'
#' `workflow`: [valRunWorkflow()] with the route's action script and the arm builder, standard
#' tables into the cell's folder, the post-cell hook; `kernel`, `custom`, `published`: the ticket's
#' function on the unit's prepared inputs; `baseline`: [valBaselineArm()] on the declared input. A
#' refusal is signalled as `val_arm_refused` (step and class recorded); an unexpected end of a
#' workflow cell is an error (`did_not_complete`).
#'
#' @param cell One row of [valPhaseCells()] as a list.
#' @param ctx List: `decl`, `dataset`, `units` (prepared), `env` (ticket code), `work_dir`,
#'   `pkg_root`, optional `helpers`, `record_messages`, `out_dir`.
#' @return A list: `da` (path of the cell's `da.tsv`), `tables`, `ledger`.
#' @export
valRunCell <- function(cell, ctx) {
    decl <- ctx$decl
    ds <- valDeclDataset(decl, ctx$dataset)
    arm <- valDeclArm(decl, cell$arm)
    base_run <- sub("@[^@]*$", "", cell$run)
    inputs <- ctx$units[[base_run]] %||% ctx$units[[cell$unit]]
    if (is.null(inputs)) {
        valAbort("prepareInputs() gave no inputs for {cell$unit}", "val_spec_invalid")
    }
    cell_dir <- valCellDir(ctx$work_dir, cell$run)
    da_path <- file.path(cell_dir, "da.tsv")
    contrasts <- unlist(inputs[["contrasts"]] %||% ds[["contrasts"]])
    kind <- arm[["kind"]]
    if (identical(kind, "baseline")) {
        tab <- utils::read.delim(inputs$files[[arm[["input_role"]]]], check.names = FALSE,
                                 stringsAsFactors = FALSE)
        id <- arm[["column"]] %||% names(tab)[[1L]]
        q <- as.matrix(tab[, setdiff(names(tab), id), drop = FALSE])
        rownames(q) <- as.character(tab[[id]])
        res <- tryCatch(valBaselineArm(q, inputs$design, contrasts,
                                       data_type = arm[["data_type"]] %||% "intensity",
                                       normalise = arm[["normalise"]] %||% "none",
                                       batch = arm[["batch"]]),
                        val_baseline_refused = function(e) valRefuseFrom(e, "baseline"))
        d <- res$da
        # The baseline's rows carry the assay its input belongs to (4.1, I02).
        valWriteDa(data.frame(feature_id = d$key, assay = arm[["assay"]] %||% "all",
                              contrast = d$contrast, log2fc = d$estimate,
                              p = d$p, s_own = d$s_own, call_statistic = d$call_statistic,
                              called_own = d$called_own, stringsAsFactors = FALSE), da_path)
        return(list(da = da_path))
    }
    if (kind %in% c("custom", "published")) {
        spec <- list(id = arm[["id"]], kind = kind, fn = ctx$env[[arm[["function"]]]],
                     anchor = if (!is.null(arm[["anchor"]])) ctx$env[[arm[["anchor"]]]],
                     anchor_on_fail = arm[["anchor_on_fail"]])
        res <- valRunCustomArm(spec, cell$unit, inputs)
        if (isTRUE(res$refused)) valArmRefused("arm", res$class, res$message_sha256)
        valWriteDa(res$da, da_path)
        return(list(da = da_path, label = res$label, primary = res$primary))
    }
    if (identical(kind, "kernel")) {
        fn <- ctx$env[[arm[["function"]]]]
        if (!is.function(fn)) valAbort("kernel arm {arm$id}: no function", "val_spec_invalid")
        res <- tryCatch(fn(cell$unit, inputs), error = function(e) valRefuseFrom(e, "kernel"))
        if (is.list(res) && !is.data.frame(res)) {
            if (isTRUE(res$refused)) {
                valArmRefused(res$step %||% "kernel", res$class %||% "unclassed",
                              res$message_sha256)
            }
            res <- res$da
        }
        valWriteDa(res, da_path)
        return(list(da = da_path))
    }
    if (!identical(kind, "workflow")) {
        valAbort("arm {arm$id}: unknown kind {kind}", "val_spec_invalid")
    }
    valRunWorkflowCell(cell, ctx, ds, arm, inputs, contrasts, cell_dir)
}

#' A workflow cell (7.3): action script, arm builder, forced settings, hook
#' @noRd
valRunWorkflowCell <- function(cell, ctx, ds, arm, inputs, contrasts, cell_dir) {
    owned <- c("multischolar.test_mode", "MULTISCHOLAR_TEST_MODE")
    if (any(vapply(arm[["forced"]] %||% list(), function(f) f[["setting"]] %in% owned,
                   logical(1)))) {
        valAbort("arm {arm$id}: test mode cannot be forced (the harness owns it)",
                 "val_spec_invalid")
    }
    builder <- get(arm[["builder"]], envir = asNamespace("msrval"))
    arm_obj <- do.call(builder, arm[intersect(names(arm), names(formals(builder)))])
    pairs <- inputs[["contrast_pairs"]] %||%
        lapply(contrasts, function(x) strsplit(x, "-", fixed = TRUE)[[1L]])
    script <- valActionScript(cell$route, assays = inputs[["assays"]],
                              import_tool = inputs[["import_tool"]],
                              import_values = inputs[["import_values"]] %||% list(),
                              contrasts = pairs, ruv = valArmRuvRun(arm_obj))
    forced <- arm[["forced"]] %||% list()
    is_option <- vapply(forced, function(f) {
        grepl("^(MultiScholaR|multischolar)\\.", f[["setting"]])
    }, logical(1))
    as_named <- function(fs) {
        stats::setNames(lapply(fs, function(f) f[["value"]]),
                        vapply(fs, function(f) f[["setting"]], character(1)))
    }
    forced_inputs <- c(as_named(forced[!is_option]), inputs[["forced_inputs"]] %||% list())
    attr(forced_inputs, "reason") <- paste(vapply(forced[!is_option], function(f) {
        f[["reason"]] %||% ""
    }, character(1)), collapse = "; ")
    hook_name <- ctx$decl[["hooks"]][["post_cell"]]
    summary <- isTRUE(ctx$decl[["hooks"]][["summary"]])
    on_complete <- if (!is.null(hook_name)) {
        function(sctx) {
            valCallHook(ctx$env[[hook_name]], sctx$manager,
                        shiny::isolate(sctx$workflow_data$da_analysis_results_list),
                        if (summary) sctx$paths[["integration_dir"]],
                        as.list(sctx$tables), unit = cell$unit,
                        arm = cell$arm, design_path = inputs[["design_path"]],
                        out_dir = ctx$out_dir)
        }
    }
    ann <- ds[["annotation"]]
    if (!is.null(ann)) {
        ann <- utils::modifyList(ann, list(accessions = inputs[["accessions"]]))
        if (!is.null(ann[["in_taxon"]])) ann$in_taxon <- ctx$env[[ann[["in_taxon"]]]]
    }
    res <- valRunWorkflow(script, inputs$files, inputs$design, arm = arm_obj,
                          helpers = ctx$helpers %||% valGuiHelpers(ctx$pkg_root),
                          forced = as_named(forced[is_option]), forced_inputs = forced_inputs,
                          seed = as.integer(cell$seed), uniprot_cache = inputs[["uniprot_cache"]],
                          annotation = ann,
                          tables_dir = cell_dir, work_dir = file.path(cell_dir, "work"),
                          roles = inputs[["roles"]] %||% list(), on_complete = on_complete,
                          pkg_root = ctx$pkg_root,
                          technical_replicates = identical(ds[["replicate_type"]], "technical"),
                          record_messages = ctx$record_messages %||% TRUE,
                          replay = valReplayDeclared(ctx$decl, cell$unit, cell$arm))
    if (identical(res$exit, "refused")) {
        s <- res$steps[res$steps$outcome == "refused", , drop = FALSE]
        s <- s[nrow(s), , drop = FALSE]
        valArmRefused(s$step %||% "unknown", s$class %||% "unclassed", s$message_sha256)
    }
    if (!identical(res$exit, "completed")) {
        stop(sprintf("did_not_complete (%s)", res$reason %||% res$exit), call. = FALSE)
    }
    list(da = res$tables[["da"]], tables = res$tables, ledger = res$ledger_path,
         fallback = res$fallback, options_at_start = res$options_at_start,
         qvalue_fallback = res$qvalue_fallback, hook = res$hook, replay = res$replay)
}

#' The kernel identity gate before a phase's kernel cells (7.5; A14, A30)
#'
#' For each kernel arm, [valKernelIdentity()] on units drawn by its `identity_seed` (plus its
#' declared `refusal_units`), against the arm named by `identity_of` run on the same unit with the
#' same seed. A failed gate refuses every cell of that kernel arm (`val_kernel_identity_failed`).
#'
#' @param cells From [valPhaseCells()].
#' @param ctx As [valRunCell()].
#' @return A list: `refused` (runs of failed kernel arms) and one gate result per kernel arm.
#' @export
valKernelGate <- function(cells, ctx) {
    out <- list(refused = character())
    for (a in unique(cells$arm[cells$kind == "kernel"])) {
        arm <- valDeclArm(ctx$decl, a)
        kc <- cells[cells$arm == a, , drop = FALSE]
        one <- function(arm_id, u) {
            cell <- as.list(kc[kc$unit == u, , drop = FALSE][1L, ])
            ref <- valDeclArm(ctx$decl, arm_id)
            cell$arm <- arm_id
            cell$kind <- ref[["kind"]]
            cell$run <- paste0(sub("@[^@]*$", "", cell$run), "@", arm_id, "@identity")
            r <- tryCatch(valRunCell(cell, ctx), val_arm_refused = function(e) {
                list(refused = TRUE, step = e$step, class = e$refusal_class,
                     message_sha256 = e$message_sha256)
            })
            if (isTRUE(r$refused)) return(r)
            list(refused = FALSE, da = utils::read.delim(r$da, stringsAsFactors = FALSE))
        }
        gate <- tryCatch(
            valKernelIdentity(unique(kc$unit), kernel = function(u) one(a, u),  # nolint
                              workflow = function(u) one(arm[["identity_of"]], u),
                              identity_seed = as.integer(arm[["identity_seed"]]),
                              refusal_units = unlist(arm[["refusal_units"]]) %||% character()),
            val_kernel_identity_failed = function(e) {  # nolint
                list(passed = FALSE, units = character(), note = conditionMessage(e))
            })
        if (!isTRUE(gate$passed)) out$refused <- c(out$refused, kc$run)
        out[[a]] <- gate
    }
    out
}

#' The view label of a metric row: `intersection` and `union` depend on the comparison's arms
#' @noRd
valViewKey <- function(view, arms) {
    if (identical(view, "own")) view else sprintf("%s(%s)", view, paste(arms, collapse = "+"))
}

#' A metric function: from msrval or the ticket's code
#' @noRd
valMetricFunction <- function(name, env) {
    fn <- if (!is.null(env) && is.function(env[[name]])) env[[name]] else
        get0(name, envir = asNamespace("msrval"), mode = "function")
    if (!is.function(fn)) valAbort("metric function {name} is not defined", "val_spec_invalid")
    fn
}

#' Per-unit metric rows of one unit (7.8 step 8; `metrics.tsv`)
#'
#' Each arm's DA table is keyed by the ticket's `featureKey`, scored against the truth over the
#' universe ([valScoreTable()]), restricted to each declared view of each comparison, and the
#' metric function applied per arm. With `blind_seed` the truth is blinded first (phase
#' `lifecycle`, A04). An arm with no DA table (refused, not completed) gives NA rows.
#'
#' @param decl The declaration.
#' @param dataset Data set id.
#' @param unit Unit id.
#' @param da Named list (by arm) of DA table paths or data frames.
#' @param truth Named list of truth tables.
#' @param universe The key universe.
#' @param env The ticket's code environment (`featureKey`, ticket metric functions).
#' @param blind_seed The declaration's blind seed in phase `lifecycle`, else `NULL`.
#' @return Rows in the columns of `metrics.tsv`.
#' @export
valUnitMetrics <- function(decl, dataset, unit, da, truth, universe, env, blind_seed = NULL) {
    arms_tab <- valArmTables(da, env)
    ds <- valDeclDataset(decl, dataset)
    expected <- ds[["expected_not_in_universe"]] %||% list()
    outside <- vapply(arms_tab, function(t1) valNotInUniverse(t1, universe), integer(1))
    for (a in names(outside)) {
        n_out <- outside[[a]]
        want <- expected[[a]]
        if (n_out > 0L && (is.null(want) || !identical(as.integer(want), n_out))) {
            valAbort(paste("arm {a} of unit {unit}: {n_out} rows outside the key universe; an",
                           "instrument stop (declare expected_not_in_universe to admit them)"),
                     "val_spec_invalid")
        }
    }
    rows <- list()
    for (m in decl[["metrics"]]) {
        tid <- m[["truth"]] %||% "primary"
        tr <- truth[[tid]]
        if (is.null(tr)) valAbort("metric {m$id}: no truth table {tid}", "val_spec_invalid")
        if (!is.null(blind_seed)) tr <- valBlindTruth(tr, as.integer(blind_seed))
        tr <- tr[, c("key", "class", "expected", "truth_status")]
        fn <- valMetricFunction(m[["function"]], env)
        mtype <- if (identical(m[["class"]], "accuracy")) "accuracy" else "detection"
        input <- m[["metric_input"]] %||% "scored"
        for (cmp in m[["comparisons"]]) {
            cmp_arms <- unlist(cmp)
            present <- cmp_arms[cmp_arms %in% names(arms_tab)]
            scored <- if (length(present) == length(cmp_arms)) {
                s <- valScoreTable(arms_tab[cmp_arms], tr, universe, unit = unit, truth_id = tid)
                s[s$contrast == m[["contrast"]], , drop = FALSE]
            }
            for (v in unlist(m[["views"]])) {
                sv <- if (!is.null(scored)) valViews(scored, cmp_arms, v, mtype)
                for (a in cmp_arms) {
                    value <- NA_real_
                    status <- "arm_not_completed"
                    if (!is.null(sv)) {
                        args <- m[["args"]] %||% list()
                        rows_a <- sv[sv$arm == a, , drop = FALSE]
                        first <- switch(input,
                                        scored = rows_a,
                                        confusion = valConfusion(rows_a, call_rule =
                                                                     m[["call_rule"]] %||% "own"),
                                        matrix = valCellMatrix(file.path(dirname(da[[a]]),
                                                                         "matrix.tsv")))
                        if (input == "scored" && "call_rule" %in% names(formals(fn))) {
                            args$call_rule <- m[["call_rule"]]
                        }
                        if (input == "matrix" && is.null(first)) {
                            res <- list(value = NA_real_, status = "no_matrix")
                        } else {
                            res <- do.call(fn, c(list(first), args))
                        }
                        value <- if (is.list(res)) as.numeric(res[["value"]]) else as.numeric(res)
                        status <- if (is.list(res)) res[["status"]] %||% "ok" else "ok"
                    }
                    rows[[length(rows) + 1L]] <- data.frame(
                        ticket = decl[["ticket"]] %||% NA_character_, dataset = dataset,  # nolint
                        truth_id = tid, unit = unit, contrast = m[["contrast"]], arm = a,
                        assay = "all", view = valViewKey(v, cmp_arms),
                        class = m[["args"]][["class"]] %||% NA_character_,
                        call_rule = m[["call_rule"]] %||% NA_character_, metric = m[["id"]],
                        statistic_class = m[["class"]] %||% NA_character_, value = value,
                        method = NA_character_, form = NA_character_, scope = NA_character_,
                        ci_low = NA_real_, ci_high = NA_real_, achieved_level = NA_real_,
                        n_units = NA_integer_, unit_kind = "unit", B = NA_integer_,
                        seed = NA_integer_, na_share = NA_real_, status = status,
                        primary = !is.null(m[["primary_family"]]),
                        family = m[["family"]] %||% NA_character_,
                        blinded = !is.null(blind_seed), stringsAsFactors = FALSE)
                }
            }
        }
    }
    out <- if (length(rows) == 0L) {
        as.data.frame(stats::setNames(replicate(length(EVIDENCE_COLUMNS[["metrics.tsv"]]),
                                                character(), simplify = FALSE),
                                      EVIDENCE_COLUMNS[["metrics.tsv"]]))
    } else {
        do.call(rbind, rows)
    }
    kc <- lapply(names(arms_tab), function(a) {
        k <- attr(arms_tab[[a]], "key_counts")
        if (!is.null(k) && nrow(k)) cbind(data.frame(unit = unit, arm = a), k)
    })
    attr(out, "key_counts") <- do.call(rbind, kc)
    attr(out, "not_in_universe") <- if (length(outside)) {
        data.frame(unit = unit, arm = names(outside), n = unname(outside))
    }
    out
}

#' An arm's matrix (`matrix.tsv` of its cell) as a numeric matrix, or `NULL`
#' @noRd
valCellMatrix <- function(path) {
    if (!is.character(path) || !file.exists(path)) return(NULL)
    x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    num <- vapply(x, is.numeric, logical(1)) & !names(x) %in% c("feature_id", "assay")
    out <- as.matrix(x[, num, drop = FALSE])
    rownames(out) <- paste(x$assay, x$feature_id, sep = ":")
    out
}

#' Run-level intervals and classification over the units (7.8 step 8; 4.4 to 4.7)
#'
#' Per metric, comparison and view: the per-unit differences (two arms) or values (one arm against
#' the declared `fixed`), then the declared method: `jackknife` over the `full` and `deletion`
#' units, or `paired_t`, `exact_median` or `clopper_pearson` over `independent` and `draw` units
#' with the refusal rule ([valRefusalFill()]) for one-arm and two-arm comparisons alike. A unit is
#' refused only when a cell of the comparison's arms on it is not `completed` (`cells`); a completed
#' unit whose metric is NA for a data reason is dropped and counted in `flags` (`data_na=<k>`), its
#' reason in `metrics.tsv`'s `status`. A jackknife whose `full` or `deletion` unit is refused is
#' `not_evaluable`. Statuses go into `flags` (`not_evaluable`, `refusal_substitute_required`,
#' `n_refused=<k>`); a refused classification is recorded in `regime` and never aborts the run.
#'
#' The widening (4.4, A38; O8 ruled in VAL-109): a jackknife interval (forms (ii), (iii)) is
#' widened to the `full` unit's within-run interval of the same view; over independent units each
#' unit's within-run interval is reported in `flags` beside the verdict, without widening.
#'
#' A metric with `family_adjustment` (`"bonferroni"`, or `list(method = "bonferroni", m = <m>)`)
#' is also classified at 1 - 0.05/m (m: the declared count, else the metric's comparisons) into
#' `outcome_family` (4.6), the only outcome that may feed M22. Compound members of each primary
#' family are classified after their components; their family is written in `flags`.
#'
#' @param decl The declaration.
#' @param metrics Per-unit rows (`unit`, `metric`, `arm`, `view`, `truth_id`, `value`).
#' @param units The unit table.
#' @param regimes The regime table for the classifier (default the harness's).
#' @param phase The phase.
#' @param within The run's within-run rows ([valWithinRunRows()]); a metric whose method is
#'   `bootstrap_within_run` is classified on them, one row per unit.
#' @param cells The cells' exits (`unit`, `arm`, `exit`), or `NULL` (then a non-finite value is
#'   taken as a refused unit).
#' @return Rows in the columns of `comparisons.tsv`, with attribute `verdicts`: per classified
#'   component and member its regime, measured coverage, Monte Carlo SE, statistic and family (I05).
#' @export
valRunLevel <- function(decl, metrics, units, regimes = NULL, phase = "counted", within = NULL,
                        cells = NULL) {
    units <- valUnits(units)
    blind <- identical(phase, "lifecycle")
    rows <- list()
    outcomes <- list()
    verdicts <- list()
    fam <- valFamilyOfComponents(decl)
    for (m in decl[["metrics"]]) {
        method <- (m[["interval"]] %||% list())[["method"]] %||% "paired_t"
        conf <- m[["conf"]] %||% 0.95
        for (cmp in m[["comparisons"]]) {
            cmp_arms <- unlist(cmp)
            for (v in unlist(m[["views"]])) {
                vk <- valViewKey(v, cmp_arms)
                component <- paste0(m[["id"]], ":", paste(cmp_arms, collapse = "-"))
                if (identical(method, "bootstrap_within_run")) {
                    got <- valClassifyWithin(m, cmp_arms, vk, within, regimes, blind)
                    rows <- c(rows, got$rows)
                    outcomes[[component]] <- got$outcome
                    next
                }
                one <- valLevelOne(decl, m, cmp_arms, v, vk, metrics, units, cells, regimes,
                                   blind, conf, within)
                cl <- one$cl
                flags <- one$flags
                outcome_family <- NA_character_
                adj <- m[["family_adjustment"]]
                if (!blind && !is.null(adj)) {
                    fm <- as.integer(if (is.list(adj)) adj[["m"]] else length(m[["comparisons"]]))
                    adj_method <- if (is.list(adj)) adj[["method"]] else adj
                    if (!identical(adj_method, "bonferroni") || length(fm) != 1L || fm < 1L) {
                        valAbort("family_adjustment must be bonferroni with m >= 1",
                                 "val_spec_invalid")
                    }
                    conf_adj <- 1 - 0.05 / fm
                    adj_one <- valLevelOne(decl, m, cmp_arms, v, vk, metrics, units, cells,
                                           regimes, blind, conf_adj, within = NULL)
                    outcome_family <- adj_one$cl$outcome
                    flags <- c(flags, sprintf("family_adjustment=bonferroni(m=%d;level=%s)", fm,
                                              format(round(conf_adj, 6))))
                    if (!is.na(adj_one$cl$regime %||% NA) &&
                            startsWith(adj_one$cl$regime, "refused")) {
                        flags <- c(flags, "family_level_refused")
                    }
                }
                outcomes[[component]] <- cl
                iv <- one$iv
                rows[[length(rows) + 1L]] <- data.frame(
                    arm = cmp_arms[[1L]],  # nolint
                    comparator = if (length(cmp_arms) == 2L) cmp_arms[[2L]] else "fixed",
                    truth_id = m[["truth"]] %||% "primary", metric = m[["id"]], view = vk,
                    call_rule = m[["call_rule"]] %||% NA_character_,
                    diff = iv[["estimate"]] %||% NA_real_, method = method,
                    scope = iv[["scope"]] %||% NA_character_, ci_low = iv[["ci_low"]],
                    ci_high = iv[["ci_high"]], margin = m[["margin"]] %||% NA_real_,
                    conf = conf, outcome = cl$outcome, no_worse = cl$no_worse,
                    outcome_family = outcome_family, n_refused = as.integer(one$n_ref),
                    regime = cl$regime %||% NA_character_,
                    flags = paste(unique(flags), collapse = ","),
                    member = NA_character_, component = component, not_sensitive = NA,
                    amended = NA_character_, headline = isTRUE(m[["headline"]]), blinded = blind,
                    stringsAsFactors = FALSE)
                note <- sub("^within_run_widening=", "",
                            grep("^within_run_widening=", flags, value = TRUE))
                verdicts[[length(verdicts) + 1L]] <- valVerdictRow(
                    component, NA_character_, m[["id"]], cmp_arms, vk, cl, conf,
                    fam$family[component] %||% NA_character_, decl,
                    within_note = if (length(note)) note[[1L]] else NA_character_
                )
            }
        }
    }
    if (!blind) {
        for (fam_name in names(decl[["primary_family"]])) {
            members <- decl[["primary_family"]][[fam_name]][["members"]]
            for (nm in names(members)) {
                mem <- members[[nm]]
                mem$all_of <- unlist(mem$all_of)
                res <- valClassifyMember(mem, outcomes)
                flags <- c(res$guard_text %||% character(), paste0("family=", fam_name),
                           if (length(res$missing)) {
                               c("not_evaluable",
                                 paste0("no_verdict=", paste(res$missing, collapse = "+")))
                           })
                rows[[length(rows) + 1L]] <- data.frame(
                    arm = nm, comparator = NA_character_, truth_id = NA_character_,  # nolint
                    metric = NA_character_, view = NA_character_, call_rule = NA_character_,
                    diff = NA_real_, method = NA_character_, scope = NA_character_,
                    ci_low = NA_real_, ci_high = NA_real_, margin = NA_real_, conf = NA_real_,
                    outcome = res$outcome, no_worse = res$no_worse,
                    outcome_family = NA_character_, n_refused = NA_integer_,
                    regime = NA_character_, flags = paste(flags, collapse = "; "),
                    member = nm, component = NA_character_, not_sensitive = NA,
                    amended = NA_character_, headline = FALSE, blinded = FALSE,
                    stringsAsFactors = FALSE)
                verdicts[[length(verdicts) + 1L]] <- valVerdictRow(
                    NA_character_, nm, NA_character_, character(), NA_character_,
                    list(outcome = res$outcome), NA_real_, fam_name, decl
                )
            }
        }
    }
    out <- do.call(rbind, rows)
    out <- out[, EVIDENCE_COLUMNS[["comparisons.tsv"]]]
    attr(out, "verdicts") <- if (length(verdicts)) do.call(rbind, verdicts)
    out
}

#' Components of each primary family's members: `family` (component -> family) and `m`
#' @noRd
valFamilyOfComponents <- function(decl) {
    family <- character()
    m <- integer()
    for (f in names(decl[["primary_family"]])) {
        pf <- decl[["primary_family"]][[f]]
        comps <- unlist(lapply(pf[["members"]], function(mem) {
            c(unlist(mem[["all_of"]]),
              vapply(mem[["guards"]] %||% list(), function(g) g[["component"]], character(1)))
        }))
        family[comps] <- f
        m[[f]] <- as.integer(pf[["m"]] %||% length(pf[["members"]]))
    }
    list(family = family, m = m)
}

#' One row of the verdict side table (I05)
#' @noRd
valVerdictRow <- function(component, member, metric, cmp_arms, view, cl, conf, family, decl,
                          within_note = NA_character_) {
    fm <- if (is.na(family)) NA_integer_ else
        as.integer(decl[["primary_family"]][[family]][["m"]] %||%
                       length(decl[["primary_family"]][[family]][["members"]]))
    data.frame(component = component, member = member, metric = metric,
               arm = if (length(cmp_arms)) cmp_arms[[1L]] else member,
               comparator = if (length(cmp_arms) == 2L) cmp_arms[[2L]] else
                   if (length(cmp_arms) == 1L) "fixed" else NA_character_,
               view = view, outcome = cl$outcome %||% NA_character_, conf = conf,
               regime = cl$regime %||% NA_character_,
               regime_status = cl$regime_status %||% NA_character_,
               regime_status_r113 = cl$regime_status_r113 %||% NA_character_,
               coverage = as.numeric(cl$coverage %||% NA_real_),
               mc_se = as.numeric(cl$mc_se %||% NA_real_),
               nsim = as.numeric(cl$nsim %||% NA_real_),
               statistic = as.character(cl$statistic %||% NA_character_),
               coverage_min = as.numeric(cl$coverage_min %||% NA_real_),
               miss_twice_nominal = isTRUE(cl$miss_twice_nominal),
               family = unname(family), m = fm, within_note = within_note,
               stringsAsFactors = FALSE)
}

#' Units refused for a comparison: a cell of one of its arms on the unit is not `completed`
#' @noRd
valRefusedUnits <- function(cells, cmp_arms, unit_ids, d) {
    if (is.null(cells)) return(stats::setNames(!is.finite(d), unit_ids))
    out <- vapply(unit_ids, function(u) {
        ex <- cells$exit[cells$unit == u & cells$arm %in% cmp_arms]
        length(ex) < length(cmp_arms) || any(ex != "completed")
    }, logical(1))
    stats::setNames(out, unit_ids)
}

#' An interval over per-unit values by a run-level method
#' @noRd
valLevelInterval <- function(values, method, conf, unit_role, class) {
    if (identical(method, "clopper_pearson")) {
        if (length(values) == 0L) {
            return(valNewInterval(NA_real_, NA_real_, method, 0L, conf, status = "too_few_units"))
        }
        if (!all(values %in% c(0, 1))) {
            valAbort("clopper_pearson needs per-unit values 0 or 1", "val_spec_invalid")
        }
        return(valInterval(as.integer(sum(values)), method = "clopper_pearson",
                           n = length(values), conf = conf, unit_role = unit_role, class = class))
    }
    valInterval(values, method = method, conf = conf, unit_role = unit_role, class = class)
}

#' The interval, flags and classification of one comparison and view at one level
#' @noRd
valLevelOne <- function(decl, m, cmp_arms, v, vk, metrics, units, cells, regimes, blind, conf,
                        within) {
    iv_spec <- m[["interval"]] %||% list()
    method <- iv_spec[["method"]] %||% "paired_t"
    hib <- isTRUE(m[["higher_is_better"]])
    cls <- m[["class"]] %||% NA_character_
    value_of <- function(a) {
        r <- metrics[metrics$metric == m[["id"]] & metrics$arm == a &
                         metrics$view %in% c(vk, v), , drop = FALSE]
        stats::setNames(as.numeric(r$value), r$unit)[units$unit]
    }
    d <- value_of(cmp_arms[[1L]])
    if (length(cmp_arms) == 2L) d <- d - value_of(cmp_arms[[2L]])
    names(d) <- units$unit
    refused_u <- valRefusedUnits(cells, cmp_arms, units$unit, d)
    data_na <- !refused_u & !is.finite(d)
    flags <- if (any(data_na)) sprintf("data_na=%d", sum(data_na)) else character()
    rf <- NULL
    withheld <- NULL
    classify <- function(iv) {
        iv[["statistic"]] <- iv_spec[["statistic"]]
        tryCatch(valClassify(iv, m[["margin"]], hib, fixed = m[["fixed"]], view = v,
                             feature_interval_decides = isTRUE(m[["feature_interval_decides"]]),
                             regimes = regimes, accept_flagged = iv_spec[["accept_flagged"]]),
                 val_spec_invalid = function(e) {
                     list(outcome = NA_character_, no_worse = NA,
                          regime = paste("refused:", conditionMessage(e)))
                 })
    }
    if (identical(method, "jackknife")) {
        full <- units$unit[units$unit_role == "full"]
        del <- units$unit[units$unit_role == "deletion"]
        n_ref <- sum(refused_u[c(full, del)])
        iv <- valJackknife(full = unname(d[full][1L]), deletions = unname(d[del]), conf = conf,
                           form = iv_spec[["form"]] %||% "unit", critical = iv_spec[["critical"]],
                           transform = iv_spec[["transform"]], margin = m[["margin"]],
                           class = cls,
                           unequal_deletion = any(units$unequal_deletion %in% c("TRUE", TRUE)),
                           design = iv_spec[["design"]] %||% NA_character_)
        if (n_ref > 0L) {
            # 4.7: a refused full or deletion run leaves the jackknife not evaluable; the status
            # `undefined_deletion` is kept for a statistic undefined on completed data.
            iv$status <- "not_evaluable"
            iv$ci_low <- NA_real_
            iv$ci_high <- NA_real_
            flags <- c(flags, "not_evaluable", sprintf("n_refused=%d", n_ref))
        }
    } else {
        if (identical(method, "clopper_pearson") && length(cmp_arms) != 1L) {
            valAbort("clopper_pearson is a one-arm share over units (against `fixed`)",
                     "val_spec_invalid")
        }
        role <- units$unit_role
        ind <- units$unit[role %in% c("independent", "draw")]
        ind_roles <- role[role %in% c("independent", "draw")]
        unit_role <- if (length(ind_roles) && all(ind_roles == "draw")) "draw" else "independent"
        keep <- ind[!data_na[ind]]
        dd <- unname(d[keep])
        dd[refused_u[keep]] <- NA_real_
        n_ref <- sum(refused_u[keep])
        if (!blind && n_ref > 0L) {
            tol <- as.integer(decl[["refusal_tolerance"]][[m[["id"]]]] %||% 0L)
            rf <- valRefusalFill(
                dd, tolerance = tol,
                higher_is_better = hib, margin = m[["margin"]], method = method, conf = conf,
                regimes = regimes, accept_flagged = iv_spec[["accept_flagged"]],
                fixed = m[["fixed"]], unit_role = unit_role, class = cls, view = v
            )
            iv <- rf$interval_dropped %||%
                valNewInterval(NA_real_, NA_real_, method, length(dd) - n_ref, conf,
                               status = rf$status)
            flags <- c(flags, rf$flags)
        } else {
            iv <- valLevelInterval(dd[is.finite(dd)], method, conf, unit_role, cls)
        }
        if (!is.null(within) && length(cmp_arms) == 2L && is.null(rf)) {
            wr <- valWithinRowsOf(within, m[["id"]], cmp_arms, vk)
            if (nrow(wr) > 0L) {
                # O8 (VAL-109): over independent units nothing is widened at any unit count; each
                # unit's within-run interval is reported beside the verdict. Below 6 units this is
                # the owner's ruling; at 6 or more the register leaves the rule open, so the
                # harness reports only and says so (`widen_to_within_run` has no effect here).
                count <- if (length(keep) < 6L) "below_6" else "6_or_more_open"
                reported <- sprintf("within_run_reported=%s:[%s,%s]",
                                    sub("^.*unit=([^;]+);.*$", "\\1", wr$flags),
                                    signif(as.numeric(wr$ci_low), 6),
                                    signif(as.numeric(wr$ci_high), 6))
                flags <- c(flags, sprintf("within_run_widening=report_only(O8:%s)", count),
                           reported)
            }
        }
    }
    if (method == "jackknife" && !is.null(within) && length(cmp_arms) == 2L &&
            identical(iv$status, "ok")) {
        full_u <- units$unit[units$unit_role == "full"][1L]
        wr <- valWithinRowsOf(within, m[["id"]], cmp_arms, vk)
        if (nrow(wr) > 0L) {
            iv <- valWidenFromRequests(iv, wr, cmp_arms[[1L]], cmp_arms[[2L]], m[["id"]],
                                       widen = m[["widen_to_within_run"]], view = vk,
                                       unit = full_u)
            flags <- c(flags, iv[["flags"]][grepl("^within_run_reported", iv[["flags"]])])
        }
    }
    flags <- c(flags, iv[["flags"]] %||% character())
    cl <- list(outcome = "blinded", no_worse = NA, regime = NA_character_)
    if (!blind) {
        cl <- if (!is.null(withheld)) {
            list(outcome = NA_character_, no_worse = NA, regime = withheld)
        } else if (!is.null(rf)) {
            r <- rf$classification %||% list(regime = NA_character_)
            r$outcome <- rf$outcome
            r$no_worse <- rf$no_worse
            r
        } else {
            classify(iv)
        }
    }
    list(iv = iv, n_ref = n_ref, flags = flags, cl = cl)
}

#' A comparison's within-run rows (scope `within_run`) of one view
#' @noRd
valWithinRowsOf <- function(within, metric, cmp_arms, vk) {
    within[within$scope == "within_run" & within$metric == metric &
               within$arm == cmp_arms[[1L]] & within$comparator == cmp_arms[[2L]] &
               within$view == vk, , drop = FALSE]
}

#' Within-run rows of one unit: a paired percentile bootstrap over keys (4.4; contract row 15)
#'
#' For each metric that declares `within_run` (`B`, `seed`, optional `strata` = `"class"` and
#' `cluster` = a function of the keys in the ticket's code) or whose interval method is
#' `bootstrap_within_run`, and each two-arm comparison and view: one set of resampled keys per
#' replicate applied to both arms ([valBootstrap()]), the statistic the metric's difference. Rows
#' are in `comparisons.tsv`'s columns with scope `within_run`; `flags` carries `unit=<id>;n=<keys>`.
#'
#' @inheritParams valUnitMetrics
#' @return Rows in the columns of `comparisons.tsv` (none when no metric asks for them).
#' @export
valWithinRunRows <- function(decl, unit, da, truth, universe, env, blind_seed = NULL) {
    rows <- list()
    arms_tab <- valArmTables(da, env)
    for (m in decl[["metrics"]]) {
        spec <- m[["within_run"]]
        if (is.null(spec) && !identical(m[["interval"]][["method"]], "bootstrap_within_run")) next
        spec <- spec %||% list()
        tid <- m[["truth"]] %||% "primary"
        tr <- truth[[tid]]
        if (!is.null(blind_seed)) tr <- valBlindTruth(tr, as.integer(blind_seed))
        tr <- tr[, c("key", "class", "expected", "truth_status")]
        fn <- valMetricFunction(m[["function"]], env)
        mtype <- if (identical(m[["class"]], "accuracy")) "accuracy" else "detection"
        conf <- m[["conf"]] %||% 0.95
        for (cmp in m[["comparisons"]]) {
            cmp_arms <- unlist(cmp)
            if (length(cmp_arms) != 2L || !all(cmp_arms %in% names(arms_tab))) next
            s <- valScoreTable(arms_tab[cmp_arms], tr, universe, unit = unit, truth_id = tid)
            s <- s[s$contrast == m[["contrast"]], , drop = FALSE]
            for (v in unlist(m[["views"]])) {
                sv <- valViews(s, cmp_arms, v, mtype)
                a1 <- sv[sv$arm == cmp_arms[[1L]], , drop = FALSE]
                a2 <- sv[sv$arm == cmp_arms[[2L]], , drop = FALSE]
                keys <- intersect(paste(a1$assay, a1$key), paste(a2$assay, a2$key))
                i1 <- match(keys, paste(a1$assay, a1$key))
                i2 <- match(keys, paste(a2$assay, a2$key))
                args <- m[["args"]] %||% list()
                if ("call_rule" %in% names(formals(fn))) args$call_rule <- m[["call_rule"]]
                val <- function(rows_) {
                    r <- do.call(fn, c(list(rows_), args))
                    if (is.list(r)) as.numeric(r[["value"]]) else as.numeric(r)
                }
                data <- data.frame(j = seq_along(keys), class = a1$class[i1],
                                   stringsAsFactors = FALSE)
                cl_name <- spec[["cluster"]]
                if (!is.null(cl_name) && is.function(env[[cl_name]])) {
                    data$cluster <- env[[spec[["cluster"]]]](a1$key[i1])
                }
                stat <- function(df) {
                    val(a1[i1[df$j], , drop = FALSE]) - val(a2[i2[df$j], , drop = FALSE])
                }
                if (is.null(spec[["seed"]])) {
                    valAbort("metric {m$id}: within_run declares no seed", "val_spec_invalid")
                }
                if (!is.null(spec[["events"]])) {
                    data$events <- as.logical(a1[[spec[["events"]]]][i1])
                }
                iv <- valBootstrap(data, stat, B = as.integer(spec[["B"]] %||%
                                                                  valBootstrapMinB(conf)),
                                   conf = conf, seed = as.integer(spec[["seed"]]),
                                   events = if (!is.null(spec[["events"]])) "events",
                                   min_events = as.integer(spec[["min_events"]] %||% 5L),
                                   strata = if (identical(spec[["strata"]], "class")) "class",
                                   cluster = if (!is.null(data$cluster)) "cluster",
                                   class = m[["class"]] %||% NA_character_)
                rows[[length(rows) + 1L]] <- data.frame(
                    arm = cmp_arms[[1L]], comparator = cmp_arms[[2L]], truth_id = tid,  # nolint
                    metric = m[["id"]], view = valViewKey(v, cmp_arms),
                    call_rule = m[["call_rule"]] %||% NA_character_,
                    diff = iv[["estimate"]] %||% NA_real_, method = "bootstrap_within_run",
                    scope = "within_run", ci_low = iv[["ci_low"]], ci_high = iv[["ci_high"]],
                    margin = m[["margin"]] %||% NA_real_, conf = conf, outcome = NA_character_,
                    no_worse = NA, outcome_family = NA_character_, n_refused = 0L,
                    regime = NA_character_,
                    flags = sprintf("unit=%s;n=%d;status=%s", unit, nrow(data), iv[["status"]]),
                    member = NA_character_,
                    component = paste0(m[["id"]], ":", paste(cmp_arms, collapse = "-")),
                    not_sensitive = NA, amended = NA_character_, headline = FALSE,
                    blinded = !is.null(blind_seed), stringsAsFactors = FALSE)
            }
        }
    }
    if (length(rows) == 0L) return(NULL)
    do.call(rbind, rows)
}

#' Arm DA tables keyed for scoring by the ticket's `featureKey(feature_id, arm)` (4.1)
#'
#' Each arm's table carries the key and its status from the ticket's function (a feature
#' `unmapped` or `ambiguous` in one arm is excluded for all the comparison's arms by
#' [valViews()]); attribute `key_counts` holds the arm's counts by status and reason.
#' @noRd
valArmTables <- function(da, env) {
    Map(function(p, arm) {
        d <- if (is.character(p)) {
            utils::read.delim(p, stringsAsFactors = FALSE, colClasses = c(feature_id = "character"))
        } else {
            p
        }
        fid <- as.character(d$feature_id)
        fk <- valFeatureKeys(env, fid, arm)
        out <- data.frame(feature_id = fid, key = fk$key, key_status = fk$status,
                          assay = as.character(d$assay %||% rep("all", nrow(d))),
                          contrast = as.character(d$contrast), estimate = as.numeric(d$log2fc),
                          p = as.numeric(d$p), s_own = as.numeric(d$s_own),
                          stringsAsFactors = FALSE)
        attr(out, "key_counts") <- valKeyCounts(fk)
        out
    }, da, names(da))
}

#' Classify the within-run rows of a `bootstrap_within_run` metric, one per unit
#'
#' The component's outcome across units is the units' common outcome; where the units disagree
#' it is NA (TODO(KIT-159 audit O8, I16): how a within-run component's outcome is formed across
#' units is the owner's rule; until then no verdict is formed from disagreeing units).
#' @noRd
valClassifyWithin <- function(m, cmp_arms, vk, within, regimes, blind) {
    r <- if (is.null(within)) within else
        within[within$metric == m[["id"]] & within$arm == cmp_arms[[1L]] &
                   within$comparator == (if (length(cmp_arms) == 2L) cmp_arms[[2L]] else "") &
                   within$view == vk, , drop = FALSE]
    if (is.null(r) || nrow(r) == 0L) {
        valAbort("metric {m$id}: no within-run rows for a bootstrap_within_run verdict",
                 "val_spec_invalid")
    }
    out <- list()
    all_cl <- list()
    for (i in seq_len(nrow(r))) {
        x <- r[i, ]
        n <- as.integer(sub("^.*n=([0-9]+).*$", "\\1", x$flags))
        st <- sub("^.*status=([a-z_]+).*$", "\\1", x$flags)
        iv <- valNewInterval(as.numeric(x$ci_low), as.numeric(x$ci_high), "bootstrap_within_run",
                             n, as.numeric(x$conf), "within_run", m[["class"]] %||% NA_character_,
                             status = st, estimate = as.numeric(x$diff))
        cl <- list(outcome = "blinded", no_worse = NA, regime = NA_character_)
        if (!blind) {
            cl <- tryCatch(valClassify(iv, m[["margin"]], isTRUE(m[["higher_is_better"]]),
                                       fixed = m[["fixed"]], view = "intersection",
                                       feature_interval_decides = isTRUE(
                                           m[["feature_interval_decides"]]),  # nolint
                                       regimes = regimes,
                                       accept_flagged = m[["interval"]][["accept_flagged"]]),
                           val_spec_invalid = function(e) {  # nolint
                               list(outcome = NA_character_, no_worse = NA,
                                    regime = paste("refused:", conditionMessage(e)))
                           })  # nolint
        }
        x$outcome <- cl$outcome
        x$no_worse <- cl$no_worse
        x$regime <- cl$regime %||% NA_character_
        x$blinded <- blind
        out[[length(out) + 1L]] <- x
        all_cl[[length(all_cl) + 1L]] <- cl
    }
    oc <- unique(vapply(all_cl, function(cl) cl$outcome %||% NA_character_, character(1)))
    nw <- unique(vapply(all_cl, function(cl) isTRUE(cl$no_worse), logical(1)))
    outcome <- if (length(oc) == 1L) {
        list(outcome = oc, no_worse = if (length(nw) == 1L) nw else NA,
             regime = all_cl[[1L]]$regime)
    } else {
        list(outcome = NA_character_, no_worse = NA,
             regime = "withheld: the units' within-run outcomes differ; rule pending (O8)")
    }
    list(rows = out, outcome = outcome)
}

#' Run a phase's cells and score them into `runs/<label>/` (7.8 steps 5 to 9)
#'
#' The caller has checked the declaration and fetched the inputs (steps 1 to 4). Here: the
#' declaration's statistical schema, the ticket's preparation, the kernel gate, the cells in
#' concurrent chunks (each chunk with its backstop: cells x `cell_minutes` + 10), blind or unblinded
#' scoring, the run folder of 5.7 (`cells.tsv`, `calls.tsv`, `metrics.tsv`, `comparisons.tsv`,
#' `coverage.tsv`, one ledger per workflow cell, `tables.sha256`, `not-committed.sha256`,
#' `scored.tsv.gz` where permitted, `results.md` outside the lifecycle, `attribution.tsv`,
#' `allowed-paths.txt`), provenance with the fields of 5.6, the evidence check over everything
#' the job wrote under `val-out/` and `outputs.sha256` last.
#'
#' @param phase `lifecycle`, `counted`, `correction`, `determinism` or `amended`.
#' @param decl The declaration (a list).
#' @param dataset Data set id.
#' @param units The unit table.
#' @param label The request's label (the run folder's name).
#' @param code_files The declaration's `code` files.
#' @param data_dir The fetched inputs.
#' @param work_dir The job's private folder.
#' @param out_dir `val-out/`.
#' @param chunks Concurrent chunks.
#' @param regimes Regime table for the classifier (default the harness's; with ticket regimes from
#'   [valClassifierRegimes()]).
#' @param arms,unit_ids Restrict the cells.
#' @param manifest The input manifest (data frame), for `attribution.tsv` and the evidence check.
#' @param code The provenance `code` block.
#' @param pkg_root The package checkout.
#' @param ticket_coverage The committed `ticket-coverage.tsv` rows (provenance), or `NULL`.
#' @param declaration_diff The counted declaration's differences from the lifecycle one, or `NULL`.
#' @return The run folder, invisibly.
#' @export
valRunPhase <- function(phase, decl, dataset, units, label, code_files, data_dir, work_dir,
                        out_dir, chunks = 4L, regimes = NULL, arms = NULL, unit_ids = NULL,
                        manifest = NULL, code = NULL, pkg_root = ".", ticket_coverage = NULL,
                        declaration_diff = NULL) {
    started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    units <- valUnits(units)
    valCheckDeclarationSchema(decl, units)
    env <- valTicketEnv(code_files)
    prep <- valPrepare(env, data_dir, file.path(work_dir, "prepare"))
    run_dir <- file.path(out_dir, "runs", label)
    dir.create(file.path(run_dir, "ticket"), recursive = TRUE, showWarnings = FALSE)
    cells <- valPhaseCells(decl, dataset, units, phase, unit_ids, arms)
    universe_ids <- unique(c(prep$universe$key, as.character(prep$features$feature_id)))
    if (length(prep$result_tables)) {
        valCopyResultTables(prep$result_tables, decl[["outputs"]], run_dir,
                            arms = unique(cells$arm), units = unique(cells$unit),
                            universe = universe_ids)
    }
    no_text <- !is.null(manifest) && valPrivateLogs(manifest)
    # The hook writes only under runs/<label>/ticket/ (5.2, I08).
    ctx <- list(decl = decl, dataset = dataset, units = prep$units, code_files = code_files,
                work_dir = file.path(work_dir, "cells"), pkg_root = normalizePath(pkg_root),
                record_messages = !no_text,
                out_dir = normalizePath(file.path(run_dir, "ticket")))
    gate <- valKernelGate(cells, c(ctx, list(env = env)))
    ctx_path <- file.path(work_dir, "cell-context.rds")
    saveRDS(ctx, ctx_path)
    runner <- file.path(work_dir, "cell-runner.R")
    writeLines(c(sprintf("ctx <- readRDS(%s)", deparse(normalizePath(ctx_path))),
                 "ctx$env <- msrval::valTicketEnv(ctx$code_files)",
                 "runCell <- function(cell) msrval::valRunCell(cell, ctx)"), runner)
    to_run <- cells[!cells$run %in% gate$refused, , drop = FALSE]
    minutes <- as.numeric(decl[["ceilings"]][["cell_minutes"]] %||% 120)
    chunk_root <- file.path(work_dir, "chunks")
    ran <- if (nrow(to_run) > 0L) {
        valRunCells(to_run, runner, k = chunks, cell_minutes = minutes, work_dir = chunk_root,
                    chunk_minutes = valChunkMinutes(length(unique(to_run$run)), chunks, minutes))
    }
    res <- valCellResults(cells, ran)
    values <- valCellValues(res, chunk_root)
    res <- valCellSettings(res, values)
    cells_out <- res
    cells_out$run <- sub("@[^@]*$", "", cells_out$run)
    utils::write.table(cells_out[, EVIDENCE_COLUMNS[["cells.tsv"]]],
                       file.path(run_dir, "cells.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
    blind_seed <- if (identical(phase, "lifecycle")) decl[["blind"]][["seed"]] else NULL
    done <- res[res$exit == "completed", , drop = FALSE]
    da_of <- function(u) {
        r <- done[done$unit == u, , drop = FALSE]
        stats::setNames(as.list(file.path(valCellDir(ctx$work_dir, r$run), "da.tsv")), r$arm)
    }
    truths <- lapply(stats::setNames(unique(cells$unit), unique(cells$unit)), function(u) {
        valUnitTruth(decl, dataset, env, prep, u)
    })
    key_counts <- list()
    outside <- list()
    metrics <- do.call(rbind, lapply(unique(cells$unit), function(u) {
        mt <- valUnitMetrics(decl, dataset, u, da_of(u), truths[[u]], prep$universe, env,
                             blind_seed = blind_seed)
        key_counts[[u]] <<- attr(mt, "key_counts")
        outside[[u]] <<- attr(mt, "not_in_universe")
        mt
    }))
    utils::write.table(metrics, file.path(run_dir, "metrics.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "NA")
    within <- do.call(rbind, lapply(unique(cells$unit), function(u) {
        valWithinRunRows(decl, u, da_of(u), truths[[u]], prep$universe, env,
                         blind_seed = blind_seed)
    }))
    scored_units <- units[units$unit %in% cells$unit, , drop = FALSE]
    comparisons <- valRunLevel(decl, metrics, scored_units, regimes = regimes, phase = phase,
                               within = within, cells = res[, c("unit", "arm", "exit")])
    verdicts <- attr(comparisons, "verdicts")
    if (!is.null(within)) {
        key <- function(x) paste(x$component, x$scope, x$flags)
        comparisons <- rbind(comparisons, within[!key(within) %in% key(comparisons), ])
    }
    utils::write.table(comparisons, file.path(run_dir, "comparisons.tsv"), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = "NA")
    primary_units <- valPrimaryUnits(decl, dataset, unique(cells$unit))
    valWriteCalls(res, values, ctx$work_dir, run_dir)
    valWriteCoverageTable(decl, done, ctx$work_dir, truths, prep$universe, env, run_dir,
                          blind_seed, primary_units)
    valCopyLedgers(res, cells, ctx$work_dir, run_dir)
    private <- !is.null(manifest) &&
        !valPublishable(manifest, manifest$role[manifest$kind == "data"])
    withheld <- valWriteScored(done, ctx$work_dir, truths, prep$universe, env, run_dir,
                               lifecycle = identical(phase, "lifecycle"), private = private)
    tables <- valWriteTableHashes(res, values, run_dir, withheld,
                                  if (private) prep$prepared else NULL)
    attribution <- if (!is.null(manifest)) {
        manifest[, c("role", "kind", "licence", "redistribution", "citation")]
    } else {
        data.frame(role = character(), kind = character(), licence = character(),
                   redistribution = character(), citation = character())
    }
    if (!is.null(manifest)) {
        utils::write.table(attribution, file.path(run_dir, "attribution.tsv"), sep = "\t",
                           quote = FALSE, row.names = FALSE)
    }
    if (phase %in% c("counted", "amended", "correction") && !is.null(verdicts)) {
        valWriteResults(file.path(run_dir, "results.md"), verdicts, attribution)
    }
    allowed <- unlist(lapply(decl[["outputs"]] %||% list(), function(o) {
        valExpandOutputPaths(o, arms = unique(cells$arm), units = unique(cells$unit))
    }))
    writeLines(sort(unique(as.character(allowed))), file.path(run_dir, "allowed-paths.txt"))
    code <- code %||% valNoCode(decl[["stage"]] %||% phase)
    code$code_files <- code$code_files %||% lapply(code_files, function(f) {
        list(path = f, sha256 = valSha256(f))
    })
    if (!length(code$helpers_extracted)) {
        code$helpers_extracted <- tryCatch(valHelpersExtracted(cells, decl, ctx$pkg_root),
                                           error = function(e) list())
    }
    gates <- gate[setdiff(names(gate), "refused")]
    chunk_tab <- if (!is.null(ran)) attr(ran, "chunks")
    if (!is.null(chunk_tab)) {
        chunk_tab$peak_rss_mib <- vapply(chunk_tab$chunk, function(ch) {
            v <- res$vmhwm_mib[res$chunk %in% ch]
            if (any(is.finite(v))) max(v, na.rm = TRUE) else NA_real_
        }, numeric(1))
    }
    prov <- valProvenance(
        ticket = decl[["ticket"]] %||% "unknown", stage = decl[["stage"]] %||% phase,
        label = label, phase = phase, status = "completed", blinded = !is.null(blind_seed),
        started = started, code = code,
        inputs = list(prepared = prep$prepared,
                      universe = list(sha256 = prep$universe_sha256,
                                      rows = nrow(prep$universe),
                                      by_status = prep$universe_counts)),
        truth = c(valTruthProvenance(prep$truth, work_dir, blind_seed),
                  list(per_unit = valPerUnitTruthProvenance(truths, prep$truth, work_dir))),
        cells = valCellsProvenance(res, cells, values, run_dir, ctx$work_dir,
                                   !is.null(decl[["hooks"]][["post_cell"]]), decl = decl),
        chunks = chunk_tab %||% list(),
        scoring = list(
            metrics_sha256 = valSha256(file.path(run_dir, "metrics.tsv")),
            comparisons_sha256 = valSha256(file.path(run_dir, "comparisons.tsv")),
            intervals = verdicts %||% list(), ticket_coverage = ticket_coverage %||% list(),
            declaration_diff = declaration_diff %||% list(),
            key_status = do.call(rbind, key_counts) %||% list(),
            not_in_universe = do.call(rbind, outside) %||% list(),
            tables = tables,
            cell_counts = list(n = nrow(cells), completed = sum(res$exit == "completed"),
                               refused = sum(res$exit == "refused")),
            kernel_gates = lapply(gates, function(g) {
                list(passed = isTRUE(g$passed), units = g$units, note = g$note)
            })
        ),
        outputs = lapply(sort(list.files(run_dir, recursive = TRUE)), function(f) {
            list(path = f, sha256 = valSha256(file.path(run_dir, f)))
        })
    )
    valWriteProvenance(prov, file.path(run_dir, "provenance.json"))
    valCheckProvenance(file.path(run_dir, "provenance.json"))
    valCheckOutsideRuns(out_dir)
    if (!is.null(manifest)) {
        classes <- unique(unlist(lapply(prep$truth, function(t1) as.character(t1$class))))
        valCheckEvidence(run_dir, manifest, outputs = decl[["outputs"]] %||% list(),
                         arms = unique(cells$arm), units = unique(cells$unit),
                         universe = universe_ids, classes = classes,
                         row_bounds = list(metrics.tsv = valMetricsRowBound(decl, cells)),
                         phase = phase, primary_units = primary_units,
                         forbidden_sha256 = if (private) prep$prepared$sha256 else character())
    }
    valWriteOutputsSha(run_dir)
    invisible(run_dir)
}
