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
#' Writes `universe.tsv` (the union of the truth tables' keys, assay `all`) into `work_dir` (it is
#' feature-level and stays in the job) and hashes every prepared file.
#'
#' @param env From [valTicketEnv()].
#' @param data_dir The fetched inputs.
#' @param work_dir The preparation's private folder.
#' @return A list: `units` (per unit: `files`, `design`, `contrasts`, ...), `truth` (named list),
#'   `universe`, `prepared` (path, sha256), `result_tables`.
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
    truth <- env$buildTruth(data_dir)
    if (is.data.frame(truth)) truth <- list(primary = truth)
    for (id in names(truth)) {
        t1 <- truth[[id]]
        if (!all(c("key", "class", "expected", "truth_status") %in% names(t1))) {
            valAbort("truth table {id} lacks key, class, expected or truth_status",
                     "val_truth_invalid")
        }
        if (anyDuplicated(t1$key)) {
            valAbort("truth table {id} holds a key more than once", "val_truth_invalid")
        }
    }
    keys <- unique(unlist(lapply(truth, function(t1) as.character(t1$key))))
    universe <- data.frame(assay = "all", key = sort(keys, method = "radix"),
                           stringsAsFactors = FALSE)
    utils::write.table(universe, file.path(work_dir, "universe.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    files <- unique(unlist(lapply(prep[["units"]], function(u) unlist(u[["files"]]))))
    files <- files[file.exists(files)]
    prepared <- data.frame(path = files, sha256 = vapply(files, valSha256, character(1)),
                           stringsAsFactors = FALSE, row.names = NULL)
    list(units = prep[["units"]], truth = truth, universe = universe, prepared = prepared,
         result_tables = prep[["result_tables"]] %||% list())
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
        valWriteDa(data.frame(feature_id = d$key, contrast = d$contrast, log2fc = d$estimate,
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
                              contrasts = pairs)
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
                          record_messages = ctx$record_messages %||% TRUE)
    if (identical(res$exit, "refused")) {
        s <- res$steps[res$steps$outcome == "refused", , drop = FALSE]
        s <- s[nrow(s), , drop = FALSE]
        valArmRefused(s$step %||% "unknown", s$class %||% "unclassed", s$message_sha256)
    }
    if (!identical(res$exit, "completed")) {
        stop(sprintf("did_not_complete (%s)", res$reason %||% res$exit), call. = FALSE)
    }
    list(da = res$tables[["da"]], tables = res$tables, ledger = res$ledger_path,
         fallback = res$fallback, options_at_start = res$options_at_start)
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
            valKernelIdentity(unique(kc$unit), kernel = function(u) one(a, u),
                              workflow = function(u) one(arm[["identity_of"]], u),
                              identity_seed = as.integer(arm[["identity_seed"]]),
                              refusal_units = unlist(arm[["refusal_units"]]) %||% character()),
            val_kernel_identity_failed = function(e) {
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
    rows <- list()
    for (m in decl[["metrics"]]) {
        tid <- m[["truth"]] %||% "primary"
        tr <- truth[[tid]]
        if (!is.null(blind_seed)) tr <- valBlindTruth(tr, as.integer(blind_seed))
        tr <- tr[, c("key", "class", "expected", "truth_status")]
        fn <- valMetricFunction(m[["function"]], env)
        mtype <- if (identical(m[["class"]], "accuracy")) "accuracy" else "detection"
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
                        if ("call_rule" %in% names(formals(fn))) args$call_rule <- m[["call_rule"]]
                        res <- do.call(fn, c(list(sv[sv$arm == a, , drop = FALSE]), args))
                        value <- if (is.list(res)) as.numeric(res[["value"]]) else as.numeric(res)
                        status <- if (is.list(res)) res[["status"]] %||% "ok" else "ok"
                    }
                    rows[[length(rows) + 1L]] <- data.frame(
                        ticket = decl[["ticket"]] %||% NA_character_, dataset = dataset,
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
    if (length(rows) == 0L) {
        return(as.data.frame(stats::setNames(replicate(length(EVIDENCE_COLUMNS[["metrics.tsv"]]),
                                                       character(), simplify = FALSE),
                                             EVIDENCE_COLUMNS[["metrics.tsv"]])))
    }
    do.call(rbind, rows)
}

#' Run-level intervals and classification over the units (7.8 step 8; 4.4 to 4.7)
#'
#' Per metric, comparison and view: the per-unit differences (two arms) or values (one arm against
#' the declared `fixed`), then the declared method: `jackknife` over the `full` and `deletion`
#' units, or `paired_t` / `exact_median` over `independent` and `draw` units with the refusal rule
#' ([valRefusalFill()]) when units are missing. In phase `lifecycle` the outcome is `blinded`.
#' Compound members of each primary family are classified after their components.
#'
#' @param decl The declaration.
#' @param metrics Per-unit rows (`unit`, `metric`, `arm`, `view`, `truth_id`, `value`).
#' @param units The unit table.
#' @param regimes The regime table for the classifier (default the harness's).
#' @param phase The phase.
#' @param within The run's within-run rows ([valWithinRunRows()]): the run-level interval is widened
#'   to them bound by bound (4.4, A38; mandatory below 6 units), and a metric whose method is
#'   `bootstrap_within_run` is classified on them, one row per unit.
#' @return Rows in the columns of `comparisons.tsv`.
#' @export
valRunLevel <- function(decl, metrics, units, regimes = NULL, phase = "counted", within = NULL) {
    units <- valUnits(units)
    blind <- identical(phase, "lifecycle")
    rows <- list()
    outcomes <- list()
    for (m in decl[["metrics"]]) {
        iv_spec <- m[["interval"]] %||% list()
        method <- iv_spec[["method"]] %||% "paired_t"
        hib <- isTRUE(m[["higher_is_better"]])
        for (cmp in m[["comparisons"]]) {
            cmp_arms <- unlist(cmp)
            for (v in unlist(m[["views"]])) {
                vk <- valViewKey(v, cmp_arms)
                value_of <- function(a) {
                    r <- metrics[metrics$metric == m[["id"]] & metrics$arm == a &
                                     metrics$view %in% c(vk, v), , drop = FALSE]
                    stats::setNames(as.numeric(r$value), r$unit)[units$unit]
                }
                component <- paste0(m[["id"]], ":", paste(cmp_arms, collapse = "-"))
                if (identical(method, "bootstrap_within_run")) {
                    got <- valClassifyWithin(m, cmp_arms, vk, within, regimes, blind)
                    rows <- c(rows, got$rows)
                    outcomes[[component]] <- got$outcome
                    next
                }
                d <- value_of(cmp_arms[[1L]])
                if (length(cmp_arms) == 2L) d <- d - value_of(cmp_arms[[2L]])
                names(d) <- units$unit
                rf <- NULL
                if (identical(method, "jackknife")) {
                    full <- units$unit[units$unit_role == "full"]
                    del <- units$unit[units$unit_role == "deletion"]
                    n_ref <- sum(!is.finite(d[c(full, del)]))
                    iv <- valJackknife(full = unname(d[full][1L]), deletions = unname(d[del]),
                                       conf = m[["conf"]] %||% 0.95,
                                       form = iv_spec[["form"]] %||% "unit",
                                       critical = iv_spec[["critical"]],
                                       transform = iv_spec[["transform"]], margin = m[["margin"]],
                                       class = m[["class"]] %||% NA_character_,
                                       unequal_deletion = any(units$unequal_deletion %in%
                                                                  c("TRUE", TRUE)),
                                       design = iv_spec[["design"]] %||% NA_character_)
                } else if (identical(method, "clopper_pearson")) {
                    if (length(cmp_arms) != 1L) {
                        valAbort("clopper_pearson is a one-arm share over units (against `fixed`)",
                                 "val_spec_invalid")
                    }
                    ind <- units$unit[units$unit_role %in% c("independent", "draw")]
                    dd <- unname(d[ind])
                    n_ref <- sum(!is.finite(dd))
                    ok <- dd[is.finite(dd)]
                    if (!all(ok %in% c(0, 1))) {
                        valAbort("clopper_pearson needs per-unit values 0 or 1", "val_spec_invalid")
                    }
                    draw <- all(units$unit_role[units$unit %in% ind] == "draw")
                    iv <- valInterval(as.integer(sum(ok)), method = "clopper_pearson",
                                      n = length(ok), conf = m[["conf"]] %||% 0.95,
                                      unit_role = if (draw) "draw" else "independent",
                                      class = m[["class"]] %||% NA_character_)
                } else {
                    role <- units$unit_role
                    ind <- units$unit[role %in% c("independent", "draw")]
                    dd <- unname(d[ind])
                    n_ref <- sum(!is.finite(dd))
                    unit_role <- if (all(role[role %in% c("independent", "draw")] == "draw")) {
                        "draw"
                    } else {
                        "independent"
                    }
                    if (!blind && n_ref > 0L && length(cmp_arms) == 2L) {
                        rf <- valRefusalFill(dd, tolerance = as.integer(
                            decl[["refusal_tolerance"]][[m[["id"]]]] %||% 0L),
                            higher_is_better = hib, margin = m[["margin"]], method = method,
                            conf = m[["conf"]] %||% 0.95, regimes = regimes,
                            accept_flagged = iv_spec[["accept_flagged"]])
                        iv <- rf$interval_dropped %||%
                            valNewInterval(NA_real_, NA_real_, method, length(dd) - n_ref,
                                           m[["conf"]] %||% 0.95, status = rf$status)
                    } else {
                        iv <- valInterval(dd[is.finite(dd)], method = method,
                                          conf = m[["conf"]] %||% 0.95, unit_role = unit_role,
                                          class = m[["class"]] %||% NA_character_)
                    }
                }
                if (length(cmp_arms) == 2L && !is.null(within) && is.null(rf)) {
                    iv <- valWidenFromRequests(iv, within, cmp_arms[[1L]], cmp_arms[[2L]],
                                               m[["id"]], widen = m[["widen_to_within_run"]])
                }
                cl <- list(outcome = "blinded", no_worse = NA, regime = NA_character_)
                if (!blind) {
                    cl <- if (!is.null(rf)) {
                        list(outcome = rf$outcome %||% rf$status, no_worse = rf$no_worse,
                             regime = NA_character_)
                    } else {
                        tryCatch(valClassify(iv, m[["margin"]], hib, fixed = m[["fixed"]],
                                             view = v,
                                             feature_interval_decides = isTRUE(
                                                 m[["feature_interval_decides"]]),
                                             regimes = regimes,
                                             accept_flagged = iv_spec[["accept_flagged"]]),
                                 val_spec_invalid = function(e) {
                                     list(outcome = NA_character_, no_worse = NA,
                                          regime = paste("refused:", conditionMessage(e)))
                                 })
                    }
                }
                outcomes[[component]] <- cl
                rows[[length(rows) + 1L]] <- data.frame(
                    arm = cmp_arms[[1L]],
                    comparator = if (length(cmp_arms) == 2L) cmp_arms[[2L]] else "fixed",
                    truth_id = m[["truth"]] %||% "primary", metric = m[["id"]], view = vk,
                    call_rule = m[["call_rule"]] %||% NA_character_,
                    diff = iv[["estimate"]] %||% NA_real_, method = method,
                    scope = iv[["scope"]] %||% NA_character_, ci_low = iv[["ci_low"]],
                    ci_high = iv[["ci_high"]], margin = m[["margin"]] %||% NA_real_,
                    conf = m[["conf"]] %||% 0.95, outcome = cl$outcome, no_worse = cl$no_worse,
                    outcome_family = NA_character_, n_refused = as.integer(n_ref),
                    regime = cl$regime %||% NA_character_,
                    flags = paste(iv[["flags"]] %||% character(), collapse = ","),
                    member = NA_character_, component = component, not_sensitive = NA,
                    amended = NA_character_, headline = isTRUE(m[["headline"]]), blinded = blind,
                    stringsAsFactors = FALSE)
            }
        }
    }
    if (!blind) {
        for (fam in names(decl[["primary_family"]])) {
            members <- decl[["primary_family"]][[fam]][["members"]]
            for (nm in names(members)) {
                mem <- members[[nm]]
                mem$all_of <- unlist(mem$all_of)
                res <- valClassifyMember(mem, outcomes)
                rows[[length(rows) + 1L]] <- data.frame(
                    arm = nm, comparator = NA_character_, truth_id = NA_character_,
                    metric = NA_character_, view = NA_character_, call_rule = NA_character_,
                    diff = NA_real_, method = NA_character_, scope = NA_character_,
                    ci_low = NA_real_, ci_high = NA_real_, margin = NA_real_, conf = NA_real_,
                    outcome = res$outcome, no_worse = res$no_worse, outcome_family = fam,
                    n_refused = NA_integer_, regime = NA_character_,
                    flags = paste(res$guard_text %||% character(), collapse = "; "),
                    member = nm, component = NA_character_, not_sensitive = NA,
                    amended = NA_character_, headline = FALSE, blinded = FALSE,
                    stringsAsFactors = FALSE)
            }
        }
    }
    out <- do.call(rbind, rows)
    out[, EVIDENCE_COLUMNS[["comparisons.tsv"]]]
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
                if (is.function(env[[spec[["cluster"]] %||% ""]])) {
                    data$cluster <- env[[spec[["cluster"]]]](a1$key[i1])
                }
                stat <- function(df) val(a1[i1[df$j], , drop = FALSE]) -
                    val(a2[i2[df$j], , drop = FALSE])
                iv <- valBootstrap(data, stat, B = as.integer(spec[["B"]] %||%
                                                                  valBootstrapMinB(conf)),
                                   conf = conf, seed = as.integer(spec[["seed"]] %||% 1L),
                                   strata = if (identical(spec[["strata"]], "class")) "class",
                                   cluster = if (!is.null(data$cluster)) "cluster",
                                   class = m[["class"]] %||% NA_character_)
                rows[[length(rows) + 1L]] <- data.frame(
                    arm = cmp_arms[[1L]], comparator = cmp_arms[[2L]], truth_id = tid,
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

#' Arm DA tables keyed for scoring
#' @noRd
valArmTables <- function(da, env) {
    lapply(da, function(p) {
        d <- if (is.character(p)) utils::read.delim(p, stringsAsFactors = FALSE) else p
        data.frame(feature_id = as.character(d$feature_id),
                   key = as.character(env$featureKey(as.character(d$feature_id))),
                   assay = as.character(d$assay %||% "all"), contrast = as.character(d$contrast),
                   estimate = as.numeric(d$log2fc), p = as.numeric(d$p),
                   s_own = as.numeric(d$s_own), stringsAsFactors = FALSE)
    })
}

#' Classify the within-run rows of a `bootstrap_within_run` metric, one per unit
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
    last <- NULL
    for (i in seq_len(nrow(r))) {
        x <- r[i, ]
        n <- as.integer(sub("^.*n=([0-9]+).*$", "\\1", x$flags))
        st <- sub("^.*status=([a-z_]+).*$", "\\1", x$flags)
        iv <- valNewInterval(x$ci_low, x$ci_high, "bootstrap_within_run", n, x$conf,
                             "within_run", m[["class"]] %||% NA_character_, status = st,
                             estimate = x$diff)
        cl <- list(outcome = "blinded", no_worse = NA, regime = NA_character_)
        if (!blind) {
            cl <- tryCatch(valClassify(iv, m[["margin"]], isTRUE(m[["higher_is_better"]]),
                                       fixed = m[["fixed"]], view = "intersection",
                                       feature_interval_decides = isTRUE(
                                           m[["feature_interval_decides"]]),
                                       regimes = regimes,
                                       accept_flagged = m[["interval"]][["accept_flagged"]]),
                           val_spec_invalid = function(e) {
                               list(outcome = NA_character_, no_worse = NA,
                                    regime = paste("refused:", conditionMessage(e)))
                           })
        }
        x$outcome <- cl$outcome
        x$no_worse <- cl$no_worse
        x$regime <- cl$regime %||% NA_character_
        x$blinded <- blind
        out[[length(out) + 1L]] <- x
        last <- cl
    }
    list(rows = out, outcome = last)
}

#' Run a phase's cells and score them into `runs/<label>/` (7.8 steps 5 to 9)
#'
#' The caller has checked the declaration and fetched the inputs (steps 1 to 4). Here: the ticket's
#' preparation, the kernel gate, the cells in concurrent chunks, blind or unblinded scoring, the
#' fixed tables, provenance, the evidence check (with a manifest) and `outputs.sha256` last.
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
#' @param regimes Regime table for the classifier (default the harness's).
#' @param arms,unit_ids Restrict the cells.
#' @param manifest The input manifest (data frame), for `attribution.tsv` and the evidence check.
#' @param code The provenance `code` block.
#' @param pkg_root The package checkout.
#' @return The run folder, invisibly.
#' @export
valRunPhase <- function(phase, decl, dataset, units, label, code_files, data_dir, work_dir,
                        out_dir, chunks = 4L, regimes = NULL, arms = NULL, unit_ids = NULL,
                        manifest = NULL, code = NULL, pkg_root = ".") {
    started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    units <- valUnits(units)
    env <- valTicketEnv(code_files)
    prep <- valPrepare(env, data_dir, file.path(work_dir, "prepare"))
    run_dir <- file.path(out_dir, "runs", label)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    cells <- valPhaseCells(decl, dataset, units, phase, unit_ids, arms)
    if (length(prep$result_tables)) {
        valCopyResultTables(prep$result_tables, decl[["outputs"]], run_dir,
                            arms = unique(cells$arm), units = unique(cells$unit),
                            universe = prep$universe)
    }
    no_text <- !is.null(manifest) &&
        any(manifest$kind == "data" & manifest$redistribution == "not_permitted")
    ctx <- list(decl = decl, dataset = dataset, units = prep$units, code_files = code_files,
                work_dir = file.path(work_dir, "cells"), pkg_root = normalizePath(pkg_root),
                record_messages = !no_text, out_dir = normalizePath(out_dir))
    gate <- valKernelGate(cells, c(ctx, list(env = env)))
    ctx_path <- file.path(work_dir, "cell-context.rds")
    saveRDS(ctx, ctx_path)
    runner <- file.path(work_dir, "cell-runner.R")
    writeLines(c(sprintf("ctx <- readRDS(%s)", deparse(normalizePath(ctx_path))),
                 "ctx$env <- msrval::valTicketEnv(ctx$code_files)",
                 "runCell <- function(cell) msrval::valRunCell(cell, ctx)"), runner)
    to_run <- cells[!cells$run %in% gate$refused, , drop = FALSE]
    minutes <- as.numeric(decl[["ceilings"]][["cell_minutes"]] %||% 120)
    ran <- if (nrow(to_run) > 0L) {
        valRunCells(to_run, runner, k = chunks, cell_minutes = minutes,
                    work_dir = file.path(work_dir, "chunks"))
    }
    res <- cells[, c("unit", "run", "arm")]
    res$exit <- "refused"
    res$stop_reason <- NA_character_
    res$step <- "identity_gate"
    res$class <- "val_kernel_identity_failed"
    res$message_sha256 <- NA_character_
    res$seconds <- NA_real_
    res$chunk <- NA_character_
    if (!is.null(ran)) {
        i <- match(ran$run, res$run)
        for (col in c("exit", "stop_reason", "step", "class", "message_sha256", "seconds",
                      "chunk")) {
            res[[col]][i] <- ran[[col]]
        }
    }
    res$fallback <- NA
    res$options_at_start <- NA_character_
    res$ruv_k <- NA_integer_
    res$ruv_controls <- NA_integer_
    cells_out <- res
    cells_out$run <- sub("@[^@]*$", "", cells_out$run)
    utils::write.table(cells_out[, EVIDENCE_COLUMNS[["cells.tsv"]]],
                       file.path(run_dir, "cells.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    blind_seed <- if (identical(phase, "lifecycle")) decl[["blind"]][["seed"]] else NULL
    done <- res[res$exit == "completed", , drop = FALSE]
    metrics <- do.call(rbind, lapply(unique(cells$unit), function(u) {
        r <- done[done$unit == u, , drop = FALSE]
        da <- stats::setNames(as.list(file.path(valCellDir(ctx$work_dir, r$run), "da.tsv")),
                              r$arm)
        valUnitMetrics(decl, dataset, u, da, prep$truth, prep$universe, env,
                       blind_seed = blind_seed)
    }))
    utils::write.table(metrics, file.path(run_dir, "metrics.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    within <- do.call(rbind, lapply(unique(cells$unit), function(u) {
        r <- done[done$unit == u, , drop = FALSE]
        da <- stats::setNames(as.list(file.path(valCellDir(ctx$work_dir, r$run), "da.tsv")),
                              r$arm)
        valWithinRunRows(decl, u, da, prep$truth, prep$universe, env, blind_seed = blind_seed)
    }))
    scored_units <- units[units$unit %in% cells$unit, , drop = FALSE]
    comparisons <- valRunLevel(decl, metrics, scored_units, regimes = regimes, phase = phase,
                               within = within)
    if (!is.null(within)) {
        key <- function(x) paste(x$component, x$scope, x$flags)
        comparisons <- rbind(comparisons, within[!key(within) %in% key(comparisons), ])
    }
    utils::write.table(comparisons, file.path(run_dir, "comparisons.tsv"), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    if (!is.null(manifest)) {
        utils::write.table(manifest[, c("role", "kind", "licence", "redistribution", "citation")],
                           file.path(run_dir, "attribution.tsv"), sep = "\t", quote = FALSE,
                           row.names = FALSE)
    }
    code <- code %||% list(package_commit = "none", release_commit = "none",
                           harness = list(commit = "none", tree = "none"), declaration = list())
    gates <- gate[setdiff(names(gate), "refused")]
    prov <- valProvenance(
        ticket = decl[["ticket"]] %||% "unknown", stage = decl[["stage"]] %||% phase,
        label = label, phase = phase, status = "completed", blinded = !is.null(blind_seed),
        started = started, code = code,
        inputs = list(prepared = prep$prepared),
        truth = list(ids = names(prep$truth), blinded = !is.null(blind_seed),
                     universe_keys = nrow(prep$universe)),
        cells = list(n = nrow(cells), completed = sum(res$exit == "completed"),
                     refused = sum(res$exit == "refused")),
        chunks = if (!is.null(ran)) attr(ran, "chunks") else list(),
        scoring = list(kernel_gates = lapply(gates, function(g) {
            list(passed = isTRUE(g$passed), units = g$units, note = g$note)
        }))
    )
    valWriteProvenance(prov, file.path(run_dir, "provenance.json"))
    valCheckProvenance(file.path(run_dir, "provenance.json"))
    if (!is.null(manifest)) {
        valCheckEvidence(run_dir, manifest, outputs = decl[["outputs"]] %||% list(),
                         arms = unique(cells$arm), units = unique(cells$unit),
                         universe = prep$universe, phase = phase)
    }
    valWriteOutputsSha(run_dir)
    invisible(run_dir)
}
