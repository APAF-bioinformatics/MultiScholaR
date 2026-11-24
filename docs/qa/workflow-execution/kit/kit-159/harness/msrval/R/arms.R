# Arm kinds' gates and hooks (KIT-159 sections 7.3, 7.5, 7.8).

#' Are two DA tables the same (digest, then values to 1e-12)?
#' @noRd
valSameDaTable <- function(a, b, tolerance = 1e-12) {
    if (identical(digest::digest(a, algo = "sha256"), digest::digest(b, algo = "sha256"))) {
        return(TRUE)
    }
    if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) return(FALSE)
    for (col in names(a)) {
        x <- a[[col]]
        y <- b[[col]]
        if (is.numeric(x) && is.numeric(y)) {
            if (!identical(is.na(x), is.na(y))) return(FALSE)
            if (any(abs(x - y) > tolerance, na.rm = TRUE)) return(FALSE)
        } else if (!identical(as.character(x), as.character(y))) {
            return(FALSE)
        }
    }
    TRUE
}

#' The kernel identity gate (section 7.5; A21, A30)
#'
#' On units drawn by the declared `identity_seed` (never by hand) plus one unit per refusal step
#' seen in the blind lifecycle, the kernel's DA table equals the workflow arm's (digest, then
#' values to 1e-12) and its refusals occur at the same step with the same class (or message
#' digest).
#'
#' @param units The unit ids.
#' @param kernel,workflow Functions of a unit returning `list(refused, step, class, da)`.
#' @param identity_seed The declared seed.
#' @param n Units drawn (at least three).
#' @param refusal_units Units of the refusal steps seen in the lifecycle.
#' @return A list: `passed`, `units`, `note` (`refusal_path_not_exercised` when no refusal unit).
#' @export
valKernelIdentity <- function(units, kernel, workflow, identity_seed, n = 3L,
                              refusal_units = character()) {
    pool <- setdiff(units, refusal_units)
    drawn <- valWithSeed(identity_seed, pool[sample.int(length(pool), min(n, length(pool)))])
    chosen <- sort(unique(c(drawn, refusal_units)))
    for (u in chosen) {
        w <- workflow(u)
        k <- kernel(u)
        same <- if (isTRUE(w$refused) || isTRUE(k$refused)) {
            isTRUE(w$refused) && isTRUE(k$refused) && identical(w$step, k$step) &&
                identical(w$class %||% w$message_sha256, k$class %||% k$message_sha256)
        } else {
            valSameDaTable(w$da, k$da)
        }
        if (!same) valAbort("the kernel differs from the workflow arm on unit {u}",
                            "val_kernel_identity_failed")
    }
    list(passed = TRUE, units = chosen,
         note = if (length(refusal_units)) "" else "refusal_path_not_exercised")
}

#' Call a post-cell hook outside test mode (section 7.3; revision 1.8 item 2, C22)
#'
#' The positional arguments are as in revision 1.6; `unit`, `arm`, `design_path` and `out_dir` are
#' passed only where the hook's formals name them explicitly (`...` does not receive them).
#'
#' @param hook The ticket's hook function.
#' @param manager,results,summary_dir,tables The positional arguments.
#' @param unit,arm,design_path,out_dir The named arguments.
#' @return The hook's value.
#' @export
valCallHook <- function(hook, manager, results, summary_dir, tables, unit = NULL, arm = NULL,
                        design_path = NULL, out_dir = NULL) {
    named <- list(unit = unit, arm = arm, design_path = design_path, out_dir = out_dir)
    named <- named[names(named) %in% names(formals(hook))]
    args <- c(list(manager, results, summary_dir, tables), named)
    withr::with_options(list(multischolar.test_mode = FALSE), {
        withr::with_envvar(c(MULTISCHOLAR_TEST_MODE = "false"), do.call(hook, args))
    })
}

#' Copy `prepareInputs()$result_tables` into the run folder and check them (revision 1.8 item 3)
#'
#' Each table must be a declared output with `key: "none"`; it is written under `run_dir` and
#' checked as a result table; an undeclared one is `val_publish_forbidden`.
#'
#' @param tables Named list: declared path = data frame.
#' @param outputs The declaration's `outputs`.
#' @param run_dir `runs/<label>/` inside `val-out/`.
#' @param arms,units Ids for placeholder expansion.
#' @param universe Keys of `universe.tsv` for a `not_permitted` lineage, or `NULL`.
#' @return The paths written.
#' @export
valCopyResultTables <- function(tables, outputs, run_dir, arms = character(), units = character(),
                                universe = NULL) {
    declared <- list()
    for (o in outputs) {
        if (!identical(o[["key"]] %||% "none", "none") || isTRUE(o[["handoff"]])) next
        for (p in valExpandOutputPaths(o, arms = arms, units = units)) declared[[p]] <- o
    }
    for (nm in names(tables)) {
        o <- declared[[nm]]
        if (is.null(o)) {
            valAbort("result table {nm} is not a declared output", "val_publish_forbidden")
        }
        path <- file.path(run_dir, nm)
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        utils::write.table(tables[[nm]], path, sep = "\t", quote = FALSE, row.names = FALSE)
        cols <- vapply(o[["columns"]], function(x) x[["type"]], character(1))
        names(cols) <- vapply(o[["columns"]], function(x) x[["name"]], character(1))
        labels <- lapply(o[["labels"]] %||% list(), function(x) as.character(unlist(x)))
        valCheckResultTable(path, nm, cols, labels, o[["max_rows"]] %||% 1000L, universe)
    }
    names(tables)
}

# ---------------------------------------------------------------------------------------------
# Workflow arm builders (section 7.5) and their inputs per route; the custom and published kinds.
# ---------------------------------------------------------------------------------------------

RUV_DESIGNS <- c("groups_all_features", "technical_replicates", "first_pass_nulls", "legacy",
                 "none", "legacy_correction", "supplied_controls")
RUV_FACTORS_IN_MODEL <- c("groups_all_features", "technical_replicates", "supplied_controls",
                          "first_pass_nulls")

#' A workflow arm
#' @noRd
valWorkflowArm <- function(id, builder, ...) {
    list(id = id, kind = "workflow", builder = builder, ...)
}

#' The release as a default user runs it: no option, no override
#' @return An arm.
#' @export
valArmShipped <- function() valWorkflowArm("shipped", "shipped")

#' The shipped legacy default, pinned: gate off, the mode and its parameters written explicitly
#' @return An arm.
#' @export
valArmRuvLegacy <- function() valWorkflowArm("ruv_legacy", "legacy")

#' The calibrated design (VAL-3, VAL-13): gate on, the design the panel selects (recorded as
#' expected), the k rule and seed at the panel's defaults
#' @return An arm.
#' @export
valArmRuvCalibrated <- function() valWorkflowArm("ruv_calibrated", "calibrated")

#' RUV off: `skip` with the gate off; transcriptomics design `none` with the gate on
#' @return An arm.
#' @export
valArmRuvOff <- function() valWorkflowArm("ruv_off", "off")

#' Any RUV design the route offers, with the gate on
#' @param design One of the designs of section 7.5.
#' @param k_rule The k rule, or `NULL` for the panel's default.
#' @return An arm.
#' @export
valArmRuvDesign <- function(design, k_rule = NULL) {
    if (!design %in% RUV_DESIGNS) valAbort("RUV design {design} is not offered", "val_spec_invalid")
    valWorkflowArm(paste0("ruv_", design), "design", design = design, k_rule = k_rule)
}

#' The supplied-controls design: the list from a declared or prepared role, one per line (A42)
#' @param controls_role The role holding the control identifiers.
#' @return An arm.
#' @export
valArmRuvSupplied <- function(controls_role) {
    valWorkflowArm("ruv_supplied_controls", "supplied", design = "supplied_controls",
                   controls_role = controls_role)
}

#' The input-id prefix of a route's normalisation module
#' @noRd
valNormPrefix <- function(route) {
    if (identical(valRoute(route)$omic, "proteomics")) "normalization-" else "norm-"
}

#' An arm's options, overrides and expected server values on a route
#'
#' Legacy values are those of the panels at a60c4e9a (`R/utils_omic_gui_config.R`: proteomics
#' `ruv_iii` default `automatic`, k NULL; lipidomics and metabolomics `automatic`, 5 %, k 2;
#' transcriptomics grid `automatic`; DIA_LIMPA `covariates`, set by the server); the calibrated
#' design is the table's `selected` (`groups_all_features`, `technical_replicates` when the design
#' has technical replicates). G3 re-reads them in the lifecycle.
#'
#' @param arm A workflow arm.
#' @param route Route id.
#' @param roles Identifier lists by role (supplied controls).
#' @param technical_replicates The design has technical replicates.
#' @return A list: `options`, `overrides`, `expect`, `mask`.
#' @export
valArmInputs <- function(arm, route, roles = list(), technical_replicates = FALSE) {
    out <- list(options = list(), overrides = list(), expect = list(), mask = character())
    if (identical(arm$builder, "shipped")) return(out)
    r <- valRoute(route)
    p <- valNormPrefix(route)
    id <- function(x) paste0(p, x)
    gate <- function(on) list(MultiScholaR.ruv_design_panel = on)
    switch(arm$builder,
           legacy = {
               out$options <- gate(FALSE)
               out$overrides <- if (identical(r$workflow, "DIA_LIMPA")) {
                   stats::setNames(list("covariates"), id("ruv_mode"))
               } else if (r$omic %in% c("lipidomics", "metabolomics")) {
                   stats::setNames(list("automatic", 5, 2), id(c("ruv_mode", "ruv_percentage",
                                                                 "ruv_k")))
               } else {
                   stats::setNames(list("automatic"), id("ruv_mode"))
               }
           },
           calibrated = {
               out$options <- gate(TRUE)
               cal_design <- if (technical_replicates) "technical_replicates" else
                   "groups_all_features"
               out$expect <- stats::setNames(list(cal_design), id("ruv_design"))
           },
           off = {
               if (identical(r$omic, "transcriptomics")) {
                   out$options <- gate(TRUE)
                   out$overrides <- stats::setNames(list("none"), id("ruv_design"))
               } else {
                   out$options <- gate(FALSE)
                   out$overrides <- stats::setNames(list("skip"), id("ruv_mode"))
               }
           },
           design = {
               out$options <- gate(TRUE)
               out$overrides <- stats::setNames(list(arm$design), id("ruv_design"))
               if (!is.null(arm$k_rule)) out$overrides[[id("ruv_k_rule")]] <- arm$k_rule
           },
           supplied = {
               ids <- roles[[arm$controls_role]]
               if (!length(ids)) {
                   valAbort("no identifiers for the controls role {arm$controls_role}",
                            "val_spec_invalid")
               }
               out$options <- gate(TRUE)
               out$overrides <- stats::setNames(list("supplied_controls",
                                                     paste(ids, collapse = "\n")),
                                                id(c("ruv_design", "ruv_control_list")))
               out$mask <- id("ruv_control_list")
           },
           valAbort("unknown arm builder {arm$builder}", "val_spec_invalid"))
    out
}

#' Where a cell's RUV puts its factors (the stage rule of 7.3): `off`, `correction` or
#' `factors_in_model`
#'
#' @param route Route id.
#' @param arm A workflow arm.
#' @return One of the three.
#' @export
valCellRuv <- function(route, arm) {
    r <- valRoute(route)
    b <- arm$builder %||% "shipped"
    design <- arm$design
    if (identical(b, "off") || identical(design, "none")) return("off")
    if (identical(r$omic, "transcriptomics")) return("factors_in_model")
    if (identical(r$workflow, "DIA_LIMPA")) {
        return(if (identical(design, "legacy_correction")) "correction" else "factors_in_model")
    }
    if (r$omic %in% c("lipidomics", "metabolomics")) {
        if (identical(b, "calibrated")) return("factors_in_model")
        if (!is.null(design) && design %in% RUV_FACTORS_IN_MODEL) return("factors_in_model")
        return("correction")
    }
    "correction"
}

#' Whether an arm skips RUV, always applies it, or may give either state (I09, I18)
#'
#' @param arm A workflow arm.
#' @return `skip` (RUV off, or design `none`), `applied` (the pinned legacy mode) or `either`.
#' @export
valArmRuvRun <- function(arm) {
    b <- arm$builder %||% "shipped"
    if (identical(b, "off") || identical(arm$design, "none")) return("skip")
    if (identical(b, "legacy")) return("applied")
    "either"
}

#' Run an arm of kind `custom` or `published` on one unit (section 7.5; 4.7)
#'
#' A `custom` arm is labelled "not the harness baseline". A `published` arm (a published
#' processing re-run) must declare its `anchor` check; on failure `anchor_on_fail = "stop"` (the
#' default) refuses the arm (`val_arm_refused`) and `"context"` keeps it, labelled "anchor not
#' reproduced: context only", never primary (KIT-163 N6). A ticket function that signals an error
#' is a refusal with the condition's class and the SHA-256 of its masked message.
#'
#' @param arm List: `id`, `kind`, `fn(unit, inputs)`, and for `published` `anchor(da)` and
#'   `anchor_on_fail`.
#' @param unit Unit id.
#' @param inputs The unit's inputs.
#' @return A list: `refused`, `class`, `message_sha256`, `da`, `label`, `anchor`, `primary`.
#' @export
valRunCustomArm <- function(arm, unit, inputs) {
    kind <- arm$kind %||% ""
    if (!kind %in% c("custom", "published") || !is.function(arm$fn)) {
        valAbort("arm {arm$id} is not a custom or published arm with a function",
                 "val_spec_invalid")
    }
    if (identical(kind, "published") && !is.function(arm$anchor)) {
        valAbort("published arm {arm$id} declares no anchor check", "val_spec_invalid")
    }
    label <- if (identical(kind, "custom")) "not the harness baseline" else
        "published processing re-run"
    da <- tryCatch(arm$fn(unit, inputs), error = function(e) e)
    if (inherits(da, "error")) {
        return(list(refused = TRUE, class = class(da)[[1L]],
                    message_sha256 = digest::digest(valMaskMessage(conditionMessage(da)),
                                                    algo = "sha256", serialize = FALSE),
                    da = NULL, label = label, anchor = NA_character_, primary = FALSE))
    }
    anchor <- NA_character_
    primary <- TRUE
    if (identical(kind, "published")) {
        if (isTRUE(arm$anchor(da))) {
            anchor <- "reproduced"
        } else if (identical(arm$anchor_on_fail %||% "stop", "context")) {
            anchor <- "not_reproduced"
            label <- "anchor not reproduced: context only"
            primary <- FALSE
        } else {
            # I20: refused at step `anchor`, so cells.tsv records the step and class.
            valArmRefused("anchor", "val_anchor_not_reproduced")
        }
    }
    list(refused = FALSE, class = NA_character_, message_sha256 = NA_character_, da = da,
         label = label, anchor = anchor, primary = primary)
}
