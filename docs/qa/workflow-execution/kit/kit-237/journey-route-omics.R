# KIT-237 (r1.2, OD-1 (b), owner 2026-11-25, VAL-225): the three further Windows routes for KIT-233's journey machinery,
# one route table each, on the fixtures of the existing e2e lanes (tests/testdata/e2e/manifest.json):
#   lipid_canonical  lipidomics, LipidSearch LCMS_Pos + LCMS_Neg, 6 samples (WT_1..3, KO_1..3), contrast KO_vs_WT;
#                    entry :335-370, SHA-256 design.tsv :365, lipidsearch_lcms_neg.txt :367, lipidsearch_lcms_pos.txt :368
#   metab_lc         metabolomics, custom LCMS_Pos + LCMS_Neg feature tables, 6 samples, contrast KO_vs_WT;
#                    entry :221-260, SHA-256 design.tsv :251, lcms_neg_features.tsv :253, lcms_pos_features.tsv :254
#   trx_screen       transcriptomics, counts.tsv, 12 samples (RPMI, Sera), contrast Sera_vs_RPMI;
#                    entry :372-410, SHA-256 counts.tsv :407, design.tsv :409
# Sourced by journey-win.R after journey.R and journey-route-lfq.R (which defines k237StepHome).
#
# journey.R's step functions are the proteomics shell's (namespace `proteomics_workflow-`, the proteomics digest keys
# and stepper). These routes therefore carry their own step functions, in journey.R's form (one k233Action per press,
# its classifier, captures and records), over the shared lipid/metab/trx shell: tabset `<omic>_workflow-<omic>_tabs`
# with values import, design, qc, norm, de, summary; test ids `<prefix>-...` (R/utils_omic_gui_config.R:1079-1121,
# :1275-1300, :1953-1960; R/mod_omic_summary.R:60). Import settings are the lane's (the user's own choices of assay
# names, format and columns: tests/testthat/helper-e2e-browser.R:809-866; test-e2e-transcriptomics-screen.R:31-33); every
# later step runs at the rendered defaults, and a default refused on the fixture gets the lane's value for that step
# only, pressed again (WL-1's rule). Completion: the digest's step status or R6 state in test mode; the DOM signal is
# the application's own message notification (R/mod_omic_design.R:260, R/mod_omic_qc_tabs.R:753, R/mod_omic_norm.R:349,
# :482, :549, :650, R/mod_omic_da.R:514, R/mod_omic_summary.R:259-332, R/utils_workflow_summary.R:439, :468). These DOM
# signals are by reading only; WL-5 on these routes is their first check.
#
# r1.3 (the script check before WL, AUDIT-KIT-237-plan-20261125.md "Script check before WL"):
#   SC-01  lipid and metab follow the owner's default-user scripts (docs/qa/workflow-execution/kit/kit-159/
#          action-scripts.md:51-52, :10-12): QC presses the intensity filter, duplicate detection and resolution and the
#          internal-standard analysis before Finalize, each at its rendered defaults; `sanitize_names` stays at its
#          rendered TRUE (the lane's FALSE is dropped). The ids and states are the default-user script's
#          (kit-159/harness/msrval/R/scripts.R:316-340). The fixtures hold no internal standards (no row of either
#          lipid assay or metab assay matches; known.tsv K12, KIT-242), so the analysis is expected to refuse with
#          "No internal standards found" (R/func_omics_qc_shared.R:279); that refusal is recorded and the route goes on.
#   SC-02  user mode: each input chosen with the shinyFiles chooser (`<omic>_workflow-import-<slot>_file`,
#          R/mod_omic_import_feature_table.R:105, :453-460), as the label-free route does.
#   SC-04  the DOM signal for import, design, QC, export, DA (and normalisation on trx) is the omic's own stepper
#          (`<omic>_workflow-workflow_progress`, R/mod_omic_workflow.R:40, :119), read by journey.R's state script
#          rebound per session in journey-win.R (k237SetOmic); the trx filter's is its exact text
#          ("Expression filter applied:", R/mod_trx_qc_filter.R:57); the others keep their exact texts.
#   Ruling (coordinator, 2026-11-25): trx_screen's design assigns its technical-replicate sets as the default-user
#          script does, and a named check records the builder's sample keys (known KIT-222, known.tsv K13).
#   SC-05  screen outputs named (import status; the DA status, summary statistics and results table); N2 computed on
#          metab_lc; WL-4 (a) planted once on metab_lc (journey-win.R).

K237_ROUTES <- if (exists("K237_ROUTES", inherits = FALSE)) K237_ROUTES else list()
K237_RUNS <- if (exists("K237_RUNS", inherits = FALSE)) K237_RUNS else list()

k237OmicRoute <- function(id, omic, prefix, qc_tabset, trx = FALSE) {
    steps <- list(
        list(id = "S1", name = "home and project (the rendered default base directory kept)", fn = "k237StepHome",
            tile = paste0("tile-", omic), omic = omic, tab_testid = paste0(prefix, "-tab-import")),
        list(id = "S2", name = "import (the lane's assay names, format and files)", fn = "k237StepImportOmic"),
        list(id = "S3", name = "design", fn = "k237StepDesignOmic"),
        list(id = "S5", name = "quality control", fn = "k237StepQcOmic"),
        list(id = "S6", name = "normalisation", fn = "k237StepNormOmic"))
    if (!trx) steps <- c(steps, list(list(id = "S7", name = "export", fn = "k237StepExportOmic")))
    steps <- c(steps, list(
        list(id = "S8", name = "differential analysis", fn = "k237StepDaOmic"),
        list(id = "S10", name = "summary", fn = "k237StepSummaryOmic")))
    list(id = id, omic = omic, prefix = prefix, trx = trx, tabs_id = sprintf("%s_workflow-%s_tabs", omic, omic),
        qc_tabset = qc_tabset, steps = steps)
}

K237_ROUTES$lipid_canonical <- k237OmicRoute("LIPIDOMICS", "lipidomics", "lipid", "lipid_qc_tabs")
K237_ROUTES$metab_lc <- k237OmicRoute("METABOLOMICS", "metabolomics", "metab", "metab_qc_tabs")
K237_ROUTES$trx_screen <- k237OmicRoute("TRANSCRIPTOMICS", "transcriptomics", "trx", "trx_qc_tabs", trx = TRUE)

#' Named files (5.2) by reading: the design files (test-e2e-lipidomics-canonical.R:18-21), the DA tables
#' (:51-55), study parameters, the publication workbook and the session export (:91-113), the report (the template
#' copied to scripts/<omic>, as on proteomics), and the lane's `expected_exports` matched by name anywhere.
k237OmicFiles <- function(omic, exports, trx = FALSE) {
    src <- sprintf("^scripts/%s/", omic)
    prefix <- c(lipidomics = "lipid", metabolomics = "metab", transcriptomics = "trx")[[omic]]
    files <- list(
        S3 = paste0(src, c("design_matrix\\.tab$", "contrast_strings\\.tab$")),
        S7 = paste0(src, prefix, c("_filtered_session_data_latest\\.rds$", "_filtered_session_summary\\.txt$")),
        S8.2 = sprintf("^results/%s/[^/]*/?(de|da)_[^/]*_long_annot\\.(tsv|xlsx)$", omic),
        S10.1 = paste0(src, "study_parameters\\.txt$"),
        S10.2 = sprintf("(^|/)Publication_tables/DA_results_%s\\.xlsx$", omic),
        S10.3 = sprintf("%s%s_report\\.rmd$", src, omic),
        S10.5 = c(paste0(src, "session_state_[^/]*\\.RDS$"), sprintf("(^|/)%s$", gsub(".", "\\.", exports, fixed = TRUE))))
    if (trx) {
        # WL 45605975875: the trx DA writes results/transcriptomics/da_genes/de_results_<contrast>.csv; trx has no S7.
        files$S8.2 <- NULL
        files$S7 <- NULL
        files$S8 <- "^results/transcriptomics/da_genes/de_results_[^/]*\\.csv$"
    }
    files
}

k237OmicRun <- function(id, omic, dataset, contrast, import, exports, trx = FALSE) {
    force(import)
    function(G, dir, label) {
        design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE)
        groups <- unique(design$group)
        list(dataset = dataset, omic = omic, groups = groups,
            members = lapply(stats::setNames(groups, groups), function(g) design$sample[design$group == g]),
            contrasts = list(contrast), label = label, description = sprintf("KIT-237 %s journey on %s", omic, dataset),
            permissive = TRUE, project_dir = NULL, import = import(dir), files = k237OmicFiles(omic, exports, trx),
            last_step = "S10.5", summary_step = "S10", fixture_dir = dir,
            # Coordinator's ruling (2026-11-25): where the design carries technical-replicate sets (trx_screen's
            # tech_rep_group), the design step assigns them as the default-user script does (kit-159/action-scripts.md,
            # Design; tests/testthat/helper-gui-entry.R:452-462): each set of two or more runs, mode "manual", the
            # replicate number of the set's first run in design.tsv.
            tech_reps = if (!is.null(design$tech_rep_group)) Filter(function(x) length(x$samples) >= 2L,
                lapply(unique(design$tech_rep_group), function(set) {
                    samples <- design$sample[design$tech_rep_group == set]
                    list(set = set, samples = samples, replicate = as.integer(design$replicates[match(samples[[1L]],
                        design$sample)]))
                })),
            biological = if (!is.null(design$tech_rep_group)) stats::setNames(design$tech_rep_group, design$sample))
    }
}

# Import actions: list(kind = "set" | "upload", input, value); ids under `<omic>_workflow-import-`.
K237_RUNS$lipid_canonical <- k237OmicRun("lipid_canonical", "lipidomics", "lipid_canonical", c("KO", "WT"),
    function(dir) list(
        list(kind = "set", input = "assay1_name", value = "LCMS_Pos"),
        list(kind = "set", input = "vendor_format", value = "lipidsearch"),
        list(kind = "upload", input = "assay1_file_std", value = file.path(dir, "lipidsearch_lcms_pos.txt")),
        list(kind = "set", input = "assay2_name", value = "LCMS_Neg"),
        list(kind = "upload", input = "assay2_file_std", value = file.path(dir, "lipidsearch_lcms_neg.txt"))),
    c("lipid_summary.xlsx", "differential_abundance.xlsx"))

K237_RUNS$metab_lc <- k237OmicRun("metab_lc", "metabolomics", "metab_lc", c("KO", "WT"),
    function(dir) list(
        list(kind = "set", input = "assay1_name", value = "LCMS_Pos"),
        list(kind = "set", input = "vendor_format", value = "custom"),
        list(kind = "upload", input = "assay1_file_std", value = file.path(dir, "lcms_pos_features.tsv")),
        list(kind = "set", input = "assay2_name", value = "LCMS_Neg"),
        list(kind = "upload", input = "assay2_file_std", value = file.path(dir, "lcms_neg_features.tsv")),
        list(kind = "set", input = "metabolite_id_col_custom", value = "Feature.Name"),
        list(kind = "set", input = "annotation_col_custom", value = ""),
        list(kind = "set", input = "sample_cols_pattern", value = "^(WT|KO)_"),
        list(kind = "set", input = "is_pattern", value = "^ITSD|_d[0-9]+$")),
    c("feature_summary.xlsx", "differential_abundance.xlsx"))

K237_RUNS$trx_screen <- k237OmicRun("trx_screen", "transcriptomics", "trx_screen", c("Sera", "RPMI"),
    function(dir) list(list(kind = "upload", input = "assay1_file_std", value = file.path(dir, "counts.tsv"))),
    c("DA_results_transcriptomics.xlsx", "transcriptomics_normalized_data.tsv", "ruviii_before_and_after.png"),
    trx = TRUE)

# --- helpers ------------------------------------------------------------------------------------------------------

k237Id <- function(route, module, input) sprintf("%s_workflow-%s-%s", route$omic, module, input)
k237Tab <- function(route, value) list(route$tabs_id, value)

k237Status <- function(st, omic, key) {
    value <- st$dg$step_status_per_omic[[omic]][[key]]
    if (is.null(value)) NA_character_ else as.character(value)
}

k237HasState <- function(st, omic, states) {
    head <- st$dg$r6_current_state_per_omic[[omic]]
    any(c(as.character(head), as.character(unlist(st$dg$r6_state_history_per_omic[[omic]]))) %in% states)
}

#' One press with journey.R's classifier; on a refusal and only where the route gives the lane's values, those values
#' are set and the same press made once more (WL-1's rule).
k237Press <- function(S, spec, permissive = NULL) {
    res <- k233Action(S, spec)
    # The application asks for the press again when another step is still running (WL 45605975875: the normalisation
    # tab's own entry QC was running at S6.1): a user waits and presses again with the same settings; the lane's values
    # are only for a default refused on the fixture.
    if (identical(res$outcome, "refused") && grepl("still running", res$message %||% "", fixed = TRUE)) {
        Sys.sleep(30)
        again <- spec
        again$step <- paste0(spec$step, "w")
        again$name <- paste(spec$name, "(pressed again, as the application asks, same settings)")
        again$prepare <- NULL
        res <- k233Action(S, again)
    }
    if (identical(res$outcome, "refused") && isTRUE(S$R$permissive) && is.function(permissive)) {
        retry <- spec
        retry$step <- paste0(spec$step, "r")
        retry$name <- paste(spec$name, "(the lane's value, WL-1's rule)")
        retry$prepare <- permissive
        res <- k233Action(S, retry)
    }
    res
}

# --- steps --------------------------------------------------------------------------------------------------------

k237StepImportOmic <- function(S, R, step, route) {
    prepare <- function(S) {
        for (a in R$import) {
            id <- k237Id(route, "import", a$input)
            if (identical(a$kind, "upload") && !S$test_mode) {
                # User mode (SC-02): the slot's shinyFiles button, the file copied where the chooser reaches it.
                chosen <- k233ChooseFile(S, "S2", sub("_std$", "", id), k237UserInput(R, a$value))
                S$chooser <- c(S$chooser, isTRUE(chosen$ok))
                k233Observation(S$G, S$label, paste("S2 chooser", a$input), chosen$how)
                if (!chosen$ok) return(list(refuse = paste("not reachable as a user: chooser", a$input, ":", chosen$how)))
            } else if (identical(a$kind, "upload")) {
                up <- k233Upload(S, "S2", id, a$value)
                if (!up$ok) return(list(refuse = paste("upload", basename(a$value), ":", up$text)))
            } else {
                k233Set(S, "S2", id, a$value)
            }
        }
        NULL
    }
    k233Action(S, list(step = "S2", name = step$name, target = paste0(route$prefix, "-import-process"), by = "testid",
        tabs = list(k237Tab(route, "import")), prepare = prepare,
        digest = function(st) identical(k237Status(st, route$omic, "setup_import"), "complete"),
        dom = function(st) k233Stepper(st, 1L),
        outputs = if (!route$trx) k237Id(route, "import", "import_status")))
}

#' The design builder as journey.R's S3 does it (factors, metadata per group, the contrast, Save), on the omic's
#' builder ids.
k237StepDesignOmic <- function(S, R, step, route) {
    builder <- function(x) k237Id(route, "design", paste0("builder-", x))
    tabs <- builder("main_tabs")
    prepare <- function(S) {
        if (!k233WaitFor(S, sprintf("document.getElementById(%s) !== null", k233Json(builder("selected_runs"))), 60)) {
            return(list(refuse = "the design builder was not shown within 60 s"))
        }
        k233ActivateTab(S, tabs, "Factors")
        for (g in R$groups) {
            k233Set(S, "S3", builder("new_factor"), g)
            p <- k233Press(S, "S3", builder("add_factor"))
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        k233ActivateTab(S, tabs, "Assign Metadata")
        options <- k233Options(S, builder("selected_runs"))
        for (g in R$groups) {
            runs <- options[k233Norm(options) %in% k233Norm(R$members[[g]])]
            if (length(runs) != length(R$members[[g]])) {
                return(list(refuse = sprintf("driver: %d of %d samples of %s found among the builder's %d", length(runs),
                    length(R$members[[g]]), g, length(options))))
            }
            k233Selectize(S, "S3", builder("selected_runs"), runs)
            k233WaitFor(S, sprintf("document.getElementById(%s) !== null", k233Json(builder("replicate_start"))), 30)
            tryCatch(S$G$helpers$e2e_wait_for_select_options(S$driver, builder("factor1_select"), g, timeout = 30000),
                error = function(e) NULL)
            k233Selectize(S, "S3", builder("factor1_select"), g, multiple = FALSE)
            p <- k233Press(S, "S3", builder("assign_metadata"))
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        if (length(R$tech_reps)) {
            done <- k237AssignTechReps(S, R, builder, tabs)
            if (!is.null(done)) return(done)
        }
        k233ActivateTab(S, tabs, "Contrasts")
        for (pair in R$contrasts) {
            choices <- k233Options(S, builder("contrast_group1"))
            matched <- S$G$helpers$e2e_match_available_values(pair, choices)
            k233Selectize(S, "S3", builder("contrast_group1"), matched[[1L]], multiple = FALSE)
            k233Selectize(S, "S3", builder("contrast_group2"), matched[[2L]], multiple = FALSE)
            p <- k233Press(S, "S3", builder("add_contrast"))
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        NULL
    }
    k233Action(S, list(step = "S3", name = step$name, target = builder("save_results"),
        tabs = list(k237Tab(route, "design")), prepare = prepare,
        digest = function(st) identical(k237Status(st, route$omic, "design_matrix"), "complete"),
        dom = function(st) k233Stepper(st, 2L),
        after = if (length(R$tech_reps)) function(S, res, cap) k237SampleKeyCheck(S, R, route)))
}

#' The Technical Replicates tab of the builder (R/mod_omic_design_builder.R:127-142): per set, the runs selected in
#' `builder-tech_rep_samples`, `builder-tech_rep_assignment_mode` = "manual", `builder-manual_replicate_number` = the
#' set's replicate number, then `builder-assign_tech_reps`. Completion of each assignment: the builder's
#' `builder-tech_rep_summary` text names the set's first run within 30 s (its message "Assigned N samples as technical
#' replicates", R/func_omics_design_shared.R:652-653, is recorded with the step's notifications). Returns NULL, or a
#' refusal for the step.
k237AssignTechReps <- function(S, R, builder, tabs) {
    k233ActivateTab(S, tabs, "Technical Replicates")
    options <- k233Options(S, builder("tech_rep_samples"))
    summary_text <- function() as.character(k233Js(S, sprintf("(document.getElementById(%s)||{textContent:''}).textContent",
        k233Json(builder("tech_rep_summary")))) %||% "")
    for (set in R$tech_reps) {
        runs <- options[k233Norm(options) %in% k233Norm(set$samples)]
        if (length(runs) != length(set$samples)) {
            return(list(refuse = sprintf("driver: %d of %d runs of technical-replicate set %s found among the builder's %d",
                length(runs), length(set$samples), set$set, length(options))))
        }
        k233Selectize(S, "S3", builder("tech_rep_samples"), runs)
        k233Set(S, "S3", builder("tech_rep_assignment_mode"), "manual")
        k233Set(S, "S3", builder("manual_replicate_number"), set$replicate)
        p <- k233Press(S, "S3", builder("assign_tech_reps"))
        if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
        if (!k233WaitFor(S, sprintf("(document.getElementById(%s)||{textContent:''}).textContent.indexOf(%s) >= 0",
            k233Json(builder("tech_rep_summary")), k233Json(runs[[1L]])), 30)) {
            return(list(refuse = sprintf("the technical-replicate set %s was not shown as assigned within 30 s: %s",
                set$set, summary_text())))
        }
        k233Idle(S, 10)
    }
    k233Observation(S$G, S$label, "S3 technical-replicate assignments (builder summary)", summary_text())
    NULL
}

#' The named check after S3 (coordinator's ruling; KIT-222): the saved design's sample keys
#' (scripts/<omic>/design_matrix.tab, `tech_rep_group` = group_replicates, R/func_omics_design_shared.R:680-687) against
#' the biological samples of design.tsv (its tech_rep_group). Written to `<label>-S3-sample-keys.tsv` (run, biological
#' sample, builder key); one row in numbers.tsv (K222: distinct keys against distinct biological samples); a key shared
#' by two biological samples is a kind-2 finding, known as KIT-222 (known.tsv K13).
k237SampleKeyCheck <- function(S, R, route) {
    file <- file.path(k233ProjectPaths(R)$source, "design_matrix.tab")
    if (!file.exists(file)) {
        k233Finding(S, "S3", 3L, paste("the saved design is missing:", file))
        return(invisible(NULL))
    }
    saved <- utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
    file.copy(file, file.path(S$G$out, paste0(S$label, "-design_matrix.tab.txt")), overwrite = TRUE)
    run_col <- intersect(c("Run", "run", "sample"), names(saved))[1L]
    if (is.na(run_col) || !"tech_rep_group" %in% names(saved)) {
        k233Observation(S$G, S$label, "S3 sample keys", paste("the saved design has no Run or tech_rep_group column:",
            paste(names(saved), collapse = ", ")))
        return(invisible(NULL))
    }
    bio <- R$biological[match(k233Norm(saved[[run_col]]), k233Norm(names(R$biological)))]
    keys <- as.character(saved$tech_rep_group)
    k233WriteTsv(data.frame(run = saved[[run_col]], biological_sample = bio, builder_key = keys,
        stringsAsFactors = FALSE), file.path(S$G$out, paste0(S$label, "-S3-sample-keys.tsv")))
    shared <- Filter(function(k) length(unique(stats::na.omit(bio[keys == k]))) > 1L, unique(keys))
    n_keys <- length(unique(keys))
    n_bio <- length(unique(stats::na.omit(bio)))
    S$G$numbers[[length(S$G$numbers) + 1L]] <- list(session = S$label, id = "K222", contrast = "sample keys",
        shown = n_keys, saved = n_bio, rule = "no key shared by two biological samples",
        result = if (length(shared)) "differ" else if (anyNA(bio)) "not_computed" else "equal",
        where = sprintf("%s tech_rep_group vs design.tsv tech_rep_group", basename(file)), planted = FALSE)
    for (k in shared) {
        k233Finding(S, "S3", 2L, sprintf("builder key %s holds two different biological samples (%s): KIT-222", k,
            paste(unique(stats::na.omit(bio[keys == k])), collapse = ", ")))
    }
    invisible(shared)
}

#' QC (SC-01): the rendered sub-tabs recorded; on trx the expression filter's Apply, then Finalize QC; on lipid and
#' metab the intensity filter, duplicate detection and resolution, the internal-standard analysis, then Finalize QC,
#' each at its rendered defaults in rendered order (the default-user script, action-scripts.md:51-52).
k237StepQcOmic <- function(S, R, step, route) {
    qc <- function(x) k237Id(route, "qc", x)
    tabset <- qc(route$qc_tabset)
    k233ActivateTab(S, route$tabs_id, "qc")
    titles <- character()
    started <- Sys.time()
    while (!length(titles) && k233Seconds(started) < 60) {
        titles <- as.character(unlist(k233Js(S, sprintf(K233_TITLES_JS, k233Json(tabset)))))
        if (!length(titles)) Sys.sleep(1)
    }
    k233Observation(S$G, S$label, "S5 rendered QC sub-tabs", paste(titles, collapse = " > "))
    state <- function(x) paste0(route$prefix, "_", x)
    sub_tab <- function(title) list(k237Tab(route, "qc"), list(tabset, title))
    message_seen <- function(pattern) function(st) k233Noted(st, pattern, "message")
    actions <- if (route$trx) list(
        list(step = "S5.1", name = "expression filter at the defaults", target = "trx-qc-apply-filter", by = "testid",
            tabs = list(k237Tab(route, "qc")), digest = function(st) k237HasState(st, route$omic, "trx_expression_filtered"),
            dom = message_seen("Expression filter applied:")),
        list(step = "S5.2", name = "Finalize QC", target = qc("s4_finalize-finalize_qc"), tabs = sub_tab("Finalize QC"),
            digest = function(st) k237HasState(st, route$omic, "trx_qc_complete"), dom = function(st) k233Stepper(st, 3L))
    ) else list(
        list(step = "S5.1", name = "intensity filter at the defaults", target = qc("intensity-apply_filter"),
            tabs = sub_tab("Intensity Filter"), digest = function(st) k237HasState(st, route$omic, state("intensity_filtered")),
            dom = message_seen("Intensity filter applied:")),
        list(step = "S5.2", name = "detect duplicates", target = qc("duplicates-detect_duplicates"),
            tabs = sub_tab("Duplicate Resolution"), dom = message_seen("Detection complete:")),
        list(step = "S5.3", name = "resolve duplicates", target = qc("duplicates-resolve_duplicates"),
            tabs = sub_tab("Duplicate Resolution"),
            digest = function(st) k237HasState(st, route$omic, state("duplicates_resolved")),
            dom = message_seen("Duplicates resolved:")),
        list(step = "S5.4", name = "analyse internal standards at the rendered pattern", target = qc("itsd-analyze_is"),
            tabs = sub_tab("Internal Standards"), dom = message_seen("Found [0-9]+ internal standards"), itsd = TRUE),
        list(step = "S5.5", name = "Finalize QC", target = qc("s4_finalize-finalize_qc"), tabs = sub_tab("Finalize QC"),
            digest = function(st) k237HasState(st, route$omic, state("qc_complete")), dom = function(st) k233Stepper(st, 3L)))
    res <- list(outcome = "complete")
    for (a in actions) {
        itsd <- isTRUE(a$itsd)
        a$itsd <- NULL
        res <- k233Action(S, a)
        # The fixtures hold no internal standards (known.tsv K12, KIT-242): that refusal is recorded and QC goes on.
        if (itsd && identical(res$outcome, "refused") && grepl("internal standards", res$message %||% "", ignore.case = TRUE)) {
            k233Observation(S$G, S$label, "S5.4 internal standards (known K12, KIT-242)", res$message)
            next
        }
        if (!identical(res$outcome, "complete")) return(res)
    }
    res
}

#' Normalisation at the rendered defaults (method, RUV mode, ITSD as rendered, recorded); then, except on trx, the
#' correlation filter at the rendered threshold. The lane's values (ITSD off, method none, RUV skip; skip the filter)
#' only after a refusal (WL-1's rule).
k237StepNormOmic <- function(S, R, step, route) {
    norm <- function(x) k237Id(route, "norm", x)
    tabs <- list(k237Tab(route, "norm"))
    prepare <- function(S) {
        k233Observation(S$G, S$label, "S6 rendered defaults", sprintf("apply_itsd=%s norm_method=%s ruv_mode=%s",
            k233Clean(k233InputValue(S, norm("apply_itsd"))), k233Clean(k233InputValue(S, norm("norm_method"))),
            k233Clean(k233InputValue(S, norm("ruv_mode")))))
        NULL
    }
    lane <- function(S) {
        if (!route$trx) k233Set(S, "S6.1r", norm("apply_itsd"), FALSE)
        k233Set(S, "S6.1r", norm("norm_method"), "none")
        k233Set(S, "S6.1r", norm("ruv_mode"), "skip")
        NULL
    }
    first <- k237Press(S, list(step = "S6.1", name = "normalisation at the defaults", target = paste0(route$prefix, "-norm-run"),
        by = "testid", tabs = tabs, prepare = prepare,
        digest = function(st) if (route$trx) identical(k237Status(st, route$omic, "normalization"), "complete") else
            k237HasState(st, route$omic, paste0(route$prefix, c("_norm_complete", "_ruv_corrected", "_normalized"))),
        dom = function(st) if (route$trx) k233Stepper(st, 4L) else
            k233Noted(st, "Normali[sz]ation (pipeline )?complete", "message")), lane)
    if (route$trx || !identical(first$outcome, "complete")) return(first)
    corr_tabs <- c(tabs, list(list(norm("norm_qc_tabs"), "Correlation Filtering")))
    second <- k233Action(S, list(step = "S6.2", name = "correlation filter at the rendered threshold",
        target = norm("apply_correlation_filter"), tabs = corr_tabs,
        digest = function(st) k237HasState(st, route$omic, paste0(route$prefix, "_correlation_filtered")),
        dom = function(st) k233Noted(st, "Correlation filtering complete", "message")))
    if (!identical(second$outcome, "refused") || !isTRUE(R$permissive)) return(second)
    # The lane's value for this step: Skip instead of Apply (helper-e2e-browser.R:1464-1465), pressed once.
    k233Action(S, list(step = "S6.2r", name = "skip the correlation filter (the lane's value, WL-1's rule)",
        target = norm("skip_correlation_filter"), tabs = corr_tabs,
        digest = function(st) k237HasState(st, route$omic, paste0(route$prefix, c("_correlation_filtered", "_norm_complete"))),
        dom = function(st) any(vapply(st$notes, function(n) identical(n$type, "message"), logical(1)))))
}

k237StepExportOmic <- function(S, R, step, route) {
    # WL 45605975875: the normalisation tab's status and stepper were already complete before the press (both signals
    # at 0 s), so the export's own message is the signal (R/mod_omic_norm.R:650), with no digest signal.
    k233Action(S, list(step = "S7", name = step$name, target = paste0(route$prefix, "-norm-export-session"), by = "testid",
        tabs = list(k237Tab(route, "norm")),
        dom = function(st) k233Noted(st, "Session data exported successfully", "message")))
}

k237StepDaOmic <- function(S, R, step, route) {
    de <- function(x) k237Id(route, "de", x)
    tabs <- list(k237Tab(route, "de"))
    run <- list(name = "run the differential analysis at the defaults", target = paste0(route$prefix, "-da-run"),
        by = "testid", tabs = tabs, texts = de("da_status"), status_out = de("da_status"), refusal = "\\[ERROR\\]|^Error",
        digest = function(st) identical(k237Status(st, route$omic, "differential_analysis"), "complete"),
        dom = function(st) k233Stepper(st, 5L), outputs = de("da_status"))
    if (route$trx) return(k233Action(S, c(list(step = "S8"), run)))
    loaded <- k233Action(S, list(step = "S8.1", name = "load the filtered session",
        target = paste0(route$prefix, "-da-load-session"), by = "testid", tabs = tabs,
        # WL 45605975875: loading commits no state (the head does not change, kit-159/action-scripts.md "Declared
        # states"), so the page's message is the signal (R/mod_omic_da.R:514), with no digest signal.
        dom = function(st) k233Noted(st, "Session loaded successfully", "message")))
    if (!identical(loaded$outcome, "complete")) return(loaded)
    ran <- k233Action(S, c(list(step = "S8.2"), run))
    if (!identical(ran$outcome, "complete")) return(ran)
    # SC-05: the results table opened as a user does; its summary statistics and rows are named outputs; N2 on metab.
    k233Action(S, list(step = "S8.3", name = "open the DA results table", by = "css",
        target = sprintf("#%s a[data-value='DA Results Table']", de("da_results_tabs")), tabs = tabs,
        dom = function(st) isTRUE(st$vis[[1L]]), sels = paste0("#", de("da_results_table"), " tbody tr"),
        outputs = de(c("da_summary_stats", "da_results_table")),
        after = function(S, res, cap) if (identical(route$omic, "metabolomics")) k237NumberN2Omic(S, R, route)))
}

#' N2 on an omics route (SC-05): the DA tab's significant count per contrast, assay "All", against the saved
#' `de_*_long_annot.tsv` tables of that contrast (every assay), counted with the module's rule: fdr_qvalue below the
#' tab's threshold and |logFC| at least its cut-off (R/func_metab_da_orchestration.R:80-81; R/mod_omic_da_panels.R:283-287).
k237NumberN2Omic <- function(S, R, route) {
    de <- function(x) k237Id(route, "de", x)
    threshold <- suppressWarnings(as.numeric(k233InputValue(S, de("da_q_val_thresh"))))
    lfc <- suppressWarnings(as.numeric(k233InputValue(S, de("treat_lfc_cutoff"))))
    options <- as.character(unlist(k233Js(S, sprintf(
        "Array.prototype.map.call((document.getElementById(%s)||{options:[]}).options,function(o){return o.value;})",
        k233Json(de("table_contrast"))))))
    files <- list.files(k233ProjectPaths(R)$base, "^de_.*_metabolites_.*_long_annot\\.tsv$", full.names = TRUE,
        recursive = TRUE)
    files <- files[!duplicated(basename(files))]
    S$n2 <- list()
    for (contrast in options) {
        k233Set(S, "S8.3", de("table_contrast"), contrast)
        Sys.sleep(2)
        text <- as.character(k233Js(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText",
            k233Json(de("da_summary_stats")))) %||% "")
        shown <- suppressWarnings(as.integer(sub("(?s).*Significant \\([^)]*\\): ([0-9]+).*", "\\1", text, perl = TRUE)))
        # RA-03: an exact two-sided match of the contrast's group names against each file's contrast part
        # (de_<mode>_metabolites_<contrast>_long_annot.tsv, R/func_metab_da_results_output.R:166), one file per assay.
        want <- k237ContrastPair(contrast)
        parts <- sub("^de_[^_]+_metabolites_(.*)_long_annot\\.tsv$", "\\1", basename(files))
        mine <- files[vapply(parts, function(x) {
            got <- k237ContrastPair(x)
            !is.null(want) && !is.null(got) && identical(k233Norm(got), k233Norm(want))
        }, logical(1))]
        modes <- sub("^de_([^_]+)_metabolites_.*$", "\\1", basename(mine))
        if (anyDuplicated(modes)) mine <- character()
        k233Observation(S$G, S$label, sprintf("N2 tables for %s (RA-03, RA-05)", contrast), sprintf(
            "%d file(s), one per assay required; read: %s", length(mine),
            paste(substring(k237Slash(mine), nchar(k237Slash(k233ProjectPaths(R)$base)) + 2L), collapse = " + ")))
        saved <- if (!length(mine)) NA_integer_ else sum(vapply(mine, function(f) {
            t <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
            q <- suppressWarnings(as.numeric(t$fdr_qvalue))
            keep <- !is.na(q) & q < threshold
            if (isTRUE(lfc > 0) && "logFC" %in% names(t)) keep <- keep & abs(as.numeric(t$logFC)) >= lfc
            sum(keep)
        }, numeric(1)))
        for (f in mine) file.copy(f, file.path(S$G$out, paste0(S$label, "-", basename(f))), overwrite = TRUE)
        k233Number(S, "N2", contrast, shown, saved, sprintf("da_summary_stats (All) vs %s (fdr_qvalue < %s%s)",
            paste(basename(mine), collapse = " + "), threshold, if (isTRUE(lfc > 0)) sprintf(", |logFC| >= %s", lfc) else ""))
        S$n2[[contrast]] <- list(shown = shown, saved = saved, files = mine, threshold = threshold, lfc = lfc)
    }
    invisible(S$n2)
}

#' The summary as journey.R's S10 does it, on the omic's ids (the shared summary module, R/mod_omic_summary.R); trx
#' runs its RUVIII export first (S10.0). N3 against study_parameters.txt; the report's N4 is a proteomics sentence and
#' is not computed here.
k237StepSummaryOmic <- function(S, R, step, route) {
    sm <- function(x) k237Id(route, "summary", x)
    p <- route$prefix
    tabs <- list(k237Tab(route, "summary"))
    summary_id <- sm("session_summary")
    copy_id <- sm("copy_status")
    actions <- list(
        list(step = "S10.1", name = "Save Workflow Arguments", target = paste0(p, "-summary-save-args"),
            prepare = function(S) {
                k233Set(S, "S10.1", sm("experiment_label"), R$label)
                k233Set(S, "S10.1", sm("description"), R$description)
                NULL
            },
            dom = function(st) grepl("Parameters saved [OK]", k233Text(st, summary_id), fixed = TRUE),
            outputs = summary_id, after = function(S, res, cap) k237NumberN3(S, R, route)),
        list(step = "S10.2", name = "Copy to Publication Directory", target = paste0(p, "-summary-copy-publication"),
            dom = function(st) grepl("Files copied", k233Text(st, summary_id), fixed = TRUE),
            outputs = c(copy_id, summary_id)),
        list(step = "S10.3", name = "Generate Report", target = paste0(p, "-summary-generate-report"),
            sels = sprintf("[data-testid='%s-summary-download-report']", p),
            dom = function(st) isTRUE(st$vis[[1L]]) || k233Noted(st, "Report generated successfully", "message"),
            outputs = summary_id),
        list(step = "S10.4", name = "Download Report", custom = function(S) k237DownloadReport(S, R, route)),
        list(step = "S10.5", name = "Export Session State", target = paste0(p, "-summary-export-state"),
            dom = function(st) k233Noted(st, "Session state exported to:", "message")))
    if (route$trx) {
        actions <- c(list(list(step = "S10.0", name = "RUVIII export", target = "trx-summary-ruviii-export",
            texts = sm("ruviii_status"),
            dom = function(st) grepl("M = group", k233Text(st, sm("ruviii_status")), fixed = TRUE),
            outputs = sm("ruviii_status"))), actions)
    }
    res <- list(outcome = "complete")
    for (a in actions) {
        a$by <- a$by %||% "testid"
        a$tabs <- tabs
        a$texts <- unique(c(a$texts, summary_id, copy_id))
        res <- k233Action(S, a)
        if (!identical(res$outcome, "complete")) return(res)
    }
    res
}

k237DownloadReport <- function(S, R, route) {
    p <- k233Press(S, "S10.4", paste0(route$prefix, "-summary-download-report"), "testid", click = FALSE)
    k233LogAction(S, "S10.4", paste0(route$prefix, "-summary-download-report"), "download (the link's file fetched)")
    if (!isTRUE(p$ok)) return(list(outcome = "refused", message = paste("not reachable as a user:", p$reason)))
    path <- tryCatch(S$G$helpers$e2e_get_download(S$driver, k237Id(route, "summary", "download_report"),
        artifact_dir = S$G$work, filename = paste0(S$label, "-report-download")), error = function(e) conditionMessage(e))
    if (!file.exists(path) || file.size(path) == 0) {
        return(list(outcome = "refused", message = paste("the download gave no file:", path)))
    }
    note <- sprintf("report %d bytes, sha256 %s", file.size(path), k233Sha(path))
    list(outcome = "complete", message = "", dom_t = NA, dig_t = NA, agreement = "dom_only", note = note)
}

k237NumberN3 <- function(S, R, route) {
    sm <- function(x) k237Id(route, "summary", x)
    file <- file.path(k233ProjectPaths(R)$source, "study_parameters.txt")
    lines <- k233ReadLog(file)
    field <- function(label) {
        hit <- grep(paste0("^", label, ": "), lines, value = TRUE)
        if (length(hit)) sub(paste0("^", label, ": "), "", hit[[1L]]) else NA_character_
    }
    if (file.exists(file)) file.copy(file, file.path(S$G$out, paste0(S$label, "-study_parameters.txt")), overwrite = TRUE)
    where <- sprintf("summary inputs vs scripts/%s/study_parameters.txt", route$omic)
    k233Number(S, "N3", "label", k233Clean(k233InputValue(S, sm("experiment_label"))), field("Workflow Name"), where)
    k233Number(S, "N3", "description", k233Clean(k233InputValue(S, sm("description"))), field("Description"), where)
}
