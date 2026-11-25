# KIT-237 (r1.2) 5.2: the label-free (MaxQuant) route for KIT-233's journey machinery, lane `prot_lfq`, fixture
# tests/testdata/e2e/prot_lfq_maxquant (MaxQuant proteinGroups.txt, 400 proteins plus 5 CON__ and 5 REV__ rows the
# reader drops, 10 samples WT_1..WT_5 and KO_1..KO_5, proteins.fasta, design.tsv; manifest.json:134-175, SHA-256
# :172-174). Sourced by journey-win.R after journey.R. Steps S3, S6, S7, S8 and S10 are journey.R's own; S1, S2, S5 and
# S9 are this route's variants (5.2). r1.3: the chooser's selections are counted for WL-0u (SC-02); the S9 download
# uses journey.R's click-free form (journey_commit at or after 628545be). The ids, named outputs and named files below are fixed before WL and not changed
# after WL passes.

K237_ROUTES <- if (exists("K237_ROUTES", inherits = FALSE)) K237_ROUTES else list()
K237_RUNS <- if (exists("K237_RUNS", inherits = FALSE)) K237_RUNS else list()

# Named files per action (5.2, "Files each step must leave"): regular expressions on the path relative to the project
# base (<project folder>/<label>), matched anywhere below it. By reading at dd769671:
#   S6.2: the normalisation figures (R/mod_prot_norm_qc_generation_helpers.R:98-215, saveArtifactFn names);
#   S7:   scripts/proteomics/filtered_session_summary.txt (journey.R N1 reads it);
#   S8.2: the per-contrast DA tables (*_long_annot.tsv in results/proteomics/da_proteins; journey.R N2 reads them);
#   S10.1: scripts/proteomics/study_parameters.txt (N3);
#   S10.2: Publication_tables/DA_results_proteomics.xlsx (copyToResultsSummary,
#          R/func_general_filemgmt_results_summary_helpers.R:643, :712);
#   S10.3: the report template and the rendered report (tests/testthat/test-e2e-proteomics-tmt-lfq.R:99, :108-115);
#   S10.5: the session state export (test-e2e-proteomics-tmt-lfq.R:117-123);
#   S10.5 also: the lane's `expected_exports` (manifest.json:146-150), matched by name anywhere. By reading no code
#          under R/ writes these three names for proteomics (no match for them in R/ at dd769671), so a kind-3 finding
#          is expected for each; known.tsv lists them.
K237_LFQ_FILES <- list(
    S6.2 = paste0("(^|/)", c("pre_norm_pca", "pre_norm_rle", "pre_norm_density", "pre_norm_correlation", "post_norm_pca",
        "post_norm_rle", "post_norm_density", "post_norm_correlation"), "\\.png$"),
    S7 = "^scripts/proteomics/filtered_session_summary\\.txt$",
    S8.2 = "^results/proteomics/da_proteins/.*_long_annot\\.tsv$",
    S10.1 = "^scripts/proteomics/study_parameters\\.txt$",
    S10.2 = "(^|/)Publication_tables/DA_results_proteomics\\.xlsx$",
    S10.3 = c("^scripts/proteomics/LFQ_report\\.rmd$", "(^|/)LFQ_report_proteomics_[^/]*\\.(docx|html|pdf)$"),
    # SC-09: the lane's placeholder `expected_exports` are not named (no product output by those names; known.tsv K2).
    S10.5 = "^scripts/proteomics/session_state_[^/]*\\.RDS$"
)

K237_ROUTES$prot_lfq <- list(
    id = "LFQ",
    omic = "proteomics",
    tabs_id = "proteomics_workflow-workflow_tabs",
    steps = list(
        list(id = "S1", name = "home and project (the rendered default base directory kept)", fn = "k237StepHome",
            tile = "tile-proteomics", omic = "proteomics", tab_testid = "prot-tab-setup"),
        list(id = "S2", name = "import (format left at Auto-detect)", fn = "k237StepImportLfq"),
        list(id = "S3", name = "design", fn = "k233StepDesign"),
        # The label-free route has one QC group, Protein QC (R/utils_omic_gui_config.R:2234-2259: tabset qc_tabs_tmt;
        # the peptide tabs carry the DIA routes only). S4 is therefore absent; the rendered group titles are recorded.
        list(id = "S5", name = "protein QC", fn = "k237StepQcLfq", group = "protein_qc", title = "Protein QC",
            tabset = "protein_qc-protein_filter_tabs"),
        list(id = "S6", name = "normalisation", fn = "k233StepNormalisation"),
        list(id = "S7", name = "export", fn = "k233StepExport"),
        list(id = "S8", name = "differential abundance", fn = "k233StepDa"),
        list(id = "S9", name = "enrichment, then the results download", fn = "k237StepEnrichmentLfq",
            continue_on_refusal = TRUE),
        list(id = "S10", name = "summary", fn = "k233StepSummary")
    ),
    qc_tabset = "qc_tabs_tmt",
    qc = list(
        protein_qc = list(
            list(id = "cleanup", title = "Accession Cleanup", button = "apply_accession_cleanup",
                state = "protein_accession_cleaned", outputs = c("accession_cleanup_results", "accession_cleanup_plot")),
            # Permissive values: the lane's own (tests/testthat/helper-e2e-browser.R:1315-1324), used only by WL-1's
            # rule (a default refused on 10 samples).
            list(id = "intensity_filter", title = "Protein Intensity Filter",
                button = "apply_protein_intensity_filter", state = "protein_intensity_filtered",
                outputs = c("protein_intensity_filter_results", "protein_intensity_filter_plot"),
                permissive = list(use_strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1,
                    proteins_intensity_cutoff_percentile = 0.1)),
            list(id = "duplicate_removal", title = "Duplicate Removal", button = "apply_duplicate_removal",
                state = "duplicates_removed", outputs = c("duplicate_removal_results", "duplicate_removal_plot")),
            list(id = "replicate_filter", title = "Protein Replicate Filter",
                button = "apply_protein_replicate_filter", state = "protein_replicate_filtered",
                outputs = c("protein_replicate_filter_results", "protein_replicate_filter_plot"),
                permissive = list(parallel_cores = 1))
        )
    )
)

K237_RUNS$prot_lfq <- function(G, dir, label) {
    design <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE)
    groups <- unique(design$group)
    list(dataset = "prot_lfq_maxquant", omic = "proteomics", report = file.path(dir, "proteinGroups.txt"),
        fasta = file.path(dir, "proteins.fasta"), groups = groups,
        members = lapply(stats::setNames(groups, groups), function(g) design$sample[design$group == g]),
        contrasts = list(c("KO", "WT")), label = label,
        description = "KIT-237 label-free journey on prot_lfq_maxquant", permissive = TRUE, project_dir = NULL,
        files = K237_LFQ_FILES, last_step = "S10.5", summary_step = "S10", fixture_dir = dir)
}

# --- route steps ------------------------------------------------------------------------------------------------

#' S1 for every route: the tile, Start, the Experiment Setup dialog with the label typed and the rendered default base
#' directory kept (R/app_server.R:164, `normalizePath("~", winslash = "/")`), unless the run names a folder (WL-2).
#' The folder used is recorded with its length and kept in `R$live` for the later steps.
k237StepHome <- function(S, R, step, route) {
    prepare <- function(S) {
        for (t in c(step$tile, "btn-start-analysis")) {
            p <- k233Press(S, "S1", t, "testid")
            if (!isTRUE(p$ok)) return(list(refuse = paste("not reachable as a user:", p$reason)))
            k233Idle(S, 10)
        }
        if (!k233WaitFor(S, "document.getElementById('experiment_label') !== null", 30)) {
            return(list(refuse = "the Experiment Setup dialog did not open within 30 s"))
        }
        k233Set(S, "S1", "experiment_label", R$label)
        rendered <- k233Clean(k233InputValue(S, "project_dir"))
        k233Observation(S$G, S$label, "S1 rendered base directory", sprintf("%s (%d characters)", rendered,
            nchar(rendered)))
        if (!is.null(R$project_dir_override)) {
            k233Set(S, "S1", "project_dir", R$project_dir_override)
            R$live$project_dir <- R$project_dir_override
        } else {
            R$live$project_dir <- rendered
        }
        NULL
    }
    omic <- step$omic
    k233Action(S, list(step = "S1", name = step$name, target = "btn-confirm-setup", by = "testid", prepare = prepare,
        tab_note = "none (home page and the Experiment Setup dialog)", sels = sprintf("[data-testid='%s']", step$tab_testid),
        dom = function(st) isTRUE(st$vis[[1L]]), digest = function(st) omic %in% unlist(st$dg$project_dir_keys),
        modals = "Experiment Setup"))
}

#' S2, label-free: proteinGroups.txt and the FASTA (test mode: the upload boxes; user mode: the shinyFiles chooser,
#' 5.2), the format left at Auto-detect and the detected format recorded (MaxQuant expected; otherwise a kind-2 finding
#' and MaxQuant chosen by hand).
k237StepImportLfq <- function(S, R, step, route) {
    ns <- paste0(K233_NS, "setup_import-")
    prepare <- function(S) {
        if (S$test_mode) {
            up <- k233Upload(S, "S2", paste0(ns, "search_results_standard"), R$report)
            if (!up$ok) return(list(refuse = paste("search results upload:", up$text)))
        } else {
            chosen <- k233ChooseFile(S, "S2", paste0(ns, "search_results"), k237UserInput(R, R$report))
            S$chooser <- c(S$chooser, isTRUE(chosen$ok))
            k233Observation(S$G, S$label, "S2 search results chooser", chosen$how)
            if (!chosen$ok) return(list(refuse = paste("not reachable as a user: search results chooser:", chosen$how)))
        }
        k233WaitFor(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText.trim().length > 0",
            k233Json(paste0(ns, "format_detection"))), 60)
        detected <- as.character(k233Js(S, sprintf("(document.getElementById(%s)||{innerText:''}).innerText.trim()",
            k233Json(paste0(ns, "format_detection")))) %||% "")
        k233Observation(S$G, S$label, "S2 detected format", detected)
        if (!grepl("MaxQuant", detected, ignore.case = TRUE)) {
            k233Finding(S, "S2", 2L, paste("detected format is not MaxQuant; MaxQuant chosen by hand. Shown:", detected))
            k233Set(S, "S2", paste0(ns, "format_override"), "maxquant")
        }
        if (S$test_mode) {
            up <- k233Upload(S, "S2", paste0(ns, "fasta_file_standard"), R$fasta)
            if (!up$ok) return(list(refuse = paste("FASTA upload:", up$text)))
        } else {
            chosen <- k233ChooseFile(S, "S2", paste0(ns, "fasta_file"), k237UserInput(R, R$fasta))
            S$chooser <- c(S$chooser, isTRUE(chosen$ok))
            k233Observation(S$G, S$label, "S2 FASTA chooser", chosen$how)
            if (!chosen$ok) return(list(refuse = paste("not reachable as a user: FASTA chooser:", chosen$how)))
        }
        NULL
    }
    on_modal <- function(S, modal) {
        p <- k233Press(S, "S2", paste0(ns, "confirm_organism"))
        k233Observation(S$G, S$label, "S2 organism modal", sprintf("%s | answered: %s", modal$text, isTRUE(p$ok)))
    }
    k233Action(S, list(step = "S2", name = step$name, target = "prot-import-process", by = "testid",
        tabs = list(k233Workflow("setup")), prepare = prepare,
        digest = function(st) identical(k233DigestStatus(st, "setup_import"), "complete"),
        dom = function(st) k233Stepper(st, 1L), modals = "Select Primary Organism", on_modal = on_modal,
        outputs = paste0(ns, c("format_detection", "import_status", "data_summary"))))
}

#' User mode: the input copied to <the kept base directory>/k237-input/, where the chooser reaches it (5.2).
k237UserInput <- function(R, path) {
    home <- file.path(R$live$project_dir %||% normalizePath("~", winslash = "/"), "k237-input")
    dir.create(home, recursive = TRUE, showWarnings = FALSE)
    file.copy(path, home, overwrite = TRUE)
    k237Slash(file.path(home, basename(path)))
}

#' S5, label-free: the rendered QC group titles recorded, then journey.R's QC step on the Protein QC group.
k237StepQcLfq <- function(S, R, step, route) {
    k233ActivateTab(S, k233Workflow("qc")[[1L]], "qc")
    titles <- as.character(unlist(k233Js(S, sprintf(K233_TITLES_JS, k233Json(paste0(K233_NS, "quality_control-",
        route$qc_tabset))))))
    k233Observation(S$G, S$label, "S5 rendered QC groups", paste(titles, collapse = " > "))
    k233StepQc(S, R, step, route)
}

#' S9: journey.R's enrichment step (test mode: the substitute; user mode offline: the expected refusal), then the
#' "Download All Results" button pressed once in either case and the download recorded (bytes; whether it opens as a
#' zip; kind 3 if not). In WL, before KIT-241 merges, the kind-3 finding is expected (known.tsv, KIT-241).
k237StepEnrichmentLfq <- function(S, R, step, route) {
    res <- k233StepEnrichment(S, R, step, route)
    if (res$outcome %in% c("complete", "refused")) {
        k233Action(S, list(step = "S9.d", name = "Download All Results (enrichment)",
            tabs = list(k233Workflow("enrichment")), custom = function(S) k237DownloadEnrichment(S, R)))
    }
    res
}

k237DownloadEnrichment <- function(S, R) {
    id <- paste0(K233_NS, "enrichment_analysis-download_enrichment_results")
    # As journey.R's report download since 628545be: the link's reachability checked without a click, then the link's
    # file fetched (one request: the user's press); needs a journey_commit at or after 628545be.
    p <- k233Press(S, "S9.d", id, "id", click = FALSE)
    k233LogAction(S, "S9.d", id, "download (the link's file fetched)")
    if (!isTRUE(p$ok)) {
        S$enrich_download <- list(bytes = 0, opens = FALSE, note = paste("not reachable:", p$reason))
        return(list(outcome = "refused", message = paste("not reachable as a user:", p$reason)))
    }
    error <- ""
    path <- tryCatch(S$G$helpers$e2e_get_download(S$driver, id, artifact_dir = S$G$work,
        filename = paste0(S$label, "-enrichment-download.zip")), error = function(e) {
            error <<- conditionMessage(e)
            ""
        })
    bytes <- if (nzchar(path) && file.exists(path)) file.size(path) else 0
    entries <- if (bytes > 0) tryCatch(utils::unzip(path, list = TRUE)$Name, error = function(e) NULL) else NULL
    opens <- !is.null(entries)
    note <- sprintf("download %s bytes, sha256 %s; opens as a zip: %s; entries: %s%s", bytes,
        if (bytes > 0) k233Sha(path) else "-", opens, paste(entries, collapse = ", "),
        if (nzchar(error)) paste("; the download call said:", error) else "")
    S$enrich_download <- list(bytes = bytes, opens = opens, note = note)
    if (!opens) k233Finding(S, "S9.d", 3L, paste("the enrichment results download is not a zip that opens:", note))
    list(outcome = "complete", message = "", dom_t = NA, dig_t = NA, agreement = "dom_only", note = note)
}
