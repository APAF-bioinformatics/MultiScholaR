# KIT-106 W1: the DA tab through its public entry points, written before the shared DA frame
# (R/mod_omic_da.R, R/mod_omic_da_panels.R; W3 to W5) and the core moves (W2) exist, and frozen on the
# start commit (b75561de; docs/qa/workflow-execution/kit/kit-106/w1/test-first.txt).
#   (a) Identity: the four DA UIs' HTML digests (the `checked` attribute of `volcano_interactive` masked:
#       D3 pins it separately) and their namespaced input and output ids, test ids, tab titles, labels,
#       help texts, choice lists and defaults (`volcano_interactive`'s left out).
#   (b) Server captures, in process, through `mod_{lipid,metab,prot}_da_server` and the trx mount
#       (`mod_lipid_da_server(profile = trxOmicProfile())`) in `testServer`, on the MCI fixtures: lipid
#       MCI-022 (helper-module-ci-lipid-da.R), metab MCI-016, trx MCI-030/031
#       (`module_ci_trx_project(module_ci_trx_corpus())`, the normalisation run first on the same
#       session) and proteomics MCI-009 (helper-module-ci-prot-da.R; its object on a DIA owner, opened
#       through the tab selection, and its session payload as the legacy export that owner adopts on the
#       load). Cases: the load (success, a missing file, an assay mismatch or reorder); the run
#       (success, no data, no contrasts, an invalid preflight); every selector after completion
#       (contrast, assay, the heatmap's controls with each tree-cut method, the table's significance and
#       row limit); Save Heatmap; the download. Proteomics' renders are also read with the interactive
#       volcano box unticked and ticked (its render reads no box today), with its label inputs and
#       dendrogram switch.
#   The buses are the MCI fixtures': an environment for lipid, metab and trx (an output that reads only
#   the bus, such as `contrasts_display` or trx's `formula_display`, is computed once per session), and
#   reactive values for proteomics. The entry lanes hold the reactive bus.
#   (b-managed) Lipid, metab and trx with `omicsManagedChainMocks`: the worker's operations in order,
#       the notifications, `tab_status` and the rows' digest. The volcano page is never requested.
#   (c) The per-omic differences (design D7) recorded as today's behaviour: the refusal order, the
#       assay check's rule and text, the table filter on "Combined", a long table without
#       `friendly_name`, and (its own block: D2 changes it) the value each server returns.
# Captured after each step: the registered output names and each value as the entry suite reads it
# (`guiEntryOutputValue`; `volcano_glimma` is not read for lipid, metab and trx), each plot and table
# the renders built (the kernel's name, its scalar arguments and a digest of its data), the
# notifications (type, duration, id and message; each removal), `tab_status`, the history and the
# state digests, the tab's reactive values (names and a reference description, never content), the
# rows the tab holds (a digest of each long table; proteomics' manifest without its store source),
# the files under the project (table digests) and the selector updates.
# MULTISCHOLAR_KIT106_W1_DUMP=<file> writes the captures as R code (how the frozen values were made;
# two dumps in separate processes were byte-identical before they were pasted).

omicDaW1Sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)

#' Where MULTISCHOLAR_KIT106_W1_DUMP points, or NULL.
omicDaW1DumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT106_W1_DUMP")
    if (nzchar(path)) path else NULL
}

#' R code for a capture (lists, named vectors and scalars): a list's elements one per line, a long
#' vector's packed into lines, indented (test-trx-norm-module.R's `trxNormW1Code`).
omicDaW1Code <- function(x, indent = 0L) {
    pad <- strrep(" ", indent + 4L)
    label <- function(names, i) {
        if (is.null(names) || !nzchar(names[[i]])) return("")
        paste0(if (make.names(names[[i]]) == names[[i]]) names[[i]] else paste0("`", names[[i]], "`"), " = ")
    }
    if (is.list(x)) {
        if (length(x) == 0L) return("list()")
        items <- vapply(seq_along(x), function(i) paste0(pad, label(names(x), i), omicDaW1Code(x[[i]], indent + 4L)), "")
        return(paste0("list(\n", paste(items, collapse = ",\n"), ")"))
    }
    if (is.null(x)) return("NULL")
    scalar <- function(value) deparse(value, control = c("keepNA", "keepInteger"), width.cutoff = 500L)
    if (length(x) == 0L) return(if (is.null(names(x))) scalar(x) else sprintf("stats::setNames(%s, character())", scalar(unname(x))))
    if (length(x) == 1L && is.null(names(x))) return(scalar(x))
    items <- paste0(vapply(seq_along(x), function(i) paste0(label(names(x), i), scalar(unname(x[[i]]))), ""), ",")
    items[[length(items)]] <- sub(",$", ")", items[[length(items)]])
    if (sum(nchar(items) + 1L) + indent < 90L) return(paste0("c(", paste(items, collapse = " ")))
    lines <- character()
    for (item in items) {
        last <- length(lines)
        if (last && nchar(lines[[last]]) + nchar(item) < 110L) {
            lines[[last]] <- paste(lines[[last]], item)
        } else {
            lines <- c(lines, paste0(pad, item))
        }
    }
    paste0("c(\n", paste(lines, collapse = "\n"))
}

#' Append `name <- <value>` to the dump file.
omicDaW1Dump <- function(name, value) {
    cat(sprintf("%s <- %s\n\n", name, omicDaW1Code(value)), file = omicDaW1DumpPath(), append = TRUE)
}

# ------------------------------------------------------------------------------------------------
# (a) Identity: the four DA UIs, frozen at b75561de.
# ------------------------------------------------------------------------------------------------

#' The four DA UIs, as the shell mounts them (mod_omic_workflow.R:221-246).
omicDaW1Uis <- function() {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    list(lipid = mod_lipid_da_ui("de"), metab = mod_metab_da_ui("de"),
        trx = mod_lipid_da_ui("de", profile = trxOmicProfile()), prot = mod_prot_da_ui("differential_expression"))
}

#' The rendered HTML with shiny's random tabset ids masked, and the `checked` attribute of the
#' interactive volcano's box removed (D3: its default changes and is pinned by W3 to W5).
omicDaW1Html <- function(ui) {
    html <- as.character(htmltools::renderTags(ui)$html)
    html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"ID\"", html)
    html <- gsub("tab-[0-9]+-([0-9]+)", "tab-ID-\\1", html)
    gsub("(<input id=\"[^\"]*-volcano_interactive\"[^>]*?) checked=\"checked\"", "\\1", html)
}

#' One DA UI's HTML digest and its parts: inputs (every namespaced `id` that is not an output, a
#' label or a tab pane), outputs, `data-testid` values, tab titles in order, each label's text by its
#' input, the help texts in order, each select's choices (`value=text`, the selected one marked) and
#' the inputs' defaults (the textarea's value too; the interactive volcano's box left out).
omicDaW1UiParts <- function(ui) {
    html <- omicDaW1Html(ui)
    tags <- regmatches(html, gregexpr("<[a-zA-Z][^>]*>", html))[[1L]]
    attr_of <- function(name) vapply(tags, function(tag) {
        m <- regmatches(tag, regexec(sprintf("[[:space:]]%s=\"([^\"]*)\"", name), tag))[[1L]]
        if (length(m)) m[[2L]] else NA_character_
    }, "", USE.NAMES = FALSE)
    ids <- attr_of("id")
    classes <- attr_of("class")
    is_output <- !is.na(ids) & grepl("(^| )(shiny-[a-z]+-output|html-widget-output)( |$)", classes)
    inputs <- ids[!is.na(ids) & !is_output & !grepl("-label$", ids) & !startsWith(ids, "tab-ID-")]
    labels <- regmatches(html, gregexpr("<label[^>]*for=\"[^\"]+\"[^>]*>[^<]*</label>", html))[[1L]]
    label_text <- stats::setNames(sub("^<label[^>]*>([^<]*)</label>$", "\\1", labels),
        sub("^.*for=\"([^\"]+)\".*$", "\\1", labels))
    boxes <- regmatches(html, gregexpr("<input id=\"[^\"]+\" type=\"checkbox\"[^>]*/>[[:space:]]*<span>[^<]*</span>", html))[[1L]]
    box_text <- stats::setNames(sub("^.*<span>([^<]*)</span>$", "\\1", boxes), sub("^<input id=\"([^\"]+)\".*$", "\\1", boxes))
    help <- regmatches(html, gregexpr("<span class=\"help-block\">[^<]*</span>", html))[[1L]]
    selects <- regmatches(html, gregexpr("(?s)<select id=\"[^\"]+\"[^>]*>.*?</select>", html, perl = TRUE))[[1L]]
    choices <- lapply(stats::setNames(selects, sub("^<select id=\"([^\"]+)\".*$", "\\1", selects)), function(select) {
        options <- regmatches(select, gregexpr("<option value=\"[^\"]*\"( selected)?>[^<]*</option>", select))[[1L]]
        sub("^<option value=\"([^\"]*)\"( selected)?>([^<]*)</option>$", "\\1=\\3\\2", options)
    })
    defaults <- omicsUiInputDefaults(ui)
    areas <- regmatches(html, gregexpr("<textarea id=\"[^\"]+\"[^>]*>[^<]*</textarea>", html))[[1L]]
    defaults[sub("^<textarea id=\"([^\"]+)\".*$", "\\1", areas)] <- as.list(sub("^.*>([^<]*)</textarea>$", "\\1", areas))
    defaults <- defaults[!endsWith(names(defaults), "-volcano_interactive")]
    defaults <- vapply(defaults[order(names(defaults))], function(value) paste(format(value), collapse = ","), "")
    list(html = digest::digest(html, algo = "sha256", serialize = FALSE), inputs = sort(unique(inputs)),
        outputs = sort(unique(ids[is_output])), testids = sort(unique(stats::na.omit(attr_of("data-testid")))),
        tabs = attr_of("data-value")[grepl("data-toggle=\"tab\"", tags)],
        labels = c(label_text, box_text)[order(names(c(label_text, box_text)))],
        help = sub("^<span class=\"help-block\">([^<]*)</span>$", "\\1", help),
        choices = if (length(choices)) choices[order(names(choices))] else list(),
        defaults = defaults)
}

# OMIC_DA_FROZEN_UI_START
omicDaW1FrozenUi <- list(
    lipid = list(
        html = "9b8de0acc639a0aab8889b190ef7e66c7a3e395a4fa5b99cbb1c5b0944ec4c76",
        inputs = c(
            "de-da_q_val_thresh", "de-da_results_tabs", "de-download_da_results", "de-formula_string",
            "de-heatmap_assay", "de-heatmap_cluster_method", "de-heatmap_clustering",
            "de-heatmap_color_scheme", "de-heatmap_container", "de-heatmap_contrast", "de-heatmap_cut_height",
            "de-heatmap_distance_method", "de-heatmap_min_cluster_size", "de-heatmap_n_clusters",
            "de-heatmap_scaling", "de-heatmap_show_labels", "de-heatmap_top_n", "de-heatmap_tree_cut_method",
            "de-load_filtered_session", "de-run_da_analysis", "de-save_heatmap", "de-table_assay",
            "de-table_contrast", "de-table_max_rows", "de-table_significance", "de-treat_lfc_cutoff",
            "de-volcano_assay", "de-volcano_contrast", "de-volcano_glimma_container",
            "de-volcano_interactive"),
        outputs = c(
            "de-cluster_summary", "de-contrasts_display", "de-da_results_table", "de-da_status",
            "de-da_summary_stats", "de-heatmap_manual_save_warning", "de-heatmap_plot", "de-volcano_glimma",
            "de-volcano_static"),
        testids = c(
            "lipid-da-download", "lipid-da-load-session", "lipid-da-results-tabs", "lipid-da-run"),
        tabs = c("Volcano Plot", "Heatmap", "DA Results Table"),
        labels = c(
            `de-da_q_val_thresh` = "Q-value Threshold:", `de-formula_string` = "Model Formula:",
            `de-heatmap_assay` = "Select Assay:", `de-heatmap_cluster_method` = "Clustering Method:",
            `de-heatmap_clustering` = "Apply Clustering:", `de-heatmap_color_scheme` = "Color Scheme:",
            `de-heatmap_contrast` = "Select Contrast:", `de-heatmap_cut_height` = "Cut Height:",
            `de-heatmap_distance_method` = "Distance Metric:",
            `de-heatmap_min_cluster_size` = "Min Cluster Size:",
            `de-heatmap_n_clusters` = "Number of Clusters:", `de-heatmap_scaling` = "Data Scaling:",
            `de-heatmap_show_labels` = "Show Lipid Labels", `de-heatmap_top_n` = "Top N Lipids:",
            `de-heatmap_tree_cut_method` = "Tree Cutting:", `de-table_assay` = "Select Assay:",
            `de-table_contrast` = "Select Contrast:", `de-table_max_rows` = "Max Rows:",
            `de-table_significance` = "Show:", `de-treat_lfc_cutoff` = "Log Fold-Change Cutoff:",
            `de-volcano_assay` = "Select Assay:", `de-volcano_contrast` = "Select Contrast:",
            `de-volcano_interactive` = "Interactive Plot (Glimma)"),
        help = c(
            "Formula fetched from S4 object @args if available", "Contrasts defined in design matrix",
            "Load previously exported filtered data for DA analysis"),
        choices = list(
            `de-heatmap_assay` = "Combined=Combined selected",
            `de-heatmap_cluster_method` = c(
                "ward.D2=Ward (minimum variance) selected", "ward.D=Ward (original)",
                "complete=Complete linkage", "single=Single linkage", "average=Average linkage",
                "mcquitty=McQuitty (WPGMA)"),
            `de-heatmap_clustering` = c(
                "both=Both rows &amp; columns selected", "row=Rows only", "column=Columns only", "none=None"),
            `de-heatmap_color_scheme` = c(
                "RdBu=Red-Blue selected", "RdYlBu=Red-Yellow-Blue", "coolwarm=Blue-White-Red",
                "viridis=Viridis", "plasma=Plasma", "inferno=Inferno"),
            `de-heatmap_contrast` = character(0),
            `de-heatmap_distance_method` = c(
                "euclidean=Euclidean selected", "manhattan=Manhattan", "pearson=Pearson correlation",
                "spearman=Spearman correlation", "maximum=Maximum"),
            `de-heatmap_scaling` = c(
                "row=Row (lipid) scaling selected", "column=Column (sample) scaling", "both=Both",
                "none=None"),
            `de-heatmap_tree_cut_method` = c(
                "none=None selected", "k_clusters=K Clusters", "height_cutoff=Height Cutoff",
                "dynamic=Dynamic"),
            `de-table_assay` = "All=All selected",
            `de-table_contrast` = character(0),
            `de-table_significance` = c(
                "all=All Results", "significant=Significant Only selected", "up=Up-regulated",
                "down=Down-regulated"),
            `de-volcano_assay` = "Combined=Combined selected",
            `de-volcano_contrast` = character(0)),
        defaults = c(
            `de-da_q_val_thresh` = "0.05", `de-formula_string` = "~ 0 + group",
            `de-heatmap_assay` = "Combined", `de-heatmap_cluster_method` = "ward.D2",
            `de-heatmap_clustering` = "both", `de-heatmap_color_scheme` = "RdBu",
            `de-heatmap_cut_height` = "10", `de-heatmap_distance_method` = "euclidean",
            `de-heatmap_min_cluster_size` = "5", `de-heatmap_n_clusters` = "4", `de-heatmap_scaling` = "row",
            `de-heatmap_show_labels` = "FALSE", `de-heatmap_top_n` = "50",
            `de-heatmap_tree_cut_method` = "none", `de-table_assay` = "All", `de-table_max_rows` = "1000",
            `de-table_significance` = "significant", `de-treat_lfc_cutoff` = "0",
            `de-volcano_assay` = "Combined")),
    metab = list(
        html = "93ae770af34c04eaf88782a58645b08fbfe6927e9d020d708b4f58e3ee99dfa5",
        inputs = c(
            "de-da_q_val_thresh", "de-da_results_tabs", "de-download_da_results", "de-formula_string",
            "de-heatmap_assay", "de-heatmap_cluster_method", "de-heatmap_clustering",
            "de-heatmap_color_scheme", "de-heatmap_container", "de-heatmap_contrast", "de-heatmap_cut_height",
            "de-heatmap_distance_method", "de-heatmap_min_cluster_size", "de-heatmap_n_clusters",
            "de-heatmap_scaling", "de-heatmap_show_labels", "de-heatmap_top_n", "de-heatmap_tree_cut_method",
            "de-load_filtered_session", "de-run_da_analysis", "de-save_heatmap", "de-table_assay",
            "de-table_contrast", "de-table_max_rows", "de-table_significance", "de-treat_lfc_cutoff",
            "de-volcano_assay", "de-volcano_contrast", "de-volcano_glimma_container",
            "de-volcano_interactive"),
        outputs = c(
            "de-cluster_summary", "de-contrasts_display", "de-da_results_table", "de-da_status",
            "de-da_summary_stats", "de-heatmap_manual_save_warning", "de-heatmap_plot", "de-volcano_glimma",
            "de-volcano_static"),
        testids = c(
            "metab-da-download", "metab-da-load-session", "metab-da-results-tabs", "metab-da-run"),
        tabs = c("Volcano Plot", "Heatmap", "DA Results Table"),
        labels = c(
            `de-da_q_val_thresh` = "Q-value Threshold:", `de-formula_string` = "Model Formula:",
            `de-heatmap_assay` = "Select Assay:", `de-heatmap_cluster_method` = "Clustering Method:",
            `de-heatmap_clustering` = "Apply Clustering:", `de-heatmap_color_scheme` = "Color Scheme:",
            `de-heatmap_contrast` = "Select Contrast:", `de-heatmap_cut_height` = "Cut Height:",
            `de-heatmap_distance_method` = "Distance Metric:",
            `de-heatmap_min_cluster_size` = "Min Cluster Size:",
            `de-heatmap_n_clusters` = "Number of Clusters:", `de-heatmap_scaling` = "Data Scaling:",
            `de-heatmap_show_labels` = "Show Metabolite Labels", `de-heatmap_top_n` = "Top N Metabolites:",
            `de-heatmap_tree_cut_method` = "Tree Cutting:", `de-table_assay` = "Select Assay:",
            `de-table_contrast` = "Select Contrast:", `de-table_max_rows` = "Max Rows:",
            `de-table_significance` = "Show:", `de-treat_lfc_cutoff` = "Log Fold-Change Cutoff:",
            `de-volcano_assay` = "Select Assay:", `de-volcano_contrast` = "Select Contrast:",
            `de-volcano_interactive` = "Interactive Plot (Glimma)"),
        help = c(
            "Formula fetched from S4 object @args if available", "Contrasts defined in design matrix",
            "Load previously exported filtered data for DA analysis"),
        choices = list(
            `de-heatmap_assay` = "Combined=Combined selected",
            `de-heatmap_cluster_method` = c(
                "ward.D2=Ward (minimum variance) selected", "ward.D=Ward (original)",
                "complete=Complete linkage", "single=Single linkage", "average=Average linkage",
                "mcquitty=McQuitty (WPGMA)"),
            `de-heatmap_clustering` = c(
                "both=Both rows &amp; columns selected", "row=Rows only", "column=Columns only", "none=None"),
            `de-heatmap_color_scheme` = c(
                "RdBu=Red-Blue selected", "RdYlBu=Red-Yellow-Blue", "coolwarm=Blue-White-Red",
                "viridis=Viridis", "plasma=Plasma", "inferno=Inferno"),
            `de-heatmap_contrast` = character(0),
            `de-heatmap_distance_method` = c(
                "euclidean=Euclidean selected", "manhattan=Manhattan", "pearson=Pearson correlation",
                "spearman=Spearman correlation", "maximum=Maximum"),
            `de-heatmap_scaling` = c(
                "row=Row (metabolite) scaling selected", "column=Column (sample) scaling", "both=Both",
                "none=None"),
            `de-heatmap_tree_cut_method` = c(
                "none=None selected", "k_clusters=K Clusters", "height_cutoff=Height Cutoff",
                "dynamic=Dynamic"),
            `de-table_assay` = "All=All selected",
            `de-table_contrast` = character(0),
            `de-table_significance` = c(
                "all=All Results", "significant=Significant Only selected", "up=Up-regulated",
                "down=Down-regulated"),
            `de-volcano_assay` = "Combined=Combined selected",
            `de-volcano_contrast` = character(0)),
        defaults = c(
            `de-da_q_val_thresh` = "0.05", `de-formula_string` = "~ 0 + group",
            `de-heatmap_assay` = "Combined", `de-heatmap_cluster_method` = "ward.D2",
            `de-heatmap_clustering` = "both", `de-heatmap_color_scheme` = "RdBu",
            `de-heatmap_cut_height` = "10", `de-heatmap_distance_method` = "euclidean",
            `de-heatmap_min_cluster_size` = "5", `de-heatmap_n_clusters` = "4", `de-heatmap_scaling` = "row",
            `de-heatmap_show_labels` = "FALSE", `de-heatmap_top_n` = "50",
            `de-heatmap_tree_cut_method` = "none", `de-table_assay` = "All", `de-table_max_rows` = "1000",
            `de-table_significance` = "significant", `de-treat_lfc_cutoff` = "0",
            `de-volcano_assay` = "Combined")),
    trx = list(
        html = "001abe86f81a8a3479caac17a5752f3ff4e054d31a1bb38993f44417dd56c7ac",
        inputs = c(
            "de-da_q_val_thresh", "de-da_results_tabs", "de-download_da_results", "de-run_da_analysis",
            "de-table_assay", "de-table_contrast", "de-table_max_rows", "de-table_significance",
            "de-treat_lfc_cutoff", "de-volcano_assay", "de-volcano_contrast", "de-volcano_glimma_container",
            "de-volcano_interactive"),
        outputs = c(
            "de-contrasts_display", "de-da_results_table", "de-da_status", "de-da_summary_stats",
            "de-formula_display", "de-volcano_glimma", "de-volcano_static"),
        testids = c("trx-da-download", "trx-da-results-tabs", "trx-da-run"),
        tabs = c("Volcano Plot", "DA Results Table"),
        labels = c(
            `de-da_q_val_thresh` = "FDR (BH) threshold:", `de-table_assay` = "Select Assay:",
            `de-table_contrast` = "Select Contrast:", `de-table_max_rows` = "Max Rows:",
            `de-table_significance` = "Show:",
            `de-treat_lfc_cutoff` = "Fold-change line (display only; no test):",
            `de-volcano_assay` = "Select Assay:", `de-volcano_contrast` = "Select Contrast:",
            `de-volcano_interactive` = "Interactive Plot (Glimma)"),
        help = c(
            "\n        edgeR quasi-likelihood: estimateDisp and glmQLFit(robust = TRUE), then glmQLFTest per contrast on the\n        one fit; Benjamini-Hochberg FDR. The formula is realised from the design: ~ 0 + group plus the RUVs factors\n        W_1..W_k as covariates. No fold-change test is applied.\n      ",
            "Contrasts defined in design matrix"),
        choices = list(
            `de-table_assay` = "All=All selected",
            `de-table_contrast` = character(0),
            `de-table_significance` = c(
                "all=All Results", "significant=Significant Only selected", "up=Up-regulated",
                "down=Down-regulated"),
            `de-volcano_assay` = "Combined=Combined selected",
            `de-volcano_contrast` = character(0)),
        defaults = c(
            `de-da_q_val_thresh` = "0.05", `de-table_assay` = "All", `de-table_max_rows` = "1000",
            `de-table_significance` = "significant", `de-treat_lfc_cutoff` = "0",
            `de-volcano_assay` = "Combined")),
    prot = list(
        html = "469abf274e873a5d9b607a8d4f7635d318cca6d8048d2d6cb4543356b450a765",
        inputs = c(
            "differential_expression-da_q_val_thresh", "differential_expression-da_results_tabs",
            "differential_expression-download_da_results", "differential_expression-formula_string",
            "differential_expression-heatmap_cluster_method", "differential_expression-heatmap_clustering",
            "differential_expression-heatmap_color_scheme", "differential_expression-heatmap_container",
            "differential_expression-heatmap_contrast", "differential_expression-heatmap_cut_height",
            "differential_expression-heatmap_distance_method",
            "differential_expression-heatmap_min_cluster_size", "differential_expression-heatmap_n_clusters",
            "differential_expression-heatmap_scaling", "differential_expression-heatmap_show_dendro",
            "differential_expression-heatmap_show_labels", "differential_expression-heatmap_top_n",
            "differential_expression-heatmap_tree_cut_method",
            "differential_expression-load_filtered_session", "differential_expression-run_da_analysis",
            "differential_expression-save_heatmap", "differential_expression-session_archive_standard",
            "differential_expression-session_archive_standard_progress",
            "differential_expression-table_contrast", "differential_expression-table_max_rows",
            "differential_expression-table_significance", "differential_expression-treat_lfc_cutoff",
            "differential_expression-volcano_contrast", "differential_expression-volcano_glimma_container",
            "differential_expression-volcano_interactive", "differential_expression-volcano_label_top_n",
            "differential_expression-volcano_show_labels"),
        outputs = c(
            "differential_expression-cluster_summary", "differential_expression-contrasts_display",
            "differential_expression-da_results_table", "differential_expression-da_status",
            "differential_expression-da_summary_stats", "differential_expression-heatmap_manual_save_warning",
            "differential_expression-heatmap_plot", "differential_expression-volcano_glimma",
            "differential_expression-volcano_plot_static"),
        testids = c(
            "prot-da-download", "prot-da-load-session", "prot-da-results-tabs", "prot-da-run"),
        tabs = c("Volcano Plot", "Heatmap", "DA Results Table"),
        labels = c(
            `differential_expression-da_q_val_thresh` = "Q-value Threshold:",
            `differential_expression-formula_string` = "Model Formula:",
            `differential_expression-heatmap_cluster_method` = "Clustering Method:",
            `differential_expression-heatmap_clustering` = "Apply Clustering:",
            `differential_expression-heatmap_color_scheme` = "Color Scheme:",
            `differential_expression-heatmap_contrast` = "Select Contrast:",
            `differential_expression-heatmap_cut_height` = "Cut Height:",
            `differential_expression-heatmap_distance_method` = "Distance Metric:",
            `differential_expression-heatmap_min_cluster_size` = "Min Cluster Size:",
            `differential_expression-heatmap_n_clusters` = "Number of Clusters:",
            `differential_expression-heatmap_scaling` = "Data Scaling:",
            `differential_expression-heatmap_show_dendro` = "Show Dendrogram",
            `differential_expression-heatmap_show_labels` = "Show Gene Labels",
            `differential_expression-heatmap_top_n` = "Top N Genes:",
            `differential_expression-heatmap_tree_cut_method` = "Tree Cutting:",
            `differential_expression-session_archive_standard` = "Saved session (optional):",
            `differential_expression-table_contrast` = "Select Contrast:",
            `differential_expression-table_max_rows` = "Max Rows:",
            `differential_expression-table_significance` = "Show:",
            `differential_expression-treat_lfc_cutoff` = "Log Fold-Change Cutoff:",
            `differential_expression-volcano_contrast` = "Select Contrast:",
            `differential_expression-volcano_interactive` = "Interactive Plot (Glimma)",
            `differential_expression-volcano_label_top_n` = "Number of Labels:",
            `differential_expression-volcano_show_labels` = "Label Top Proteins"),
        help = c(
            "Formula fetched from S4 object @args", "Contrasts defined in design matrix",
            "Load the latest session, a portable native DIA archive, or a historical RDS session."),
        choices = list(
            `differential_expression-heatmap_cluster_method` = c(
                "ward.D2=Ward (minimum variance) selected", "ward.D=Ward (original)",
                "complete=Complete linkage", "single=Single linkage", "average=Average linkage",
                "mcquitty=McQuitty (WPGMA)", "median=Median (WPGMC)", "centroid=Centroid (UPGMC)"),
            `differential_expression-heatmap_clustering` = c(
                "both=Both rows &amp; columns selected", "row=Rows only", "column=Columns only", "none=None"),
            `differential_expression-heatmap_color_scheme` = c(
                "RdBu=Red-Blue selected", "RdYlBu=Red-Yellow-Blue", "coolwarm=Blue-White-Red",
                "viridis=Viridis", "plasma=Plasma", "inferno=Inferno"),
            `differential_expression-heatmap_contrast` = character(0),
            `differential_expression-heatmap_distance_method` = c(
                "euclidean=Euclidean selected", "manhattan=Manhattan", "pearson=Pearson correlation",
                "spearman=Spearman correlation", "maximum=Maximum", "canberra=Canberra", "binary=Binary",
                "minkowski=Minkowski"),
            `differential_expression-heatmap_scaling` = c(
                "row=Row (gene) scaling selected", "column=Column (sample) scaling", "both=Both", "none=None"),
            `differential_expression-heatmap_tree_cut_method` = c(
                "k_clusters=Number of clusters selected", "height_cutoff=Height cutoff",
                "dynamic=Dynamic tree cutting", "none=No cutting"),
            `differential_expression-table_contrast` = character(0),
            `differential_expression-table_significance` = c(
                "all=All Results", "significant=Significant Only selected", "up=Up-regulated",
                "down=Down-regulated"),
            `differential_expression-volcano_contrast` = character(0)),
        defaults = c(
            `differential_expression-da_q_val_thresh` = "0.05",
            `differential_expression-formula_string` = "~ 0 + group",
            `differential_expression-heatmap_cluster_method` = "ward.D2",
            `differential_expression-heatmap_clustering` = "both",
            `differential_expression-heatmap_color_scheme` = "RdBu",
            `differential_expression-heatmap_cut_height` = "0.5",
            `differential_expression-heatmap_distance_method` = "euclidean",
            `differential_expression-heatmap_min_cluster_size` = "3",
            `differential_expression-heatmap_n_clusters` = "4",
            `differential_expression-heatmap_scaling` = "row",
            `differential_expression-heatmap_show_dendro` = "TRUE",
            `differential_expression-heatmap_show_labels` = "FALSE",
            `differential_expression-heatmap_top_n` = "50",
            `differential_expression-heatmap_tree_cut_method` = "k_clusters",
            `differential_expression-table_max_rows` = "1000",
            `differential_expression-table_significance` = "significant",
            `differential_expression-treat_lfc_cutoff` = "0",
            `differential_expression-volcano_label_top_n` = "10",
            `differential_expression-volcano_show_labels` = "TRUE")))
# OMIC_DA_FROZEN_UI_END

test_that("(a) the lipid, metab, trx-profile and proteomics DA UIs keep their HTML, ids, test ids, tab titles, labels, choices and defaults", {
    actual <- lapply(omicDaW1Uis(), omicDaW1UiParts)
    if (!is.null(omicDaW1DumpPath())) omicDaW1Dump("omicDaW1FrozenUi", actual)
    expect_named(actual, names(omicDaW1FrozenUi))
    for (omic in names(omicDaW1FrozenUi)) {
        for (part in names(omicDaW1FrozenUi[[omic]])) {
            expect_identical(actual[[omic]][[part]], omicDaW1FrozenUi[[omic]][[part]], info = paste(omic, part))
        }
    }
})

# ------------------------------------------------------------------------------------------------
# (b), (b-managed) and (c): the servers, driven in testServer.
# ------------------------------------------------------------------------------------------------

#' The browser and notification doubles, scoped to `env`. A notification is kept as
#' "[type duration id] message" (shiny's default duration 5 when the caller gives none), a removal as
#' "remove <id>", a modal as "modal"; a selector update as "<id>: <choices> -> <selected or value>",
#' and an update with a value is queued for the browser double, which sends it back as the input.
omicDaW1InstallDoubles <- function(log, env = parent.frame()) {
    note <- function(text) log$notes <- c(log$notes, text)
    notify <- function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
        type = c("default", "message", "warning", "error"), session = NULL) {
        note(sprintf("[%s %s %s] %s", type[[1L]], if (is.null(duration)) "NULL" else format(duration), id %||% "-",
            paste(as.character(ui), collapse = " ")))
        invisible(id %||% "notification")
    }
    queue <- function(session, input_id, value = NULL, choices = NULL) {
        if (is.null(session)) return(invisible(NULL))
        id <- session$ns(input_id)
        log$updates <- c(log$updates, sprintf("%s: %s -> %s", id,
            if (is.null(choices)) "<same>" else paste(as.character(unlist(choices, use.names = FALSE)), collapse = "|"),
            if (is.null(value)) "<none>" else paste(format(value), collapse = "|")))
        if (!is.null(value)) log$pending[[length(log$pending) + 1L]] <- list(id = id, value = value)
        invisible(NULL)
    }
    domain <- function() shiny::getDefaultReactiveDomain()
    bindings <- list(
        updateSelectInput = function(session = domain(), inputId, label = NULL, choices = NULL, selected = NULL) {
            queue(session, inputId, selected, choices)
        },
        updateSelectizeInput = function(session = domain(), inputId, label = NULL, choices = NULL, selected = NULL,
            options = list(), server = FALSE) queue(session, inputId, selected, choices),
        updateTextAreaInput = function(session = domain(), inputId, label = NULL, value = NULL, placeholder = NULL) {
            queue(session, inputId, value)
        },
        updateTextInput = function(session = domain(), inputId, label = NULL, value = NULL, placeholder = NULL) {
            queue(session, inputId, value)
        },
        updateNumericInput = function(session = domain(), inputId, label = NULL, value = NULL, min = NULL, max = NULL,
            step = NULL) queue(session, inputId, value),
        updateCheckboxInput = function(session = domain(), inputId, label = NULL, value = NULL) queue(session, inputId, value),
        updateTabsetPanel = function(session = domain(), inputId, selected = NULL) queue(session, inputId, selected),
        showNotification = notify,
        removeNotification = function(id, session = domain()) {
            note(paste("remove", id))
            invisible(NULL)
        },
        showModal = function(ui, session = domain()) {
            note("modal")
            invisible(NULL)
        },
        removeModal = function(session = domain()) invisible(NULL))
    do.call(testthat::local_mocked_bindings, c(bindings, list(.package = "shiny", .env = env)))
    imported <- intersect(names(bindings), unlist(getNamespaceImports("MultiScholaR")$shiny, use.names = FALSE))
    if (length(imported)) do.call(testthat::local_mocked_bindings, c(bindings[imported], list(.package = "MultiScholaR", .env = env)))
    invisible(log)
}

#' A digest of what a plot or table kernel returned: a ggplot's data, a heatmap's matrix (and its row
#' clusters), a table's rows.
omicDaW1ValueDigest <- function(value) {
    short <- function(x) substr(guiEntryDigest(x), 1L, 16L)
    if (is.null(value)) return("<NULL>")
    if (inherits(value, "ggplot")) return(paste("ggplot", short(value$data)))
    if (methods::is(value, "Heatmap")) return(paste("Heatmap", short(value@matrix)))
    if (is.list(value) && !is.null(value$plot)) {
        return(sprintf("list(plot = %s, row_clusters = %s)", omicDaW1ValueDigest(value$plot), short(value$row_clusters)))
    }
    if (is.data.frame(value)) return(sprintf("data.frame %dx%d %s", nrow(value), ncol(value), short(value)))
    paste0("<", class(value)[[1L]], ">")
}

#' One kernel call as text: its name, its scalar arguments (and the save's parameters) and a digest of
#' what it returned (of its data, for the table widget).
omicDaW1CallText <- function(name, args, value) {
    keep <- vapply(args, function(arg) is.atomic(arg) && !is.null(arg) && length(arg) <= 3L && is.null(dim(arg)), logical(1))
    named <- names(args) %||% rep("", length(args))
    shown <- vapply(which(keep), function(i) paste0(if (nzchar(named[[i]])) paste0(named[[i]], "=") else "",
        paste(format(args[[i]]), collapse = "|")), "")
    if (is.list(args$params_list)) {
        params <- Filter(function(arg) is.atomic(arg) && length(arg) <= 3L, args$params_list)
        shown <- c(shown, paste0("params_list$", names(params), "=",
            vapply(params, function(arg) paste(format(arg), collapse = "|"), "")))
    }
    data <- if (identical(name, "DT::datatable")) omicDaW1ValueDigest(args[[1L]]) else omicDaW1ValueDigest(value)
    sprintf("%s(%s) -> %s", name, paste(shown, collapse = ", "), data)
}

#' Pass-through recorders on the plot, heatmap, save and table kernels the renders call, scoped to
#' `env`; each call is appended to `acc$calls`. The kernels are the core's (`func_*`, the trx static
#' volcano), so they keep their names through W2 to W5 and are looked up by name.
omicDaW1InstallRecorders <- function(acc, env = parent.frame()) {
    ns <- asNamespace("MultiScholaR")
    kernels <- c("generateLipidDAVolcanoStatic", "generateMetabDAVolcanoStatic", "generateProtDAVolcanoStatic",
        "trxDaStaticVolcano", "generateLipidDAHeatmap", "generateMetabDAHeatmap", "generateProtDAHeatmap",
        "save_heatmap_products")
    recorder <- function(name, original) {
        force(name)
        force(original)
        function(...) {
            value <- original(...)
            acc$calls <- c(acc$calls, omicDaW1CallText(name, list(...), value))
            value
        }
    }
    mocks <- lapply(stats::setNames(kernels, kernels), function(name) recorder(name, get(name, envir = ns)))
    do.call(testthat::local_mocked_bindings, c(mocks, list(.package = "MultiScholaR", .env = env)))
    table <- recorder("DT::datatable", DT::datatable)
    testthat::local_mocked_bindings(datatable = table, .package = "DT", .env = env)
    if ("datatable" %in% unlist(getNamespaceImports("MultiScholaR")$DT, use.names = FALSE)) {
        testthat::local_mocked_bindings(datatable = table, .package = "MultiScholaR", .env = env)
    }
    invisible(acc)
}

#' The tab's reactive values: each field's name and a description of what it references (a short
#' atomic value as itself; a list by its names; a data frame by its shape; any other object by its
#' class), never its content.
omicDaW1TabValues <- function(values) {
    if (is.null(values)) return("<none>")
    fields <- shiny::isolate(shiny::reactiveValuesToList(values))
    fields <- fields[order(names(fields))]
    vapply(fields, function(value) {
        if (is.null(value)) return("<NULL>")
        if (is.data.frame(value)) return(sprintf("<data.frame %dx%d>", nrow(value), ncol(value)))
        if (isS4(value)) return(paste0("<S4 ", class(value)[[1L]], ">"))
        if (is.list(value) && !is.object(value)) return(paste0("<list: ", paste(names(value), collapse = ","), ">"))
        if (is.atomic(value) && is.null(dim(value)) && length(value) <= 10L && !is.object(value)) {
            return(paste(as.character(value), collapse = "|"))
        }
        paste0("<", paste(class(value), collapse = "|"), ">")
    }, "")
}

#' The rows the tab holds: a digest of each long table in its results list (lipid's, metab's and
#' trx's; proteomics' results are a store manifest, digested without its references).
omicDaW1Rows <- function(values) {
    if (is.null(values)) return("<none>")
    results <- shiny::isolate(values$da_results_list)
    if (is.null(results)) return("<NULL>")
    long <- results[grepl("_long$", names(results)) & vapply(results, is.data.frame, logical(1))]
    # A manifest's store references (its `source`, its `results`) carry process-local entity ids: its other fields only.
    if (!length(long)) {
        kept <- intersect(c("schema", "protein_id_column", "sample_id", "contrasts", "parameters", "model"), names(results))
        return(paste("<manifest>", paste(kept, collapse = ","), substr(guiEntryDigest(results[kept]), 1L, 16L)))
    }
    vapply(long, function(table) sprintf("%dx%d %s", nrow(table), ncol(table), substr(guiEntryDigest(table), 1L, 16L)), "")
}

#' A message as captured: the entry suite's masks, then without terminal styling or non-ASCII
#' characters (cli's bullets differ between a UTF-8 console and check's).
omicDaW1Text <- function(x, roots) {
    x <- gsub("\033\\[[0-9;]*m", "", guiEntryMask(x, roots))
    iconv(x, "UTF-8", "ASCII", sub = "")
}

#' What a server returned, by kind (D2 changes lipid's and metab's; (c) records it apart).
omicDaW1Returned <- function(value) {
    if (shiny::is.reactivevalues(value)) return(paste("reactivevalues:", paste(sort(names(value)), collapse = ",")))
    if (inherits(value, "shinyoutput")) return("shinyoutput")
    if (is.list(value)) return(paste("list:", paste(names(value), collapse = ",")))
    paste0("<", class(value)[[1L]], ">")
}

#' Drive one case: mount its server in testServer, set the DA UI's defaults (as the browser sends
#' them), run the steps (named functions of the context) and capture after each step.
#' @param setup The case's fixture: `omic`, `server(acc)` (a server function), `ns` (the DA module
#'   id), `ui` (the DA UI, for its defaults), `bus()` (the workflow bus), `roots` (masked),
#'   `files_root`, and for a managed case `probe`; for proteomics `tab` (the selected tab).
#' @return `captures` (by step) and `returned` (what the DA server returned).
omicDaW1Drive <- function(setup, steps) {
    log <- new.env(parent = emptyenv())
    log$notes <- character()
    log$updates <- character()
    log$pending <- list()
    omicDaW1InstallDoubles(log)
    acc <- new.env(parent = emptyenv())
    acc$calls <- character()
    acc$captures <- list()
    acc$previous <- list()
    acc$extra <- list()
    acc$operations <- 0L
    acc$errors <- character()
    omicDaW1InstallRecorders(acc)
    # The tab's reactive values where the server does not return them (lipid's and metab's servers at
    # the start commit): their state factories, recorded while they exist.
    for (factory in c("initializeLipidDaServerState", "createMetabDaServerState")) {
        if (exists(factory, envir = asNamespace("MultiScholaR"), inherits = FALSE)) {
            hook <- local({
                make <- get(factory, envir = asNamespace("MultiScholaR"))
                function(...) {
                    values <- make(...)
                    acc$tab <- values
                    values
                }
            })
            do.call(testthat::local_mocked_bindings, c(stats::setNames(list(hook), factory), list(.package = "MultiScholaR")))
        }
    }
    withr::local_seed(20261109L)
    prefix <- paste0(setup$ns, "-")
    # An observer's error reaches no notification; the session's unhandled-error callback records it.
    server <- local({
        mount <- setup$server(acc)
        function(input, output, session) {
            session$onUnhandledError(function(error) acc$errors <- c(acc$errors, paste("<unhandled>", conditionMessage(error))))
            mount(input, output, session)
        }
    })
    text <- function(x) omicDaW1Text(x, setup$roots)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(server, {
        apply_updates <- function() {
            for (i in seq_len(20L)) {
                if (!length(log$pending)) return(invisible())
                pending <- log$pending
                log$pending <- list()
                values <- stats::setNames(lapply(pending, `[[`, "value"), vapply(pending, `[[`, "", "id"))
                do.call(session$setInputs, values[!duplicated(names(values), fromLast = TRUE)])
            }
        }
        set <- function(...) {
            values <- list(...)
            names(values) <- ifelse(startsWith(names(values), "norm-"), names(values), paste0(prefix, names(values)))
            tryCatch({
                do.call(session$setInputs, values)
                apply_updates()
            }, error = function(e) acc$errors <- c(acc$errors, conditionMessage(e)))
            invisible(NULL)
        }
        read <- function(name) tryCatch(output[[paste0(prefix, name)]], error = function(e) {
            if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", text(conditionMessage(e)))
        })
        defaults <- omicsUiInputDefaults(setup$ui("x"))
        names(defaults) <- sub("^x-", "", names(defaults))
        do.call(set, defaults)
        capture <- function(step) {
            session$flushReact()
            apply_updates()
            bus <- setup$bus()
            manager <- shiny::isolate(bus$state_manager)
            out_names <- sort(names(session$.__enclos_env__$private$outs))
            out_names <- sub(paste0("^", prefix), "", out_names[startsWith(out_names, prefix)])
            # The interactive volcano is D3's: lipid's, metab's and trx's are not read (a managed page would be
            # requested); proteomics' DA-time page is.
            shown <- setdiff(out_names, if (!identical(setup$omic, "prot")) "volcano_glimma")
            # A download output renders its file when read: it is recorded by name (the download step reads it).
            values <- vapply(shown, function(name) {
                if (startsWith(name, "download_")) return("<download>")
                guiEntryOutputValue(read(name), setup$roots)
            }, "")
            returned <- session$userData$da
            acc$returned <- omicDaW1Returned(returned)
            operations <- NULL
            if (!is.null(setup$probe)) {
                operations <- setup$probe$operations[seq_along(setup$probe$operations) > acc$operations]
                acc$operations <- length(setup$probe$operations)
            }
            files <- unlist(guiEntryArtifacts(setup$files_root, setup$roots))
            record <- c(list(
                outputs = out_names,
                values = values,
                calls = vapply(acc$calls, text, "", USE.NAMES = FALSE),
                notes = vapply(log$notes, text, "", USE.NAMES = FALSE),
                updates = vapply(log$updates, text, "", USE.NAMES = FALSE),
                tab_status = vapply(shiny::isolate(bus$tab_status), as.character, ""),
                history = guiEntryHistory(manager),
                states = tryCatch(omicDaW1Sha(guiEntryStateDigests(manager)), error = function(e) "<none>"),
                tab = omicDaW1TabValues(if (shiny::is.reactivevalues(returned)) returned else acc$tab),
                rows = omicDaW1Rows(if (shiny::is.reactivevalues(returned)) returned else acc$tab),
                files = stats::setNames(as.character(files), names(files)),
                operations = operations,
                errors = vapply(acc$errors, text, "", USE.NAMES = FALSE)), acc$extra)
            # A field the step left as the previous step had it is held as "<unchanged>".
            for (field in c("tab_status", "history", "states", "tab", "rows", "files")) {
                value <- record[[field]]
                if (identical(value, acc$previous[[field]])) record[[field]] <- "<unchanged>"
                acc$previous[[field]] <- value
            }
            log$notes <- character()
            log$updates <- character()
            acc$calls <- character()
            acc$errors <- character()
            acc$extra <- list()
            acc$captures[[step]] <- record
        }
        ctx <- list(set = set, read = read, bus = setup$bus, roots = setup$roots, tab = setup$tab,
            keep = function(name, value) acc$extra[[name]] <- value)
        for (step in names(steps)) {
            steps[[step]](ctx)
            capture(step)
        }
    }))))
    list(captures = acc$captures, returned = acc$returned)
}

#' The download as the browser saves it: the file's name (its date masked) and its bytes' digest.
omicDaW1Download <- function(ctx) {
    path <- ctx$read("download_da_results")
    ctx$keep("download", if (is.character(path) && length(path) == 1L && file.exists(path)) {
        c(file = guiEntryMask(basename(path)), sha256 = digest::digest(file = path, algo = "sha256"))
    } else {
        paste("<no file>", paste(path, collapse = " "))
    })
}

#' Lipid's or metab's fixture: MCI-022's / MCI-016's object (LCMS_Pos and GCMS, groups WT, KO and RES),
#' the case's contrasts and the MCI legacy session file, on an owner of the omic with the project's
#' paths. The bus is the MCI fixtures' (an environment, `module_ci_{lipid,metab}_da_workflow`).
#' @param state Whether the bus has an owner holding the object (its head); FALSE: no owner on the
#'   bus, so nothing is loaded (the omic's "no data" refusal).
#' @param session_assays The session file's `assay_names` (NULL: the object's); `session_file = FALSE`
#'   writes none.
omicDaW1FeatureSetup <- function(omic, contrasts = "multi_group", state = TRUE, session_file = TRUE,
    session_assays = NULL, managed = FALSE, env = parent.frame()) {
    spec <- switch(omic,
        lipid = list(omic_type = "lipidomics", object = module_ci_lipid_da_object, contrasts = module_ci_lipid_da_contrasts,
            session = module_ci_lipid_da_session_data, state_name = "lipid_correlation_filtered",
            file = "lipid_filtered_session_data_latest.rds", server = mod_lipid_da_server, ui = mod_lipid_da_ui),
        metab = list(omic_type = "metabolomics", object = module_ci_metab_da_object, contrasts = module_ci_metab_da_contrasts,
            session = module_ci_metab_da_session_data, state_name = "metab_correlation_filtered",
            file = "metab_filtered_session_data_latest.rds", server = mod_metab_da_server, ui = mod_metab_da_ui))
    root <- normalizePath(withr::local_tempdir(pattern = paste0("kit106-", omic, "-"), .local_envir = env))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da_output"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE)
    object <- spec$object()
    bus <- list2env(list(contrasts_tbl = spec$contrasts(contrasts), design_matrix = object@design_matrix,
        config_list = list(da = list(mode = "module-ci")), tab_status = list(setup_import = "complete",
            design_matrix = "complete", quality_control = "complete", normalization = "complete",
            differential_analysis = "pending")), parent = emptyenv())
    if (state) {
        bus <- omicsTestOwner(object, spec$state_name, workflow = bus, omic_type = spec$omic_type, experiment_paths = paths,
            env = env)
    }
    if (session_file) saveRDS(spec$session(assay_names = session_assays), file.path(paths$source_dir, spec$file))
    probe <- NULL
    if (managed) {
        probe <- new.env()
        omicsManagedChainMocks(list(owner = bus$state_manager, workflow = bus, paths = paths, source_dir = paths$source_dir),
            probe, env = env)
    }
    server <- function(acc) {
        function(input, output, session) {
            if (managed) {
                runner <- installWorkflowOmicsActions(session, bus, spec$omic_type, "kit106", paths)
                runner$progress <- module_ci_trx_quiet_progress
            }
            session$userData$da <- spec$server("de", bus, paths, spec$omic_type, "kit106")
        }
    }
    list(omic = omic, server = server, ns = "de", ui = spec$ui, bus = function() bus, roots = root, files_root = root,
        probe = probe)
}

#' Trx's fixture: MCI-030's project (`module_ci_trx_project(module_ci_trx_corpus())`, at
#' `trx_qc_complete`) with the normalisation and DA tabs under the transcriptomics profile, as the
#' shell mounts them (the normalisation first, on the same session, makes the RUVs head the DA runs on).
omicDaW1TrxSetup <- function(managed = FALSE, env = parent.frame()) {
    testthat::skip_if_not_installed("edgeR")
    fixture <- module_ci_trx_project(module_ci_trx_corpus(), env = env)
    paths <- fixture$paths
    project <- dirname(dirname(paths$source_dir))
    probe <- NULL
    if (managed) {
        probe <- new.env()
        omicsManagedChainMocks(fixture, probe, env = env)
    }
    server <- function(acc) {
        function(input, output, session) {
            if (managed) {
                runner <- installWorkflowOmicsActions(session, fixture$workflow, "transcriptomics", "kit106", paths)
                runner$progress <- module_ci_trx_quiet_progress
            }
            profile <- trxOmicProfile()
            mod_lipid_norm_server("norm", fixture$workflow, paths, "transcriptomics", "kit106", profile = profile)
            session$userData$da <- mod_lipid_da_server("de", fixture$workflow, paths, "transcriptomics", "kit106",
                profile = profile)
        }
    }
    list(omic = "trx", server = server, ns = "de", ui = function(id) mod_lipid_da_ui(id, profile = trxOmicProfile()),
        bus = function() fixture$workflow, roots = unique(c(normalizePath(project, mustWork = FALSE), project)),
        files_root = project, probe = probe)
}

#' Proteomics' fixture: a DIA owner with the project's paths (as test-workflow-saved-da.R makes one),
#' the bus as the proteomics shell keeps it, and the tab selection (`selected_tab`) on `qc`. MCI-009's
#' object has its row names dropped: the DIA reference's tables refuse a data frame with row names,
#' and an application's tables have none.
#' @param state Whether the owner holds the object as `correlation_filtered` (FALSE: an empty owner).
#' @param session_file Whether MCI-009's session payload is written as the legacy export the owner adopts
#'   on the load (its states as records, its workflow type), with MCI-009's annotation file.
omicDaW1ProtSetup <- function(contrasts = "two_group", state = TRUE, session_file = FALSE, env = parent.frame()) {
    root <- normalizePath(withr::local_tempdir(pattern = "kit106-prot-", .local_envir = env))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da_output"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE)
    design <- module_ci_prot_da_design()
    rownames(design) <- NULL
    quant <- module_ci_prot_da_table()
    rownames(quant) <- NULL
    object <- module_ci_prot_da_object(design = design, protein_quant_table = quant)
    if (session_file) {
        payload <- module_ci_prot_da_session_payload()
        record <- list(data = object, config = list(), description = "MCI-009 fixture state")
        payload$r6_complete_states <- list(raw_data_s4 = record, correlation_filtered = record)
        payload$current_s4_object <- object
        payload$design_matrix <- design
        payload$workflow_type <- "LFQ"
        module_ci_prot_da_write_session(paths$source_dir, payload)
    }
    tab <- shiny::reactiveVal("qc")
    holder <- new.env(parent = emptyenv())
    server <- function(acc) {
        function(input, output, session) {
            manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit106-prot", "DIA"),
                WorkflowState$new(audit_enabled = FALSE))
            session$onSessionEnded(function() closeDIAStateAdapter(manager))
            manager$.context$reference_states <- TRUE
            manager$.context$experiment_paths <- paths
            session$userData$dia_capacity_manager <- manager
            if (state) manager$saveState("correlation_filtered", object, list(), "MCI-009 fixture state")
            holder$bus <- shiny::reactiveValues(state_manager = manager,
                contrasts_tbl = if (!is.null(contrasts)) module_ci_prot_da_contrasts(contrasts),
                tab_status = list(normalization = "complete", differential_expression = "pending",
                    enrichment_analysis = "disabled"),
                da_analysis_results_list = NULL, da_ui_params = NULL, uniprot_dat_cln = NULL, state_update_trigger = NULL)
            session$userData$da <- mod_prot_da_server("differential_expression", holder$bus, paths, "proteomics", "kit106",
                selected_tab = function() tab())
        }
    }
    list(omic = "prot", server = server, ns = "differential_expression", ui = mod_prot_da_ui,
        bus = function() holder$bus, roots = root, files_root = root, tab = tab)
}

#' The steps of lipid's and metab's cases (both fixtures hold LCMS_Pos and GCMS; contrasts KO_vs_WT and
#' RES_vs_WT).
omicDaW1FeatureSteps <- function(case) {
    run <- function(n) {
        function(ctx) {
            ctx$set(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
            ctx$set(run_da_analysis = n)
        }
    }
    load <- function(ctx) ctx$set(load_filtered_session = 1L)
    switch(case,
        refusals = list(load_missing = load, run_no_data = run(1L)),
        mismatch = list(load_mismatch = load),
        reordered = list(load_reordered = load),
        no_contrasts = list(run_no_contrasts = run(1L)),
        invalid_preflight = list(run_invalid = run(1L)),
        main = list(
            start = function(ctx) NULL,
            load = load,
            run = run(1L),
            select = function(ctx) ctx$set(volcano_contrast = "RES_vs_WT", volcano_assay = "LCMS_Pos",
                heatmap_contrast = "RES_vs_WT", heatmap_assay = "GCMS", heatmap_top_n = 20, heatmap_clustering = "row",
                heatmap_scaling = "column", heatmap_cluster_method = "average", heatmap_distance_method = "manhattan",
                heatmap_color_scheme = "viridis", heatmap_show_labels = TRUE, table_contrast = "RES_vs_WT",
                table_assay = "LCMS_Pos", table_significance = "up", table_max_rows = 100),
            tree_k = function(ctx) ctx$set(heatmap_clustering = "both", heatmap_tree_cut_method = "k_clusters",
                heatmap_n_clusters = 3),
            tree_height = function(ctx) ctx$set(heatmap_tree_cut_method = "height_cutoff", heatmap_cut_height = 5),
            tree_dynamic = function(ctx) ctx$set(heatmap_tree_cut_method = "dynamic", heatmap_min_cluster_size = 2),
            save_heatmap = function(ctx) ctx$set(save_heatmap = 1L),
            table_down = function(ctx) ctx$set(table_assay = "All", table_significance = "down"),
            table_limit = function(ctx) ctx$set(table_significance = "all", table_max_rows = 5),
            table_combined = function(ctx) ctx$set(table_assay = "Combined"),
            download = omicDaW1Download),
        no_friendly = list(
            run = run(1L),
            table = function(ctx) ctx$set(volcano_contrast = "groupKO-groupWT", volcano_assay = "LCMS_Pos",
                table_contrast = "groupKO-groupWT", table_assay = "All", table_significance = "all")),
        managed = list(
            run = run(1L),
            select = function(ctx) ctx$set(volcano_contrast = "RES_vs_WT", volcano_assay = "GCMS",
                table_contrast = "RES_vs_WT", table_assay = "GCMS", table_significance = "all")))
}

#' The steps of trx's cases (the corpus's one contrast; the assay `counts`).
omicDaW1TrxSteps <- function(case) {
    norm <- function(ctx) {
        ctx$set(`norm-ruv_mode` = "automatic")
        ctx$set(`norm-run_normalization` = 1L)
    }
    run <- function(n) {
        function(ctx) {
            ctx$set(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
            ctx$set(run_da_analysis = n)
        }
    }
    switch(case,
        refusals = list(
            run_no_data = run(1L),
            norm = norm,
            run_no_contrasts = function(ctx) {
                bus <- ctx$bus()
                bus$contrasts_tbl <- bus$contrasts_tbl[0, , drop = FALSE]
                run(2L)(ctx)
            },
            run_invalid = function(ctx) {
                bus <- ctx$bus()
                bus$contrasts_tbl <- module_ci_trx_builder_contrasts("Missing_vs_RPMI=groupMissing-groupRPMI")
                run(3L)(ctx)
            }),
        main = list(
            norm = norm,
            run = run(1L),
            # The table offers "All" only (the long table has no assay column).
            select = function(ctx) ctx$set(volcano_assay = "counts", table_significance = "up", table_max_rows = 100),
            table_down = function(ctx) ctx$set(table_assay = "All", table_significance = "down"),
            table_limit = function(ctx) ctx$set(table_significance = "all", table_max_rows = 5),
            table_combined = function(ctx) ctx$set(table_assay = "Combined"),
            download = omicDaW1Download),
        managed = list(
            norm = norm,
            run = run(1L),
            select = function(ctx) ctx$set(volcano_assay = "counts", table_significance = "all")))
}

#' The steps of proteomics' cases (MCI-009's contrast B_vs_A).
omicDaW1ProtSteps <- function(case) {
    run <- function(n) {
        function(ctx) {
            ctx$set(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
            ctx$set(run_da_analysis = n)
        }
    }
    open <- function(ctx) ctx$tab("da")
    load <- function(ctx) ctx$set(load_filtered_session = 1L)
    switch(case,
        refusals = list(load_missing = load, run_no_data = run(1L)),
        load = list(load = load, run = run(1L)),
        main = list(
            start = function(ctx) NULL,
            open = open,
            run = run(1L),
            select = function(ctx) ctx$set(volcano_contrast = "B_vs_A", heatmap_contrast = "B_vs_A", heatmap_top_n = 20,
                heatmap_clustering = "row", heatmap_scaling = "column", heatmap_cluster_method = "average",
                heatmap_distance_method = "manhattan", heatmap_color_scheme = "viridis", heatmap_show_labels = TRUE,
                table_contrast = "B_vs_A", table_significance = "up", table_max_rows = 100),
            volcano_unticked = function(ctx) ctx$set(volcano_interactive = FALSE),
            volcano_ticked = function(ctx) ctx$set(volcano_interactive = TRUE),
            volcano_labels = function(ctx) ctx$set(volcano_show_labels = FALSE, volcano_label_top_n = 3),
            tree_k = function(ctx) ctx$set(heatmap_clustering = "both", heatmap_tree_cut_method = "k_clusters",
                heatmap_n_clusters = 3),
            tree_height = function(ctx) ctx$set(heatmap_tree_cut_method = "height_cutoff", heatmap_cut_height = 5),
            tree_dynamic = function(ctx) ctx$set(heatmap_tree_cut_method = "dynamic", heatmap_min_cluster_size = 2),
            tree_none = function(ctx) ctx$set(heatmap_tree_cut_method = "none"),
            dendrogram_off = function(ctx) ctx$set(heatmap_tree_cut_method = "k_clusters", heatmap_show_dendro = FALSE),
            save_heatmap = function(ctx) ctx$set(save_heatmap = 1L),
            table_down = function(ctx) ctx$set(table_significance = "down"),
            table_limit = function(ctx) ctx$set(table_significance = "all", table_max_rows = 5),
            download = omicDaW1Download),
        invalid_preflight = list(open = open, run_invalid = run(1L)),
        no_contrasts = list(open = open, run_no_contrasts = run(1L)))
}

#' Each case: its fixture (made in `env`) and its steps.
omicDaW1CaseSpec <- function(case, env) {
    feature <- function(omic, name, ...) list(setup = omicDaW1FeatureSetup(omic, ..., env = env),
        steps = omicDaW1FeatureSteps(name))
    switch(case,
        lipid_refusals = feature("lipid", "refusals", state = FALSE, session_file = FALSE),
        lipid_mismatch = feature("lipid", "mismatch", session_assays = c("LCMS_Pos", "LCMS_Neg")),
        lipid_reordered = feature("lipid", "reordered", session_assays = c("GCMS", "LCMS_Pos")),
        lipid_no_contrasts = feature("lipid", "no_contrasts", contrasts = "no_contrasts", session_file = FALSE),
        lipid_invalid_preflight = feature("lipid", "invalid_preflight", contrasts = "invalid_term", session_file = FALSE),
        lipid_main = feature("lipid", "main"),
        lipid_no_friendly = feature("lipid", "no_friendly", session_file = FALSE),
        lipid_managed = feature("lipid", "managed", session_file = FALSE, managed = TRUE),
        metab_refusals = feature("metab", "refusals", state = FALSE, session_file = FALSE),
        metab_mismatch = feature("metab", "mismatch", session_assays = c("LCMS_Pos", "LCMS_Neg")),
        metab_reordered = feature("metab", "reordered", session_assays = c("GCMS", "LCMS_Pos")),
        metab_no_contrasts = feature("metab", "no_contrasts", contrasts = "no_contrasts", session_file = FALSE),
        metab_invalid_preflight = feature("metab", "invalid_preflight", contrasts = "invalid_term", session_file = FALSE),
        metab_main = feature("metab", "main"),
        metab_no_friendly = feature("metab", "no_friendly", session_file = FALSE),
        metab_managed = feature("metab", "managed", session_file = FALSE, managed = TRUE),
        trx_refusals = list(setup = omicDaW1TrxSetup(env = env), steps = omicDaW1TrxSteps("refusals")),
        trx_main = list(setup = omicDaW1TrxSetup(env = env), steps = omicDaW1TrxSteps("main")),
        trx_managed = list(setup = omicDaW1TrxSetup(managed = TRUE, env = env), steps = omicDaW1TrxSteps("managed")),
        prot_refusals = list(setup = omicDaW1ProtSetup(state = FALSE, env = env), steps = omicDaW1ProtSteps("refusals")),
        prot_load = list(setup = omicDaW1ProtSetup(state = FALSE, session_file = TRUE, env = env),
            steps = omicDaW1ProtSteps("load")),
        prot_main = list(setup = omicDaW1ProtSetup(env = env), steps = omicDaW1ProtSteps("main")),
        prot_invalid_preflight = list(setup = omicDaW1ProtSetup(contrasts = "invalid_term", env = env),
            steps = omicDaW1ProtSteps("invalid_preflight")),
        prot_no_contrasts = list(setup = omicDaW1ProtSetup(contrasts = NULL, env = env),
            steps = omicDaW1ProtSteps("no_contrasts")))
}

#' The long table of a reference DA's results without its `friendly_name` column (D7: the table
#' filter and selectors on a table that lacks it), scoped to `env`.
omicDaW1DropFriendlyName <- function(omic, env = parent.frame()) {
    name <- if (identical(omic, "lipid")) "runLipidReferenceDaAnalysis" else "runMetabReferenceDaAnalysis"
    long <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"
    original <- get(name, envir = asNamespace("MultiScholaR"))
    mock <- function(...) {
        out <- original(...)
        out$results[[long]]$friendly_name <- NULL
        out
    }
    do.call(testthat::local_mocked_bindings, c(stats::setNames(list(mock), name), list(.package = "MultiScholaR", .env = env)))
}

omicDaW1Runs <- new.env(parent = emptyenv())

#' One case's captures, run once per process (each case has its own fixture and session).
omicDaW1Case <- function(case) {
    if (is.null(omicDaW1Runs[[case]])) {
        env <- environment()
        if (case %in% c("lipid_no_friendly", "metab_no_friendly")) omicDaW1DropFriendlyName(sub("_no_friendly$", "", case))
        spec <- omicDaW1CaseSpec(case, env)
        omicDaW1Runs[[case]] <- omicDaW1Drive(spec$setup, spec$steps)
        if (!is.null(omicDaW1DumpPath())) {
            cat(sprintf("omicDaW1Frozen[[\"%s\"]] <- %s\n\n", case, omicDaW1Code(omicDaW1Runs[[case]]$captures)),
                file = omicDaW1DumpPath(), append = TRUE)
        }
    }
    omicDaW1Runs[[case]]
}

# The frozen captures, per case (the dump's entries, pasted).
# OMIC_DA_FROZEN_START
omicDaW1Frozen <- list()
omicDaW1Frozen[["lipid_refusals"]] <- list(
    load_missing = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] Session file not found: <ROOT>/source/lipid_filtered_session_data_latest.rds",
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = character(0),
        states = "743a1d2c832c72b269ed1dc4181d24657bf9c7296c6f9b1c4934328d5583605e",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = character(0),
        operations = NULL,
        errors = character(0)),
    run_no_data = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] No lipidomics data loaded. Please load a filtered session first.",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["lipid_mismatch"]] <- list(
    load_mismatch = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Data loaded. Ready to run analysis.", da_summary_stats = "<req>",
            download_da_results = "<download>", heatmap_manual_save_warning = "<NULL>",
            heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL loading_session] Loading filtered session...", "remove loading_session",
            "[error 10 -] Error loading session: Loaded lipid session assay_names do not match S4 lipid_data assays: session [LCMS_Pos, LCMS_Neg], S4 [LCMS_Pos, GCMS]"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<NULL>",
            current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "lipid_correlation_filtered", da_results_list = "<NULL>",
            formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/lipid_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["lipid_reordered"]] <- list(
    load_reordered = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Data loaded. Ready to run analysis.", da_summary_stats = "<req>",
            download_da_results = "<download>", heatmap_manual_save_warning = "<NULL>",
            heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL loading_session] Loading filtered session...", "remove loading_session",
            "[message 3 -] Session loaded successfully!"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|GCMS|LCMS_Pos -> <none>",
            "de-heatmap_assay: Combined|GCMS|LCMS_Pos -> <none>",
            "de-table_assay: All|GCMS|LCMS_Pos -> <none>", "de-formula_string: <same> -> ~ 0 + group"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "FALSE", assays_available = "GCMS|LCMS_Pos",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<NULL>",
            current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "lipid_correlation_filtered", da_results_list = "<NULL>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/lipid_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["lipid_no_contrasts"]] <- list(
    run_no_contrasts = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>",
            contrasts_display = "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab.",
            da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] No contrasts defined. Please define contrasts in the design tab.",
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["lipid_invalid_preflight"]] <- list(
    run_invalid = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "Missing_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[error 5 -] Lipidomics DA contrast 'groupMISSING-groupWT' references terms absent from the model matrix: groupMISSING"),
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["lipid_main"]] <- list(
    start = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/lipid_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    load = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Data loaded. Ready to run analysis.", da_summary_stats = "<req>",
            download_da_results = "<download>", heatmap_manual_save_warning = "<NULL>",
            heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL loading_session] Loading filtered session...", "remove loading_session",
            "[message 3 -] Session loaded successfully!"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> <none>",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> <none>",
            "de-table_assay: All|LCMS_Pos|GCMS -> <none>", "de-formula_string: <same> -> ~ 0 + group"),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "FALSE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<NULL>",
            current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "lipid_correlation_filtered", da_results_list = "<NULL>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/lipid_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 1bca5d30123f39779d9e6238fe985fcedf9631330542910d6850ffc05db0b701",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateLipidDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=KO_vs_WT, selected_assay=Combined, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=KO_vs_WT, selected_assay=Combined, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>",
            "generateLipidDAVolcanoStatic(selected_contrast=KO_vs_WT, selected_assay=Combined, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=KO_vs_WT, selected_assay=Combined, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-table_assay: All|LCMS_Pos|GCMS -> All"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "complete"),
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<hclust>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "<NULL>",
            current_state_name = "lipid_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_lipids_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = c(da_lipids_long = "24x23 276ea5100bd6ff77"),
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "45",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/lipid-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/lipid_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/lipid_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.tsv` = "074f019210bb05e4db534de339e88fba7e321a2f2e9b0d50bb503d7d5bb8eb98",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.tsv` = "8597e926b23cddc4514ddec594433c871289662093c37c5f845074adafcdc8e3",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "3ef61cc44f484ab6f9f7c6a1fef40960ffe082360f5bed58ddbff740cdbcdd23",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.tsv` = "c10aa8e923a83d86c58b8e06db8c79d9917eeb185799fdb802622f5286c9e7e8",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>",
            `source/lipid_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    select = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 09ef8be32102d8d33317cd28d2a3c85a2bd212788a83e75fad2a92fcefb21a22",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateLipidDAVolcanoStatic(selected_contrast=RES_vs_WT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=FALSE, scale_data=column, color_scheme=viridis, show_lipid_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 3x7 <HEX>"),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<NULL>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "<NULL>",
            current_state_name = "lipid_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_lipids_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_k = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 lipids):\nL_UP, L_TIED\n\nCluster 2 (1 lipids):\nL_DOWN\n\nCluster 3 (1 lipids):\nL_METALESS",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 09ef8be32102d8d33317cd28d2a3c85a2bd212788a83e75fad2a92fcefb21a22",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "generateLipidDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_lipid_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=k_clusters, n_clusters=3, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<hclust>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "1|2|1|3",
            current_state_name = "lipid_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_lipids_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_height = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 lipids):\nL_UP, L_TIED\n\nCluster 2 (1 lipids):\nL_DOWN\n\nCluster 3 (1 lipids):\nL_METALESS",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 09ef8be32102d8d33317cd28d2a3c85a2bd212788a83e75fad2a92fcefb21a22",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "generateLipidDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_lipid_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=height_cutoff, n_clusters=3, cut_height=5, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_dynamic = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 09ef8be32102d8d33317cd28d2a3c85a2bd212788a83e75fad2a92fcefb21a22",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "generateLipidDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_lipid_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=dynamic, n_clusters=3, cut_height=5, min_cluster_size=2) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<hclust>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "<NULL>",
            current_state_name = "lipid_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_lipids_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    save_heatmap = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 09ef8be32102d8d33317cd28d2a3c85a2bd212788a83e75fad2a92fcefb21a22",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "save_heatmap_products(output_dir=<ROOT>/publication_graphs, file_prefix=lipid_RES_vs_WT, params_list$contrast=RES_vs_WT, params_list$top_n=20, params_list$clustering_method=average, params_list$distance_method=manhattan, params_list$cluster_rows=TRUE, params_list$cluster_cols=TRUE, params_list$scaling=column, params_list$color_scheme=viridis, params_list$tree_cut_method=dynamic, params_list$n_clusters=3, params_list$cut_height=5, params_list$min_cluster_size=2) -> <list>",
        notes = "[message 5 -] Heatmap and cluster info saved to publication_graphs/Heatmap",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "45",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/lipid-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/lipid_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/lipid_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.tsv` = "074f019210bb05e4db534de339e88fba7e321a2f2e9b0d50bb503d7d5bb8eb98",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.tsv` = "8597e926b23cddc4514ddec594433c871289662093c37c5f845074adafcdc8e3",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "3ef61cc44f484ab6f9f7c6a1fef40960ffe082360f5bed58ddbff740cdbcdd23",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.tsv` = "c10aa8e923a83d86c58b8e06db8c79d9917eeb185799fdb802622f5286c9e7e8",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmap/lipid_RES_vs_WT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmap/lipid_RES_vs_WT_heatmap.png` = "<present>",
            `publication_graphs/Heatmap/lipid_RES_vs_WT_heatmap_methods.md` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>",
            `source/lipid_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    table_down = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 64baef3392b7ebbaf808699b6313102b0fbd74dfd7094f147e188bc2b23eba26",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 2x7 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_limit = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> ee7f7bb3dab700825030c1c66b9f2db18a8d1343b1c8b80e2cd22332335ad28e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 5x7 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_combined = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "{\"x\":null,\"evals\":[],\"jsHooks\":[],\"deps\":[{\"name\":\"jquery\",\"version\":\"3.6.0\",\"src\":{\"href\":\"jquery-3.6.0\"},\"meta\":null,\"script\":\"jquery-3.6.0.min.js\",\"stylesheet\":null,\"head\":null,\"attachment\":null,\"all_files\":true},{\"name\":\"dt-core\",\"version\":\"1.13.6\",\"src\":{\"href\":\"dt-core-1.13.6\"},\"meta\":null,\"script\":\"js/jquery.dataTables.min.js\",\"stylesheet\":[\"css/jquery.dataTables.min.css\",\"css/jquery.dataTables.extra.css\"],\"head\":null,\"attachment\":null,\"package\":null,\"all_files\":false},{\"name\":\"crosstalk\",\"version\":\"1.2.2\",\"src\":{\"href\":\"crosstalk-1.2.2\"},\"meta\":null,\"script\":\"js/crosstalk.min.js\",\"stylesheet\":\"css/crosstalk.min.css\",\"head\":null,\"attachment\":null,\"all_files\":true}]}",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 0\nSignificant (Q < 0.050): 0 (0.0%)\n  Up-regulated: 0\n  Down-regulated: 0",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable() -> <NULL>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    download = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "{\"x\":null,\"evals\":[],\"jsHooks\":[],\"deps\":[{\"name\":\"jquery\",\"version\":\"3.6.0\",\"src\":{\"href\":\"jquery-3.6.0\"},\"meta\":null,\"script\":\"jquery-3.6.0.min.js\",\"stylesheet\":null,\"head\":null,\"attachment\":null,\"all_files\":true},{\"name\":\"dt-core\",\"version\":\"1.13.6\",\"src\":{\"href\":\"dt-core-1.13.6\"},\"meta\":null,\"script\":\"js/jquery.dataTables.min.js\",\"stylesheet\":[\"css/jquery.dataTables.min.css\",\"css/jquery.dataTables.extra.css\"],\"head\":null,\"attachment\":null,\"package\":null,\"all_files\":false},{\"name\":\"crosstalk\",\"version\":\"1.2.2\",\"src\":{\"href\":\"crosstalk-1.2.2\"},\"meta\":null,\"script\":\"js/crosstalk.min.js\",\"stylesheet\":\"css/crosstalk.min.css\",\"head\":null,\"attachment\":null,\"all_files\":true}]}",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 0\nSignificant (Q < 0.050): 0 (0.0%)\n  Up-regulated: 0\n  Down-regulated: 0",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0),
        download = c(
            file = "lipidomics_da_results_<DATE>.csv",
            sha256 = "aad6f3f2fa52185a86eed5428dfef9f0c5f0548110add82a790e73486a522be4")))

omicDaW1Frozen[["lipid_no_friendly"]] <- list(
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<error> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "<error> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            volcano_static = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 16x7 <HEX>",
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: groupKO-groupWT|groupRES-groupWT -> groupKO-groupWT",
            "de-heatmap_contrast: groupKO-groupWT|groupRES-groupWT -> groupKO-groupWT",
            "de-table_contrast: groupKO-groupWT|groupRES-groupWT -> groupKO-groupWT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-table_assay: All|LCMS_Pos|GCMS -> All"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "complete"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "TRUE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>",
            da_results_list = "<list: contrasts_results,da_lipids_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "<NULL>"),
        rows = c(da_lipids_long = "24x22 7fddbdb5400fbc8c"),
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "13",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid_da_contrast-da_contrast-<ID>` = "2",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.tsv` = "d2f89fe4d910977451525faf47b8df2eb865e4f9dc672e79010ee282ce6795d2",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.tsv` = "19ab37e5e8eeb250111d0163ec7dd332bf924603ab81ebb5b847f5d734c1e1c8",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "9633c5a76d0bf66dbfb318b6e8d7be9861e5a3c88671c3e54abbc06e9f3912a9",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.tsv` = "85c142271d7aa11a681709e08a4ce344d18605a19dea87e171a19a22685bb5cb",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.xlsx` = "<present>"),
        operations = NULL,
        errors = c(
            "<unhandled> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            "<unhandled> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            "<unhandled> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.",
            "<unhandled> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.")),
    table = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<error> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "<error> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            volcano_static = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = c(
            "<unhandled> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            "<unhandled> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.",
            "<unhandled> Can't subset rows with `results$comparison == selectedContrast | ...`.\nx Logical subscript `results$comparison == selectedContrast | ...` must be size 1 or 24, not 0.")))

omicDaW1Frozen[["metab_refusals"]] <- list(
    load_missing = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] Session file not found: <ROOT>/source/metab_filtered_session_data_latest.rds",
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = character(0),
        states = "743a1d2c832c72b269ed1dc4181d24657bf9c7296c6f9b1c4934328d5583605e",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = character(0),
        operations = NULL,
        errors = character(0)),
    run_no_data = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] No metabolomics data loaded. Please load a filtered session first.",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["metab_mismatch"]] <- list(
    load_mismatch = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL loading_session] Loading filtered session...", "remove loading_session",
            "[error 10 -] Error loading session: Session assay_names do not match current_s4_object assays: session=LCMS_Pos, LCMS_Neg; object=LCMS_Pos, GCMS"),
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/metab_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["metab_reordered"]] <- list(
    load_reordered = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL loading_session] Loading filtered session...", "remove loading_session",
            "[error 10 -] Error loading session: Session assay_names do not match current_s4_object assays: session=GCMS, LCMS_Pos; object=LCMS_Pos, GCMS"),
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/metab_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["metab_no_contrasts"]] <- list(
    run_no_contrasts = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>",
            contrasts_display = "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab.",
            da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] No contrasts defined. Please define contrasts in the design tab.",
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["metab_invalid_preflight"]] <- list(
    run_invalid = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "Missing_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] Metabolomics DA contrast 'groupMISSING-groupWT' references terms absent from the model matrix: groupMISSING",
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["metab_main"]] <- list(
    start = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Waiting for data.\nClick 'Load Filtered Session' to begin.",
            da_summary_stats = "<req>", download_da_results = "<download>",
            heatmap_manual_save_warning = "<NULL>", heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "pending"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/metab_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    load = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT", da_results_table = "<req>",
            da_status = "Data loaded. Ready to run analysis.", da_summary_stats = "<req>",
            download_da_results = "<download>", heatmap_manual_save_warning = "<NULL>",
            heatmap_plot = "<req>", volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL loading_session] Loading filtered session...", "remove loading_session",
            "[message 3 -] Session loaded successfully!"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> <none>",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> <none>",
            "de-table_assay: All|LCMS_Pos|GCMS -> <none>", "de-formula_string: <same> -> ~ 0 + group"),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "FALSE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<NULL>",
            current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "metab_correlation_filtered", da_results_list = "<NULL>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/metab_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 13c95e9b7cba61911ce889b5c985e49654e6dad862f963f51ea6aa4728e0d1f8",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateMetabDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=KO_vs_WT, selected_assay=Combined, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=KO_vs_WT, selected_assay=Combined, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>",
            "generateMetabDAVolcanoStatic(selected_contrast=KO_vs_WT, selected_assay=Combined, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=KO_vs_WT, selected_assay=Combined, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-table_assay: All|LCMS_Pos|GCMS -> All"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "complete"),
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<hclust>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "<NULL>",
            current_state_name = "metab_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_metabolites_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = c(da_metabolites_long = "24x23 134e87ea555ea0a0"),
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "45",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/metab-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/metab_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/metab_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "8182dd547c4ecf1bfcd64eede3d2ca7fed28796be30db143270b4e1938b81f4a",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.tsv` = "751ac64ceb5f2d1bb8602c40f16806d39341ff1e4c25cdf5b9e31fb66604cfce",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "cc3b4d6c4323e14b01e7caa7f564239f0dccf2816bee6c3c4d191ad2d90b273b",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.tsv` = "b557e8c011c4f6809ca4d940a91686b44fdc3ed0b2af331174826e78ff3dcb7b",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>",
            `source/metab_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    select = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 7627557da9f8d3b5c631ff4cae841a41ce9951146353660acf9991c9b000f09e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateMetabDAVolcanoStatic(selected_contrast=RES_vs_WT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=FALSE, scale_data=column, color_scheme=viridis, show_metabolite_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 3x7 <HEX>"),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<NULL>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "<NULL>",
            current_state_name = "metab_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_metabolites_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_k = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 metabolites):\nM_UP, M_TIED\n\nCluster 2 (1 metabolites):\nM_DOWN\n\nCluster 3 (1 metabolites):\nM_METALESS",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 7627557da9f8d3b5c631ff4cae841a41ce9951146353660acf9991c9b000f09e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "generateMetabDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_metabolite_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=k_clusters, n_clusters=3, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<hclust>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "1|2|1|3",
            current_state_name = "metab_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_metabolites_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_height = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 metabolites):\nM_UP, M_TIED\n\nCluster 2 (1 metabolites):\nM_DOWN\n\nCluster 3 (1 metabolites):\nM_METALESS",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 7627557da9f8d3b5c631ff4cae841a41ce9951146353660acf9991c9b000f09e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "generateMetabDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_metabolite_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=height_cutoff, n_clusters=3, cut_height=5, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_dynamic = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 7627557da9f8d3b5c631ff4cae841a41ce9951146353660acf9991c9b000f09e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "generateMetabDAHeatmap(selected_contrast=RES_vs_WT, selected_assay=GCMS, top_n=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_metabolite_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=dynamic, n_clusters=3, cut_height=5, min_cluster_size=2) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_col_clusters = "<hclust>",
            current_heatmap_plot = "<S4 Heatmap>", current_row_clusters = "<NULL>",
            current_state_name = "metab_correlation_filtered",
            da_results_list = "<list: contrasts_results,da_metabolites_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "~ 0 + group"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    save_heatmap = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 7627557da9f8d3b5c631ff4cae841a41ce9951146353660acf9991c9b000f09e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "save_heatmap_products(output_dir=<ROOT>/publication_graphs, file_prefix=metab_RES_vs_WT, params_list$contrast=RES_vs_WT, params_list$top_n=20, params_list$clustering_method=average, params_list$distance_method=manhattan, params_list$cluster_rows=TRUE, params_list$cluster_cols=TRUE, params_list$scaling=column, params_list$color_scheme=viridis, params_list$tree_cut_method=dynamic, params_list$n_clusters=3, params_list$cut_height=5, params_list$min_cluster_size=2) -> <list>",
        notes = "[message 5 -] Heatmap and cluster info saved to publication_graphs/Heatmap",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "5",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "45",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/metab-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/metab_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/metab_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "8182dd547c4ecf1bfcd64eede3d2ca7fed28796be30db143270b4e1938b81f4a",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.tsv` = "751ac64ceb5f2d1bb8602c40f16806d39341ff1e4c25cdf5b9e31fb66604cfce",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "cc3b4d6c4323e14b01e7caa7f564239f0dccf2816bee6c3c4d191ad2d90b273b",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.tsv` = "b557e8c011c4f6809ca4d940a91686b44fdc3ed0b2af331174826e78ff3dcb7b",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmap/metab_RES_vs_WT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmap/metab_RES_vs_WT_heatmap.png` = "<present>",
            `publication_graphs/Heatmap/metab_RES_vs_WT_heatmap_methods.md` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>",
            `source/metab_filtered_session_data_latest.rds` = "<present>"),
        operations = NULL,
        errors = character(0)),
    table_down = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> b463ab571c88721c6cf5936e17a97f093187530f84b2d6e6e1a2e3842648b834",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 2x7 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_limit = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> f820370872e4eac44393e3b560f7615c4d5ab9002ff8c41ade4fedb859284285",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 5x7 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_combined = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> f820370872e4eac44393e3b560f7615c4d5ab9002ff8c41ade4fedb859284285",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 5x7 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    download = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.",
            contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> f820370872e4eac44393e3b560f7615c4d5ab9002ff8c41ade4fedb859284285",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0),
        download = c(
            file = "metabolomics_da_results_<DATE>.csv",
            sha256 = "7adfe999098a7dd429d3c5f97a522f719ca4e67c02d41b8a6837987c7acce284")))

omicDaW1Frozen[["metab_no_friendly"]] <- list(
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 13c95e9b7cba61911ce889b5c985e49654e6dad862f963f51ea6aa4728e0d1f8",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            volcano_static = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found"),
        calls = c(
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 16x7 <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: groupKO-groupWT|groupRES-groupWT -> groupKO-groupWT",
            "de-heatmap_contrast: groupKO-groupWT|groupRES-groupWT -> groupKO-groupWT",
            "de-table_contrast: groupKO-groupWT|groupRES-groupWT -> groupKO-groupWT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-table_assay: All|LCMS_Pos|GCMS -> All"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "complete"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "TRUE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>",
            da_results_list = "<list: contrasts_results,da_metabolites_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "<NULL>"),
        rows = c(da_metabolites_long = "24x22 aa6b404d253370a5"),
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "13",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab_da_contrast-da_contrast-<ID>` = "2",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "3e3c6fc5d8d3c979b055d3c20ac40196fa1d5e9559e3d0094564d616cb421712",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.tsv` = "a717b83eaa50cbfe1da2c6f173204e476507406d34669f63e690693d7dd6cbcb",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "5aa0397dfea01933f4c75f1baa67c7d456f3245c414e0e531a53d6a6927644a8",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.tsv` = "5197cc3c9a90f81bba25e5cf83feb93eae5ddfa6091efc97876bccea812a1bef",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.xlsx` = "<present>"),
        operations = NULL,
        errors = c(
            "<unhandled> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            "<unhandled> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found")),
    table = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 13c95e9b7cba61911ce889b5c985e49654e6dad862f963f51ea6aa4728e0d1f8",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found",
            volcano_static = "<error> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 12x7 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = "<unhandled> i In argument: `comparison == comparison_to_search | friendly_name == comparison_to_search`.\nCaused by error:\n! object 'friendly_name' not found"))

omicDaW1Frozen[["trx_refusals"]] <- list(
    run_no_data = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI", da_results_table = "<req>",
            da_status = "Ready. The DA runs on the normalised (RUVs) state.", da_summary_stats = "<req>",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] Run the normalisation first: the DA runs on the normalised (RUVs) state.",
        updates = character(0),
        tab_status = c(
            quality_control = "complete", normalization = "pending", differential_analysis = "disabled",
            session_summary = "disabled"),
        history = c("initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete"),
        states = "d6b0c4f539ac28e638197e985adea5895477821122c9dee997f57ebc72f6859c",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "13",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = NULL,
        errors = character(0)),
    norm = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI", da_results_table = "<req>",
            da_status = "Ready. The DA runs on the normalised (RUVs) state.", da_summary_stats = "<req>",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<req>"),
        calls = "DT::datatable(rownames=FALSE) -> data.frame 6x6 <HEX>",
        notes = "[message 5 -] Normalization complete! Proceeding to DE analysis.",
        updates = character(0),
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "disabled",
            session_summary = "disabled"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"),
        states = "bbcc8e91c9921cafbad19de435d2542e4dcbcaccf24addea238253fb56938ddc",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = NULL,
        errors = character(0)),
    run_no_contrasts = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI", da_results_table = "<req>",
            da_status = "Ready. The DA runs on the normalised (RUVs) state.", da_summary_stats = "<req>",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<req>"),
        calls = character(0),
        notes = "[error 5 -] No contrasts defined. Please define contrasts in the design tab.",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    run_invalid = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI", da_results_table = "<req>",
            da_status = "Ready. The DA runs on the normalised (RUVs) state.", da_summary_stats = "<req>",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[error 10 -] Analysis error: the contrast 'groupMissing-groupRPMI' is not estimable: object 'groupMissing' not found"),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "25",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["trx_main"]] <- list(
    norm = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI", da_results_table = "<req>",
            da_status = "Ready. The DA runs on the normalised (RUVs) state.", da_summary_stats = "<req>",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<req>"),
        calls = "DT::datatable(rownames=FALSE) -> data.frame 6x6 <HEX>",
        notes = "[message 5 -] Normalization complete! Proceeding to DE analysis.",
        updates = character(0),
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "disabled",
            session_summary = "disabled"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"),
        states = "bbcc8e91c9921cafbad19de435d2542e4dcbcaccf24addea238253fb56938ddc",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = NULL,
        errors = character(0)),
    run = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "<text> 5f7c2b7dc125a9b9d579fe8f2f20aef07b3767effca4ecaa5c10ab2be6ece1f8",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 187\nSignificant (FDR (BH) < 0.050): 23 (12.3%)\n  Up-regulated: 9\n  Down-regulated: 14",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 23x6 <HEX>",
            "trxDaStaticVolcano(selected_contrast=Sera_vs_RPMI, selected_assay=counts, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 23x6 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-heatmap_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-table_contrast: Sera_vs_RPMI -> Sera_vs_RPMI", "de-volcano_assay: Combined|counts -> counts",
            "de-heatmap_assay: Combined|counts -> Combined", "de-table_assay: All -> All"),
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "complete",
            session_summary = "disabled"),
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>",
            da_results_list = "<list: contrasts_results,da_genes_long,per_assay_results,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts>",
            formula_from_s4 = "<NULL>"),
        rows = c(da_genes_long = "187x11 121582784dc95c0a"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "28",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "6",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx-da-fit-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/da_genes/de_results_Sera_vs_RPMI.csv` = "c11158b3af965955f4284fad607de9462097852a94c2b074ffab0e78c0a4af6b",
            `results/transcriptomics/da_genes/ruv_optimization_summary.txt` = "<present>",
            `results/transcriptomics/da_genes/top_de_genes.tsv` = "10aa0f7b3fa91cdc9137a9e4f11dca3a1413f4bd4e236542637c471bf87daa21",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = NULL,
        errors = character(0)),
    select = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "<text> f1a394434d8e77e898a3df94550468f1a8718ea1c94756c3f519d5185a0107fb",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 187\nSignificant (FDR (BH) < 0.050): 23 (12.3%)\n  Up-regulated: 9\n  Down-regulated: 14",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "trxDaStaticVolcano(selected_contrast=Sera_vs_RPMI, selected_assay=counts, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 9x6 <HEX>"),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_down = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "<text> 65049c87c2e5f012c687d9ea09832323f15265865fc13639756c0e9f2b9d6e6e",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 187\nSignificant (FDR (BH) < 0.050): 23 (12.3%)\n  Up-regulated: 9\n  Down-regulated: 14",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 14x6 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_limit = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "<text> 6e05b7d95d04ce037425fc036dea45ce7ce912f07224dfe2b41a977b2b063365",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 187\nSignificant (FDR (BH) < 0.050): 23 (12.3%)\n  Up-regulated: 9\n  Down-regulated: 14",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 5x6 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_combined = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "{\"x\":null,\"evals\":[],\"jsHooks\":[],\"deps\":[{\"name\":\"jquery\",\"version\":\"3.6.0\",\"src\":{\"href\":\"jquery-3.6.0\"},\"meta\":null,\"script\":\"jquery-3.6.0.min.js\",\"stylesheet\":null,\"head\":null,\"attachment\":null,\"all_files\":true},{\"name\":\"dt-core\",\"version\":\"1.13.6\",\"src\":{\"href\":\"dt-core-1.13.6\"},\"meta\":null,\"script\":\"js/jquery.dataTables.min.js\",\"stylesheet\":[\"css/jquery.dataTables.min.css\",\"css/jquery.dataTables.extra.css\"],\"head\":null,\"attachment\":null,\"package\":null,\"all_files\":false},{\"name\":\"crosstalk\",\"version\":\"1.2.2\",\"src\":{\"href\":\"crosstalk-1.2.2\"},\"meta\":null,\"script\":\"js/crosstalk.min.js\",\"stylesheet\":\"css/crosstalk.min.css\",\"head\":null,\"attachment\":null,\"all_files\":true}]}",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 0\nSignificant (FDR (BH) < 0.050): 0 (0.0%)\n  Up-regulated: 0\n  Down-regulated: 0",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = "DT::datatable() -> <NULL>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    download = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "{\"x\":null,\"evals\":[],\"jsHooks\":[],\"deps\":[{\"name\":\"jquery\",\"version\":\"3.6.0\",\"src\":{\"href\":\"jquery-3.6.0\"},\"meta\":null,\"script\":\"jquery-3.6.0.min.js\",\"stylesheet\":null,\"head\":null,\"attachment\":null,\"all_files\":true},{\"name\":\"dt-core\",\"version\":\"1.13.6\",\"src\":{\"href\":\"dt-core-1.13.6\"},\"meta\":null,\"script\":\"js/jquery.dataTables.min.js\",\"stylesheet\":[\"css/jquery.dataTables.min.css\",\"css/jquery.dataTables.extra.css\"],\"head\":null,\"attachment\":null,\"package\":null,\"all_files\":false},{\"name\":\"crosstalk\",\"version\":\"1.2.2\",\"src\":{\"href\":\"crosstalk-1.2.2\"},\"meta\":null,\"script\":\"js/crosstalk.min.js\",\"stylesheet\":\"css/crosstalk.min.css\",\"head\":null,\"attachment\":null,\"all_files\":true}]}",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 0\nSignificant (FDR (BH) < 0.050): 0 (0.0%)\n  Up-regulated: 0\n  Down-regulated: 0",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0),
        download = c(
            file = "transcriptomics_da_results_<DATE>.csv",
            sha256 = "d9e7dd3dfa432147d91b1d36d3ce60ecf26ee6de44ade9f217ae792e4ee346c4")))

omicDaW1Frozen[["prot_refusals"]] <- list(
    load_missing = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "No contrasts available.\nComplete normalization and\ncorrelation filtering first.",
            da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = "[error 10 -] Error loading session: No exported session data found. Please export session data from the Normalization tab first.",
        updates = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = "initial",
        states = "6fd42083fde9c0d24d05ab41df9d5bc4835696e0abd42b463a04efc255c3d1c0",
        tab = c(
            analysis_complete = "FALSE", contrasts_available = "<NULL>", current_col_clusters = "<NULL>",
            current_row_clusters = "<NULL>", current_s4_object = "<NULL>", da_results_list = "<NULL>",
            formula_from_s4 = "<NULL>", volcano_plot_identity = "<NULL>"),
        rows = "<NULL>",
        files = c(`<store entries>..workflow-entities/entities/dia-state-<ID>` = "2"),
        operations = NULL,
        errors = character(0)),
    run_no_data = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "No contrasts available.\nComplete normalization and\ncorrelation filtering first.",
            da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = "[error 10 -] Error in DA analysis: no filtered protein data is loaded.",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["prot_load"]] <- list(
    load = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupB-groupA", da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = c(
            "[message 5 -] UniProt annotations loaded: 6 protein annotations available for enrichment analysis",
            "[message 10 -] Session loaded successfully!\n\nData Summary:\n- Proteins: 6\n- Samples: 4\n- Contrasts: 1\n- State: correlation_filtered\n- Export time: <TIME>\n\nReady for differential abundance analysis."),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>",
            "differential_expression-formula_string: <same> -> ~ 0 + group"),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("raw_data_s4", "correlation_filtered"),
        states = "482262ee3c967fcc3cb47199cafeab8db582af062118fd3327c6a4156d7f0900",
        tab = c(
            analysis_complete = "FALSE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<NULL>", formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "10",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "3",
            `source/.dia-native-legacy-<HEX>/legacy-protein-qc/ruv_normalised_results_cln_with_replicates.tsv` = "00600417dee1d305e8e7ab21e4fc6a3d52b6ae698d3b63e18a6ea21c00e4cb10",
            `source/.dia-native-legacy-<HEX>/legacy-source/accession_cleanup_results.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/fasta_metadata.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/mixed_species_analysis.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/qc_params.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/ruv_optimization_results.RDS` = "<present>",
            `source/accession_cleanup_results.RDS` = "<present>", `source/fasta_metadata.RDS` = "<present>",
            `source/filtered_session_data_latest.rds` = "<present>",
            `source/mixed_species_analysis.RDS` = "<present>", `source/qc_params.RDS` = "<present>",
            `source/ruv_optimization_results.RDS` = "<present>", `source/uniprot_dat_cln.RDS` = "<present>"),
        operations = NULL,
        errors = character(0)),
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupB-groupA", da_results_table = "<req>",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>",
            "differential_expression-formula_string: <same> -> ~ 0 + group"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = "<unchanged>",
        states = "b2aa2a422f89ad5b29d9de380489185368d3f3787891e4c76cf139e8c5cdfe8c",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<manifest> schema,protein_id_column,sample_id,contrasts,parameters,model 1cdb0864e51a3fc1",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "2",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-state-<ID>` = "21",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>.source/.dia-native-legacy-<ID>/.workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "180c330673a89c3d58f1e830639d9011ddfc16f96aa8ee5ebdcbfc489d85c7c0",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-protein-qc/ruv_normalised_results_cln_with_replicates.tsv` = "00600417dee1d305e8e7ab21e4fc6a3d52b6ae698d3b63e18a6ea21c00e4cb10",
            `source/.dia-native-legacy-<HEX>/legacy-source/accession_cleanup_results.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/fasta_metadata.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/mixed_species_analysis.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/qc_params.RDS` = "<present>",
            `source/.dia-native-legacy-<HEX>/legacy-source/ruv_optimization_results.RDS` = "<present>",
            `source/accession_cleanup_results.RDS` = "<present>", `source/fasta_metadata.RDS` = "<present>",
            `source/filtered_session_data_latest.rds` = "<present>",
            `source/mixed_species_analysis.RDS` = "<present>", `source/qc_params.RDS` = "<present>",
            `source/ruv_optimization_results.RDS` = "<present>", `source/uniprot_dat_cln.RDS` = "<present>"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["prot_main"]] <- list(
    start = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "No contrasts available.\nComplete normalization and\ncorrelation filtering first.",
            da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = "3a69181ec204208698f47019efc2cf574bfec461bf61ff3effd0820746d01416",
        tab = c(
            analysis_complete = "FALSE", contrasts_available = "<NULL>", current_col_clusters = "<NULL>",
            current_row_clusters = "<NULL>", current_s4_object = "<NULL>", da_results_list = "<NULL>",
            formula_from_s4 = "<NULL>", volcano_plot_identity = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)),
    open = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupB-groupA", da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = character(0),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupB-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupB-groupA -> <none>",
            "differential_expression-table_contrast: groupB-groupA -> <none>"),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "FALSE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<NULL>", formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupB-groupA", da_results_table = "<req>",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-volcano_contrast: B_vs_A -> <none>",
            "differential_expression-heatmap_contrast: B_vs_A -> <none>",
            "differential_expression-table_contrast: B_vs_A -> <none>"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = "<unchanged>",
        states = "ee10cc7a21277b1f1ad516b65868f9ae736e5e28c90ebb69aec239da357a02a8",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<manifest> schema,protein_id_column,sample_id,contrasts,parameters,model 1cdb0864e51a3fc1",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    select = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 4\n----------------------------------\n\nCluster 1 (1 proteins):\nP1\n\nCluster 2 (1 proteins):\nP2\n\nCluster 3 (1 proteins):\nP3\n\nCluster 4 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = c(
            "generateProtDAVolcanoStatic(selected_contrast=B_vs_A, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=10) -> ggplot <HEX>",
            "generateProtDAHeatmap(selected_contrast=B_vs_A, top_n_genes=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=FALSE, scale_data=column, color_scheme=viridis, show_gene_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=k_clusters, n_clusters=4, cut_height=0.5, min_cluster_size=3, show_dendrogram=TRUE) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE) -> data.frame 3x4 <HEX>"),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "1|2|3|4",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "18",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    volcano_unticked = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 4\n----------------------------------\n\nCluster 1 (1 proteins):\nP1\n\nCluster 2 (1 proteins):\nP2\n\nCluster 3 (1 proteins):\nP3\n\nCluster 4 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    volcano_ticked = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 4\n----------------------------------\n\nCluster 1 (1 proteins):\nP1\n\nCluster 2 (1 proteins):\nP2\n\nCluster 3 (1 proteins):\nP3\n\nCluster 4 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    volcano_labels = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 4\n----------------------------------\n\nCluster 1 (1 proteins):\nP1\n\nCluster 2 (1 proteins):\nP2\n\nCluster 3 (1 proteins):\nP3\n\nCluster 4 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "generateProtDAVolcanoStatic(selected_contrast=B_vs_A, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=FALSE, n_labels=3) -> ggplot <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    tree_k = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 proteins):\nP1, P2\n\nCluster 2 (1 proteins):\nP3\n\nCluster 3 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "generateProtDAHeatmap(selected_contrast=B_vs_A, top_n_genes=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_gene_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=k_clusters, n_clusters=3, cut_height=0.5, min_cluster_size=3, show_dendrogram=TRUE) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "1|1|2|3",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "2",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "19",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    tree_height = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 2\n----------------------------------\n\nCluster 1 (3 proteins):\nP1, P2, P3\n\nCluster 2 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "generateProtDAHeatmap(selected_contrast=B_vs_A, top_n_genes=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_gene_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=height_cutoff, n_clusters=3, cut_height=5, min_cluster_size=3, show_dendrogram=TRUE) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "1|1|1|2",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "3",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "20",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    tree_dynamic = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 4\n----------------------------------\n\nCluster 1 (1 proteins):\nP1\n\nCluster 2 (1 proteins):\nP2\n\nCluster 3 (1 proteins):\nP3\n\nCluster 4 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "generateProtDAHeatmap(selected_contrast=B_vs_A, top_n_genes=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_gene_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=dynamic, n_clusters=3, cut_height=5, min_cluster_size=2, show_dendrogram=TRUE) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "1|2|3|4",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "21",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    tree_none = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "generateProtDAHeatmap(selected_contrast=B_vs_A, top_n_genes=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_gene_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=3, cut_height=5, min_cluster_size=2, show_dendrogram=TRUE) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "5",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "22",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    dendrogram_off = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 proteins):\nP1, P2\n\nCluster 2 (1 proteins):\nP3\n\nCluster 3 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "generateProtDAHeatmap(selected_contrast=B_vs_A, top_n_genes=20, clustering_method=average, distance_method=manhattan, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=column, color_scheme=viridis, show_gene_names=TRUE, da_q_val_thresh=0.05, tree_cut_method=k_clusters, n_clusters=3, cut_height=5, min_cluster_size=2, show_dendrogram=FALSE) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "groupB-groupA",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "1|1|2|3",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "6",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    save_heatmap = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 proteins):\nP1, P2\n\nCluster 2 (1 proteins):\nP3\n\nCluster 3 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "save_heatmap_products(output_dir=<ROOT>/publication_graphs, file_prefix=prot_B_vs_A, params_list$contrast=B_vs_A, params_list$top_n=20, params_list$clustering_method=average, params_list$distance_method=manhattan, params_list$cluster_rows=TRUE, params_list$cluster_cols=TRUE, params_list$scaling=column, params_list$color_scheme=viridis, params_list$tree_cut_method=k_clusters, params_list$n_clusters=3, params_list$cut_height=5, params_list$min_cluster_size=2) -> <list>",
        notes = "[message 5 -] Heatmap and cluster info saved to publication_graphs/Heatmap",
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-<ID>` = "1",
            `<store entries>..workflow-entities/entities/da-selection-product-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-heatmap-<ID>` = "6",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "23",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupB-groupA_long_annot.tsv` = "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6",
            `da_output/da_proteins_groupB-groupA_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmap/prot_B_vs_A_clusters.csv` = "297557cb0dfdc9797f9f058e6828e3abaa3bf575fe716b3801f7fd727a4d6c81",
            `publication_graphs/Heatmap/prot_B_vs_A_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmap/prot_B_vs_A_heatmap.png` = "<present>",
            `publication_graphs/Heatmap/prot_B_vs_A_heatmap_methods.md` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/B_vs_A.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "8a7cae2958352a9afa06dadd17250eacf2c56a50583e8b0b1c17146d9991c109",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupB-groupA.png` = "<present>"),
        operations = NULL,
        errors = character(0)),
    table_down = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 proteins):\nP1, P2\n\nCluster 2 (1 proteins):\nP3\n\nCluster 3 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 29a63a244e14e7331ce5695853ad9adad0ed7877664af9802aea8f8dc4173486",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE) -> data.frame 1x4 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    table_limit = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 proteins):\nP1, P2\n\nCluster 2 (1 proteins):\nP3\n\nCluster 3 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 7bc710ab6798d25022839efcad87239fe1823853790bd0fdfc0027c8fd3ac5c8",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = "DT::datatable(extensions=Buttons, rownames=FALSE) -> data.frame 5x4 <HEX>",
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)),
    download = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "Total Clusters: 3\n----------------------------------\n\nCluster 1 (2 proteins):\nP1, P2\n\nCluster 2 (1 proteins):\nP3\n\nCluster 3 (1 proteins):\nP6",
            contrasts_display = "groupB-groupA",
            da_results_table = "<text> 7bc710ab6798d25022839efcad87239fe1823853790bd0fdfc0027c8fd3ac5c8",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "Total genes: 6\nSignificant (q < 0.050): 4\nUp-regulated: 3\nDown-regulated: 1\nFold-change cutoff: 0.00",
            download_da_results = "<download>", heatmap_plot = "<plot> alt=Plot object",
            volcano_glimma = "<ui> 6c24457f4f3ce8d348007280dfd26da93e5791de2d3131cecd367fad73baee06",
            volcano_plot_static = "<plot> alt=Plot object"),
        calls = character(0),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0),
        download = c(
            file = "DA_results_<DATE>.zip",
            sha256 = "48ae2a71a12d601ce240dceb9160e6a184a5938ca69d24578091ab8b94bc1bcb")))

omicDaW1Frozen[["prot_invalid_preflight"]] <- list(
    open = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupZ-groupA", da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = character(0),
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupZ-groupA -> <none>",
            "differential_expression-heatmap_contrast: groupZ-groupA -> <none>",
            "differential_expression-table_contrast: groupZ-groupA -> <none>"),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = "3a69181ec204208698f47019efc2cf574bfec461bf61ff3effd0820746d01416",
        tab = c(
            analysis_complete = "FALSE", contrasts_available = "groupZ-groupA",
            current_col_clusters = "<NULL>", current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<NULL>", formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)),
    run_invalid = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupZ-groupA", da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[error 10 -] Error in DA analysis: DA contrast 'groupZ-groupA' references terms absent from the model matrix: groupZ",
            "remove da_working"),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["prot_no_contrasts"]] <- list(
    open = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupA-groupB", da_results_table = "<req>",
            da_status = "[WAITING] Waiting for analysis...\nClick 'Run DA Analysis' to start",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = "[warning 5 -] No user-specified contrasts found. Auto-generating all possible contrasts.",
        updates = c(
            "differential_expression-formula_string: <same> -> ~ 0 + group",
            "differential_expression-volcano_contrast: groupA-groupB -> <none>",
            "differential_expression-heatmap_contrast: groupA-groupB -> <none>",
            "differential_expression-table_contrast: groupA-groupB -> <none>"),
        tab_status = c(
            normalization = "complete", differential_expression = "pending", enrichment_analysis = "disabled"),
        history = c("initial", "correlation_filtered"),
        states = "3a69181ec204208698f47019efc2cf574bfec461bf61ff3effd0820746d01416",
        tab = c(
            analysis_complete = "FALSE", contrasts_available = "<array>", current_col_clusters = "<NULL>",
            current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<NULL>", formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "4",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "3"),
        operations = NULL,
        errors = character(0)),
    run_no_contrasts = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_plot", "volcano_glimma", "volcano_plot_static"),
        values = c(
            cluster_summary = "No clusters defined. Enable clustering and tree cutting on the heatmap.NULL",
            contrasts_display = "groupA-groupB", da_results_table = "<req>",
            da_status = "[OK] Analysis Complete\n Contrasts analyzed: 1\n Results available in all tabs",
            da_summary_stats = "<req>", download_da_results = "<download>", heatmap_plot = "<req>",
            volcano_glimma = "<req>", volcano_plot_static = "<req>"),
        calls = character(0),
        notes = c(
            "[default NULL da_working] Running differential abundance analysis...",
            "[message 5 -] Differential abundance analysis completed successfully!"),
        updates = c(
            "differential_expression-volcano_contrast: A_vs_B -> <none>",
            "differential_expression-heatmap_contrast: A_vs_B -> <none>",
            "differential_expression-table_contrast: A_vs_B -> <none>"),
        tab_status = c(
            normalization = "complete", differential_expression = "complete", enrichment_analysis = "pending",
            differential_abundance = "complete"),
        history = "<unchanged>",
        states = "d943ed2d07378a01099b974ae1be0e283169a5fac0dd5cf54479175659b63481",
        tab = c(
            analysis_complete = "TRUE", contrasts_available = "<array>", current_col_clusters = "<NULL>",
            current_row_clusters = "<NULL>",
            current_s4_object = "<list: schema,owner,kind,identity,payload,metadata,parent_id,revision,entity_id>",
            da_results_list = "<list: schema,source,protein_id_column,sample_id,contrasts,results,parameters,model>",
            formula_from_s4 = "~ 0 + group", volcano_plot_identity = "<NULL>"),
        rows = "<manifest> schema,protein_id_column,sample_id,contrasts,parameters,model 6e817ac35b5d05b7",
        files = c(
            `<store entries>..workflow-entities/entities/DA-contrast-output-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules-<ID>` = "1",
            `<store entries>..workflow-entities/entities/DA-contrast-output-num_sig_da_molecules_first_go-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-base-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-contrasted-fit-<ID>` = "1",
            `<store entries>..workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>..workflow-entities/entities/dia-table-<ID>` = "8",
            `<store entries>..workflow-entities/entities/dia-volcano-html-<ID>` = "1",
            `da_output/da_proteins_groupA-groupB_long_annot.tsv` = "8b56ef913db712729cc4169dbcc8a98372446c986b6e4401ad17da377cc34847",
            `da_output/da_proteins_groupA-groupB_long_annot.xlsx` = "<present>",
            `publication_graphs/Interactive_Volcano_Plots/A_vs_B.html` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.pdf` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png` = "<present>",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.tab` = "34e7cdf6569b1fd322acdbbf846c9ecaef198bf58b2c6d0434c0d838a580c4e6",
            `publication_graphs/NumSigDaMolecules/da_proteins_num_sig_da_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupA-groupB.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/groupA-groupB.png` = "<present>"),
        operations = NULL,
        errors = character(0)))

omicDaW1Frozen[["lipid_managed"]] <- list(
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 1bca5d30123f39779d9e6238fe985fcedf9631330542910d6850ffc05db0b701",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateLipidDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateLipidDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 16x7 <HEX>",
            "generateLipidDAVolcanoStatic(selected_contrast=KO_vs_WT, selected_assay=Combined, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateLipidDAHeatmap(selected_contrast=KO_vs_WT, selected_assay=Combined, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_lipid_names=FALSE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...",
            "[default NULL lipidomics_action_working] Applying lipidomics step ...",
            "remove lipidomics_action_working",
            "[default NULL lipidomics_action_working] Applying lipidomics step ...",
            "remove lipidomics_action_working",
            "[default NULL lipidomics_action_working] Applying lipidomics step ...",
            "remove lipidomics_action_working",
            "[default NULL lipidomics_action_working] Applying lipidomics step ...",
            "remove lipidomics_action_working", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-table_assay: All|LCMS_Pos|GCMS -> All"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "complete"),
        history = c("initial", "lipid_correlation_filtered"),
        states = "672439f0bdf18a5a8ee3d7579dd29fee0560710d5d3fab6f415ede89199c64f6",
        tab = c(
            analysis_complete = "TRUE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<hclust>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "<NULL>", current_state_name = "<NULL>",
            da_results_list = "<list: da_lipids_long,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts,download_columns,write_download>",
            formula_from_s4 = "<NULL>"),
        rows = c(da_lipids_long = "24x17 731129c5141db751"),
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "14",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.tsv` = "074f019210bb05e4db534de339e88fba7e321a2f2e9b0d50bb503d7d5bb8eb98",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.tsv` = "8597e926b23cddc4514ddec594433c871289662093c37c5f845074adafcdc8e3",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "3ef61cc44f484ab6f9f7c6a1fef40960ffe082360f5bed58ddbff740cdbcdd23",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.tsv` = "c10aa8e923a83d86c58b8e06db8c79d9917eeb185799fdb802622f5286c9e7e8",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>"),
        operations = c(
            "lipid_da_contrast", "lipid_da_contrast", "lipid_da_tab_results", "lipid_da_publication"),
        errors = character(0)),
    select = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> ee7f7bb3dab700825030c1c66b9f2db18a8d1343b1c8b80e2cd22332335ad28e",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total lipids: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateLipidDAVolcanoStatic(selected_contrast=RES_vs_WT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 6x7 <HEX>"),
        notes = c(
            "[default NULL lipidomics_action_working] Applying lipidomics step ...",
            "remove lipidomics_action_working"),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/lipid_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-volcano-html-<ID>` = "1",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.tsv` = "074f019210bb05e4db534de339e88fba7e321a2f2e9b0d50bb503d7d5bb8eb98",
            `da_output/de_gcms_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.tsv` = "8597e926b23cddc4514ddec594433c871289662093c37c5f845074adafcdc8e3",
            `da_output/de_gcms_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.tsv` = "3ef61cc44f484ab6f9f7c6a1fef40960ffe082360f5bed58ddbff740cdbcdd23",
            `da_output/de_posmode_lipids_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.tsv` = "c10aa8e923a83d86c58b8e06db8c79d9917eeb185799fdb802622f5286c9e7e8",
            `da_output/de_posmode_lipids_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/lipids_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>"),
        operations = "lipid_volcano_page",
        errors = character(0)))

omicDaW1Frozen[["metab_managed"]] <- list(
    run = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> 13c95e9b7cba61911ce889b5c985e49654e6dad862f963f51ea6aa4728e0d1f8",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateMetabDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=LCMS_Pos, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=groupKO-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupKO-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "generateMetabDAVolcanoStatic(selected_contrast=groupRES-groupWT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=groupRES-groupWT, selected_assay=GCMS, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 16x7 <HEX>",
            "generateMetabDAVolcanoStatic(selected_contrast=KO_vs_WT, selected_assay=Combined, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "generateMetabDAHeatmap(selected_contrast=KO_vs_WT, selected_assay=Combined, top_n=50, clustering_method=ward.D2, distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, show_metabolite_names=FALSE, da_q_val_thresh=0.05, tree_cut_method=none, n_clusters=4, cut_height=10, min_cluster_size=5) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 8x7 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...",
            "[default NULL metabolomics_action_working] Applying metabolomics step ...",
            "remove metabolomics_action_working",
            "[default NULL metabolomics_action_working] Applying metabolomics step ...",
            "remove metabolomics_action_working",
            "[default NULL metabolomics_action_working] Applying metabolomics step ...",
            "remove metabolomics_action_working",
            "[default NULL metabolomics_action_working] Applying metabolomics step ...",
            "remove metabolomics_action_working", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-heatmap_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-table_contrast: KO_vs_WT|RES_vs_WT -> KO_vs_WT",
            "de-volcano_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-heatmap_assay: Combined|LCMS_Pos|GCMS -> Combined",
            "de-table_assay: All|LCMS_Pos|GCMS -> All"),
        tab_status = c(
            setup_import = "complete", design_matrix = "complete", quality_control = "complete",
            normalization = "complete", differential_analysis = "complete"),
        history = c("initial", "metab_correlation_filtered"),
        states = "6c9129cd1741dcb566b9899e1e1abc5e14149d4d16741cf1f4844656c1a0e73f",
        tab = c(
            analysis_complete = "TRUE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<hclust>", current_heatmap_plot = "<S4 Heatmap>",
            current_row_clusters = "<NULL>", current_state_name = "<NULL>",
            da_results_list = "<list: da_metabolites_long,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts,download_columns,write_download>",
            formula_from_s4 = "<NULL>"),
        rows = c(da_metabolites_long = "24x17 088acb43fc1f0c6c"),
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "14",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "8182dd547c4ecf1bfcd64eede3d2ca7fed28796be30db143270b4e1938b81f4a",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.tsv` = "751ac64ceb5f2d1bb8602c40f16806d39341ff1e4c25cdf5b9e31fb66604cfce",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "cc3b4d6c4323e14b01e7caa7f564239f0dccf2816bee6c3c4d191ad2d90b273b",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.tsv` = "b557e8c011c4f6809ca4d940a91686b44fdc3ed0b2af331174826e78ff3dcb7b",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>"),
        operations = c(
            "metab_da_contrast", "metab_da_contrast", "metab_da_tab_results", "metab_da_publication"),
        errors = character(0)),
    select = list(
        outputs = c(
            "cluster_summary", "contrasts_display", "da_results_table", "da_status", "da_summary_stats",
            "download_da_results", "heatmap_manual_save_warning", "heatmap_plot", "volcano_glimma",
            "volcano_static"),
        values = c(
            cluster_summary = "<req>", contrasts_display = "KO_vs_WT\nRES_vs_WT",
            da_results_table = "<text> f820370872e4eac44393e3b560f7615c4d5ab9002ff8c41ade4fedb859284285",
            da_status = "Analysis Complete\n\nLCMS_Pos:\n  Up: 5 | Down: 3 | NS: 4\nGCMS:\n  Up: 5 | Down: 3 | NS: 4",
            da_summary_stats = "Total metabolites: 6\nSignificant (Q < 0.050): 4 (66.7%)\n  Up-regulated: 3\n  Down-regulated: 1",
            download_da_results = "<download>",
            heatmap_manual_save_warning = "<ui> b49e16f946860932c7dd705a255211e9526ac5c4a85467ed0a5fd32b6e9e4d48",
            heatmap_plot = "<plot> alt=Plot object", volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "generateMetabDAVolcanoStatic(selected_contrast=RES_vs_WT, selected_assay=GCMS, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 6x7 <HEX>"),
        notes = c(
            "[default NULL metabolomics_action_working] Applying metabolomics step ...",
            "remove metabolomics_action_working"),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = c(
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "15",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "7",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab-da-fit-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab_da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/metab_da_contrast-da_contrast-<ID>` = "2",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-<ID>` = "11",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-figure-asset-<ID>` = "19",
            `<store entries>.source/.dia-native-<ID>/.workflow-entities/entities/omics-volcano-html-<ID>` = "1",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.tsv` = "8182dd547c4ecf1bfcd64eede3d2ca7fed28796be30db143270b4e1938b81f4a",
            `da_output/de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.tsv` = "751ac64ceb5f2d1bb8602c40f16806d39341ff1e4c25cdf5b9e31fb66604cfce",
            `da_output/de_gcms_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.tsv` = "cc3b4d6c4323e14b01e7caa7f564239f0dccf2816bee6c3c4d191ad2d90b273b",
            `da_output/de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx` = "<present>",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.tsv` = "b557e8c011c4f6809ca4d940a91686b44fdc3ed0b2af331174826e78ff3dcb7b",
            `da_output/de_posmode_metabolites_groupRES-groupWT_long_annot.xlsx` = "<present>",
            `publication_graphs/Heatmaps/all_heatmaps_combined.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/gcms_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupKO-groupWT_heatmap.png` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.pdf` = "<present>",
            `publication_graphs/Heatmaps/posmode_groupRES-groupWT_heatmap.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_barplot.png` = "<present>",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.tab` = "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713",
            `publication_graphs/NumSigDeMolecules/metabolites_num_sig_de_molecules.xlsx` = "<present>",
            `publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/gcms_groupRES-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupKO-groupWT.png` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.pdf` = "<present>",
            `publication_graphs/Volcano_Plots/posmode_groupRES-groupWT.png` = "<present>"),
        operations = "metab_volcano_page",
        errors = character(0)))

omicDaW1Frozen[["trx_managed"]] <- list(
    norm = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI", da_results_table = "<req>",
            da_status = "Ready. The DA runs on the normalised (RUVs) state.", da_summary_stats = "<req>",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<req>"),
        calls = "DT::datatable(rownames=FALSE) -> data.frame 6x6 <HEX>",
        notes = c(
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[message 5 -] Normalization complete! Proceeding to DE analysis."),
        updates = character(0),
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "disabled",
            session_summary = "disabled"),
        history = c(
            "initial", "trx_raw_data_s4", "trx_expression_filtered", "trx_qc_complete", "trx_ruvs"),
        states = "bbcc8e91c9921cafbad19de435d2542e4dcbcaccf24addea238253fb56938ddc",
        tab = c(
            analysis_complete = "FALSE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>", da_results_list = "<NULL>", formula_from_s4 = "<NULL>"),
        rows = "<NULL>",
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "36",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "5",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = c(
            "trx_qc_figures", "trx_tmm_logcpm", "trx_ruv_grid", "trx_ruvs", "trx_qc_figures",
            "trx_qc_figures", "trx_ruv_publication", "trx_composite_figure"),
        errors = character(0)),
    run = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "<text> 5f7c2b7dc125a9b9d579fe8f2f20aef07b3767effca4ecaa5c10ab2be6ece1f8",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 187\nSignificant (FDR (BH) < 0.050): 23 (12.3%)\n  Up-regulated: 9\n  Down-regulated: 14",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 23x6 <HEX>",
            "trxDaStaticVolcano(selected_contrast=Sera_vs_RPMI, selected_assay=counts, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 23x6 <HEX>"),
        notes = c(
            "[default NULL da_running] Running differential expression analysis...",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working",
            "[default NULL transcriptomics_action_working] Applying transcriptomics step ...",
            "remove transcriptomics_action_working", "remove da_running",
            "[message 5 -] Differential expression analysis complete!",
            "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)"),
        updates = c(
            "de-volcano_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-heatmap_contrast: Sera_vs_RPMI -> Sera_vs_RPMI",
            "de-table_contrast: Sera_vs_RPMI -> Sera_vs_RPMI", "de-volcano_assay: Combined|counts -> counts",
            "de-heatmap_assay: Combined|counts -> Combined", "de-table_assay: All -> All"),
        tab_status = c(
            quality_control = "complete", normalization = "complete", differential_analysis = "complete",
            session_summary = "disabled"),
        history = "<unchanged>",
        states = "<unchanged>",
        tab = c(
            analysis_complete = "TRUE", assays_available = "<NULL>", contrasts_available = "<NULL>",
            current_col_clusters = "<NULL>", current_heatmap_plot = "<NULL>", current_row_clusters = "<NULL>",
            current_state_name = "<NULL>",
            da_results_list = "<list: da_genes_long,significant_counts,qvalue_warnings,da_q_val_thresh,treat_lfc_cutoff,failed_contrasts,download_columns,write_download>",
            formula_from_s4 = "<NULL>"),
        rows = c(da_genes_long = "187x11 121582784dc95c0a"),
        files = c(
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-state-<ID>` = "44",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/dia-table-<ID>` = "6",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/partitioned-table-<ID>` = "2",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx-da-fit-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_da_contrast-da_contrast-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_expression_filter-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruv_grid-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_ruvs-<ID>` = "1",
            `<store entries>.scripts/transcriptomics/.dia-native-<ID>/.workflow-entities/entities/trx_tmm_logcpm-<ID>` = "1",
            `results/transcriptomics/da_genes/de_results_Sera_vs_RPMI.csv` = "c11158b3af965955f4284fad607de9462097852a94c2b074ffab0e78c0a4af6b",
            `results/transcriptomics/da_genes/ruv_optimization_summary.txt` = "<present>",
            `results/transcriptomics/da_genes/top_de_genes.tsv` = "10aa0f7b3fa91cdc9137a9e4f11dca3a1413f4bd4e236542637c471bf87daa21",
            `results/transcriptomics/gene_qc/ruv_normalised_log_values.tsv` = "83de1dbdd1018161c5ec6c1dac4dc179b0c1d0340fd2fe826d4b3c9b4ae93387",
            `results/transcriptomics/gene_qc/trx_raw_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_raw_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_ruvs_rle.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_pca.png` = "<present>",
            `results/transcriptomics/gene_qc/trx_tmm_rle.png` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/NumSigDeMolecules/summary_barplot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/Volcano_Plots/volcano_plot_Sera_vs_RPMI.png` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_optimization_summary.xlsx` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/canonical_correlation_panels.png` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.pdf` = "<present>",
            `results/transcriptomics/publication_graphs/normalization_comparison_grid.png` = "<present>",
            `scripts/transcriptomics/design_matrix.tab` = "e3b6f40550967e975506c582126b035f5fcc887aae4ffe97e3200f3999a3584e",
            `scripts/transcriptomics/import_sources/counts.tsv` = "2c9e3b6290daf902297e7dab06425b831bfeb6261bd5b5d1819e6fd3e21dde49"),
        operations = c("trx_da_contrast", "trx_da_tab_results", "trx_da_publication"),
        errors = character(0)),
    select = list(
        outputs = c(
            "contrasts_display", "da_results_table", "da_status", "da_summary_stats", "download_da_results",
            "formula_display", "volcano_glimma", "volcano_static"),
        values = c(
            contrasts_display = "Sera_vs_RPMI",
            da_results_table = "<text> 5f7c2b7dc125a9b9d579fe8f2f20aef07b3767effca4ecaa5c10ab2be6ece1f8",
            da_status = "Analysis Complete\n\ncounts:\n  Up: 9 | Down: 14 | NS: 164",
            da_summary_stats = "Total genes: 187\nSignificant (FDR (BH) < 0.050): 23 (12.3%)\n  Up-regulated: 9\n  Down-regulated: 14",
            download_da_results = "<download>",
            formula_display = "Run the normalisation (RUVs) first; the formula is realised from its design.",
            volcano_static = "<plot> alt=Plot object"),
        calls = c(
            "trxDaStaticVolcano(selected_contrast=Sera_vs_RPMI, selected_assay=counts, da_q_val_thresh=0.05, lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>",
            "DT::datatable(extensions=Buttons, rownames=FALSE, filter=top) -> data.frame 187x6 <HEX>"),
        notes = character(0),
        updates = character(0),
        tab_status = "<unchanged>",
        history = "<unchanged>",
        states = "<unchanged>",
        tab = "<unchanged>",
        rows = "<unchanged>",
        files = "<unchanged>",
        operations = character(0),
        errors = character(0)))
# OMIC_DA_FROZEN_END

expectOmicDaW1Case <- function(case) {
    seen <- omicDaW1Case(case)$captures
    frozen <- omicDaW1Frozen[[case]]
    testthat::expect_false(is.null(frozen), info = case)
    testthat::expect_identical(names(seen), names(frozen), info = case)
    for (step in names(frozen)) {
        for (field in union(names(frozen[[step]]), names(seen[[step]]))) {
            testthat::expect_identical(seen[[step]][[field]], frozen[[step]][[field]], info = paste(case, step, field))
        }
    }
}

for (case in c("lipid_refusals", "lipid_mismatch", "lipid_reordered", "lipid_no_contrasts", "lipid_invalid_preflight",
    "lipid_main", "lipid_no_friendly", "metab_refusals", "metab_mismatch", "metab_reordered", "metab_no_contrasts",
    "metab_invalid_preflight", "metab_main", "metab_no_friendly", "trx_refusals", "trx_main", "prot_refusals",
    "prot_load", "prot_main", "prot_invalid_preflight", "prot_no_contrasts")) {
    local({
        name <- case
        test_that(sprintf("(b) in process: %s", name), expectOmicDaW1Case(name))
    })
}

for (case in c("lipid_managed", "metab_managed", "trx_managed")) {
    local({
        name <- case
        test_that(sprintf("(b-managed) %s: the worker's operations in order, the notifications, tab_status and the rows", name),
            expectOmicDaW1Case(name))
    })
}

# ------------------------------------------------------------------------------------------------
# (c) The per-omic differences (design D7), as today's behaviour, from the (b) captures.
# ------------------------------------------------------------------------------------------------

omicDaW1Notes <- function(case, step) omicDaW1Case(case)$captures[[step]]$notes

test_that("(c) the refusal order: metab refuses an invalid preflight before the running notice, lipid after it; no data and no contrasts come first for lipid, metab and trx", {
    running <- "[default NULL da_running] Running differential expression analysis..."
    lipid <- omicDaW1Notes("lipid_invalid_preflight", "run_invalid")
    expect_identical(lipid, c(running, "remove da_running",
        "[error 5 -] Lipidomics DA contrast 'groupMISSING-groupWT' references terms absent from the model matrix: groupMISSING"))
    expect_identical(omicDaW1Notes("metab_invalid_preflight", "run_invalid"),
        "[error 5 -] Metabolomics DA contrast 'groupMISSING-groupWT' references terms absent from the model matrix: groupMISSING")
    # Trx has no app-side preflight: the declaration refuses the contrast after the running notice.
    trx <- omicDaW1Notes("trx_refusals", "run_invalid")
    expect_identical(trx[1:2], c(running, "remove da_running"))
    expect_match(trx[[3]], "^\\[error 10 -\\] Analysis error: the contrast 'groupMissing-groupRPMI' is not estimable")
    expect_identical(omicDaW1Notes("lipid_refusals", "run_no_data"),
        "[error 5 -] No lipidomics data loaded. Please load a filtered session first.")
    expect_identical(omicDaW1Notes("metab_refusals", "run_no_data"),
        "[error 5 -] No metabolomics data loaded. Please load a filtered session first.")
    expect_identical(omicDaW1Notes("trx_refusals", "run_no_data"),
        "[error 5 -] Run the normalisation first: the DA runs on the normalised (RUVs) state.")
    no_contrasts <- "[error 5 -] No contrasts defined. Please define contrasts in the design tab."
    for (case in c("lipid_no_contrasts", "metab_no_contrasts")) expect_identical(omicDaW1Notes(case, "run_no_contrasts"), no_contrasts)
    expect_identical(omicDaW1Notes("trx_refusals", "run_no_contrasts"), no_contrasts)
    # Proteomics' run plugin: its own notice id and texts; an empty contrast table is auto-generated on the tab's opening.
    expect_identical(omicDaW1Notes("prot_refusals", "run_no_data"), "[error 10 -] Error in DA analysis: no filtered protein data is loaded.")
    expect_identical(omicDaW1Notes("prot_invalid_preflight", "run_invalid")[[1L]],
        "[default NULL da_working] Running differential abundance analysis...")
    expect_identical(omicDaW1Notes("prot_no_contrasts", "open"),
        "[warning 5 -] No user-specified contrasts found. Auto-generating all possible contrasts.")
})

test_that("(c) the assay check: lipid compares the session's assay names as a set, metab in order; each with its own text", {
    loaded <- c("[default NULL loading_session] Loading filtered session...", "remove loading_session",
        "[message 3 -] Session loaded successfully!")
    expect_identical(omicDaW1Notes("lipid_reordered", "load_reordered"), loaded)
    expect_identical(omicDaW1Notes("lipid_main", "load"), loaded)
    expect_identical(omicDaW1Notes("metab_main", "load"), loaded)
    expect_identical(omicDaW1Notes("lipid_mismatch", "load_mismatch")[[3L]], paste("[error 10 -] Error loading session: Loaded lipid",
        "session assay_names do not match S4 lipid_data assays: session [LCMS_Pos, LCMS_Neg], S4 [LCMS_Pos, GCMS]"))
    expect_identical(omicDaW1Notes("metab_mismatch", "load_mismatch")[[3L]], paste("[error 10 -] Error loading session: Session",
        "assay_names do not match current_s4_object assays: session=LCMS_Pos, LCMS_Neg; object=LCMS_Pos, GCMS"))
    expect_identical(omicDaW1Notes("metab_reordered", "load_reordered")[[3L]], paste("[error 10 -] Error loading session:",
        "Session assay_names do not match current_s4_object assays: session=GCMS, LCMS_Pos; object=LCMS_Pos, GCMS"))
})

test_that("(c) the table filter: \"Combined\" selects every assay for metab and none for lipid and trx; a long table without friendly_name", {
    value <- function(case, step, output) omicDaW1Case(case)$captures[[step]]$values[[output]]
    expect_identical(value("lipid_main", "table_combined", "da_summary_stats"),
        "Total lipids: 0\nSignificant (Q < 0.050): 0 (0.0%)\n  Up-regulated: 0\n  Down-regulated: 0")
    expect_identical(value("metab_main", "table_combined", "da_summary_stats"),
        "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 6\n  Down-regulated: 2")
    expect_identical(value("trx_main", "table_combined", "da_summary_stats"),
        "Total genes: 0\nSignificant (FDR (BH) < 0.050): 0 (0.0%)\n  Up-regulated: 0\n  Down-regulated: 0")
    # Without friendly_name, lipid's filter fails (its `|` with an absent column) and metab's matches on comparison.
    expect_match(value("lipid_no_friendly", "table", "da_summary_stats"), "^<error> Can't subset rows")
    expect_identical(value("metab_no_friendly", "table", "da_summary_stats"),
        "Total metabolites: 12\nSignificant (Q < 0.050): 8 (66.7%)\n  Up-regulated: 4\n  Down-regulated: 4")
})

test_that("(c, D2) the value each DA server returns today: lipid its outputs, metab its registrations, trx and proteomics their reactive values", {
    # D2: the frame returns its reactive values for every omic. W3 and W4 change lipid's and metab's (a recorded divergence;
    # no caller reads the value) and retarget these two lines. W3: lipid's are the frame's reactive values (trx's fields).
    expect_identical(omicDaW1Case("lipid_refusals")$returned, paste0("reactivevalues: analysis_complete,assays_available,",
        "contrasts_available,current_col_clusters,current_heatmap_plot,current_row_clusters,current_state_name,da_results_list,",
        "formula_from_s4"))
    expect_identical(omicDaW1Case("metab_refusals")$returned,
        "list: loadFilteredSessionObserver,runAnalysisObserver,visualizationOutputs,saveHeatmapObserver,resultsOutputs")
    expect_identical(omicDaW1Case("trx_refusals")$returned, paste0("reactivevalues: analysis_complete,assays_available,",
        "contrasts_available,current_col_clusters,current_heatmap_plot,current_row_clusters,current_state_name,da_results_list,",
        "formula_from_s4"))
    expect_identical(omicDaW1Case("prot_refusals")$returned, paste0("reactivevalues: analysis_complete,contrasts_available,",
        "current_col_clusters,current_row_clusters,current_s4_object,da_results_list,formula_from_s4,volcano_plot_identity"))
})

# ------------------------------------------------------------------------------------------------
# W3 (D3, F106-05): the shared frame (R/mod_omic_da.R, R/mod_omic_da_panels.R) for lipid and trx, written before
# it (kit-106/w3/test-first.txt). The interactive volcano is unticked in the UI, `<req>` until ticked, then the
# binding's `<prefix>_volcano_page` (in process and managed); the tab's reactive values hold no state, widget or
# page bytes. F106-05: an owner whose head is its initial state is refused as no data.
# ------------------------------------------------------------------------------------------------

#' The paths in a tab's reactive values to a state (an omic's assay data object, or a state owner), a widget, a tag or
#' raw bytes. The heatmap's last drawing (`current_heatmap_plot`, the plot Save Heatmap writes; W1 (b) pins it) is left
#' out; so are the model fits the in-process kernel's list carries (trx's edgeR `qlf`), which are not states.
omicDaW3Heavy <- function(values) {
    found <- character()
    states <- c(lipidomicsDeclarationBinding$class, metabolomicsDeclarationBinding$class, transcriptomicsDeclarationBinding$class)
    walk <- function(x, path) {
        state <- isS4(x) && any(vapply(states, function(class) methods::is(x, class), logical(1)))
        if (state || inherits(x, c("DIAStateAdapter", "htmlwidget", "shiny.tag", "shiny.tag.list")) || is.raw(x)) {
            found <<- c(found, path)
        } else if (is.list(x) && !is.data.frame(x)) {
            for (i in seq_along(x)) walk(x[[i]], paste0(path, "$", if (is.null(names(x))) i else names(x)[[i]]))
        }
    }
    if (!shiny::is.reactivevalues(values)) return("<not reactive values>")
    fields <- shiny::isolate(shiny::reactiveValuesToList(values))
    for (name in setdiff(names(fields), "current_heatmap_plot")) walk(fields[[name]], name)
    found
}

#' Drive a DA tab to its completion and its interactive volcano: the UI's defaults (as the browser sends them),
#' `run(put)`, the selection, the volcano read unticked, then ticked (a managed page is read again after the flush that
#' delivers it). The page actions (`computeWorkflowOmicsVolcanoPageAction`, in process or in the worker) are counted.
omicDaW3Volcano <- function(setup, run, select) {
    log <- new.env(parent = emptyenv())
    log$notes <- character()
    log$updates <- character()
    log$pending <- list()
    omicDaW1InstallDoubles(log)
    seen <- new.env(parent = emptyenv())
    seen$pages <- character()
    page_action <- get("computeWorkflowOmicsVolcanoPageAction", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(computeWorkflowOmicsVolcanoPageAction = function(manager, request, variant) {
        seen$pages <- c(seen$pages, request$operation)
        page_action(manager, request, variant)
    }, .package = "MultiScholaR")
    out <- new.env(parent = emptyenv())
    mount <- setup$server(new.env(parent = emptyenv()))
    prefix <- paste0(setup$ns, "-")
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(mount, {
        flush_updates <- function() {
            for (k in seq_len(20L)) {
                if (!length(log$pending)) return(invisible())
                queued <- log$pending
                log$pending <- list()
                sent <- stats::setNames(lapply(queued, `[[`, "value"), vapply(queued, `[[`, "", "id"))
                do.call(session$setInputs, sent[!duplicated(names(sent), fromLast = TRUE)])
            }
        }
        put <- function(inputs) {
            names(inputs) <- ifelse(startsWith(names(inputs), "norm-"), names(inputs), paste0(prefix, names(inputs)))
            do.call(session$setInputs, inputs)
            session$flushReact()
            flush_updates()
        }
        glimma <- function() tryCatch(paste(as.character(output[[paste0(prefix, "volcano_glimma")]]$html), collapse = ""),
            error = function(e) if (inherits(e, c("shiny.silent.error", "shiny.output.cancel"))) "<req>" else paste("<error>", conditionMessage(e)))
        ui_defaults <- omicsUiInputDefaults(setup$ui("x"))
        names(ui_defaults) <- sub("^x-", "", names(ui_defaults))
        put(ui_defaults)
        run(put)
        put(select)
        out$unticked <- glimma()
        out$pages_unticked <- seen$pages
        put(list(volcano_interactive = TRUE))
        for (k in seq_len(3L)) {
            out$ticked <- glimma()
            if (grepl("<iframe", out$ticked, fixed = TRUE)) break
            session$flushReact()
        }
        out$pages <- seen$pages
        out$heavy <- omicDaW3Heavy(session$userData$da)
        out$notes <- log$notes
    }))))
    as.list(out)
}

omicDaW3FeatureRun <- function(put) {
    put(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
    put(list(run_da_analysis = 1L))
}

omicDaW3TrxRun <- function(put) {
    put(list(`norm-ruv_mode` = "automatic"))
    put(list(`norm-run_normalization` = 1L))
    put(list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
    put(list(run_da_analysis = 1L))
}

test_that("(W3, D3) lipid's and trx's DA UIs have the interactive volcano unticked", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    for (ui in list(lipid = mod_lipid_da_ui("de"), trx = mod_lipid_da_ui("de", profile = trxOmicProfile()))) {
        html <- as.character(htmltools::renderTags(ui)$html)
        expect_match(html, "<input id=\"de-volcano_interactive\" type=\"checkbox\"", fixed = TRUE)
        expect_false(grepl("<input id=\"de-volcano_interactive\"[^>]*checked", html))
        expect_false(omicsUiInputDefaults(ui)[["de-volcano_interactive"]])
    }
})

test_that("(W3) the shared frame's trx UI is W1 (a)'s trx UI", {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
    expect_true(exists("mod_omic_da_ui", envir = asNamespace("MultiScholaR"), inherits = FALSE))
    if (!exists("mod_omic_da_ui", envir = asNamespace("MultiScholaR"), inherits = FALSE)) return(invisible())
    parts <- omicDaW1UiParts(mod_omic_da_ui("de", omicGuiConfig("transcriptomics")$da))
    for (part in names(omicDaW1FrozenUi$trx)) expect_identical(parts[[part]], omicDaW1FrozenUi$trx[[part]], info = part)
})

test_that("(W3, D3) lipid in process: the volcano is <req> until ticked, then the lipid_volcano_page iframe; the tab holds no state, widget or page", {
    testthat::skip_if_not_installed("Glimma")
    setup <- omicDaW1FeatureSetup("lipid", session_file = FALSE)
    out <- omicDaW3Volcano(setup, omicDaW3FeatureRun, list(volcano_contrast = "KO_vs_WT", volcano_assay = "LCMS_Pos"))
    expect_identical(out$unticked, "<req>")
    expect_identical(out$pages_unticked, character())
    expect_match(out$ticked, "<iframe", fixed = TRUE)
    expect_identical(out$pages, "lipid_volcano_page")
    expect_identical(out$heavy, character())
})

test_that("(W3, D3) lipid managed: the volcano is <req> until ticked, then the worker's lipid_volcano_page iframe; the tab holds no state, widget or page", {
    testthat::skip_if_not_installed("Glimma")
    setup <- omicDaW1FeatureSetup("lipid", session_file = FALSE, managed = TRUE)
    out <- omicDaW3Volcano(setup, omicDaW3FeatureRun, list(volcano_contrast = "KO_vs_WT", volcano_assay = "LCMS_Pos"))
    expect_identical(out$unticked, "<req>")
    expect_identical(out$pages_unticked, character())
    expect_match(out$ticked, "<iframe", fixed = TRUE)
    expect_identical(out$pages, "lipid_volcano_page")
    expect_identical(out$heavy, character())
})

test_that("(W3, D3) trx in process: the volcano is <req> until ticked, then the trx_volcano_page iframe of the contrast's entry; the tab holds no state, widget or page", {
    testthat::skip_if_not_installed("Glimma")
    setup <- omicDaW1TrxSetup()
    out <- omicDaW3Volcano(setup, omicDaW3TrxRun, list(volcano_assay = "counts"))
    expect_identical(out$unticked, "<req>")
    expect_identical(out$pages_unticked, character())
    expect_match(out$ticked, "<iframe", fixed = TRUE)
    expect_match(out$ticked, "Interactive volcano Sera_vs_RPMI", fixed = TRUE)
    expect_identical(out$pages, "trx_volcano_page")
    expect_identical(out$heavy, character())
})

test_that("(W3, F106-05) lipid: Run DA on an owner whose head is its initial state is refused as no data", {
    setup <- omicDaW1FeatureSetup("lipid", state = FALSE, session_file = FALSE)
    omicsTestOwner(workflow = setup$bus(), omic_type = "lipidomics", experiment_paths = list(source_dir = file.path(setup$roots, "source")))
    run <- omicDaW1Drive(setup, list(run_no_data = function(ctx) {
        ctx$set(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
        ctx$set(run_da_analysis = 1L)
    }))$captures$run_no_data
    expect_identical(run$errors, character())
    expect_identical(run$notes, "[error 5 -] No lipidomics data loaded. Please load a filtered session first.")
    expect_identical(run$tab_status[["differential_analysis"]], "pending")
})
