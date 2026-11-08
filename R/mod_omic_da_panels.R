# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-106 W3 (D2, D5): the shared DA tab's panel library, from lipid's DA files at 2269a45f (mod_lipid_da_ui.R and
# mod_lipid_da_{display,results}_helpers.R, deleted by W3). The UI is layout data: a tree of nodes, each a Shiny tag
# function of `omicDaTagFunctions()` with its arguments and children, rendered by `omicDaRenderNode`. Lipid, metab and
# trx share one builder of those data from their `da` group (`omicDaFeatureLayout`); proteomics' are its own (W5). The
# pure views of the omics that hold rows follow: the status, contrasts, summary and cluster texts, the table filter, the
# table widget, the static volcano, the heatmap and the download. Recorded divergence from conventions/R-shiny.md: the
# wellPanel/column layout of the start commit is kept (the UI's HTML is pinned by test-omic-da-module.R (a)), under the
# existing-code scope rule of conventions/R.md.

# ------------------------------------------------------------------------------------------------
# The layout data and its renderer.
# ------------------------------------------------------------------------------------------------

#' One layout node: a tag function's name (`omicDaTagFunctions()`), its arguments and its children.
#' @param tag The tag function's name.
#' @param ... The function's arguments, in its order (positional or named).
#' @param id An input or output id, namespaced and passed first.
#' @param ns_id An element id, namespaced and passed first as `id`.
#' @param children Child nodes (NULL entries are dropped).
#' @param testid A `data-testid` set on the rendered tag.
#' @param resizable TRUE to wrap the tag with `shinyjqui::jqui_resizable`.
#' @return A plain list.
#' @noRd
omicDaNode <- function(tag, ..., id = NULL, ns_id = NULL, children = list(), testid = NULL, resizable = FALSE) {
    list(tag = tag, args = list(...), id = id, ns_id = ns_id, children = children, testid = testid, resizable = resizable)
}

#' The tag functions a layout may name.
#' @noRd
omicDaTagFunctions <- function() {
    list(fluidRow = shiny::fluidRow, column = shiny::column, wellPanel = shiny::wellPanel, div = shiny::div,
        h4 = shiny::h4, h5 = shiny::h5, hr = shiny::hr, br = shiny::br, helpText = shiny::helpText,
        textAreaInput = shiny::textAreaInput, numericInput = shiny::numericInput, selectInput = shiny::selectInput,
        checkboxInput = shiny::checkboxInput, actionButton = shiny::actionButton, downloadButton = shiny::downloadButton,
        verbatimTextOutput = shiny::verbatimTextOutput, uiOutput = shiny::uiOutput, htmlOutput = shiny::htmlOutput,
        plotOutput = shiny::plotOutput, DTOutput = DT::DTOutput, tabsetPanel = shiny::tabsetPanel,
        tabPanel = shiny::tabPanel, conditionalPanel = shiny::conditionalPanel)
}

#' Render one layout node under the namespace `ns`: an `icon` argument is a Shiny icon name; a conditional panel is
#' given `ns`.
#' @return A tag.
#' @noRd
omicDaRenderNode <- function(ns, node) {
    # W5: a plugin node is its group plugin's UI, `function(ns)`, looked up by name (proteomics' archive chooser).
    if (!is.null(node$plugin)) return(omicGuiFunction(node$plugin)(ns))
    args <-c(if (!is.null(node$ns_id)) list(id = ns(node$ns_id)), node$args)
    if (!is.null(args$icon)) args$icon <- shiny::icon(args$icon)
    if (identical(node$tag, "conditionalPanel")) args$ns <- ns
    children <- lapply(Filter(Negate(is.null), node$children), \(child) omicDaRenderNode(ns, child))
    tag <- do.call(omicDaTagFunctions()[[node$tag]], c(if (!is.null(node$id)) list(ns(node$id)), args, children))
    if (isTRUE(node$resizable)) tag <- shinyjqui::jqui_resizable(tag)
    if (!is.null(node$testid)) tag <- testid(tag, node$testid)
    tag
}

#' The layout of an omic's DA tab: the feature omics' (lipid, metab, trx) from their `da` group; an omic whose group
#' has plugins (proteomics) its own.
#' @param da The omic's `da` group.
#' @return The root node.
#' @noRd
omicDaLayout <- function(da) {
    if (!is.null(da$plugins)) return(omicDaPluginLayout(da))
    omicDaFeatureLayout(da)
}

#' The DA tab of the omics that hold rows, as layout data (mod_lipid_da_ui.R:9-444 at 2269a45f): the settings column
#' (the formula input or the model shown, the threshold and fold-change inputs, the contrasts, the session load when the
#' group has its test id, the run and download buttons and the status) and the results column (the save note with the
#' heatmap's, then the volcano, heatmap and table tabs the group names). The interactive volcano's box is unticked
#' (D3: built on demand).
#' @param da The omic's `da` group.
#' @return The root node.
#' @noRd
omicDaFeatureLayout <- function(da) {
    node <- omicDaNode
    model <- da$model
    model_nodes <- if (isTRUE(model$editable)) {
        list(node("textAreaInput", model$label, value = model$default, height = "60px", id = "formula_string"),
            do.call(node, c(list("helpText"), as.list(model$help))))
    } else {
        list(node("h5", model$heading), node("verbatimTextOutput", id = model$output),
            do.call(node, c(list("helpText"), as.list(model$help))))
    }
    numeric <- function(id, spec) {
        node("numericInput", spec$label, value = spec$value, min = spec$min, max = spec$max, step = spec$step, id = id)
    }
    load <- if ("load_session" %in% names(da$testids)) {
        list(node("actionButton", "Load Filtered Session", class = "btn-info", width = "100%", icon = "upload",
            id = "load_filtered_session", testid = da$testids[["load_session"]]),
            node("helpText", "Load previously exported filtered data for DA analysis"), node("br"), node("br"))
    }
    sidebar <- c(list(node("h4", "Differential Abundance Analysis")), model_nodes, list(node("hr"),
        node("h5", "Analysis Parameters"), numeric("da_q_val_thresh", da$threshold), numeric("treat_lfc_cutoff", da$lfc),
        node("hr"),
        node("h5", "Available Contrasts"), node("verbatimTextOutput", id = "contrasts_display"),
        node("helpText", "Contrasts defined in design matrix"), node("hr")), load, list(
        node("actionButton", "Run DA Analysis", class = "btn-primary", width = "100%", icon = "play", id = "run_da_analysis",
            testid = da$testids[["run"]]),
        node("br"), node("br"),
        node("downloadButton", "Download All Results", class = "btn-success", style = "width: 100%;",
            id = "download_da_results", testid = da$testids[["download"]]),
        node("hr"),
        node("h5", "Analysis Status"), node("verbatimTextOutput", id = "da_status")))
    tabs <- list(volcano = omicDaFeatureVolcanoTab(da), heatmap = if ("heatmap" %in% names(da$tabs)) omicDaFeatureHeatmapTab(da),
        table = omicDaFeatureTableTab(da))
    tabset <- node("tabsetPanel", ns_id = "da_results_tabs", testid = da$testids[["results_tabs"]],
        children = unname(tabs[names(da$tabs)]))
    note <- if (isTRUE(da$views$heatmap$save_note)) node("uiOutput", id = "heatmap_manual_save_warning")
    node("fluidRow", children = list(node("column", 3, children = list(node("wellPanel", children = sidebar))),
        node("column", 9, children = list(note, tabset))))
}

#' A select input's node.
#' @noRd
omicDaSelect <- function(id, label, choices, selected = NULL, width = NULL) {
    omicDaNode("selectInput", label, choices = choices, selected = selected, width = width, id = id)
}

#' A fluid row of columns, each `list(width, node)`.
#' @noRd
omicDaRow <- function(...) {
    omicDaNode("fluidRow", children = lapply(list(...), \(cell) omicDaNode("column", cell[[1L]], children = cell[-1L])))
}

#' The volcano tab (mod_lipid_da_ui.R:119-171 at 2269a45f).
#' @noRd
omicDaFeatureVolcanoTab <- function(da) {
    node <- omicDaNode
    static <- da$views$static_volcano$output
    node("tabPanel", da$tabs[["volcano"]], children = list(node("br"),
        omicDaRow(list(4, omicDaSelect("volcano_contrast", "Select Contrast:", NULL, width = "100%")),
            list(4, omicDaSelect("volcano_assay", "Select Assay:", c("Combined"), width = "100%")),
            list(4, node("checkboxInput", "Interactive Plot (Glimma)", value = FALSE, id = "volcano_interactive"))),
        node("conditionalPanel", condition = "input.volcano_interactive == true", children = list(
            node("div", style = "height: 600px; background-color: white; border-radius: 8px; padding: 10px;",
                ns_id = "volcano_glimma_container", children = list(node("htmlOutput", id = "volcano_glimma"))))),
        node("conditionalPanel", condition = "input.volcano_interactive == false", children = list(
            node("plotOutput", height = "600px", id = static, resizable = TRUE)))))
}

#' The heatmap tab with the group's nouns and tree-cut defaults (mod_lipid_da_ui.R:176-367 at 2269a45f).
#' @noRd
omicDaFeatureHeatmapTab <- function(da) {
    node <- omicDaNode
    nouns <- da$nouns
    cut <- da$views$heatmap$tree_cut
    when <- function(method, input) {
        node("conditionalPanel", condition = sprintf("input.heatmap_tree_cut_method == '%s'", method), children = list(input))
    }
    options <- node("wellPanel", children = list(node("h5", "Clustering Options"),
        omicDaRow(list(4, omicDaSelect("heatmap_cluster_method", "Clustering Method:", list("Ward (minimum variance)" = "ward.D2",
                "Ward (original)" = "ward.D", "Complete linkage" = "complete", "Single linkage" = "single",
                "Average linkage" = "average", "McQuitty (WPGMA)" = "mcquitty"), selected = "ward.D2")),
            list(4, omicDaSelect("heatmap_distance_method", "Distance Metric:", list("Euclidean" = "euclidean",
                "Manhattan" = "manhattan", "Pearson correlation" = "pearson", "Spearman correlation" = "spearman",
                "Maximum" = "maximum"), selected = "euclidean")),
            list(4, omicDaSelect("heatmap_color_scheme", "Color Scheme:", list("Red-Blue" = "RdBu", "Red-Yellow-Blue" = "RdYlBu",
                "Blue-White-Red" = "coolwarm", "Viridis" = "viridis", "Plasma" = "plasma", "Inferno" = "inferno"),
                selected = "RdBu"))),
        omicDaRow(list(4, omicDaSelect("heatmap_assay", "Select Assay:", c("Combined"), width = "100%")),
            list(4, node("checkboxInput", nouns[["labels"]], value = FALSE, id = "heatmap_show_labels"))),
        omicDaRow(list(3, omicDaSelect("heatmap_tree_cut_method", "Tree Cutting:", c("None" = "none",
                "K Clusters" = "k_clusters", "Height Cutoff" = "height_cutoff", "Dynamic" = "dynamic"), selected = cut$default)),
            list(3,
                when("k_clusters", node("numericInput", "Number of Clusters:", value = cut$n_clusters, min = 2, max = 20, step = 1,
                    id = "heatmap_n_clusters")),
                when("height_cutoff", node("numericInput", "Cut Height:", value = cut$cut_height, min = 1, max = 100, step = 1,
                    id = "heatmap_cut_height")),
                when("dynamic", node("numericInput", "Min Cluster Size:", value = cut$min_cluster_size, min = 2, max = 20,
                    step = 1, id = "heatmap_min_cluster_size"))))))
    node("tabPanel", da$tabs[["heatmap"]], children = list(node("br"),
        omicDaRow(list(3, omicDaSelect("heatmap_contrast", "Select Contrast:", NULL, width = "100%")),
            list(3, node("numericInput", nouns[["top_n"]], value = 50, min = 10, max = 500, step = 10, id = "heatmap_top_n")),
            list(3, omicDaSelect("heatmap_clustering", "Apply Clustering:", list("Both rows & columns" = "both",
                "Rows only" = "row", "Columns only" = "column", "None" = "none"), selected = "both")),
            list(3, omicDaSelect("heatmap_scaling", "Data Scaling:", stats::setNames(list("row", "column", "both", "none"),
                c(nouns[["scaling"]], "Column (sample) scaling", "Both", "None")), selected = "row"))),
        node("conditionalPanel", condition = "input.heatmap_clustering != 'none'", children = list(options)),
        node("div", style = "height: 650px;", ns_id = "heatmap_container", children = list(
            node("plotOutput", height = "600px", id = "heatmap_plot", resizable = TRUE))),
        node("div", style = "margin-top: 10px; margin-bottom: 20px;", children = list(
            node("actionButton", "Save Heatmap", icon = "download", class = "btn-primary", id = "save_heatmap"))),
        node("conditionalPanel", condition = "input.heatmap_tree_cut_method != 'none'", children = list(
            node("wellPanel", children = list(node("h4", "Cluster Information"),
                node("verbatimTextOutput", id = "cluster_summary")))))))
}

#' The table tab with the group's defaults (mod_lipid_da_ui.R:372-436 at 2269a45f).
#' @noRd
omicDaFeatureTableTab <- function(da) {
    node <- omicDaNode
    table <- da$views$table
    node("tabPanel", da$tabs[["table"]], children = list(node("br"),
        omicDaRow(list(3, omicDaSelect("table_contrast", "Select Contrast:", NULL, width = "100%")),
            list(3, omicDaSelect("table_assay", "Select Assay:", c("All"), width = "100%")),
            list(3, omicDaSelect("table_significance", "Show:", list("All Results" = "all", "Significant Only" = "significant",
                "Up-regulated" = "up", "Down-regulated" = "down"), selected = table$significance)),
            list(3, node("numericInput", "Max Rows:", value = table$max_rows, min = 100, max = 10000, step = 100,
                id = "table_max_rows"))),
        omicDaRow(list(12, node("wellPanel", children = list(node("h5", "Summary Statistics"),
            node("verbatimTextOutput", id = "da_summary_stats"))))),
        node("DTOutput", id = "da_results_table")))
}

# ------------------------------------------------------------------------------------------------
# The views of the omics that hold rows.
# ------------------------------------------------------------------------------------------------

#' The status text (mod_lipid_da_display_helpers.R:1-33 at 2269a45f): each assay's counts once the analysis has run,
#' else the group's ready text when a state is loaded and its waiting text otherwise.
#' @param texts The group's `texts`.
#' @noRd
omicDaStatusText <- function(analysisComplete, daResultsList, loaded, texts) {
    if (isTRUE(analysisComplete)) {
        counts <- daResultsList$significant_counts
        if (is.null(counts)) return("Analysis complete.")
        lines <- vapply(names(counts), \(assay) sprintf("%s:\n  Up: %d | Down: %d | NS: %d", assay, counts[[assay]]$up,
            counts[[assay]]$down, counts[[assay]]$ns), character(1))
        return(paste(c("Analysis Complete", "", lines), collapse = "\n"))
    }
    if (!is.null(loaded)) return(texts[["status_ready"]])
    texts[["status_waiting"]]
}

#' The contrasts text (mod_lipid_da_display_helpers.R:35-52 at 2269a45f).
#' @noRd
omicDaContrastsText <- function(contrastsTbl, texts) {
    if (is.null(contrastsTbl) || nrow(contrastsTbl) == 0) return(texts[["contrasts_none"]])
    if ("friendly_names" %in% colnames(contrastsTbl)) return(paste(contrastsTbl$friendly_names, collapse = "\n"))
    if ("contrasts" %in% colnames(contrastsTbl)) return(paste(contrastsTbl$contrasts, collapse = "\n"))
    paste(utils::capture.output(print(contrastsTbl)), collapse = "\n")
}

#' The rows of a contrast and an assay under the group's table filter (D7): `all_only`, lipid's (a contrast matched by
#' `comparison` or `friendly_name`; every assay for "All" only; mod_lipid_da_results_helpers.R:20, :117 at 2269a45f), or
#' `all_or_combined`, metab's (a missing `comparison` or `friendly_name` column matches nothing; every assay for "All"
#' and "Combined"; rows kept as a data frame; mod_metab_da_display_helpers.R:13-29 at 2269a45f).
#' @noRd
omicDaFilterRows <- function(results, selectedContrast, selectedAssay, tableFilter = "all_only") {
    if (identical(tableFilter, "all_or_combined")) {
        if (!is.null(selectedContrast)) {
            matches <- function(column) if (column %in% colnames(results)) results[[column]] == selectedContrast else rep(FALSE, nrow(results))
            results <- results[matches("comparison") | matches("friendly_name"), , drop = FALSE]
        }
        if (!is.null(selectedAssay) && !selectedAssay %in% c("All", "Combined")) {
            results <- results[results$assay == selectedAssay, , drop = FALSE]
        }
        return(results)
    }
    if (!is.null(selectedContrast)) {
        results <- results[results$comparison == selectedContrast | results$friendly_name == selectedContrast, ]
    }
    if (!is.null(selectedAssay) && selectedAssay != "All") results <- results[results$assay == selectedAssay, ]
    results
}

#' The summary statistics of a contrast and an assay (mod_lipid_da_results_helpers.R:1-41 at 2269a45f), in the group's
#' noun and significance label.
#' @param summary The group's `summary` (`noun`, `significance_label`).
#' @noRd
omicDaSummaryText <- function(daResults, selectedContrast, selectedAssay, daQValThresh, summary, tableFilter = "all_only") {
    if (is.null(daResults) || nrow(daResults) == 0) return("No results available.")
    results <- omicDaFilterRows(daResults, selectedContrast, selectedAssay, tableFilter)
    total <- nrow(results)
    sig <- sum(results$significant != "NS", na.rm = TRUE)
    paste(c(sprintf("Total %s: %d", summary$noun, total),
        sprintf("Significant (%s < %.3f): %d (%.1f%%)", summary$significance_label, daQValThresh, sig, 100 * sig / max(total, 1)),
        sprintf("  Up-regulated: %d", sum(results$significant == "Up", na.rm = TRUE)),
        sprintf("  Down-regulated: %d", sum(results$significant == "Down", na.rm = TRUE))), collapse = "\n")
}

#' The results table of a contrast and an assay (mod_lipid_da_results_helpers.R:87-171 at 2269a45f): filtered by
#' significance and fold change, cut to the row limit, the group's columns shown (the significance column under its
#' label when the group has one).
#' @param columns The group's `columns` (`id`, `label`, `fold_change`, `pvalue`, `significance`, `significance_label`).
#' @param datatableFactory,formatRoundFactory,formatStyleFactory,styleEqualFactory The widget's DT functions.
#' @return The widget, or NULL when no row is left.
#' @noRd
omicDaTableWidget <- function(daResults, selectedContrast, selectedAssay, tableSignificance, daQValThresh, lfcThreshold, tableMaxRows,
    columns, tableFilter = "all_only", datatableFactory = DT::datatable, formatRoundFactory = DT::formatRound,
    formatStyleFactory = DT::formatStyle, styleEqualFactory = DT::styleEqual) {
    if (is.null(daResults) || nrow(daResults) == 0) return(NULL)
    significance <- columns$significance
    fold <- columns$fold_change
    results <- omicDaFilterRows(daResults, selectedContrast, selectedAssay, tableFilter)
    if (identical(tableSignificance, "significant")) {
        results <- results[results[[significance]] < daQValThresh, ]
    } else if (identical(tableSignificance, "up")) {
        results <- results[results[[significance]] < daQValThresh & results[[fold]] > lfcThreshold, ]
    } else if (identical(tableSignificance, "down")) {
        results <- results[results[[significance]] < daQValThresh & results[[fold]] < -lfcThreshold, ]
    }
    if (nrow(results) == 0) return(NULL)
    if (nrow(results) > tableMaxRows) results <- results[seq_len(tableMaxRows), ]
    shown <- intersect(c(columns$id, columns$label, "assay", fold, columns$pvalue, significance, "significant"), colnames(results))
    results <- results[, shown, drop = FALSE]
    args <- list(results, options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
        extensions = "Buttons", rownames = FALSE, filter = "top")
    if (!is.null(columns$significance_label)) args$colnames <- ifelse(shown == significance, columns$significance_label, shown)
    do.call(datatableFactory, args) |>
        formatRoundFactory(columns = fold, digits = 3) |>
        formatRoundFactory(columns = c(columns$pvalue, significance), digits = 6) |>
        formatStyleFactory("significant", backgroundColor = styleEqualFactory(c("Up", "Down", "NS"), c("#ffcccc", "#cce5ff", "white")))
}

#' The static volcano of a contrast and an assay (mod_lipid_da_display_helpers.R:121-140 at 2269a45f).
#' @param view The group's `views$static_volcano` (`plot`, `show_labels`, `n_labels`).
#' @noRd
omicDaStaticVolcano <- function(daResultsList, selectedContrast, selectedAssay, daQValThresh, lfcThreshold, view) {
    shiny::req(daResultsList, selectedContrast)
    view$plot(da_results_list = daResultsList, selected_contrast = selectedContrast, selected_assay = selectedAssay,
        da_q_val_thresh = daQValThresh, lfc_threshold = lfcThreshold, show_labels = view$show_labels, n_labels = view$n_labels)
}

#' The interactive volcano of a contrast and an assay (mod_lipid_da_display_helpers.R:54-101 at 2269a45f): the banner for
#' "Combined" (no page is asked for), the page, or a message when the page is absent or failed.
#' @param widgetFactory The page factory (`function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh)`).
#' @noRd
omicDaVolcanoGlimmaUi <- function(daResultsList, selectedContrast, selectedAssay, daQValThresh, widgetFactory) {
    shiny::req(daResultsList, selectedContrast)
    if (is.null(selectedAssay) || identical(selectedAssay, "Combined")) {
        return(shiny::div(class = "alert alert-info", style = "margin: 20px;", shiny::icon("info-circle"),
            shiny::tags$strong(" Combined View: "), "Interactive Glimma plots require a single assay selection. ",
            "Please select a specific assay (e.g., LCMS_Pos) or uncheck 'Interactive Plot' to use the static volcano plot for combined view."))
    }
    tryCatch({
        widget <- widgetFactory(da_results_list = daResultsList, selected_contrast = selectedContrast,
            selected_assay = selectedAssay, da_q_val_thresh = daQValThresh)
        if (is.null(widget)) {
            return(shiny::div(class = "alert alert-warning", "Could not generate Glimma plot. Try the static plot option."))
        }
        widget
    }, error = function(e) {
        logger::log_error(sprintf("Glimma error: %s", conditionMessage(e)))
        shiny::div(class = "alert alert-danger", sprintf("Error generating plot: %s", conditionMessage(e)))
    })
}

#' The heatmap save note (mod_lipid_da_display_helpers.R:103-119 at 2269a45f), once the analysis has run.
#' @noRd
omicDaHeatmapSaveNote <- function(analysisComplete) {
    if (!isTRUE(analysisComplete)) return(NULL)
    shiny::div(class = "alert alert-info", style = "margin-bottom: 15px; padding: 10px; border-left: 5px solid #17a2b8;",
        shiny::icon("info-circle"), shiny::tags$b(" Note: "),
        "Heatmaps are NOT saved automatically during analysis. Please navigate to the ", shiny::tags$b("Heatmap"),
        " tab and use the ", shiny::tags$b("Save Heatmap"), " button to save your specific clustering configuration and image.")
}

#' The heatmap of a contrast and an assay (mod_lipid_da_display_helpers.R:142-200 at 2269a45f): the group's kernel, its
#' labels argument by name, and the plot with its row and column clusters. Each input is passed unevaluated, so the
#' kernel reads it when it uses it, as it did.
#' @param heatmap The group's `views$heatmap` (`kernel`, `labels_argument`).
#' @return `list(plot, rowClusters, colClusters, currentHeatmapPlot)`, or NULL when the kernel draws nothing.
#' @noRd
omicDaHeatmapState <- function(daResultsList, selectedContrast, selectedAssay, topN, heatmapClustering, scaleData, clusteringMethod,
    distanceMethod, colorScheme, showLabels, daQValThresh, treeCutMethod, nClusters, cutHeight, minClusterSize, heatmap) {
    shiny::req(daResultsList, selectedContrast)
    args <- alist(da_results_list = daResultsList, selected_contrast = selectedContrast, selected_assay = selectedAssay,
        top_n = topN, clustering_method = clusteringMethod, distance_method = distanceMethod,
        cluster_rows = heatmapClustering %in% c("both", "row"), cluster_cols = heatmapClustering %in% c("both", "column"),
        scale_data = scaleData, color_scheme = colorScheme, labels = showLabels, da_q_val_thresh = daQValThresh,
        tree_cut_method = treeCutMethod, n_clusters = nClusters, cut_height = cutHeight, min_cluster_size = minClusterSize)
    names(args)[names(args) == "labels"] <- heatmap$labels_argument
    hm <- eval(as.call(c(list(heatmap$kernel), args)))
    if (is.null(hm)) return(NULL)
    if (is.list(hm) && "plot" %in% names(hm)) {
        return(list(plot = hm$plot, rowClusters = hm$row_clusters, colClusters = hm$col_clusters, currentHeatmapPlot = hm$plot))
    }
    list(plot = hm, rowClusters = NULL, colClusters = NULL, currentHeatmapPlot = hm)
}

#' The tree cut's cluster summary (mod_lipid_da_display_helpers.R:214-256 at 2269a45f), in the group's cluster noun.
#' @noRd
omicDaClusterSummaryText <- function(rowClusters, noun) {
    if (is.null(rowClusters)) return("No clusters defined. Enable clustering and tree cutting on the heatmap.")
    clusters <- sort(unique(rowClusters))
    lines <- c(sprintf("Total Clusters: %d", length(clusters)), "----------------------------------")
    for (cluster in clusters) {
        members <- names(rowClusters)[rowClusters == cluster]
        shown <- paste(utils::head(members, 20), collapse = ", ")
        if (length(members) > 20) shown <- paste0(shown, sprintf(", ... and %d more", length(members) - 20))
        lines <- c(lines, "", sprintf("Cluster %d (%d %s):", cluster, length(members), noun), shown)
    }
    paste(lines, collapse = "\n")
}

#' The download handler (mod_lipid_da_results_helpers.R:173-210 at 2269a45f): the group's file name, and at each click
#' the results read then: a managed tab's writer (the store's tables), else its long table as CSV.
#' @param readResults `function()` giving the tab's results list now (KIT-082: a value here would be forced once and
#'   serve the first run's results for ever).
#' @param download The group's `download` (`stem`, `extension`).
#' @noRd
omicDaDownloadHandler <- function(readResults, longTable, download, downloadHandlerFactory = shiny::downloadHandler, reqFn = shiny::req) {
    downloadHandlerFactory(filename = function() paste0(download$stem, Sys.Date(), ".", download$extension),
        content = function(file) {
            results <- omicDaDownloadData(readResults(), longTable)
            reqFn(results)
            if (is.function(results)) return(results(file))
            utils::write.csv(results, file, row.names = FALSE)
        })
}

#' What the download writes: a managed tab's writer (KIT-082: it holds no per-sample columns), else its long table.
#' @noRd
omicDaDownloadData <- function(daResultsList, longTable) {
    daResultsList$write_download %||% daResultsList[[longTable]]
}

# ------------------------------------------------------------------------------------------------
# Proteomics' layout data (KIT-106 W5, D5).
# ------------------------------------------------------------------------------------------------

#' The DA tab of an omic whose group has plugins, as layout data: proteomics' (mod_prot_da_ui.R:13-493 at fe97aac1).
#' Its own layout and values: the outer well, the sidebar heading and spacing, the archive chooser (the group's
#' `archive` plugin), tab icons, column widths, the "Advanced Clustering Options" well with its choices and the
#' dendrogram switch, `h5` for "Cluster Information", and the save note's output that no server renders (P-25). The
#' interactive volcano's box is unticked (D3: built on demand).
#' @param da The omic's `da` group.
#' @return The root node.
#' @noRd
omicDaPluginLayout <- function(da) {
    node <- omicDaNode
    numeric <- function(id, spec) {
        node("numericInput", spec$label, value = spec$value, min = spec$min, max = spec$max, step = spec$step, id = id)
    }
    sidebar <- node("wellPanel", children = list(node("h4", "DA Analysis Settings"),
        node("textAreaInput", da$model$label, value = da$model$default, height = "60px", id = "formula_string"),
        node("helpText", da$model$help), node("hr"),
        node("h5", "Analysis Parameters"), numeric("da_q_val_thresh", da$threshold), numeric("treat_lfc_cutoff", da$lfc),
        node("hr"),
        node("h5", "Available Contrasts"), node("verbatimTextOutput", id = "contrasts_display"),
        node("helpText", "Contrasts defined in design matrix"), node("br"),
        list(plugin = da$plugins[["archive"]]),
        node("actionButton", "Load Filtered Session", class = "btn-info", width = "100%", icon = "upload",
            id = "load_filtered_session", testid = da$testids[["load_session"]]),
        node("helpText", "Load the latest session, a portable native DIA archive, or a historical RDS session."), node("br"),
        node("actionButton", "Run DA Analysis", class = "btn-primary", width = "100%", icon = "play", id = "run_da_analysis",
            testid = da$testids[["run"]]),
        node("br"), node("br"),
        node("downloadButton", "Download All Results", class = "btn-success", width = "100%", id = "download_da_results",
            testid = da$testids[["download"]]),
        node("br"), node("br"),
        node("h5", "Analysis Status"), node("verbatimTextOutput", id = "da_status")))
    tabset <- node("tabsetPanel", ns_id = "da_results_tabs", testid = da$testids[["results_tabs"]],
        children = list(omicDaPluginVolcanoTab(da), omicDaPluginHeatmapTab(da), omicDaPluginTableTab(da)))
    node("wellPanel", children = list(node("fluidRow", children = list(node("column", 3, children = list(sidebar)),
        node("column", 9, children = list(node("uiOutput", id = "heatmap_manual_save_warning"), tabset))))))
}

#' Proteomics' volcano tab (mod_prot_da_ui.R:136-204 at fe97aac1).
#' @noRd
omicDaPluginVolcanoTab <- function(da) {
    node <- omicDaNode
    static <- da$views$static_volcano
    node("tabPanel", da$tabs[["volcano"]], icon = "chart-area", children = list(node("br"),
        omicDaRow(list(6, omicDaSelect("volcano_contrast", "Select Contrast:", NULL, width = "100%")),
            list(6, node("checkboxInput", "Interactive Plot (Glimma)", value = FALSE, id = "volcano_interactive"))),
        node("hr"),
        node("conditionalPanel", condition = "input.volcano_interactive == true", children = list(
            node("div", style = "height: 600px; background-color: white; border-radius: 8px; padding: 10px;",
                ns_id = "volcano_glimma_container", children = list(node("htmlOutput", id = "volcano_glimma"))))),
        node("conditionalPanel", condition = "input.volcano_interactive == false", children = list(
            # KIT-101 F-02: the static renderer's label inputs, at the defaults it applies when absent.
            omicDaRow(list(6, node("checkboxInput", da$nouns[["label_top"]], value = static$show_labels, id = "volcano_show_labels")),
                list(6, node("numericInput", "Number of Labels:", value = static$n_labels, min = 0, max = 100, step = 1,
                    id = "volcano_label_top_n"))),
            node("plotOutput", height = "600px", id = static$output, resizable = TRUE)))))
}

#' Proteomics' heatmap tab with its tree-cut defaults (mod_prot_da_ui.R:206-432 at fe97aac1).
#' @noRd
omicDaPluginHeatmapTab <- function(da) {
    node <- omicDaNode
    nouns <- da$nouns
    cut <- da$views$heatmap$tree_cut
    when <- function(method, input) {
        node("conditionalPanel", condition = sprintf("input.heatmap_tree_cut_method == '%s'", method), children = list(input))
    }
    options <- node("wellPanel", children = list(node("h5", "Advanced Clustering Options"),
        omicDaRow(list(3, omicDaSelect("heatmap_cluster_method", "Clustering Method:", list("Ward (minimum variance)" = "ward.D2",
                "Ward (original)" = "ward.D", "Complete linkage" = "complete", "Single linkage" = "single",
                "Average linkage" = "average", "McQuitty (WPGMA)" = "mcquitty", "Median (WPGMC)" = "median",
                "Centroid (UPGMC)" = "centroid"), selected = "ward.D2")),
            list(3, omicDaSelect("heatmap_distance_method", "Distance Metric:", list("Euclidean" = "euclidean",
                "Manhattan" = "manhattan", "Pearson correlation" = "pearson", "Spearman correlation" = "spearman",
                "Maximum" = "maximum", "Canberra" = "canberra", "Binary" = "binary", "Minkowski" = "minkowski"),
                selected = "euclidean")),
            list(3, omicDaSelect("heatmap_tree_cut_method", "Tree Cutting:", list("Number of clusters" = "k_clusters",
                "Height cutoff" = "height_cutoff", "Dynamic tree cutting" = "dynamic", "No cutting" = "none"),
                selected = cut$default)),
            list(3,
                when("k_clusters", node("numericInput", "Number of Clusters:", value = cut$n_clusters, min = 2, max = 20, step = 1,
                    id = "heatmap_n_clusters")),
                when("height_cutoff", node("numericInput", "Cut Height:", value = cut$cut_height, min = 0.1, max = 2.0, step = 0.1,
                    id = "heatmap_cut_height")),
                when("dynamic", node("numericInput", "Min Cluster Size:", value = cut$min_cluster_size, min = 2, max = 20,
                    step = 1, id = "heatmap_min_cluster_size")))),
        omicDaRow(list(4, node("checkboxInput", "Show Dendrogram", value = TRUE, id = "heatmap_show_dendro")),
            list(4, node("checkboxInput", nouns[["labels"]], value = FALSE, id = "heatmap_show_labels")),
            list(4, omicDaSelect("heatmap_color_scheme", "Color Scheme:", list("Red-Blue" = "RdBu", "Red-Yellow-Blue" = "RdYlBu",
                "Blue-White-Red" = "coolwarm", "Viridis" = "viridis", "Plasma" = "plasma", "Inferno" = "inferno"),
                selected = "RdBu")))))
    node("tabPanel", da$tabs[["heatmap"]], icon = "th", children = list(node("br"),
        omicDaRow(list(3, omicDaSelect("heatmap_contrast", "Select Contrast:", NULL, width = "100%")),
            list(3, node("numericInput", nouns[["top_n"]], value = 50, min = 10, max = 500, step = 10, id = "heatmap_top_n")),
            list(3, omicDaSelect("heatmap_clustering", "Apply Clustering:", list("Both rows & columns" = "both",
                "Rows only" = "row", "Columns only" = "column", "None" = "none"), selected = "both")),
            list(3, omicDaSelect("heatmap_scaling", "Data Scaling:", stats::setNames(list("row", "column", "both", "none"),
                c(nouns[["scaling"]], "Column (sample) scaling", "Both", "None")), selected = "row"))),
        node("conditionalPanel", condition = "input.heatmap_clustering != 'none'", children = list(options)),
        node("hr"),
        node("div", style = "height: 650px;", ns_id = "heatmap_container", children = list(
            node("plotOutput", height = "600px", id = "heatmap_plot", resizable = TRUE))),
        node("div", style = "margin-top: 10px; margin-bottom: 20px;", children = list(
            node("actionButton", "Save Heatmap", icon = "download", class = "btn-primary", id = "save_heatmap"))),
        node("conditionalPanel", condition = "input.heatmap_tree_cut_method != 'none'", children = list(
            node("wellPanel", children = list(node("h5", "Cluster Information"),
                node("verbatimTextOutput", id = "cluster_summary")))))))
}

#' Proteomics' table tab with its defaults (mod_prot_da_ui.R:434-493 at fe97aac1).
#' @noRd
omicDaPluginTableTab <- function(da) {
    node <- omicDaNode
    table <- da$views$table
    node("tabPanel", da$tabs[["table"]], icon = "table", children = list(node("br"),
        omicDaRow(list(4, omicDaSelect("table_contrast", "Select Contrast:", NULL, width = "100%")),
            list(4, omicDaSelect("table_significance", "Show:", list("All Results" = "all", "Significant Only" = "significant",
                "Up-regulated" = "up", "Down-regulated" = "down"), selected = table$significance)),
            list(4, node("numericInput", "Max Rows:", value = table$max_rows, min = 100, max = 10000, step = 100,
                id = "table_max_rows"))),
        node("hr"),
        omicDaRow(list(12, node("wellPanel", children = list(node("h5", "Summary Statistics"),
            node("verbatimTextOutput", id = "da_summary_stats"))))),
        node("DTOutput", id = "da_results_table")))
}
