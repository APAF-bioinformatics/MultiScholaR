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

# KIT-105 W4: the normalisation tab's panel library (kit-101 design 4.7), from the lipid tab's
# helper files at f8d4c5fc (mod_lipid_norm_ui_helpers.R, mod_lipid_norm_support_helpers.R,
# mod_lipid_norm_runtime_helpers.R and mod_lipid_norm_observer_helpers.R; deleted by W4): the
# option sections of `config$norm$panels`, the QC tabset with its stage grid, the correlation-filter
# tab, and the renderers and registrars of the tab's outputs. The shared frame (mod_omic_norm.R)
# assembles them. Ids, texts and markup are the lipid tab's, byte for byte; KIT-110 added trx's
# TMM/logCPM and RUV grid sections and the reset's optional test id, trx's markup byte for byte.
# Recorded divergence from conventions/R-shiny.md: the wellPanel layout is kept (existing-code
# scope rule).

# ------------------------------------------------------------------------------------------------
# The option sections (column 3).
# ------------------------------------------------------------------------------------------------

#' The design-driven select of the aesthetics and the RUV grouping (`ns(id)`, its four default choices).
#' @noRd
omicNormDesignSelect <- function(ns, id, label) {
    shiny::selectInput(ns(id), label, choices = c("group", "factor1", "factor2", "batch"), selected = "group", width = "100%")
}

#' The option sections `norm$panels` names, in its order: each a list of tags in the well panel
#' (mod_lipid_norm_ui_helpers.R:6-261 at f8d4c5fc; KIT-110: trx's `tmm_logcpm` and `ruv_grid`, mod_trx_norm.R:26-47 at
#' 38d9d25e). Only the named sections are built. `norm` is the omic's `norm` group. KIT-120: a section shows the
#' omic's own text where its options name one (`help`, `run_label`, `label`), else the library's.
#' @noRd
omicNormOptionSections <- function(ns, norm) {
    panels <- norm$panels
    section <- function(name) switch(name,
        aesthetics = list(shiny::hr(), shiny::h5("Plot Aesthetics"), omicNormDesignSelect(ns, "color_variable", "Color by:"),
            shiny::helpText("Variable to use for plot coloring"), omicNormDesignSelect(ns, "shape_variable", "Shape by:"),
            shiny::helpText("Variable to use for plot shapes")),
        itsd = list(shiny::hr(), shiny::h4("ITSD Normalization"),
            shiny::checkboxInput(ns("apply_itsd"), "Apply ITSD normalization", value = TRUE, width = "100%"),
            shiny::helpText("Normalize using internal standards (per assay)"),
            shiny::conditionalPanel(condition = "input.apply_itsd == true", ns = ns,
                shiny::selectInput(ns("itsd_aggregation"), "ITSD Aggregation:", choices = panels$itsd$aggregation,
                    selected = unname(panels$itsd$aggregation[[1L]]), width = "100%"),
                shiny::helpText("Method to combine multiple internal standards per sample. Select standards in ITSD Selection tab."))),
        # KIT-172: plain log2, no offset offered; a project whose record used one says so (the server's output).
        log2 = list(shiny::hr(), shiny::h4("Log2 Transform"),
            shiny::helpText(omicGuiLog2Texts()[["help"]]),
            shiny::uiOutput(ns("log2_recorded_offset")),
            shiny::helpText("Applied after ITSD normalization")),
        between_sample = list(shiny::hr(), shiny::h4("Between-Sample Normalization"),
            testid(shiny::selectInput(ns("norm_method"), "Normalization Method:", choices = panels$between_sample$methods,
                selected = panels$between_sample$default, width = "100%"), norm$testids[["method"]]),
            shiny::helpText(panels$between_sample$help %||% "Method for between-sample normalization")),
        tmm_logcpm = list(shiny::hr(), shiny::h4("TMM and logCPM"),
            shiny::helpText("edgeR TMM normalisation factors of the filtered counts; log2 counts per million with a prior count of 2.")),
        # KIT-129 W2: with the gate on, today's RUV controls leave the column for the RUV design section above it
        # (`omicNormRuvDesignSection`), where they are shown under today's design; else today's.
        ruv_iii = , ruv_grid = if (!isTRUE(norm$ruv_design$enabled)) omicNormRuvToday(ns, norm),
        actions = list(shiny::hr(),
            testid(shiny::actionButton(ns("run_normalization"), panels$actions$run_label %||% "Run Normalization Pipeline",
                class = "btn-primary", width = "100%",
                icon = shiny::icon("play")), norm$testids[["run"]]),
            shiny::br(), shiny::br(),
            omicNormResetButton(ns, if ("reset" %in% names(norm$testids)) norm$testids[["reset"]])),
        export = list(shiny::br(), shiny::br(),
            testid(shiny::actionButton(ns(panels$export$button), panels$export$label %||% "Export Normalized Data",
                class = "btn-success", width = "100%", icon = shiny::icon("download")), norm$testids[["export"]]),
            shiny::helpText(panels$export$help %||% "Save normalized data for later analysis")),
        # KIT-120: a download button; its handler is the omic's server's (proteomics' S4 compatibility file).
        s4_download = list(shiny::downloadButton(ns(panels$s4_download$button), panels$s4_download$label),
            shiny::helpText(panels$s4_download$help)))
    lapply(stats::setNames(names(panels), names(panels)), section)
}

#' The reset button, with its test id when the omic has one (KIT-110: trx's; lipid and metab have none).
#' @noRd
omicNormResetButton <- function(ns, testId = NULL) {
    button <- shiny::actionButton(ns("reset_normalization"), "Reset to Pre-Normalization", class = "btn-warning", width = "100%",
        icon = shiny::icon("undo"))
    if (is.null(testId)) button else testid(button, testId)
}

#' The RUV grid section (KIT-110, trx; mod_trx_norm.R:26-47 at 38d9d25e): RUVs' modes, and the automatic grid's and the
#' manual run's inputs with the grid's values. The input bounds and help texts are trx's, literal.
#' @noRd
omicNormRuvGridSection <- function(ns, grid, testId) {
    automatic <- shiny::conditionalPanel(condition = "input.ruv_mode == 'automatic'", ns = ns,
        shiny::h5("Automatic Optimization Parameters"),
        shiny::fluidRow(shiny::column(4, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = grid$percentage_min, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = grid$percentage_max, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_step"), "Step %:", value = grid$percentage_step, min = 1, max = 25, width = "100%"))),
        shiny::numericInput(ns("k_penalty_weight"), "K Value Penalty Weight:", value = grid$k_penalty_weight, min = 0, max = 0.9, step = 0.01,
            width = "100%"),
        shiny::helpText("The notebook's grid: 5 to 30 % by 5, k capped at 3 below 15 samples and 4 otherwise, penalty 0.05; the",
            "percentage and k with the maximum composite score are chosen."))
    manual <- shiny::conditionalPanel(condition = "input.ruv_mode == 'manual'", ns = ns,
        shiny::sliderInput(ns("ruv_percentage"), "% Genes as Inferred Negative Controls:", min = 1, max = 50, value = grid$percentage, step = 1,
            width = "100%"),
        shiny::numericInput(ns("ruv_k"), "Number of Factors (k):", value = grid$k, min = 1, max = 4, width = "100%"),
        shiny::helpText("k must not exceed the cap (3 below 15 samples, 4 otherwise); a larger k is refused."))
    list(shiny::hr(), shiny::h4("RUVs (RUVSeq)"),
        testid(shiny::radioButtons(ns("ruv_mode"), "RUV Parameter Tuning:", choices = grid$modes, selected = grid$default,
            width = "100%"), testId),
        shiny::helpText("RUVs uses inferred negative controls (the least significant genes in a preliminary group GLM). Its factors",
            "W_1..W_k enter the DA model as covariates; the RUVs-corrected logCPM is for display only."),
        automatic, manual)
}

#' The RUV-III section: the modes, the grouping, and the automatic, manual and skip panels
#' (mod_lipid_norm_ui_helpers.R:88-222 at f8d4c5fc), with the omic's control noun.
#' @noRd
omicNormRuvSection <- function(ns, ruv, testId) {
    # KIT-120: the omic's own help text by name (`ruv$help`), else the library's; a text only an omic has is shown
    # when it names one.
    help <- function(name, default = NULL) {
        text <- if (name %in% names(ruv$help)) ruv$help[[name]] else default
        if (is.null(text)) NULL else shiny::helpText(text)
    }
    list(shiny::hr(), shiny::h4("RUV-III Batch Correction"),
        testid(shiny::radioButtons(ns("ruv_mode"), "RUV Parameter Tuning:", choices = ruv$modes, selected = ruv$default,
            width = "100%"), testId),
        omicNormDesignSelect(ns, "ruv_grouping_variable", "RUV Grouping Variable:"),
        help("grouping", "Experimental factor to preserve during batch correction (shared across assays)"),
        shiny::conditionalPanel(condition = "input.ruv_mode == 'automatic'", ns = ns,
            shiny::h5("Automatic Optimization Parameters"),
            help("automatic", "Auto-optimization runs independently for each assay (pos/neg)"),
            shiny::fluidRow(
                shiny::column(6, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = 1, min = 1, max = 50, width = "100%")),
                shiny::column(6, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = 20, min = 1, max = 50, width = "100%"))),
            help("range", "Range of percentages to test (1-50%)"),
            shiny::selectInput(ns("separation_metric"), "Separation Quality Metric:", choices = list(
                "Maximum Difference (recommended)" = "max_difference", "Mean Difference" = "mean_difference",
                "Area Under Curve" = "auc", "Weighted Difference (deprecated)" = "weighted_difference"),
                selected = "max_difference", width = "100%"),
            help("metric"),
            shiny::sliderInput(ns("k_penalty_weight"), "K Value Penalty Weight:", min = 0.1, max = 0.9, value = 0.5, step = 0.1,
                width = "100%"),
            help("penalty", "Controls penalty for high k values (0.1 = lenient, 0.9 = strict)"),
            shiny::checkboxInput(ns("adaptive_k_penalty"), "Adaptive K Penalty (sample size aware)", value = TRUE, width = "100%"),
            help("adaptive"),
            shiny::br()),
        shiny::conditionalPanel(condition = "input.ruv_mode == 'manual'", ns = ns,
            shiny::sliderInput(ns("ruv_percentage"), sprintf("%% %s as Negative Controls:", ruv$noun), min = 1, max = 20,
                value = ruv$percentage, step = 1, width = "100%"),
            help("manual", sprintf("Percentage of most stable %s used as negative controls", tolower(ruv$noun))),
            shiny::numericInput(ns("ruv_k"), "Number of Factors (k):", value = ruv$k, min = 1, max = 10, width = "100%"),
            help("manual_k", "Number of unwanted variation factors to remove"),
            shiny::br()),
        shiny::conditionalPanel(condition = "input.ruv_mode == 'skip'", ns = ns,
            shiny::tags$div(style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; border-radius: 4px;",
                shiny::h5("RUV-III Will Be Skipped", style = "color: #856404; margin-top: 0;"),
                shiny::helpText(if ("skip" %in% names(ruv$help)) ruv$help[["skip"]] else
                    "Batch correction will not be applied. Consider including batch as a covariate in DE model if batch effects are present.",
                    style = "margin-bottom: 0;"),
                if ("skip_alternative" %in% names(ruv$help)) {
                    shiny::helpText(ruv$help[["skip_alternative"]], style = "margin-bottom: 0; font-style: italic;")
                }),
            shiny::br()))
}

# ------------------------------------------------------------------------------------------------
# KIT-129 W2: the RUV design section (ticket sections 1 and 2), shared by the four omics. It reads the
# configuration's presets (`norm$ruv_design`, utils_omic_gui_config.R) and is shown only with their gate on:
# with it off every function below returns today's markup or nothing. Presentation only: no input here is read
# by a computation yet. New ids, all static: the inputs `ruv_design`, `ruv_k_rule` and `ruv_k_user`, the output
# slot `ruv_design_summary` under the options (the factor-to-group association and a run's notices) and, in the
# RUV QC tab, one plot slot per assay row for the estimated count's diagnostic (`ruv_k_estimate_assay<n>`).
# The section is a block of the tab's full width above the options column and the QC tabset (the column is too
# narrow for it): the designs as rows of one comparison list (tier, design, false discoveries, power, note), the
# k rules side by side. Its styles are the `ruv-design` rules of app_ui.R's style block.
# ------------------------------------------------------------------------------------------------

#' Today's RUV section of an omic: the RUVSeq grid's where its panels name one, else RUV-III's.
#' @param ns The module's namespace function.
#' @param norm The omic's `norm` group.
#' @return A list of tags.
#' @noRd
omicNormRuvToday <- function(ns, norm) {
    if ("ruv_grid" %in% names(norm$panels)) {
        return(omicNormRuvGridSection(ns, norm$panels$ruv_grid, norm$testids[["ruv_mode"]]))
    }
    omicNormRuvSection(ns, norm$panels$ruv_iii, norm$testids[["ruv_mode"]])
}

#' A measured figure in its parts, as shown: the false discovery proportion with its scenario, and the power (the
#' share of real differences found, as a percentage). The column headings say which scenarios the figures are
#' from (the worst tested; the reference scenario), so a cell names the scenario only; a figure measured in one
#' scenario says so in both cells.
#' @param option A measured preset of the table, or another figure of one (`fdp`, `fdp_where`, `power`).
#' @return A named character vector: `fdp`, `fdp_where`, `power`, `power_where` (empty in the reference scenario).
#' @noRd
omicNormRuvFigureParts <- function(option) {
    single <- isTRUE(option$single_scenario)
    c(fdp = sprintf("%.1f%%", 100 * option$fdp),
        fdp_where = if (single) paste("only scenario tested:", option$fdp_where) else option$fdp_where,
        power = paste0(format(100 * option$power), "%"), power_where = if (single) "in that scenario" else "")
}

#' A preset's measured figure in words: the false discovery proportion with its cell and the power; a pending
#' or unavailable preset shows no number.
#' @param option A preset of the table.
#' @return One string.
#' @noRd
omicNormRuvFigureText <- function(option) {
    parts <- if (identical(option$status, "measured")) omicNormRuvFigureParts(option)
    switch(option$status,
        measured = sprintf("Measured in simulation: %s false discoveries at a 5%% FDR threshold (%s); power %s %s.",
            parts[["fdp"]], if (isTRUE(option$single_scenario)) parts[["fdp_where"]] else
                paste("worst scenario tested:", parts[["fdp_where"]]), parts[["power"]],
            if (isTRUE(option$single_scenario)) parts[["power_where"]] else "in the reference scenario"),
        pending = "Not yet measured in simulation.",
        unavailable = paste("Not available:", option$unavailable))
}

#' The figures of a preset measured with another rule for k than its pre-selected one, in words, one per rule.
#' @param option A preset of the table.
#' @param k_rules The k rules' texts (the group's `k_rules`), for their names.
#' @return A character vector, empty when the preset has no such figure.
#' @noRd
omicNormRuvAltText <- function(option, k_rules = omicGuiRuvKRules()) {
    vapply(names(option$k_rule_figures), \(rule) {
        parts <- omicNormRuvFigureParts(option$k_rule_figures[[rule]])
        sprintf("With the %s for k: %s false discoveries; worst scenario: %s. Power %s.", tolower(k_rules[[rule]]$label),
            parts[["fdp"]], parts[["fdp_where"]], parts[["power"]])
    }, character(1), USE.NAMES = FALSE)
}

#' One preset as a radio choice, a row of the comparison list: its tier, its label with its description under it,
#' the figure and the note; the tooltip (the exact statement, extra detail only) is the title. A measured figure is shown in two cells (false discoveries, power); its sentence is the
#' text a screen reader is given in their place. A figure measured with the other rule for k is a line under them,
#' in the same cell group (`ruv-design-figures`, which spans the two figure columns).
#' @param option A preset of the table.
#' @param tiers The tier labels.
#' @return A tag.
#' @noRd
omicNormRuvOptionTag <- function(option, tiers) {
    cell <- \(name, ...) shiny::tags$span(class = paste0("ruv-design-", name), ...)
    figure <- if (identical(option$status, "measured")) {
        parts <- omicNormRuvFigureParts(option)
        list(shiny::tags$span(class = "sr-only", omicNormRuvFigureText(option)),
            cell("fdp", `aria-hidden` = "true", shiny::tags$strong(parts[["fdp"]]), cell("where", parts[["fdp_where"]])),
            cell("power", `aria-hidden` = "true", shiny::tags$strong(parts[["power"]]),
                if (nzchar(parts[["power_where"]])) cell("where", parts[["power_where"]])))
    } else {
        cell("status", omicNormRuvFigureText(option))
    }
    # KIT-129 fourth pass item 2: an option the run path does not take here, or whose precondition the design does
    # not meet, keeps its figure and says why it cannot be chosen.
    reason <- if (identical(option$status, "measured")) option$unavailable_now
    shiny::tags$span(title = option$tooltip,
        class = paste(c("ruv-design-option", if (!identical(option$status, "measured")) paste0("is-", option$status),
            if (!is.null(reason)) "is-unavailable"), collapse = " "),
        cell("tier", shiny::tags$strong(class = paste0("ruv-design-tier-", gsub("_", "-", option$tier, fixed = TRUE)),
            tiers[[option$tier]], shiny::tags$span(class = "sr-only", ":"))),
        cell("label", option$label, cell("description", option$description)),
        cell("figures", figure,
            if (length(option$k_rule_figures) > 0L) cell("alt", lapply(omicNormRuvAltText(option), shiny::tags$span)),
            if (!is.null(reason)) cell("status", paste("Not available:", reason))),
        cell("note", option$note))
}

#' The comparison list's column headings, above the options: the two figure columns say what is measured and in
#' which scenario, so the rows name the scenario only. They are for the eye: each option's own text carries its
#' tier and its figure in words.
#' @return A tag.
#' @noRd
omicNormRuvOptionHeader <- function() {
    cell <- \(name, ...) shiny::tags$span(class = paste0("ruv-design-", name), ...)
    shiny::tags$div(class = "ruv-design-header", `aria-hidden` = "true", cell("tier", "Tier"), cell("label", "Approach"),
        cell("fdp", "False discoveries", cell("where", "target 5%; worst scenario tested")),
        cell("power", "Power", cell("where", "real differences found; reference scenario")),
        cell("note", "When to use it, and the evidence"))
}

#' A JavaScript condition for `conditionalPanel`: whether an input's value is one of `values`.
#' @param input_id The input's id, without its namespace.
#' @param values The values.
#' @return One string.
#' @noRd
omicNormRuvJsIn <- function(input_id, values) {
    sprintf("[%s].indexOf(input.%s) >= 0", paste0("'", values, "'", collapse = ", "), input_id)
}

#' Every preset of a design group's tables, once per id.
#' @param design The norm group's `ruv_design`.
#' @return A list of presets named by id.
#' @noRd
omicNormRuvPresets <- function(design) {
    options <- unlist(lapply(unname(design$tables), `[[`, "options"), recursive = FALSE)
    ids <- vapply(options, `[[`, character(1), "id")
    stats::setNames(options, ids)[unique(ids)]
}

#' The radio group of a table's presets under the list's column headings; a preset that is not available is
#' disabled.
#' @param ns The module's namespace function.
#' @param design The norm group's `ruv_design`.
#' @param table The table shown.
#' @return A tag.
#' @noRd
omicNormRuvDesignRadio <- function(ns, design, table) {
    # KIT-129 fourth pass item 2: as the tab opens, before the session's design is read, no technical replicates.
    table <- omicNormRuvMarkTable(table, technical = FALSE)
    radio <- shiny::radioButtons(ns("ruv_design"), "How unwanted variation is removed:",
        choiceNames = lapply(table$options, omicNormRuvOptionTag, tiers = design$tiers),
        choiceValues = vapply(table$options, `[[`, character(1), "id"), selected = table$selected, width = "100%")
    for (option in table$options) {
        if (!identical(option$status, "unavailable") && is.null(option$unavailable_now)) next
        radio <- htmltools::tagQuery(radio)$find("input")$
            filter(\(input, i) identical(htmltools::tagGetAttribute(input, "value"), option$id))$
            addAttrs(disabled = NA)$allTags()
    }
    radio <- htmltools::tagQuery(radio)$find(".shiny-options-group")$before(omicNormRuvOptionHeader())$allTags()
    testid(radio, design$testid)
}

#' The choice of k: the rule's radio group (its choices follow the selected design; `registerOmicNormRuvDesign`)
#' and the number, always visible; behind a disclosure that is closed until opened (a native `details` element,
#' keyboard operable), each rule's reading, strengths, limits and diagnostic shown where the rule is defined, the
#' rules side by side with their rows aligned.
#' @param ns The module's namespace function.
#' @param design The norm group's `ruv_design`.
#' @return A list of tags.
#' @noRd
omicNormRuvKSection <- function(ns, design) {
    presets <- omicNormRuvPresets(design)
    offering <- \(rule) names(Filter(\(preset) rule %in% preset$k_rules, presets))
    with_k <- names(Filter(\(preset) length(preset$k_rules) > 0L, presets))
    # The choices rendered first are the selected design's; today's design has none, so the recommended one's.
    table <- design$tables[[1L]]
    first <- presets[[if (table$selected %in% with_k) table$selected else table$recommended]]
    first_rules <- setdiff(first$k_rules, names(table$k_rules_not_offered))
    labels <- vapply(design$k_rules, `[[`, character(1), "label")
    texts <- table$texts
    rule_text <- lapply(names(design$k_rules), \(rule) {
        text <- design$k_rules[[rule]]
        shiny::conditionalPanel(condition = omicNormRuvJsIn("ruv_design", offering(rule)), ns = ns,
            class = "ruv-design-k-rule",
            shiny::tags$p(class = "ruv-design-k-rule-name", shiny::tags$strong(text$label)), shiny::tags$p(text$reads),
            shiny::tags$p(shiny::tags$strong("Strengths: "), text$pros),
            shiny::tags$p(shiny::tags$strong("Limits: "), text$cons),
            shiny::tags$p(shiny::tags$strong("Plot to check: "), text$diagnostic))
    })
    pending <- unlist(lapply(presets, \(preset) sprintf("(input.ruv_design == '%s' && input.ruv_k_rule == '%s')",
        rep(preset$id, length(preset$k_rules_pending)), preset$k_rules_pending)))
    list(
        shiny::conditionalPanel(condition = omicNormRuvJsIn("ruv_design", with_k), ns = ns, class = "ruv-design-k",
            shiny::h5("Number of factors (k)"),
            shiny::helpText("k is the number of separate patterns of unwanted variation (factors) that are estimated",
                "and removed."),
            shiny::tags$div(class = "ruv-design-k-choice",
                shiny::radioButtons(ns("ruv_k_rule"), "How the number of factors is chosen:",
                    choices = stats::setNames(first_rules, labels[first_rules]), selected = first[["k_rule"]],
                    width = "100%"),
                shiny::conditionalPanel(condition = "input.ruv_k_rule == 'user'", ns = ns,
                    shiny::numericInput(ns("ruv_k_user"), "Number of factors (k):", value = NA, min = 1, step = 1,
                        width = "100%")),
                # KIT-129 fourth pass item 2: the seed of the estimated number (also counted beside your own number).
                shiny::conditionalPanel(condition = "input.ruv_k_rule != 'sweep_bounded'", ns = ns,
                    shiny::numericInput(ns("ruv_k_seed"), texts[["seed"]], value = RUV_FACTOR_COUNT_SEED, step = 1,
                        width = "100%"))),
            if (length(table$k_rules_not_offered) > 0L) shiny::conditionalPanel(
                condition = omicNormRuvJsIn("ruv_design", intersect(with_k, offering(names(table$k_rules_not_offered)))),
                ns = ns, shiny::helpText(paste(texts[["k_rule_not_available"]], table$k_rules_not_offered[[1L]]))),
            if (length(pending) > 0L) shiny::conditionalPanel(condition = paste(pending, collapse = " || "), ns = ns,
                shiny::helpText("No figure has been measured yet for this approach with this rule. The figure shown",
                    "for the approach in the list above was measured with its pre-selected rule.")),
            shiny::tags$details(class = "ruv-design-k-why",
                shiny::tags$summary("How is this chosen? What each rule does, its strengths and its limits"),
                shiny::tags$div(class = "ruv-design-k-rules", rule_text),
                shiny::conditionalPanel(condition = omicNormRuvJsIn("ruv_design", setdiff(with_k, offering("sweep_bounded"))),
                    ns = ns, shiny::helpText("The separation rule is not offered when all features are the controls: the",
                        "two curves it compares are then the same."))),
            shiny::helpText("Whichever rule you use, its diagnostic plot is shown in the RUV QC tab and saved with the",
                "QC outputs.")))
}

#' The RUV design section, a block of the tab's full width: every preset with its figure and note, the output slot
#' under them, the choice of k, and today's RUV controls shown under today's design. Nothing with the gate off:
#' today's section is then in the options column (`omicNormOptionSections`).
#' @param ns The module's namespace function.
#' @param norm The omic's `norm` group.
#' @return A tag, or NULL.
#' @noRd
omicNormRuvDesignSection <- function(ns, norm) {
    design <- norm$ruv_design
    if (!isTRUE(design$enabled)) return(NULL)
    presets <- omicNormRuvPresets(design)
    legacy <- names(Filter(\(preset) identical(preset$tier, "legacy"), presets))
    shiny::tags$section(class = "well ruv-design", `aria-label` = "Unwanted variation (RUV)",
        shiny::h4("Unwanted variation (RUV)"),
        shiny::helpText("Unwanted variation is variation between samples that comes from the measurement rather than",
            "the biology, such as batch, run order or sample handling. RUV (Remove Unwanted Variation) estimates it",
            "as a small number of factors so that it does not distort the differential analysis."),
        omicNormRuvDesignRadio(ns, design, design$tables[[1L]]),
        shiny::helpText(shiny::tags$strong("How to read the figures."), "They were measured in simulation, not on your",
            "data. False discoveries: the average share of the features reported as different that are not, at a 5%",
            "false discovery rate (FDR) threshold, so 5% is the target; shown for the worst scenario tested that",
            "contained real differences. Power: the share of real differences found, in the reference scenario (6",
            "samples per group, 20% of features truly different, one unwanted factor). Each figure was measured with",
            "the approach's pre-selected rule for the number of factors (k); where the other rule was also measured,",
            "its figure is on the line beneath. \"Tied to the groups\": the",
            "unwanted variation differs between the groups, for example when most of one group was run in the same",
            "batch."),
        if (!is.null(design$conditions)) shiny::helpText(design$conditions),
        shiny::helpText("Whichever approach you choose, the corrected values are still calculated and used for the",
            "plots, the exported tables and the multi-omics integration. The approach decides what the differential",
            "test uses: in the recommended approach and those valid with conditions, the test is run on the",
            "uncorrected values with the estimated factors included as covariates."),
        # KIT-129 fourth pass item 2: the supplied list, read by the run under that design.
        shiny::conditionalPanel(condition = "input.ruv_design == 'supplied_controls'", ns = ns,
            shiny::textAreaInput(ns("ruv_control_list"), design$tables[[1L]]$texts[["control_list"]], value = "",
                rows = 4, width = "100%")),
        shiny::uiOutput(ns("ruv_design_summary")),
        omicNormRuvKSection(ns, design),
        shiny::conditionalPanel(condition = omicNormRuvJsIn("ruv_design", legacy), ns = ns, class = "ruv-design-today",
            omicNormRuvToday(ns, norm)))
}

#' The RUV QC tab's static plot slots for the estimated count's diagnostic, one per assay row, shown under the
#' designs that choose a k. Nothing with the gate off.
#' @param ns The module's namespace function.
#' @param view The omic's `norm` group.
#' @return A tag, or NULL.
#' @noRd
omicNormRuvKDiagnosticSlots <- function(ns, view) {
    design <- view$ruv_design
    if (!isTRUE(design$enabled)) return(NULL)
    with_k <- names(Filter(\(preset) length(preset$k_rules) > 0L, omicNormRuvPresets(design)))
    with_separation <- names(Filter(\(preset) "sweep_bounded" %in% preset$k_rules, omicNormRuvPresets(design)))
    # KIT-129 fifth pass item 5: the check that a factor follows the groups, one column (ticket section 3 (2))
    list(shiny::tableOutput(ns("ruv_design_association_table")),
        shiny::conditionalPanel(condition = omicNormRuvJsIn("ruv_design", with_k), ns = ns, shiny::hr(),
        shiny::h5("Number of factors (k): the estimated number", style = "text-align: center;"),
        shiny::helpText(design$k_rules$num_sv_be$diagnostic),
        lapply(seq_len(view$assay_slots), \(n) shinyjqui::jqui_resizable(
            shiny::plotOutput(ns(sprintf("ruv_k_estimate_assay%d", n)), height = "400px")))),
        # KIT-129 fourth pass item 3: the separation rule's plot of the last run, where the run computed one.
        shiny::conditionalPanel(condition = paste(omicNormRuvJsIn("ruv_design", with_separation), "&& input.ruv_k_rule == 'sweep_bounded'"),
            ns = ns, shiny::hr(),
            shiny::h5(design$tables[[1L]]$texts[["separation_plot_heading"]], style = "text-align: center;"),
            shiny::helpText(design$k_rules$sweep_bounded$diagnostic),
            # KIT-129 fifth pass item 4: one plot per assay row where the omic has several (lipid/metab)
            if ((view$assay_slots %||% 1L) > 1L) lapply(seq_len(view$assay_slots), \(n) shinyjqui::jqui_resizable(
                shiny::plotOutput(ns(sprintf("ruv_k_separation_assay%d", n)), height = "400px"))) else
                shinyjqui::jqui_resizable(shiny::plotOutput(ns("ruv_k_separation_plot"), height = "400px"))))
}

#' The table of a design group for a workflow type: the route's where the group names one, else the first.
#' @param design The norm group's `ruv_design`.
#' @param manager The bus's state manager (proteomics' carries `workflow_type`); anything else has no route.
#' @return The table's name.
#' @noRd
omicNormRuvTableName <- function(design, manager) {
    route <- if (is.environment(manager)) manager$workflow_type else NULL
    if (checkmate::test_string(route) && route %in% names(design$route_tables)) {
        return(design$route_tables[[route]])
    }
    names(design$tables)[[1L]]
}

#' Replace the presets shown by another table's (a route with its own figures), with that table's selection.
#' @param session The module's session.
#' @param design The norm group's `ruv_design`.
#' @param table The table to show.
#' @return NULL, invisibly.
#' @noRd
omicNormRuvShowTable <- function(session, design, table) {
    # (KIT-129 fourth pass item 2: `table` is marked, `omicNormRuvMarkTable`; its selection is the default there.)
    shiny::updateRadioButtons(session, "ruv_design",
        choiceNames = lapply(table$options, omicNormRuvOptionTag, tiers = design$tiers),
        choiceValues = vapply(table$options, `[[`, character(1), "id"), selected = table$selected)
    invisible(NULL)
}

#' Offer the k rules of the selected design, with the rule its figure was measured with selected. A design
#' without a k of its own (today's, none) leaves the choice as it is: its panel is hidden.
#' @param session The module's session.
#' @param design The norm group's `ruv_design`.
#' @param table The table shown.
#' @param selected The selected design's id.
#' @return NULL, invisibly.
#' @noRd
omicNormRuvShowKRules <- function(session, design, table, selected) {
    option <- Find(\(option) identical(option$id, selected), table$options)
    if (length(option$k_rules) == 0L) return(invisible(NULL))
    # KIT-129 fourth pass item 2: a rule the run path does not take on this route is not offered.
    rules <- setdiff(option$k_rules, names(table$k_rules_not_offered))
    labels <- vapply(design$k_rules[rules], `[[`, character(1), "label")
    shiny::updateRadioButtons(session, "ruv_k_rule", choices = stats::setNames(rules, unname(labels)),
        selected = option[["k_rule"]])
    invisible(NULL)
}

#' The section's two observers: the table of the workflow type once it is known (an omic with tables by route),
#' and the k rules of the selected design. None with the gate off.
#' @param input,session The module's.
#' @param workflowData The bus.
#' @param design The norm group's `ruv_design`.
#' @return NULL, invisibly.
#' @noRd
registerOmicNormRuvDesign <- function(input, session, workflowData, design, output = NULL, normData = NULL, slots = 1L) {
    if (!isTRUE(design$enabled)) return(invisible(NULL))
    if (!is.null(output) && !is.null(normData)) registerOmicNormRuvDesignOutputs(output, normData, design, slots)
    # KIT-129 fourth pass item 2: the table of the route and whether the session's design declares technical
    # replicates, read again on each state change; the options are shown again, with their default selection, only
    # when what is available changes.
    shown_table <- shiny::reactive({
        # The manager's workflow type is not reactive; the bus's trigger marks a state change (proteomics' idiom).
        workflowData$state_update_trigger
        name <- omicNormRuvTableName(design, workflowData$state_manager)
        list(name = name, table = omicNormRuvMarkTable(design$tables[[name]],
            omicNormRuvHeadTechnicalReplicates(workflowData$state_manager)))
    })
    shown <- shiny::reactiveVal(omicNormRuvAvailabilityKey(names(design$tables)[[1L]],
        omicNormRuvMarkTable(design$tables[[1L]], technical = FALSE)))
    shiny::observe({
        current <- shown_table()
        key <- omicNormRuvAvailabilityKey(current$name, current$table)
        if (identical(key, shiny::isolate(shown()))) return(invisible(NULL))
        omicNormRuvShowTable(session, design, current$table)
        shown(key)
    })
    shiny::observeEvent(input$ruv_design, {
        table <- shown_table()$table
        option <- Find(\(option) identical(option$id, input$ruv_design), table$options)
        reason <- option$unavailable_now %||% (if (identical(option$status, "unavailable")) option$unavailable)
        if (!is.null(reason)) {
            # Not selectable here: said, and the default of this design selected again.
            shiny::showNotification(paste(table$texts[["not_available"]], reason), type = "warning", duration = 10)
            shiny::updateRadioButtons(session, "ruv_design", selected = table$selected)
            return(invisible(NULL))
        }
        omicNormRuvShowKRules(session, design, table, input$ruv_design)
    })
    invisible(NULL)
}

#' The section's outputs of the last run (KIT-129 fourth pass item 3): under the options, what was done
#' (`ruv_design_summary`); in the RUV QC tab the estimated number's diagnostic per assay row (`plotRuvFactorCount`) and
#' the separation rule's plot when the run computed one. They read `normData$ruv_design_result`
#' (`omicNormRuvDesignResult`), which a run of a new design sets and a run of a legacy design clears.
#' @param output The module's output.
#' @param normData The tab's state.
#' @param design The norm group's `ruv_design`.
#' @param slots The tab's assay rows (`norm$assay_slots`).
#' @return NULL, invisibly.
#' @noRd
registerOmicNormRuvDesignOutputs <- function(output, normData, design, slots = 1L) {
    presets <- omicNormRuvPresets(design)
    texts <- design$tables[[1L]]$texts
    output$ruv_design_summary <- shiny::renderUI(omicNormRuvResultTags(normData$ruv_design_result, presets, texts))
    # KIT-129 fifth pass item 5: the RUV QC tab's column of the check, from the same recorded values
    output$ruv_design_association_table <- shiny::renderTable({
        rows <- omicNormRuvAssociationRows(normData$ruv_design_result, texts)
        shiny::req(nrow(rows) > 0L)
        rows
    }, striped = TRUE, width = "100%")
    for (n in seq_len(slots %||% 1L)) local({
        slot <- n
        output[[sprintf("ruv_k_estimate_assay%d", slot)]] <- shiny::renderPlot({
            result <- normData$ruv_design_result
            estimate <- if (length(result$assays) >= slot) result$assays[[slot]]$estimate
            shiny::req(is.list(estimate))
            plotRuvFactorCount(estimate)
        })
    })
    output$ruv_k_separation_plot <- shiny::renderPlot({
        plot <- normData$ruv_design_result$separation_plot
        shiny::req(!is.null(plot))
        if (inherits(plot, "grob")) grid::grid.draw(plot) else print(plot)
    })
    # KIT-129 fifth pass item 4: per assay row, the separation rule's plot of a lipid/metab run that computed it
    if ((slots %||% 1L) > 1L) for (n in seq_len(slots)) local({
        slot <- n
        output[[sprintf("ruv_k_separation_assay%d", slot)]] <- shiny::renderPlot({
            result <- normData$ruv_design_result
            plot <- if (length(result$assays) >= slot) result$assays[[slot]]$separation_plot
            shiny::req(!is.null(plot))
            if (inherits(plot, "grob")) grid::grid.draw(plot) else print(plot)
        })
    })
    invisible(NULL)
}

#' What a run of a new design did, in the panel's words (KIT-129 fourth pass item 3): the approach, per assay k with the
#' rule that gave it (the estimated number beside a typed k, decision M33), or why no factor was removed; the replicate
#' sets; the seed.
#' @param result `omicNormRuvDesignResult`'s value, or NULL (nothing shown).
#' @param presets The design group's presets by id.
#' @param texts The table's texts.
#' @return A tag, or NULL.
#' @noRd
omicNormRuvResultTags <- function(result, presets, texts) {
    if (is.null(result)) return(NULL)
    line <- \(...) shiny::tags$p(...)
    text <- \(key) if (checkmate::test_string(key) && key %in% names(texts)) texts[[key]]
    label <- presets[[result$design]]$label %||% result$design
    assays <- result$assays
    several <- length(assays) > 1L
    k_lines <- lapply(names(assays), function(name) {
        assay <- assays[[name]]
        prefix <- if (several) paste(texts[["result_assay"]], name, "-") else NULL
        if (!isTRUE(assay$applied)) return(line(paste(c(prefix, texts[["result_not_applied"]], assay$reason), collapse = " ")))
        rule <- text(paste0("result_rule_", result$k_rule)) %||% result$k_rule
        beside <- if (identical(result$k_rule, "user")) {
            if (!is.na(assay$k_estimate %||% NA)) paste0("; ", texts[["result_k_estimate"]], " ", assay$k_estimate) else
                if (!is.null(assay$k_estimate_reason) && !is.na(assay$k_estimate_reason))
                    paste0("; ", texts[["result_k_estimate_none"]], " ", assay$k_estimate_reason)
        }
        line(paste(c(prefix, texts[["result_k"]], paste0(assay$k, " (", rule, beside, ")")), collapse = " "))
    })
    # KIT-129 fifth pass item 3: a supplied list, how much of it the data holds and the controls used per assay.
    list_lines <- if (!is.null(result$control_list)) c(list(line(sprintf(texts[["result_list"]], result$control_list$listed,
        result$control_list$listed - result$control_list$unmatched))), lapply(names(assays), function(name) {
            if (is.null(assays[[name]]$control_count)) return(NULL)
            line(paste(c(if (several) paste(texts[["result_assay"]], name, "-"), texts[["result_controls"]], assays[[name]]$control_count),
                collapse = " "))
        }))
    sets <- text(paste0("result_sets_", result$replicate_structure))
    # KIT-129 fourth pass item 4 (sixth pass, M35, M36): the check that a factor follows the groups, per assay: detectable
    # or not against the approach's detection limit, the value without a verdict where none applies, and under a legacy
    # design that it is not available.
    association_lines <- lapply(names(assays), function(name) {
        a <- assays[[name]]$association
        words <- omicNormRuvAssociationText(a, texts, result$route)
        if (is.null(words)) return(NULL)
        prefix <- if (several) paste(texts[["result_assay"]], name, "-") else NULL
        if (identical(a$ruv_group_association_outcome, "detectable")) {
            return(shiny::tags$p(class = "ruv-design-warning", role = "alert", shiny::tags$strong(paste(c(prefix, words), collapse = " "))))
        }
        line(paste(c(prefix, words), collapse = " "))
    })
    shiny::tags$div(class = "ruv-design-result", `aria-live` = "polite",
        shiny::tags$h5(texts[["result_heading"]]),
        line(texts[["result_design"]], label),
        if (!identical(result$design, "none")) k_lines,
        list_lines,
        if (!identical(result$design, "none") && !is.null(sets)) line(sets),
        if (!is.na(result$seed %||% NA) && result$k_rule %in% c("num_sv_be", "user")) line(texts[["result_seed"]], result$seed),
        association_lines,
        if (isTRUE(result$legacy)) line(texts[["association_legacy"]]))
}

#' The check that a factor follows the groups in the panel's words (KIT-129 fourth pass item 4; fifth pass item 5: shared
#' by the panel, the RUV QC tab's table and the DA tab's line; sixth pass, M35 and M36: a detection limit with the
#' approach's measured consequences): the recorded value and limit to 3 decimals, the false discovery proportions as
#' percentages to one decimal ("not measured" where the cell was not measured); NULL without a value or for another
#' outcome.
#' @param a The recorded check (`ruvGroupAssociationRecordFields`). @param texts The table's texts.
#' @param route The result's route (`dia_limpa` adds its sentence).
#' @noRd
omicNormRuvAssociationText <- function(a, texts, route = NA_character_) {
    if (is.null(a) || is.na(a$ruv_group_association %||% NA)) return(NULL)
    value <- sprintf("%.3f", a$ruv_group_association)
    limit <- sprintf("%.3f", a$ruv_group_association_threshold %||% NA_real_)
    percent <- function(field) {
        x <- a[[field]] %||% NA_real_
        if (is.na(x)) texts[["association_not_measured"]] else sprintf("%.1f%%", 100 * x)
    }
    typed <- if (identical(a$ruv_group_association_basis, "typed_k_equals_estimated_number")) texts[["association_typed_k_equals"]]
    words <- switch(a$ruv_group_association_outcome %||% "",
        detectable = sprintf(texts[["association_detectable"]], value, limit, percent("ruv_group_association_fdp_0"),
            percent("ruv_group_association_fdp_03"), percent("ruv_group_association_fdp_06"),
            percent("ruv_group_association_uncorrected_06")),
        not_detectable = sprintf(texts[["association_not_detectable"]], value, limit),
        no_limit = if (identical(route, "dia_limpa")) sprintf(texts[["association_dia_limpa"]], value) else
            sprintf(texts[["association_no_limit"]], value),
        NULL)
    if (is.null(words)) return(NULL)
    paste(c(words, typed), collapse = " ")
}

#' The RUV QC tab's table of the check (KIT-129 fifth pass item 5, ticket section 3 (2)): one row per assay with a
#' recorded value, its column the panel's words for it; the assay column only where the omic has several assays.
#' @param result `omicNormRuvDesignResult`'s value, or NULL. @param texts The table's texts.
#' @return A data frame, with no row when there is nothing to show.
#' @noRd
omicNormRuvAssociationRows <- function(result, texts) {
    assays <- result$assays %||% list()
    words <- lapply(assays, function(assay) omicNormRuvAssociationText(assay$association, texts, result$route))
    shown <- names(Filter(Negate(is.null), words))
    rows <- stats::setNames(data.frame(unlist(words[shown], use.names = FALSE) %||% character(), stringsAsFactors = FALSE),
        texts[["association_column"]])
    if (length(assays) > 1L) rows <- cbind(stats::setNames(data.frame(shown, stringsAsFactors = FALSE), texts[["result_assay"]]), rows)
    rows
}

#' One run's result for the panel (KIT-129 fourth pass item 3), in one shape for every omic.
#' @param design The design's id. @param k_rule The k rule. @param seed The count's seed, or NA.
#' @param replicate_structure `group`, `technical_replicates` or NA.
#' @param assays Per assay (one, `all`, outside lipid/metab): `k`, `k_estimate` (beside a typed k), `k_estimate_reason`,
#'   `applied`, `reason`, `estimate` (the count's record for `plotRuvFactorCount`, or NULL).
#' @param separation_plot The separation rule's plot of the run, or NULL.
#' @param control_list A supplied list (KIT-129 fifth pass item 3): `list(listed, unmatched)`, the counts; NULL otherwise.
#'   Its assays carry `control_count`, the controls used.
#' @noRd
omicNormRuvDesignResult <- function(design, k_rule, seed, replicate_structure, assays, separation_plot = NULL, route = NA_character_,
    legacy = FALSE, control_list = NULL) {
    list(design = design, k_rule = k_rule %||% NA_character_, seed = seed %||% NA_integer_,
        replicate_structure = replicate_structure %||% NA_character_, assays = assays, separation_plot = separation_plot, route = route,
        legacy = legacy, control_list = control_list)
}

#' A product's supplied list for the panel (KIT-129 fifth pass item 3): its counts, or NULL without one.
#' @param record The product's record of the list (`ruv_control_list_count`, `ruv_control_list_unmatched`).
#' @noRd
omicNormRuvControlList <- function(record) {
    if (is.null(record$ruv_control_list_count)) return(NULL)
    list(listed = as.integer(record$ruv_control_list_count), unmatched = length(record$ruv_control_list_unmatched))
}

#' A run of a legacy design, for the panel (KIT-129 fourth pass item 4): its name and that the check is not available.
#' @noRd
omicNormRuvLegacyResult <- function(design = "legacy") omicNormRuvDesignResult(design, NA, NA, NA, list(), legacy = TRUE)

#' Set the tab's last-run result (KIT-129 fourth pass item 3); nothing with the gate off.
#' @param normData The tab's state. @param result The result, or NULL (a legacy run: the summary is cleared).
#' @noRd
omicNormSetRuvDesignResult <- function(normData, result) {
    if (!omicGuiRuvDesignEnabled() || is.null(normData)) return(invisible(NULL))
    normData$ruv_design_result <- result
    invisible(result)
}

#' The result of the classic normalisation's product (KIT-129 fourth pass item 3); NULL for today's design.
#' @param saved The product's value (`diaNormalizationOperation`).
#' @noRd
protNormRuvDesignResult <- function(saved) {
    design <- saved$parameters$ruv_design
    if (is.null(design)) return(omicNormRuvLegacyResult())
    r <- saved$ruv_result
    if (identical(design, "none")) return(omicNormRuvDesignResult("none", NA, NA, NA, list()))
    applied <- isTRUE(r$ruv_applied)
    omicNormRuvDesignResult(design, r$ruv_k_rule, r$ruv_k_seed, saved$ruv_model$replicate_structure,
        list(all = list(k = if (applied) saved$best_k else 0L, k_estimate = r$ruv_k_estimate, k_estimate_reason = r$ruv_k_estimate_reason,
            applied = applied, reason = r$ruv_not_applied_reason, estimate = r$k_estimate, association = r$ruv_group_association,
            control_count = if (!is.null(r$ruv_control_list_count)) saved$ruv_model$control_count)),
        separation_plot = if (identical(r$ruv_k_rule, "sweep_bounded")) r$best_cancor_plot, control_list = omicNormRuvControlList(r))
}

#' The result of DIA with limpa's RUV covariates product (KIT-129 fourth pass item 3); NULL for today's rule.
#' @param saved The product's value. @param design The panel's design id (the route's product names its control rule only).
#' @noRd
diaRuvDesignResult <- function(saved, design) {
    if (is.null(saved$ruv_control_rule)) return(omicNormRuvLegacyResult(if (identical(saved$ruv_mode, "correction")) "legacy_correction" else "legacy"))
    r <- saved$ruv_result
    applied <- isTRUE(saved$ruv_applied)
    omicNormRuvDesignResult(design %||% saved$ruv_control_rule, saved$ruv_k_rule, saved$parameters$ruv_k_seed,
        saved$replicate_structure, list(all = list(k = saved$best_k, k_estimate = saved$ruv_k_estimate,
            k_estimate_reason = r$ruv_k_estimate_reason, applied = applied, reason = saved$ruv_not_applied_reason, estimate = r$k_estimate,
            association = r$ruv_group_association)),
        separation_plot = r$k_sweep$cancor_plot, route = "dia_limpa")
}

#' The result of the lipid/metab new design's RUV product (KIT-129 fourth pass item 3), per assay; NULL for today's.
#' @param value The product's value (`omicsRuvDesignOperation`).
#' @noRd
omicsRuvDesignResult <- function(value) {
    design <- value$ruv_design
    if (is.null(design)) return(omicNormRuvLegacyResult())
    if (identical(design, "none")) return(omicNormRuvDesignResult("none", NA, NA, NA, list()))
    assays <- lapply(stats::setNames(names(value$factors), names(value$factors)), function(assay) {
        record <- value$factors[[assay]]
        list(k = record$k, k_estimate = record$ruv_k_estimate, k_estimate_reason = record$k_estimate_reason,
            applied = isTRUE(value$ruv_applied[[assay]]), reason = record$reason, estimate = value$k_estimate[[assay]],
            association = record$ruv_group_association,
            control_count = if (!is.null(value$ruv_control_list_count)) value$control_count[[assay]],
            separation_plot = value$separation_plots[[assay]])
    })
    omicNormRuvDesignResult(design, value$ruv_k_rule, value$ruv_k_seed, value$ruv_replicate_structure, assays,
        control_list = omicNormRuvControlList(value))
}

#' The result of the transcriptomics declared choice (KIT-129 fourth pass item 3).
#' @param value The choice's value (`trxChooseRuvDesign`, the declared `trx_ruv_grid` product).
#' @param outcomes The RUVs step's outcomes (`trx_ruvs`; item 4: the factor-tied-to-groups check), or NULL.
#' @noRd
trxRuvDesignResult <- function(value, outcomes = NULL) {
    d <- value[["ruv_design", exact = TRUE]]
    if (is.null(d)) return(omicNormRuvLegacyResult())
    if (identical(d$id, "none")) return(omicNormRuvDesignResult("none", NA, NA, NA, list()))
    omicNormRuvDesignResult(d$id, d$k_rule, value$k_record$k_seed %||% value$settings$ruv_k_seed,
        if (identical(d$method, "RUVs")) "group", list(all = list(k = d$k, k_estimate = value$k_record$k_estimate,
            k_estimate_reason = value$k_record$k_estimate_reason, applied = isTRUE(d$applied), reason = d$reason,
            estimate = value$k_estimate, association = if (!is.null(outcomes$ruv_group_association_outcome))
                outcomes[ruvGroupAssociationRecordFields])),
        # KIT-129 fifth pass item 4: the separation rule's plot, where the run computed the rule
        separation_plot = value[["separation_plot", exact = TRUE]])
}

#' Why an option cannot be chosen on this session (KIT-129 fourth pass item 2): the reason the run path does not take
#' it on this route, or its precondition on the design (technical replicates declared, or not); NULL when it can.
#' An option the table marks unavailable carries its own reason (`unavailable`) and is disabled as rendered.
#' @param option A preset of the table.
#' @param table Its table (its `texts`).
#' @param technical Whether the session's design declares technical replicates.
#' @return One string, or NULL.
#' @noRd
omicNormRuvUnavailableReason <- function(option, table, technical) {
    if (identical(option$status, "unavailable")) return(NULL)
    if (!is.null(option$not_offered)) return(option$not_offered)
    needs <- option$needs
    if (identical(needs, "technical_replicates") && !isTRUE(technical)) return(table$texts[["no_technical_replicates"]])
    if (identical(needs, "no_technical_replicates") && isTRUE(technical)) return(table$texts[["technical_replicates_declared"]])
    NULL
}

#' A table as shown on this session (KIT-129 fourth pass item 2): each option that cannot be chosen carries its
#' reason (`unavailable_now`), and the selection is the default of this design: the design for technical replicates
#' when the design declares them (decision M6), else the table's (the recommended design, decision M2).
#' @param table A preset table.
#' @param technical Whether the session's design declares technical replicates.
#' @return The table.
#' @noRd
omicNormRuvMarkTable <- function(table, technical) {
    table$options <- lapply(table$options, function(option) {
        option$unavailable_now <- omicNormRuvUnavailableReason(option, table, technical)
        option
    })
    preferred <- if (isTRUE(technical)) table$selected_with_technical_replicates
    open <- vapply(table$options, \(option) is.null(option$unavailable_now) && !identical(option$status, "unavailable"), NA)
    ids <- vapply(table$options, `[[`, character(1), "id")
    table$selected <- if (!is.null(preferred) && preferred %in% ids[open]) preferred else
        if (table$selected %in% ids[open]) table$selected else table$recommended
    table
}

#' What decides whether the options are shown again: the table's name, the options that cannot be chosen and the
#' selection.
#' @noRd
omicNormRuvAvailabilityKey <- function(name, table) {
    list(name = name, closed = vapply(table$options, \(option) option$unavailable_now %||% "", character(1)),
        selected = table$selected)
}

#' Whether the head state's design declares technical replicates (the kernels' rule,
#' `resolveProtTechnicalReplicateBlock`), from the state's metadata and its design table only; FALSE without a
#' session state or when it cannot be read.
#' @param manager The bus's state manager.
#' @return One logical.
#' @noRd
omicNormRuvHeadTechnicalReplicates <- function(manager) {
    if (!inherits(manager, "DIAStateAdapter")) return(FALSE)
    tryCatch(withWorkflowStateScope(manager, function() {
        reference <- getWorkflowStateReference(manager)
        if (is.null(reference)) return(FALSE)
        spec <- readWorkflowStateMetadata(manager, reference)$data
        design <- as.data.frame(materializeWorkflowTable(manager$.context, spec$components$design_matrix))
        isTRUE(resolveProtTechnicalReplicateBlock(design, spec$slots$group_id %||% NA_character_,
            spec$slots$technical_replicate_id %||% NA_character_, check = FALSE)$technical)
    }), error = function(error) FALSE)
}

#' The panel's choice for the run (KIT-129 fourth pass item 2), read when the run starts: the design, its k rule, the
#' typed k, the seed and the supplied list, in the run path's names. NULL with the gate off or without the section's
#' input; `legacy` TRUE under a legacy design (the omic then builds today's request from today's controls).
#' Refuses, with the table's words, an option the run path does not take on this route, a rule it does not take, a
#' missing or fractional typed k, a seed that is not a whole number and an empty list; the preconditions on the data
#' (technical replicates, k against the limit and the controls) are the kernels' refusals, with their messages.
#' @param input The module's input (or a list with its values).
#' @param table The route's preset table.
#' @return NULL, or `list(ruv_design, legacy, ruv_k_rule, ruv_k, ruv_k_seed, ruv_control_ids)`.
#' @noRd
omicNormRuvDesignChoice <- function(input, table) {
    design <- input[["ruv_design"]]
    if (!omicGuiRuvDesignEnabled() || is.null(design)) return(NULL)
    refuse <- function(message) rlang::abort(message, class = "omic_norm_ruv_design_invalid")
    texts <- table$texts
    option <- Find(\(option) identical(option$id, design), table$options)
    if (is.null(option)) refuse(paste(texts[["not_available"]], design))
    reason <- if (identical(option$status, "unavailable")) option$unavailable else option$not_offered
    if (!is.null(reason)) refuse(paste(texts[["not_available"]], reason))
    if (identical(option$tier, "legacy")) return(list(ruv_design = design, legacy = TRUE))
    if (identical(design, "none")) {
        return(list(ruv_design = "none", legacy = FALSE, ruv_k_rule = NULL, ruv_k = NULL, ruv_k_seed = NULL, ruv_control_ids = NULL))
    }
    rule <- input[["ruv_k_rule"]]
    if (!checkmate::test_choice(rule, option$k_rules)) refuse(paste(texts[["k_rule_not_available"]], rule %||% ""))
    if (rule %in% names(table$k_rules_not_offered)) refuse(paste(texts[["k_rule_not_available"]], table$k_rules_not_offered[[rule]]))
    k <- NULL
    if (identical(rule, "user")) {
        k <- input[["ruv_k_user"]]
        if (!checkmate::test_number(k) || k < 1 || k != round(k)) refuse(texts[["k_missing"]])
        k <- as.integer(k)
    }
    seed <- NULL
    if (rule %in% c("num_sv_be", "user")) {
        seed <- input[["ruv_k_seed"]] %||% RUV_FACTOR_COUNT_SEED
        if (!checkmate::test_number(seed) || seed != round(seed) || abs(seed) > .Machine$integer.max) refuse(texts[["seed_invalid"]])
        seed <- as.integer(seed)
    }
    ids <- NULL
    if (identical(design, "supplied_controls")) {
        ids <- unique(strsplit(trimws(input[["ruv_control_list"]] %||% ""), "[,;[:space:]]+")[[1L]])
        ids <- ids[nzchar(ids)]
        if (!length(ids)) refuse(texts[["list_missing"]])
    }
    list(ruv_design = design, legacy = FALSE, ruv_k_rule = rule, ruv_k = k, ruv_k_seed = seed, ruv_control_ids = ids)
}

#' The route's preset table for the run (the configuration's, without building the whole configuration).
#' @param omic `proteomics`, `lipidomics`, `metabolomics` or `transcriptomics`.
#' @param route The proteomics workflow type (`DIA_LIMPA` takes its own table).
#' @noRd
omicNormRuvRunTable <- function(omic, route = NULL) {
    switch(omic,
        proteomics = if (identical(route, "DIA_LIMPA")) omicGuiRuvTableDiaLimpa() else omicGuiRuvTableClassic(),
        lipidomics = omicGuiRuvTableFeature("lipids"), metabolomics = omicGuiRuvTableFeature("metabolites"),
        transcriptomics = omicGuiRuvTableTrx())
}

#' The options column: `norm$panels`' sections in their order.
#' @noRd
omicNormOptionsColumn <- function(ns, norm) {
    sections <- omicNormOptionSections(ns, norm)
    shiny::column(3, do.call(shiny::wellPanel, c(list(shiny::h4("Normalization Options")), do.call(c, unname(sections)))))
}

# ------------------------------------------------------------------------------------------------
# The QC tabset (column 9).
# ------------------------------------------------------------------------------------------------

#' One figure tab's rows: per assay its label and the three stages' images (`<kind>_<stage>_assay<n>`;
#' mod_lipid_norm_ui_helpers.R:270-279 at f8d4c5fc).
#' @param labels The three stages' captions.
#' @param assays The assay rows.
#' @noRd
omicNormStageRows <- function(ns, kind, labels, assays, assayLabels = TRUE) {
    stages <- c("post_filter", "post_norm", "ruv_corrected")
    rows <- lapply(assays, function(n) shiny::tagList(
        if (assayLabels) shiny::fluidRow(shiny::column(12, shiny::h5(shiny::textOutput(ns(sprintf("assay%d_label_%s", n, kind))), style = "margin-bottom: 5px;"))),
        shiny::fluidRow(lapply(seq_along(stages), function(i) shiny::column(4,
            shiny::tags$small(labels[[i]], style = "display: block; text-align: center;"),
            shinyjqui::jqui_resizable(shiny::imageOutput(ns(sprintf("%s_%s_assay%d", kind, stages[[i]], n)), height = "300px")))))))
    if (length(rows) > 1L) rows <- c(rows[1L], unlist(lapply(rows[-1L], function(row) list(shiny::hr(), row)), recursive = FALSE))
    shiny::tagList(rows)
}

#' The QC tabset `norm_qc_tabs`: the tabs of `view$plot_tabs` in their order, the stage grid of `view$assay_slots` rows
#' with `view$stage_labels`, and the RUV tab's `view$ruv_heading` and `view$ruv_help` (mod_lipid_norm_ui_helpers.R:285-307
#' at f8d4c5fc). `view` is the omic's `norm` group.
#' @noRd
omicNormQcTabset <- function(ns, view) {
    assays <- seq_len(view$assay_slots)
    figure <- function(title, icon, kind) shiny::tabPanel(title, icon = shiny::icon(icon), shiny::br(),
        omicNormStageRows(ns, kind, view$stage_labels, assays, !isFALSE(view$assay_labels)))
    centred <- function(title, icon, heading, help, ...) shiny::tabPanel(title, icon = shiny::icon(icon), shiny::br(),
        shiny::fluidRow(shiny::column(12, shiny::h5(heading, style = "text-align: center;"), shiny::helpText(help), shiny::hr(), ...)))
    tabs <- list(
        itsd = centred("ITSD Selection", "flask", "Internal Standard Selection per Assay",
            "Select internal standards for each assay. Pre-selected candidates are auto-detected from naming patterns (ITSD, IS, labeled compounds).",
            shiny::uiOutput(ns("itsd_selection_ui"))),
        ruv = centred("RUV QC", "chart-line", view$ruv_heading, view$ruv_help,
            if (isTRUE(view$ruv_single_result)) omicNormRuvSingleResult(ns) else shiny::uiOutput(ns("ruv_qc_ui")),
            omicNormRuvKDiagnosticSlots(ns, view)),
        pca = figure("PCA", "project-diagram", "pca"),
        density = figure("Density", "chart-area", "density"),
        rle = figure("RLE", "chart-bar", "rle"),
        correlation = figure("Correlation", "th", "correlation"),
        # KIT-120: the filtering-progression figure and text, filled by the omic's server (proteomics).
        post_norm_qc = shiny::tabPanel("Post-Normalisation QC", icon = shiny::icon("chart-bar"), shiny::br(),
            shiny::fluidRow(shiny::column(12, shiny::h5(view$post_norm_qc$heading, style = "text-align: center;"),
                shiny::helpText(view$post_norm_qc$help),
                shinyjqui::jqui_resizable(shiny::plotOutput(ns("post_norm_filtering_summary"), height = "600px")), shiny::br(),
                shiny::verbatimTextOutput(ns("filtering_summary_text"))))),
        correlation_filter = omicNormCorrelationFilterTab(ns, view$correlation_slot),
        log = shiny::tabPanel("Log", icon = shiny::icon("list-alt"), shiny::br(), shiny::fluidRow(shiny::column(12,
            shiny::h5("Normalization Pipeline Log", style = "text-align: center;"), shiny::verbatimTextOutput(ns("norm_log"))))))
    shiny::column(9, do.call(shiny::tabsetPanel, c(list(id = ns("norm_qc_tabs")), unname(tabs[view$plot_tabs]))))
}

#' The RUV QC tab's content for an omic with one RUV result (`norm$ruv_single_result`; KIT-122, audit F5; proteomics):
#' the canonical-correlation plot, the optimisation summary and the results table as static outputs, which the omic's
#' server fills (`registerProtNormRenderOutputs`). They were created through a `renderUI` on `ruv_qc_ui`.
#' @noRd
omicNormRuvSingleResult <- function(ns) {
    shiny::fluidRow(
        shiny::column(8, shinyjqui::jqui_resizable(shiny::plotOutput(ns("ruv_canonical_correlation_plot"), height = "500px"))),
        shiny::column(4, shiny::wellPanel(
            shiny::h5("Optimization Summary"), shiny::verbatimTextOutput(ns("ruv_optimization_summary")), shiny::br(),
            shiny::h5("Optimization Results"), shiny::helpText("Detailed results for all tested percentages:"),
            DT::dataTableOutput(ns("ruv_optimization_table")))))
}

#' The correlation filtering tab, the normalisation's final step (mod_lipid_norm_ui_helpers.R:311-371 at f8d4c5fc).
#' @noRd
omicNormCorrelationFilterTab <- function(ns, slot = NULL) {
    shiny::tabPanel("Correlation Filtering", icon = shiny::icon("filter"), shiny::br(), shiny::fluidRow(
        shiny::column(4, shiny::wellPanel(
            shiny::h5("Sample Correlation Filtering"),
            shiny::helpText("FINAL STEP: Filter samples based on their correlation with other samples in the same group. Low correlation indicates potential technical issues."),
            shiny::sliderInput(ns("min_pearson_correlation_threshold"), "Min. Pearson Correlation Threshold:", min = 0.0, max = 1.0,
                value = 0.5, step = 0.05, width = "100%"),
            shiny::helpText("Samples with correlation below this threshold will be flagged"),
            shiny::br(),
            shiny::actionButton(ns("apply_correlation_filter"), "Apply Correlation Filter", class = "btn-primary", width = "100%",
                icon = shiny::icon("play")),
            shiny::br(), shiny::br(),
            shiny::actionButton(ns("skip_correlation_filter"), "Skip Filtering & Proceed to DE", class = "btn-info", width = "100%",
                icon = shiny::icon("forward")),
            shiny::helpText("Bypass sample filtering and proceed directly to DE"),
            shiny::br(), shiny::br(),
            shiny::h5("Filter Results"),
            shiny::verbatimTextOutput(ns("correlation_filter_summary")),
            # KIT-120: an omic's own panel under the summary (`norm$correlation_slot`; proteomics' DIA_LIMPA covariates).
            if (!is.null(slot)) shiny::uiOutput(ns(slot)))),
        shiny::column(8,
            shiny::h5("Final Filtered Data QC", style = "text-align: center;"),
            shiny::helpText("QC plots for the final dataset after correlation filtering"),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("final_qc_plot"), height = "500px")))))
}

# ------------------------------------------------------------------------------------------------
# The tab's state, log and aesthetics.
# ------------------------------------------------------------------------------------------------

#' The tab's reactive values (mod_lipid_norm_runtime_helpers.R:1-23 at f8d4c5fc). KIT-075: the `*_obj` slots hold state
#' references, not objects.
#' @noRd
omicNormReactiveState <- function() {
    shiny::reactiveValues(assay_names = NULL, itsd_selections = list(), pre_norm_qc_generated = FALSE,
        normalization_complete = FALSE, ruv_complete = FALSE, correlation_filtering_complete = FALSE, qc_plot_paths = list(),
        post_filter_obj = NULL, post_itsd_obj = NULL, post_log2_obj = NULL, post_norm_obj = NULL, ruv_corrected_obj = NULL,
        correlation_filtered_obj = NULL, ruv_optimization_results = list(), correlation_results = list(),
        plot_refresh_trigger = 0, normalization_log = character(0))
}

#' The tab's log: `addLog(message)` appends `[HH:MM:SS] message` to `normData$normalization_log`.
#' @noRd
omicNormAddLog <- function(normData) {
    function(message) {
        normData$normalization_log <- c(normData$normalization_log, sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), message))
        invisible(normData$normalization_log)
    }
}

#' The plot aesthetics: the colour and shape variables, `group` when unset.
#' @noRd
omicNormPlotAestheticsGetter <- function(input) {
    function() {
        list(color_var = if (is.null(input$color_variable) || identical(input$color_variable, "")) "group" else input$color_variable,
            shape_var = if (is.null(input$shape_variable) || identical(input$shape_variable, "")) "group" else input$shape_variable)
    }
}

# ------------------------------------------------------------------------------------------------
# The stage grid's renderers and registrars (with an omic's file stem and stage prefixes: trx's, KIT-110).
# ------------------------------------------------------------------------------------------------

#' An image slot's renderer factory: `function(assaySlot, plotType, stagePrefix)` gives the renderImage of
#' `<sanitised assay>_<stagePrefix>_<plotType>.png` in `experimentPaths[[qcDirKey]]`.
#' @param sanitizeAssayNameFn The file stem of an assay name (trx: `trx`, its `norm$image_stem`).
#' @param qcDirKey The QC directory's key in `experimentPaths`.
#' @noRd
omicNormQcImageRenderer <- function(normData, experimentPaths,
    sanitizeAssayNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)), qcDirKey) {
    function(assaySlot, plotType, stagePrefix) {
        # Bind this output's own spec now: the registrars call this in a loop (F105-01).
        force(assaySlot)
        force(plotType)
        force(stagePrefix)
        shiny::renderImage({
            normData$plot_refresh_trigger
            assay_names <- normData$assay_names
            if (is.null(assay_names) || length(assay_names) < assaySlot) return(list(src = "", alt = "Assay not detected yet"))
            assay_name <- assay_names[[assaySlot]]
            # KIT-120: an empty stem gives `<stagePrefix>_<plotType>.png` (`norm$image_stemless`).
            stem <- sanitizeAssayNameFn(assay_name)
            filename <- paste0(if (nzchar(stem)) paste0(stem, "_"), stagePrefix, "_", plotType, ".png")
            qc_dir <- experimentPaths[[qcDirKey]]
            if (is.null(qc_dir)) return(list(src = "", alt = "QC directory not configured"))
            img_path <- file.path(qc_dir, filename)
            if (!file.exists(img_path)) return(list(src = "", alt = paste("Plot not generated yet:", filename)))
            list(src = img_path, contentType = "image/png", width = "100%", height = "auto", alt = paste(plotType, "-", assay_name))
        }, deleteFile = FALSE)
    }
}

#' The stage grid's image slots `<kind>_<post_filter|post_norm|ruv_corrected>_assay<n>`, for `assaySlots` assays, each the
#' renderer of its assay, kind and stage file prefix.
#' @param stagePrefixes The file prefix of each stage (KIT-097: trx's `raw`, `tmm`, `ruvs`, its `norm$stage_prefixes`).
#' @noRd
registerOmicNormQcImageOutputs <- function(output, renderQcImageForAssay,
    stagePrefixes = c(pre_norm = "pre_norm", post_norm = "post_norm", ruv_corrected = "ruv_corrected"), assaySlots = 2L) {
    stages <- c(post_filter = "pre_norm", post_norm = "post_norm", ruv_corrected = "ruv_corrected")
    for (kind in c("pca", "density", "rle", "correlation")) {
        for (n in seq_len(assaySlots)) {
            for (stage in names(stages)) {
                output[[sprintf("%s_%s_assay%d", kind, stage, n)]] <- renderQcImageForAssay(n, kind, stagePrefixes[[stages[[stage]]]])
            }
        }
    }
    invisible(output)
}

#' An assay label's renderer factory: `function(assaySlot)` gives `Assay: <name>`, or `Assay <n>: (detecting...)`.
#' @noRd
omicNormAssayLabelRenderer <- function(normData) {
    function(assaySlot) {
        # Bind this output's own slot now: the registrars call this in a loop (F105-01).
        force(assaySlot)
        shiny::renderText({
            if (!is.null(normData$assay_names) && length(normData$assay_names) >= assaySlot) {
                paste("Assay:", normData$assay_names[[assaySlot]])
            } else {
                paste0("Assay ", assaySlot, ": (detecting...)")
            }
        })
    }
}

#' The stage grid's labels `assay<n>_label_<kind>`.
#' @noRd
registerOmicNormAssayLabelOutputs <- function(output, renderAssayLabel, assaySlots = 2L) {
    for (kind in c("pca", "density", "rle", "correlation")) {
        for (n in seq_len(assaySlots)) output[[sprintf("assay%d_label_%s", n, kind)]] <- renderAssayLabel(n)
    }
    invisible(output)
}

#' The log tab's renderer factory.
#' @noRd
omicNormLogRenderer <- function(normData) {
    function() {
        shiny::renderText({
            if (length(normData$normalization_log) == 0) return("Normalization log will appear here as you apply steps...")
            paste(normData$normalization_log, collapse = "\n")
        })
    }
}

#' @noRd
registerOmicNormLogOutput <- function(output, renderNormLog) {
    output$norm_log <- renderNormLog()
    invisible(output)
}

#' The RUV QC tab's renderer factory: per assay its canonical-correlation plot, summary and results table.
#' @noRd
omicNormRuvQcUiRenderer <- function(normData, ns) {
    function() {
        shiny::renderUI({
            shiny::req(normData$assay_names)
            # KIT-129 fourth pass item 5: after a new design's run (set under the gate only) there is no sweep or grid
            # to show; the RUV design section's outputs show what the run did.
            result <- normData$ruv_design_result
            if (!is.null(result) && !isTRUE(result$legacy)) return(NULL)
            shiny::tagList(purrr::map(normData$assay_names, function(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                shiny::tagList(
                    shiny::fluidRow(shiny::column(12, shiny::h5(paste("Assay:", assay_name),
                        style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"))),
                    shiny::fluidRow(
                        shiny::column(8, shinyjqui::jqui_resizable(shiny::plotOutput(ns(paste0("cancor_plot_", safe_name)), height = "400px"))),
                        shiny::column(4, shiny::wellPanel(shiny::h6("Optimization Summary"),
                            shiny::verbatimTextOutput(ns(paste0("ruv_summary_", safe_name))), shiny::br(),
                            shiny::h6("Results Table"), DT::dataTableOutput(ns(paste0("ruv_table_", safe_name)))))),
                    shiny::hr())
            }))
        })
    }
}

#' @noRd
registerOmicNormRuvQcOutput <- function(output, renderRuvQcUi) {
    output$ruv_qc_ui <- renderRuvQcUi()
    invisible(output)
}

#' Each assay's RUV outputs once the optimisation has results: the canonical-correlation plot, the summary and the table
#' (mod_lipid_norm_observer_helpers.R:87-148 at f8d4c5fc).
#' @noRd
registerOmicNormRuvCancorOutputs <- function(output, normData) {
    shiny::observe({
        shiny::req(normData$assay_names)
        shiny::req(length(normData$ruv_optimization_results) > 0)
        purrr::walk(normData$assay_names, function(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            output[[paste0("cancor_plot_", safe_name)]] <- shiny::renderPlot({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success) && !is.null(result$cancor_plot)) return(result$cancor_plot)
                omicNormTextPlot("No data", 6)
            })
            output[[paste0("ruv_summary_", safe_name)]] <- shiny::renderText(omicNormRuvSummary(normData$ruv_optimization_results[[assay_name]]))
            output[[paste0("ruv_table_", safe_name)]] <- DT::renderDataTable({
                result <- normData$ruv_optimization_results[[assay_name]]
                if (!is.null(result) && isTRUE(result$success) && !is.null(result$optimization_results)) {
                    DT::datatable(result$optimization_results, options = list(pageLength = 5, dom = "t"), rownames = FALSE)
                } else {
                    NULL
                }
            })
        })
    })
    invisible(output)
}

#' One assay's RUV optimisation summary text.
#' @noRd
omicNormRuvSummary <- function(result) {
    if (is.null(result)) return("Not yet computed")
    if (!isTRUE(result$success)) return(paste0("Failed: ", result$error))
    paste0("Best k: ", result$best_k, "\n", "Best %: ", result$best_percentage, "\n", "Separation: ",
        round(result$separation_score, 4), "\n", "Controls: ", sum(result$control_genes_index, na.rm = TRUE))
}

#' A text placeholder plot.
#' @noRd
omicNormTextPlot <- function(label, size) {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = size) + ggplot2::theme_void()
}

# ------------------------------------------------------------------------------------------------
# The ITSD selection, the correlation summary and the final QC plot.
# ------------------------------------------------------------------------------------------------

#' The ITSD selection tab's UI: one table per assay.
#' @noRd
omicNormItsdSelectionUi <- function(normData, ns) {
    shiny::renderUI({
        shiny::req(normData$assay_names)
        shiny::tagList(purrr::map(normData$assay_names, function(assay_name) {
            safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
            shiny::wellPanel(shiny::h5(paste("Assay:", assay_name)), DT::dataTableOutput(ns(paste0("itsd_table_", safe_name))), shiny::br())
        }))
    })
}

#' The ITSD selection tables, built from the head's state a page at a time (KIT-082; nothing whole kept), and the shown
#' rows' feature ids (`itsd_candidate_ids`, KIT-078: a managed run resolves a selection without the matrix). Synchronous,
#' so that Run cannot race the selection (mod_lipid_norm_observer_helpers.R:150-230 at f8d4c5fc).
#' @param config The omic's GUI configuration (its binding and `norm$itsd_table`).
#' @noRd
omicNormItsdTables <- function(output, workflowData, normData, config) {
    shiny::req(normData$assay_names)
    shiny::req(workflowData$state_manager)
    binding <- config$binding$binding
    shell <- tryCatch(readWorkflowOmicsShell(binding, workflowData$state_manager), error = function(e) NULL)
    if (is.null(shell) || !inherits(shell, binding$class)) return(invisible(NULL))
    id_col <- methods::slot(shell, binding$id_slot)
    annotation_cols <- shell@annotation_id_column
    design_samples <- shell@design_matrix[[shell@sample_id]]
    rm(shell)
    build <- config$norm$itsd_table
    # KIT-214 (E-8): each candidate's mean and CV over the design's samples that hold numbers (the builder's own rule only
    # when there are none).
    build_page <- function(page) {
        sample_columns <- omicsDesignNumericColumnsOf(page, design_samples, id_col)
        if (length(sample_columns)) build(page, id_col, annotation_cols, sample_columns = sample_columns) else {
            build(page, id_col, annotation_cols)
        }
    }
    stored_tables <- tryCatch(readWorkflowOmicsItsdTables(binding, workflowData$state_manager, normData$assay_names,
        build_page), error = function(e) {
        logger::log_warn(paste("Could not read the ITSD candidate tables:", conditionMessage(e)))
        NULL
    })
    if (is.null(stored_tables)) return(invisible(NULL))
    selection_tables <- lapply(stats::setNames(normData$assay_names, normData$assay_names), function(assay_name) stored_tables[[assay_name]])
    normData$itsd_candidate_ids <- lapply(selection_tables, function(table) table$feature_id)
    purrr::walk(normData$assay_names, function(assay_name) {
        selection_table <- selection_tables[[assay_name]]
        output[[paste0("itsd_table_", gsub("[^A-Za-z0-9]", "_", tolower(assay_name)))]] <- DT::renderDataTable({
            if (is.null(selection_table)) return(NULL)
            DT::datatable(selection_table, selection = list(mode = "multiple", selected = which(selection_table$is_candidate)),
                filter = "top", options = list(pageLength = 10, scrollX = TRUE, order = list(list(4, "desc"), list(3, "asc"))),
                rownames = FALSE) |>
                DT::formatStyle("is_candidate", backgroundColor = DT::styleEqual(TRUE, "#d4edda")) |>
                DT::formatRound(columns = c("mean_intensity", "cv_percent"), digits = 2)
        })
    })
    invisible(NULL)
}

#' The ITSD selection's tracking: each assay table's selected rows into `normData$itsd_selections`.
#' @noRd
registerOmicNormItsdSelectionTracking <- function(input, normData) {
    shiny::observe({
        shiny::req(normData$assay_names)
        purrr::walk(normData$assay_names, function(assay_name) {
            input_id <- paste0("itsd_table_", gsub("[^A-Za-z0-9]", "_", tolower(assay_name)), "_rows_selected")
            shiny::observeEvent(input[[input_id]], {
                normData$itsd_selections[[assay_name]] <- input[[input_id]]
                logger::log_info(paste("ITSD selection updated for", assay_name, ":", length(input[[input_id]]), "features selected"))
            }, ignoreNULL = FALSE)
        })
    })
    invisible(normData)
}

#' The correlation filter's summary: per assay the pairs and their correlations, and the samples kept (the designs from
#' the states' metadata, KIT-075).
#' @param readMetadata The omic's metadata reader (`config$binding$read_metadata`).
#' @noRd
omicNormCorrelationSummary <- function(normData, workflowData, readMetadata) {
    shiny::renderText({
        if (!normData$correlation_filtering_complete) return("Apply correlation filter to see results...")
        corr_results <- normData$correlation_results
        filtered_obj <- normData$correlation_filtered_obj
        original_obj <- if (!is.null(normData$ruv_corrected_obj)) normData$ruv_corrected_obj else normData$post_norm_obj
        if (is.null(corr_results) || length(corr_results) == 0) return("No correlation results available.")
        summary_lines <- c("=== Correlation Filtering Summary ===\n")
        for (assay_name in names(corr_results)) {
            assay_corr <- corr_results[[assay_name]]
            if (!is.null(assay_corr) && nrow(assay_corr) > 0) {
                values <- assay_corr$pearson_correlation
                summary_lines <- c(summary_lines, sprintf("\n[%s]\n  Sample pairs: %d\n  Correlation: mean=%.3f, min=%.3f, max=%.3f",
                    assay_name, nrow(assay_corr), round(mean(values, na.rm = TRUE), 3), round(min(values, na.rm = TRUE), 3),
                    round(max(values, na.rm = TRUE), 3)))
            }
        }
        if (!is.null(original_obj) && !is.null(filtered_obj)) {
            original_samples <- nrow(readMetadata(workflowData$state_manager, original_obj)$design_matrix)
            filtered_samples <- nrow(readMetadata(workflowData$state_manager, filtered_obj)$design_matrix)
            summary_lines <- c(summary_lines, sprintf(
                "\n\n[Sample Filtering]\n  Original: %d samples\n  After filtering: %d samples\n  Removed: %d samples",
                original_samples, filtered_samples, original_samples - filtered_samples))
        }
        paste(summary_lines, collapse = "")
    })
}

#' The final QC plot: the PCA of the correlation-filtered state (else the RUV-corrected, else the normalised). KIT-082: a
#' managed session's worker draws it (`workflowOmicsFinalPca`), with the omic's placeholder sizes (D7); in process the
#' state is read for this render only (KIT-075) and the binding's `final_pca` draws it.
#' @param getPlotAesthetics The tab's aesthetics getter.
#' @param config The omic's GUI configuration.
#' @noRd
omicNormFinalQcPlot <- function(normData, getPlotAesthetics, workflowData, config) {
    binding <- config$binding$binding
    sizes <- config$norm$final_pca
    shiny::renderPlot({
        shiny::req(normData$correlation_filtering_complete || normData$ruv_complete)
        current_reference <- normData$correlation_filtered_obj %||% normData$ruv_corrected_obj %||% normData$post_norm_obj
        if (is.null(current_reference)) return(omicNormTextPlot("No data", 6))
        aesthetics <- getPlotAesthetics()
        runner <- workflowOmicsManagedActions(binding$omic_type)
        if (!is.null(runner)) {
            entry <- workflowOmicsFinalPca(binding, runner, workflowData$state_manager, current_reference, aesthetics$color_var,
                aesthetics$shape_var)
            if (is.null(entry)) return(omicNormTextPlot("Preparing the final QC plot...", sizes$preparing_size))
            if (!is.null(entry$message)) return(omicNormTextPlot(paste("Error:", entry$message), sizes$error_size))
            return(grid::grid.draw(entry$plot))
        }
        current_s4 <- config$binding$read_state(workflowData$state_manager, current_reference)
        tryCatch(binding$final_pca(current_s4, aesthetics$color_var, aesthetics$shape_var),
            error = function(e) omicNormTextPlot(paste("Error:", e$message), 4))
    })
}
