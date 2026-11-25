# KIT-129 W2: the RUV design presets of the configuration (`norm$ruv_design`), their schema, and the shared panel
# section of the normalisation tab. Presentation and configuration only: no option is wired to a computation, and the
# section is not shown until the gate (`options(MultiScholaR.ruv_design_panel = TRUE)`) is on. The expected values
# below are written from the ticket (sections 1a to 1d and 2) and the design review of 2026-11-18 (sections 3 to 5),
# independently of the implementation's table.

ruvDesignOmics <- c("proteomics", "lipidomics", "metabolomics", "transcriptomics")
ruvDesignIds <- c("groups_all_features", "technical_replicates", "supplied_controls", "first_pass_nulls", "legacy", "none")

# Each table's measured figures: false discovery proportion in the worst cell, and power (ticket 1a to 1d; since
# decision M15 of 2026-11-18 the first-pass nulls' is the one measured with the estimated number). Since the KIT-129
# final figures (2026-11-20, decisions M40, M41) every figure is the shipped method's: KIT-124's confirmation round
# (scored/confirm-gha), or the simulator arm the round found identical (docs/qa/workflow-execution/kit/kit-129/
# final-figures.tsv).
ruvDesignFigures <- list(
    classic = list(groups_all_features = c(0.067, 0.92), technical_replicates = c(0.055, 0.90),
        # KIT-201 (VAL-131, D2): the Legacy figure is KIT-124 S1 arm 28's, the legacy test with TREAT's trend and robust
        # settings (scored/exact-gha: worst cell C5 0.389, C0 0.101); it was arm 3's 0.385.
        supplied_controls = c(0.054, 0.91), first_pass_nulls = c(0.079, 0.92), legacy = c(0.389, 0.92),
        none = c(0.052, 0.49)),
    dia_limpa = list(groups_all_features = c(0.056, 0.91), technical_replicates = c(0.056, 0.77),
        supplied_controls = c(0.056, 0.91), first_pass_nulls = c(0.147, 0.92), legacy = c(0.435, 0.90),
        legacy_correction = c(0.466, 0.91), none = c(0.191, 0.56)),
    transcriptomics = list(groups_all_features = c(0.055, 0.77), supplied_controls = c(0.056, 0.76),
        first_pass_nulls = c(0.071, 0.77), legacy = c(0.299, 0.73), none = c(0.044, 0.005)),
    feature = list(groups_all_features = c(0.085, 0.47), technical_replicates = c(0.063, 0.34),
        supplied_controls = c(0.079, 0.40), first_pass_nulls = c(0.203, 0.48), legacy = c(0.557, 0.42),
        none = c(0.059, 0.003)))

# The figures measured with the separation rule, where it is not the pre-selected rule (the same definition).
ruvDesignSeparationFigures <- list(
    classic = list(supplied_controls = c(0.060, 0.91), first_pass_nulls = c(0.079, 0.91)),
    dia_limpa = list(supplied_controls = c(0.057, 0.90), first_pass_nulls = c(0.151, 0.92)),
    transcriptomics = list(supplied_controls = c(0.053, 0.61), first_pass_nulls = c(0.064, 0.60)),
    feature = list(supplied_controls = c(0.095, 0.26), first_pass_nulls = c(0.261, 0.53)))

# The k rules each design offers and the one pre-selected (ticket section 2, the user's direction U1).
ruvDesignKRules <- list(
    groups_all_features = list(rules = c("num_sv_be", "user"), default = "num_sv_be", pending = NULL),
    technical_replicates = list(rules = c("num_sv_be", "user"), default = "num_sv_be", pending = NULL),
    supplied_controls = list(rules = c("num_sv_be", "sweep_bounded", "user"), default = "num_sv_be", pending = NULL),
    first_pass_nulls = list(rules = c("num_sv_be", "sweep_bounded", "user"), default = "num_sv_be", pending = NULL),
    legacy = list(rules = character(), default = NULL, pending = NULL),
    legacy_correction = list(rules = character(), default = NULL, pending = NULL),
    none = list(rules = character(), default = NULL, pending = NULL))

#' Every table of the four omics, named `<omic>.<table>`.
ruvDesignTables <- function() {
    tables <- lapply(stats::setNames(ruvDesignOmics, ruvDesignOmics), \(omic) omicGuiConfig(omic)$norm$ruv_design$tables)
    unlist(tables, recursive = FALSE)
}

#' The figures expected of one table.
ruvDesignExpectedFigures <- function(name, figures = ruvDesignFigures) {
    key <- sub("^(lipidomics|metabolomics)\\..*$", "feature", name)
    key <- sub("^transcriptomics\\..*$", "transcriptomics", key)
    figures[[sub("^proteomics\\.", "", key)]]
}

ruvDesignOptions <- function(table) stats::setNames(table$options, vapply(table$options, `[[`, "", "id"))

#' The norm UI's HTML with the gate on or off.
ruvDesignHtml <- function(omic, enabled) {
    withr::local_options(MultiScholaR.ruv_design_panel = enabled)
    as.character(htmltools::renderTags(mod_omic_norm_ui("norm", omicGuiConfig(omic)$norm))$html)
}

ruvDesignMaskTabs <- function(html) {
    gsub("tab-[0-9]+-([0-9]+)", "tab-ID-\\1", gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"ID\"", html))
}

expectRuvDesignRefused <- function(config, field) {
    expect_error(validateOmicGuiConfig(config), regexp = field, fixed = TRUE, class = "omic_gui_config_invalid")
}

# ------------------------------------------------------------------------------------------------
# The table: completeness per omic, the figures, the pending cells.
# ------------------------------------------------------------------------------------------------

test_that("every omic carries the preset table: every design once, in four tiers, one recommended", {
    tables <- ruvDesignTables()
    expect_setequal(names(tables), c("proteomics.classic", "proteomics.dia_limpa", "lipidomics.features",
        "metabolomics.features", "transcriptomics.genes"))
    expect_identical(omicGuiConfig("proteomics")$norm$ruv_design$route_tables, c(DIA_LIMPA = "dia_limpa"))
    for (name in names(tables)) {
        options <- ruvDesignOptions(tables[[name]])
        expected <- if (identical(name, "proteomics.dia_limpa")) append(ruvDesignIds, "legacy_correction", 5L) else ruvDesignIds
        expect_identical(names(options), expected, info = name)
        tiers <- vapply(options, `[[`, "", "tier")
        expect_identical(unname(tiers[ruvDesignIds]), c("recommended", "valid_with_conditions", "valid_with_conditions",
            "valid_with_conditions", "legacy", "none"), info = name)
        # KIT-129 fourth pass item 2 (decision M2): the recommended design is the selected one (the section is shown only
        # with its gate on).
        expect_identical(tables[[name]]$selected, "groups_all_features", info = name)
        expect_identical(tables[[name]]$recommended, "groups_all_features", info = name)
        expect_identical(tables[[name]]$selected_with_technical_replicates,
            if (startsWith(name, "transcriptomics")) NULL else "technical_replicates", info = name)
        for (option in options) {
            for (field in c("label", "description", "tooltip", "note", "source")) {
                expect_true(nzchar(option[[field]]), info = paste(name, option$id, field))
            }
        }
    }
})

test_that("each measured figure is the ticket's; the cells still running are pending, not numbers", {
    tables <- ruvDesignTables()
    for (name in names(tables)) {
        options <- ruvDesignOptions(tables[[name]])
        expected <- ruvDesignExpectedFigures(name)
        measured <- names(Filter(\(option) identical(option$status, "measured"), options))
        expect_setequal(measured, names(expected))
        for (id in names(expected)) {
            expect_identical(c(options[[id]]$fdp, options[[id]]$power), expected[[id]], info = paste(name, id))
            expect_true(nzchar(options[[id]]$fdp_where), info = paste(name, id))
        }
        for (id in setdiff(names(options), measured)) {
            expect_false(any(c("fdp", "power") %in% names(options[[id]])), info = paste(name, id))
        }
        # The separation rule's own figure where it is offered and not pre-selected, and nowhere else.
        separation <- ruvDesignExpectedFigures(name, ruvDesignSeparationFigures)
        with_other <- names(Filter(\(option) length(option$k_rule_figures) > 0L, options))
        expect_setequal(with_other, names(separation))
        for (id in names(separation)) {
            figure <- options[[id]]$k_rule_figures
            expect_identical(names(figure), "sweep_bounded", info = paste(name, id))
            expect_identical(c(figure$sweep_bounded$fdp, figure$sweep_bounded$power), separation[[id]], info = paste(name, id))
            expect_true(all(nzchar(unlist(figure$sweep_bounded[c("fdp_where", "source")]))), info = paste(name, id))
        }
    }
    # A figure measured in one scenario says so: the technical-replicate designs, and no other.
    single <- unlist(lapply(tables, \(table) vapply(table$options, \(option) isTRUE(option$single_scenario), NA)))
    ids <- unlist(lapply(tables, \(table) vapply(table$options, `[[`, "", "id")))
    expect_identical(unname(ids[single]), rep("technical_replicates", 4L))
    trx <- ruvDesignOptions(tables[["transcriptomics.genes"]])$technical_replicates
    expect_identical(trx$status, "unavailable")
    expect_match(trx$unavailable, "added together into one set of counts before the analysis", fixed = TRUE)
    statuses <- unlist(lapply(tables, \(table) vapply(table$options, `[[`, "", "status")))
    expect_identical(sum(statuses == "pending"), 0L)
    expect_identical(sum(statuses == "unavailable"), 1L)
})

test_that("the wording rulings of 2026-11-19 are the shipped notes and labels (R1 to R3)", {
    tables <- ruvDesignTables()
    # R3: the caution leads the note of controls from the data, with the omic's own figure under association.
    # KIT-129 final figures (2026-11-20): the shipped method's figure at a correlation of 0.3; since decision M42 also
    # its figure at 0.6, measured with the separation rule (scored/confirm06-gha).
    at06 <- "at a correlation of 0.6, with the separation rule for the number of factors)."
    caution <- c(proteomics.classic = paste("reported more false discoveries than the 5% target (7.9% at a correlation of 0.3,",
            "with the estimated number of factors; 28.0%", at06),
        proteomics.dia_limpa = paste("far more false discoveries than the 5% target (14.7% at a correlation of 0.3, with the",
            "estimated number of factors; 30.2%", at06),
        lipidomics.features = paste("far more false discoveries than the 5% target (20.3% at a correlation of 0.3, with the",
            "estimated number of factors; 56.2%", at06),
        metabolomics.features = paste("far more false discoveries than the 5% target (20.3% at a correlation of 0.3, with the",
            "estimated number of factors; 56.2%", at06),
        transcriptomics.genes = paste("more false discoveries than the 5% target (6.4% at a correlation of 0.3 and 9.6%", at06))
    # The last sentence points at the check's verdict only where the check has a detection limit for this approach
    # (classic proteomics today); elsewhere it says the check gives none.
    verdict <- "If the check below reports a detectable association, choose the recommended approach."
    no_verdict <- paste("The check below shows how strongly the estimated unwanted factors are associated with your groups,",
        "but gives no verdict for this approach. Unless you know the unwanted variation does not differ between your groups,",
        "choose the recommended approach.")
    for (name in names(caution)) {
        note <- ruvDesignOptions(tables[[name]])$first_pass_nulls$note
        expect_true(startsWith(note, "Use this approach only when the unwanted variation is unlikely to differ"), info = name)
        expect_match(note, caution[[name]], fixed = TRUE, info = name)
        with_limit <- identical(name, "proteomics.classic")
        expect_identical(grepl(verdict, note, fixed = TRUE), with_limit, info = name)
        expect_identical(grepl(no_verdict, note, fixed = TRUE), !with_limit, info = name)
    }
    # The sentence follows the thresholds table, not a list of omics.
    expect_match(omicGuiRuvFirstPassCaution("proteomics", "more", "x"), verdict, fixed = TRUE)
    expect_match(omicGuiRuvFirstPassCaution("dia_limpa", "more", "x"), no_verdict, fixed = TRUE)
    testthat::local_mocked_bindings(ruvGroupAssociationThresholds = function() data.frame(omic = "lipidomics",
        ruv_design = "first_pass_nulls"))
    expect_match(omicGuiRuvFirstPassCaution("lipidomics", "more", "x"), verdict, fixed = TRUE)
    expect_match(omicGuiRuvFirstPassCaution("proteomics", "more", "x"), no_verdict, fixed = TRUE)
    # M42 (the user, 2026-11-21): DIA with limpa's recommended note states its result at a correlation of 0.6.
    expect_identical(ruvDesignOptions(tables[["proteomics.dia_limpa"]])$groups_all_features$note, paste("The general-purpose",
        "choice. In every scenario tested up to a correlation of 0.3 with the groups it did as well as an analysis that was",
        "given the true unwanted factor (5.3% false discoveries with unwanted variation tied to the groups, correlation 0.3).",
        "With unwanted variation tied to the groups at a correlation of 0.6 it had 6.6% false discoveries, against 4.7% for",
        "an analysis given the true unwanted factor, and 9.4% of simulated data sets had more than 10% false discoveries,",
        "against 0.6%. Measured with the fastlimpadpc engine; expected to hold for the limpa engine, not measured there."))
    # R1: transcriptomics quotes the separation rule that runs (power 25% tied to the groups, 5% with three groups).
    expect_match(ruvDesignOptions(tables[["transcriptomics.genes"]])$first_pass_nulls$note,
        "tied to the groups (power 25%) or with three groups (5%).", fixed = TRUE)
    # R2: no default list promised for lipids and metabolites.
    for (name in c("lipidomics.features", "metabolomics.features")) {
        features <- tables[[name]]$features
        supplied <- ruvDesignOptions(tables[[name]])$supplied_controls
        expect_identical(supplied$label, sprintf(paste("Your own list of control %s, for example internal standards that",
            "were not used for normalisation"), features))
        expect_false(grepl("by default|unless you supply", paste(supplied$label, supplied$description)), info = features)
    }
})

test_that("each design offers the k rules defined for it, with the rule its figure was measured with selected", {
    design <- omicGuiConfig("lipidomics")$norm$ruv_design
    expect_identical(names(design$k_rules), c("num_sv_be", "sweep_bounded", "user"))
    for (rule in design$k_rules) expect_true(all(nzchar(unlist(rule[c("label", "reads", "pros", "cons", "diagnostic")]))))
    tables <- ruvDesignTables()
    for (name in names(tables)) {
        for (option in tables[[name]]$options) {
            expected <- ruvDesignKRules[[option$id]]
            expect_identical(option$k_rules, expected$rules, info = paste(name, option$id))
            # (`[[`: `$k_rule` would partially match `k_rules`)
            expect_identical(option[["k_rule"]], expected$default, info = paste(name, option$id))
            expect_identical(option$k_rules_pending, expected$pending, info = paste(name, option$id))
        }
    }
})

test_that("the texts shown name no internal arm, cell, ticket or question", {
    internal <- "\\b(arm|arms|KIT-[0-9]+|C[0-9]|S[0-9]{1,2}|U[0-9]{1,2}|W[0-9]{1,2}|R0|B[0-9])\\b"
    for (omic in ruvDesignOmics) {
        design <- omicGuiConfig(omic)$norm$ruv_design
        shown <- c(design$tiers, design$conditions, unlist(design$k_rules), unlist(lapply(design$tables, \(table) c(
            unlist(lapply(table$options, \(option) c(option[c("label", "description", "tooltip", "note", "fdp_where",
                "unavailable")], omicNormRuvAltText(option))))))))
        expect_identical(grep(internal, unname(shown), value = TRUE), character(), info = omic)
    }
})

# ------------------------------------------------------------------------------------------------
# The schema.
# ------------------------------------------------------------------------------------------------

test_that("the schema refuses an incomplete or inconsistent preset table by class, naming the field", {
    config <- omicGuiConfig("proteomics")
    path <- "norm$ruv_design"
    change <- function(edit) {
        bad <- config
        bad$norm$ruv_design <- edit(bad$norm$ruv_design)
        bad
    }
    for (field in c("enabled", "testid", "tiers", "k_rules", "tables")) {
        expectRuvDesignRefused(change(\(d) d[setdiff(names(d), field)]), paste0(path, "$", field))
    }
    expectRuvDesignRefused(change(\(d) utils::modifyList(d, list(enabled = "yes"))), paste0(path, "$enabled"))
    expectRuvDesignRefused(change(\(d) {
        d$tiers <- d$tiers[-1L]
        d
    }), paste0(path, "$tiers"))
    expectRuvDesignRefused(change(\(d) {
        d$k_rules$sweep_bounded$cons <- NULL
        d
    }), paste0(path, "$k_rules$sweep_bounded$cons"))
    expectRuvDesignRefused(change(\(d) {
        d$route_tables <- c(DIA_LIMPA = "absent")
        d
    }), paste0(path, "$route_tables"))
    # One table: every design once; the selected and the recommended among them.
    table_path <- paste0(path, "$tables$classic")
    table <- function(edit) change(\(d) {
        d$tables$classic <- edit(d$tables$classic)
        d
    })
    expectRuvDesignRefused(table(\(t) {
        t$options <- t$options[-2L]
        t
    }), paste0(table_path, "$options"))
    expectRuvDesignRefused(table(\(t) {
        t$options <- c(t$options, t$options[1L])
        t
    }), paste0(table_path, "$options"))
    expectRuvDesignRefused(table(\(t) utils::modifyList(t, list(selected = "automatic"))), paste0(table_path, "$selected"))
    expectRuvDesignRefused(table(\(t) utils::modifyList(t, list(recommended = "legacy"))), paste0(table_path, "$recommended"))
    expectRuvDesignRefused(table(\(t) t[setdiff(names(t), "selected")]), paste0(table_path, "$selected"))
    # One option (the first: the recommended design, measured, all features as controls).
    option_path <- paste0(table_path, "$options[[1]]")
    option <- function(edit) table(\(t) {
        t$options[[1L]] <- edit(t$options[[1L]])
        t
    })
    for (field in c("id", "tier", "label", "description", "tooltip", "note", "status", "source", "test_form", "k_rules")) {
        expectRuvDesignRefused(option(\(o) o[setdiff(names(o), field)]), paste0(option_path, "$", field))
    }
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(tier = "preferred"))), paste0(option_path, "$tier"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(status = "estimated"))), paste0(option_path, "$status"))
    # A measured option carries its figure; a pending or unavailable one carries no number.
    expectRuvDesignRefused(option(\(o) o[setdiff(names(o), "fdp")]), paste0(option_path, "$fdp"))
    expectRuvDesignRefused(option(\(o) o[setdiff(names(o), "power")]), paste0(option_path, "$power"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(fdp = 6.6))), paste0(option_path, "$fdp"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(status = "pending"))), paste0(option_path, "$fdp"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o[setdiff(names(o), c("fdp", "fdp_where", "power"))],
        list(status = "unavailable"))), paste0(option_path, "$unavailable"))
    # The k rules: known ones, the selected among them, none pending that is not offered, and no separation rule
    # where all features are the controls.
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(k_rules = c("num_sv_be", "elbow")))),
        paste0(option_path, "$k_rules"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(k_rule = "sweep_bounded"))), paste0(option_path, "$k_rule"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(k_rules_pending = "sweep_bounded"))),
        paste0(option_path, "$k_rules_pending"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(k_rules = c("num_sv_be", "sweep_bounded", "user")))),
        paste0(option_path, "$k_rules"))
    # Another rule's figure: of a rule the option offers besides the selected one, complete, proportions.
    other <- list(fdp = 0.06, fdp_where = "x", power = 0.9, source = "s")
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(k_rule_figures = list(num_sv_be = other)))),
        paste0(option_path, "$k_rule_figures"))
    expectRuvDesignRefused(option(\(o) utils::modifyList(o, list(k_rule_figures = list(sweep_bounded = other)))),
        paste0(option_path, "$k_rule_figures"))
    supplied <- function(edit) table(\(t) {
        t$options[[3L]]$k_rule_figures$sweep_bounded <- edit(t$options[[3L]]$k_rule_figures$sweep_bounded)
        t
    })
    figure_path <- paste0(table_path, "$options[[3]]$k_rule_figures$sweep_bounded")
    expectRuvDesignRefused(supplied(\(f) f[setdiff(names(f), "source")]), paste0(figure_path, "$source"))
    expectRuvDesignRefused(supplied(\(f) utils::modifyList(f, list(power = 90))), figure_path)
    # A group without the table still validates (the field is checked when present).
    without <- config
    without$norm$ruv_design <- NULL
    expect_identical(validateOmicGuiConfig(without), without)
})

# ------------------------------------------------------------------------------------------------
# The gate.
# ------------------------------------------------------------------------------------------------

test_that("the gate is off by default: no omic's tab shows the section, and its HTML is today's", {
    expect_null(getOption("MultiScholaR.ruv_design_panel"))
    for (omic in ruvDesignOmics) {
        norm <- omicGuiConfig(omic)$norm
        expect_false(norm$ruv_design$enabled, info = omic)
        html <- as.character(htmltools::renderTags(mod_omic_norm_ui("norm", norm))$html)
        expect_false(grepl("ruv_design|ruv_k_rule|ruv_k_user|ruv_k_estimate", html), info = omic)
        # The same bytes as the tab of a group that has no table at all.
        without <- norm[setdiff(names(norm), "ruv_design")]
        expect_identical(ruvDesignMaskTabs(html),
            ruvDesignMaskTabs(as.character(htmltools::renderTags(mod_omic_norm_ui("norm", without))$html)), info = omic)
    }
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    for (omic in ruvDesignOmics) expect_true(omicGuiConfig(omic)$norm$ruv_design$enabled, info = omic)
})

test_that("with the gate off the server registers nothing for the section", {
    updates <- list()
    testthat::local_mocked_bindings(updateRadioButtons = function(session, inputId, ...) {
        updates[[length(updates) + 1L]] <<- inputId
    }, .package = "shiny")
    bus <- shiny::reactiveValues(state_manager = NULL, design_matrix = NULL, config_list = list(),
        tab_status = list(normalization = "pending"))
    suppressWarnings(suppressMessages(shiny::testServer(mod_omic_norm_server, args = list(id = "norm", workflow_data = bus,
        experiment_paths = list(), omic_type = "lipidomics", experiment_label = "kit129",
        config = omicGuiConfig("lipidomics")$norm), {
        session$setInputs(ruv_design = "first_pass_nulls")
        session$flushReact()
    })))
    expect_length(updates, 0L)
})

# ------------------------------------------------------------------------------------------------
# The panel section (gate on).
# ------------------------------------------------------------------------------------------------

test_that("with the gate on each omic's tab shows every design with its figure and note, today's controls under today's design", {
    today <- list(
        proteomics = c("ruv_mode", "ruv_grouping_variable", "auto_percentage_min", "auto_percentage_max", "separation_metric",
            "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k"),
        transcriptomics = c("ruv_mode", "auto_percentage_min", "auto_percentage_max", "auto_percentage_step", "k_penalty_weight",
            "ruv_percentage", "ruv_k"))
    today$lipidomics <- today$metabolomics <- today$proteomics
    for (omic in ruvDesignOmics) {
        design <- omicGuiConfig(omic)$norm$ruv_design
        # (fourth pass item 2) as the tab opens: no technical replicates read yet
        table <- omicNormRuvMarkTable(design$tables[[1L]], technical = FALSE)
        html <- ruvDesignHtml(omic, TRUE)
        off <- ruvDesignHtml(omic, FALSE)
        # One radio group under a static id, one input per design in the table's order, today's design checked.
        radios <- regmatches(html, gregexpr('<input type="radio" name="norm-ruv_design"[^>]*>', html))[[1L]]
        expect_identical(sub('.*value="([^"]*)".*', "\\1", radios), vapply(table$options, `[[`, "", "id"), info = omic)
        expect_identical(sub('.*value="([^"]*)".*', "\\1", grep("checked", radios, value = TRUE)), table$selected, info = omic)
        expect_true(grepl(sprintf('data-testid="%s"', design$testid), html, fixed = TRUE), info = omic)
        # A design that is not available is shown disabled, with its reason; no other is disabled.
        unavailable <- vapply(table$options, \(option) identical(option$status, "unavailable") || !is.null(option$unavailable_now), NA)
        expect_identical(grepl("disabled", radios), unavailable, info = omic)
        for (option in table$options) {
            where <- paste(omic, option$id)
            # The tier, the label, the description, the figure and the note are text in the panel; the tooltip is
            # the option's title.
            for (text in c(design$tiers[[option$tier]], option$label, option$description, option$note,
                omicNormRuvFigureText(option), omicNormRuvAltText(option))) {
                expect_true(grepl(htmltools::htmlEscape(text), html, fixed = TRUE), info = paste(where, text))
            }
            expect_true(grepl(sprintf('title="%s"', htmltools::htmlEscape(option$tooltip, attribute = TRUE)), html,
                fixed = TRUE), info = where)
            expect_false(grepl(htmltools::htmlEscape(option$note), off, fixed = TRUE), info = where)
        }
        # Today's controls keep their ids and are shown under today's design only.
        ids <- sub('^id="norm-', "", sub('"$', "", regmatches(html, gregexpr('id="norm-[A-Za-z0-9_]+"', html))[[1L]]))
        expect_identical(setdiff(today[[omic]], ids), character(), info = omic)
        legacy <- regmatches(html, regexpr('data-display-if="[^"]*"[^>]*>\\s*<hr/>\\s*<h4>RUV', html))
        expect_length(legacy, 1L)
        expect_match(legacy, "legacy", fixed = TRUE)
        expect_false(grepl("groups_all_features|none", legacy), info = omic)
        # The k rule, the number and the static slots (the summary under the option; one diagnostic per assay row).
        expect_true(all(c("ruv_design", "ruv_k_rule", "ruv_k_user", "ruv_k_seed", "ruv_control_list", "ruv_design_summary") %in% ids),
            info = omic)
        slots <- sprintf("ruv_k_estimate_assay%d", seq_len(omicGuiConfig(omic)$norm$assay_slots))
        expect_identical(sort(grep("^ruv_k_estimate", ids, value = TRUE)), slots, info = omic)
        for (rule in design$k_rules) {
            for (text in rule[c("label", "reads", "pros", "cons", "diagnostic")]) {
                expect_true(grepl(htmltools::htmlEscape(text), html, fixed = TRUE), info = paste(omic, text))
            }
        }
        # Every other id of today's tab is still there, and no new output beyond the static slots.
        off_ids <- sub('^id="norm-', "", sub('"$', "", regmatches(off, gregexpr('id="norm-[A-Za-z0-9_]+"', off))[[1L]]))
        expect_identical(setdiff(off_ids, ids), character(), info = omic)
        expect_identical(sort(setdiff(grep("-label$", ids, invert = TRUE, value = TRUE), off_ids)),
            # (fifth pass items 4, 5) the separation plot per assay row where the omic has several; the check's table
            sort(c("ruv_design", "ruv_k_rule", "ruv_k_user", "ruv_k_seed", "ruv_control_list", "ruv_design_summary",
                "ruv_design_association_table", if (length(slots) > 1L) sub("estimate", "separation", slots) else "ruv_k_separation_plot",
                slots)), info = omic)
    }
})

test_that("with the gate on the section is a block above the options column, its designs a comparison list", {
    for (omic in ruvDesignOmics) {
        design <- omicGuiConfig(omic)$norm$ruv_design
        html <- ruvDesignHtml(omic, TRUE)
        # One labelled block, before the options column; today's RUV section is in it and no longer in the column.
        block <- regexpr('<section class="well ruv-design" aria-label="Unwanted variation (RUV)">', html, fixed = TRUE)
        column <- regexpr('<div class="col-sm-3">', html, fixed = TRUE)
        expect_true(block > 0L && block < column, info = omic)
        expect_length(gregexpr("RUV Parameter Tuning:", html, fixed = TRUE)[[1L]], 1L)
        expect_lt(regexpr("RUV Parameter Tuning:", html, fixed = TRUE), column)
        # The column headings once, between the group's label and its options; they are not read out.
        header <- regexpr('<div class="ruv-design-header" aria-hidden="true">', html, fixed = TRUE)
        expect_true(header > regexpr("How unwanted variation is removed:", html, fixed = TRUE) &&
            header < regexpr('class="shiny-options-group"', html, fixed = TRUE), info = omic)
        for (option in omicNormRuvMarkTable(design$tables[[1L]], technical = FALSE)$options) {
            row <- as.character(omicNormRuvOptionTag(option, design$tiers))
            expect_true(grepl(row, html, fixed = TRUE), info = paste(omic, option$id))
            # A measured figure in two cells for the eye and as its sentence for a screen reader; else its status.
            if (identical(option$status, "measured")) {
                parts <- omicNormRuvFigureParts(option)
                expect_match(row, sprintf('<span class="ruv-design-fdp" aria-hidden="true">\\s*<strong>%s</strong>',
                    parts[["fdp"]]), info = option$id)
                expect_match(row, sprintf('<span class="ruv-design-power" aria-hidden="true">\\s*<strong>%s</strong>',
                    parts[["power"]]), info = option$id)
                expect_match(row, paste0('<span class="sr-only">', htmltools::htmlEscape(omicNormRuvFigureText(option)),
                    "</span>"), fixed = TRUE, info = option$id)
            } else {
                expect_match(row, sprintf('class="ruv-design-option is-%s"', option$status), fixed = TRUE)
                expect_match(row, paste0('<span class="ruv-design-status">',
                    htmltools::htmlEscape(omicNormRuvFigureText(option)), "</span>"), fixed = TRUE, info = option$id)
                expect_false(grepl("aria-hidden", row, fixed = TRUE), info = option$id)
            }
        }
        # One panel per k rule, side by side, behind one disclosure that is closed until opened; the choice of the
        # rule and the number are before it, outside.
        expect_length(gregexpr('ruv-design-k-rule"', html, fixed = TRUE)[[1L]], length(design$k_rules))
        opened <- regexpr('<details class="ruv-design-k-why">\\s*<summary>How is this chosen\\?', html)
        closed <- regexpr("</details>", html, fixed = TRUE)
        expect_length(gregexpr("<details", html, fixed = TRUE)[[1L]], 1L)
        expect_false(grepl("<details[^>]* open", html), info = omic)
        rules <- gregexpr('ruv-design-k-rule"', html, fixed = TRUE)[[1L]]
        expect_true(all(rules > opened & rules < closed), info = omic)
        for (id in c('id="norm-ruv_k_rule"', 'id="norm-ruv_k_user"')) {
            expect_true(regexpr(id, html, fixed = TRUE) < opened, info = paste(omic, id))
        }
        # The headings say which scenario each figure is from; no row repeats it.
        expect_match(html, "target 5%; worst scenario tested", fixed = TRUE)
        expect_match(html, "real differences found; reference scenario", fixed = TRUE)
        expect_false(grepl('<span class="ruv-design-where">(worst scenario tested:|in the reference scenario)', html), info = omic)
    }
})

test_that("the figure text is the table's numbers in words; a pending or unavailable cell shows no number", {
    expect_identical(omicNormRuvFigureText(list(status = "measured", fdp = 0.066,
        fdp_where = "12 samples per group", power = 0.91)),
        paste("Measured in simulation: 6.6% false discoveries at a 5% FDR threshold (worst scenario tested: 12 samples",
            "per group); power 91% in the reference scenario."))
    expect_identical(omicNormRuvFigureText(list(status = "measured", fdp = 0.044, fdp_where = "x", power = 0.005,
        single_scenario = TRUE)),
        "Measured in simulation: 4.4% false discoveries at a 5% FDR threshold (only scenario tested: x); power 0.5% in that scenario.")
    # A row shows the scenario only (the heading says it is the worst tested); a single scenario says so.
    expect_identical(omicNormRuvFigureParts(list(fdp = 0.066, fdp_where = "x", power = 0.91))[c("fdp_where", "power_where")],
        c(fdp_where = "x", power_where = ""))
    expect_identical(omicNormRuvFigureParts(list(fdp = 0.066, fdp_where = "x", power = 0.91, single_scenario = TRUE))[c("fdp_where",
        "power_where")], c(fdp_where = "only scenario tested: x", power_where = "in that scenario"))
    # The power of every table as a percentage, exactly.
    powers <- c(0.91, 0.90, 0.92, 0.49, 0.77, 0.56, 0.76, 0.60, 0.61, 0.73, 0.005, 0.47, 0.34, 0.40, 0.26, 0.48, 0.53, 0.42,
        0.003)
    expect_identical(vapply(powers, \(power) omicNormRuvFigureParts(list(fdp = 0, fdp_where = "", power = power))[["power"]], ""),
        c("91%", "90%", "92%", "49%", "77%", "56%", "76%", "60%", "61%", "73%", "0.5%", "47%", "34%", "40%", "26%", "48%", "53%",
            "42%", "0.3%"))
    expect_identical(omicNormRuvFigureText(list(status = "pending")), "Not yet measured in simulation.")
    # The figure with the other rule for k, in words; none where there is none.
    expect_identical(omicNormRuvAltText(list(k_rule_figures = list(sweep_bounded = list(fdp = 0.118, fdp_where = "x",
        power = 0.26)))), "With the separation rule for k: 11.8% false discoveries; worst scenario: x. Power 26%.")
    expect_identical(omicNormRuvAltText(list(status = "measured")), character())
    expect_identical(omicNormRuvFigureText(list(status = "unavailable", unavailable = "No replicates.")),
        "Not available: No replicates.")
})

test_that("the k rule's choices follow the selected design, and the route's table replaces the options", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    updates <- list()
    testthat::local_mocked_bindings(updateRadioButtons = function(session, inputId, label = NULL, choices = NULL,
        selected = NULL, inline = FALSE, choiceNames = NULL, choiceValues = NULL) {
        updates[[length(updates) + 1L]] <<- list(id = inputId, choices = choices, selected = selected,
            values = choiceValues, names = choiceNames)
    }, .package = "shiny")
    bus <- shiny::reactiveValues(state_manager = NULL, design_matrix = NULL, config_list = list(),
        tab_status = list(normalization = "pending"), state_update_trigger = NULL)
    suppressWarnings(suppressMessages(shiny::testServer(mod_omic_norm_server, args = list(id = "norm", workflow_data = bus,
        experiment_paths = list(), omic_type = "lipidomics", experiment_label = "kit129",
        config = omicGuiConfig("lipidomics")$norm), {
        session$setInputs(ruv_design = "first_pass_nulls")
        session$setInputs(ruv_design = "groups_all_features")
        session$setInputs(ruv_design = "legacy")
    })))
    expect_identical(lapply(updates, `[[`, "id"), list("ruv_k_rule", "ruv_k_rule"))
    expect_identical(unname(updates[[1L]]$choices), c("num_sv_be", "sweep_bounded", "user"))
    expect_identical(updates[[1L]]$selected, "num_sv_be")
    expect_identical(unname(updates[[2L]]$choices), c("num_sv_be", "user"))
    expect_identical(updates[[2L]]$selected, "num_sv_be")
    expect_identical(names(updates[[2L]]$choices), c("Estimated number", "Your own number"))

    # Proteomics: the DIA with limpa route has its own table; the other routes keep the first.
    updates <- list()
    norm <- omicGuiConfig("proteomics")$norm
    norm$server <- function(input, output, session, context) NULL
    manager <- new.env()
    manager$workflow_type <- "LFQ"
    bus <- shiny::reactiveValues(state_manager = manager, design_matrix = NULL, config_list = list(),
        tab_status = list(normalization = "pending"), state_update_trigger = NULL)
    suppressWarnings(suppressMessages(shiny::testServer(mod_omic_norm_server, args = list(id = "normalization",
        workflow_data = bus, experiment_paths = list(), omic_type = "proteomics", experiment_label = "kit129", config = norm), {
        session$flushReact()
        expect_length(updates, 0L)
        manager$workflow_type <- "DIA_LIMPA"
        bus$state_update_trigger <- Sys.time()
        session$flushReact()
    })))
    expect_identical(lapply(updates, `[[`, "id"), list("ruv_design"))
    limpa <- norm$ruv_design$tables$dia_limpa
    expect_identical(updates[[1L]]$values, vapply(limpa$options, `[[`, "", "id"))
    expect_identical(updates[[1L]]$selected, limpa$selected)
    shown <- vapply(updates[[1L]]$names, as.character, "")
    for (i in seq_along(limpa$options)) {
        expect_true(grepl(htmltools::htmlEscape(limpa$options[[i]]$note), shown[[i]], fixed = TRUE), info = limpa$options[[i]]$id)
        expect_true(grepl(htmltools::htmlEscape(omicNormRuvFigureText(limpa$options[[i]])), shown[[i]], fixed = TRUE))
    }
})

# ------------------------------------------------------------------------------------------------
# KIT-129 fourth pass item 2: the panel's choices build the run path's request.
# ------------------------------------------------------------------------------------------------

# The options a session cannot choose, by table and by whether its design declares technical replicates, with the
# text of the table that says why (written from the ticket's answers: M6, U2, U3; and the run path as it stands).
# KIT-129 fifth pass item 3: supplied lists are on the classic and lipid/metab run paths; only the preconditions close.
ruvDesignClosed <- list(
    proteomics.classic = list(no = c(technical_replicates = "no_technical_replicates"), yes = character()),
    proteomics.dia_limpa = list(no = c(technical_replicates = "no_technical_replicates"),
        yes = c(groups_all_features = "technical_replicates_declared")),
    transcriptomics.genes = list(no = character(), yes = character()),
    lipidomics.features = list(no = c(technical_replicates = "no_technical_replicates"),
        yes = c(groups_all_features = "technical_replicates_declared")))
ruvDesignClosed$metabolomics.features <- ruvDesignClosed$lipidomics.features

test_that("(fourth pass item 2) only the options a route's run path takes can be chosen; the others say why; the defaults are M2 and M6", {
    tables <- ruvDesignTables()
    for (name in names(tables)) {
        for (technical in c(FALSE, TRUE)) {
            marked <- omicNormRuvMarkTable(tables[[name]], technical)
            closed <- Filter(Negate(is.null), stats::setNames(lapply(marked$options, `[[`, "unavailable_now"),
                vapply(marked$options, `[[`, "", "id")))
            expected <- ruvDesignClosed[[name]][[if (technical) "yes" else "no"]]
            expect_identical(names(closed) %||% character(), names(expected) %||% character(), info = paste(name, technical))
            expect_identical(unname(unlist(closed)) %||% character(), unname(marked$texts[expected]), info = paste(name, technical))
            # the recommended design, or the technical-replicate design where the design declares them (not trx)
            expect_identical(marked$selected, if (technical && !startsWith(name, "transcriptomics")) "technical_replicates" else
                "groups_all_features", info = paste(name, technical))
        }
    }
    # transcriptomics: technical replicates disabled with the table's reason (M6, U3); (fifth pass item 4) the separation
    # rule is on every route's run path where it is defined
    for (name in names(tables)) expect_null(tables[[name]]$k_rules_not_offered, info = name)
    # every new text is plain words: no internal names
    internal <- "\\b(arm|arms|KIT-[0-9]+|C[0-9]|S[0-9]{1,2}|U[0-9]{1,2}|W[0-9]{1,2}|R0|B[0-9]|num_sv_be|sweep_bounded)\\b"
    for (name in names(tables)) expect_identical(grep(internal, unname(tables[[name]]$texts), value = TRUE), character(), info = name)
    # as rendered: a closed option is disabled and its row says why, and keeps its figure; (fifth pass item 3) the
    # supplied list is open in classic proteomics
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    html <- ruvDesignHtml("proteomics", TRUE)
    radio <- function(id) regmatches(html, regexpr(sprintf('<input type="radio" name="norm-ruv_design" value="%s"[^>]*>', id), html))
    expect_match(radio("technical_replicates"), "disabled", fixed = TRUE)
    expect_false(grepl("disabled", radio("supplied_controls"), fixed = TRUE))
    expect_match(html, htmltools::htmlEscape(paste("Not available:", tables[["proteomics.classic"]]$texts[["no_technical_replicates"]])),
        fixed = TRUE)
    expect_match(html, '<span class="ruv-design-fdp" aria-hidden="true">\\s*<strong>5.4%</strong>')
})

ruvDesignInput <- function(...) {
    utils::modifyList(list(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_user = NA, ruv_k_seed = 7,
        ruv_control_list = "", norm_method = "quantile", ruv_mode = "automatic", auto_percentage_min = 1, auto_percentage_max = 20,
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_percentage = 5, ruv_k = 2,
        color_variable = "group", shape_variable = "group", ruv_grouping_variable = "group", apply_itsd = FALSE,
        itsd_aggregation = "sum", log_offset = 1, auto_percentage_step = 5), list(...))
}

test_that("(fourth pass item 2) the panel's choice reads the design, the k rule, the typed k, the seed and the list, and refuses in words", {
    table <- omicNormRuvRunTable("proteomics", "DIA_LIMPA")
    expect_null(omicNormRuvDesignChoice(ruvDesignInput(), table))
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(), table), list(ruv_design = "groups_all_features", legacy = FALSE,
        ruv_k_rule = "num_sv_be", ruv_k = NULL, ruv_k_seed = 7L, ruv_control_ids = NULL))
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(ruv_k_rule = "user", ruv_k_user = 2), table)[c("ruv_k_rule", "ruv_k", "ruv_k_seed")],
        list(ruv_k_rule = "user", ruv_k = 2L, ruv_k_seed = 7L))
    supplied <- omicNormRuvDesignChoice(ruvDesignInput(ruv_design = "supplied_controls", ruv_k_rule = "sweep_bounded",
        ruv_control_list = " P2, P1\nP3;P1 "), table)
    expect_identical(supplied[c("ruv_k_rule", "ruv_k_seed", "ruv_control_ids")], list(ruv_k_rule = "sweep_bounded", ruv_k_seed = NULL,
        ruv_control_ids = c("P2", "P1", "P3")))
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(ruv_design = "legacy_correction"), table),
        list(ruv_design = "legacy_correction", legacy = TRUE))
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(ruv_design = "none"), table)[c("ruv_design", "ruv_k_rule", "ruv_k")],
        list(ruv_design = "none", ruv_k_rule = NULL, ruv_k = NULL))
    texts <- table$texts
    refused <- function(input, table, text) expect_error(omicNormRuvDesignChoice(input, table), text, fixed = TRUE,
        class = "omic_norm_ruv_design_invalid")
    refused(ruvDesignInput(ruv_k_rule = "user", ruv_k_user = NA), table, texts[["k_missing"]])
    refused(ruvDesignInput(ruv_k_rule = "user", ruv_k_user = 1.5), table, texts[["k_missing"]])
    refused(ruvDesignInput(ruv_k_seed = 1.5), table, texts[["seed_invalid"]])
    refused(ruvDesignInput(ruv_design = "supplied_controls", ruv_control_list = " \n "), table, texts[["list_missing"]])
    refused(ruvDesignInput(ruv_k_rule = "sweep_bounded"), table, texts[["k_rule_not_available"]])
    classic <- omicNormRuvRunTable("proteomics", "LFQ")
    # (fifth pass item 3) the supplied list is read in classic proteomics and lipid/metab, the identifiers in entry order
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(ruv_design = "supplied_controls", ruv_control_list = "P2 P1"),
        classic)$ruv_control_ids, c("P2", "P1"))
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(ruv_design = "supplied_controls", ruv_control_list = "L2;L1"),
        omicNormRuvRunTable("lipidomics"))$ruv_control_ids, c("L2", "L1"))
    refused(ruvDesignInput(ruv_design = "legacy_correction"), classic, texts[["not_available"]])
    trx <- omicNormRuvRunTable("transcriptomics")
    refused(ruvDesignInput(ruv_design = "technical_replicates"), trx, trx$options[[2L]]$unavailable)
    expect_identical(omicNormRuvDesignChoice(ruvDesignInput(ruv_design = "first_pass_nulls", ruv_k_rule = "sweep_bounded"), trx)$ruv_k_rule,
        "sweep_bounded")
    refused(ruvDesignInput(ruv_k_rule = "sweep_bounded"), trx, texts[["k_rule_not_available"]])
})

test_that("(fourth pass item 2) each omic's tab builds the run path's request from the choice; legacy builds today's; the gate off changes nothing", {
    today <- list()
    withr::with_options(list(MultiScholaR.ruv_design_panel = NULL), {
        input <- ruvDesignInput()
        today$classic <- applyProtNormRuvDesign(resolveProtNormReferenceSettings(input),
            omicNormRuvDesignChoice(input, omicNormRuvRunTable("proteomics", "LFQ")))
        today$lipid <- omicNormRunSettings(input, list(color_var = "group", shape_var = "group"), list(assay_names = NULL), list(),
            lipidomicsDeclarationBinding, ruvDesign = omicNormRuvDesignChoice(input, omicNormRuvRunTable("lipidomics")))
        today$trx <- trxNormGridParameters(input)
    })
    expect_identical(today$classic, resolveProtNormReferenceSettings(ruvDesignInput()))
    expect_identical(today$trx, trxRuvGridParameters("automatic", 1, 20, 5, 0.5))
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    norm_data <- new.env()
    norm_data$normalization_input_reference <- list(entity_id = "source")
    classic <- function(...) {
        input <- ruvDesignInput(...)
        settings <- applyProtNormRuvDesign(resolveProtNormReferenceSettings(input), omicNormRuvDesignChoice(input,
            omicNormRuvRunTable("proteomics", "LFQ")))
        diaNormalizationRequest(NULL, settings, norm_data, "group")$parameters
    }
    # classic proteomics: legacy is today's request; the new designs carry the design, rule and seed, the mode from the rule
    expect_identical(classic(ruv_design = "legacy"), diaNormalizationRequest(NULL, today$classic, norm_data, "group")$parameters)
    counted <- classic()
    expect_identical(counted[c("ruv_mode", "ruv_design", "ruv_k_rule", "ruv_k_seed")], list(ruv_mode = "automatic",
        ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = 7L))
    expect_identical(workflowRequestOperation(NULL, workflowOperationEntry("normalization-ruv"), NULL, counted)$declaration$parameters,
        diaNormalizationOperation("declared")$parameters)
    typed <- classic(ruv_design = "technical_replicates", ruv_k_rule = "user", ruv_k_user = 3)
    expect_identical(typed[c("ruv_mode", "ruv_k", "ruv_design", "ruv_k_rule")], list(ruv_mode = "manual", ruv_k = 3L,
        ruv_design = "technical_replicates", ruv_k_rule = "user"))
    none <- classic(ruv_design = "none")
    expect_identical(none[c("ruv_mode", "ruv_design")], list(ruv_mode = "skip", ruv_design = "none"))
    expect_null(none$ruv_k_rule)
    expect_null(none$ruv_control_ids)
    # (fifth pass item 3) the supplied list travels beside the design; the request's identity is its hash and count
    supplied <- classic(ruv_design = "supplied_controls", ruv_control_list = "P2 P1 P2")
    expect_identical(supplied[c("ruv_design", "ruv_k_rule", "ruv_control_ids")], list(ruv_design = "supplied_controls",
        ruv_k_rule = "num_sv_be", ruv_control_ids = c("P2", "P1")))
    expect_null(counted$ruv_control_ids)
    chosen <- workflowRequestOperation(NULL, workflowOperationEntry("normalization-ruv"), NULL, supplied)
    expect_identical(chosen$declaration$parameters, c(diaNormalizationOperation("declared")$parameters, RUV_SUPPLIED_LIST_PARAMETERS))
    expect_false("ruv_control_ids" %in% chosen$declaration$parameters)
    expect_identical(chosen$parameters[RUV_SUPPLIED_LIST_PARAMETERS], list(ruv_control_list_sha256 = trxControlListSha256(c("P1", "P2")),
        ruv_control_count = 2L))
    # DIA with limpa: the mode and the declared rules
    dia <- function(...) {
        choice <- omicNormRuvDesignChoice(ruvDesignInput(...), omicNormRuvRunTable("proteomics", "DIA_LIMPA"))
        settings <- c(list(ruv_mode = protNormDiaRuvDesignMode(choice)), protNormDiaRuvDesignRules(choice))
        # (the skip mode records a skip on the tab and sends no request)
        list(mode = settings$ruv_mode, parameters = if (!identical(settings$ruv_mode, "skip"))
            diaRuvCovariatesRequestParameters(list(), "group", settings))
    }
    expect_identical(dia(ruv_design = "legacy"), list(mode = "covariates", parameters = resolveDIARuvCovariatesParameters(list(), "group")))
    expect_identical(dia(ruv_design = "legacy_correction")$mode, "correction")
    expect_identical(dia(ruv_design = "none")$mode, "skip")
    expect_identical(dia()$parameters, resolveDIARuvCovariatesParameters(list(), "group", "covariates", ruv_control_rule = "all_features",
        ruv_k_rule = "num_sv_be", ruv_k_seed = 7L))
    expect_identical(dia(ruv_design = "technical_replicates")$parameters, dia()$parameters)
    expect_identical(dia(ruv_design = "supplied_controls", ruv_k_rule = "user", ruv_k_user = 1, ruv_control_list = "P2 P1")$parameters,
        resolveDIARuvCovariatesParameters(list(), "group", "covariates", ruv_control_rule = "supplied_list", ruv_k_rule = "user", ruv_k = 1L,
            ruv_control_ids = c("P2", "P1"), ruv_k_seed = 7L))
    expect_identical(dia(ruv_design = "first_pass_nulls")$parameters$ruv_control_rule, "first_pass_drop_top_25")
    # lipidomics and metabolomics: the RUV step's request names the new design (its mode never skip); legacy is today's
    for (binding in list(lipidomicsDeclarationBinding, metabolomicsDeclarationBinding)) {
        omics <- function(...) {
            input <- ruvDesignInput(...)
            omicNormRunSettings(input, list(color_var = "group", shape_var = "group"), list(assay_names = NULL), list(), binding,
                ruvDesign = omicNormRuvDesignChoice(input, omicNormRuvRunTable(binding$omic_type)))
        }
        legacy <- omics(ruv_design = "legacy")
        expect_identical(legacy, today$lipid, info = binding$prefix)
        new <- omics()
        expect_identical(new$ruv_mode, "groups_all_features", info = binding$prefix)
        step <- omicsRuvStep(binding, new$ruv_parameters)
        expect_identical(step$parameters, list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = NULL,
            ruv_k_rule = "num_sv_be", ruv_k_seed = 7L, ruv_control_ids = NULL), info = binding$prefix)
        expect_identical(workflowRequestOperation(NULL, workflowOperationEntry(step$operation), NULL, step$parameters)$declaration$parameters,
            omicsRuvDesignParameterNames, info = binding$prefix)
        technical <- omicsRuvStep(binding, omics(ruv_design = "technical_replicates", ruv_k_rule = "user", ruv_k_user = 1)$ruv_parameters)
        expect_identical(technical$parameters[c("ruv_design", "ruv_grouping_variable", "ruv_k", "ruv_k_rule")],
            list(ruv_design = "technical_replicates", ruv_grouping_variable = NULL, ruv_k = 1L, ruv_k_rule = "user"), info = binding$prefix)
        expect_identical(omics(ruv_design = "none")$ruv_mode, "none", info = binding$prefix)
        # (fifth pass item 3) the one supplied list, by the hash and the count in the identity
        listed <- omicsRuvStep(binding, omics(ruv_design = "supplied_controls", ruv_control_list = "f2 f1")$ruv_parameters)
        expect_identical(listed$parameters$ruv_control_ids, c("f2", "f1"), info = binding$prefix)
        chosen <- workflowRequestOperation(NULL, workflowOperationEntry(listed$operation), NULL, listed$parameters)
        expect_identical(chosen$declaration$parameters, c(setdiff(omicsRuvDesignParameterNames, "ruv_control_ids"),
            RUV_SUPPLIED_LIST_PARAMETERS), info = binding$prefix)
        expect_identical(chosen$parameters$ruv_control_count, 2L, info = binding$prefix)
    }
    # transcriptomics: the declared design alone, in its canonical form; legacy is the grid's request
    trx <- function(...) trxNormRuvRequest(trxNormGridParameters(ruvDesignInput(...)))
    expect_identical(trx(ruv_design = "legacy"), list(grid = today$trx, design = NULL))
    expect_identical(trx()$design, trxRuvDesignParameters("groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = 7L))
    expect_identical(trx(ruv_k_rule = "user", ruv_k_user = 2)$design, trxRuvDesignParameters("groups_all_features", 2L, "user"))
    expect_identical(trx(ruv_design = "supplied_controls", ruv_control_list = "g2,g1")$design,
        trxRuvDesignParameters("supplied_controls", ruv_k_rule = "num_sv_be", ruv_control_ids = c("g1", "g2"), ruv_k_seed = 7L))
    expect_identical(trx(ruv_design = "none")$design, trxRuvDesignParameters("none"))
})

test_that("(fourth pass item 2) on a design with technical replicates the panel selects that design and refuses a closed one", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    updates <- list()
    notices <- list()
    testthat::local_mocked_bindings(updateRadioButtons = function(session, inputId, label = NULL, choices = NULL,
        selected = NULL, inline = FALSE, choiceNames = NULL, choiceValues = NULL) {
        updates[[length(updates) + 1L]] <<- list(id = inputId, selected = selected, names = choiceNames)
    }, showNotification = function(ui, ...) notices[[length(notices) + 1L]] <<- ui, .package = "shiny")
    head <- new.env()
    head$technical <- FALSE
    testthat::local_mocked_bindings(omicNormRuvHeadTechnicalReplicates = function(manager) head$technical, .package = "MultiScholaR")
    bus <- shiny::reactiveValues(state_manager = NULL, design_matrix = NULL, config_list = list(),
        tab_status = list(normalization = "pending"), state_update_trigger = NULL)
    table <- omicGuiConfig("lipidomics")$norm$ruv_design$tables$features
    suppressWarnings(suppressMessages(shiny::testServer(mod_omic_norm_server, args = list(id = "norm", workflow_data = bus,
        experiment_paths = list(), omic_type = "lipidomics", experiment_label = "kit129",
        config = omicGuiConfig("lipidomics")$norm), {
        session$flushReact()
        expect_length(updates, 0L)
        head$technical <- TRUE
        bus$state_update_trigger <- 1
        session$flushReact()
        expect_identical(updates[[1L]][c("id", "selected")], list(id = "ruv_design", selected = "technical_replicates"))
        shown <- vapply(updates[[1L]]$names, as.character, "")
        expect_match(shown[[1L]], htmltools::htmlEscape(paste("Not available:", table$texts[["technical_replicates_declared"]])),
            fixed = TRUE)
        expect_false(grepl("Not available", shown[[2L]], fixed = TRUE))
        updates <<- list()
        session$setInputs(ruv_design = "groups_all_features")
        expect_identical(updates[[1L]][c("id", "selected")], list(id = "ruv_design", selected = "technical_replicates"))
        expect_identical(notices[[1L]], paste(table$texts[["not_available"]], table$texts[["technical_replicates_declared"]]))
        # the same state again: nothing is shown again
        updates <<- list()
        bus$state_update_trigger <- 2
        session$flushReact()
        expect_length(updates, 0L)
    })))
})

test_that("(fourth pass item 2) a refused choice is the tab's error notice and the run does not start", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    notices <- list()
    lines <- character()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        notices[[length(notices) + 1L]] <<- list(ui = ui, type = type)
    }, .package = "shiny")
    read <- omicNormReadRuvDesign(ruvDesignInput(ruv_k_rule = "user", ruv_k_user = NA), "lipidomics", function(line) lines <<- c(lines, line))
    texts <- omicNormRuvRunTable("lipidomics")$texts
    expect_identical(read, list(ok = FALSE, choice = NULL))
    expect_identical(notices, list(list(ui = texts[["k_missing"]], type = "error")))
    expect_identical(lines, paste("ERROR:", texts[["k_missing"]]))
    expect_true(omicNormReadRuvDesign(ruvDesignInput(), "lipidomics", function(line) NULL)$ok)
})

# ------------------------------------------------------------------------------------------------
# KIT-129 fourth pass item 3: after the run the panel shows what was done.
# ------------------------------------------------------------------------------------------------

test_that("(fourth pass item 3) each omic's product gives one result: the design, k and its rule, the count beside a typed k, the sets", {
    classic <- protNormRuvDesignResult(list(parameters = list(ruv_design = "groups_all_features"), best_k = 2L,
        ruv_model = list(replicate_structure = "group"), ruv_result = list(ruv_k_rule = "user", ruv_k_estimate = 1L,
            ruv_applied = TRUE, ruv_k_seed = 7L, k_estimate = list(k = 1L), ruv_not_applied_reason = NA_character_)))
    expect_identical(classic[c("design", "k_rule", "seed", "replicate_structure")], list(design = "groups_all_features",
        k_rule = "user", seed = 7L, replicate_structure = "group"))
    expect_identical(classic$assays$all[c("k", "k_estimate", "applied", "estimate")], list(k = 2L, k_estimate = 1L, applied = TRUE,
        estimate = list(k = 1L)))
    expect_null(classic$separation_plot)
    # (item 4) a legacy run is shown as such, with the check not available
    expect_true(protNormRuvDesignResult(list(parameters = list(norm_method = "quantile"), ruv_result = list()))$legacy)
    expect_identical(protNormRuvDesignResult(list(parameters = list(ruv_design = "none")))$design, "none")
    dia <- diaRuvDesignResult(list(ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be", ruv_k_estimate = 1L, best_k = 1L,
        ruv_applied = TRUE, replicate_structure = "technical_replicates", parameters = list(ruv_k_seed = 7L),
        ruv_result = list(k_estimate = list(k = 1L), k_sweep = list(cancor_plot = "plot"))), "technical_replicates")
    expect_identical(dia[c("design", "replicate_structure", "separation_plot")], list(design = "technical_replicates",
        replicate_structure = "technical_replicates", separation_plot = "plot"))
    expect_identical(diaRuvDesignResult(list(ruv_result = list(), ruv_mode = "correction"), NULL)[c("design", "legacy")],
        list(design = "legacy_correction", legacy = TRUE))
    omics <- omicsRuvDesignResult(list(ruv_design = "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = 7L,
        ruv_replicate_structure = "group", ruv_applied = list(pos = TRUE, neg = FALSE), k_estimate = list(pos = list(k = 1L), neg = NULL),
        factors = list(pos = list(k = 1L, ruv_k_estimate = 1L, reason = NA_character_), neg = list(k = 0L, ruv_k_estimate = 0L,
            reason = "the count found no unwanted factor"))))
    expect_identical(names(omics$assays), c("pos", "neg"))
    expect_identical(omics$assays$neg[c("applied", "reason")], list(applied = FALSE, reason = "the count found no unwanted factor"))
    expect_true(omicsRuvDesignResult(list(best_k = list(pos = 1L)))$legacy)
    trx <- trxRuvDesignResult(list(ruv_design = list(id = "supplied_controls", method = "RUVg", k = 2L, k_rule = "user", applied = TRUE,
        reason = NA_character_), k_record = list(k_estimate = 1L, k_seed = 7L), k_estimate = list(k = 1L)))
    expect_identical(trx[c("design", "k_rule", "seed", "replicate_structure")], list(design = "supplied_controls", k_rule = "user",
        seed = 7L, replicate_structure = NA_character_))
    expect_true(trxRuvDesignResult(list(chosen = list(k = 1L)))$legacy)
})

test_that("(fourth pass item 3) the summary says what was done in the table's words", {
    design <- omicGuiConfig("lipidomics")$norm$ruv_design
    presets <- omicNormRuvPresets(design)
    texts <- design$tables[[1L]]$texts
    html <- function(result) as.character(omicNormRuvResultTags(result, presets, texts))
    expect_identical(html(NULL), character(0))
    typed <- html(omicNormRuvDesignResult("groups_all_features", "user", 7L, "group", list(pos = list(k = 2L, k_estimate = 1L,
        applied = TRUE), neg = list(k = 0L, applied = FALSE, reason = "k = 3 is above the limit of 2"))))
    expect_match(typed, htmltools::htmlEscape(texts[["result_heading"]]), fixed = TRUE)
    expect_match(typed, htmltools::htmlEscape(presets$groups_all_features$label), fixed = TRUE)
    expect_match(typed, paste("Assay pos -", texts[["result_k"]], "2 (your own number; estimated number for comparison: 1)"), fixed = TRUE)
    expect_match(typed, paste("Assay neg -", texts[["result_not_applied"]], "k = 3 is above the limit of 2"), fixed = TRUE)
    expect_match(typed, htmltools::htmlEscape(texts[["result_sets_group"]]), fixed = TRUE)
    expect_match(typed, texts[["result_seed"]], fixed = TRUE)
    counted <- html(omicNormRuvDesignResult("technical_replicates", "num_sv_be", 7L, "technical_replicates",
        list(all = list(k = 1L, applied = TRUE))))
    expect_match(counted, paste(texts[["result_k"]], "1 (estimated number)"), fixed = TRUE)
    expect_false(grepl("Assay", counted, fixed = TRUE))
    expect_match(counted, htmltools::htmlEscape(texts[["result_sets_technical_replicates"]]), fixed = TRUE)
    none <- html(omicNormRuvDesignResult("none", NA, NA, NA, list()))
    expect_match(none, htmltools::htmlEscape(presets$none$label), fixed = TRUE)
    expect_false(grepl(texts[["result_k"]], none, fixed = TRUE))
    # (fifth pass item 3) a supplied list: the identifiers listed and found, the controls used per assay
    listed <- html(omicNormRuvDesignResult("supplied_controls", "num_sv_be", 7L, "group", list(pos = list(k = 1L, applied = TRUE,
        control_count = 8L), neg = list(k = 1L, applied = TRUE, control_count = 5L)), control_list = list(listed = 15L, unmatched = 2L)))
    expect_match(listed, sprintf(texts[["result_list"]], 15L, 13L), fixed = TRUE)
    expect_match(listed, paste("Assay neg -", htmltools::htmlEscape(texts[["result_controls"]]), "5"), fixed = TRUE)
    expect_false(grepl(sprintf(texts[["result_list"]], 15L, 13L), typed, fixed = TRUE))
})

test_that("(fourth pass item 3) the tab renders the last run's summary and the estimated number's plot; the gate off sets nothing", {
    state <- new.env()
    omicNormSetRuvDesignResult(state, list(design = "x"))
    expect_null(state$ruv_design_result)
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    estimate <- NULL
    if (requireNamespace("sva", quietly = TRUE)) {
        set.seed(1)
        y <- matrix(stats::rnorm(60L * 12L), 60L) + outer(stats::rnorm(60L), stats::rnorm(12L))
        estimate <- suppressMessages(estimateRuvFactorCount(y, rep(c("A", "B"), each = 6L), seed = 7L))
    }
    bus <- shiny::reactiveValues(state_manager = NULL, design_matrix = NULL, config_list = list(),
        tab_status = list(normalization = "pending"), state_update_trigger = NULL)
    suppressWarnings(suppressMessages(shiny::testServer(mod_omic_norm_server, args = list(id = "norm", workflow_data = bus,
        experiment_paths = list(), omic_type = "lipidomics", experiment_label = "kit129",
        config = omicGuiConfig("lipidomics")$norm), {
        omicNormSetRuvDesignResult(session$returned, omicNormRuvDesignResult("groups_all_features", "num_sv_be", 7L, "group",
            list(pos = list(k = 1L, applied = TRUE, estimate = estimate))))
        session$flushReact()
        expect_match(as.character(output$ruv_design_summary$html), "1 (estimated number)", fixed = TRUE)
        if (!is.null(estimate)) expect_true(nzchar(output$ruv_k_estimate_assay1$src))
    })))
})

test_that("(fourth pass item 3) refusals reach the tab's error notice with the kernel's message; DIA's None is shown as done", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    captured <- new.env(parent = emptyenv())
    refusal <- "the design technical_replicates needs declared technical replicates; this design has none (no_duplicates)"
    testthat::local_mocked_bindings(runDIAReferenceRuvCovariates = function(settings, ...) {
        captured$settings <- settings
        rlang::abort(refusal, class = "workflow_execution_invalid")
    }, updateProtNormRuvAuditTrail = function(...) invisible(NULL))
    workflow <- list2env(list(state_manager = list(workflow_type = "DIA_LIMPA")))
    norm_data <- list2env(list(correlation_filtering_complete = TRUE))
    notes <- character()
    observe <- function(input) runProtNormRuvCovariatesObserver(input, list2env(list()), workflow, norm_data, list(), "proteomics", "x",
        function(reference) NULL, function() "group", renderTextFn = function(text) text,
        showNotificationFn = function(message, ...) notes <<- c(notes, message), messageFn = function(...) NULL)
    observe(c(ruvDesignInput(ruv_design = "technical_replicates"), list(color_variable = NULL, shape_variable = NULL)))
    expect_identical(captured$settings[c("ruv_mode", "ruv_control_rule", "ruv_k_rule", "ruv_design")], list(ruv_mode = "covariates",
        ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be", ruv_design = "technical_replicates"))
    expect_identical(notes, paste("Error computing the RUV-III covariates:", refusal))
    notes <- character()
    observe(c(ruvDesignInput(ruv_design = "none"), list(color_variable = NULL, shape_variable = NULL)))
    expect_identical(norm_data$ruv_design_result$design, "none")
    expect_true(norm_data$ruv_covariates_complete)
    # classic proteomics: the kernel's refusal is the normalisation's error handler's
    errors <- list()
    testthat::local_mocked_bindings(runDIAReferenceNormalization = function(settings, ...) {
        captured$classic <- settings
        rlang::abort(refusal, class = "workflow_execution_invalid")
    }, recordProtNormMethod = function(...) invisible(NULL))
    runProtNormNormalizationWorkflow(ruvDesignInput(ruv_design = "technical_replicates"), list2env(list(state_manager = list(workflow_type = "LFQ"))),
        norm_data, list(), "proteomics", "x", function(...) NULL, function(...) NULL, function() "group",
        handleErrorFn = function(error) errors[[length(errors) + 1L]] <<- conditionMessage(error))
    expect_identical(captured$classic[c("ruv_design", "ruv_mode", "ruv_k_rule")], list(ruv_design = "technical_replicates",
        ruv_mode = "automatic", ruv_k_rule = "num_sv_be"))
    expect_identical(errors, list(refusal))
})

# KIT-129 sixth pass: a recorded check of the given outcome (M35, M36 fields).
ruvDesignCheck <- function(outcome, value = 0.21, threshold = 0.058077853) {
    limit <- outcome %in% c("detectable", "not_detectable")
    list(ruv_group_association = value, ruv_group_association_threshold = if (limit) threshold else NA_real_,
        ruv_group_association_flag = if (limit) identical(outcome, "detectable") else NA,
        ruv_group_association_outcome = outcome, ruv_group_association_source = "a source",
        ruv_group_association_basis = if (limit) "estimated_number" else NA_character_,
        # classic proteomics, the recommended design (arm 7) and arm 1 at 0.6 (ruvGroupAssociationConsequences)
        ruv_group_association_fdp_0 = if (limit) 0.0535118 else NA_real_, ruv_group_association_fdp_03 = if (limit) 0.0560489 else NA_real_,
        ruv_group_association_fdp_06 = if (limit) 0.0616284 else NA_real_, ruv_group_association_uncorrected_06 = if (limit) 0.593446 else NA_real_)
}

# the detectable text of `ruvDesignCheck("detectable")`, value 0.210
ruvDesignDetectableText <- function(texts, limit) sprintf(texts[["association_detectable"]], "0.210", limit, "5.4%", "5.6%", "6.2%", "59.3%")

test_that("(fourth pass item 4; sixth pass M36) the panel shows the check: detectable, not detectable, no limit, not available", {
    design <- omicGuiConfig("proteomics")$norm$ruv_design
    presets <- omicNormRuvPresets(design)
    texts <- design$tables[[1L]]$texts
    html <- function(association, route = NA_character_) as.character(omicNormRuvResultTags(omicNormRuvDesignResult("groups_all_features",
        "num_sv_be", 7L, "group", list(all = list(k = 1L, applied = TRUE, association = association)), route = route), presets, texts))
    detectable <- html(ruvDesignCheck("detectable"))
    expect_match(detectable, htmltools::htmlEscape(ruvDesignDetectableText(texts, "0.058")), fixed = TRUE)
    expect_match(detectable, 'role="alert"', fixed = TRUE)
    expect_match(html(ruvDesignCheck("not_detectable", 0.01)), sprintf(texts[["association_not_detectable"]], "0.010", "0.058"), fixed = TRUE)
    plain <- html(ruvDesignCheck("no_limit"))
    expect_match(plain, sprintf(texts[["association_no_limit"]], "0.210"), fixed = TRUE)
    expect_match(html(ruvDesignCheck("no_limit"), "dia_limpa"), sprintf(texts[["association_dia_limpa"]], "0.210"), fixed = TRUE)
    # M36: no "Warning" and no "false discoveries can rise above the target" in any of the check's texts
    for (key in grep("^association_", names(texts), value = TRUE)) {
        expect_false(grepl("Warning|warning|can rise above the target", texts[[key]]), info = key)
    }
    # a cell not measured for the approach reads "not measured", never another figure
    unmeasured <- html(utils::modifyList(ruvDesignCheck("detectable"), list(ruv_group_association_fdp_06 = NA_real_)))
    expect_match(unmeasured, "5.6% at 0.3 and not measured at 0.6", fixed = TRUE)
    # M35: a typed k equal to the estimated number says that the estimated number's limit applies
    typed <- html(utils::modifyList(ruvDesignCheck("not_detectable", 0.01), list(ruv_group_association_basis = "typed_k_equals_estimated_number")))
    expect_match(typed, texts[["association_typed_k_equals"]], fixed = TRUE)
    legacy <- as.character(omicNormRuvResultTags(omicNormRuvLegacyResult(), presets, texts))
    expect_match(legacy, htmltools::htmlEscape(presets$legacy$label), fixed = TRUE)
    expect_match(legacy, texts[["association_legacy"]], fixed = TRUE)
    # the texts with a number to fill take exactly that many
    keys <- c("association_detectable", "association_not_detectable", "association_no_limit", "association_dia_limpa")
    expect_identical(lengths(regmatches(texts[keys], gregexpr("%s", texts[keys], fixed = TRUE))),
        c(association_detectable = 6L, association_not_detectable = 2L, association_no_limit = 1L, association_dia_limpa = 1L))
})

test_that("(M42) every approach with a detection limit shows its measured figures at 0.6, never \"not measured\"", {
    texts <- omicGuiConfig("proteomics")$norm$ruv_design$tables[[1L]]$texts
    group <- rep(c("A", "B"), each = 6L)
    W <- cbind(W1 = (group == "B") + seq(0, 0.11, by = 0.01))
    thresholds <- ruvGroupAssociationThresholds()
    for (i in seq_len(nrow(thresholds))) {
        limit <- thresholds[i, ]
        check <- ruvGroupAssociationCheck(W, group, limit$omic, limit$ruv_design, limit$ruv_k_rule)
        expect_identical(check$ruv_group_association_outcome, "detectable", info = paste(limit$omic, limit$ruv_design))
        words <- omicNormRuvAssociationText(check, texts)
        expect_false(grepl(texts[["association_not_measured"]], words, fixed = TRUE), info = paste(limit$omic, limit$ruv_design))
    }
    # classic proteomics, the recommended design: 6.3% at 0.6, 59.3% without correction (the shipped No RUV)
    check <- ruvGroupAssociationCheck(W, group, "proteomics", "groups_all_features", "num_sv_be")
    expect_identical(omicNormRuvAssociationText(check, texts), sprintf(texts[["association_detectable"]],
        sprintf("%.3f", check$ruv_group_association), "0.058", "5.5%", "5.8%", "6.3%", "59.3%"))
})

# ------------------------------------------------------------------------------------------------
# KIT-129 fifth pass item 5: the warning's other two places (ticket section 3 (2) and (3)).
# ------------------------------------------------------------------------------------------------

test_that("(fifth pass item 5) the RUV QC tab's column holds the panel's words for the recorded check", {
    texts <- omicNormRuvRunTable("lipidomics")$texts
    two <- omicNormRuvDesignResult("groups_all_features", "num_sv_be", 7L, "group", list(
        pos = list(k = 1L, applied = TRUE, association = ruvDesignCheck("detectable")),
        neg = list(k = 0L, applied = FALSE, association = ruvDesignCheck("no_factors", NA))))
    rows <- omicNormRuvAssociationRows(two, texts)
    expect_identical(names(rows), c(texts[["result_assay"]], texts[["association_column"]]))
    expect_identical(rows[[1L]], "pos")
    expect_identical(rows[[2L]], ruvDesignDetectableText(texts, "0.058"))
    # the panel shows the same words
    shown <- as.character(omicNormRuvResultTags(two, omicNormRuvPresets(omicGuiConfig("lipidomics")$norm$ruv_design), texts))
    expect_match(shown, htmltools::htmlEscape(paste("Assay pos -", ruvDesignDetectableText(texts, "0.058"))), fixed = TRUE)
    one <- omicNormRuvAssociationRows(omicNormRuvDesignResult("groups_all_features", "user", 7L, "group",
        list(all = list(k = 1L, applied = TRUE, association = ruvDesignCheck("no_limit")))), texts)
    expect_identical(names(one), texts[["association_column"]])
    expect_identical(one[[1L]], sprintf(texts[["association_no_limit"]], "0.210"))
    expect_identical(nrow(omicNormRuvAssociationRows(NULL, texts)), 0L)
    expect_identical(nrow(omicNormRuvAssociationRows(omicNormRuvLegacyResult(), texts)), 0L)
    # the slot is in the RUV QC tab with the gate on only
    expect_match(ruvDesignHtml("lipidomics", TRUE), 'id="norm-ruv_design_association_table"', fixed = TRUE)
    expect_false(grepl("ruv_design_association_table", ruvDesignHtml("lipidomics", FALSE), fixed = TRUE))
})

test_that("(fifth pass item 5) the DA tab shows the warning above the run button when the run that feeds it was flagged", {
    texts <- omicGuiRuvTexts("features")
    steps <- list()
    testthat::local_mocked_bindings(buildWorkflowAnalysisRecord = function(manager, now = Sys.time()) list(steps = steps),
        .package = "MultiScholaR")
    manager <- structure(list(), class = "DIAStateAdapter")
    single <- function(outcome) c(list(rows = 10L), ruvDesignCheck(outcome))
    steps <- list(list(operation = "normalization-ruv", outcomes = single("detectable")))
    expect_identical(omicDaRuvAssociationLines(manager), character())   # the gate off: nothing
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    expect_identical(omicDaRuvAssociationLines(manager), ruvDesignDetectableText(texts, "0.058"))
    # the last RUV step is the one that feeds the test
    steps <- list(list(operation = "normalization-ruv", outcomes = single("detectable")),
        list(operation = "normalization-ruv", outcomes = single("not_detectable")))
    expect_identical(omicDaRuvAssociationLines(manager), character())
    steps <- list(list(operation = "trx_ruvs", outcomes = single("no_limit")))
    expect_identical(omicDaRuvAssociationLines(manager), character())
    # lipid/metab: per assay
    per <- function(field, values) stats::setNames(list(values), paste0(field, "_per_assay"))
    steps <- list(list(operation = "lipid_ruv", outcomes = c(per("ruv_group_association", list(pos = 0.21, neg = 0.01)),
        per("ruv_group_association_threshold", list(pos = 0.088611623, neg = 0.088611623)),
        per("ruv_group_association_flag", list(pos = TRUE, neg = FALSE)),
        per("ruv_group_association_outcome", list(pos = "detectable", neg = "not_detectable")),
        per("ruv_group_association_source", list(pos = "s", neg = "s")))))
    expect_identical(omicDaRuvAssociationLines(manager), paste("Assay pos -", sprintf(texts[["association_detectable"]], "0.210", "0.089",
        "not measured", "not measured", "not measured", "not measured")))
    expect_identical(omicDaRuvAssociationLines(NULL), character())
    tags <- as.character(omicDaRuvAssociationTags(omicDaRuvAssociationLines(manager)))
    expect_match(tags, 'role="alert"', fixed = TRUE)
    expect_null(omicDaRuvAssociationTags(character()))
    # the slot sits above the run button with the gate on, in both layouts; nothing with the gate off
    for (omic in c("lipidomics", "proteomics")) {
        html <- function(enabled) withr::with_options(list(MultiScholaR.ruv_design_panel = enabled),
            ruvDesignMaskTabs(as.character(htmltools::renderTags(mod_omic_da_ui("da", omicGuiConfig(omic)$da))$html)))
        on <- html(TRUE)
        expect_lt(regexpr('id="da-ruv_association_warning"', on, fixed = TRUE), regexpr('id="da-run_da_analysis"', on, fixed = TRUE),
            label = omic)
        expect_gt(regexpr('id="da-ruv_association_warning"', on, fixed = TRUE), 0L, label = omic)
        expect_identical(html(NULL), gsub('<div id="da-ruv_association_warning" class="shiny-html-output"></div>\\s*', "", on), info = omic)
    }
})
