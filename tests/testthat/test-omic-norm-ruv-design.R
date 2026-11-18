# KIT-129 W2: the RUV design presets of the configuration (`norm$ruv_design`), their schema, and the shared panel
# section of the normalisation tab. Presentation and configuration only: no option is wired to a computation, and the
# section is not shown until the gate (`options(MultiScholaR.ruv_design_panel = TRUE)`) is on. The expected values
# below are written from the ticket (sections 1a to 1d and 2) and the design review of 2026-11-18 (sections 3 to 5),
# independently of the implementation's table.

ruvDesignOmics <- c("proteomics", "lipidomics", "metabolomics", "transcriptomics")
ruvDesignIds <- c("groups_all_features", "technical_replicates", "supplied_controls", "first_pass_nulls", "legacy", "none")

# Each table's measured figures: false discovery proportion in the worst cell, and power (ticket 1a to 1d).
ruvDesignFigures <- list(
    classic = list(groups_all_features = c(0.066, 0.91), technical_replicates = c(0.054, 0.90),
        supplied_controls = c(0.053, 0.91), first_pass_nulls = c(0.078, 0.91), legacy = c(0.385, 0.92),
        none = c(0.052, 0.49)),
    dia_limpa = list(groups_all_features = c(0.055, 0.91), technical_replicates = c(0.055, 0.77),
        supplied_controls = c(0.056, 0.91), first_pass_nulls = c(0.149, 0.92), legacy = c(0.435, 0.90),
        legacy_correction = c(0.466, 0.91), none = c(0.191, 0.56)),
    transcriptomics = list(groups_all_features = c(0.055, 0.77), supplied_controls = c(0.056, 0.76),
        first_pass_nulls = c(0.064, 0.60), legacy = c(0.299, 0.73), none = c(0.044, 0.005)),
    feature = list(groups_all_features = c(0.072, 0.43), supplied_controls = c(0.068, 0.37),
        first_pass_nulls = c(0.239, 0.50), legacy = c(0.556, 0.42), none = c(0.059, 0.003)))

# The k rules each design offers and the one pre-selected (ticket section 2, the user's direction U1).
ruvDesignKRules <- list(
    groups_all_features = list(rules = c("num_sv_be", "user"), default = "num_sv_be", pending = NULL),
    technical_replicates = list(rules = c("num_sv_be", "user"), default = "num_sv_be", pending = NULL),
    supplied_controls = list(rules = c("num_sv_be", "sweep_bounded", "user"), default = "num_sv_be",
        pending = "sweep_bounded"),
    first_pass_nulls = list(rules = c("num_sv_be", "sweep_bounded", "user"), default = "sweep_bounded",
        pending = "num_sv_be"),
    legacy = list(rules = character(), default = NULL, pending = NULL),
    legacy_correction = list(rules = character(), default = NULL, pending = NULL),
    none = list(rules = character(), default = NULL, pending = NULL))

#' Every table of the four omics, named `<omic>.<table>`.
ruvDesignTables <- function() {
    tables <- lapply(stats::setNames(ruvDesignOmics, ruvDesignOmics), \(omic) omicGuiConfig(omic)$norm$ruv_design$tables)
    unlist(tables, recursive = FALSE)
}

#' The figures expected of one table.
ruvDesignExpectedFigures <- function(name) {
    key <- sub("^(lipidomics|metabolomics)\\..*$", "feature", name)
    key <- sub("^transcriptomics\\..*$", "transcriptomics", key)
    ruvDesignFigures[[sub("^proteomics\\.", "", key)]]
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
        # Until the switch commit today's design is the selected one; the recommended one is recorded for it.
        expect_identical(tables[[name]]$selected, "legacy", info = name)
        expect_identical(tables[[name]]$recommended, "groups_all_features", info = name)
        expect_identical(tables[[name]]$selected_with_technical_replicates,
            if (startsWith(name, "transcriptomics")) NULL else "technical_replicates", info = name)
        for (option in options) {
            for (field in c("label", "tooltip", "note", "source")) {
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
    }
    # Lipid and metab: the technical-replicate design is being re-measured. Transcriptomics: not available.
    for (name in c("lipidomics.features", "metabolomics.features")) {
        expect_identical(ruvDesignOptions(tables[[name]])$technical_replicates$status, "pending", info = name)
    }
    trx <- ruvDesignOptions(tables[["transcriptomics.genes"]])$technical_replicates
    expect_identical(trx$status, "unavailable")
    expect_match(trx$unavailable, "summed into one sample per biological unit", fixed = TRUE)
    statuses <- unlist(lapply(tables, \(table) vapply(table$options, `[[`, "", "status")))
    expect_identical(sum(statuses == "pending"), 2L)
    expect_identical(sum(statuses == "unavailable"), 1L)
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
            unlist(lapply(table$options, \(option) option[c("label", "tooltip", "note", "fdp_where", "power_where",
                "unavailable")]))))))
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
    for (field in c("id", "tier", "label", "tooltip", "note", "status", "source", "test_form", "k_rules")) {
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
        table <- design$tables[[1L]]
        html <- ruvDesignHtml(omic, TRUE)
        off <- ruvDesignHtml(omic, FALSE)
        # One radio group under a static id, one input per design in the table's order, today's design checked.
        radios <- regmatches(html, gregexpr('<input type="radio" name="norm-ruv_design"[^>]*>', html))[[1L]]
        expect_identical(sub('.*value="([^"]*)".*', "\\1", radios), vapply(table$options, `[[`, "", "id"), info = omic)
        expect_identical(sub('.*value="([^"]*)".*', "\\1", grep("checked", radios, value = TRUE)), table$selected, info = omic)
        expect_true(grepl(sprintf('data-testid="%s"', design$testid), html, fixed = TRUE), info = omic)
        # A design that is not available is shown disabled, with its reason; no other is disabled.
        unavailable <- vapply(table$options, \(option) identical(option$status, "unavailable"), NA)
        expect_identical(grepl("disabled", radios), unavailable, info = omic)
        for (option in table$options) {
            where <- paste(omic, option$id)
            # The tier, the label, the figure and the note are text in the panel; the tooltip is the option's title.
            for (text in c(design$tiers[[option$tier]], option$label, option$note, omicNormRuvFigureText(option))) {
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
        expect_true(all(c("ruv_design", "ruv_k_rule", "ruv_k_user", "ruv_design_summary") %in% ids), info = omic)
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
            sort(c("ruv_design", "ruv_k_rule", "ruv_k_user", "ruv_design_summary", slots)), info = omic)
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
        expect_true(header > regexpr("RUV design:", html, fixed = TRUE) &&
            header < regexpr('class="shiny-options-group"', html, fixed = TRUE), info = omic)
        for (option in design$tables[[1L]]$options) {
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
        # One panel per k rule, side by side.
        expect_length(gregexpr('ruv-design-k-rule"', html, fixed = TRUE)[[1L]], length(design$k_rules))
    }
})

test_that("the figure text is the table's numbers in words; a pending or unavailable cell shows no number", {
    expect_identical(omicNormRuvFigureText(list(status = "measured", fdp = 0.066,
        fdp_where = "worst design tested: 12 per group", power = 0.91)),
        "Measured: 6.6% false discoveries at a nominal 5% (worst design tested: 12 per group); power 0.91 in the standard design.")
    expect_identical(omicNormRuvFigureText(list(status = "measured", fdp = 0.044, fdp_where = "worst design tested: x",
        power = 0.005, power_where = "there")),
        "Measured: 4.4% false discoveries at a nominal 5% (worst design tested: x); power 0.005 there.")
    expect_identical(omicNormRuvFigureText(list(status = "pending")), "Not yet measured: the simulation is still running.")
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
    expect_identical(updates[[1L]]$selected, "sweep_bounded")
    expect_identical(unname(updates[[2L]]$choices), c("num_sv_be", "user"))
    expect_identical(updates[[2L]]$selected, "num_sv_be")
    expect_identical(names(updates[[2L]]$choices), c("Estimated count", "A number"))

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
