# KIT-193 target 2 (ticket section 6.2, 5.5): the shipped replay-rule selection against the Lean model.
# PENDING: moves to tests/testthat/ when KIT-187 and KIT-190 have merged (it calls workflowReplayRule,
# workflowFixedInForce, workflowReplayRuleOptions and reads workflowReplayRuleMarkers with KIT-190's row). Before it
# first runs, section 8.2: diff the merged code against the lines Spec/ReplayRule.lean cites. Still to write at merge:
# the replay-build.tsv comparison through buildWorkflowAnalysisRecord on a minimal store, built as KIT-187 T8 (iii)
# builds one (read that test at merge; not on staging today). Kind: package (bindings.tsv rows T2).
# Tables: tests/testdata/kit193/replay-{rule,inforce,options}.tsv (Oracle/ReplayRule.lean).

kit193ReplayTable <- function(name) {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", name), sep = "\t", quote = "",
        stringsAsFactors = FALSE, colClasses = "character", na.strings = character())
}

#' A record step from the table's text (`name|source|value;…`, `-` when empty)
kit193Step <- function(text) {
    if (identical(text, "-")) return(list(parameters = list()))
    entries <- strsplit(strsplit(text, ";", fixed = TRUE)[[1L]], "|", fixed = TRUE)
    list(parameters = lapply(entries, function(e) list(name = e[[1L]], source = e[[2L]],
        value = encodeWorkflowRecordValue(e[[3L]]), provenance = "fixed")))
}

#' Fixed constants from the table's text (`name|value;…`, `-` when empty)
kit193Fixed <- function(text) {
    if (identical(text, "-")) return(list())
    pairs <- strsplit(strsplit(text, ";", fixed = TRUE)[[1L]], "|", fixed = TRUE)
    stats::setNames(lapply(pairs, `[[`, 2L), vapply(pairs, `[[`, character(1), 1L))
}

kit193Declarations <- list(
    none = list(fixed = list(other = "x")),
    missing = list(fixed = list(missing_value_comparison = "counts", other = "x")),
    ruv = list(fixed = list(ruv_control_grouping = "kept", other = "x")),
    both = list(fixed = list(missing_value_comparison = "counts", ruv_control_grouping = "kept", other = "x")))

test_that("KIT-193 T2 pair: the marker table is the model's (two rows, distinct names, options and tokens)", {
    expect_identical(sort(names(workflowReplayRuleMarkers)), c("missing_value_comparison", "ruv_control_grouping"))
    expect_identical(workflowReplayRuleMarkers$missing_value_comparison,
        list(option = "MultiScholaR.missing_filter_rule", current = "counts", earlier = "rounded_percentage"))
    expect_identical(workflowReplayRuleMarkers$ruv_control_grouping,
        list(option = "MultiScholaR.ruv_control_grouping", current = "kept", earlier = "dropped"))
})

test_that("KIT-193 T2 pair: a token survives the record's encode-decode round trip", {
    for (token in c("counts", "rounded_percentage", "kept", "dropped", "other", "x")) {
        expect_identical(decodeWorkflowRecordValue(encodeWorkflowRecordValue(token)), token)
    }
})

test_that("KIT-193 T2 pair: workflowReplayRule equals the model on every listed step (replay-rule.tsv)", {
    table <- kit193ReplayTable("replay-rule.tsv")
    expect_identical(nrow(table), 1057L)
    for (i in seq_len(nrow(table))) {
        step <- kit193Step(table$step[[i]])
        expect_identical(workflowReplayRule(step, "missing_value_comparison"),
            table$rule_missing_value_comparison[[i]], info = table$step[[i]])
        expect_identical(workflowReplayRule(step, "ruv_control_grouping"),
            table$rule_ruv_control_grouping[[i]], info = table$step[[i]])
    }
})

test_that("KIT-193 T2 pair: workflowFixedInForce equals the model under every option state (replay-inforce.tsv)", {
    table <- kit193ReplayTable("replay-inforce.tsv")
    expect_identical(nrow(table), 128L)
    state <- function(marker, label) switch(label, unset = NULL, current = workflowReplayRuleMarkers[[marker]]$current,
        earlier = workflowReplayRuleMarkers[[marker]]$earlier, another = "another")
    for (i in seq_len(nrow(table))) {
        options <- list(state("missing_value_comparison", table$opt_missing_value_comparison[[i]]),
            state("ruv_control_grouping", table$opt_ruv_control_grouping[[i]]))
        names(options) <- c("MultiScholaR.missing_filter_rule", "MultiScholaR.ruv_control_grouping")
        got <- withr::with_options(options, workflowFixedInForce(kit193Fixed(table$fixed[[i]])))
        expect_identical(got %||% list(), kit193Fixed(table$in_force[[i]]), info = paste(table[i, 1:3], collapse = " "))
    }
})

test_that("KIT-193 T2 pair: workflowReplayRuleOptions equals the model on every step and declaration shape", {
    table <- kit193ReplayTable("replay-options.tsv")
    steps <- kit193ReplayTable("replay-rule.tsv")
    for (i in seq_len(nrow(table))) {
        step <- kit193Step(steps$step[[as.integer(table$step_id[[i]])]])
        got <- workflowReplayRuleOptions(step, kit193Declarations[[table$declaration[[i]]]])
        want <- if (identical(table$options[[i]], "-")) list() else {
            kv <- strsplit(strsplit(table$options[[i]], ";", fixed = TRUE)[[1L]], "=", fixed = TRUE)
            stats::setNames(lapply(kv, `[[`, 2L), vapply(kv, `[[`, character(1), 1L))
        }
        expect_identical(got[order(names(got))], want[order(names(want))], info = paste(table[i, ], collapse = " "))
    }
})

test_that("KIT-193 T2 pair, planted defect (ticket 5.5): a rule that takes source registry as current fails", {
    step <- kit193Step("missing_value_comparison|registry|counts")
    testthat::local_mocked_bindings(workflowReplayRule = function(step, marker) {
        spec <- workflowReplayRuleMarkers[[marker]]
        for (parameter in step$parameters) {
            if (identical(parameter$name, marker) && parameter$source %in% c("fixed", "registry")) return(spec$current)
        }
        spec$earlier
    })
    testthat::expect_failure(expect_identical(workflowReplayRule(step, "missing_value_comparison"), "rounded_percentage"))
})
