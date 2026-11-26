# KIT-247 correction pass (R2 45612382879: test-workflow-native-session.R's KIT-026 test regressed). A re-save is a design
# saved over an EARLIER saved design: with no earlier design on the bus there is nothing to compare and the save is
# committed as before. A real change over later states with no session to ask cannot be confirmed: it is refused with a
# classed condition, nothing committed, the later states left as they are (coordinator's ruling, ticket r1.4 section 19).

kit247Manager <- function(history = list("initial", "lipid_raw_data_s4", "lipid_normalized")) {
    manager <- structure(new.env(parent = emptyenv()), class = "DIAStateAdapter")
    manager$.context <- new.env(parent = emptyenv())
    manager$.context$snapshot <- list(state_history = history)
    manager
}

kit247Design <- function(group = c("x", "x", "y", "y")) {
    data.frame(Run = c("A", "B", "C", "D"), group = group, replicates = c(1L, 2L, 1L, 2L), stringsAsFactors = FALSE)
}

test_that("KIT-247 N1: with no earlier design on the bus there is nothing to compare: 'first' (committed as before)", {
    manager <- kit247Manager()
    expect_identical(compareOmicsDesignResult(manager, list(design_matrix = NULL, data = NULL), "lipid_raw_data_s4",
        list(design_matrix = NULL, data = NULL)), "first")
    expect_identical(compareOmicsDesignResult(manager, list(design_matrix = kit247Design(), data = NULL), "lipid_raw_data_s4",
        list(design_matrix = NULL, data = NULL)), "first")
})

test_that("KIT-247 N2: a real change over later states with no session is refused with a classed condition, nothing committed", {
    committed <- FALSE
    expect_error(runOmicsDesignResaveGate(NULL, "changed", TRUE, commit = function(cut) committed <<- TRUE,
        unchanged = function() NULL, showModal = function(...) NULL), class = "omics_design_resave_unconfirmed")
    expect_false(committed)
    condition <- tryCatch(runOmicsDesignResaveGate(NULL, "changed", TRUE, commit = function(cut) NULL,
        unchanged = function() NULL, showModal = function(...) NULL), error = function(e) e)
    expect_match(conditionMessage(condition), "differs from the saved design", fixed = TRUE)
})

test_that("KIT-247 N3: a real change over later states with a session still asks first (unchanged behaviour)", {
    session <- shiny::MockShinySession$new()
    asked <- 0L
    committed <- FALSE
    shiny::isolate(runOmicsDesignResaveGate(session, "changed", TRUE, commit = function(cut) committed <<- TRUE,
        unchanged = function() NULL, showModal = function(...) asked <<- asked + 1L))
    expect_identical(asked, 1L)
    expect_false(committed)
})
