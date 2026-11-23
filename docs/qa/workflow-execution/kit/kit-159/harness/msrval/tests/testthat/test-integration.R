# Part B (KIT-159 7.10): valRunIntegration drives the exported workbook functions in a user's order,
# records every setting, departures as forced settings with their reasons, tool warnings verbatim,
# the maths-library pin, and calls the post-integration hook. Functions are passed in so the order
# and records are tested without MOFA2 (the lane's lifecycle G9 runs the real ones).

fakeWorkbook <- function(log) {
    step <- function(name, value) {
        force(value)
        function(...) {
            log$calls <- c(log$calls, name)
            value
        }
    }
    list(
        integrationSettings = step("settings", list(
            model = data.frame(setting = c("factors", "seed"), default = c("15", "42")),
            results = data.frame(setting = "top_n", default = "20"))),
        integrationFindInputs = step("find", list(paths = "p", scan = "s")),
        integrationOpenProject = step("open", list(project = TRUE)),
        integrationMatchSamples = step("match", list(preview = 1, object = 2)),
        integrationConfirmMatch = step("confirm", list(parameters = 1, samples = 2)),
        integrationTrainModel = function(project, confirmed, input_matrix = NULL,
                                         overrides = list()) {
            log$calls <- c(log$calls, "train")
            log$overrides <- overrides
            warning("MOFA2: factor 3 explains little variance")
            list(value = 1, reference = 2, parameters = 3)
        },
        integrationGetResults = step("results", list(factors = data.frame(f = 1))),
        integrationReplayProject = step("replay", list(states = "s")),
        integrationPinMathsLibrary = step("pin", list(library = "openblas", threads = 1L))
    )
}

test_that("the workbook functions run in a user's order with every setting recorded", {
    log <- new.env()
    log$calls <- character()
    spec <- list(project_dir = "proj", sample_match_file = "match.tsv",
                 forced = list(factors = list(value = 10, reason = "small fixture")))
    res <- valRunIntegration(spec, fns = fakeWorkbook(log))
    expect_identical(log$calls, c("settings", "pin", "find", "open", "match", "confirm", "train",
                                  "results", "replay"))
    expect_identical(log$overrides, list(factors = 10))
    expect_named(res$settings, c("table", "setting", "default", "value", "source", "reason"))
    f <- res$settings[res$settings$setting == "factors", ]
    expect_identical(f$value, "10")
    expect_identical(f$source, "forced")
    expect_identical(f$reason, "small fixture")
    expect_identical(res$settings$source[res$settings$setting == "seed"], "default")
    expect_identical(res$warnings, "MOFA2: factor 3 explains little variance")
    expect_identical(res$pin, list(library = "openblas", threads = 1L))
})

test_that("a departure from a setting that does not exist is refused; the hook gets the results", {
    log <- new.env()
    log$calls <- character()
    bad <- list(project_dir = "proj", forced = list(nonsense = list(value = 1, reason = "r")))
    expect_error(valRunIntegration(bad, fns = fakeWorkbook(log)), class = "val_spec_invalid")
    no_reason <- list(project_dir = "proj", forced = list(factors = list(value = 1)))
    expect_error(valRunIntegration(no_reason, fns = fakeWorkbook(log)), class = "val_spec_invalid")
    got <- new.env()
    hook <- function(project, results) {
        got$results <- results
        got$test_mode <- getOption("multischolar.test_mode")
        "hooked"
    }
    res <- valRunIntegration(list(project_dir = "proj"), fns = fakeWorkbook(log),
                             post_integration = hook)
    expect_identical(got$results, list(factors = data.frame(f = 1)))
    expect_false(got$test_mode)
    expect_identical(res$hook, "hooked")
})
