# KIT-156 smoke (run-k156-rc instrumentation, not a package test): the test setting as the packaged entry applies it.
# The block of packaging/electron/entry.R that reads MULTISCHOLAR_TEST_FEATURES is evaluated here, in the bundle's R
# against its installed MultiScholaR, and the pages the two options switch on are rendered with the options that block
# set. The smoke job sets MULTISCHOLAR_TEST_FEATURES=1 for every file it runs.

k156EntryBlock <- function() {
    dir <- normalizePath(getwd(), winslash = "/")
    repeat {
        entry <- file.path(dir, "packaging", "electron", "entry.R")
        if (file.exists(entry)) break
        parent <- dirname(dir)
        if (identical(parent, dir)) stop("packaging/electron/entry.R not found above ", getwd())
        dir <- parent
    }
    found <- Filter(\(e) is.call(e) && identical(e[[1L]], as.name("if")) &&
        grepl("MULTISCHOLAR_TEST_FEATURES", paste(deparse(e[[2L]]), collapse = ""), fixed = TRUE),
        as.list(parse(entry, keep.source = FALSE)))
    stopifnot(length(found) == 1L)
    found[[1L]]
}

# The entry's block with the variable and the shell's probe as given; the options it sets stay set (callers reset them).
k156Apply <- function(value, probe) {
    withr::local_envvar(MULTISCHOLAR_TEST_FEATURES = value, MSR_LAUNCHER_PROBE = probe)
    output <- utils::capture.output(eval(k156EntryBlock(), envir = new.env(parent = globalenv())))
    list(output = output, probe = Sys.getenv("MSR_LAUNCHER_PROBE"),
        ruv = getOption("MultiScholaR.ruv_design_panel"), integration = getOption("MultiScholaR.integration_gui"))
}

k156Gates <- function(env = parent.frame()) {
    withr::local_options(list(MultiScholaR.ruv_design_panel = NULL, MultiScholaR.integration_gui = NULL), .local_envir = env)
}

k156Omics <- c("proteomics", "lipidomics", "metabolomics", "transcriptomics")
k156Html <- function(tag) as.character(htmltools::renderTags(tag)$html)
systemd <- '{"launcher":"systemd","reason":"every requirement holds"}'
process <- '{"launcher":"process","reason":"kernel 6.8 is older than 6.12"}'

test_that("KIT-156: the smoke runs with the test setting", {
    expect_identical(Sys.getenv("MULTISCHOLAR_TEST_FEATURES"), "1")
})

test_that("KIT-156: without the setting the entry's block changes nothing", {
    k156Gates()
    for (value in c("", "0", "true", "yes")) {
        applied <- k156Apply(value, systemd)
        expect_identical(applied$output, character(), info = value)
        expect_null(applied$ruv)
        expect_null(applied$integration)
        expect_identical(applied$probe, systemd, info = value)
    }
})

test_that("KIT-156: with the setting both options are TRUE, said on one line, and the systemd choice becomes the plain launcher", {
    k156Gates()
    applied <- k156Apply("1", systemd)
    expect_true(applied$ruv)
    expect_true(applied$integration)
    expect_length(applied$output, 1L)
    expect_match(applied$output, "^MSR-TEST-FEATURES ")
    line <- jsonlite::fromJSON(sub("^MSR-TEST-FEATURES ", "", applied$output))
    expect_identical(line, list(variable = "MULTISCHOLAR_TEST_FEATURES", ruv_design_panel = TRUE, integration_gui = TRUE))
    probe <- jsonlite::fromJSON(applied$probe)
    expect_identical(probe$launcher, "process")
    expect_match(probe$reason, "MULTISCHOLAR_TEST_FEATURES", fixed = TRUE)
    expect_match(probe$reason, "every requirement holds", fixed = TRUE)
    k156Gates()
    expect_identical(k156Apply("1", process)$probe, process)
})

test_that("KIT-156: with the options the setting sets, the Integration tile is live and every omic shows the RUV design panel", {
    k156Gates()
    off_page <- k156Html(app_ui(NULL))
    off_norm <- vapply(k156Omics, \(omic) k156Html(mod_omic_norm_ui("norm", omicGuiConfig(omic)$norm)), "")
    k156Apply("1", process)
    on_page <- k156Html(app_ui(NULL))
    expect_match(off_page, '<div id="integration_box" class="omic-selection-box disabled"', fixed = TRUE)
    expect_match(on_page, '<div id="integration_box" class="omic-selection-box" onclick=', fixed = TRUE)
    for (omic in k156Omics) {
        on_norm <- k156Html(mod_omic_norm_ui("norm", omicGuiConfig(omic)$norm))
        expect_match(on_norm, 'id="norm-ruv_design"', fixed = TRUE, info = omic)
        expect_false(grepl('id="norm-ruv_design"', off_norm[[omic]], fixed = TRUE), info = omic)
    }
    ui <- k156Html(mod_integration_ui("integration_workflow"))
    for (id in c("integration-import-commit", "integration-alignment-confirm", "integration-model-run",
        "integration-results-save", "integration-record-save")) {
        expect_match(ui, sprintf('data-testid="%s"', id), fixed = TRUE)
    }
})
