# KIT-245 product part, option (b) (the owner's decision of 2026-11-25, VAL-244; ticket section 3.4 and 10):
# "Push to GitHub" is not offered in the installed application, recognised by MSR_BUNDLE (set by the Windows shell,
# the Linux launcher and passed to systemd units; E-1). Declared before the run, the base value beside each case:
#   E1  with MSR_BUNDLE set, the Summary UI of the four omics creates no enable_github, github_org, github_email,
#       github_username, project_id or push_to_github input, keeps the heading "Version Control (Optional)" and shows
#       the sentence saying pushing is not available there (base: every input present, no sentence);
#   G1  with MSR_BUNDLE unset, the UI's HTML is identical to the one rendered by the same call with MSR_BUNDLE
#       explicitly unset (passes on the base too; the frozen digests of test-omic-summary-module.R guard the rest).

k245PanelUi <- function(omic) {
    module_id <- if (identical(omic, "proteomics")) "session_summary" else "summary"
    as.character(mod_omic_summary_ui(module_id, omicGuiConfig(omic)$summary))
}

k245PanelOmics <- c("proteomics", "lipidomics", "metabolomics", "transcriptomics")

test_that("KIT-245 E1: the installed application shows no GitHub push controls, and says why", {
    withr::local_envvar(MSR_BUNDLE = withr::local_tempdir("kit245-bundle-"))
    for (omic in k245PanelOmics) {
        module_id <- if (identical(omic, "proteomics")) "session_summary" else "summary"
        html <- k245PanelUi(omic)
        for (input in c("enable_github", "github_org", "github_email", "github_username", "project_id",
            "push_to_github")) {
            expect_false(grepl(paste0("id=\"", module_id, "-", input, "\""), html, fixed = TRUE),
                info = paste(omic, input))
        }
        expect_true(grepl("Version Control (Optional)", html, fixed = TRUE), info = omic)
        expect_true(grepl("Pushing to GitHub is not available in the installed application.", html, fixed = TRUE),
            info = omic)
    }
})

test_that("KIT-245 G1: outside the installed application the Summary UI is unchanged", {
    for (omic in k245PanelOmics) {
        unset <- withr::with_envvar(c(MSR_BUNDLE = NA), k245PanelUi(omic))
        module_id <- if (identical(omic, "proteomics")) "session_summary" else "summary"
        expect_true(grepl(paste0("id=\"", module_id, "-push_to_github\""), unset, fixed = TRUE), info = omic)
        expect_false(grepl("Pushing to GitHub is not available", unset, fixed = TRUE), info = omic)
    }
})

# Correction pass (AUDIT-KIT-245-impl-20261125, I04): the push observer itself refuses in the installed application.
#   E2  with MSR_BUNDLE set, a client that sends the five inputs and the press reaches no push (base: one push);
#   G2  without MSR_BUNDLE the same press reaches the push once (passes on the base too).
k245PanelPresses <- function() {
    pressed <- 0L
    testthat::local_mocked_bindings(omicSummaryGithub = function(...) {
        pressed <<- pressed + 1L
        invisible(TRUE)
    }, .package = "MultiScholaR")
    root <- withr::local_tempdir("kit245-panel-")
    dirs <- list(lipidomics = list(base_dir = root, source_dir = file.path(root, "scripts", "lipidomics")))
    shiny::testServer(mod_omic_summary_server, args = list(project_dirs = dirs, omic_type = "lipidomics",
        experiment_label = NULL, workflow_data = NULL, config = omicGuiConfig("lipidomics")$summary), {
        values <- session$getReturned()
        values$report_generated <- TRUE
        session$setInputs(enable_github = TRUE, github_org = "kit245-org", github_email = "kit245@example.invalid",
            github_username = "kit245", project_id = "kit245-panel")
        session$setInputs(push_to_github = 1)
    })
    pressed
}

test_that("KIT-245 E2: in the installed application the push observer refuses inputs a client sends itself", {
    withr::local_envvar(MSR_BUNDLE = withr::local_tempdir("kit245-bundle-"))
    expect_identical(k245PanelPresses(), 0L)
})

test_that("KIT-245 G2: outside the installed application the push observer still pushes", {
    withr::local_envvar(MSR_BUNDLE = NA)
    expect_identical(k245PanelPresses(), 1L)
})
