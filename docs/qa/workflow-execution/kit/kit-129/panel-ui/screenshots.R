# KIT-129 panel UI: screenshots of the normalisation tab's RUV design section, gate on, in headless Chromium.
# Usage: Rscript screenshots.R <worktree> <output dir> <before|after>
# The page is the app's own shell (app_ui(): shinydashboard, its styles and scripts, dark mode by default) with the
# home tab's content replaced by the real mod_omic_norm_ui() of one omic; the server runs only the section's own
# observers (registerOmicNormRuvDesign). Nothing else of the tab is served: plots and other outputs stay empty.
args <- commandArgs(TRUE)
worktree <- normalizePath(args[[1L]]); out_dir <- args[[2L]]; tag <- args[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(K129_WORKTREE = worktree)

harnessApp <- function() {
    options(MultiScholaR.ruv_design_panel = TRUE)
    pkgload::load_all(Sys.getenv("K129_WORKTREE"), quiet = TRUE, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
    omic <- Sys.getenv("K129_OMIC"); route <- Sys.getenv("K129_ROUTE")
    norm <- omicGuiConfig(omic)$norm
    ui <- function(request) {
        htmltools::tagQuery(app_ui(request))$find("#shiny-tab-home")$empty()$
            append(mod_omic_norm_ui("norm", norm))$allTags()
    }
    server <- function(input, output, session) {
        manager <- new.env(); manager$workflow_type <- if (nzchar(route)) route else NULL
        bus <- shiny::reactiveValues(state_manager = manager, state_update_trigger = NULL)
        shiny::moduleServer("norm", function(input, output, session) {
            registerOmicNormRuvDesign(input, session, bus, norm$ruv_design)
        })
    }
    shiny::shinyApp(ui, server)
}

Sys.setenv(CHROMOTE_CHROME = Sys.glob(path.expand(
    "~/.cache/ms-playwright/chromium_headless_shell-1223/*/chrome-headless-shell"))[[1L]])
cases <- list(c("proteomics", "", "prot-classic"), c("proteomics", "DIA_LIMPA", "prot-dialimpa"),
    c("transcriptomics", "", "trx"), c("lipidomics", "", "lipid"))
# design, k rule (NA: leave the pre-selected one), name
states <- list(c("legacy", NA, "1-legacy-default"), c("groups_all_features", NA, "2-recommended-estimated"),
    c("supplied_controls", "sweep_bounded", "3-supplied-separation-pending"),
    c("first_pass_nulls", "user", "4-firstpass-number"), c("none", NA, "5-none"))
section <- if (identical(tag, "before")) ".col-sm-3 > .well" else ".ruv-design"
for (case in cases) {
    Sys.setenv(K129_OMIC = case[[1L]], K129_ROUTE = case[[2L]])
    app <- shinytest2::AppDriver$new(harnessApp, name = case[[3L]], load_timeout = 180000, timeout = 30000,
        width = 1440, height = 900)
    cat("pid", case[[3L]], app$get_chromote_session()$parent$get_browser()$get_process()$get_pid(), "\n")
    shot <- function(width, name, selector = section) {
        file <- file.path(out_dir, sprintf("%s_%s_%d_%s.png", tag, case[[3L]], width, name))
        app$get_screenshot(file, selector = "scrollable_area")
        if (!identical(selector, "scrollable_area")) {
            # The section cut from the page's image by its box (a selector's own screenshot came out shifted).
            box <- app$get_js(sprintf("(function() { var r = document.querySelector('%s').getBoundingClientRect(); return [r.left + window.scrollX, r.top + window.scrollY, r.width, r.height]; })()", selector))
            box <- round(unlist(box))
            magick::image_write(magick::image_crop(magick::image_read(file),
                sprintf("%dx%d+%d+%d", box[[3L]], box[[4L]], box[[1L]], box[[2L]])), file)
        }
        cat("wrote", file, "\n")
    }
    for (state in states) {
        if (identical(state[[1L]], "none") && !identical(case[[3L]], "prot-classic")) next
        app$set_inputs(`norm-ruv_design` = state[[1L]], wait_ = FALSE); app$wait_for_idle(500)
        if (!is.na(state[[2L]])) { app$set_inputs(`norm-ruv_k_rule` = state[[2L]], wait_ = FALSE); app$wait_for_idle(500) }
        for (width in c(1440L, 1024L)) {
            app$set_window_size(width, 900); Sys.sleep(0.4)
            shot(width, state[[3L]])
            # The whole page for the alignment with the rest of the tab: three cases only.
            if ((identical(case[[3L]], "prot-classic") && identical(state[[3L]], "1-legacy-default")) ||
                (case[[3L]] %in% c("prot-classic", "trx") && identical(state[[3L]], "2-recommended-estimated"))) shot(width, paste0(state[[3L]], "_page"), "scrollable_area")
        }
    }
    # The k rules' explanation: closed in every state above; opened here from the keyboard (focus on its summary,
    # Enter), as a user without a mouse opens it.
    app$set_inputs(`norm-ruv_design` = "supplied_controls", wait_ = FALSE); app$wait_for_idle(500)
    if (!identical(tag, "before") && isTRUE(app$get_js("document.querySelector('.ruv-design-k-why') !== null"))) {
        app$run_js("document.querySelector('.ruv-design-k-why > summary').focus();")
        chrome <- app$get_chromote_session()
        for (type in c("keyDown", "keyUp")) chrome$Input$dispatchKeyEvent(type = type, key = "Enter", code = "Enter",
            windowsVirtualKeyCode = 13L, nativeVirtualKeyCode = 13L, text = if (identical(type, "keyDown")) "\r" else NULL)
        Sys.sleep(0.3)
        cat("disclosure", case[[3L]], "summary focusable:", app$get_js("document.activeElement.tagName"),
            "open after Enter:", app$get_js("document.querySelector('.ruv-design-k-why').open"), "\n")
        for (width in c(1440L, 1024L)) { app$set_window_size(width, 900); Sys.sleep(0.4); shot(width, "7-supplied-rules-open") }
        app$run_js("document.querySelector('.ruv-design-k-why').open = false;")
    }
    if (identical(case[[3L]], "prot-classic")) {
        app$set_inputs(`norm-ruv_design` = "groups_all_features", wait_ = FALSE); app$wait_for_idle(500)
        app$run_js("toggleDarkMode();"); Sys.sleep(0.6)
        for (width in c(1440L, 1024L)) { app$set_window_size(width, 900); Sys.sleep(0.4); shot(width, "6-recommended-light") }
    }
    app$stop()
}
