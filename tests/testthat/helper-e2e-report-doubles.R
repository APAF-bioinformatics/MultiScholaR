# E2E test doubles for report template downloading.
# Replaces downloadReportTemplate() with a deterministic stub that avoids
# live GitHub fetches during CI/E2E runs.

.e2e_report_template_dir <- function() {
    testthat::test_path("..", "testdata", "e2e", "report_templates")
}

.e2e_stub_template_path <- function(rmd_filename) {
    template_dir <- .e2e_report_template_dir()
    stub_path <- file.path(template_dir, rmd_filename)

    if (!file.exists(stub_path)) {
        available <- if (dir.exists(template_dir)) {
            paste(
                list.files(template_dir, pattern = "\\.rmd$", ignore.case = TRUE)
                , collapse = ", "
            )
        } else {
            "(template dir not found)"
        }
        rlang::abort(paste0(
            "E2E report template stub not found: "
            , stub_path
            , "\n"
            , "Available stubs: "
            , available
        ))
    }

    stub_path
}

#' Install a deterministic double for downloadReportTemplate()
#'
#' Replaces `downloadReportTemplate()` in the MultiScholaR namespace with a
#' function that returns a path to the matching stub template from
#' `tests/testdata/e2e/report_templates/`. No network access occurs.
#'
#' Automatically restored when `.local_envir` exits (withr::defer).
#'
#' @param .local_envir Environment to scope the teardown into (default: caller frame).
install_report_template_double <- function(.local_envir = parent.frame()) {
    double_fn <- function(omic_type, rmd_filename) {
        .e2e_stub_template_path(rmd_filename)
    }

    pkg_env <- asNamespace("MultiScholaR")
    binding_name <- "downloadReportTemplate"

    had_binding <- exists(binding_name, envir = pkg_env, inherits = FALSE)
    old_value <- if (had_binding) get(binding_name, envir = pkg_env, inherits = FALSE) else NULL
    was_locked <- had_binding && bindingIsLocked(binding_name, pkg_env)

    if (was_locked) unlockBinding(binding_name, pkg_env)
    assign(binding_name, double_fn, envir = pkg_env)
    if (was_locked) lockBinding(binding_name, pkg_env)

    withr::defer({
        env <- asNamespace("MultiScholaR")
        if (exists(binding_name, envir = env, inherits = FALSE) &&
                bindingIsLocked(binding_name, env)) {
            unlockBinding(binding_name, env)
        }
        if (had_binding) {
            assign(binding_name, old_value, envir = env)
        } else if (exists(binding_name, envir = env, inherits = FALSE)) {
            rm(list = binding_name, envir = env)
        }
        if (was_locked && exists(binding_name, envir = env, inherits = FALSE)) {
            lockBinding(binding_name, env)
        }
    }, envir = .local_envir)

    invisible(double_fn)
}
