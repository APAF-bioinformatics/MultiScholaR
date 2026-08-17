# Shared browser-runtime contract for local E2E helpers and CI runners.

e2e_flag_enabled <- function(value) {
  isTRUE(tolower(trimws(as.character(value))) %in% c("1", "true", "yes", "y"))
}

e2e_browser_required <- function() {
  e2e_flag_enabled(Sys.getenv("MULTISCHOLAR_E2E_BROWSER_REQUIRED", unset = "false"))
}

e2e_resolve_executable <- function(candidate) {
  candidate <- trimws(as.character(candidate))
  if (!nzchar(candidate)) {
    return("")
  }

  expanded <- path.expand(candidate)
  if (file.exists(expanded)) {
    return(normalizePath(expanded, winslash = "/", mustWork = TRUE))
  }
  unname(Sys.which(candidate))
}

e2e_browser_executable_usable <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(FALSE)
  }
  if (!identical(.Platform$OS.type, "windows") && file.access(path, mode = 1L) != 0L) {
    return(FALSE)
  }

  version_output <- tryCatch(
    suppressWarnings(system2(path, "--version", stdout = TRUE, stderr = TRUE)),
    error = function(condition) structure(character(), status = 1L)
  )
  status <- attr(version_output, "status")
  is.null(status) || identical(as.integer(status), 0L)
}

e2e_browser_launch_probe <- function(path) {
  tryCatch(
    {
      chromote::with_chromote_chrome(path, {
        session <- chromote::ChromoteSession$new()
        on.exit(try(session$close(), silent = TRUE), add = TRUE)
      })
      list(usable = TRUE, error = NULL)
    },
    error = function(condition) list(
      usable = FALSE,
      error = conditionMessage(condition)
    )
  )
}

e2e_browser_preflight <- function(
    package_available = function(package) requireNamespace(package, quietly = TRUE),
    chrome_finder = function() chromote::find_chrome(),
    executable_check = e2e_browser_executable_usable,
    launch_probe = e2e_browser_launch_probe) {
  shinytest2_available <- isTRUE(package_available("shinytest2"))
  chromote_available <- isTRUE(package_available("chromote"))
  failures <- character()

  if (!shinytest2_available) {
    failures <- c(failures, "missing R package 'shinytest2'")
  }
  if (!chromote_available) {
    failures <- c(failures, "missing R package 'chromote'")
  }

  chrome_path <- ""
  chrome_error <- NULL
  if (chromote_available) {
    candidates <- tryCatch(
      suppressMessages(chrome_finder()),
      error = function(condition) {
        chrome_error <<- conditionMessage(condition)
        character()
      }
    )
    candidates <- unique(as.character(candidates))
    candidates <- candidates[!is.na(candidates) & nzchar(trimws(candidates))]
    for (candidate in candidates) {
      resolved <- e2e_resolve_executable(candidate)
      if (nzchar(resolved) && isTRUE(executable_check(resolved))) {
        launch_result <- launch_probe(resolved)
        if (isTRUE(launch_result$usable)) {
          chrome_path <- resolved
          break
        }
        chrome_error <- launch_result$error
      }
    }
  }

  if (chromote_available && !nzchar(chrome_path)) {
    detail <- if (is.null(chrome_error) || !nzchar(chrome_error)) {
      "Chromium executable was not found or did not pass its executable probe"
    } else {
      paste("Chromium launch probe failed:", chrome_error)
    }
    failures <- c(failures, detail)
  }

  list(
    available = length(failures) == 0L,
    shinytest2_available = shinytest2_available,
    chromote_available = chromote_available,
    chrome_path = if (nzchar(chrome_path)) chrome_path else NA_character_,
    failures = as.list(failures)
  )
}

e2e_browser_preflight_message <- function(preflight) {
  failures <- unlist(preflight$failures, use.names = FALSE)
  if (!length(failures)) {
    return("E2E browser runtime is available")
  }
  paste(failures, collapse = "; ")
}

e2e_require_browser_preflight <- function(preflight = e2e_browser_preflight()) {
  if (!isTRUE(preflight$available)) {
    stop(
      paste("Required E2E browser preflight failed:", e2e_browser_preflight_message(preflight)),
      call. = FALSE
    )
  }
  invisible(preflight)
}
