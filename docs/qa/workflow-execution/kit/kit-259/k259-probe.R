# KIT-259 diagnostic (triage of KIT-247 G0-4 run 45604291093): shared condition capture and the figure renderer.
# Sourced by g0-diag.R (the plain R process) and, inside the application's own process, by the probe app that
# g0-diag.R starts through shinytest2 (part 3). No package code; observation only.

k259Text <- function(value, n = 3000L) {
    if (is.null(value)) return("<null>")
    value <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(value)), collapse = " | "))
    if (!nzchar(value)) "<none>" else substr(value, 1L, n)
}

#' The condition chain of `cnd`: class, message and call of the condition and of every `parent`, outermost first.
k259Chain <- function(cnd) {
    out <- character()
    depth <- 0L
    while (!is.null(cnd) && inherits(cnd, "condition") && depth < 10L) {
        call <- tryCatch(paste(deparse(conditionCall(cnd), nlines = 2L), collapse = " "), error = function(e) "<call?>")
        out <- c(out, sprintf("[%d] <%s> %s {call: %s}", depth, paste(class(cnd), collapse = "/"),
            conditionMessage(cnd), call))
        cnd <- cnd$parent
        depth <- depth + 1L
    }
    paste(out, collapse = " ==> ")
}

#' The rlang traces stored on the conditions of the chain (rlang errors carry one), formatted.
k259StoredTraces <- function(cnd) {
    out <- character()
    depth <- 0L
    while (!is.null(cnd) && inherits(cnd, "condition") && depth < 10L) {
        if (!is.null(cnd$trace)) out <- c(out, sprintf("[%d] %s", depth, tryCatch(paste(format(cnd$trace,
            simplify = "none"), collapse = " / "), error = function(e) paste("<trace format failed>", conditionMessage(e)))))
        cnd <- cnd$parent
        depth <- depth + 1L
    }
    if (length(out)) paste(out, collapse = " ### ") else "<none>"
}

#' Evaluate `expr`; return status, the full chain, the warnings seen (not muffled) and the backtrace at the first
#' error a handler outside the expression sees.
k259Capture <- function(expr) {
    warnings <- character()
    trace <- NULL
    result <- tryCatch(
        withCallingHandlers({
            force(expr)
            list(status = "ok")
        }, warning = function(w) {
            warnings <<- c(warnings, paste0("<", paste(class(w), collapse = "/"), "> ", conditionMessage(w)))
        }, error = function(e) {
            if (is.null(trace)) trace <<- tryCatch(paste(format(rlang::trace_back(), simplify = "none"), collapse = " / "),
                error = function(x) paste("<trace failed>", conditionMessage(x)))
        }),
        error = function(e) list(status = "error", error = e))
    list(status = result$status,
        chain = if (identical(result$status, "error")) k259Chain(result$error) else "<none>",
        stored_traces = if (identical(result$status, "error")) k259StoredTraces(result$error) else "<none>",
        backtrace = trace %||% "<none>",
        warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else "<none>",
        error = result$error)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' The device ggsave picks for a .png file in this process (deparsed), if the internal can be read.
k259GgsaveDevice <- function() {
    ns <- asNamespace("ggplot2")
    for (name in c("validate_device", "plot_dev")) {
        if (exists(name, envir = ns, inherits = FALSE)) {
            dev <- tryCatch(get(name, envir = ns)("png", "x.png", 150), error = function(e) NULL)
            if (!is.null(dev)) return(paste0(name, ": ", k259Text(deparse(dev), 400L)))
        }
    }
    "<not read>"
}

k259Environment <- function() {
    version <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "<absent>")
    list(R = R.version.string, ggplot2 = version("ggplot2"), ragg = version("ragg"), scales = version("scales"),
        rlang = version("rlang"), cli = version("cli"), shiny = version("shiny"), gtable = version("gtable"),
        systemfonts = version("systemfonts"), textshaping = version("textshaping"),
        warn = as.character(getOption("warn")), shiny_running = as.character(tryCatch(shiny::isRunning(),
            error = function(e) NA)), display = Sys.getenv("DISPLAY", "<unset>"),
        capabilities = paste(names(capabilities()), capabilities(), sep = "=", collapse = ","),
        ggsave_png_device = k259GgsaveDevice(), device_option = k259Text(deparse(getOption("device")), 200L),
        current_device = paste(names(grDevices::dev.cur()), grDevices::dev.cur()),
        global_handlers = tryCatch(k259Text(names(globalCallingHandlers()), 300L), error = function(e) "<not read>"))
}

#' Render each plot RDS with ggsave (as saveProtNormQcPlotArtifact does: width 10, height 8, dpi 150) and write one
#' TSV row per plot plus the environment, to `<out>/p3-<where>.tsv`. Used inside the application's process.
k259ProbeRun <- function(where, plots, out) {
    rows <- list()
    env <- k259Environment()
    for (name in names(env)) rows[[length(rows) + 1L]] <- data.frame(key = paste0("p3.", where, ".env.", name),
        value = k259Text(env[[name]]), stringsAsFactors = FALSE)
    for (path in plots) {
        plot <- readRDS(path)
        target <- tempfile(fileext = ".png")
        got <- k259Capture(ggplot2::ggsave(target, plot, width = 10, height = 8, dpi = 150))
        id <- sub("[.]rds$", "", basename(path))
        rows[[length(rows) + 1L]] <- data.frame(key = paste0("p3.", where, ".", id, ".status"),
            value = if (identical(got$status, "ok")) "ok" else paste("error:", k259Text(got$chain)), stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(key = paste0("p3.", where, ".", id, ".warnings"), value = k259Text(got$warnings),
            stringsAsFactors = FALSE)
        rows[[length(rows) + 1L]] <- data.frame(key = paste0("p3.", where, ".", id, ".backtrace"),
            value = k259Text(got$backtrace), stringsAsFactors = FALSE)
    }
    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out, paste0("p3-", where, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    invisible(table)
}
