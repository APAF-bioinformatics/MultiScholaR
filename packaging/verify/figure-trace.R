# KIT-134 / KIT-052 (Windows figure paths): record every figure a run asks a graphics device to write,
# every file copy, and every check of expected figures, with the path's length and whether the file is
# there afterwards. Build-time verification, not package code: it traces functions of the loaded
# namespaces (grDevices::png/pdf/svg/dev.off, ragg::agg_png (what ggplot2::ggsave draws with), base::file.copy,
# MultiScholaR's checkWorkflowFigureFiles and setupDirectories) and changes no result.
#
# Usage (from another verification script, after library(MultiScholaR)):
#   source(".../figure-trace.R"); kit134FigureTraceInstall(); .kit134FigureTrace$item <- "label"; ...
#   kit134FigureTraceReport(<output tsv>, <case label>)
#
# "exists" is read from the folder's listing (the folder is short), so it is right even where the
# file's own path is too long for the file functions; "bytes" is the file's size, NA where the
# size cannot be read at that path.

.kit134FigureTrace <- new.env()
local({
    e <- .kit134FigureTrace
    e$rows <- list()
    e$roots <- character()
    e$open <- list()
    e$pending <- NULL
    e$closing <- NULL
    e$item <- ""
    e$absolute <- function(path) {
        path <- as.character(path)[[1L]]
        if (grepl("^([A-Za-z]:|/|\\\\)", path)) path else file.path(getwd(), path)
    }
    e$listed <- function(path) {
        listing <- tryCatch(list.files(dirname(path), all.files = TRUE, no.. = TRUE), error = function(error) character())
        basename(path) %in% listing
    }
    e$add <- function(op, path, note = "") {
        if (!length(path) || is.na(path[[1L]]) || !nzchar(path[[1L]])) return(invisible())
        path <- gsub("\\\\", "/", e$absolute(path))
        bytes <- suppressWarnings(file.size(path))
        e$rows[[length(e$rows) + 1L]] <- data.frame(item = e$item, op = op, length = nchar(path),
            exists = e$listed(path), bytes = bytes, note = note, path = path, stringsAsFactors = FALSE)
        invisible()
    }
    e$deviceOpened <- function(op, path) e$pending <- list(op = op, path = path, before = grDevices::dev.list())
    # Runs when the device function returns or fails: a device that did not open is recorded at once.
    e$deviceReady <- function() {
        if (!is.null(e$pending)) {
            if (identical(grDevices::dev.list(), e$pending$before)) {
                e$add(e$pending$op, e$pending$path, "the device did not open")
            } else {
                e$open[[as.character(grDevices::dev.cur())]] <- e$pending
            }
        }
        e$pending <- NULL
    }
    e$deviceClosed <- function() {
        key <- as.character(e$closing)
        if (length(key) == 1L && !is.null(e$open[[key]])) {
            e$add(e$open[[key]]$op, e$open[[key]]$path, "at dev.off")
            e$open[[key]] <- NULL
        }
        e$closing <- NULL
    }
    e$copied <- function(from, to, result) {
        if (!length(from)) return(invisible())
        targets <- if (length(to) == 1L && dir.exists(to)) file.path(to, basename(from)) else rep_len(to, length(from))
        for (i in seq_along(from)) e$add("file.copy", targets[[i]], sprintf("returned %s", if (length(result) >= i) result[[i]] else NA))
    }
    e$checked <- function(files, message) for (file in files) e$add("checkWorkflowFigureFiles", file, message)
})

kit134FigureTraceInstall <- function() {
    # trace() rebuilds its own call, so its arguments are named here, not passed through `...`.
    quietly <- function(what, where, tracer = NULL, exit = NULL) {
        if (is.null(tracer)) {
            invisible(suppressMessages(trace(what, exit = exit, where = where, print = FALSE)))
        } else if (is.null(exit)) {
            invisible(suppressMessages(trace(what, tracer = tracer, where = where, print = FALSE)))
        } else {
            invisible(suppressMessages(trace(what, tracer = tracer, exit = exit, where = where, print = FALSE)))
        }
    }
    quietly("png", tracer = quote(.kit134FigureTrace$deviceOpened(paste0("grDevices::png type=", paste(type, collapse = "|")), filename)),
        exit = quote(.kit134FigureTrace$deviceReady()), where = asNamespace("grDevices"))
    quietly("pdf", tracer = quote(.kit134FigureTrace$deviceOpened("grDevices::pdf", file)),
        exit = quote(.kit134FigureTrace$deviceReady()), where = asNamespace("grDevices"))
    quietly("svg", tracer = quote(.kit134FigureTrace$deviceOpened("grDevices::svg", filename)),
        exit = quote(.kit134FigureTrace$deviceReady()), where = asNamespace("grDevices"))
    if (requireNamespace("ragg", quietly = TRUE)) {
        quietly("agg_png", tracer = quote(.kit134FigureTrace$deviceOpened("ragg::agg_png", filename)),
            exit = quote(.kit134FigureTrace$deviceReady()), where = asNamespace("ragg"))
    }
    quietly("dev.off", tracer = quote(.kit134FigureTrace$closing <- which), exit = quote(.kit134FigureTrace$deviceClosed()),
        where = asNamespace("grDevices"))
    quietly("file.copy", exit = quote(.kit134FigureTrace$copied(from, to, returnValue(logical()))), where = baseenv())
    namespace <- asNamespace("MultiScholaR")
    if (exists("checkWorkflowFigureFiles", envir = namespace, inherits = FALSE)) {
        quietly("checkWorkflowFigureFiles", tracer = quote(.kit134FigureTrace$checked(files, message)), where = namespace)
    }
    quietly("setupDirectories", tracer = quote(.kit134FigureTrace$roots <- unique(c(.kit134FigureTrace$roots,
        gsub("\\\\", "/", .kit134FigureTrace$absolute(base_dir))))), where = namespace)
    invisible(TRUE)
}

#' Write the recorded rows (with each path's project root and its length) and print them: every device,
#' ggsave and figure check in full; copies only where the target is 200 characters or longer or the copy
#' or the file is missing, with a count of the rest.
kit134FigureTraceReport <- function(path, case) {
    e <- .kit134FigureTrace
    rows <- if (length(e$rows)) do.call(rbind, e$rows) else data.frame(item = character(), op = character(), length = integer(),
        exists = logical(), bytes = numeric(), note = character(), path = character())
    temp <- gsub("\\\\", "/", tempdir())
    rootOf <- function(p) {
        known <- e$roots[!is.na(e$roots) & nzchar(e$roots)]
        roots <- if (length(known)) known[startsWith(p, paste0(known, "/"))] else character()
        if (length(roots)) return(roots[[which.max(nchar(roots))]])
        store <- regmatches(p, regexpr("^.*?(?=/\\.workflow-entities/)", p, perl = TRUE))
        if (length(store)) return(store)
        if (startsWith(p, paste0(temp, "/"))) return(paste0(temp, "/", sub("/.*$", "", substring(p, nchar(temp) + 2L))))
        NA_character_
    }
    rows$root <- vapply(rows$path, rootOf, character(1), USE.NAMES = FALSE)
    rows$root_length <- nchar(rows$root)
    rows <- cbind(case = rep(case, nrow(rows)), rows)
    utils::write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
    copies <- rows$op == "file.copy"
    shown <- !copies | rows$length >= 200L | !rows$exists | !grepl("TRUE", rows$note)
    cat(sprintf("FIGURES %s: tempdir %s (%d characters); project roots seen: %s\n", case, temp, nchar(temp),
        if (length(e$roots)) paste(sprintf("%s (%d)", e$roots, nchar(e$roots)), collapse = ", ") else "none"))
    cat(sprintf("FIGURES %s: %d rows (%d device or ggsave, %d checks, %d copies; %d copies under 200 characters and present not shown)\n",
        case, nrow(rows), sum(rows$op %in% c("ggsave") | grepl("::", rows$op)), sum(rows$op == "checkWorkflowFigureFiles"),
        sum(copies), sum(copies & !shown)))
    if (any(shown)) {
        out <- rows[shown, c("item", "op", "length", "root_length", "exists", "bytes", "note", "path")]
        cat(paste("FIGROW", apply(out, 1L, function(r) paste(trimws(r), collapse = "\t"))), sep = "\n")
    }
    for (item in unique(rows$item)) {
        r <- rows[rows$item == item, ]
        figures <- r[r$op != "file.copy", ]
        cat(sprintf("FIGSUM %s %s: longest figure path %s (root %s), longest copy target %s; figures missing %d of %d; copies missing %d of %d\n",
            case, item, if (nrow(figures)) max(figures$length) else NA, if (nrow(figures)) figures$root_length[which.max(figures$length)] else NA,
            if (any(r$op == "file.copy")) max(r$length[r$op == "file.copy"]) else NA, sum(!figures$exists), nrow(figures),
            sum(!r$exists[r$op == "file.copy"]), sum(r$op == "file.copy")))
    }
    invisible(rows)
}
