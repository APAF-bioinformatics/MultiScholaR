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
    # KIT-052 second confirmation: the final name of a published file (publishWorkflowActionOutputs renames its
    # `.publish-` copy into place) and every classed refusal of unwritten figures, with its message.
    e$renamed <- function(from, to, result) {
        for (i in seq_along(to)) e$add("file.rename", to[[i]], sprintf("returned %s", if (length(result) >= i) result[[i]] else NA))
    }
    e$raised <- list()
    e$abort <- function(message, class) {
        if (!"workflow_figure_unwritten" %in% class) return(invisible())
        labels <- names(message)
        if (is.null(labels)) labels <- rep("", length(message))
        text <- paste(ifelse(nzchar(labels), paste0("[", labels, "] "), ""), message, collapse = " | ")
        e$raised[[length(e$raised) + 1L]] <- data.frame(item = e$item, class = paste(class, collapse = "/"),
            message = cli::ansi_strip(text), stringsAsFactors = FALSE)
        invisible()
    }
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
    quietly("file.rename", exit = quote(.kit134FigureTrace$renamed(from, to, returnValue(logical()))), where = baseenv())
    quietly("abort", tracer = quote(.kit134FigureTrace$abort(message, class)), where = asNamespace("rlang"))
    namespace <- asNamespace("MultiScholaR")
    if (exists("checkWorkflowFigureFiles", envir = namespace, inherits = FALSE)) {
        quietly("checkWorkflowFigureFiles", tracer = quote(.kit134FigureTrace$checked(files, message)), where = namespace)
    }
    quietly("setupDirectories", tracer = quote(.kit134FigureTrace$roots <- unique(c(.kit134FigureTrace$roots,
        gsub("\\\\", "/", .kit134FigureTrace$absolute(base_dir))))), where = namespace)
    # KIT-052 Summary confirmation (767ad7d6): every notification the Summary tab's save, copy and report give the user
    # (SUMMARY lines, as shown), and every rmarkdown render with the file it wrote (RENDER lines).
    if (exists("omicSummaryNotify", envir = namespace, inherits = FALSE)) {
        quietly("omicSummaryNotify", tracer = quote(cat(sprintf("SUMMARY %s %s: %s\n", .kit134FigureTrace$item, type,
            substr(gsub("[\r\n\t]+", " ", paste(message, collapse = " ")), 1L, 1200L)))), where = namespace)
    }
    if (requireNamespace("rmarkdown", quietly = TRUE)) {
        quietly("render", exit = quote(try({
            written <- returnValue(NULL)
            if (is.character(written) && length(written) == 1L && file.exists(written)) {
                cat(sprintf("RENDER %s: %s, %.0f bytes, %d characters\n", .kit134FigureTrace$item, basename(written), file.size(written),
                    nchar(written)))
            } else cat(sprintf("RENDER %s: returned no file\n", .kit134FigureTrace$item))
        }, silent = TRUE)), where = asNamespace("rmarkdown"))
    }
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
        figures <- r[!r$op %in% c("file.copy", "file.rename"), ]
        cat(sprintf("FIGSUM %s %s: longest figure path %s (root %s), longest copy target %s; figures missing %d of %d; copies missing %d of %d\n",
            case, item, if (nrow(figures)) max(figures$length) else NA, if (nrow(figures)) figures$root_length[which.max(figures$length)] else NA,
            if (any(r$op == "file.copy")) max(r$length[r$op == "file.copy"]) else NA, sum(!figures$exists), nrow(figures),
            sum(!r$exists[r$op == "file.copy"]), sum(r$op == "file.copy")))
    }
    kit134FigurePairs(rows, case, sub("\\.tsv$", "-pairs.tsv", path))
    for (x in e$raised) cat(sprintf("RAISED %s %s: class %s: %s\n", case, x$item, x$class, substr(gsub("[\r\n\t]+", " ", x$message), 1L, 900L)))
    invisible(rows)
}

#' Each figure a device drew, paired with where it ended: the copies and renames after it of a file with the same
#' name, up to the next drawing of that name in the same item, outside the short scratch folders (`wfx-<hex>` under
#' tempdir()). A figure drawn in place is its own final file. "ok": listed in its folder and not empty (a size that
#' cannot be read counts as listed only).
kit134FigurePairs <- function(rows, case, path) {
    device <- grepl("::", rows$op)
    temp <- gsub("\\\\", "/", tempdir())
    scratch <- startsWith(rows$path, paste0(temp, "/wfx-"))
    moved <- rows$op %in% c("file.copy", "file.rename") & !scratch
    name <- basename(rows$path)
    ok <- function(i) rows$exists[i] & (is.na(rows$bytes[i]) | rows$bytes[i] > 0)
    pairs <- list()
    for (i in which(device)) {
        later <- which(seq_len(nrow(rows)) > i & rows$item == rows$item[[i]] & name == name[[i]])
        stop_at <- later[device[later]]
        if (length(stop_at)) later <- later[later < min(stop_at)]
        finals <- later[moved[later]]
        if (!length(finals)) finals <- i
        longest <- finals[which.max(rows$length[finals])]
        native <- regmatches(rows$path[[longest]], regexpr("/\\.dia-native-[0-9a-f]+/.*$", rows$path[[longest]]))
        pairs[[length(pairs) + 1L]] <- data.frame(case = case, item = rows$item[[i]], device = rows$op[[i]],
            drawn_length = rows$length[[i]], drawn_ok = ok(i), drawn_in_scratch = scratch[[i]], finals = length(finals),
            final_length = rows$length[[longest]], final_root_length = rows$root_length[[longest]], final_ok = all(ok(finals)),
            after_native = if (length(native)) nchar(sub("^/\\.dia-native-[0-9a-f]+", "", native)) else NA_integer_,
            drawn_path = rows$path[[i]], final_path = rows$path[[longest]], stringsAsFactors = FALSE)
    }
    if (!length(pairs)) return(invisible())
    pairs <- do.call(rbind, pairs)
    utils::write.table(pairs, path, sep = "\t", quote = FALSE, row.names = FALSE)
    short <- function(p) paste(utils::tail(strsplit(p, "/", fixed = TRUE)[[1L]], 3L), collapse = "/")
    for (j in seq_len(nrow(pairs))) {
        p <- pairs[j, ]
        cat(sprintf("FIGPAIR %s\t%s\t%s\tdrawn %d %s%s\tfinal %d (root %s, after .dia-native %s) %s x%d\t%s\n", case, p$item, p$device,
            p$drawn_length, if (p$drawn_ok) "ok" else "MISSING", if (p$drawn_in_scratch) " scratch" else "", p$final_length,
            p$final_root_length, p$after_native, if (p$final_ok) "ok" else "MISSING", p$finals, short(p$final_path)))
    }
    for (item in unique(pairs$item)) {
        p <- pairs[pairs$item == item, ]
        cat(sprintf(paste("FIGPAIRSUM %s %s: figures drawn %d (in scratch %d; PNG %d, PDF %d), drawn ok %d, final ok %d;",
            "longest drawn %d, longest final %d (root %s); longest after .dia-native %s\n"),
            case, item, nrow(p), sum(p$drawn_in_scratch), sum(grepl("\\.png$", p$final_path, ignore.case = TRUE)),
            sum(grepl("\\.pdf$", p$final_path, ignore.case = TRUE)), sum(p$drawn_ok), sum(p$final_ok), max(p$drawn_length),
            max(p$final_length), p$final_root_length[which.max(p$final_length)],
            if (all(is.na(p$after_native))) NA else max(p$after_native, na.rm = TRUE)))
    }
    invisible(pairs)
}
