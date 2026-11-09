# KIT-108 W2, instrument I1: the controller trace by action (D7).
#
# Splits the harness's existing 50 ms trace (memory-trace.tsv, uptime clock) into the report step's windows using
# the driver's marks (driver-phases.jsonl, epoch clock) and the step marks (journey-phases.tsv, uptime clock):
#   pre      the 5 s before the report step starts
#   da_step, ruviii_export_step   the preceding worker steps, for reference (journey-phases.tsv)
#   save     report start to the copy's enter (open tab, label, description, the save click)
#   copy     enter/exit:action:<ns>summary-copy_to_publication
#   render   the copy's exit to report complete (the generate_report click)
#   post     report complete to the journey teardown mark
# Read-only. No package code is loaded.
#
# Clock alignment: epoch = uptime + offset. Every driver mark belongs to one journey step, so for the true offset
# o: o <= epoch(first mark of step k) - uptime(start k) and o >= epoch(last mark of step k) - uptime(complete k).
# The bracket [max lower, min upper] over all steps is the alignment interval; the midpoint is used and the
# half-width is reported as the alignment error. Marks are assigned to steps with a first estimate from the first
# step and refined once.
#
# Usage: Rscript analyse-traces.R <out.tsv> <member_dir> [<member_dir> ...]
#   member_dir is .../kit/runtime/kit-journey-<key>/candidate-1 (or /member)

args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
out <- args[[1L]]
dirs <- args[-1L]
GiB <- 1024^3

readMarks <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(lines)]
    rows <- lapply(lines, jsonlite::fromJSON)
    data.frame(epoch = vapply(rows, function(r) as.numeric(r$time), numeric(1)),
        stage = vapply(rows, function(r) as.character(r$stage), character(1)), stringsAsFactors = FALSE)
}

alignClocks <- function(steps, marks) {
    assign <- function(offset) {
        u <- marks$epoch - offset
        vapply(u, function(x) {
            inside <- which(steps$start <= x & steps$complete >= x)
            if (length(inside)) return(inside[[1L]])
            d <- pmin(abs(steps$start - x), abs(steps$complete - x))
            which.min(d)
        }, integer(1))
    }
    o <- marks$epoch[[1L]] - steps$start[[1L]]
    for (pass in 1:2) {
        k <- assign(o)
        lower <- upper <- numeric(0)
        for (s in unique(k)) {
            e <- marks$epoch[k == s]
            upper <- c(upper, min(e) - steps$start[[s]])
            lower <- c(lower, max(e) - steps$complete[[s]])
        }
        lo <- max(lower); hi <- min(upper)
        o <- (lo + hi) / 2
    }
    list(offset = o, lower = lo, upper = hi, error = (hi - lo) / 2, consistent = lo <= hi)
}

valueAt <- function(trace, column, t, side = c("before", "after")) {
    side <- match.arg(side)
    i <- if (side == "before") max(which(trace$uptime <= t)) else min(which(trace$uptime >= t))
    if (!is.finite(i)) return(NA_real_)
    trace[[column]][[i]]
}

windowRow <- function(member, corpus, name, from, to, trace) {
    rows <- trace$uptime >= from & trace$uptime <= to
    gap <- if (sum(rows) > 1L) max(diff(trace$uptime[rows])) else NA_real_
    data.frame(member = member, corpus = corpus, window = name, from = from, to = to, seconds = to - from,
        app_entry = valueAt(trace, "controller_anon", from, "before") / GiB,
        app_peak = if (any(rows)) max(trace$controller_anon[rows], na.rm = TRUE) / GiB else NA_real_,
        app_exit = valueAt(trace, "controller_anon", to, "before") / GiB,
        slice_entry = valueAt(trace, "anon", from, "before") / GiB,
        slice_peak = if (any(rows)) max(trace$anon[rows], na.rm = TRUE) / GiB else NA_real_,
        slice_exit = valueAt(trace, "anon", to, "before") / GiB,
        app_min = if (any(rows)) min(trace$controller_anon[rows], na.rm = TRUE) / GiB else NA_real_,
        slice_min = if (any(rows)) min(trace$anon[rows], na.rm = TRUE) / GiB else NA_real_,
        samples = sum(rows), max_gap = gap, stringsAsFactors = FALSE)
}

analyseMember <- function(dir) {
    member <- sub("^kit-journey-", "", basename(dirname(dir)))
    timed <- tryCatch(readRDS(file.path(dir, "timed-member.rds")), error = function(e) list())
    corpus <- if (!is.null(timed$corpus)) timed$corpus else NA_character_
    trace <- utils::read.delim(file.path(dir, "memory-trace.tsv"))
    for (col in c("uptime", "anon", "controller_anon")) trace[[col]] <- suppressWarnings(as.numeric(trace[[col]]))
    trace <- trace[is.finite(trace$uptime), , drop = FALSE]
    phases <- utils::read.delim(file.path(dir, "journey-phases.tsv"), stringsAsFactors = FALSE)
    starts <- phases[phases$event == "start", ]
    ends <- phases[phases$event == "complete", ]
    steps <- data.frame(step = starts$step, start = starts$uptime,
        complete = ends$uptime[match(starts$step, ends$step)], stringsAsFactors = FALSE)
    steps <- steps[is.finite(steps$complete), ]
    marks <- readMarks(file.path(dir, "driver-phases.jsonl"))
    clock <- alignClocks(steps, marks)
    # /proc/uptime has 10 ms resolution: a bracket inverted by less than that is consistent at resolution.
    if (!clock$consistent && clock$lower - clock$upper > 0.02)
        warning(member, ": the clock bracket is inverted by ", round(clock$lower - clock$upper, 3), " s")
    epochOf <- function(stage) {
        e <- marks$epoch[marks$stage == stage]
        if (length(e) != 1L) stop(member, ": expected one mark '", stage, "', found ", length(e))
        e
    }
    copyStage <- grep("summary-copy_to_publication$", grep("^enter:action:", marks$stage, value = TRUE), value = TRUE)
    ns <- sub("summary-copy_to_publication$", "", sub("^enter:action:", "", copyStage))
    copyInE <- epochOf(paste0("enter:action:", ns, "summary-copy_to_publication"))
    copyOutE <- epochOf(paste0("exit:action:", ns, "summary-copy_to_publication"))
    report <- steps[steps$step == "report", ]
    teardown <- phases$uptime[phases$event == "teardown"]
    build <- function(offset) {
        copyIn <- copyInE - offset; copyOut <- copyOutE - offset
        w <- list(windowRow(member, corpus, "pre", report$start - 5, report$start, trace))
        for (s in intersect(c("da", "ruviii_export"), steps$step)) {
            st <- steps[steps$step == s, ]
            w <- c(w, list(windowRow(member, corpus, paste0(s, "_step"), st$start, st$complete, trace)))
        }
        w <- c(w, list(windowRow(member, corpus, "report_step", report$start, report$complete, trace),
            windowRow(member, corpus, "save", report$start, copyIn, trace),
            windowRow(member, corpus, "copy", copyIn, copyOut, trace),
            windowRow(member, corpus, "render", copyOut, report$complete, trace)))
        if (length(teardown) == 1L)
            w <- c(w, list(windowRow(member, corpus, "post", report$complete, teardown, trace),
                windowRow(member, corpus, "teardown_5s", teardown - 5, teardown, trace)))
        do.call(rbind, w)
    }
    mid <- build(clock$offset)
    # The bracket's ends: the copy's enter and exit move by up to the alignment error; report the app's anon at a
    # window's entry and exit under both ends, so a conclusion that depends on the alignment shows it.
    early <- build(clock$upper); late <- build(clock$lower)
    rng <- function(a, b) sprintf("%.3f-%.3f", pmin(a, b) / 1, pmax(a, b))
    mid$app_entry_range <- rng(early$app_entry, late$app_entry)
    mid$app_exit_range <- rng(early$app_exit, late$app_exit)
    mid$seconds_range <- rng(early$seconds, late$seconds)
    mid$app_entry_early <- early$app_entry; mid$app_entry_late <- late$app_entry
    mid$app_exit_early <- early$app_exit; mid$app_exit_late <- late$app_exit
    mid$offset_epoch_minus_uptime <- sprintf("%.3f", clock$offset)
    mid$align_error_s <- max(clock$error, 0.01)
    mid
}

all <- do.call(rbind, lapply(dirs, analyseMember))
num <- vapply(all, is.numeric, logical(1))
all[num] <- lapply(all[num], function(x) round(x, 3))
utils::write.table(all, out, sep = "\t", quote = FALSE, row.names = FALSE)
# The per-member summary: the application's anon at each boundary of the report step (GiB). A boundary set by a
# driver mark carries its range over the clock bracket ("early" = the copy marks at their earliest uptime).
summary <- do.call(rbind, lapply(split(all, all$member), function(m) {
    w <- function(name) m[m$window == name, ]
    pre <- w("pre"); save <- w("save"); copy <- w("copy"); render <- w("render"); td <- w("teardown_5s")
    rep <- w("report_step"); da <- w("da_step")
    rise <- function(x) c(x$app_exit_early - x$app_entry_early, x$app_exit_late - x$app_entry_late, x$app_exit - x$app_entry)
    data.frame(member = m$member[[1]], corpus = m$corpus[[1]], align_error_s = m$align_error_s[[1]],
        da_peak = if (nrow(da)) da$app_peak else NA, before_save = pre$app_exit,
        after_save = save$app_exit, after_save_range = save$app_exit_range,
        after_copy = copy$app_exit, after_copy_range = copy$app_exit_range,
        after_render = render$app_exit, report_peak = rep$app_peak,
        teardown = if (nrow(td)) td$app_exit else NA, teardown_min5s = if (nrow(td)) td$app_min else NA,
        slice_resident_end = if (nrow(td)) td$slice_min else NA,
        save_s = save$seconds, copy_s = copy$seconds, render_s = render$seconds, render_s_range = render$seconds_range,
        save_rise_min = min(rise(save)), save_rise_max = max(rise(save)),
        copy_rise_min = min(rise(copy)), copy_rise_max = max(rise(copy)),
        render_rise_min = min(rise(render)), render_rise_max = max(rise(render)),
        report_rise = rep$app_exit - rep$app_entry, stringsAsFactors = FALSE)
}))
num <- vapply(summary, is.numeric, logical(1))
summary[num] <- lapply(summary[num], function(x) round(x, 3))
utils::write.table(summary, sub("[.]tsv$", "-summary.tsv", out), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250)
print(summary, row.names = FALSE)
print(all[, c("member", "window", "seconds", "seconds_range", "app_entry", "app_peak", "app_exit", "app_exit_range", "slice_peak",
    "max_gap", "align_error_s")], row.names = FALSE)
