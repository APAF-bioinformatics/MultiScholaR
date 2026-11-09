# KIT-108 W2: reads I2's per-process rows (sampler.sh) against a GUI member's harness files (D7, M1/M2).
#
# 1. The lifecycle checks (D7 M1): the tick gaps (overall and inside the report step), the controller's R PID, and
#    the R PID's RssAnon against the harness's 50 ms controller_anon (memory-trace.tsv, nearest sample by uptime) and
#    against the controller cgroup's anon on the same tick.
# 2. The report step by window, on I2's own clocks: the sampler's rows carry epoch and uptime, so the driver's epoch
#    marks map without the I1 bracket (the offset is the median of epoch - uptime over the rows; its spread is shown).
#    Per window: the R PID's RssAnon, the controller cgroup's anon, and the peak of every non-R process in the
#    controller's cgroup (a pandoc child shows here), plus the other units (monitor, workers) for context.
#
# Usage: Rscript analyse-sampler.R <sampler.tsv> <member_dir> <out-prefix> [<controller unit regex>]
#   member_dir = .../kit/runtime/kit-journey-<key>/candidate-1; regex default "-controller\\.service$"

args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
samplerPath <- args[[1L]]; dir <- args[[2L]]; prefix <- args[[3L]]
controllerRegex <- if (length(args) >= 4L) args[[4L]] else "-controller\\.service$"
GiB <- 1024^3

s <- utils::read.delim(samplerPath, stringsAsFactors = FALSE, colClasses = "character")
num <- c("epoch", "uptime", "pid", "vm_rss_kb", "rss_anon_kb", "rss_file_kb", "rss_shmem_kb", "vm_hwm_kb", "smaps_rss_kb",
    "smaps_pss_anon_kb", "smaps_anonymous_kb", "cg_anon", "cg_file", "cg_shmem")
for (n in num) s[[n]] <- suppressWarnings(as.numeric(s[[n]]))
events <- s[!s$event %in% c("process", "cgroup"), c("epoch", "uptime", "event", "unit", "note")]
proc <- s[s$event == "process", ]
cg <- s[s$event == "cgroup", ]
ctl <- unique(cg$unit[grepl(controllerRegex, cg$unit)])
if (length(ctl) != 1L) stop("expected one controller unit matching ", controllerRegex, ", found: ", paste(ctl, collapse = ","))
cgc <- cg[cg$unit == ctl, ]
pc <- proc[proc$unit == ctl, ]
rp <- pc[pc$comm == "R", ]
rPids <- unique(rp$pid)
offset <- stats::median(s$epoch - s$uptime, na.rm = TRUE)
spread <- diff(range(s$epoch - s$uptime, na.rm = TRUE))

trace <- utils::read.delim(file.path(dir, "memory-trace.tsv"))
for (col in c("uptime", "anon", "controller_anon")) trace[[col]] <- suppressWarnings(as.numeric(trace[[col]]))
trace <- trace[is.finite(trace$uptime), ]
phases <- utils::read.delim(file.path(dir, "journey-phases.tsv"), stringsAsFactors = FALSE)
report <- list(start = phases$uptime[phases$step == "report" & phases$event == "start"],
    complete = phases$uptime[phases$step == "report" & phases$event == "complete"])
teardown <- phases$uptime[phases$event == "teardown"]
marks <- do.call(rbind, lapply(readLines(file.path(dir, "driver-phases.jsonl"), warn = FALSE), function(l) {
    r <- jsonlite::fromJSON(l); data.frame(epoch = as.numeric(r$time), stage = r$stage, stringsAsFactors = FALSE)
}))
markUptime <- function(pattern) {
    e <- marks$epoch[grepl(pattern, marks$stage)]
    if (length(e) != 1L) stop("expected one mark matching ", pattern, ", found ", length(e))
    e - offset
}

# 1. Lifecycle checks.
ticks <- sort(unique(cgc$uptime))
gapAll <- if (length(ticks) > 1L) max(diff(ticks)) else NA_real_
inReport <- ticks[ticks >= report$start & ticks <= report$complete]
gapReport <- if (length(inReport) > 1L) max(diff(c(report$start, inReport, report$complete))) else NA_real_
nearest <- function(u) {
    i <- findInterval(u, trace$uptime)
    i[i < 1L] <- 1L
    j <- pmin(i + 1L, nrow(trace))
    ifelse(abs(trace$uptime[j] - u) < abs(trace$uptime[i] - u), j, i)
}
rTick <- stats::aggregate(rss_anon_kb ~ uptime, data = rp, FUN = sum)
k <- nearest(rTick$uptime)
rTick$harness_controller_anon <- trace$controller_anon[k]
rTick$harness_dt <- abs(trace$uptime[k] - rTick$uptime)
# The same tick's cgroup row: written just before the tick's process rows (correction 1, 2026-11-09: the first
# version interpolated between ticks, which showed up to 15% at a step change that was not in the data).
cgc <- cgc[order(cgc$uptime), ]
rTick$cg_anon <- cgc$cg_anon[pmax(findInterval(rTick$uptime, cgc$uptime), 1L)]
rTick$rel_err_harness <- (rTick$rss_anon_kb * 1024 - rTick$harness_controller_anon) / rTick$harness_controller_anon
rTick$rel_err_cgroup <- (rTick$rss_anon_kb * 1024 - rTick$cg_anon) / rTick$cg_anon
check <- data.frame(
    item = c("sampler rows", "controller unit", "controller R PIDs", "ticks (controller cgroup rows)", "max tick gap (s)",
        "max tick gap inside the report step (s)", "epoch - uptime offset (median)", "offset spread over rows (s)",
        "R RssAnon vs harness controller_anon: max |relative error|", "rows within 10% of controller_anon",
        "R RssAnon vs controller cgroup anon: max |relative error|", "max |dt| to the nearest harness sample (s)",
        "events"),
    value = c(nrow(s), ctl, paste(rPids, collapse = ","), length(ticks), round(gapAll, 3), round(gapReport, 3),
        sprintf("%.3f", offset), round(spread, 3), round(max(abs(rTick$rel_err_harness), na.rm = TRUE), 4),
        sprintf("%d / %d", sum(abs(rTick$rel_err_harness) <= 0.10, na.rm = TRUE), nrow(rTick)),
        round(max(abs(rTick$rel_err_cgroup), na.rm = TRUE), 4), round(max(rTick$harness_dt), 3),
        paste(events$event, collapse = ",")), stringsAsFactors = FALSE)
utils::write.table(check, paste0(prefix, "-check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(rTick, paste0(prefix, "-tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 2. Windows on the sampler's clocks.
copyIn <- markUptime("^enter:action:.*summary-copy_to_publication$")
copyOut <- markUptime("^exit:action:.*summary-copy_to_publication$")
windows <- list(pre = c(report$start - 5, report$start), save = c(report$start, copyIn), copy = c(copyIn, copyOut),
    render = c(copyOut, report$complete))
if (length(teardown) == 1L) windows$post <- c(report$complete, teardown)
at <- function(x, u, value) { i <- which(x$uptime <= u); if (length(i)) value[[max(i)]] else NA_real_ }
rows <- lapply(names(windows), function(w) {
    from <- windows[[w]][[1L]]; to <- windows[[w]][[2L]]
    r <- rTick[rTick$uptime >= from & rTick$uptime <= to, ]
    c0 <- cgc[cgc$uptime >= from & cgc$uptime <= to, ]
    kids <- pc[pc$comm != "R" & pc$uptime >= from & pc$uptime <= to, ]
    others <- cg[cg$unit != ctl & !grepl("[.]slice$", cg$unit) & cg$uptime >= from & cg$uptime <= to, ]
    otherTick <- if (nrow(others)) stats::aggregate(cg_anon ~ uptime, data = others, FUN = sum) else NULL
    data.frame(window = w, from = round(from, 2), to = round(to, 2), seconds = round(to - from, 2), ticks = nrow(r),
        r_anon_entry = at(rTick, from, rTick$rss_anon_kb) * 1024 / GiB,
        r_anon_peak = if (nrow(r)) max(r$rss_anon_kb) * 1024 / GiB else NA,
        r_anon_exit = at(rTick, to, rTick$rss_anon_kb) * 1024 / GiB,
        cg_anon_peak = if (nrow(c0)) max(c0$cg_anon) / GiB else NA,
        children = if (nrow(kids)) paste(unique(kids$comm), collapse = ",") else "",
        children_anon_peak = if (nrow(kids)) max(stats::aggregate(rss_anon_kb ~ uptime, data = kids, FUN = sum)$rss_anon_kb) * 1024 / GiB else 0,
        other_units_anon_peak = if (!is.null(otherTick)) max(otherTick$cg_anon) / GiB else 0, stringsAsFactors = FALSE)
})
win <- do.call(rbind, rows)
num <- vapply(win, is.numeric, logical(1)); win[num] <- lapply(win[num], function(x) round(x, 3))
utils::write.table(win, paste0(prefix, "-windows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

kids <- pc[pc$comm != "R", ]
if (nrow(kids)) {
    ch <- do.call(rbind, lapply(split(kids, kids$pid), function(x) data.frame(pid = x$pid[[1]], comm = x$comm[[1]],
        first_uptime = min(x$uptime), last_uptime = max(x$uptime), max_rss_anon_kb = max(x$rss_anon_kb), rows = nrow(x))))
    utils::write.table(ch, paste0(prefix, "-children.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
options(width = 200)
print(check, row.names = FALSE)
print(win, row.names = FALSE)
