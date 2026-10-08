# CORE-013 premise probe (2026-10-08, in-process; not a counted run, no performance claim).
# Times one systemctl property read (two per controller poll tick) and one durable journal write with sync.
# Usage: R_LIBS=<installed library> Rscript overhead-probe.R <directory on the project filesystem>
d <- normalizePath(commandArgs(TRUE)[1]); ns <- asNamespace("MultiScholaR")
props <- get("workflowServiceProperties", ns); sync <- get("syncWorkflowExportPaths", ns)
unit <- "multischolarwfx0123456789abcdef-w0123456789ab.service"; invisible(props(unit))
t1 <- system.time(for (i in 1:100) props(unit))[["elapsed"]] / 100
p <- file.path(d, "j.json")
t2 <- system.time(for (i in 1:30) { writeLines(as.character(i), p); sync(c(p, d), d) })[["elapsed"]] / 30
unlink(p)
cat(sprintf("RESULT systemctl_show_ms=%.1f per_tick_ms=%.1f journal_write_sync_ms=%.1f thirty_writes_s=%.2f\n",
    1000 * t1, 2000 * t1, 1000 * t2, 30 * t2))
