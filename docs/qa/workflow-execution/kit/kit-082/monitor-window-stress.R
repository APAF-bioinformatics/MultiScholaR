# KIT-082 T01: the monitor's worker read against real transient units starting and stopping, the two windows: read before the
# unit exists, and read as it exits. The read of d52303f4 (a second look at the file after a failed read) must abort here and
# the one-read version must not. Times and counts only. Usage: R_LIBS=<library> Rscript .../monitor-window-stress.R [units]
suppressMessages(library(MultiScholaR)); current <- get("workflowWorkerMemory", asNamespace("MultiScholaR"))
counter <- get("workflowMemoryCounter", asNamespace("MultiScholaR"))
second_look <- function(path) tryCatch(counter(path), workflow_execution_unsupported = function(error) { if (file.exists(path)) stop(error); NULL })
check_then_read <- function(path) if (file.exists(path)) counter(path)
manager <- system2("systemctl", c("--user", "show", "-p", "ControlGroup", "--value"), stdout = TRUE)
units <- { a <- commandArgs(TRUE); if (length(a)) as.integer(a[[1L]]) else 40L }
stress <- function(label, read) {
    aborted <- numbers <- nones <- 0L
    for (index in seq_len(units)) {
        unit <- sprintf("multischolar-kit082-stress-%d-%s-%d", Sys.getpid(), label, index)
        path <- file.path("/sys/fs/cgroup", sub("^/", "", manager), "app.slice", paste0(unit, ".service"), "memory.current")
        system2("systemd-run", c("--user", "--quiet", "--collect", "--no-block", paste0("--unit=", unit), "sleep", "0.15"), stdout = FALSE, stderr = FALSE)
        until <- Sys.time() + 0.6
        while (Sys.time() < until) {
            value <- tryCatch(suppressWarnings(read(path)), error = function(error) "aborted")
            if (identical(value, "aborted")) aborted <- aborted + 1L else if (is.null(value)) nones <- nones + 1L else numbers <- numbers + 1L
        }
    }
    cat(sprintf("%-16s units %d reads with a number %d, none %d, ABORTS %d\n", label, units, numbers, nones, aborted))
}
stress("second-look", second_look)
stress("check-then-read", check_then_read)
stress("one-read", current)
