# KIT-134 (the Linux library path): the bundled R on one host, started through the bundle's launcher.
# Shows that the bundled shared libraries serve R and its child R processes without reaching the
# host programs R starts. Build-time check, not package code.
#
# Usage: <runtime>/bin/Rscript --vanilla host-check.R <runtime directory> <entry.R> <label>
#
# Checks: (1) the environment R's children inherit has no path into the bundle's syslib/;
# (2) `bash -c true` and `sh -c true` started from R (system, system2, processx) exit 0;
# (3) a child R started as the application starts its workers (R's own Rscript) loads packages
# that need bundled libraries (ICU, libxml2, libcurl); (4) the application's entry starts with the
# plain launcher, whose trial unit is a child R process, and answers HTTP 200. Printed, not
# decided: whether the path the launcher exported before this repair breaks this host's bash (the
# defect, where the host has it), which host libraries beside the C library R has mapped, and
# whose libreadline and ICU are in use.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("usage: Rscript --vanilla host-check.R <runtime directory> <entry.R> <label>")
runtime <- normalizePath(arguments[[1L]], mustWork = TRUE)
entry <- normalizePath(arguments[[2L]], mustWork = TRUE)
label <- arguments[[3L]]
syslib <- file.path(runtime, "syslib")
failures <- character()
check <- function(name, passed, detail = "") {
    cat(sprintf("%s %s%s\n", if (isTRUE(passed)) "ok  " else "FAIL", name, if (nzchar(detail)) paste0(": ", detail) else ""))
    if (!isTRUE(passed)) failures <<- c(failures, name)
}
os <- tryCatch(sub('^PRETTY_NAME="?([^"]*)"?$', "\\1", grep("^PRETTY_NAME=", readLines("/etc/os-release"), value = TRUE)[[1L]]),
    error = function(e) "unknown")
cat(sprintf("HOST %s: %s; kernel %s; %s; /bin/sh -> %s\n", label, os, Sys.info()[["release"]],
    tryCatch(system2("getconf", "GNU_LIBC_VERSION", stdout = TRUE), error = function(e) "glibc unknown"),
    Sys.readlink("/bin/sh")))

inherited <- Sys.getenv("LD_LIBRARY_PATH")
cat(sprintf("LD_LIBRARY_PATH inherited by R's children: %s\n", inherited))
check("children inherit no path into the bundle's syslib", !grepl(syslib, inherited, fixed = TRUE))

old <- suppressWarnings(system(sprintf("LD_LIBRARY_PATH=%s bash -c true 2>/dev/null", shQuote(syslib))))
cat(sprintf("for the record, bash under the path exported before the repair: %s\n",
    if (identical(old, 0L)) "runs on this host (the defect does not show here)" else sprintf("fails on this host (exit %d): the defect", old)))

check("system: bash -c true", identical(system("bash -c true"), 0L))
check("system: sh -c true", identical(system("sh -c true"), 0L))
check("system2: sh -c with a pipe of host programs", identical(system2("sh", c("-c", shQuote("ls / | wc -l > /dev/null"))), 0L))
shelled <- tryCatch(processx::run("bash", c("-c", "echo shelled"))$stdout, error = function(e) conditionMessage(e))
check("processx: bash", identical(trimws(shelled), "shelled"), if (identical(trimws(shelled), "shelled")) "" else shelled)

# The child also uses LAPACK and reports the shared objects it has mapped (for the table below).
child <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e",
    shQuote(paste0('for (p in c("stringi", "xml2", "curl", "MultiScholaR")) loadNamespace(p); invisible(solve(matrix(c(2, 1, 1, 3), 2))); ',
        'm <- unique(sub("^\\\\S+ +\\\\S+ +\\\\S+ +\\\\S+ +\\\\S+ +", "", grep(" /.*\\\\.so[.0-9]*$", readLines("/proc/self/maps"), value = TRUE))); ',
        'cat(paste0("CHILDMAP ", m), sep = "\\n"); cat("CHILD-LOADED\\n")'))),
    stdout = TRUE, stderr = TRUE))
check("child R loads stringi, xml2, curl and MultiScholaR", any(child == "CHILD-LOADED"),
    if (any(child == "CHILD-LOADED")) "" else paste(utils::tail(child, 3L), collapse = " | "))

for (package in c("stringi", "xml2", "curl")) loadNamespace(package)
invisible(solve(matrix(c(2, 1, 1, 3), 2))) # loads LAPACK
cat(sprintf("La_library(): %s; extSoftVersion BLAS: %s\n", La_library(), extSoftVersion()[["BLAS"]]))
# /proc/self/maps: address, permissions, offset, device, inode, then the path (which may hold spaces)
mapped <- sub("^\\S+ +\\S+ +\\S+ +\\S+ +\\S+ +", "", grep(" /.*\\.so[.0-9]*$", readLines("/proc/self/maps"), value = TRUE))
mapped <- unique(mapped)
inside <- startsWith(mapped, paste0(runtime, "/"))
cat(sprintf("shared objects mapped in this R: %d from the bundle, %d from the host: %s\n", sum(inside), sum(!inside),
    paste(sort(basename(mapped[!inside])), collapse = " ")))
# Printed, not decided: R's own etc/ldpaths puts /usr/local/lib and /usr/lib/x86_64-linux-gnu in front
# of the run path, so on Debian and Ubuntu a host library of the same name is taken first.
cat(sprintf("libreadline and ICU in use: %s\n", paste(grep("/(libreadline|libicuuc)", mapped, value = TRUE), collapse = " ")))
# Decision M26's measurement: which copy of readline, BLAS, LAPACK, curl and the Fortran runtime each R uses.
roles <- c(readline = "/libreadline[.]so", BLAS = "/lib(R|open|flexi)?blas[^/]*[.]so", LAPACK = "/lib(R)?lapack[^/]*[.]so",
    curl = "/libcurl[^/]*[.]so", gfortran = "/libgfortran[.]so", quadmath = "/libquadmath[.]so")
childMapped <- sub("^CHILDMAP ", "", grep("^CHILDMAP ", child, value = TRUE))
for (process in c("R", "child R")) {
    paths <- if (process == "R") mapped else childMapped
    for (role in names(roles)) {
        found <- grep(roles[[role]], paths, value = TRUE)
        cat(sprintf("LIBS\t%s\t%s\t%s\t%s\t%s\n", label, process, role,
            if (!length(found)) "none" else paste(ifelse(startsWith(found, paste0(runtime, "/")), "bundle", "HOST"), collapse = ","),
            if (!length(found)) "-" else paste(gsub(runtime, "<bundle>", found, fixed = TRUE), collapse = " ")))
    }
}

# The application's entry, as the shell starts it; this process stands for the shell (the entry's
# watchdog ends the application when it is gone).
run <- tempfile("msr-host-check-run-")
dir.create(run)
port <- httpuv::randomPort()
application <- processx::process$new(file.path(runtime, "bin", "Rscript"),
    c("--vanilla", entry, port, run, Sys.getpid()), stdout = file.path(run, "out.log"), stderr = file.path(run, "err.log"))
answers <- function() {
    connection <- url(sprintf("http://127.0.0.1:%d/", port))
    on.exit(try(close(connection), silent = TRUE), add = TRUE)
    tryCatch(length(suppressWarnings(readLines(connection, n = 1L, warn = FALSE))) > 0L, error = function(e) FALSE)
}
answered <- FALSE
deadline <- Sys.time() + 240
while (!answered && application$is_alive() && Sys.time() < deadline) {
    Sys.sleep(2)
    answered <- answers()
}
launcher <- grep("^MSR-LAUNCHER ", readLines(file.path(run, "out.log"), warn = FALSE), value = TRUE)
cat(if (length(launcher)) launcher[[1L]] else "no MSR-LAUNCHER line", "\n")
check("the entry chose the plain launcher (its trial unit, a child R, ran)", length(launcher) > 0L && grepl('"launcher":"process"', launcher[[1L]], fixed = TRUE))
check("the application answers HTTP 200", answered,
    if (answered) "" else paste(utils::tail(readLines(file.path(run, "err.log"), warn = FALSE), 5L), collapse = " | "))
try(application$kill_tree(), silent = TRUE)
try(application$kill(), silent = TRUE)

if (length(failures)) {
    cat(sprintf("KIT-134 HOST CHECK FAILED: %s: %s\n", label, paste(failures, collapse = "; ")))
    quit(save = "no", status = 1L)
}
cat(sprintf("KIT-134 HOST CHECK PASSED: %s\n", label))
