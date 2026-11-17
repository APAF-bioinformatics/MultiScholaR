# KIT-134 P2: the bundle's R and library, checked from the folder it was moved to, with the
# bundle's own R started the way the shell starts it (the Linux launcher or the Windows
# environment). Build-time check, not package code.
#
# Usage: Rscript --vanilla check-library.R <bundle directory> <online|offline>
#
# Checks: R 4.6.1 with its home and every library path inside the bundle; the network as stated;
# every package in the bundle's library loads (loadNamespace, in child R processes of 10 packages
# each so the DLL limit is never reached) and every namespace loaded on the way comes from the
# bundle; MultiScholaR attaches and has run_app; every version MultiScholaR's DESCRIPTION declares is
# met by the bundled package (U4). Prints "KIT-134 CHECK PASSED: <name>" lines.

arguments <- commandArgs(trailingOnly = TRUE)
is_windows <- .Platform$OS.type == "windows"
canonical <- function(paths) {
    paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
    if (is_windows) tolower(paths) else paths
}
insideBundle <- function(paths, bundle) startsWith(canonical(paths), paste0(canonical(bundle), "/"))
Sys.setenv(RGL_USE_NULL = "TRUE") # rgl: no display

if (length(arguments) >= 2L && identical(arguments[[2L]], "--child")) {
    # A child: load the named packages, then report where every loaded namespace came from.
    bundle <- arguments[[1L]]
    for (package in arguments[-(1:2)]) {
        result <- tryCatch({
            loadNamespace(package)
            sprintf("LOADED %s %s", package, as.character(utils::packageVersion(package)))
        }, error = function(e) sprintf("FAILED %s: %s", package, gsub("\n", " ", conditionMessage(e))))
        cat(result, "\n", sep = "")
    }
    for (namespace in setdiff(loadedNamespaces(), "base")) {
        path <- getNamespaceInfo(namespace, "path")
        if (!insideBundle(path, bundle)) cat(sprintf("OUTSIDE %s %s\n", namespace, path))
    }
    quit(save = "no", status = 0L)
}

if (length(arguments) != 2L || !arguments[[2L]] %in% c("online", "offline")) {
    stop("usage: Rscript --vanilla check-library.R <bundle directory> <online|offline>")
}
bundle <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
failures <- character()
check <- function(name, passed, detail = "") {
    cat(sprintf("%s %s%s\n", if (passed) "ok  " else "FAIL", name, if (nzchar(detail)) paste0(": ", detail) else ""))
    if (!passed) failures <<- c(failures, name)
}

check("R is 4.6.1", identical(as.character(getRversion()), "4.6.1"), R.version.string)
check("R home is inside the bundle", insideBundle(R.home(), bundle), R.home())
check("every library path is inside the bundle", all(insideBundle(.libPaths(), bundle)),
    paste(.libPaths(), collapse = "; "))
reached <- tryCatch({
    connection <- socketConnection("cran.r-project.org", 80L, open = "r+", blocking = TRUE, timeout = 10)
    close(connection)
    TRUE
}, error = function(e) FALSE, warning = function(w) FALSE)
check(sprintf("the network is %s", arguments[[2L]]), reached == identical(arguments[[2L]], "online"),
    sprintf("a connection to cran.r-project.org %s", if (reached) "opened" else "failed"))

packages <- utils::installed.packages(lib.loc = .libPaths())
# R for Windows keeps its message translations in library/translations, a folder with a DESCRIPTION
# but no namespace: not a package (KIT-134 run 43942816664).
packages <- packages[!packages[, "Priority"] %in% "base" &
    file.exists(file.path(packages[, "LibPath"], rownames(packages), "NAMESPACE")), , drop = FALSE]
packages <- sort(rownames(packages))
rscript <- file.path(R.home("bin"), if (is_windows) "Rscript.exe" else "Rscript")
self <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), mustWork = TRUE)
quote <- function(x) shQuote(x, type = if (is_windows) "cmd" else "sh")
started <- Sys.time()
output <- character()
for (chunk in split(packages, ceiling(seq_along(packages) / 10))) {
    output <- c(output, suppressWarnings(system2(rscript, c("--vanilla", quote(self), quote(bundle), "--child", chunk),
        stdout = TRUE, stderr = FALSE)))
}
loaded <- sub("^LOADED ([^ ]+) .*$", "\\1", grep("^LOADED ", output, value = TRUE))
failed <- grep("^FAILED ", output, value = TRUE)
outside <- unique(grep("^OUTSIDE ", output, value = TRUE))
cat(failed, outside, sep = "\n")
check(sprintf("every bundled package loads (%d packages)", length(packages)),
    length(failed) == 0L && setequal(loaded, packages),
    sprintf("%d loaded, %d failed, %d not reported, %.0f s", length(loaded), length(failed),
        length(setdiff(packages, c(loaded, sub("^FAILED ([^:]+):.*$", "\\1", failed)))),
        as.numeric(difftime(Sys.time(), started, units = "secs"))))
check("every loaded namespace comes from the bundle", length(outside) == 0L, sprintf("%d outside", length(outside)))

attach_output <- suppressWarnings(system2(rscript, c("--vanilla", "-e",
    quote("suppressPackageStartupMessages(library(MultiScholaR)); cat('MultiScholaR', as.character(packageVersion('MultiScholaR')), is.function(run_app), find.package('MultiScholaR'), '\\n')")),
    stdout = TRUE, stderr = TRUE))
attach_line <- grep("^MultiScholaR ", attach_output, value = TRUE)
check("MultiScholaR attaches from the bundle and has run_app",
    length(attach_line) == 1L && grepl(" TRUE ", attach_line) &&
        insideBundle(trimws(sub("^.* TRUE ", "", attach_line)), bundle),
    paste(utils::tail(attach_output, 3L), collapse = " | "))

# U4: every version DESCRIPTION declares (Depends, Imports, LinkingTo, Suggests) is met by the bundle.
# A declared package that is not bundled (never-installed Suggests; the internal packages in a build
# without the token) is listed, not failed: the version table accounts for it.
fields <- utils::packageDescription("MultiScholaR", lib.loc = .libPaths())[c("Depends", "Imports", "LinkingTo", "Suggests")]
entries <- trimws(unlist(strsplit(gsub("\n", " ", unlist(fields[!vapply(fields, is.null, logical(1))])), ",")))
parts <- regmatches(entries, regexec("^([A-Za-z0-9.]+)\\s*\\(\\s*(>=|>|==|<=|<)\\s*([^) ]+)\\s*\\)$", entries))
parts <- do.call(rbind, parts[lengths(parts) == 4L])
minimums <- data.frame(package = parts[, 2L], operator = parts[, 3L], declared = parts[, 4L], stringsAsFactors = FALSE)
minimums$bundled <- vapply(minimums$package, function(package) {
    if (identical(package, "R")) return(as.character(getRversion()))
    path <- find.package(package, lib.loc = .libPaths(), quiet = TRUE)
    if (length(path)) as.character(utils::packageVersion(package, lib.loc = .libPaths())) else NA_character_
}, character(1))
minimums$met <- vapply(seq_len(nrow(minimums)), function(i) {
    if (is.na(minimums$bundled[[i]])) return(NA)
    do.call(minimums$operator[[i]], list(package_version(minimums$bundled[[i]]), package_version(minimums$declared[[i]])))
}, logical(1))
cat(sprintf("MINIMUM %s %s %s: bundled %s, %s\n", minimums$package, minimums$operator, minimums$declared,
    ifelse(is.na(minimums$bundled), "none", minimums$bundled),
    ifelse(is.na(minimums$met), "not bundled", ifelse(minimums$met, "met", "NOT MET"))), sep = "")
check(sprintf("every bundled package meets DESCRIPTION's declared version (%d declared)", nrow(minimums)),
    !any(minimums$met %in% FALSE),
    sprintf("%d met, %d not met (%s), %d not bundled (%s)", sum(minimums$met %in% TRUE), sum(minimums$met %in% FALSE),
        paste(minimums$package[minimums$met %in% FALSE], collapse = ", "), sum(is.na(minimums$met)),
        paste(minimums$package[is.na(minimums$met)], collapse = ", ")))

if (length(failures)) {
    cat(sprintf("KIT-134 CHECK FAILED: %s\n", paste(failures, collapse = "; ")))
    quit(save = "no", status = 1L)
}
cat("KIT-134 CHECK PASSED: the bundle's R and library from the moved folder\n")
