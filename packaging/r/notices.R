# KIT-134 P4 (revision 2.0, sections 2.0.2 and 2.0.6): the third-party notices, the list for the
# sources archive, and the version stamp. Build time only; run with the bundle's own R.
#
# Usage: Rscript --vanilla notices.R <bundle directory> <checkout> <bundle-sources.tsv> <output directory> <electron version>
#
# Writes into <output directory>:
#   THIRD-PARTY-NOTICES.txt  every component shipped (R, every R package, the bundled shared
#                            libraries on Linux, Python and its packages, pandoc, Electron): version,
#                            licence and where its source is; and where the release's sources archive is (U6)
#   sources-list.tsv         the same components, one row each, for the sources archive
#   VERSION.json             release, lockfile sha256, commit, package commit, snapshot date, build date, runner image
#
# KIT-176: MSR_PACKAGE_COMMIT (40 hexadecimal characters) and MSR_PACKAGE_ON_RELEASE_LINE (true or false), set by the
# workflow's package-commit step (packaging/package-commit.sh), are required; nothing falls back to GITHUB_SHA.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 5L) {
    stop("usage: Rscript --vanilla notices.R <bundle> <checkout> <bundle-sources.tsv> <output directory> <electron version>")
}
bundle <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
checkout <- normalizePath(arguments[[2L]], winslash = "/", mustWork = TRUE)
sources <- utils::read.delim(arguments[[3L]], stringsAsFactors = FALSE, colClasses = "character")
out_dir <- arguments[[4L]]
electron <- arguments[[5L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
is_windows <- .Platform$OS.type == "windows"
PPM <- "https://packagemanager.posit.co"
BIOC_VERSION <- "3.23"

lockfile_path <- file.path(checkout, "packaging", "multischolar.lock.json")
lockfile <- jsonlite::fromJSON(lockfile_path, simplifyVector = FALSE)
snapshot <- lockfile$Snapshot$Date
python <- lockfile$Python
description <- read.dcf(file.path(checkout, "DESCRIPTION"))
release <- unname(description[1L, "Version"])
commit <- Sys.getenv("GITHUB_SHA", "unknown")
package_commit <- Sys.getenv("MSR_PACKAGE_COMMIT")
package_on_release_line <- Sys.getenv("MSR_PACKAGE_ON_RELEASE_LINE")
if (!grepl("^[0-9a-f]{40}$", package_commit)) {
    stop(sprintf("MSR_PACKAGE_COMMIT must be the 40-character package commit (packaging/package-commit.sh), not '%s'", package_commit))
}
if (!package_on_release_line %in% c("true", "false")) {
    stop(sprintf("MSR_PACKAGE_ON_RELEASE_LINE must be true or false, not '%s'", package_on_release_line))
}
repository <- sprintf("https://github.com/%s", Sys.getenv("GITHUB_REPOSITORY", "APAF-bioinformatics/MultiScholaR"))

rows <- list()
component <- function(kind, name, version, licence, source) {
    rows[[length(rows) + 1L]] <<- data.frame(kind = kind, name = name, version = version,
        licence = gsub("[\r\n]+", " ", licence), source = source, stringsAsFactors = FALSE)
}

r_version <- as.character(getRversion())
component("R", "R", r_version, "GPL-2 | GPL-3",
    sprintf("https://cran.r-project.org/src/base/R-%s/R-%s.tar.gz", substr(r_version, 1L, 1L), r_version))

installed <- utils::installed.packages(lib.loc = .Library, fields = c("License", "RemoteSha", "RemoteUsername", "RemoteRepo"))
# Not packages: R for Windows' library/translations (part of R) has a DESCRIPTION but no namespace.
installed <- installed[file.exists(file.path(installed[, "LibPath"], rownames(installed), "NAMESPACE")), , drop = FALSE]
for (package in sort(rownames(installed))) {
    version <- unname(installed[package, "Version"])
    licence <- unname(installed[package, "License"])
    from <- sources$origin[match(package, sources$package)]
    source <- if (installed[package, "Priority"] %in% "base") {
        sprintf("part of R %s (above)", r_version)
    } else if (is.na(from)) {
        sprintf("installed with R %s (recommended package); %s/cran/%s/src/contrib/%s_%s.tar.gz", r_version, PPM,
            snapshot, package, version)
    } else if (startsWith(from, "GitHub ")) {
        sprintf("https://github.com/%s/%s/archive/%s.tar.gz", installed[package, "RemoteUsername"],
            installed[package, "RemoteRepo"], installed[package, "RemoteSha"])
    } else if (startsWith(from, "checkout ")) {
        sprintf("%s/archive/%s.tar.gz", repository, commit)
    } else if (identical(from, "CRAN")) {
        sprintf("%s/cran/%s/src/contrib/%s_%s.tar.gz", PPM, snapshot, package, version)
    } else {
        area <- if (identical(package, "GO.db")) "data/annotation" else "bioc"
        sprintf("%s/bioconductor/%s/packages/%s/%s/src/contrib/%s_%s.tar.gz", PPM, snapshot, BIOC_VERSION, area,
            package, version)
    }
    component(if (package == "MultiScholaR") "MultiScholaR" else "R package", package, version, licence, source)
}

syslib <- file.path(bundle, "syslib", "PACKAGES.tsv")
if (file.exists(syslib)) {
    libraries <- utils::read.delim(syslib, stringsAsFactors = FALSE, colClasses = "character")
    for (i in seq_len(nrow(libraries))) {
        source_package <- if (nzchar(libraries$source_package[[i]])) libraries$source_package[[i]] else libraries$package[[i]]
        component("shared library (Ubuntu 22.04)", libraries$file[[i]],
            sprintf("%s %s", libraries$package[[i]], libraries$version[[i]]),
            sprintf("see /usr/share/doc/%s/copyright in Ubuntu 22.04", libraries$package[[i]]),
            sprintf("https://launchpad.net/ubuntu/+source/%s/%s", source_package, libraries$version[[i]]))
    }
}

if (!is.null(python$CPython)) {
    component("Python", "CPython", python$CPython, "PSF-2.0",
        sprintf("https://www.python.org/ftp/python/%s/Python-%s.tgz", python$CPython, python$CPython))
    component("Python", "python-build-standalone", python$PythonBuildStandalone,
        "BSD-3-Clause (build scripts); the libraries it links are listed in the release's licence files",
        sprintf("https://github.com/astral-sh/python-build-standalone/releases/tag/%s", python$PythonBuildStandalone))
}
site <- if (is_windows) file.path(bundle, "python", "Lib", "site-packages") else
    Sys.glob(file.path(bundle, "python", "lib", "python3*", "site-packages"))[1L]
for (metadata in Sys.glob(file.path(site, "*.dist-info", "METADATA"))) {
    lines <- readLines(metadata, warn = FALSE, encoding = "UTF-8")
    lines <- lines[seq_len(match("", lines, nomatch = length(lines) + 1L) - 1L)] # the header only
    field <- function(name) {
        value <- sub(sprintf("^%s: *", name), "", grep(sprintf("^%s:", name), lines, value = TRUE))
        if (length(value)) value[[1L]] else ""
    }
    classifiers <- sub("^Classifier: License :: ", "", grep("^Classifier: License :: ", lines, value = TRUE))
    licence <- field("License-Expression")
    if (!nzchar(licence)) licence <- field("License")
    if (!nzchar(licence) || nchar(licence) > 120L) licence <- paste(c(classifiers, if (!length(classifiers)) licence), collapse = "; ")
    name <- field("Name")
    component("Python package", name, field("Version"), licence,
        sprintf("https://pypi.org/project/%s/%s/#files", name, field("Version")))
}

pandoc <- Sys.getenv("PANDOC_VERSION", "")
component("pandoc", "pandoc", pandoc, "GPL-2.0-or-later",
    sprintf("https://github.com/jgm/pandoc/archive/refs/tags/%s.tar.gz", pandoc))
component("Electron", "Electron", electron,
    "MIT (LICENSE.electron.txt); Chromium and its components: LICENSES.chromium.html (both beside the program)",
    sprintf("https://github.com/electron/electron/archive/refs/tags/v%s.tar.gz", electron))

notices <- do.call(rbind, rows)
utils::write.table(notices, file.path(out_dir, "sources-list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# U6 (the user, 2026-11-17): no written offer and no named contact; each release publishes one sources
# archive beside its installers.
text <- c(
    sprintf("MultiScholaR %s: third-party notices", release),
    "",
    sprintf("MultiScholaR is licensed under %s (the LICENSE file of %s).", description[1L, "License"], repository),
    "This distribution also contains the components listed below, each under its own licence. The",
    "components are distributed together; none of them is modified, except that R's front-end script",
    if (is_windows) "is changed at build time so that R finds its own folder." else c(
        "is changed at build time so that R finds its own folder, and R's etc/ldpaths so that R uses only",
        "its own libraries."),
    "",
    # Decision M25 (the user, 2026-11-18): the Linux host requirements, as install.sh checks them.
    if (!is_windows) c(
        "Requirements (Linux): the GNU C library 2.35 or later (Ubuntu 22.04, Debian 12, Fedora 36 and",
        "later) and the program /usr/bin/which (Fedora, RHEL: dnf install which; Arch Linux: pacman -S",
        "which; Debian, Ubuntu: apt install debianutils; openSUSE: zypper install which).",
        ""),
    "Source code. Each component's source, at the version listed, is at the address beside it. The",
    "complete corresponding source of every component in this distribution is also published as one",
    sprintf("archive beside the installers of this release, at %s/releases. Questions:", repository),
    sprintf("%s/issues.", repository),
    "",
    if (!identical(package_commit, commit)) sprintf(paste("Analysis code from commit %s; the source address below names the",
        "commit that was built."), package_commit),
    sprintf("Built from commit %s with packages from the Posit Package Manager snapshot of %s.", commit, snapshot),
    "",
    sprintf("%-28s %-22s %-45s %s", "Component", "Version", "Licence", "Source"),
    sprintf("%-28s %-22s %-45s %s", notices$name, notices$version, notices$licence, notices$source)
)
writeLines(text, file.path(out_dir, "THIRD-PARTY-NOTICES.txt"), useBytes = TRUE)

stamp <- list(
    release = release,
    lockfile_sha256 = digest::digest(file = lockfile_path, algo = "sha256"),
    commit = commit,
    package_commit = package_commit,
    package_on_release_line = identical(package_on_release_line, "true"),
    snapshot = snapshot,
    build_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    runner_image = trimws(paste(Sys.getenv("ImageOS"), Sys.getenv("ImageVersion"))),
    platform = if (is_windows) "windows-x86_64" else "linux-x86_64",
    r = r_version,
    electron = electron,
    internal_packages = paste(vapply(c("fastduckdb", "fastlimpadpc"), function(package) if (package %in% rownames(installed))
        sprintf("%s %s", package, installed[package, "Version"]) else paste(package, "not included"), character(1)), collapse = "; ")
)
writeLines(jsonlite::toJSON(stamp, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "VERSION.json"))
cat(sprintf("notices: %d components; VERSION.json: release %s, commit %s, package commit %s (on the release line: %s), snapshot %s\n",
    nrow(notices), release, commit, package_commit, package_on_release_line, snapshot))
