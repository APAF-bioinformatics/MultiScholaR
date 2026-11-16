# KIT-134 P2 (revision 2.0, section 2.0.3): the version table. One row per package in the lockfile
# or in the bundle: locked version and source against bundled version and source. Build time only.
#
# Usage: Rscript --vanilla version-table.R <checkout> <bundle library> <bundle-sources.tsv> <install log> <output directory>
#
# <install log> is the output of install-library.R; it tells which packages R installed as binaries
# and which it built from source. Writes version-table.tsv and version-table.md (the summary and the
# rows that differ, for the job summary). Fails when a hard-closure package of the lockfile is
# missing from the bundle, or when duckdb is not the version the code pins.
# Needs jsonlite on the library path.

PINNED <- c(duckdb = "1.5.5")
INTERNAL <- c("fastduckdb", "fastlimpadpc")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 5L) {
    stop("usage: Rscript --vanilla version-table.R <checkout> <bundle library> <bundle-sources.tsv> <install log> <output directory>")
}
checkout <- arguments[[1L]]
library_path <- arguments[[2L]]
sources <- utils::read.delim(arguments[[3L]], stringsAsFactors = FALSE, colClasses = "character")
log_lines <- readLines(arguments[[4L]], warn = FALSE)
out_dir <- arguments[[5L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lockfile <- jsonlite::fromJSON(file.path(checkout, "packaging", "multischolar.lock.json"), simplifyVector = FALSE)
excluded <- readLines(file.path(checkout, "packaging", "r", "library-exclude.txt"), warn = FALSE)
excluded <- trimws(excluded[!grepl("^\\s*(#|$)", excluded)])
`%||%` <- function(a, b) if (is.null(a)) b else a

locked <- do.call(rbind, lapply(lockfile$Packages, function(p) {
    sha <- p$RemoteSha %||% ""
    source <- switch(p$Source,
        GitHub = sprintf("GitHub %s/%s@%s", p$RemoteUsername, p$RemoteRepo, substr(sha, 1L, 8L)),
        Bioconductor = sprintf("Bioconductor %s", p$Bioconductor %||% ""),
        local = sprintf("local build; Remotes pin %s", sub("^.*@", "", p$DescriptionRemotes %||% "none")),
        p$Source)
    data.frame(package = p$Package, locked_version = p$Version %||% NA_character_, locked_source = source,
        locked_sha = sha, closure = p$Closure, locked_compiled = isTRUE(p$NeedsCompilation),
        stringsAsFactors = FALSE)
}))

installed <- utils::installed.packages(lib.loc = library_path,
    fields = c("RemoteSha", "RemoteUsername", "RemoteRepo", "NeedsCompilation"))
bundled <- data.frame(package = rownames(installed), bundled_version = unname(installed[, "Version"]),
    priority = unname(installed[, "Priority"]), bundled_sha = unname(installed[, "RemoteSha"]),
    bundled_compiled = unname(installed[, "NeedsCompilation"]) %in% "yes", stringsAsFactors = FALSE)

# Which packages R built from source ("installing *source* package"); the rest were binaries.
built <- regmatches(log_lines, regexpr("installing \\*source\\* package .[^ ]+. ", log_lines))
built <- unique(gsub("^installing \\*source\\* package .|. $", "", built))

table <- merge(locked, bundled, by = "package", all = TRUE)
from <- match(table$package, sources$package)
table$bundled_source <- ifelse(is.na(table$bundled_version), "",
    ifelse(!is.na(table$priority) & table$priority == "base", sprintf("R %s", getRversion()),
        ifelse(is.na(from), "R installer (recommended, not reinstalled)",
            ifelse(grepl("^(GitHub|checkout)", sources$origin[from]), sources$origin[from],
                sprintf("snapshot %s %s, %s", sources$snapshot[from], sources$origin[from],
                    ifelse(table$package %in% built, "built from source", "binary"))))))
table$same <- !is.na(table$locked_version) & !is.na(table$bundled_version) &
    table$locked_version == table$bundled_version &
    (table$locked_sha == "" | (!is.na(table$bundled_sha) & table$locked_sha == table$bundled_sha))
table$same <- ifelse(table$same, "yes", "no")
table$closure[is.na(table$closure)] <- "not in the lockfile"
table$compiled <- ifelse(is.na(table$bundled_version), ifelse(table$locked_compiled %in% TRUE, "yes", "no"),
    ifelse(table$bundled_compiled, "yes", "no"))

note <- character(nrow(table))
add <- function(rows, text) {
    rows <- rows %in% TRUE
    if (length(text) > 1L) text <- text[rows]
    note[rows] <<- ifelse(nzchar(note[rows]), paste0(note[rows], "; ", text), text)
}
missing <- is.na(table$bundled_version)
add(missing & table$package %in% excluded, "excluded (library-exclude.txt)")
add(missing & table$package %in% INTERNAL, "internal: needs the read token; not in this build")
add(missing & table$locked_source == "not installed", "never installed on the laptop; not added")
add(missing & !table$package %in% c(excluded, INTERNAL) & table$locked_source != "not installed",
    "left out: needed only by excluded tools")
add(table$package %in% names(PINNED), sprintf("pinned by the code to %s", PINNED[table$package]))
outside <- !missing & table$closure == "not in the lockfile"
add(outside & table$priority %in% c("base", "recommended"), "part of R, outside the closure")
add(outside & !table$priority %in% c("base", "recommended"), "not in the verified library (D7: MOFA2 and what it adds)")
add(table$package == "MOFA2", "not in the verified library (D7)")
add(table$package == "MultiScholaR", "this build's source")
table$note <- note

columns <- c("package", "locked_version", "locked_source", "bundled_version", "bundled_source", "same",
    "closure", "compiled", "note")
table <- table[order(tolower(table$package)), columns]
utils::write.table(table, file.path(out_dir, "version-table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

in_both <- !is.na(table$locked_version) & !is.na(table$bundled_version)
differs <- in_both & table$same == "no"
summary <- c(
    sprintf("rows: %d (lockfile %d, bundle %d)", nrow(table), sum(!is.na(table$locked_source)),
        sum(!is.na(table$bundled_version))),
    sprintf("in both: %d; same version: %d; different version: %d (hard closure %d, compiled %d)", sum(in_both),
        sum(in_both & table$same == "yes"), sum(differs), sum(differs & table$closure == "hard"),
        sum(differs & table$compiled == "yes")),
    sprintf("in the lockfile, not bundled: %d (%s)", sum(missing & !is.na(table$locked_version)),
        paste(table$package[missing & !is.na(table$locked_version)], collapse = ", ")),
    sprintf("bundled, not in the lockfile: %d (%s)", sum(!is.na(table$bundled_version) & table$closure == "not in the lockfile"),
        paste(table$package[!is.na(table$bundled_version) & table$closure == "not in the lockfile"], collapse = ", ")),
    sprintf("bundled packages built from source: %s", paste(sort(intersect(built, table$package)), collapse = ", "))
)
cat(summary, sep = "\n")
rows <- table[differs, , drop = FALSE]
markdown <- c("### Version table: the bundle against the laptop's lockfile", "", paste0("- ", summary), "",
    "| package | locked | bundled | closure | compiled | bundled from |", "|---|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s | %s | %s |", rows$package, rows$locked_version, rows$bundled_version,
        rows$closure, rows$compiled, rows$bundled_source))
writeLines(markdown, file.path(out_dir, "version-table.md"))

problems <- character()
hard_missing <- table$package[missing & table$closure == "hard"]
if (length(hard_missing)) problems <- c(problems, sprintf("hard-closure packages missing from the bundle: %s",
    paste(hard_missing, collapse = ", ")))
for (package in names(PINNED)) {
    version <- table$bundled_version[table$package == package]
    if (!identical(version, PINNED[[package]])) problems <- c(problems, sprintf("%s is %s, the code pins %s",
        package, if (length(version)) version else "missing", PINNED[[package]]))
}
if (length(problems)) stop(paste(problems, collapse = "; "))
cat("KIT-134 CHECK PASSED: version table complete; no hard-closure package missing; duckdb", PINNED[["duckdb"]], "\n")
