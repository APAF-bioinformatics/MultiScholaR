# KIT-130 rev 1.1, the dependency closure: resolve MultiScholaR's whole DESCRIPTION (Depends,
# Imports, LinkingTo, Suggests; without the two internal Suggests the workflow cannot read) into
# a pak lockfile and list every dependency that has no binary for this R on this platform.
# Resolve only: nothing is installed.
#
# Usage (from the package root): Rscript dependency-closure.R <output prefix>
# Writes <prefix>.lock, <prefix>-all.csv and <prefix>-no-binary.csv; needs pak and jsonlite.

out_prefix <- commandArgs(trailingOnly = TRUE)[[1L]]
refs <- c("deps::.", "fastduckdb=?ignore", "fastlimpadpc=?ignore")
lockfile <- paste0(out_prefix, ".lock")
summary_file <- Sys.getenv("GITHUB_STEP_SUMMARY", unset = NA)
platform_label <- sprintf("%s, %s", R.version.string, R.version[["platform"]])
# setup-r leaves CRAN unset on Windows; binaries for Windows are CRAN's own.
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
}
cat(sprintf("repositories: %s\n", paste(names(getOption("repos")), getOption("repos"), sep = "=", collapse = ", ")))

appendSummary <- function(lines) {
    cat(lines, sep = "\n")
    if (!is.na(summary_file)) cat(lines, file = summary_file, sep = "\n", append = TRUE)
}

tryCatch(
    pak::lockfile_create(refs, lockfile = lockfile, dependencies = TRUE),
    error = \(e) {
        appendSummary(c(
            sprintf("### Dependency closure (%s): resolution failed", platform_label),
            "```",
            conditionMessage(e),
            "```"
        ))
        stop(e)
    }
)

lock <- jsonlite::fromJSON(lockfile, simplifyVector = FALSE)
packages <- Filter(\(p) !identical(p[["type"]], "deps"), lock[["packages"]])
cat("fields of a lockfile entry:", paste(names(packages[[1L]]), collapse = ", "), "\n")
field <- \(name) {
    vapply(packages, \(p) {
        value <- unlist(p[[name]])
        if (length(value) == 0L) NA_character_ else as.character(value[[1L]])
    }, character(1L))
}
closure <- data.frame(
    package = field("package"),
    version = field("version"),
    source = field("type"),
    platform = field("platform"),
    binary = field("binary"),
    needs_compilation = field("needscompilation"),
    ref = field("ref")
)
closure <- closure[order(tolower(closure$package)), ]
# No binary: pak marks the entry as not binary, or resolves it to the source platform.
no_binary <- closure[!(closure$binary %in% "TRUE") | closure$platform %in% "source", ]

utils::write.csv(closure, paste0(out_prefix, "-all.csv"), row.names = FALSE)
utils::write.csv(no_binary, paste0(out_prefix, "-no-binary.csv"), row.names = FALSE)

appendSummary(c(
    sprintf("### Dependency closure (%s)", platform_label),
    "",
    sprintf(
        "%d dependencies resolved; %d have a binary; %d have none, of which %d need compilation.",
        nrow(closure), nrow(closure) - nrow(no_binary), nrow(no_binary),
        sum(!(no_binary$needs_compilation %in% "FALSE"))
    ),
    "",
    "| Package | Version | Source | Needs compilation | Reference |",
    "|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s | %s |", no_binary$package, no_binary$version, no_binary$source,
        no_binary$needs_compilation, no_binary$ref)
))
