# run-tests: a directory checked before it is uploaded (the patterns of .github/k124-sim/check-output.R).
#
# Usage: Rscript --vanilla scan.R <dir> <strict|check|final>
#   strict  every pattern (test logs, failures, script outputs): compiler or build output, source locations and
#           names of the internal packages, anything shaped like a token; a file that matches is removed (withheld)
#           and named, with the pattern's name, in <dir>/withheld.txt; never its text.
#   check   MultiScholaR's own R CMD check log: its compiler lines are allowed; the internal packages' names and
#           sources, Rust build output and tokens are not (a match is withheld as above).
#   final   the assembled artifact: the strict patterns minus the compiler ones (the check log is in it); any match
#           fails (exit status 1) and nothing is uploaded.
# Compiled objects and archives (.so .o .a .dll .rlib .tar.gz .tgz .zip) are always withheld (final: fail).

PATTERNS <- c(
    rust_location = "\\.rs:[0-9]", rust_tree = "src/rust|\\.cargo|Cargo\\.(toml|lock)", rust_build = "rustc|cargo |panicked at",
    rust_messages = "Compiling [a-z_]+ v[0-9]|error\\[E[0-9]+\\]|warning: unused",
    c_compiler = "(^|[ /])(gcc|g\\+\\+|clang|cc1|ld)( |:)|-fPIC|-std=gnu", make = "make(\\[[0-9]+\\])?: \\*\\*\\*",
    c_location = "\\.(c|cc|cpp|h|hpp|o):[0-9]+", internal_repository = "fast_limpa_dpc|fast_duckdb|internal-build-logs",
    token = "gh[pousr]_[A-Za-z0-9]{20,}|github_pat_|Bearer "
)
COMPILER <- c("c_compiler", "make", "c_location")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !arguments[[2L]] %in% c("strict", "check", "final")) {
    stop("usage: Rscript --vanilla scan.R <dir> <strict|check|final>")
}
dir <- normalizePath(arguments[[1L]])
mode <- arguments[[2L]]
patterns <- if (mode == "strict") PATTERNS else PATTERNS[setdiff(names(PATTERNS), COMPILER)]

values <- function(path) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        object <- tryCatch(readRDS(path), error = function(e) NULL)
        return(unlist(rapply(list(object), function(x) as.character(x), classes = c("character", "factor"), how = "list"),
            use.names = FALSE))
    }
    suppressWarnings(readLines(path, warn = FALSE))
}

files <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
files <- files[basename(files) != "withheld.txt"]
withheld <- character()
for (path in files) {
    relative <- substring(path, nchar(dir) + 2L)
    hit <- if (grepl("\\.(so|o|a|dll|rlib|tar\\.gz|tgz|zip)$", path, ignore.case = TRUE)) "compiled_or_archive" else {
        text <- values(path)
        names(patterns)[vapply(patterns, function(p) any(grepl(p, text, useBytes = TRUE)), logical(1L))]
    }
    if (length(hit)) {
        withheld <- c(withheld, sprintf("%s\t%s", relative, paste(hit, collapse = ",")))
        if (mode != "final") unlink(path)
    }
}
cat(sprintf("%d files scanned (%s)\n", length(files), mode))
if (length(withheld)) {
    if (mode == "final") {
        cat(sprintf("::error::%s matches a forbidden pattern\n", withheld), sep = "")
        quit(save = "no", status = 1L)
    }
    cat(sprintf("withheld (not uploaded): %s\n", withheld), sep = "")
    cat(withheld, file = file.path(dir, "withheld.txt"), sep = "\n", append = TRUE)
}
