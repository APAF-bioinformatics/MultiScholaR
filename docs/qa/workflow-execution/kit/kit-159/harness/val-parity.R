# KIT-159 section 7.9: the parity step inside a bundle job of distributable.yml (stage parity).
# Run with the installed bundle's Rscript, after the test-tool library step, before the uninstall:
#
#   Rscript --vanilla <harness root>/docs/qa/workflow-execution/kit/kit-159/harness/val-parity.R \
#       <request.json> <installed application folder> <test-tool library> <harness root> <out dir>
#
# The harness comes from the request's harness_commit (git archive, no token). Per route of the
# request: (a) the driver in process on the declared generated input; (b) the installed application
# started as the shell starts it (entry.R, launcher probe `process`, test mode in that process only);
# (c) shinytest2::AppDriver performs the same actions and records what the browser sent; (d) the
# ledgers compared (parity-<route>.tsv); (e) the bundle driver's DA table against L1's committed
# digest and fingerprint (parity-da.tsv). valCheckEvidence runs here; check-evidence.sh after.
# Environment: MSR_PKG_PATHS (KIT-176's PKG_PATHS, read from packaging/package-commit.sh).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("usage: val-parity.R <request> <app> <test tools> <harness root> <out>",
                             call. = FALSE)
request <- args[[1L]]
app <- normalizePath(args[[2L]])
tools_lib <- args[[3L]]
harness_root <- args[[4L]]
out_dir <- args[[5L]]
harness <- file.path(harness_root, "docs/qa/workflow-execution/kit/kit-159/harness")
started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

val_lib <- file.path(tempdir(), "val-lib")
dir.create(val_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(val_lib, tools_lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
status <- system2(file.path(R.home("bin"), "R"),
                  c("CMD", "INSTALL", paste0("--library=", val_lib), file.path(harness, "msrval")),
                  stdout = FALSE, stderr = FALSE)
if (!identical(as.integer(status), 0L)) stop("msrval did not install", call. = FALSE)

req <- msrval::valParityRequest(request)
# r1.12 (7.9 item 1): the request's copied L1 digests and fingerprints equal L1's committed file.
msrval::valCheckParityL1(req, repo = ".", at = req$evidence_commit %||% "HEAD")
version_file <- c(file.path(app, "VERSION.json"), file.path(app, "resources", "VERSION.json"))
version_file <- version_file[file.exists(version_file)]
if (length(version_file) == 0L) {
    version_file <- list.files(app, "^VERSION\\.json$", recursive = TRUE, full.names = TRUE)
}
version <- jsonlite::read_json(version_file[[1L]], simplifyVector = FALSE)
pkg_paths <- strsplit(trimws(Sys.getenv("MSR_PKG_PATHS")), "[[:space:]]+")[[1L]]
if (length(pkg_paths) == 0L || !nzchar(pkg_paths[[1L]])) stop("MSR_PKG_PATHS is empty", call. = FALSE)
bundle_check <- msrval::valCheckBundleCommit(version, req$release_commit, pkg_paths, repo = ".")
msrval::valStatus(phase = "parity", status = paste("bundle commit:", bundle_check))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (x in req$inputs) {
    if (!identical(digest::digest(file = x$path, algo = "sha256"), x$sha256)) {
        msrval:::valAbort("parity input {x$path} differs from its declared SHA-256",
                          "val_input_hash_mismatch")
    }
}
# The committed exclusion list of client-state inputs (on the run branch, outside the harness).
excl_path <- req$exclusions %||% "docs/qa/workflow-execution/kit/kit-159/parity/exclusions.tsv"
exclusions <- if (file.exists(excl_path)) {
    utils::read.delim(excl_path, colClasses = "character")
} else {
    data.frame(input_id = character(), class = character())
}
helpers <- msrval::valGuiHelpers(".")
entry <- file.path(app, c("resources", "app/resources"), "entry.R")
entry <- entry[file.exists(entry)][[1L]]
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
da_rows <- list()
for (route in unlist(req$routes)) {
    spec <- req$route_inputs[[route]]
    files <- lapply(spec$files, function(p) normalizePath(p))
    design <- utils::read.delim(spec$design, stringsAsFactors = FALSE, check.names = FALSE)
    script <- msrval::valActionScript(sub("_two_species$", "", route),
                                      assays = unlist(spec$assays), import_tool = spec$import_tool,
                                      import_values = spec$import_values %||% list(),
                                      contrasts = lapply(unlist(spec$contrasts), function(x) {
                                          strsplit(x, "_vs_")[[1L]]
                                      }))
    work <- file.path(tempdir(), "parity", route)
    drv <- msrval::valRunWorkflow(script, files, design, helpers = helpers, seed = 20261107L,
                                  uniprot_cache = spec$uniprot_cache,
                                  tables_dir = file.path(work, "tables"),
                                  work_dir = file.path(work, "driver"), pkg_root = ".")
    if (!identical(drv$exit, "completed")) stop(sprintf("driver on %s: %s", route, drv$exit))
    port <- 18000L + match(route, unlist(req$routes))
    run_folder <- file.path(work, "app-run")
    dir.create(run_folder, recursive = TRUE, showWarnings = FALSE)
    log <- file.path(work, "app.log")
    proc <- processx::process$new(
        rscript, c("--vanilla", entry, as.character(port), run_folder, as.character(Sys.getpid())),
        env = c("current", MSR_LAUNCHER_PROBE = '{"launcher":"process","reason":"hosted runner"}',
                MULTISCHOLAR_TEST_MODE = "true"),
        stdout = log, stderr = log)
    url <- sprintf("http://127.0.0.1:%d", port)
    for (i in seq_len(120L)) {
        ok <- tryCatch({
            curl::curl_fetch_memory(url)$status_code == 200L
        }, error = function(e) FALSE)
        if (ok) break
        Sys.sleep(1)
    }
    launcher <- grep("^MSR-LAUNCHER ", readLines(log, warn = FALSE), value = TRUE)
    msrval::valStatus(phase = "parity", status = paste(route, "launcher line seen:",
                                                       length(launcher) > 0L))
    browser <- shinytest2::AppDriver$new(url, load_timeout = 120000)
    prefix <- paste0(script$omic, "_workflow-")
    sent <- msrval::valParityBrowser(browser, drv$actions, prefix)
    browser$stop()
    proc$kill()
    cmp <- msrval::valCompareLedgers(drv$ledger, sent, exclusions)
    utils::write.table(cmp, file.path(out_dir, sprintf("parity-%s.tsv", route)), sep = "\t",
                       quote = FALSE, row.names = FALSE)
    da <- utils::read.delim(drv$tables[["da"]], stringsAsFactors = FALSE)
    l1 <- req$l1[[route]]
    if (!is.null(l1)) {
        pd <- msrval::valCompareParityDa(da, l1)
        da_rows[[length(da_rows) + 1L]] <- data.frame(route = route, comparison = "bundle_vs_l1",
                                                      status = pd$status, digest = pd$digest,
                                                      sha256_10 = pd$sha256_10)
    }
}
utils::write.table(if (length(da_rows)) do.call(rbind, da_rows) else
                       data.frame(route = character(), comparison = character(),
                                  status = character(), digest = character(),
                                  sha256_10 = character()),
                   file.path(out_dir, "parity-da.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)
manifest <- data.frame(role = vapply(req$inputs, function(x) x$role, character(1)), kind = "data",
                       licence = "generated", redistribution = "permitted",
                       citation = "KIT-159 generated parity input", stringsAsFactors = FALSE)
manifest$sha256 <- vapply(req$inputs, function(x) x$sha256, character(1))
utils::write.table(manifest, file.path(out_dir, "attribution.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE)
prov <- msrval::valProvenance(
    ticket = "KIT-159", stage = "parity", label = basename(request), phase = "parity",
    status = "completed", started = started,
    code = list(package_commit = version$package_commit %||% "unknown",
                release_commit = req$release_commit,
                harness = list(commit = req$harness_commit, tree = req$harness_tree),
                declaration = list(), bundle_check = bundle_check),
    inputs = list(files = req$inputs)
)
msrval::valWriteProvenance(prov, file.path(out_dir, "provenance.json"))
msrval::valCheckProvenance(file.path(out_dir, "provenance.json"))
msrval::valCheckEvidence(out_dir, manifest, phase = "parity")
msrval::valWriteOutputsSha(out_dir)
