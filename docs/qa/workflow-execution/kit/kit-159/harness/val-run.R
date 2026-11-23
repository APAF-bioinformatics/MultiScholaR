# KIT-159 entry script (section 7.8). From the package root, in a `run-tests` scripts request:
#
#   Rscript docs/qa/workflow-execution/kit/kit-159/harness/val-run.R record <evidence dir> <label>
#   Rscript .../val-run.R pre <evidence dir> <label> --script <path> --pre-outputs <path>
#   Rscript .../val-run.R check <folder> --manifest <inputs.tsv>
#
# record: draft-inputs.tsv of <evidence dir> -> val-out/record/<label>/ (inputs.tsv, members-*.tsv,
# library.tsv, attribution.tsv, provenance.json, outputs.sha256); no declaration check.
# pre: the ticket's pre-declaration script on the recorded inputs (val-out/record/<label>/inputs.tsv
# of the same job) -> val-out/pre/<label>/ (the listed typed tables only).
# check: the allow-list and hashes on a committed folder (valCheckEvidence, runs.tsv corrections).
# The phases that run cells (lifecycle, counted, score, determinism, correction, amended, parity)
# are not yet wired to the driver core (msrval's valRunWorkflow and valRunCells); they refuse.

harness_dir <- "docs/qa/workflow-execution/kit/kit-159/harness"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: val-run.R <phase> <evidence dir> ...", call. = FALSE)
phase <- args[[1L]]
evidence_dir <- args[[2L]]
opts <- args[-(1:2)]
optValue <- function(name, default = NULL) {
    i <- match(name, opts)
    if (is.na(i) || i == length(opts)) default else opts[[i + 1L]]
}
label <- if (length(opts) && !startsWith(opts[[1L]], "--")) opts[[1L]] else "run"
tmp_root <- Sys.getenv("RUNNER_TEMP", tempdir())

# --- Step 1: install msrval from this checkout (7.8) ---
val_lib <- file.path(tmp_root, "val-lib")
if (!requireNamespace("msrval", lib.loc = val_lib, quietly = TRUE)) {
    dir.create(val_lib, recursive = TRUE, showWarnings = FALSE)
    status <- system2(file.path(R.home("bin"), "R"),
                      c("CMD", "INSTALL", paste0("--library=", val_lib),
                        file.path(harness_dir, "msrval")),
                      stdout = FALSE, stderr = FALSE)
    if (!identical(as.integer(status), 0L)) stop("msrval did not install", call. = FALSE)
}
.libPaths(c(val_lib, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

# --- Step 2: private logs when any data input is not permitted (5.2, A27) ---
draft_path <- file.path(evidence_dir, "draft-inputs.tsv")
draft <- if (file.exists(draft_path)) {
    utils::read.delim(draft_path, colClasses = "character", na.strings = character())
} else {
    NULL
}
private_needed <- !is.null(draft) && any(draft$kind == "data" & draft$redistribution != "permitted")
if (private_needed && !nzchar(Sys.getenv("VAL_TEMPLATE_FILE"))) {
    lines <- msrval::valRunAsChild(file.path(harness_dir, "val-run.R"), args,
                                   file.path(tmp_root, "val-private"))
    quit(status = attr(lines, "status"))
}

finish <- function(out, manifest, prov_phase, outputs = list()) {
    utils::write.table(manifest[, c("role", "kind", "licence", "redistribution", "citation")],
                       file.path(out, "attribution.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    prov <- msrval::valProvenance(
        ticket = basename(evidence_dir), stage = "record", label = label, phase = prov_phase,
        status = "completed", started = started,
        code = list(package_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                    release_commit = "none",
                    harness = list(commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                                   tree = msrval::valHarnessTree(".")),
                    declaration = list())
    )
    msrval::valWriteProvenance(prov, file.path(out, "provenance.json"))
    msrval::valCheckProvenance(file.path(out, "provenance.json"))
    msrval::valCheckEvidence(out, manifest, outputs = outputs, phase = prov_phase)
    msrval::valWriteOutputsSha(out)
}

if (phase == "record") {
    out <- file.path("val-out", "record", label)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    manifest <- msrval::valRecordInputs(draft, file.path(tmp_root, "val-data"), out, repo_root = ".")
    lib <- msrval::valLibraryTable()
    utils::write.table(lib, file.path(out, "library.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    msrval::valStatus(phase = "record", status = "completed", count = nrow(manifest))
    finish(out, manifest, "record")
} else if (phase == "pre") {
    out <- file.path("val-out", "pre", label)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    manifest_path <- optValue("--manifest", file.path("val-out", "record", label, "inputs.tsv"))
    truth_roles <- if (is.null(draft$truth)) character() else draft$role[toupper(draft$truth) == "TRUE"]
    msrval::valFetchInputs(manifest_path, file.path(tmp_root, "val-data-pre"), repo_root = ".",
                           truth_roles = truth_roles)
    work <- file.path(tmp_root, "val-pre-work")
    msrval::valRunPre(optValue("--script"), optValue("--pre-outputs"), work)
    file.copy(list.files(work, full.names = TRUE), out, recursive = TRUE)
    file.copy(optValue("--pre-outputs"), file.path(out, "pre-outputs.tsv"))
    manifest <- msrval::valReadManifest(manifest_path)
    spec <- utils::read.delim(optValue("--pre-outputs"), colClasses = "character",
                              na.strings = character())
    outputs <- lapply(unique(spec$file), function(f) {
        s <- spec[spec$file == f, ]
        lab <- s[s$type == "label", ]
        list(path = f, key = "none", max_rows = 1000L,
             lineage = manifest$role[manifest$kind == "data"],
             columns = lapply(seq_len(nrow(s)), function(i) list(name = s$column[[i]],
                                                                 type = s$type[[i]])),
             labels = stats::setNames(lapply(lab$labels, function(x) strsplit(x, ",")[[1L]]),
                                      lab$column))
    })
    msrval::valStatus(phase = "pre", status = "completed", count = length(list.files(work)))
    finish(out, manifest, "pre", outputs)
} else if (phase == "check") {
    manifest <- msrval::valReadManifest(optValue("--manifest"))
    runs <- file.path(evidence_dir, "runs.tsv")
    if (file.exists(runs)) msrval::valCheckRuns(utils::read.delim(runs, stringsAsFactors = FALSE))
    msrval::valCheckEvidence(evidence_dir, manifest, phase = optValue("--phase", "counted"))
    msrval::valStatus(phase = "check", status = "completed")
} else {
    stop(sprintf("phase %s is not yet wired to the driver core", phase), call. = FALSE)
}
