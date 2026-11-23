# KIT-159 entry script (section 7.8). From the package root, in a `run-tests` scripts request:
#
#   Rscript docs/qa/workflow-execution/kit/kit-159/harness/val-run.R record <evidence dir> <label>
#       [--uniprot <accession list>]
#   Rscript .../val-run.R pre <evidence dir> <label> --script <path> --pre-outputs <path>
#   Rscript .../val-run.R check <folder> --manifest <inputs.tsv>
#   Rscript .../val-run.R s0 <evidence dir> <label> --tests <harness-tests.tsv>
#       --test-files <a,b> [--cells <c1,c2>]
#
# s0 (10.4, r1.14): the evidence run on generated inputs -> val-out/s0/<label>/ (corpus.tsv,
# routes.tsv, steps.tsv, modal.tsv, alternatives.tsv, ledgers/, run.tsv, outputs.sha256); refused
# (exit 3) before any route unless the named test files of the job's first script all passed.
# record: draft-inputs.tsv of <evidence dir> -> val-out/record/<label>/ (inputs.tsv, members-*.tsv,
# library.tsv, attribution.tsv, provenance.json, outputs.sha256); no declaration check.
# pre: the ticket's pre-declaration script on the recorded inputs (val-out/record/<label>/inputs.tsv
# of the same job) -> val-out/pre/<label>/ (the listed typed tables only).
# check: the allow-list and hashes on a committed folder (valCheckEvidence, runs.tsv corrections).
#   Rscript .../val-run.R <lifecycle|counted|correction|determinism|amended> <evidence dir> <stage>
#       <label> --dataset <id> [--units <first:last>] [--arms <a,b>] [--chunks <k>]
#   Rscript .../val-run.R score <evidence dir> <stage> <label> --dataset <id> --runs <l1,l2,...>
# The cell phases (7.8 steps 4 to 9): declaration check for the stage, environment check (G15),
# retrieval, then msrval::valRunPhase (prepare.R, kernel gate, cells in chunks, blind scoring in
# the lifecycle) into val-out/runs/<label>/ (amended: <label>-amended). score: the committed
# runs/<label>/ folders of the named requests combined (valScoreRuns) and scored at run level.
# parity runs in the bundle jobs (val-parity.R, 7.9), not here.

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
cell_phases <- c("lifecycle", "counted", "correction", "determinism", "amended", "score")
draft_path <- if (phase %in% cell_phases && length(opts) > 0L) {
    file.path(evidence_dir, opts[[1L]], "inputs.tsv")
} else {
    file.path(evidence_dir, "draft-inputs.tsv")
}
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
    # R1 (r1.12): a UniProt cache recorded from a committed accession list; only uniprot.tsv
    # (release and hash) enters the record folder, the cache stays private
    if (!is.null(optValue("--uniprot"))) {
        msrval::valRecordUniprot(optValue("--uniprot"), out, file.path(tmp_root, "val-uniprot"))
    }
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
} else if (phase == "s0") {
    # S0 (10.4, r1.14): the evidence run on generated inputs, gating nothing (VAL-41). The job's
    # first script ran the new pieces' tests; unless they all passed, S0 ends here before any
    # corpus or route (VAL-54 instrument fault).
    out <- file.path("val-out", "s0", label)
    cells <- optValue("--cells")
    res <- tryCatch(
        msrval::valRunS0(out, ".", list(path = optValue("--tests"),
                                         files = strsplit(optValue("--test-files", ""), ",")[[1L]]),
                         file.path(tmp_root, "val-s0"),
                         cells = if (is.null(cells)) NULL else strsplit(cells, ",")[[1L]]),
        val_error = function(e) e
    )
    if (inherits(res, "val_error")) {
        msrval::valStatus(phase = "s0", status = "refused", class = class(res)[[1L]])
        quit(status = 3L)
    }
    run <- data.frame(field = c("package_commit", "harness_tree", "started", "finished"),
                      value = c(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                                msrval::valHarnessTree("."), started,
                                format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
    utils::write.table(run, file.path(out, "run.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    msrval::valWriteOutputsSha(out)
    msrval::valStatus(phase = "s0", status = "completed", count = nrow(res))
} else if (phase == "check") {
    manifest <- msrval::valReadManifest(optValue("--manifest"))
    runs <- file.path(evidence_dir, "runs.tsv")
    if (file.exists(runs)) msrval::valCheckRuns(utils::read.delim(runs, stringsAsFactors = FALSE))
    msrval::valCheckEvidence(evidence_dir, manifest, phase = optValue("--phase", "counted"))
    msrval::valStatus(phase = "check", status = "completed")
} else if (phase %in% cell_phases) {
    stage <- opts[[1L]]
    label <- opts[[2L]]
    dataset <- optValue("--dataset")
    if (is.null(dataset)) stop("--dataset is required", call. = FALSE)
    sd <- file.path(evidence_dir, stage)
    decl <- jsonlite::read_json(file.path(sd, "declaration.json"), simplifyVector = FALSE)
    ds <- Filter(function(d) identical(d$id, dataset), decl$datasets)[[1L]]
    units <- msrval::valUnits(file.path(sd, ds$units$table))
    lib <- msrval::valLibraryTable()
    # Step 4: the declaration check and the environment check (G15), before any download.
    msrval::valCheckDeclaration(".", evidence_dir, stage, at = "HEAD", phase = phase,
                                library_table = lib, units = units$unit)
    env_check <- msrval::valCheckEnvironment(lib, decl$environment$library_sha256,
                                             extra_packages = unlist(decl$environment$extra_packages))
    manifest_path <- file.path(sd, decl$inputs$path)
    manifest <- msrval::valReadManifest(manifest_path)
    commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    decl_commit <- utils::tail(system2("git", c("log", "--format=%H", "--diff-filter=A", "HEAD",
                                                "--", file.path(sd, "declaration.json")),
                                       stdout = TRUE), 1L)
    code <- list(package_commit = commit, release_commit = decl$release_commit,
                 harness = list(commit = commit, tree = msrval::valHarnessTree(".")),
                 declaration = list(stage = stage, commit = decl_commit),
                 environment_check = env_check)
    if (phase == "score") {
        labels <- strsplit(optValue("--runs", ""), ",")[[1L]]
        cells <- msrval::valPhaseCells(decl, dataset, units, "counted")
        got <- msrval::valScoreRuns(evidence_dir, labels, cells[, c("unit", "arm")])
        metrics <- got$metrics
        metrics$value <- as.numeric(metrics$value)
        out <- file.path("val-out", "runs", label)
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        cmp <- msrval::valRunLevel(decl, metrics, units, phase = stage)
        utils::write.table(cmp, file.path(out, "comparisons.tsv"), sep = "\t", quote = FALSE,
                           row.names = FALSE)
        utils::write.table(manifest[, c("role", "kind", "licence", "redistribution", "citation")],
                           file.path(out, "attribution.tsv"), sep = "\t", quote = FALSE,
                           row.names = FALSE)
        prov <- msrval::valProvenance(ticket = decl$ticket, stage = stage, label = label,
                                      phase = "score", status = "completed",
                                      blinded = identical(stage, "lifecycle"), started = started,
                                      code = code, scoring = list(requests = labels))
        msrval::valWriteProvenance(prov, file.path(out, "provenance.json"))
        msrval::valCheckProvenance(file.path(out, "provenance.json"))
        msrval::valCheckEvidence(out, manifest, phase = stage)
        msrval::valWriteOutputsSha(out)
    } else {
        # Step 5: retrieval into the private folder, then the phase (steps 5 to 9).
        data_dir <- file.path(tmp_root, "val-data")
        msrval::valFetchInputs(manifest_path, data_dir, repo_root = ".")
        range <- optValue("--units")
        unit_ids <- if (is.null(range)) NULL else {
            ends <- strsplit(range, ":", fixed = TRUE)[[1L]]
            i <- match(ends, units$unit)
            if (anyNA(i)) stop("--units names a unit not in the unit table", call. = FALSE)
            units$unit[i[[1L]]:i[[length(i)]]]
        }
        arms <- optValue("--arms")
        work <- file.path(tmp_root, "val-work")
        dir.create(work, recursive = TRUE, showWarnings = FALSE)
        msrval::valRunPhase(
            phase = phase, decl = decl, dataset = dataset, units = units,
            label = if (phase == "amended") paste0(label, "-amended") else label,
            code_files = file.path(sd, vapply(decl$code, function(x) x$path, character(1))),
            data_dir = data_dir, work_dir = normalizePath(work), out_dir = "val-out",
            chunks = as.integer(optValue("--chunks", "4")),
            arms = if (is.null(arms)) NULL else strsplit(arms, ",")[[1L]], unit_ids = unit_ids,
            manifest = manifest, code = code, pkg_root = "."
        )
    }
    msrval::valStatus(phase = phase, status = "completed")
} else {
    stop(sprintf("unknown phase %s", phase), call. = FALSE)
}
