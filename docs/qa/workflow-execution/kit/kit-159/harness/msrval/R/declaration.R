# Declarations in two stages, their check, the harness pin and the environment (KIT-159 5.3, G15).

HARNESS_PATH <- "docs/qa/workflow-execution/kit/kit-159/harness"
DECLARATION_FIELDS <- c("schema", "ticket", "stage", "revision", "date", "release_commit",
                        "harness_tree", "environment", "inputs", "code", "datasets", "arms",
                        "metrics")
STAGE_FILE_PATTERNS <- c("^declaration\\.(md|json)$", "^inputs\\.tsv$", "^members-.*\\.tsv$",
                         "^design.*\\.tsv$", "^units-.*\\.tsv$", "^truth/", "^pre-outputs\\.tsv$")

#' Run git in a repository
#' @noRd
valGit <- function(repo, ...) {
    res <- processx::run("git", c("-C", repo, ...), error_on_status = FALSE)
    list(status = res$status, out = strsplit(trimws(res$stdout), "\n", fixed = TRUE)[[1L]])
}

#' Is commit `a` an ancestor of commit `b`?
#' @noRd
valIsAncestor <- function(repo, a, b) {
    valGit(repo, "merge-base", "--is-ancestor", a, b)$status == 0L
}

#' Text of a file at a commit, or NULL
#' @noRd
valGitShow <- function(repo, at, path) {
    res <- processx::run("git", c("-C", repo, "show", paste0(at, ":", path)),
                         error_on_status = FALSE)
    if (res$status != 0L) NULL else res$stdout
}

#' SHA-256 of a file at a commit
#' @noRd
valGitSha <- function(repo, at, path) {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    res <- processx::run("git", c("-C", repo, "cat-file", "blob", paste0(at, ":", path)),
                         stdout = tmp, error_on_status = FALSE)
    if (res$status != 0L) NA_character_ else valSha256(tmp)
}

#' Read a declaration and check its required fields and types (section 5.3)
#'
#' @param path Path of `declaration.json`.
#' @param text The declaration's text (in place of `path`).
#' @return The declaration as a list.
#' @export
valReadDeclaration <- function(path = NULL, text = NULL) {
    txt <- text %||% paste(readLines(path, warn = FALSE), collapse = "\n")
    decl <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                     error = function(e) {
                         valAbort("the declaration is not valid JSON", "val_spec_invalid")
                     })
    missing <- setdiff(DECLARATION_FIELDS, names(decl))
    if (length(missing) > 0L) {
        valAbort("the declaration lacks {paste(missing, collapse = ', ')}", "val_spec_invalid")
    }
    if (!identical(decl$schema, "msrval.declaration.v1")) {
        valAbort("the declaration schema is not msrval.declaration.v1", "val_spec_invalid")
    }
    if (!decl$stage %in% c("lifecycle", "counted")) {
        valAbort("the declaration stage is not lifecycle or counted", "val_spec_invalid")
    }
    hex40 <- function(x) is.character(x) && grepl("^[0-9a-f]{40}$", x)
    if (!hex40(decl$release_commit) || !hex40(decl$harness_tree)) {
        valAbort("release_commit and harness_tree must be 40 hexadecimal characters",
                 "val_spec_invalid")
    }
    valCheckDeclarationSchema(decl)
    decl
}

#' The pinned harness tree at a commit
#'
#' @param repo Repository path.
#' @param at Commit.
#' @return The tree hash of `docs/qa/workflow-execution/kit/kit-159/harness`.
#' @export
valHarnessTree <- function(repo, at = "HEAD") {
    valGit(repo, "rev-parse", paste0(at, ":", HARNESS_PATH))$out[[1L]]
}

#' Is the product tree at `at` the release commit's (docs, tickets and .github excluded)?
#'
#' @param repo Repository path.
#' @param release_commit The declared release commit.
#' @param at Commit.
#' @return `TRUE`, or a `val_product_tree_mismatch` refusal naming the files.
#' @export
valCheckProductTree <- function(repo, release_commit, at = "HEAD") {
    res <- valGit(repo, "diff", "--name-only", release_commit, at, "--", ".",
                  ":(exclude)docs", ":(exclude)tickets", ":(exclude).github")
    files <- res$out[nzchar(res$out)]
    if (res$status != 0L || length(files) > 0L) {
        valAbort("the product differs from the release commit: {paste(files, collapse = ', ')}",
                 "val_product_tree_mismatch")
    }
    TRUE
}

#' The job's library table (package, version; sorted; msrval left out)
#'
#' @return A data frame.
#' @export
valLibraryTable <- function() {
    ip <- utils::installed.packages(lib.loc = .libPaths())
    tab <- data.frame(package = ip[, "Package"], version = ip[, "Version"], row.names = NULL,
                      stringsAsFactors = FALSE)
    tab <- tab[!duplicated(tab$package) & tab$package != "msrval", , drop = FALSE]
    tab[order(tab$package, method = "radix"), , drop = FALSE]
}

#' SHA-256 of a library table (sorted by package)
#'
#' @param table Data frame with `package` and `version`.
#' @return The hash.
#' @export
valLibrarySha <- function(table) {
    t2 <- table[table$package != "msrval", , drop = FALSE]
    t2 <- t2[order(t2$package, method = "radix"), , drop = FALSE]
    digest::digest(paste(t2$package, t2$version, sep = "\t", collapse = "\n"), algo = "sha256",
                   serialize = FALSE)
}

#' G15: is the job's library the declared one?
#'
#' The table with the declared `extra_packages` closure removed must have the declared SHA-256, and
#' the job's additions must be the whole closure or none. Packages whose version differs from the
#' lockfile are listed (information only).
#'
#' @param library_table The job's library table.
#' @param declared_sha256 `environment.library_sha256`.
#' @param extra_packages The declared closure.
#' @param lockfile Path of `packaging/multischolar.lock.json` at the release commit, or `NULL`.
#' @return A list: `ok`, `sha256`, `extras` (`whole_closure` or `none`), `differs_from_lockfile`.
#' @export
valCheckEnvironment <- function(library_table, declared_sha256, extra_packages = character(),
                                lockfile = NULL) {
    lib <- library_table[library_table$package != "msrval", , drop = FALSE]
    present <- intersect(extra_packages, lib$package)
    if (length(present) > 0L && !setequal(present, extra_packages)) {
        valAbort("the job holds part of the declared extra closure", "val_environment_mismatch")
    }
    base <- lib[!lib$package %in% extra_packages, , drop = FALSE]
    sha <- valLibrarySha(base)
    if (!identical(sha, declared_sha256)) {
        valAbort("the library table differs from the declared one", "val_environment_mismatch")
    }
    differs <- NULL
    if (!is.null(lockfile)) {
        lock <- jsonlite::fromJSON(lockfile, simplifyVector = FALSE)$Packages
        lv <- vapply(lock, function(p) as.character(p$Version %||% NA), character(1))
        names(lv) <- vapply(lock, function(p) as.character(p$Package %||% NA), character(1))
        m <- base[base$package %in% names(lv), , drop = FALSE]
        m$lockfile <- lv[m$package]
        differs <- m[m$version != m$lockfile, , drop = FALSE]
    }
    list(ok = TRUE, sha256 = sha,
         extras = if (length(present) > 0L) "whole_closure" else "none",
         differs_from_lockfile = differs)
}

#' Functions defined at the top level of R files (for name resolution)
#' @noRd
valDefinedFunctions <- function(texts) {
    unlist(lapply(texts, function(txt) {
        exprs <- tryCatch(parse(text = txt, keep.source = FALSE), error = function(e) NULL)
        unlist(lapply(exprs, function(e) {
            if (is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") && is.name(e[[2L]]) &&
                    is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) {
                as.character(e[[2L]])
            }
        }))
    }))
}

#' The paths a declaration commit may touch
#' @noRd
valStageAllowed <- function(touched, stage_dir, code_paths) {
    inside <- startsWith(touched, paste0(stage_dir, "/"))
    rel <- substring(touched, nchar(stage_dir) + 2L)
    pattern_ok <- vapply(rel, function(r) {
        any(vapply(STAGE_FILE_PATTERNS, grepl, logical(1), x = r))
    }, logical(1))
    inside & (pattern_ok | rel %in% code_paths)
}

#' Check a declaration stage at a commit (section 5.3)
#'
#' Files exist; committed alone; frozen; hashes; order; harness pinned; product is the release;
#' environment (with a library table); attribution and declared outputs; for `counted`, every path
#' whose hash differs from the lifecycle declaration listed in `changes_from_lifecycle`. Phase
#' `amended` accepts exactly the replacements `declaration-amendments.json` names.
#'
#' @param repo Repository path (full history).
#' @param evidence_dir The ticket's evidence folder, relative to the repository root.
#' @param stage `lifecycle` or `counted`.
#' @param at Commit under test.
#' @param phase The phase (`amended` reads the amendments).
#' @param library_table The job's library table, or `NULL` to skip G15.
#' @param units Unit ids for output placeholder expansion.
#' @return A list: `ok`, `declaration`, `declaration_commit`, `at`, `amended`, `declaration_diff`.
#' @export
valCheckDeclaration <- function(repo, evidence_dir, stage, at = "HEAD", phase = stage,
                                library_table = NULL, units = character()) {
    sd <- file.path(evidence_dir, stage)
    json_path <- file.path(sd, "declaration.json")
    for (p in c(json_path, file.path(sd, "declaration.md"))) {
        if (valGit(repo, "cat-file", "-e", paste0(at, ":", p))$status != 0L) {
            valAbort("{p} is not committed at {at}", "val_declaration_missing")
        }
    }
    at_commit <- valGit(repo, "rev-parse", at)$out[[1L]]
    decl <- valReadDeclaration(text = valGitShow(repo, at, json_path))
    dc <- utils::tail(valGit(repo, "log", "--format=%H", "--diff-filter=A", at, "--",
                             json_path)$out, 1L)
    code_paths <- vapply(decl$code, function(x) x$path, character(1))
    touched <- valGit(repo, "diff-tree", "--no-commit-id", "--name-only", "-r", "--root", dc)$out
    touched <- touched[nzchar(touched)]
    if (!all(valStageAllowed(touched, sd, code_paths))) {
        valAbort("the declaration commit touches other files",
                 "val_declaration_not_committed_alone")
    }
    amend <- if (identical(phase, "amended")) valReadAmendments(repo, at, sd) else NULL
    replaced <- amend[["code"]] %||% list()
    replaced_sha <- stats::setNames(vapply(replaced, function(x) x$sha256, character(1)),
                                    vapply(replaced, function(x) x$path, character(1)))
    valCheckFrozen(repo, dc, at, touched, sd, replaced_sha)
    valCheckHashes(repo, at, sd, evidence_dir, decl, dc, replaced_sha)
    if (!valIsAncestor(repo, dc, at)) {
        valAbort("the declaration commit is not an ancestor of {at}", "val_declaration_after_run")
    }
    diff <- if (stage == "counted") valCheckCountedOrder(repo, evidence_dir, decl, dc) else NULL
    expected_tree <- amend[["harness_tree"]] %||% decl[["harness_tree"]]
    if (!identical(valHarnessTree(repo, at), expected_tree)) {
        valAbort("the harness tree at {at} is not the declared one", "val_harness_mismatch")
    }
    valCheckProductTree(repo, decl$release_commit, at)
    if (!is.null(library_table)) {
        valCheckEnvironment(library_table, decl$environment$library_sha256,
                            as.character(unlist(decl$environment$extra_packages)))
    }
    valCheckDeclaredTables(repo, at, sd, decl, units)
    if (length(units)) valCheckDeclarationSchema(decl, units)
    list(ok = TRUE, declaration = decl, declaration_commit = dc, at = at_commit,
         amended = amend[["label"]], declaration_diff = diff)
}

#' Read `declaration-amendments.json` at a commit
#' @noRd
valReadAmendments <- function(repo, at, sd) {
    txt <- valGitShow(repo, at, file.path(sd, "declaration-amendments.json"))
    if (is.null(txt)) valAbort("phase amended without declaration-amendments.json",
                               "val_declaration_changed")
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

#' Frozen: every path the declaration commit touched is byte-identical at `at`, except listed
#' amendment replacements with their listed hash
#' @noRd
valCheckFrozen <- function(repo, dc, at, touched, sd, replaced_sha) {
    for (p in touched) {
        if (valGit(repo, "diff", "--quiet", dc, at, "--", p)$status == 0L) next
        rel <- substring(p, nchar(sd) + 2L)
        listed <- rel %in% names(replaced_sha)
        if (listed && identical(valGitSha(repo, at, p), replaced_sha[[rel]])) next
        valAbort("{p} changed after the declaration", "val_declaration_changed")
    }
    invisible(TRUE)
}

#' Hashes of inputs.tsv, code and preparation files; names resolve
#' @noRd
valCheckHashes <- function(repo, at, sd, evidence_dir, decl, dc, replaced_sha) {
    changed <- function(what) {
        valAbort("{what} differs from the declaration", "val_declaration_changed")
    }
    if (!identical(valGitSha(repo, at, file.path(sd, decl$inputs$path)), decl$inputs$sha256)) {
        changed("inputs.tsv")
    }
    for (x in decl$code) {
        want <- if (x$path %in% names(replaced_sha)) replaced_sha[[x$path]] else x$sha256
        if (!identical(valGitSha(repo, at, file.path(sd, x$path)), want)) changed(x$path)
    }
    prep <- decl[["preparation"]]
    if (!is.null(prep)) {
        if (!valIsAncestor(repo, prep$commit, dc)) changed("the preparation commit's order")
        for (x in prep$files) {
            if (!identical(valGitSha(repo, at, file.path(evidence_dir, x$path)), x$sha256)) {
                changed(x$path)
            }
        }
    }
    r_files <- grep("\\.R$", vapply(decl$code, function(x) x$path, character(1)), value = TRUE)
    code_text <- lapply(r_files, function(p) valGitShow(repo, at, file.path(sd, p)))
    defined <- c(valDefinedFunctions(code_text), getNamespaceExports("msrval"))
    wanted <- c(vapply(decl$metrics, function(m) m[["function"]] %||% NA_character_, character(1)),
                vapply(decl$arms, function(a) {
                    a[["function"]] %||% a[["builder"]] %||% NA_character_
                }, character(1)),
                unlist(decl[["hooks"]] %||% list()))
    wanted <- wanted[!is.na(wanted) & is.character(wanted) & !wanted %in% c("TRUE", "FALSE")]
    if (!all(wanted %in% defined)) {
        valAbort("a declared function does not resolve: {setdiff(wanted, defined)}",
                 "val_spec_invalid")
    }
    invisible(TRUE)
}

#' Counted stage: the lifecycle declaration and its provenance precede it; no counted provenance
#' precedes it; every changed path is listed
#' @noRd
valCheckCountedOrder <- function(repo, evidence_dir, decl, dc) {
    after <- function(what) valAbort(what, "val_declaration_after_run")
    life_json <- valGitShow(repo, dc, file.path(evidence_dir, "lifecycle", "declaration.json"))
    if (is.null(life_json)) after("the lifecycle declaration does not precede the counted one")
    runs <- valGit(repo, "ls-tree", "-r", "--name-only", dc, "--",
                   file.path(evidence_dir, "runs"))$out
    provs <- grep("/provenance\\.json$", runs, value = TRUE)
    stages <- vapply(provs, function(p) {
        x <- jsonlite::fromJSON(valGitShow(repo, dc, p), simplifyVector = FALSE)
        as.character(x[["stage"]] %||% "")
    }, character(1))
    if (!any(stages == "lifecycle")) {
        after("no lifecycle provenance precedes the counted declaration")
    }
    if (any(stages == "counted")) after("a counted provenance precedes the counted declaration")
    life <- jsonlite::fromJSON(life_json, simplifyVector = FALSE)
    hashes <- function(d) {
        stats::setNames(vapply(d$code, function(x) x$sha256, character(1)),
                        vapply(d$code, function(x) x$path, character(1)))
    }
    lh <- c(hashes(life), inputs.tsv = life$inputs$sha256)
    ch <- c(hashes(decl), inputs.tsv = decl$inputs$sha256)
    paths <- union(names(lh), names(ch))
    differ <- paths[is.na(lh[paths]) | is.na(ch[paths]) | lh[paths] != ch[paths]]
    listed <- vapply(decl[["changes_from_lifecycle"]] %||% list(), function(x) x[["path"]],
                     character(1))
    if (!all(differ %in% listed)) {
        valAbort("changed from the lifecycle and not listed: {setdiff(differ, listed)}",
                 "val_declaration_changed")
    }
    data.frame(path = differ, lifecycle = unname(lh[differ]), counted = unname(ch[differ]),
               listed = differ %in% listed, stringsAsFactors = FALSE)
}

#' Attribution (every inputs.tsv row cited) and declared result tables
#' @noRd
valCheckDeclaredTables <- function(repo, at, sd, decl, units) {
    txt <- valGitShow(repo, at, file.path(sd, decl$inputs$path))
    tab <- utils::read.delim(text = txt, colClasses = "character", na.strings = character(),
                             stringsAsFactors = FALSE)
    if (nrow(tab) > 0L && (!"citation" %in% names(tab) || any(!nzchar(trimws(tab$citation))))) {
        valAbort("an inputs.tsv row has no citation", "val_spec_invalid")
    }
    arms <- vapply(decl[["arms"]], function(a) a[["id"]] %||% NA_character_, character(1))
    for (o in decl[["outputs"]] %||% list()) {
        valExpandOutputPaths(o, arms = arms, units = units)
        if (identical(o[["key"]] %||% "none", "none")) {
            types <- vapply(o$columns, function(x) x$type, character(1))
            nm <- vapply(o$columns, function(x) x$name, character(1))
            if (!all(types %in% RESULT_TYPES) || (o[["max_rows"]] %||% 1001L) > 1000L ||
                    !all(nm[types == "label"] %in% names(o[["labels"]] %||% list()))) {
                valAbort(paste("output {o$path}: a result table needs typed columns, labels",
                               "and at most 1,000 rows"), "val_spec_invalid")
            }
        }
    }
    invisible(TRUE)
}

#' A correction runs at the counted run's package commit (section 5.3)
#'
#' @param repo Repository path.
#' @param evidence_dir Evidence folder.
#' @param label The counted run's label.
#' @param at The correction's commit.
#' @return `TRUE`, or `val_declaration_changed`.
#' @export
valCheckCorrection <- function(repo, evidence_dir, label, at = "HEAD") {
    prov_path <- file.path(repo, evidence_dir, "runs", label, "provenance.json")
    if (!file.exists(prov_path)) {
        valAbort("no committed provenance for the counted run {label}", "val_declaration_changed")
    }
    prov <- jsonlite::fromJSON(prov_path, simplifyVector = FALSE)
    at_commit <- valGit(repo, "rev-parse", at)$out[[1L]]
    if (!identical(prov$code$package_commit, at_commit)) {
        valAbort("a correction at another commit than the counted run's", "val_declaration_changed")
    }
    TRUE
}

#' The annotation exemption (section 5.1): `expected = "none"` only when at most 1 % of the
#' identifiers match the package's UniProt accession pattern
#'
#' @param ids The data set's protein identifiers.
#' @param annotation The declaration's `annotation` block.
#' @param pattern The accession pattern (default: `MultiScholaR::getUniprotRegexPatterns()$entry`).
#' @return `TRUE`, or `val_spec_invalid`.
#' @export
valCheckAnnotationExemption <- function(ids, annotation, pattern = NULL) {
    if (!identical(annotation$expected, "none")) return(TRUE)
    if (is.null(annotation$reason) || !nzchar(annotation$reason)) {
        valAbort("annotation 'none' needs a reason", "val_spec_invalid")
    }
    pattern <- pattern %||% getExportedValue("MultiScholaR", "getUniprotRegexPatterns")()$entry
    share <- mean(grepl(pattern, ids))
    if (share > 0.01) {
        valAbort("annotation 'none' on identifiers {round(100 * share, 1)} % UniProt-like",
                 "val_spec_invalid")
    }
    TRUE
}

#' At most one correction per counted run in `runs.tsv` (A23)
#'
#' @param runs `runs.tsv` as a data frame (`label`, `phase`).
#' @return `TRUE`, or `val_spec_invalid`.
#' @export
valCheckRuns <- function(runs) {
    corr <- runs[runs$phase == "correction", , drop = FALSE]
    if (anyDuplicated(corr$label)) {
        valAbort("runs.tsv holds more than one correction of a counted run", "val_spec_invalid")
    }
    TRUE
}
