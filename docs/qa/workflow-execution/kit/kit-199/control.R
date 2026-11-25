# KIT-199 script S1: the control and the freeze of the GUI-entry reference prot_dia_limpa (ticket section 6.2).
#
# Usage (run-tests scripts job, from the package root of the script commit F1, full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-199/control.R <A> <B> <out>     control mode (R1, R2)
#   Rscript docs/qa/workflow-execution/kit/kit-199/control.R measure <out>     measurement mode (R4)
# Internal (one child process per side, working directory at that side's tree root):
#   Rscript --vanilla control.R child <side> <private dir> <run|freeze> <steps file>
#
# Control mode. Three sides, one after another, each in a detached worktree $RUNNER_TEMP/k199/<side>:
#   anchor = A, base = B: inside testthat::test_that(), guiEntryRunLane("prot_dia_limpa") with the DA copy hook
#     (G1D's: at summary_refusal the DA table is copied as step `da` left it, its masked digest computed as the lane
#     computes it, then the step's own actions run), MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=<out>/<side>.
#   head = this script's commit (F1): inside testthat::test_that(), guiEntryExpectLane("prot_dia_limpa") with
#     MULTISCHOLAR_GUI_ENTRY_FREEZE=1, MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR=<out>/oracles and the record directory
#     <out>/head; the same DA copy is made by rebinding guiEntryRunLane, in the environment the helpers were sourced
#     into, to a wrapper that calls the sourced original with the hook. No file of the package is edited.
# Each child loads its tree as the test jobs do (pkgload::load_all with helpers, then the helpers sourced as
# testthat::test_file sources them) and writes <out>/steps-<side>.json (guiEntryComparable(record) written by
# guiEntryWriteJson). The parent then writes, under <out>: steps-frozen.json (the frozen file's `steps` member written
# alone), sha256.tsv, gui-entry-control.tsv (compareLane's logic, parsed from packaging/verify/run-lanes.R),
# verdict.tsv, da-runner.tsv, environment.tsv, negative-control.tsv and status.tsv.
# Exit status 0 only when every verdict check holds, the negative control behaves and no stage is `mechanism`.
#
# Measurement mode (ticket 6.2 R4). One side, `measure`, on this script's own commit in a detached worktree, run as
# the anchor and base sides run (inside test_that(), with the DA copy; no freeze, no comparison with the stored file).
# Writes da-runner.tsv, environment.tsv, measure.tsv (the masked DA digest), sha256.tsv, verdict.tsv and status.tsv.
#
# Private material (the raw DA table copies, child consoles, worktrees) stays under $RUNNER_TEMP; nothing is written
# into tests/. Lines written for upload avoid the upload scan's patterns (.github/run-tests/scan.R).

LANE <- "prot_dia_limpa"
DA_KEY <- "results/proteomics/da_proteins/da_proteins_groupKO-groupWT_long_annot.tsv"
STORED_DIGEST <- "80cfedaae858a2e93b87185583156c17b106b157ac02e3eacfc068edf7679b47"
RUNNER_DIGEST <- "8a4c2dd50ad521daf9b407f9d3beb2ba5bff2b3bd393b3a48fdf13e7bcd112b9"
DECLARED_CLASS <- "declared_change (KIT-199)"
INTEGRATION_CLASS <- "integration_artifacts"
CHANGED_STEP <- "da"
HOOK_STEP <- "summary_refusal"
STEP_COUNT <- 24L
CHILD_TIMEOUT_SECONDS <- 1200
RECORDS_DIR <- "docs/qa/workflow-execution/kit/kit-199"
NUMERIC_EXACT <- c("fdr_qvalue", "raw_pvalue", "log2FC", "annotation_score")
NUMERIC_PREFIXES <- c("^log2norm\\.", "^raw\\.")
UPLOAD_PATTERNS <- c(
    "\\.rs:[0-9]", "src/rust|\\.cargo|Cargo\\.(toml|lock)", "rustc|cargo |panicked at",
    "Compiling [a-z_]+ v[0-9]|error\\[E[0-9]+\\]|warning: unused",
    "(^|[ /])(gcc|g\\+\\+|clang|cc1|ld)( |:)|-fPIC|-std=gnu", "make(\\[[0-9]+\\])?: \\*\\*\\*",
    "\\.(c|cc|cpp|h|hpp|o):[0-9]+", "fast_limpa_dpc|fast_duckdb|internal-build-logs",
    "gh[pousr]_[A-Za-z0-9]{20,}|github_pat_|Bearer "
)

oneLine <- function(text) gsub("[\t\r\n]+", " ", paste(as.character(text), collapse = " "))

cleanDetail <- function(text) {
    text <- substr(oneLine(text), 1L, 300L)
    hit <- vapply(UPLOAD_PATTERNS, \(pattern) grepl(pattern, text, useBytes = TRUE), logical(1L))
    if (any(hit)) "<detail not written: it matches an upload pattern>" else text
}

formatSeconds <- function(seconds) {
    if (length(seconds) != 1L || is.na(seconds)) "" else sprintf("%.1f", seconds)
}

errorDetail <- function(error) paste0(class(error)[[1L]], ": ", conditionMessage(error))

attempt <- function(expr) {
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(list(ok = TRUE, value = expr), error = \(e) list(ok = FALSE, error = e))
    result$seconds <- proc.time()[["elapsed"]] - started
    result
}

writeTsv <- function(table, path) {
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

sha256File <- function(path) {
    if (file.exists(path)) digest::digest(file = path, algo = "sha256") else "<absent>"
}

# --------------------------------------------------------------------------------------------------------------------
# Child: one side, working directory at that side's tree root.
# --------------------------------------------------------------------------------------------------------------------

childMain <- function(side, private_dir, mode, steps_path) {
    stages_path <- file.path(private_dir, "child-stages.tsv")
    writeStage <- function(stage, outcome, detail, seconds) {
        cat(paste(stage, outcome, oneLine(detail), formatSeconds(seconds), sep = "\t"), "\n",
            sep = "", file = stages_path, append = TRUE)
    }
    loaded <- attempt({
        suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
        suppressPackageStartupMessages(library(testthat))
        # testthat::test_file() sources the helpers again into the test environment (a child of the global
        # environment), with TESTTHAT set and the working directory at the test directory.
        helper_env <- new.env(parent = globalenv())
        test_dir <- normalizePath("tests/testthat")
        withr::with_dir(test_dir, withr::with_envvar(c(TESTTHAT = "true"),
            suppressMessages(testthat::source_test_helpers(test_dir, env = helper_env))))
        helper_env
    })
    if (!loaded$ok) {
        writeStage("load", "error", errorDetail(loaded$error), loaded$seconds)
        return(1L)
    }
    writeStage("load", "ok", "", loaded$seconds)
    helper_env <- loaded$value
    run_lane <- get("guiEntryRunLane", envir = helper_env, mode = "function")
    expect_lane <- get("guiEntryExpectLane", envir = helper_env, mode = "function")
    mask <- get("guiEntryMask", envir = helper_env, mode = "function")
    comparable <- get("guiEntryComparable", envir = helper_env, mode = "function")
    write_json <- get("guiEntryWriteJson", envir = helper_env, mode = "function")

    # G1D's copy hook (kit-159/g1d/g1d.R:119-142): run in place of summary_refusal's actions.
    copy_path <- file.path(private_dir, "da-table.tsv")
    hookFn <- function() {
        caller <- parent.frame()
        step <- get("step", envir = caller, mode = "list")
        api <- get("api", envir = caller, mode = "list")
        if (!identical(step$id, HOOK_STEP)) stop("the copy hook ran at an unexpected step")
        roots <- NULL
        for (frame in rev(sys.frames())) {
            if (exists("roots", envir = frame, inherits = FALSE) &&
                exists("acc", envir = frame, inherits = FALSE)) {
                roots <- get("roots", envir = frame, inherits = FALSE)
                break
            }
        }
        if (is.null(roots)) stop("the lane's masking roots were not found")
        base_dir <- get("project_dirs", envir = globalenv())[["proteomics"]][["base_dir"]]
        source_path <- file.path(base_dir, DA_KEY)
        if (!file.exists(source_path)) stop("the DA table is absent after step da")
        if (!file.copy(source_path, copy_path, overwrite = TRUE)) stop("the DA table could not be copied")
        masked <- mask(readLines(source_path, warn = FALSE), roots)
        writeLines(digest::digest(masked, algo = "sha256", serialize = FALSE),
            file.path(private_dir, "da-table-digest.txt"))
        step$run(api)
    }
    hook <- list(step = HOOK_STEP, fn = hookFn)

    captured <- NULL
    test_error <- ""
    started <- proc.time()[["elapsed"]]
    if (identical(mode, "freeze")) {
        # guiEntryExpectLane has no `fault` argument; it resolves guiEntryRunLane in helper_env.
        assign("guiEntryRunLane", function(lane_id, fault = NULL) {
            record <- run_lane(lane_id, fault = hook)
            captured <<- record
            record
        }, envir = helper_env)
        passed <- tryCatch(
            testthat::test_that(paste("KIT-199 control", side, LANE), expect_lane(LANE)),
            error = \(e) {
                test_error <<- errorDetail(e)
                FALSE
            }
        )
    } else {
        passed <- tryCatch(
            testthat::test_that(paste("KIT-199 control", side, LANE), {
                captured <<- run_lane(LANE, fault = hook)
                testthat::succeed()
            }),
            error = \(e) {
                test_error <<- errorDetail(e)
                FALSE
            }
        )
    }
    seconds <- proc.time()[["elapsed"]] - started
    if (is.null(captured)) {
        writeStage("lane", "error", paste("the lane did not complete;", test_error), seconds)
        return(1L)
    }
    write_json(comparable(captured), steps_path)
    if (!isTRUE(passed)) {
        writeStage("lane", "test_failed", paste("test_that() did not pass;", test_error), seconds)
        return(1L)
    }
    writeStage("lane", "ok", sprintf("lane record %s s", formatSeconds(as.numeric(captured$seconds))), seconds)
    0L
}

# --------------------------------------------------------------------------------------------------------------------
# Parent: sides.
# --------------------------------------------------------------------------------------------------------------------

newStatus <- function(path) {
    writeLines(paste("stage", "side", "outcome", "detail", "seconds", sep = "\t"), path)
    state <- new.env(parent = emptyenv())
    state$mechanism <- FALSE
    state$emit <- function(stage, side, outcome, detail = "", seconds = NA_real_) {
        if (identical(outcome, "mechanism")) state$mechanism <- TRUE
        cat(paste(stage, side, outcome, cleanDetail(detail), formatSeconds(seconds), sep = "\t"), "\n",
            sep = "", file = path, append = TRUE)
        cat(sprintf("%s %s: %s\n", side, stage, outcome))
    }
    state
}

runGit <- function(args, wd = NULL) {
    processx::run("git", args, wd = wd, error_on_status = FALSE)
}

prepareTree <- function(side, commit, root, tree, private_dir, status) {
    made <- attempt(runGit(c("-C", root, "worktree", "add", "--detach", tree, commit)))
    if (!made$ok || made$value$status != 0L) {
        status$emit("worktree", side, "mechanism", "git worktree add failed", made$seconds)
        return(FALSE)
    }
    status$emit("worktree", side, "ok", "", made$seconds)
    head <- attempt(trimws(runGit(c("-C", tree, "rev-parse", "HEAD"))$stdout))
    if (!head$ok || !identical(head$value, commit)) {
        status$emit("head_commit", side, "mechanism", "HEAD differs from the commit", head$seconds)
        return(FALSE)
    }
    status$emit("head_commit", side, "ok", commit, head$seconds)
    # The workflow's native-code step for a package tree (run-tests.yml, job scripts, "Install package_commit").
    manifest <- file.path(tree, "inst", "duckdb", "linux_amd64", "manifest.rds")
    if (!file.exists(file.path(tree, "src", "build_dia_rsum.R")) || file.exists(manifest)) {
        status$emit("native", side, "ok", "not needed")
        return(TRUE)
    }
    native <- attempt(processx::run(file.path(R.home("bin"), "Rscript"), c("--vanilla", "src/build_dia_rsum.R"),
        wd = tree, error_on_status = FALSE, stdout = file.path(private_dir, "native.log"), stderr = "2>&1"))
    if (!native$ok || native$value$status != 0L || !file.exists(manifest)) {
        status$emit("native", side, "mechanism", "the native step did not complete", native$seconds)
        return(FALSE)
    }
    status$emit("native", side, "ok", "", native$seconds)
    TRUE
}

# One side: worktree, native step, child. Returns what the parent needs (paths and outcomes).
runSide <- function(side, commit, mode, ctx, status) {
    private_dir <- file.path(ctx$private_root, side)
    dir.create(private_dir, recursive = TRUE, showWarnings = FALSE)
    tree <- file.path(ctx$trees_root, side)
    result <- list(side = side, commit = commit, tree = tree, lane_ok = FALSE, copy_ok = FALSE,
        copy_path = file.path(private_dir, "da-table.tsv"),
        steps_path = file.path(ctx$out, paste0("steps-", side, ".json")))
    if (!prepareTree(side, commit, ctx$root, tree, private_dir, status)) return(result)
    env <- c("current", MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = file.path(ctx$out, side))
    if (identical(mode, "freeze")) {
        env <- c(env, MULTISCHOLAR_GUI_ENTRY_FREEZE = "1",
            MULTISCHOLAR_GUI_ENTRY_ORACLE_DIR = file.path(ctx$out, "oracles"))
    }
    child <- attempt(processx::run(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", ctx$script_path, "child", side, private_dir, mode, result$steps_path), wd = tree,
        env = env, error_on_status = FALSE, stdout = file.path(private_dir, "child.log"), stderr = "2>&1",
        timeout = CHILD_TIMEOUT_SECONDS, cleanup_tree = TRUE))
    stages_path <- file.path(private_dir, "child-stages.tsv")
    if (file.exists(stages_path)) {
        for (line in readLines(stages_path, warn = FALSE)) {
            fields <- c(strsplit(line, "\t", fixed = TRUE)[[1L]], "", "", "", "")[1:4]
            status$emit(fields[[1L]], side, fields[[2L]], fields[[3L]], suppressWarnings(as.numeric(fields[[4L]])))
        }
    }
    if (!child$ok || isTRUE(child$value$timeout)) {
        status$emit("child", side, "mechanism", if (child$ok) "timeout" else errorDetail(child$error), child$seconds)
        return(result)
    }
    status$emit("child", side, if (child$value$status == 0L) "ok" else "result",
        sprintf("the child exited %s", child$value$status), child$seconds)
    result$lane_ok <- identical(as.integer(child$value$status), 0L)
    if (!result$lane_ok) {
        status$emit("record", side, "not_checked", "the lane did not pass")
        return(result)
    }

    record_path <- file.path(ctx$out, side, paste0(LANE, ".json"))
    record <- attempt(jsonlite::read_json(record_path, simplifyVector = FALSE))
    if (!record$ok) {
        status$emit("record", side, "mechanism", "the lane's record is not readable")
        return(result)
    }
    ids <- names(record$value$steps)
    status$emit("record", side, if (length(ids) == STEP_COUNT) "ok" else "result",
        sprintf("%d step ids", length(ids)))
    result$record_digest <- record$value$steps$da$artifacts[[DA_KEY]] %||% "<absent>"
    digest_path <- file.path(private_dir, "da-table-digest.txt")
    copied <- if (file.exists(digest_path)) readLines(digest_path, warn = FALSE)[[1L]] else "<absent>"
    result$copy_ok <- file.exists(result$copy_path) && identical(copied, result$record_digest)
    status$emit("copy", side, if (result$copy_ok) "ok" else "mechanism",
        sprintf("copy present: %s; its masked digest equals the record's at da: %s", file.exists(result$copy_path),
            identical(copied, result$record_digest)))
    result
}

# --------------------------------------------------------------------------------------------------------------------
# Parent: functions read from the package's own files (parsed, not copied).
# --------------------------------------------------------------------------------------------------------------------

definitionsFrom <- function(path, names) {
    env <- new.env(parent = globalenv())
    for (expression in parse(path, keep.source = FALSE)) {
        if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
            is.name(expression[[2L]]) && as.character(expression[[2L]]) %in% names) {
            eval(expression, env)
        }
    }
    missing <- setdiff(names, ls(env))
    if (length(missing)) stop("not defined in ", basename(path), ": ", paste(missing, collapse = ", "))
    env
}

# One row per differing (step, field, key), compareLane's form; `a` is the first side, `b` the second.
compareSides <- function(compare_lane, comparison, first, second) {
    rows <- compare_lane(LANE, second, first)
    if (is.null(rows) || !nrow(rows)) {
        return(data.frame(comparison = character(), lane = character(), step = character(), field = character(),
            key = character(), kind = character(), a = character(), b = character(), stringsAsFactors = FALSE))
    }
    data.frame(comparison = comparison, lane = rows$lane, step = rows$step, field = rows$field, key = rows$key,
        kind = rows$kind, a = rows$oracle, b = rows$runner, stringsAsFactors = FALSE)
}

# The classes of link (i) (ticket 6.2 R0, ruling C1 (a)); any other difference is `outside`.
classifyFrozenBase <- function(rows, declared) {
    declared_key <- paste(declared$step, declared$key, sep = "\r")
    row_key <- paste(rows$step, rows$key, sep = "\r")
    position <- match(row_key, declared_key)
    is_declared <- rows$field == "artifacts" & !is.na(position) &
        rows$a == declared$frozen[position] & rows$b == declared$expected[position]
    is_integration <- rows$field == "artifacts" & startsWith(rows$key, "integration/proteomics/") &
        (rows$a == "<absent>" | rows$b == "<absent>")
    ifelse(is_declared, DECLARED_CLASS, ifelse(is_integration, INTEGRATION_CLASS, "outside"))
}

# The (step, field, key, absent side) rows the reference allows, from g1a-reference-dia-limpa.tsv.
expectedFromReference <- function(reference) {
    rows <- list()
    for (i in seq_len(nrow(reference))) {
        text <- reference$differing_fields_and_keys[[i]]
        step <- reference$step[[i]]
        if (identical(reference$class[[i]], DECLARED_CLASS)) {
            rows[[length(rows) + 1L]] <- data.frame(step = step, key = DA_KEY, absent = "none")
            next
        }
        live_only <- sub("^artifacts: head only ([^;]+);.*$", "\\1", text)
        stored_only <- sub("^.*; frozen only (.+)$", "\\1", text)
        rows[[length(rows) + 1L]] <- data.frame(step = step, key = strsplit(live_only, ",", fixed = TRUE)[[1L]],
            absent = "a")
        rows[[length(rows) + 1L]] <- data.frame(step = step, key = strsplit(stored_only, ",", fixed = TRUE)[[1L]],
            absent = "b")
    }
    do.call(rbind, rows)
}

# --------------------------------------------------------------------------------------------------------------------
# Parent: the DA table reduced to its keys and numeric columns (text copied, not re-formatted).
# --------------------------------------------------------------------------------------------------------------------

writeDaRunner <- function(copy_path, path) {
    lines <- readLines(copy_path, warn = FALSE)
    # A sentinel keeps a trailing empty field (strsplit drops it).
    fields <- lapply(strsplit(paste0(lines, "\tEND"), "\t", fixed = TRUE), \(row) row[-length(row)])
    header <- fields[[1L]]
    widths <- lengths(fields)
    if (any(widths != length(header))) stop("the DA table's rows do not all have the header's field count")
    numeric <- header %in% NUMERIC_EXACT | Reduce(`|`, lapply(NUMERIC_PREFIXES, \(p) grepl(p, header)))
    keep <- c(match(c("Protein.Ids", "comparison"), header), which(numeric))
    if (anyNA(keep) || sum(numeric) != 20L) stop(sprintf("expected 2 key and 20 numeric columns, found %d numeric",
        sum(numeric)))
    writeLines(vapply(fields, \(row) paste(row[keep], collapse = "\t"), ""), path)
    c(rows = length(lines) - 1L, columns = length(keep))
}

# --------------------------------------------------------------------------------------------------------------------
# Parent: the environment (bare values; no upload pattern).
# --------------------------------------------------------------------------------------------------------------------

commandValue <- function(command, args, pattern = NULL) {
    path <- Sys.which(command)
    if (!nzchar(path)) return("not found")
    output <- tryCatch(processx::run(path, args, error_on_status = FALSE)$stdout, error = \(e) "")
    output <- trimws(strsplit(output, "\n", fixed = TRUE)[[1L]][1L] %||% "")
    if (!is.null(pattern)) output <- regmatches(output, regexpr(pattern, output))
    if (!length(output) || !nzchar(output)) "not found" else output
}

cpuModel <- function() {
    info <- readLines("/proc/cpuinfo", warn = FALSE)
    model <- sub("^model name\\s*:\\s*", "", grep("^model name", info, value = TRUE))
    if (length(model)) model[[1L]] else "not found"
}

internalPin <- function(root, package) {
    text <- readLines(file.path(root, "DESCRIPTION"), warn = FALSE)
    line <- grep(paste0("^\\s*", package, "="), text, value = TRUE)
    pin <- regmatches(line, regexpr("@[0-9a-f]{40}", line))
    if (length(pin)) substring(pin[[1L]], 2L) else "not found"
}

writeEnvironment <- function(path, root, commits, mode) {
    rust <- Sys.which("rustc")
    if (!nzchar(rust)) rust <- file.path(Sys.getenv("HOME"), ".cargo", "bin", "rustc")
    rust_version <- if (file.exists(rust)) {
        output <- tryCatch(processx::run(rust, "--version", error_on_status = FALSE)$stdout, error = \(e) "")
        found <- regmatches(output, regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", output))
        if (length(found)) found else "not found"
    } else "not found"
    envOr <- function(name) {
        value <- Sys.getenv(name)
        if (nzchar(value)) value else "<unset>"
    }
    values <- c(
        mode = mode,
        commits,
        cpu_model = cpuModel(),
        cpu_count = as.character(parallel::detectCores()),
        image_os = envOr("ImageOS"),
        image_version = envOr("ImageVersion"),
        c_library_version = commandValue("getconf", "GNU_LIBC_VERSION", "[0-9]+\\.[0-9]+"),
        libblas3_package = commandValue("dpkg-query", c("-W", "-f=${Version}", "libblas3")),
        liblapack3_package = commandValue("dpkg-query", c("-W", "-f=${Version}", "liblapack3")),
        blas_library = basename(normalizePath(extSoftVersion()[["BLAS"]], mustWork = FALSE)),
        lapack_library = basename(normalizePath(La_library(), mustWork = FALSE)),
        lapack_version = La_version(),
        r_version = R.version.string,
        rust_toolchain_version = rust_version,
        pin_fastduckdb = internalPin(root, "fastduckdb"),
        pin_fastlimpadpc = internalPin(root, "fastlimpadpc"),
        openblas_num_threads = envOr("OPENBLAS_NUM_THREADS"),
        omp_num_threads = envOr("OMP_NUM_THREADS")
    )
    table <- data.frame(key = names(values), value = vapply(values, cleanDetail, ""), stringsAsFactors = FALSE)
    writeTsv(table, path)
    all(nzchar(table$value)) && !any(table$value == "not found")
}

# --------------------------------------------------------------------------------------------------------------------
# Parent: the two modes.
# --------------------------------------------------------------------------------------------------------------------

setupContext <- function(out_arg, script_path) {
    root <- normalizePath(".")
    out <- file.path(root, out_arg)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    temp <- Sys.getenv("RUNNER_TEMP")
    if (!nzchar(temp)) temp <- tempdir()
    ctx <- list(root = root, out = normalizePath(out), script_path = script_path,
        trees_root = file.path(temp, "k199"), private_root = file.path(temp, "k199-private"))
    dir.create(ctx$trees_root, recursive = TRUE, showWarnings = FALSE)
    ctx$head_commit <- trimws(runGit(c("-C", root, "rev-parse", "HEAD"))$stdout)
    ctx
}

controlMain <- function(anchor_commit, base_commit, out_arg, script_path) {
    ctx <- setupContext(out_arg, script_path)
    for (sub in c("anchor", "base", "head", "oracles")) {
        dir.create(file.path(ctx$out, sub), recursive = TRUE, showWarnings = FALSE)
    }
    status <- newStatus(file.path(ctx$out, "status.tsv"))
    verdict <- list()
    check <- function(name, value, holds) {
        verdict[[length(verdict) + 1L]] <<- data.frame(check = name, value = cleanDetail(value),
            holds = isTRUE(holds), stringsAsFactors = FALSE)
    }
    commits <- c(anchor = anchor_commit, base = base_commit, head = ctx$head_commit)
    if (!all(grepl("^[0-9a-f]{40}$", commits))) {
        status$emit("arguments", "all", "mechanism", "the commits must be full 40-digit hashes")
        return(1L)
    }
    sides <- list()
    for (side in names(commits)) {
        sides[[side]] <- runSide(side, commits[[side]], if (side == "head") "freeze" else "run", ctx, status)
        check(paste0(side, "_lane_passed"), sides[[side]]$lane_ok, sides[[side]]$lane_ok)
        check(paste0(side, "_copy_digest_equals_record_at_da"), sides[[side]]$copy_ok, sides[[side]]$copy_ok)
    }

    # Byte and parsed views of the three sides' steps and of the frozen file's steps member written alone.
    helper_defs <- definitionsFrom(file.path(ctx$root, "tests/testthat/helper-gui-entry.R"), "guiEntryWriteJson")
    oracle_path <- file.path(ctx$out, "oracles", paste0(LANE, ".json"))
    frozen_steps_path <- file.path(ctx$out, "steps-frozen.json")
    if (file.exists(oracle_path)) {
        helper_defs$guiEntryWriteJson(jsonlite::read_json(oracle_path, simplifyVector = FALSE)$steps,
            frozen_steps_path)
    }
    paths <- c(vapply(sides, `[[`, "", "steps_path"), frozen = frozen_steps_path)
    hashes <- vapply(paths, sha256File, "")
    parsed <- lapply(paths, \(p) if (file.exists(p)) jsonlite::read_json(p, simplifyVector = FALSE) else NULL)
    copies <- vapply(sides, \(s) sha256File(s$copy_path), "")
    sha_table <- data.frame(
        side = c(names(sides), "frozen", names(sides), "oracle"),
        file = c(paste0("steps-", names(sides), ".json"), "steps-frozen.json", rep("da-table-raw", 3L),
            paste0("oracles/", LANE, ".json")),
        sha256 = c(hashes, copies, sha256File(oracle_path)), stringsAsFactors = FALSE)
    writeTsv(sha_table, file.path(ctx$out, "sha256.tsv"))

    present <- !vapply(parsed, is.null, logical(1L))
    bytes_equal <- function(x, y) present[[x]] && present[[y]] && identical(hashes[[x]], hashes[[y]])
    parsed_equal <- function(x, y) present[[x]] && present[[y]] && identical(parsed[[x]], parsed[[y]])
    for (pair in list(c("anchor", "base"), c("base", "head"), c("head", "frozen"))) {
        label <- paste0("steps_", pair[[1L]], "_", pair[[2L]])
        b <- bytes_equal(pair[[1L]], pair[[2L]])
        p <- parsed_equal(pair[[1L]], pair[[2L]])
        check(paste0(label, "_bytes_identical"), b, b)
        check(paste0(label, "_parsed_identical"), p, p)
        if (p && !b) status$emit("bytes", paste(pair, collapse = "+"), "mechanism", "parsed identical, bytes differ")
    }
    raw_equal <- all(copies != "<absent>") && length(unique(copies)) == 1L
    check("raw_da_anchor_base_head_identical", paste(substr(copies, 1L, 12L), collapse = ","), raw_equal)

    # gui-entry-control.tsv.
    lane_defs <- definitionsFrom(file.path(ctx$root, "packaging/verify/run-lanes.R"), c("truncate", "compareLane"))
    compare_lane <- lane_defs$compareLane
    declared <- utils::read.delim(file.path(ctx$root, RECORDS_DIR, "declared-keys.tsv"), stringsAsFactors = FALSE,
        colClasses = "character")
    reference <- utils::read.delim(file.path(ctx$root, RECORDS_DIR, "g1a-reference-dia-limpa.tsv"),
        stringsAsFactors = FALSE, colClasses = "character")
    stored <- jsonlite::read_json(file.path(sides$base$tree, "tests", "testdata", "gui-entry",
        paste0(LANE, ".json")), simplifyVector = FALSE)$steps
    none <- list()
    side_steps <- lapply(c(anchor = "anchor", base = "base", head = "head", frozen = "frozen"),
        \(s) parsed[[s]] %||% none)
    tables <- list(
        compareSides(compare_lane, "frozen_vs_base", stored, side_steps$base),
        compareSides(compare_lane, "anchor_vs_base", side_steps$anchor, side_steps$base),
        compareSides(compare_lane, "base_vs_head", side_steps$base, side_steps$head),
        compareSides(compare_lane, "head_vs_frozen_new", side_steps$head, side_steps$frozen))
    control <- do.call(rbind, tables)
    control$class <- if (nrow(control)) "outside" else character()
    link <- control$comparison == "frozen_vs_base"
    control$class[link] <- classifyFrozenBase(control[link, ], declared)
    control$changed_step <- if (nrow(control)) CHANGED_STEP else character()
    writeTsv(control, file.path(ctx$out, "gui-entry-control.tsv"))

    # Agreement 2: link (i) against the ten reference rows.
    link_rows <- control[link, ]
    expected <- expectedFromReference(reference)
    observed_keys <- sort(paste(link_rows$step, link_rows$field, link_rows$key, sep = " | "))
    expected_keys <- sort(paste(expected$step, "artifacts", expected$key, sep = " | "))
    check("frozen_vs_base_rows", sprintf("%d observed, %d expected", nrow(link_rows), nrow(expected)),
        identical(observed_keys, expected_keys))
    absent_side <- ifelse(link_rows$a == "<absent>", "a", ifelse(link_rows$b == "<absent>", "b", "none"))
    expected_absent <- expected$absent[match(paste(link_rows$step, link_rows$key), paste(expected$step, expected$key))]
    check("frozen_vs_base_absent_sides_as_reference", sum(absent_side != expected_absent, na.rm = TRUE),
        identical(absent_side, expected_absent))
    check("frozen_vs_base_outside_rows", sum(link_rows$class == "outside"), !any(link_rows$class == "outside"))
    classes <- table(factor(link_rows$class, levels = c(DECLARED_CLASS, INTEGRATION_CLASS)))
    reference_classes <- table(factor(reference$class, levels = c(DECLARED_CLASS, INTEGRATION_CLASS)))
    check("frozen_vs_base_steps_per_class",
        sprintf("%s: %d steps; %s: %d steps", DECLARED_CLASS, length(unique(link_rows$step[link_rows$class ==
            DECLARED_CLASS])), INTEGRATION_CLASS, length(unique(link_rows$step[link_rows$class == INTEGRATION_CLASS]))),
        length(unique(link_rows$step[link_rows$class == DECLARED_CLASS])) == reference_classes[[DECLARED_CLASS]] &&
            length(unique(link_rows$step[link_rows$class == INTEGRATION_CLASS])) ==
                reference_classes[[INTEGRATION_CLASS]])
    base_values <- vapply(declared$step, \(s) side_steps$base[[s]]$artifacts[[DA_KEY]] %||% "<absent>", "")
    stored_values <- vapply(declared$step, \(s) stored[[s]]$artifacts[[DA_KEY]] %||% "<absent>", "")
    check("base_da_digest_at_six_steps", paste(unique(base_values), collapse = ","), all(base_values == RUNNER_DIGEST))
    check("stored_da_digest_at_six_steps", paste(unique(stored_values), collapse = ","),
        all(stored_values == STORED_DIGEST))
    # Agreements 1 and 3: no difference between the three sides, none between the head and its frozen file.
    for (name in c("anchor_vs_base", "base_vs_head", "head_vs_frozen_new")) {
        n <- sum(control$comparison == name)
        check(paste0(name, "_rows"), n, n == 0L)
    }

    # Negative control: two planted edits on the base record.
    negative_ok <- negativeControl(compare_lane, side_steps$base, declared, file.path(ctx$out,
        "negative-control.tsv"))
    check("negative_control_two_planted_rows_outside", negative_ok, negative_ok)

    # da-runner.tsv (the head side's table) and environment.tsv.
    written <- attempt(writeDaRunner(sides$head$copy_path, file.path(ctx$out, "da-runner.tsv")))
    if (!written$ok) status$emit("da_runner", "head", "mechanism", errorDetail(written$error))
    check("da_runner_written", if (written$ok) paste(written$value, collapse = " x ") else "no", written$ok)
    complete <- writeEnvironment(file.path(ctx$out, "environment.tsv"), ctx$root,
        c(commit_anchor = anchor_commit, commit_base = base_commit, commit_head = ctx$head_commit), "control")
    check("environment_complete", complete, complete)

    verdict_table <- do.call(rbind, verdict)
    writeTsv(verdict_table, file.path(ctx$out, "verdict.tsv"))
    print(verdict_table, row.names = FALSE)
    if (all(verdict_table$holds) && negative_ok && !status$mechanism) 0L else 1L
}

negativeControl <- function(compare_lane, base_steps, declared, path) {
    if (is.null(base_steps$da) || is.null(base_steps$norm) || !length(base_steps$norm$states)) {
        writeLines("the base record lacks the steps the negative control edits", path)
        return(FALSE)
    }
    planted <- base_steps
    planted$da$artifacts[[DA_KEY]] <- digest::digest("KIT-199 planted DA value", algo = "sha256",
        serialize = FALSE)
    state <- names(planted$norm$states)[[1L]]
    planted$norm$states[[state]] <- digest::digest("KIT-199 planted state", algo = "sha256", serialize = FALSE)
    rows <- compareSides(compare_lane, "negative_control", base_steps, planted)
    rows$class <- if (nrow(rows)) classifyFrozenBase(rows, declared) else character()
    writeTsv(rows, path)
    expected <- sort(c(paste("da", "artifacts", DA_KEY, sep = " | "), paste("norm", "states", state, sep = " | ")))
    identical(sort(paste(rows$step, rows$field, rows$key, sep = " | ")), expected) && all(rows$class == "outside")
}

measureMain <- function(out_arg, script_path) {
    ctx <- setupContext(out_arg, script_path)
    dir.create(file.path(ctx$out, "measure"), recursive = TRUE, showWarnings = FALSE)
    status <- newStatus(file.path(ctx$out, "status.tsv"))
    side <- runSide("measure", ctx$head_commit, "run", ctx, status)
    verdict <- data.frame(check = c("measure_lane_passed", "measure_copy_digest_equals_record_at_da"),
        value = c(side$lane_ok, side$copy_ok), holds = c(side$lane_ok, side$copy_ok))
    written <- attempt(writeDaRunner(side$copy_path, file.path(ctx$out, "da-runner.tsv")))
    if (!written$ok) status$emit("da_runner", "measure", "mechanism", errorDetail(written$error))
    complete <- writeEnvironment(file.path(ctx$out, "environment.tsv"), ctx$root,
        c(commit_measure = ctx$head_commit), "measure")
    verdict <- rbind(verdict, data.frame(check = c("da_runner_written", "environment_complete"),
        value = c(written$ok, complete), holds = c(written$ok, complete)))
    writeTsv(verdict, file.path(ctx$out, "verdict.tsv"))
    digest_value <- side$record_digest %||% "<absent>"
    writeTsv(data.frame(quantity = c("da_masked_digest", "da_masked_digest_equals_runner_reference"),
        value = c(digest_value, identical(digest_value, RUNNER_DIGEST))), file.path(ctx$out, "measure.tsv"))
    writeTsv(data.frame(side = "measure", file = "da-table-raw", sha256 = sha256File(side$copy_path)),
        file.path(ctx$out, "sha256.tsv"))
    print(verdict, row.names = FALSE)
    if (all(verdict$holds) && !status$mechanism) 0L else 1L
}

# --------------------------------------------------------------------------------------------------------------------
# Entry.
# --------------------------------------------------------------------------------------------------------------------

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) == 5L && identical(arguments[[1L]], "child")) {
    code <- tryCatch(childMain(arguments[[2L]], normalizePath(arguments[[3L]]), arguments[[4L]], arguments[[5L]]),
        error = \(e) {
            cat(paste("child", "mechanism", oneLine(errorDetail(e)), "", sep = "\t"), "\n", sep = "",
                file = file.path(arguments[[3L]], "child-stages.tsv"), append = TRUE)
            1L
        })
    quit(save = "no", status = code)
}
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
    value = TRUE)[[1L]]))
out_arg <- arguments[[length(arguments)]]
code <- tryCatch({
    if (length(arguments) == 2L && identical(arguments[[1L]], "measure")) {
        measureMain(out_arg, script_path)
    } else if (length(arguments) == 3L) {
        controlMain(arguments[[1L]], arguments[[2L]], out_arg, script_path)
    } else {
        stop("usage: control.R <A> <B> <out> | control.R measure <out>")
    }
}, error = \(e) {
    out <- file.path(normalizePath("."), out_arg)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    cat(paste("script", "all", "mechanism", cleanDetail(errorDetail(e)), "", sep = "\t"), "\n", sep = "",
        file = file.path(out, "status.tsv"), append = TRUE)
    cat("script: mechanism\n")
    1L
})
quit(save = "no", status = code)
